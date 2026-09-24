//! The clearance index itself: a sum of activity-to-limit ratios.
//!
//! Every regulation here uses the same arithmetic. Each radionuclide's activity
//! is divided by its tabulated limit and the ratios are summed; a total below
//! one means the material meets the limits. The German regulation calls it the
//! Summenformel, the UK calls it the summation rule, the NRC calls it the sum of
//! fractions rule and Fetter calls the result a waste disposal rating. This
//! module implements it once.
//!
//! What differs between regulations, and what [`ClearanceResult`] therefore
//! records, is what happens to a nuclide that is *not* simply looked up: one
//! whose parent already accounts for it, one the table does not list, and one
//! the regulation places outside its scope entirely.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt;

use crate::error::{Error, Result};
use crate::fmt::g;
use crate::limits::{get_limit_set, limit_sets, DynamicRule, LimitSet};
use crate::material::{ActivityUnit, Material};

/// How far below its parent a daughter may sit and still count as being in
/// secular equilibrium with it. Equilibrium means equal activities, so this is
/// generous; it exists to tolerate a decayed or freshly separated inventory,
/// not to admit a trace.
pub const EQUILIBRIUM_TOLERANCE: f64 = 0.1;

/// Five years in seconds, the NRC's short-lived and long-lived boundary.
const FIVE_YEARS: f64 = 5.0 * 365.25 * 86400.0;

/// The choices [`clearance_index`] leaves to the caller.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClearanceOptions {
    /// Whether the material is activated metal, which changes some NRC limits
    /// and adds others. Default `false`.
    pub metal: bool,
    /// Whether to apply the set's catch-all limit to nuclides the table does
    /// not list. The UK regulations define one, so leaving this on is what the
    /// regulation says; turning it off shows what the listed nuclides alone
    /// contribute. Default `true`.
    pub apply_default_limit: bool,
    /// Whether to leave out daughters whose parent's limit already covers them,
    /// per the regulation's own table. Turning this off double counts them,
    /// which is conservative but not what the regulation intends. Default
    /// `true`.
    pub exclude_daughters: bool,
}

impl Default for ClearanceOptions {
    fn default() -> Self {
        ClearanceOptions {
            metal: false,
            apply_default_limit: true,
            exclude_daughters: true,
        }
    }
}

/// The outcome of assessing one material against one limit set.
///
/// The per-nuclide lists are ordered as a report would read them: `by_nuclide`,
/// `activities`, `credited` and `uncovered` largest first.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct ClearanceResult {
    /// Name of the set assessed against.
    pub limit_set: String,
    /// The sum of activity-to-limit ratios.
    pub index: f64,
    /// The value the index must stay below, normally 1.
    pub threshold: f64,
    /// Activity units the comparison was made in.
    pub units: ActivityUnit,
    /// Each nuclide's contribution to the index, largest first.
    pub by_nuclide: Vec<(String, f64)>,
    /// Each nuclide's activity in `units`, largest first.
    pub activities: Vec<(String, f64)>,
    /// The limit applied to each nuclide, after any metal, per-gram, dynamic or
    /// secular equilibrium adjustment.
    pub limits_used: Vec<(String, f64)>,
    /// Nuclides that took the set's catch-all limit because the table does not
    /// list them, sorted.
    pub defaulted: Vec<String>,
    /// Nuclide to the reason its activity was left out, which is always that a
    /// parent's limit already accounts for it in full.
    pub excluded: Vec<(String, String)>,
    /// Nuclide to the activity a parent accounted for, where the parent could
    /// only support part of it. The remainder was assessed against the
    /// nuclide's own limit, so for these the ratio in `by_nuclide` is computed
    /// from the activity less the credit. [`ClearanceResult::assessed_activity`]
    /// returns that residual directly.
    pub credited: Vec<(String, f64)>,
    /// Activity present with no limit and no catch-all, so absent from the
    /// index entirely. The number to check before trusting a comfortable index.
    pub uncovered: Vec<(String, f64)>,
    /// Nuclides present that the source explicitly places no limit on, sorted.
    pub unlimited: Vec<String>,
    /// Whether the regulation excludes this material outright, which for the
    /// UK sets means every radionuclide present is shorter lived than 100
    /// seconds.
    pub out_of_scope: bool,
    /// The material's label, carried through for reporting.
    pub material_name: String,
}

fn lookup(pairs: &[(String, f64)], name: &str) -> Option<f64> {
    pairs.iter().find(|(n, _)| n == name).map(|&(_, v)| v)
}

impl ClearanceResult {
    /// Whether the material meets this set's limits.
    pub fn clearable(&self) -> bool {
        self.out_of_scope || self.index < self.threshold
    }

    /// Total activity with no limit, in this result's units.
    pub fn uncovered_activity(&self) -> f64 {
        crate::total(self.uncovered.iter().map(|(_, v)| v))
    }

    /// Share of total activity that falls outside the sum, from 0 to 1.
    ///
    /// A large value means the index understates the inventory: activity is
    /// present that no limit in this set applies to.
    pub fn uncovered_fraction(&self) -> f64 {
        let total = crate::total(self.activities.iter().map(|(_, v)| v));
        if total > 0.0 {
            self.uncovered_activity() / total
        } else {
            0.0
        }
    }

    /// The activity actually charged against the limit for one nuclide.
    ///
    /// The same as its total activity, except where a parent accounted for part
    /// of it, in which case the remainder is what the ratio was computed from.
    pub fn assessed_activity(&self, nuclide: &str) -> f64 {
        lookup(&self.activities, nuclide).unwrap_or(0.0)
            - lookup(&self.credited, nuclide).unwrap_or(0.0)
    }

    /// The limit applied to one nuclide, if it was assessed against one.
    pub fn limit_used(&self, nuclide: &str) -> Option<f64> {
        lookup(&self.limits_used, nuclide)
    }

    /// The `count` nuclides contributing most to the index.
    pub fn dominant(&self, count: usize) -> &[(String, f64)] {
        &self.by_nuclide[..count.min(self.by_nuclide.len())]
    }
}

impl fmt::Display for ClearanceResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut header = format!(
            "{}: index {} of {}",
            self.limit_set,
            g(self.index, 4),
            g(self.threshold, 6)
        );
        if !self.material_name.is_empty() {
            header = format!("{}  {}", self.material_name, header);
        }
        let verdict = if self.out_of_scope {
            "OUT OF SCOPE (every radionuclide is short lived)"
        } else if self.clearable() {
            "CLEARABLE"
        } else {
            "NOT CLEARABLE"
        };
        write!(f, "{header}  [{verdict}]")?;
        if !self.by_nuclide.is_empty() {
            write!(
                f,
                "\n  {:<10} {:>12} {:>12} {:>8}",
                "nuclide", "activity", "limit", "share"
            )?;
            for (name, value) in self.dominant(10) {
                let share = if self.index != 0.0 {
                    value / self.index
                } else {
                    0.0
                };
                write!(
                    f,
                    "\n  {:<10} {:>12} {:>12} {:>7}",
                    name,
                    g(self.assessed_activity(name), 4),
                    g(self.limit_used(name).unwrap_or(f64::NAN), 4),
                    format!("{:.1}%", share * 100.0),
                )?;
            }
        }
        if !self.uncovered.is_empty() {
            write!(
                f,
                "\n  {} nuclide(s) with no limit, {:.2}% of total activity",
                self.uncovered.len(),
                self.uncovered_fraction() * 100.0
            )?;
        }
        Ok(())
    }
}

/// The limit actually applied to each nuclide for one material.
struct EffectiveLimits<'a> {
    set: &'a LimitSet,
    metal: bool,
    /// nCi/g to Ci/m3, when the set has per-gram limits.
    per_gram_factor: f64,
    dynamic: HashMap<String, f64>,
    overrides: HashMap<String, f64>,
}

impl<'a> EffectiveLimits<'a> {
    fn new(material: &Material, set: &'a LimitSet, metal: bool) -> Result<Self> {
        let mut dynamic = HashMap::new();
        if let Some(DynamicRule::NrcShortLivedClassA) = set.dynamic_rule() {
            for name in material.nuclides() {
                if let Some(t) = material.decay_data().half_life_canonical(name)? {
                    if t < FIVE_YEARS {
                        dynamic.insert(name.to_string(), 700.0);
                    }
                }
            }
        }
        // 1 nCi/g is 1e-9 Ci/g, and multiplying by g/cm3 then by 1e6 cm3/m3
        // gives Ci/m3.
        let per_gram_factor = if set.limits_per_gram().is_empty() {
            0.0
        } else {
            1e-9 * material.mass_density()? * 1e6
        };
        Ok(EffectiveLimits {
            set,
            metal,
            per_gram_factor,
            dynamic,
            overrides: HashMap::new(),
        })
    }

    /// The same precedence as building the map in order: table, then the
    /// dynamic rule where the table is silent, then metal rows, then per-gram
    /// rows, then any secular equilibrium variant.
    fn get(&self, name: &str) -> Option<f64> {
        if let Some(&v) = self.overrides.get(name) {
            return Some(v);
        }
        if let Some(&v) = self.set.limits_per_gram().get(name) {
            return Some(v * self.per_gram_factor);
        }
        if self.metal {
            if let Some(&v) = self.set.metal_overrides().get(name) {
                return Some(v);
            }
        }
        if let Some(&v) = self.set.limits().get(name) {
            return Some(v);
        }
        self.dynamic.get(name).copied()
    }
}

/// Whether the daughters present are near enough equilibrium with the parent.
///
/// In secular equilibrium a daughter's activity equals its parent's. This asks
/// whether every tabulated daughter present is within [`EQUILIBRIUM_TOLERANCE`]
/// of that, which is what earns a parent the limit that assumes equilibrium.
fn in_equilibrium(activities: &BTreeMap<String, f64>, parent: &str, daughters: &[&str]) -> bool {
    let parent_activity = activities.get(parent).copied().unwrap_or(0.0);
    if parent_activity <= 0.0 {
        return false;
    }
    daughters.iter().all(|d| {
        activities.get(*d).copied().unwrap_or(0.0) >= EQUILIBRIUM_TOLERANCE * parent_activity
    })
}

/// Apply a whole-material scope test such as the UK 100 second rule.
///
/// The regulation says a substance is not radioactive material where *none* of
/// the radionuclides it contains has a half-life exceeding the threshold. That
/// is a test on the material, so a short-lived nuclide cannot be dropped from
/// the sum individually.
fn out_of_scope(material: &Material, set: &LimitSet) -> Result<bool> {
    let Some(scope) = set.min_half_life_scope() else {
        return Ok(false);
    };
    let mut any = false;
    for name in material.nuclides() {
        if let Some(t) = material.decay_data().half_life_canonical(name)? {
            any = true;
            if t > scope {
                return Ok(false);
            }
        }
    }
    // With no radionuclides at all there is nothing for the rule to be true
    // of, and saying out of scope would force clearable regardless of the
    // index, which matters when activities were supplied directly for a
    // nuclide the tables call stable.
    Ok(any)
}

/// Work out what a parent's limit already accounts for.
///
/// The regulations mark a parent whose tabulated value "already takes into
/// account the daughter radionuclides present". That claim rests on two things
/// being true, and this checks both rather than assuming them.
///
/// The parent must actually be limited by this set. A parent with no row here
/// contributes nothing to the index, so crediting its daughters against it
/// would remove them from the sum on the strength of a limit that does not
/// exist. Nine of the shipped sets contain such a pair, StrlSchV_soil Np-237
/// and Pa-233 among them.
///
/// The credit is also bounded by how much daughter activity the parent can
/// support. Secular equilibrium means the daughter's activity equals its
/// parent's, so only that much is accounted for, and any excess stays in the
/// sum. Without this a trace of Sr-90 would delete an arbitrarily large Y-90
/// activity from the index.
///
/// Fills `limits.overrides` and returns `(excluded, credited)`.
fn resolve_equilibrium(
    set: &LimitSet,
    activities: &BTreeMap<String, f64>,
    limits: &mut EffectiveLimits,
    default_applies: bool,
) -> (Vec<(String, String)>, HashMap<String, f64>) {
    let present = |name: &str| activities.get(name).is_some_and(|&v| v > 0.0);
    let mut excluded: Vec<(String, String)> = Vec::new();
    let mut credited: HashMap<String, f64> = HashMap::new();

    let parents: BTreeSet<&String> = set
        .secular_equilibrium()
        .keys()
        .chain(set.secular_equilibrium_sec().keys())
        .collect();
    for parent in parents {
        let parent = parent.as_str();
        if !present(parent) {
            continue;
        }
        // A nuclide already accounted for by someone else cannot also account
        // for others. Without this a mutually referencing pair excludes both of
        // them and the whole activity disappears.
        if excluded.iter().any(|(n, _)| n == parent) {
            continue;
        }

        // A parent can carry two published values with two daughter lists: a
        // "+" value covering a few short lived progeny, and a much stricter
        // whole chain "sec" value covering the lot. Applying the "+" value while
        // excluding the "sec" list charges the parent against a limit that
        // accounts for only part of what was removed, which is how natural
        // uranium came to clear at fifty times the set's own limit. Pick the
        // value whose own list matches what is actually there.
        let plus: &[String] = set
            .secular_equilibrium()
            .get(parent)
            .map_or(&[], Vec::as_slice);
        let chain: &[String] = set
            .secular_equilibrium_sec()
            .get(parent)
            .map_or(&[], Vec::as_slice);
        let chain_limit = set.limits().get(&format!("{parent}_sec")).copied();
        let beyond_plus = chain.iter().any(|d| present(d) && !plus.contains(d));
        let daughters = match chain_limit {
            Some(value) if beyond_plus => {
                limits.overrides.insert(parent.to_string(), value);
                chain
            }
            _ => plus,
        };

        // An explicitly unlimited parent is NOT a parent that can account for
        // anything. The source placing no limit on it says nothing about its
        // daughters, and crediting them against a limit that does not exist
        // deletes their activity outright.
        let parent_limited = limits.get(parent).is_some()
            || set.limits_secular_equilibrium().contains_key(parent)
            || default_applies;
        if !parent_limited {
            limits.overrides.remove(parent);
            continue;
        }

        let found: Vec<&str> = daughters
            .iter()
            .map(String::as_str)
            .filter(|d| present(d) && *d != parent)
            .collect();
        if found.is_empty() {
            limits.overrides.remove(parent);
            continue;
        }
        if !limits.overrides.contains_key(parent) {
            if let Some(&variant) = set.limits_secular_equilibrium().get(parent) {
                let plain = limits.get(parent);
                // A variant stricter than the plain row can be applied as soon
                // as any tabulated daughter appears, since that is the
                // conservative direction. A more LENIENT variant is different:
                // IRR 2017 gives U-240 as 0.01 Bq/g plain and 100 marked, so
                // switching on the mere presence of a daughter buys a ten
                // thousand fold relief for a trace twelve orders of magnitude
                // below equilibrium. The lenient value is only earned when the
                // daughters really are in secular equilibrium, and otherwise the
                // parent is held to its own row with nothing excluded.
                if plain.map_or(true, |p| variant <= p)
                    || in_equilibrium(activities, parent, &found)
                {
                    limits.overrides.insert(parent.to_string(), variant);
                } else {
                    continue;
                }
            }
        }

        let parent_activity = activities.get(parent).copied().unwrap_or(0.0);
        for daughter in found {
            let activity = activities[daughter];
            let supported = activity.min(parent_activity);
            if supported >= activity {
                let reason = format!(
                    "in secular equilibrium with {parent}, whose limit already accounts for it"
                );
                match excluded.iter_mut().find(|(n, _)| n == daughter) {
                    Some(entry) => entry.1 = reason,
                    None => excluded.push((daughter.to_string(), reason)),
                }
            } else if supported > 0.0 {
                // A daughter can have more than one parent present. Each parent
                // supports at most its own activity, and taking the largest
                // rather than whichever came last keeps the result independent
                // of iteration order.
                let entry = credited.entry(daughter.to_string()).or_insert(0.0);
                *entry = entry.max(supported);
            }
        }
    }
    (excluded, credited)
}

fn largest_first(map: impl IntoIterator<Item = (String, f64)>) -> Vec<(String, f64)> {
    let mut pairs: Vec<(String, f64)> = map.into_iter().collect();
    // Stable, so ties stay in name order.
    pairs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
    pairs
}

/// Assess a material against one set of clearance limits.
///
/// Returns a [`ClearanceResult`] carrying the index and everything needed to
/// judge it, including any activity that fell outside the sum. An error if a
/// volumetric set is used and the material has no density, or if a nuclide the
/// set needs a half-life for is not in the decay data.
pub fn clearance_index(
    material: &Material,
    limit_set: &LimitSet,
    options: ClearanceOptions,
) -> Result<ClearanceResult> {
    let activities = material.activities(limit_set.units())?;
    let mut limits = EffectiveLimits::new(material, limit_set, options.metal)?;
    let default_applies = options.apply_default_limit && limit_set.default_limit().is_some();
    let (excluded, credited) = if options.exclude_daughters {
        resolve_equilibrium(limit_set, &activities, &mut limits, default_applies)
    } else {
        (Vec::new(), HashMap::new())
    };

    let mut ratios: Vec<(String, f64)> = Vec::new();
    let mut used: Vec<(String, f64)> = Vec::new();
    let mut uncovered: Vec<(String, f64)> = Vec::new();
    let mut defaulted: Vec<String> = Vec::new();
    let mut unlimited: Vec<String> = Vec::new();

    for (name, &total) in &activities {
        if total <= 0.0 || excluded.iter().any(|(n, _)| n == name) {
            continue;
        }
        // Only the part of a daughter's activity its parent can support is
        // accounted for by the parent's limit. The rest is assessed normally.
        let activity = total - credited.get(name).copied().unwrap_or(0.0);
        if activity <= 0.0 {
            continue;
        }
        if limit_set.unlimited().contains(name) {
            unlimited.push(name.clone());
            continue;
        }
        let limit = match limits.get(name) {
            Some(limit) => limit,
            None if default_applies => {
                defaulted.push(name.clone());
                limit_set
                    .default_limit()
                    .expect("default_applies checks it")
            }
            None => {
                uncovered.push((name.clone(), activity));
                continue;
            }
        };
        if limit.is_nan() || limit <= 0.0 {
            // The set validates positivity, but a per-gram limit is computed at
            // evaluation time. Reporting the nuclide as uncovered is the honest
            // outcome; dividing by it would give infinity, and skipping it
            // silently would drop the activity.
            uncovered.push((name.clone(), activity));
            continue;
        }
        ratios.push((name.clone(), activity / limit));
        used.push((name.clone(), limit));
    }

    Ok(ClearanceResult {
        limit_set: limit_set.name().to_string(),
        index: crate::total(ratios.iter().map(|(_, v)| v)),
        threshold: limit_set.threshold(),
        units: limit_set.units(),
        by_nuclide: largest_first(ratios),
        activities: largest_first(activities),
        limits_used: used,
        defaulted,
        excluded,
        credited: largest_first(credited),
        uncovered: largest_first(uncovered),
        unlimited,
        out_of_scope: out_of_scope(material, limit_set)?,
        material_name: material.name.clone(),
    })
}

/// Assess a material against many limit sets at once.
///
/// `names` defaults to every registered set. Sets whose units the material
/// cannot supply, such as a volumetric set for a material with no density, are
/// skipped rather than failing, so one missing density does not hide every
/// result that does not need it. Results are in the order of `names`.
pub fn clearance_indices(
    material: &Material,
    names: Option<&[&str]>,
    options: ClearanceOptions,
) -> Result<Vec<ClearanceResult>> {
    let all;
    let names: Vec<&str> = match names {
        Some(names) => names.to_vec(),
        None => {
            all = limit_sets(None);
            all.iter().map(String::as_str).collect()
        }
    };
    let mut results = Vec::with_capacity(names.len());
    for name in names {
        match clearance_index(material, &*get_limit_set(name)?, options) {
            Ok(result) => results.push(result),
            Err(Error::InsufficientData(_)) => continue,
            Err(e) => return Err(e),
        }
    }
    Ok(results)
}

/// The limit sets this material already meets, best margin first.
///
/// `names` defaults to every registered set.
pub fn clearable_routes(
    material: &Material,
    names: Option<&[&str]>,
    options: ClearanceOptions,
) -> Result<Vec<String>> {
    let mut passing: Vec<ClearanceResult> = clearance_indices(material, names, options)?
        .into_iter()
        .filter(|r| r.clearable())
        .collect();
    passing.sort_by(|a, b| {
        (a.index / a.threshold)
            .partial_cmp(&(b.index / b.threshold))
            .unwrap_or(Ordering::Equal)
    });
    Ok(passing.into_iter().map(|r| r.limit_set).collect())
}
