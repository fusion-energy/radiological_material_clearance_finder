//! Limit sets: the regulatory tables a material is measured against.
//!
//! Every limit set is data, compiled in from the JSON under `data/limits`.
//! Adding a jurisdiction means adding a file and a build script, never a branch
//! in the calculation. Only two things in these regulations cannot be expressed
//! as a table, and both are named explicitly rather than worked around:
//!
//! * limits given per gram of material (the NRC's nCi/g transuranic entries),
//!   which need the material's density before they can be compared against
//!   Ci/m3, and
//! * limits that depend on the material's own nuclides (the NRC Class A rule
//!   that anything with a half-life under five years takes a 700 Ci/m3 limit),
//!   held in [`DynamicRule`].

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fmt;
use std::str::FromStr;
use std::sync::{Arc, OnceLock, RwLock};

use serde::Deserialize;

use crate::error::{Error, Result};
use crate::material::ActivityUnit;
use crate::nuclide;

/// The shipped tables, in the order they are registered.
const SHIPPED: [(&str, &str); 6] = [
    (
        "de_strlschv.json",
        include_str!("../data/limits/de_strlschv.json"),
    ),
    ("eu_bss.json", include_str!("../data/limits/eu_bss.json")),
    ("iaea.json", include_str!("../data/limits/iaea.json")),
    (
        "uk_epr16.json",
        include_str!("../data/limits/uk_epr16.json"),
    ),
    (
        "uk_irr17.json",
        include_str!("../data/limits/uk_irr17.json"),
    ),
    ("us.json", include_str!("../data/limits/us.json")),
];

/// Units a limit set may be expressed in.
pub const LIMIT_UNITS: [ActivityUnit; 3] = [
    ActivityUnit::BqPerG,
    ActivityUnit::CiPerM3,
    ActivityUnit::Bq,
];

/// A rule that depends on the material rather than only on the table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DynamicRule {
    /// 10 CFR 61.55 Table 2 row 1: all nuclides with a half-life under five
    /// years. The Class A column limits their total to 700 Ci/m3, so every
    /// such nuclide that does not have its own row takes that limit.
    NrcShortLivedClassA,
}

impl DynamicRule {
    /// Every rule.
    pub const ALL: [DynamicRule; 1] = [DynamicRule::NrcShortLivedClassA];

    /// The name a limit set's `dynamic_rule` field refers to the rule by.
    pub fn as_str(self) -> &'static str {
        match self {
            DynamicRule::NrcShortLivedClassA => "nrc_short_lived_class_a",
        }
    }
}

impl FromStr for DynamicRule {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        DynamicRule::ALL
            .into_iter()
            .find(|r| r.as_str() == text)
            .ok_or_else(|| {
                let known: Vec<&str> = DynamicRule::ALL.iter().map(|r| r.as_str()).collect();
                Error::Invalid(format!(
                    "unknown dynamic rule '{text}', expected one of {}",
                    known.join(", ")
                ))
            })
    }
}

fn default_threshold() -> f64 {
    1.0
}

/// The fields of a limit set as written, before validation.
///
/// This is the JSON form of the shipped tables, and the way to build a
/// site-specific or draft set: fill one in and pass it to [`LimitSet::new`],
/// which canonicalises every nuclide name and rejects values that cannot mean
/// anything. See [`LimitSet`] for what each field means.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LimitSetData {
    pub name: String,
    pub label: String,
    pub units: String,
    pub limits: BTreeMap<String, f64>,
    #[serde(default)]
    pub jurisdiction: String,
    #[serde(default)]
    pub default_limit: Option<f64>,
    #[serde(default)]
    pub unlimited: Vec<String>,
    #[serde(default)]
    pub secular_equilibrium: BTreeMap<String, Vec<String>>,
    #[serde(default)]
    pub secular_equilibrium_sec: BTreeMap<String, Vec<String>>,
    #[serde(default)]
    pub limits_secular_equilibrium: BTreeMap<String, f64>,
    #[serde(default)]
    pub metal_overrides: BTreeMap<String, f64>,
    #[serde(default)]
    pub limits_per_gram: BTreeMap<String, f64>,
    #[serde(default)]
    pub limits_upper: BTreeMap<String, f64>,
    #[serde(default)]
    pub dynamic_rule: Option<String>,
    #[serde(default)]
    pub min_half_life_scope: Option<f64>,
    #[serde(default = "default_threshold")]
    pub threshold: f64,
    #[serde(default)]
    pub source: String,
    #[serde(default)]
    pub url: String,
    #[serde(default)]
    pub retrieved: String,
    #[serde(default)]
    pub notes: String,
}

impl LimitSetData {
    /// A set with only the required fields filled in, and the rest empty.
    pub fn new(
        name: impl Into<String>,
        label: impl Into<String>,
        units: impl Into<String>,
        limits: BTreeMap<String, f64>,
    ) -> LimitSetData {
        LimitSetData {
            name: name.into(),
            label: label.into(),
            units: units.into(),
            limits,
            jurisdiction: String::new(),
            default_limit: None,
            unlimited: Vec::new(),
            secular_equilibrium: BTreeMap::new(),
            secular_equilibrium_sec: BTreeMap::new(),
            limits_secular_equilibrium: BTreeMap::new(),
            metal_overrides: BTreeMap::new(),
            limits_per_gram: BTreeMap::new(),
            limits_upper: BTreeMap::new(),
            dynamic_rule: None,
            min_half_life_scope: None,
            threshold: 1.0,
            source: String::new(),
            url: String::new(),
            retrieved: String::new(),
            notes: String::new(),
        }
    }
}

/// A named, validated table of per-nuclide activity limits.
///
/// * `name`: stable identifier, such as `UK_EPR16_out_of_scope`.
/// * `label`: human readable description for reports.
/// * `jurisdiction`: issuing authority, such as `UK` or `Germany`.
/// * `units`: `Bq/g`, `Ci/m3`, or `Bq` for a total activity limit, which needs
///   the material's mass.
/// * `limits`: limit per canonical nuclide name.
/// * `default_limit`: limit applied to a nuclide absent from `limits`, or none
///   where the regulation has no catch-all. Three UK sets define one: 0.01 Bq/g
///   for `UK_EPR16_out_of_scope` and `UK_IRR17_notification`, and 0.1 Bq/g for
///   `UK_IRR17_registration`.
/// * `unlimited`: nuclides the source explicitly places no limit on. They are
///   covered, and contribute nothing, which is different from being absent.
/// * `secular_equilibrium`: parent to the daughters whose contribution the
///   parent's limit already includes, the "+" rows.
/// * `secular_equilibrium_sec`: parent to the daughters covered by its whole
///   chain "sec" value, usually a much longer list with a much stricter limit,
///   stored in `limits` under a `_sec` key. UK EPR 2016 gives U-238 three
///   progeny at 1 Bq/g under "U-238+" and fourteen at 0.01 Bq/g under
///   "U-238sec".
/// * `limits_secular_equilibrium`: the limit for a parent listed twice, plain
///   and marked "+", when its daughters are actually present. StrlSchV gives
///   Th-232 as 10 Bq/g plain and 0.01 Bq/g marked.
/// * `metal_overrides`: limits replacing or adding to `limits` for activated
///   metal.
/// * `limits_per_gram`: limits in nCi/g, converted using the material density.
/// * `limits_upper`: upper end of a limit given as a range, kept for
///   reference. `limits` holds the conservative lower end.
/// * `dynamic_rule`: a rule that depends on the material's own nuclides.
/// * `min_half_life_scope`: half-life in seconds below which, if *every*
///   radionuclide present falls under it, the material is outside the
///   regulation altogether. A whole-material test, not a per-nuclide filter.
/// * `threshold`: index at or above which the material fails, normally 1.
/// * `source`, `url`, `retrieved`, `notes`: provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct LimitSet {
    name: String,
    label: String,
    units: ActivityUnit,
    limits: BTreeMap<String, f64>,
    jurisdiction: String,
    default_limit: Option<f64>,
    unlimited: Vec<String>,
    secular_equilibrium: BTreeMap<String, Vec<String>>,
    secular_equilibrium_sec: BTreeMap<String, Vec<String>>,
    limits_secular_equilibrium: BTreeMap<String, f64>,
    metal_overrides: BTreeMap<String, f64>,
    limits_per_gram: BTreeMap<String, f64>,
    limits_upper: BTreeMap<String, f64>,
    dynamic_rule: Option<DynamicRule>,
    min_half_life_scope: Option<f64>,
    threshold: f64,
    source: String,
    url: String,
    retrieved: String,
    notes: String,
}

/// Normalise nuclide keys, refusing to let two of them collapse silently.
///
/// `U240` and `U-240+` both normalise to `U240`, and Schedule 7 really does
/// publish both, ten thousand apart. Taking whichever happened to come last
/// would make the set that lenient depending on nothing but map order, so a
/// collision on different values is an error. The "+" and "sec" variants
/// belong in `limits_secular_equilibrium` and under a `_sec` key.
fn canonical_limits(raw: &BTreeMap<String, f64>, label: &str) -> Result<BTreeMap<String, f64>> {
    let mut out: BTreeMap<String, f64> = BTreeMap::new();
    for (name, &value) in raw {
        // "sec" entries keep an explicit key so that the whole-chain value and
        // the tabulated-daughters value stay distinguishable.
        let key = match name.strip_suffix("_sec") {
            Some(base) => format!("{}_sec", nuclide::normalise(base)?),
            None => nuclide::normalise(name)?,
        };
        if let Some(&earlier) = out.get(&key) {
            if earlier != value && !(earlier.is_nan() && value.is_nan()) {
                return Err(Error::Invalid(format!(
                    "{label}: '{name}' and an earlier key both normalise to '{key}' but \
                     give different values, {earlier:?} and {value:?}. Whichever came \
                     last would silently win. Put the secular equilibrium variant in \
                     limits_secular_equilibrium, or under a _sec key."
                )));
            }
        }
        out.insert(key, value);
    }
    Ok(out)
}

fn canonical_daughters(
    raw: &BTreeMap<String, Vec<String>>,
) -> Result<BTreeMap<String, Vec<String>>> {
    raw.iter()
        .map(|(parent, daughters)| {
            Ok((
                nuclide::normalise(parent)?,
                daughters
                    .iter()
                    .map(|d| nuclide::normalise(d))
                    .collect::<Result<_>>()?,
            ))
        })
        .collect()
}

impl LimitSet {
    /// Canonicalise and validate a set.
    ///
    /// Limits are required to be positive and finite. A limit of zero means
    /// "no activity of this is permitted", but the sum of fractions cannot
    /// express that, and an infinite limit divides every activity to nothing;
    /// both would read as no limit at all. A nuclide the source places no limit
    /// on belongs in `unlimited`.
    pub fn new(data: LimitSetData) -> Result<LimitSet> {
        let name = data.name;
        let label_of = |field: &str| format!("{name}.{field}");
        let mut limits = canonical_limits(&data.limits, &label_of("limits"))?;
        let metal_overrides =
            canonical_limits(&data.metal_overrides, &label_of("metal_overrides"))?;
        let limits_per_gram =
            canonical_limits(&data.limits_per_gram, &label_of("limits_per_gram"))?;
        let limits_upper = canonical_limits(&data.limits_upper, &label_of("limits_upper"))?;
        let limits_secular_equilibrium = canonical_limits(
            &data.limits_secular_equilibrium,
            &label_of("limits_secular_equilibrium"),
        )?;

        // A "sec" row is the value for a parent taken with its whole chain in
        // secular equilibrium. Its "_sec" key can never match a nuclide in a
        // material, so for a nuclide whose only row is the "sec" one it has to
        // become the plain limit too, or UK_IRR17_natural would give natural
        // uranium no limit at all rather than the published 1 Bq/g.
        let promoted: Vec<(String, f64)> = limits
            .iter()
            .filter_map(|(k, &v)| k.strip_suffix("_sec").map(|base| (base.to_string(), v)))
            .filter(|(base, _)| !limits.contains_key(base))
            .collect();
        limits.extend(promoted);

        let units: ActivityUnit = data
            .units
            .parse()
            .ok()
            .filter(|u| LIMIT_UNITS.contains(u))
            .ok_or_else(|| {
                Error::Invalid(format!(
                    "{name}: units '{}' is not one of Bq, Bq/g, Ci/m3",
                    data.units
                ))
            })?;
        if !(data.threshold > 0.0 && data.threshold.is_finite()) {
            return Err(Error::Invalid(format!(
                "{name}: threshold must be positive, got {:?}",
                data.threshold
            )));
        }
        if let Some(scope) = data.min_half_life_scope {
            if !(scope > 0.0 && scope.is_finite()) {
                return Err(Error::Invalid(format!(
                    "{name}: min_half_life_scope is {scope:?}, but it must be positive \
                     and finite. It is not a ratio, it decides whether the material is \
                     in scope at all, so an infinite value would put every material out \
                     of scope and report it clearable whatever its index."
                )));
            }
        }
        if let Some(default) = data.default_limit {
            if !(default > 0.0 && default.is_finite()) {
                return Err(Error::Invalid(format!(
                    "{name}: default_limit must be positive, got {default:?}"
                )));
            }
        }
        for (field, map) in [
            ("limits", &limits),
            ("metal_overrides", &metal_overrides),
            ("limits_per_gram", &limits_per_gram),
            ("limits_secular_equilibrium", &limits_secular_equilibrium),
        ] {
            for (nuclide_name, &value) in map {
                if !(value > 0.0 && value.is_finite()) {
                    return Err(Error::Invalid(format!(
                        "{name}: {field}[{nuclide_name}] is {value:?}, but a limit must be \
                         positive and finite. Zero cannot be expressed as a ratio, and an \
                         infinite limit divides every activity to nothing, so both would be \
                         read as no limit at all. A nuclide the source places no limit on \
                         belongs in unlimited."
                    )));
                }
            }
        }

        Ok(LimitSet {
            label: data.label,
            units,
            limits,
            jurisdiction: data.jurisdiction,
            default_limit: data.default_limit,
            unlimited: data
                .unlimited
                .iter()
                .map(|n| nuclide::normalise(n))
                .collect::<Result<_>>()?,
            secular_equilibrium: canonical_daughters(&data.secular_equilibrium)?,
            secular_equilibrium_sec: canonical_daughters(&data.secular_equilibrium_sec)?,
            limits_secular_equilibrium,
            metal_overrides,
            limits_per_gram,
            limits_upper,
            dynamic_rule: data.dynamic_rule.as_deref().map(str::parse).transpose()?,
            min_half_life_scope: data.min_half_life_scope,
            threshold: data.threshold,
            source: data.source,
            url: data.url,
            retrieved: data.retrieved,
            notes: data.notes,
            name,
        })
    }

    /// Stable identifier, such as `UK_EPR16_out_of_scope`.
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Human readable description.
    pub fn label(&self) -> &str {
        &self.label
    }
    /// Activity units the limits are in.
    pub fn units(&self) -> ActivityUnit {
        self.units
    }
    /// Limit per canonical nuclide name, with whole chain values under `_sec`.
    pub fn limits(&self) -> &BTreeMap<String, f64> {
        &self.limits
    }
    /// Issuing authority.
    pub fn jurisdiction(&self) -> &str {
        &self.jurisdiction
    }
    /// The catch-all limit for unlisted nuclides, if the regulation has one.
    pub fn default_limit(&self) -> Option<f64> {
        self.default_limit
    }
    /// Nuclides the source explicitly places no limit on.
    pub fn unlimited(&self) -> &[String] {
        &self.unlimited
    }
    /// The "+" rows: parent to the daughters its limit covers.
    pub fn secular_equilibrium(&self) -> &BTreeMap<String, Vec<String>> {
        &self.secular_equilibrium
    }
    /// The "sec" rows: parent to the daughters its whole chain value covers.
    pub fn secular_equilibrium_sec(&self) -> &BTreeMap<String, Vec<String>> {
        &self.secular_equilibrium_sec
    }
    /// The limit for a parent listed twice, applied when its daughters are present.
    pub fn limits_secular_equilibrium(&self) -> &BTreeMap<String, f64> {
        &self.limits_secular_equilibrium
    }
    /// Limits replacing or adding to `limits` for activated metal.
    pub fn metal_overrides(&self) -> &BTreeMap<String, f64> {
        &self.metal_overrides
    }
    /// Limits in nCi/g.
    pub fn limits_per_gram(&self) -> &BTreeMap<String, f64> {
        &self.limits_per_gram
    }
    /// Upper end of limits the source gives as a range.
    pub fn limits_upper(&self) -> &BTreeMap<String, f64> {
        &self.limits_upper
    }
    /// The rule that depends on the material's own nuclides, if any.
    pub fn dynamic_rule(&self) -> Option<DynamicRule> {
        self.dynamic_rule
    }
    /// The whole-material short half-life scope test, in seconds.
    pub fn min_half_life_scope(&self) -> Option<f64> {
        self.min_half_life_scope
    }
    /// The index value at or above which a material fails.
    pub fn threshold(&self) -> f64 {
        self.threshold
    }
    /// The regulation this set is taken from.
    pub fn source(&self) -> &str {
        &self.source
    }
    /// Where the source was retrieved from.
    pub fn url(&self) -> &str {
        &self.url
    }
    /// When the source was retrieved.
    pub fn retrieved(&self) -> &str {
        &self.retrieved
    }
    /// Anything else worth knowing about the set.
    pub fn notes(&self) -> &str {
        &self.notes
    }

    /// Daughters whose activity this set's limit for `parent` already covers.
    pub fn daughters_of(&self, parent: &str) -> &[String] {
        self.secular_equilibrium
            .get(parent)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    /// Every nuclide this set names, whether limited or explicitly unlimited.
    ///
    /// The `_sec` keys are whole chain variants of a nuclide already counted,
    /// not nuclides of their own, so they are not included.
    pub fn covered_nuclides(&self) -> BTreeSet<&str> {
        self.limits
            .keys()
            .chain(self.unlimited.iter())
            .chain(self.limits_per_gram.keys())
            .map(String::as_str)
            .filter(|n| !n.ends_with("_sec"))
            .collect()
    }

    /// The number of nuclides limited, not counting whole chain variants.
    pub fn len(&self) -> usize {
        self.limits.keys().filter(|n| !n.ends_with("_sec")).count()
    }

    /// Whether the set limits no nuclides at all.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl fmt::Display for LimitSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "LimitSet('{}', {} nuclides, {})",
            self.name,
            self.len(),
            self.units
        )
    }
}

// ----------------------------------------------------------------------
// registry
// ----------------------------------------------------------------------

struct Registry {
    sets: BTreeMap<String, Arc<LimitSet>>,
    /// Names that came from the shipped data, so replacing one can be told
    /// apart from a caller re-registering their own set.
    shipped: HashSet<String>,
}

/// Build limit sets from the JSON form the shipped tables use.
///
/// Keys beginning with an underscore are provenance for the file itself.
/// Any other unrecognised key is an error: dropping it silently would hide a
/// build script typo or a field renamed on one side only.
pub fn parse_limit_sets_json(text: &str) -> Result<Vec<LimitSet>> {
    #[derive(Deserialize)]
    struct File {
        #[serde(default)]
        sets: Vec<LimitSetData>,
    }
    let file: File = serde_json::from_str(text)
        .map_err(|e| Error::Invalid(format!("limit set data is not valid: {e}")))?;
    file.sets.into_iter().map(LimitSet::new).collect()
}

fn registry() -> &'static RwLock<Registry> {
    static REGISTRY: OnceLock<RwLock<Registry>> = OnceLock::new();
    REGISTRY.get_or_init(|| {
        let mut sets = BTreeMap::new();
        let mut shipped = HashSet::new();
        for (file, text) in SHIPPED {
            // Compiled in and covered by the tests, so a failure is a broken
            // build rather than something a caller could handle.
            let parsed = parse_limit_sets_json(text)
                .unwrap_or_else(|e| panic!("embedded limits/{file} does not load: {e}"));
            for set in parsed {
                shipped.insert(set.name.clone());
                sets.insert(set.name.clone(), Arc::new(set));
            }
        }
        RwLock::new(Registry { sets, shipped })
    })
}

/// Names of the available limit sets, sorted.
///
/// `jurisdiction` restricts to one issuing authority, such as `UK`, matched
/// case-insensitively.
pub fn limit_sets(jurisdiction: Option<&str>) -> Vec<String> {
    let registry = registry().read().expect("limit set registry poisoned");
    registry
        .sets
        .values()
        .filter(|s| jurisdiction.map_or(true, |j| s.jurisdiction.eq_ignore_ascii_case(j)))
        .map(|s| s.name.clone())
        .collect()
}

/// Look up a limit set by name.
///
/// An unknown name is an error listing the ones that are registered.
pub fn get_limit_set(name: &str) -> Result<Arc<LimitSet>> {
    let registry = registry().read().expect("limit set registry poisoned");
    registry.sets.get(name).cloned().ok_or_else(|| {
        let known: Vec<&str> = registry.sets.keys().map(String::as_str).collect();
        Error::UnknownLimitSet(format!(
            "unknown limit set '{name}'. Available: {}",
            known.join(", ")
        ))
    })
}

/// Add a limit set at runtime, for a site-specific or draft table.
///
/// Replaces any set of the same name. Returns `true` when the replaced set came
/// from the shipped regulatory data, which callers should warn about: results
/// record only the name, so shadowing a published table by accident would be
/// hard to spot.
pub fn register_limit_set(limit_set: LimitSet) -> bool {
    let mut registry = registry().write().expect("limit set registry poisoned");
    let shadowed = registry.shipped.contains(&limit_set.name);
    registry
        .sets
        .insert(limit_set.name.clone(), Arc::new(limit_set));
    shadowed
}
