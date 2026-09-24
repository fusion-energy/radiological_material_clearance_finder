//! The material a clearance index is computed for.
//!
//! A material here is just a nuclide inventory plus enough information to put
//! it on a per-gram or per-volume basis. It deliberately knows nothing about
//! geometry, temperature or cross sections, so it can be built from a map typed
//! by hand as readily as from a transport or depletion code's own material:
//! yamc and yani hand over their atom densities with
//! [`Material::from_atom_densities`].

use std::collections::BTreeMap;
use std::fmt;
use std::sync::Arc;

use crate::decay::{DecayData, AVOGADRO, BECQUEREL_PER_CURIE};
use crate::error::{Error, Result};
use crate::fmt::g;
use crate::nuclide;

/// 1 atom/barn-cm is this many atoms per cubic centimetre.
const ATOMS_PER_BARN_CM: f64 = 1e24;

/// Activity units a material can report.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize)]
pub enum ActivityUnit {
    /// Total activity in becquerel. Needs the material's mass.
    #[serde(rename = "Bq")]
    Bq,
    /// Specific activity, the unit most clearance tables use.
    #[serde(rename = "Bq/g")]
    BqPerG,
    /// Specific activity per kilogram.
    #[serde(rename = "Bq/kg")]
    BqPerKg,
    /// Volumetric activity. Needs a density.
    #[serde(rename = "Bq/cm3")]
    BqPerCm3,
    /// Volumetric activity. Needs a density.
    #[serde(rename = "Bq/m3")]
    BqPerM3,
    /// Total activity in curies. Needs the material's mass.
    #[serde(rename = "Ci")]
    Ci,
    /// Volumetric activity, the unit the US tables use. Needs a density.
    #[serde(rename = "Ci/m3")]
    CiPerM3,
}

impl ActivityUnit {
    /// Every unit, in the order [`ActivityUnit::as_str`] lists them.
    pub const ALL: [ActivityUnit; 7] = [
        ActivityUnit::Bq,
        ActivityUnit::BqPerG,
        ActivityUnit::BqPerKg,
        ActivityUnit::BqPerCm3,
        ActivityUnit::BqPerM3,
        ActivityUnit::Ci,
        ActivityUnit::CiPerM3,
    ];

    /// The unit as written, such as `"Bq/g"`.
    pub fn as_str(self) -> &'static str {
        match self {
            ActivityUnit::Bq => "Bq",
            ActivityUnit::BqPerG => "Bq/g",
            ActivityUnit::BqPerKg => "Bq/kg",
            ActivityUnit::BqPerCm3 => "Bq/cm3",
            ActivityUnit::BqPerM3 => "Bq/m3",
            ActivityUnit::Ci => "Ci",
            ActivityUnit::CiPerM3 => "Ci/m3",
        }
    }
}

impl std::str::FromStr for ActivityUnit {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self> {
        ActivityUnit::ALL
            .into_iter()
            .find(|u| u.as_str() == text)
            .ok_or_else(|| {
                let known: Vec<&str> = ActivityUnit::ALL.iter().map(|u| u.as_str()).collect();
                Error::Invalid(format!(
                    "unknown activity units '{text}', expected one of {}",
                    known.join(", ")
                ))
            })
    }
}

impl fmt::Display for ActivityUnit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// What the inventory's numbers are, which decides what can be derived.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Basis {
    /// Absolute atom counts.
    AtomCounts,
    /// Atoms per barn-cm, which fixes the mass density.
    AtomDensities,
    /// Absolute grams per nuclide.
    Masses,
    /// Relative mass fractions.
    MassFractions,
    /// Specific activities in Bq/g, as an assay reports them.
    SpecificActivities,
}

/// A nuclide inventory to assess for clearance.
///
/// Specific activity in Bq/g is scale invariant, so atom counts, atom
/// densities and atom fractions all give the same answer for the Bq/g limit
/// sets and no density is needed. Volumetric limit sets (Ci/m3) do need one,
/// supplied with [`Material::with_density`] or derived from atom densities.
///
/// The Bq/g denominator is the mass of **everything** in the inventory, so
/// stable isotopes must be included. Passing only the radioactive nuclides of
/// an activated steel gives a mass thousands of times too small and a specific
/// activity thousands of times too high. [`Material::looks_truncated`] flags
/// the common form of that mistake.
///
/// ```
/// use radiological_material_clearance_finder::Material;
///
/// let steel = Material::from_atom_counts([("Fe56", 8.4e22), ("Co60", 2.1e9)])?;
/// assert!(steel.specific_activity()? > 0.0);
/// # Ok::<(), radiological_material_clearance_finder::Error>(())
/// ```
#[derive(Debug, Clone)]
pub struct Material {
    /// A label carried through to results, for reporting.
    pub name: String,
    basis: Basis,
    amounts: BTreeMap<String, f64>,
    density: Option<f64>,
    volume: Option<f64>,
    decay_data: Arc<DecayData>,
}

/// Reject a non-positive or non-finite density or volume.
///
/// A density of zero would make every volumetric activity zero, and so every
/// Ci/m3 index zero, which reports as clearable. That is the most dangerous way
/// for this to fail, so it is rejected up front.
fn positive(what: &str, value: f64) -> Result<f64> {
    if !(value > 0.0 && value.is_finite()) {
        return Err(Error::Invalid(format!(
            "{what} must be greater than zero, got {value:?}. A non-positive {what} \
             would make every volumetric activity zero and so report the material as \
             clearable."
        )));
    }
    Ok(value)
}

impl Material {
    fn build<I, K>(basis: Basis, entries: I) -> Result<Material>
    where
        I: IntoIterator<Item = (K, f64)>,
        K: AsRef<str>,
    {
        let mut amounts: BTreeMap<String, f64> = BTreeMap::new();
        for (raw, value) in entries {
            let name = nuclide::normalise(raw.as_ref())?;
            // NaN fails every comparison, so it must be tested for directly or
            // it propagates through the whole index as a NaN verdict.
            if value.is_nan() {
                return Err(Error::Invalid(format!("amount for {name} is NaN")));
            }
            if value < 0.0 {
                return Err(Error::Invalid(format!(
                    "negative amount {value:?} for {name}"
                )));
            }
            if value.is_infinite() {
                return Err(Error::Invalid(format!("amount for {name} is infinite")));
            }
            *amounts.entry(name).or_insert(0.0) += value;
        }
        Ok(Material {
            name: String::new(),
            basis,
            amounts,
            density: None,
            volume: None,
            decay_data: DecayData::shared_default(),
        })
    }

    /// Build from absolute atom counts, what an activation or depletion
    /// calculation produces.
    pub fn from_atom_counts<I, K>(atoms: I) -> Result<Material>
    where
        I: IntoIterator<Item = (K, f64)>,
        K: AsRef<str>,
    {
        Material::build(Basis::AtomCounts, atoms)
    }

    /// Build from atom densities in atoms per barn-cm, OpenMC's and yamc's unit.
    ///
    /// The mass density follows from the atom densities and the atomic masses,
    /// so [`Material::with_density`] is refused: it could only contradict the
    /// inventory. Add a volume for total activity.
    pub fn from_atom_densities<I, K>(densities: I) -> Result<Material>
    where
        I: IntoIterator<Item = (K, f64)>,
        K: AsRef<str>,
    {
        Material::build(Basis::AtomDensities, densities)
    }

    /// Build from a mass in grams per nuclide.
    pub fn from_masses<I, K>(masses: I) -> Result<Material>
    where
        I: IntoIterator<Item = (K, f64)>,
        K: AsRef<str>,
    {
        Material::build(Basis::Masses, masses)
    }

    /// Build from mass fractions, which need not sum to one.
    pub fn from_mass_fractions<I, K>(fractions: I) -> Result<Material>
    where
        I: IntoIterator<Item = (K, f64)>,
        K: AsRef<str>,
    {
        Material::build(Basis::MassFractions, fractions)
    }

    /// Build from specific activities in Bq/g, as an assay reports them.
    ///
    /// No half-life or atomic mass data is used, since the specific activity is
    /// the direct input to the index. Volumetric limit sets then need an
    /// explicit density, because activity alone does not imply a mass.
    pub fn from_specific_activities<I, K>(activities: I) -> Result<Material>
    where
        I: IntoIterator<Item = (K, f64)>,
        K: AsRef<str>,
    {
        Material::build(Basis::SpecificActivities, activities)
    }

    /// Set the mass density in g/cm3, needed for volumetric limit sets.
    pub fn with_density(mut self, density: f64) -> Result<Material> {
        if self.basis == Basis::AtomDensities {
            return Err(Error::Invalid(
                "a material built from atom densities derives its mass density from \
                 them and the atomic masses, so setting a density could only contradict \
                 the inventory. Build it from specific activities if the density is the \
                 value you trust."
                    .into(),
            ));
        }
        self.density = Some(positive("density", density)?);
        Ok(self)
    }

    /// Set the volume in cm3, needed for total activity in Bq or Ci.
    pub fn with_volume(mut self, volume: f64) -> Result<Material> {
        self.set_volume(Some(volume))?;
        Ok(self)
    }

    /// Set the label carried through to results.
    pub fn with_name(mut self, name: impl Into<String>) -> Material {
        self.name = name.into();
        self
    }

    /// Use these half-life and mass tables instead of the vendored ones.
    pub fn with_decay_data(mut self, decay_data: Arc<DecayData>) -> Material {
        self.decay_data = decay_data;
        self
    }

    /// The decay data this material's activities are computed from.
    pub fn decay_data(&self) -> &Arc<DecayData> {
        &self.decay_data
    }

    /// The volume in cm3, if known.
    pub fn volume(&self) -> Option<f64> {
        self.volume
    }

    /// Set or clear the volume in cm3.
    pub fn set_volume(&mut self, volume: Option<f64>) -> Result<()> {
        self.volume = volume.map(|v| positive("volume", v)).transpose()?;
        Ok(())
    }

    /// Whether the inventory holds no stable nuclides, which usually means it
    /// was filtered down to its radioactive part.
    ///
    /// An activated material is overwhelmingly stable by mass. An inventory
    /// that is entirely radioactive makes every specific activity too high by
    /// the ratio of true mass to retained mass. A pure source is legitimate, so
    /// a single nuclide never counts, and this is advice rather than an error.
    pub fn looks_truncated(&self) -> bool {
        if self.amounts.len() < 2 || self.basis == Basis::SpecificActivities {
            return false;
        }
        self.amounts
            .keys()
            .all(|n| matches!(self.decay_data.half_life_canonical(n), Ok(Some(_))))
    }

    /// The nuclides present, in canonical form, sorted.
    pub fn nuclides(&self) -> impl ExactSizeIterator<Item = &str> {
        self.amounts.keys().map(String::as_str)
    }

    /// Atom amounts, keyed by canonical name.
    ///
    /// Counts or densities as supplied. Masses are converted to counts, and
    /// mass fractions to relative amounts on the same scale. A material built
    /// from specific activities has none.
    pub fn atoms(&self) -> Result<BTreeMap<String, f64>> {
        match self.basis {
            Basis::AtomCounts | Basis::AtomDensities => Ok(self.amounts.clone()),
            Basis::Masses | Basis::MassFractions => self
                .amounts
                .iter()
                .map(|(n, &grams)| {
                    Ok((
                        n.clone(),
                        grams * AVOGADRO / self.decay_data.atomic_mass_canonical(n)?,
                    ))
                })
                .collect(),
            Basis::SpecificActivities => Err(Error::InsufficientData(
                "this material was built from specific activities, which do not \
                 determine atom counts"
                    .into(),
            )),
        }
    }

    /// Mass in grams of the amounts as supplied, whatever their scale.
    fn relative_mass(&self) -> Result<f64> {
        match self.basis {
            Basis::Masses | Basis::MassFractions => Ok(crate::total(self.amounts.values())),
            Basis::AtomCounts | Basis::AtomDensities => {
                let mut total = 0.0;
                for (n, &amount) in &self.amounts {
                    total += amount * self.decay_data.atomic_mass_canonical(n)? / AVOGADRO;
                }
                Ok(total)
            }
            Basis::SpecificActivities => Err(Error::InsufficientData(
                "a material built from specific activities has no intrinsic mass".into(),
            )),
        }
    }

    /// Whether the amounts are absolute, so their mass is the material's mass.
    fn absolute(&self) -> bool {
        matches!(self.basis, Basis::AtomCounts | Basis::Masses)
    }

    /// Total mass in grams.
    ///
    /// An error if the material carries only relative amounts, such as atom
    /// densities or mass fractions, with no volume to scale them by, or if its
    /// absolute amounts and its density times volume disagree by more than one
    /// percent.
    pub fn mass(&self) -> Result<f64> {
        let from_atoms = if self.absolute() {
            Some(self.relative_mass()?)
        } else {
            None
        };
        let from_volume = match (self.density()?, self.volume) {
            (Some(d), Some(v)) => Some(d * v),
            _ => None,
        };
        match (from_atoms, from_volume) {
            (Some(a), Some(v)) => {
                // Over-determined. Disagreement means the density and the
                // inventory describe different materials, and silently
                // preferring either would put the total activity out by
                // whatever the ratio happens to be.
                if (a - v).abs() > 0.01 * a.max(v) {
                    return Err(Error::Invalid(format!(
                        "this material is over-determined and inconsistent: the atom \
                         amounts weigh {} g, but density times volume gives {} g. Drop \
                         whichever of the two is not meant to describe this material.",
                        g(a, 6),
                        g(v, 6)
                    )));
                }
                Ok(a)
            }
            (Some(a), None) => Ok(a),
            (None, Some(v)) => Ok(v),
            (None, None) => Err(Error::InsufficientData(
                "total mass is unknown. This material holds relative amounts, so give \
                 it a volume, or build it from atom counts or masses."
                    .into(),
            )),
        }
    }

    /// The mass density in g/cm3 if it is known or derivable, else `None`.
    fn density(&self) -> Result<Option<f64>> {
        if self.density.is_some() {
            return Ok(self.density);
        }
        match self.basis {
            Basis::AtomDensities => Ok(Some(self.relative_mass()? * ATOMS_PER_BARN_CM)),
            _ if self.absolute() => match self.volume {
                Some(v) => Ok(Some(self.relative_mass()? / v)),
                None => Ok(None),
            },
            _ => Ok(None),
        }
    }

    /// Mass density in g/cm3.
    ///
    /// An error, naming what to supply, if none was given and none can be
    /// derived.
    pub fn mass_density(&self) -> Result<f64> {
        self.density()?.ok_or_else(|| {
            Error::InsufficientData(
                "mass density is unknown, and volumetric limits (Ci/m3) need it. Give \
                 the material a density in g/cm3, or build it from atom densities, which \
                 derives it."
                    .into(),
            )
        })
    }

    /// Specific activity in Bq/g per nuclide.
    ///
    /// Scale invariant, so this works from atom counts, atom densities or atom
    /// fractions alike.
    pub fn specific_activities(&self) -> Result<BTreeMap<String, f64>> {
        if self.basis == Basis::SpecificActivities {
            return Ok(self.amounts.clone());
        }
        let mass = self.relative_mass()?;
        if mass <= 0.0 {
            return Err(Error::InsufficientData("material has zero mass".into()));
        }
        let mut out = BTreeMap::new();
        for (n, &amount) in &self.amounts {
            let atoms = match self.basis {
                Basis::Masses | Basis::MassFractions => {
                    amount * AVOGADRO / self.decay_data.atomic_mass_canonical(n)?
                }
                _ => amount,
            };
            out.insert(
                n.clone(),
                self.decay_data.decay_constant_canonical(n)? * atoms / mass,
            );
        }
        Ok(out)
    }

    /// Total specific activity in Bq/g.
    pub fn specific_activity(&self) -> Result<f64> {
        Ok(crate::total(self.specific_activities()?.values()))
    }

    /// Activity per nuclide in the requested units.
    pub fn activities(&self, units: ActivityUnit) -> Result<BTreeMap<String, f64>> {
        let factor = match units {
            ActivityUnit::BqPerG => 1.0,
            ActivityUnit::BqPerKg => 1000.0,
            ActivityUnit::BqPerCm3 => self.mass_density()?,
            ActivityUnit::BqPerM3 => self.mass_density()? * 1e6,
            ActivityUnit::CiPerM3 => self.mass_density()? * 1e6 / BECQUEREL_PER_CURIE,
            ActivityUnit::Bq => self.mass()?,
            ActivityUnit::Ci => self.mass()? / BECQUEREL_PER_CURIE,
        };
        let mut per_gram = self.specific_activities()?;
        if factor != 1.0 {
            per_gram.values_mut().for_each(|v| *v *= factor);
        }
        Ok(per_gram)
    }

    /// Total activity in the requested units.
    pub fn activity(&self, units: ActivityUnit) -> Result<f64> {
        Ok(crate::total(self.activities(units)?.values()))
    }
}
