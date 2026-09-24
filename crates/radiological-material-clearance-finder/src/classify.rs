//! Waste classification, where the answer is a category rather than an index.
//!
//! The US and UK both sort waste into named classes rather than reporting a
//! single ratio. The NRC's classes come out of the same sum-of-fractions
//! machinery as a clearance index, applied to several tables and combined by
//! rules in 10 CFR 61.55. The UK's categories instead compare gross alpha and
//! beta or gamma activity against concentration thresholds, so they need no
//! per-nuclide table at all and are computed directly.

use std::fmt;

use crate::error::Result;
use crate::fmt::g;
use crate::index::{clearance_index, ClearanceOptions};
use crate::limits::get_limit_set;
use crate::material::{ActivityUnit, Material};

// The 2007 policy states these per tonne. One GBq/tonne is 1e9 Bq per 1e6 g,
// which is 1e3 Bq/g, and one MBq/tonne is 1 Bq/g.
/// LLW upper bound on alpha activity, from 4 GBq/te.
pub const UK_LLW_ALPHA_BQ_PER_G: f64 = 4.0e3;
/// LLW upper bound on beta and gamma activity, from 12 GBq/te.
pub const UK_LLW_BETA_GAMMA_BQ_PER_G: f64 = 12.0e3;
/// High volume VLLW upper bound on total activity, from 4 MBq/te.
pub const UK_HIGH_VOLUME_VLLW_BQ_PER_G: f64 = 4.0;
/// High volume VLLW upper bound on tritium and carbon-14 together, from
/// 40 MBq/te. The 2007 policy gives them a shared allowance in both the low
/// volume and the high volume categories, so they are summed rather than
/// tritium being singled out.
pub const UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G: f64 = 40.0;

/// A near-surface disposal class under 10 CFR 61.55.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum NrcWasteClass {
    /// Class A.
    A,
    /// Class B.
    B,
    /// Class C.
    C,
    /// Greater than Class C, not generally acceptable for near-surface disposal.
    Gtcc,
}

impl NrcWasteClass {
    /// `"Class A"`, `"Class B"`, `"Class C"` or `"GTCC"`.
    pub fn as_str(self) -> &'static str {
        match self {
            NrcWasteClass::A => "Class A",
            NrcWasteClass::B => "Class B",
            NrcWasteClass::C => "Class C",
            NrcWasteClass::Gtcc => "GTCC",
        }
    }
}

impl fmt::Display for NrcWasteClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Classify a material for near-surface disposal under 10 CFR 61.55.
///
/// Applies the sum of fractions rule to Table 1 and Table 2 and combines the
/// results with the rules in paragraphs 61.55(a)(3) to (a)(7). The material
/// needs a density, since the NRC limits are volumetric. `metal` selects the
/// activated metal rows.
pub fn nrc_waste_class(material: &Material, metal: bool) -> Result<NrcWasteClass> {
    let options = ClearanceOptions {
        metal,
        ..ClearanceOptions::default()
    };
    let index = |name: &str| -> Result<f64> {
        Ok(clearance_index(material, &*get_limit_set(name)?, options)?.index)
    };
    let table1 = index("NRC_long")?;
    let table2 = [
        index("NRC_short_A")?,
        index("NRC_short_B")?,
        index("NRC_short_C")?,
    ];

    let by_table2 = |[a, b, c]: [f64; 3]| {
        if a < 1.0 {
            NrcWasteClass::A
        } else if b < 1.0 {
            NrcWasteClass::B
        } else if c < 1.0 {
            NrcWasteClass::C
        } else {
            NrcWasteClass::Gtcc
        }
    };
    let table1_present = table1 > 0.0;
    let table2_present = table2.iter().any(|&v| v > 0.0);

    Ok(match (table1_present, table2_present) {
        // 61.55(a)(5)
        (true, true) => {
            if table1 < 0.1 {
                by_table2(table2)
            } else if table1 < 1.0 && table2[2] < 1.0 {
                NrcWasteClass::C
            } else {
                NrcWasteClass::Gtcc
            }
        }
        // 61.55(a)(3)
        (true, false) => {
            if table1 < 0.1 {
                NrcWasteClass::A
            } else if table1 < 1.0 {
                NrcWasteClass::C
            } else {
                NrcWasteClass::Gtcc
            }
        }
        // 61.55(a)(4)
        (false, true) => by_table2(table2),
        // 61.55(a)(6)
        (false, false) => NrcWasteClass::A,
    })
}

/// Activity from alpha emission, weighted by each nuclide's alpha branch.
///
/// Bi-212 branches 35.94 percent alpha, so it contributes that share of its
/// activity here and the rest to [`beta_gamma_activity`].
pub fn alpha_activity(material: &Material, units: ActivityUnit) -> Result<f64> {
    let data = material.decay_data();
    Ok(crate::total(
        material
            .activities(units)?
            .iter()
            .map(|(n, v)| v * data.alpha_fraction_canonical(n)),
    ))
}

/// Activity from everything that is not alpha emission.
pub fn beta_gamma_activity(material: &Material, units: ActivityUnit) -> Result<f64> {
    let data = material.decay_data();
    Ok(crate::total(material.activities(units)?.iter().map(
        |(n, v)| v * (1.0 - data.alpha_fraction_canonical(n)),
    )))
}

/// A UK solid radioactive waste category.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize)]
pub enum UkCategory {
    /// Very low level waste (high volume).
    #[serde(rename = "VLLW")]
    Vllw,
    /// Low level waste.
    #[serde(rename = "LLW")]
    Llw,
    /// Intermediate level waste.
    #[serde(rename = "ILW")]
    Ilw,
}

impl UkCategory {
    /// `"VLLW"`, `"LLW"` or `"ILW"`.
    pub fn as_str(self) -> &'static str {
        match self {
            UkCategory::Vllw => "VLLW",
            UkCategory::Llw => "LLW",
            UkCategory::Ilw => "ILW",
        }
    }
}

impl fmt::Display for UkCategory {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// The UK waste category of a material, with the numbers behind it.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct UkWasteCategory {
    /// The category.
    pub category: UkCategory,
    /// Alpha activity in Bq/g.
    pub alpha: f64,
    /// Beta and gamma activity in Bq/g.
    pub beta_gamma: f64,
    /// Tritium plus carbon-14 activity in Bq/g, which share their own VLLW
    /// allowance.
    pub tritium_and_c14: f64,
    /// Total activity in Bq/g.
    pub total: f64,
    /// Why this category and not the one below it.
    pub reason: String,
}

impl fmt::Display for UkWasteCategory {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.category, self.reason)
    }
}

/// Classify a material against the UK solid low level waste categories.
///
/// Thresholds are from *Policy for the Long Term Management of Solid Low Level
/// Radioactive Waste in the United Kingdom* (Defra and the Devolved
/// Administrations, 2007), which defines LLW as not exceeding 4 GBq/te alpha or
/// 12 GBq/te beta and gamma, and high volume VLLW as not exceeding 4 MBq/te
/// total activity, with tritium allowed to 40 MBq/te.
///
/// Two limits of this classification are worth stating rather than hiding. The
/// low volume VLLW category ("dustbin disposal") is defined per 0.1 cubic metre
/// and per item, so it is a property of a consignment rather than of a material
/// and is not decided here. The boundary above which waste becomes high level
/// rather than intermediate is a thermal one, about 2 kW/m3, and needs decay
/// heat rather than activity, so this returns ILW for everything above LLW.
pub fn uk_waste_category(material: &Material) -> Result<UkWasteCategory> {
    let per_nuclide = material.activities(ActivityUnit::BqPerG)?;
    let data = material.decay_data();
    let alpha = crate::total(
        per_nuclide
            .iter()
            .map(|(n, v)| v * data.alpha_fraction_canonical(n)),
    );
    let beta_gamma = crate::total(
        per_nuclide
            .iter()
            .map(|(n, v)| v * (1.0 - data.alpha_fraction_canonical(n))),
    );
    let tritium_and_c14 = per_nuclide.get("H3").copied().unwrap_or(0.0)
        + per_nuclide.get("C14").copied().unwrap_or(0.0);
    let total = crate::total(per_nuclide.values());
    let remainder = total - tritium_and_c14;

    let (category, reason) = if alpha > UK_LLW_ALPHA_BQ_PER_G {
        (
            UkCategory::Ilw,
            format!(
                "alpha activity {} Bq/g exceeds the LLW limit of {} Bq/g (4 GBq/te)",
                g(alpha, 3),
                g(UK_LLW_ALPHA_BQ_PER_G, 6)
            ),
        )
    } else if beta_gamma > UK_LLW_BETA_GAMMA_BQ_PER_G {
        (
            UkCategory::Ilw,
            format!(
                "beta and gamma activity {} Bq/g exceeds the LLW limit of {} Bq/g (12 GBq/te)",
                g(beta_gamma, 3),
                g(UK_LLW_BETA_GAMMA_BQ_PER_G, 6)
            ),
        )
    } else if remainder > UK_HIGH_VOLUME_VLLW_BQ_PER_G {
        (
            UkCategory::Llw,
            format!(
                "total activity {} Bq/g excluding tritium and carbon-14 exceeds the high \
                 volume VLLW limit of {} Bq/g (4 MBq/te), and it is within both LLW limits",
                g(remainder, 3),
                g(UK_HIGH_VOLUME_VLLW_BQ_PER_G, 6)
            ),
        )
    } else if tritium_and_c14 > UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G {
        (
            UkCategory::Llw,
            format!(
                "tritium and carbon-14 activity {} Bq/g exceeds their shared high volume \
                 VLLW allowance of {} Bq/g (40 MBq/te)",
                g(tritium_and_c14, 3),
                g(UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G, 6)
            ),
        )
    } else {
        (
            UkCategory::Vllw,
            format!(
                "total activity {} Bq/g excluding tritium and carbon-14 is within the high \
                 volume VLLW limit of {} Bq/g (4 MBq/te)",
                g(remainder, 3),
                g(UK_HIGH_VOLUME_VLLW_BQ_PER_G, 6)
            ),
        )
    };
    Ok(UkWasteCategory {
        category,
        alpha,
        beta_gamma,
        tritium_and_c14,
        total,
        reason,
    })
}
