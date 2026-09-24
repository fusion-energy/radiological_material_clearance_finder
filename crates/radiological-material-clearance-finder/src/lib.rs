//! Clearance indexes for a radiological material.
//!
//! Give it a nuclide inventory and it tells you whether the material meets the
//! clearance, exemption or disposal limits of the UK, German, US, EU and IAEA
//! regulations, and by how much.
//!
//! ```
//! use radiological_material_clearance_finder::{
//!     clearance_index, get_limit_set, ClearanceOptions, Material,
//! };
//!
//! let steel = Material::from_atom_counts([("Fe56", 8.4e22), ("Co60", 2.1e9), ("Cs137", 3.1e8)])?;
//! let set = get_limit_set("UK_EPR16_out_of_scope")?;
//! let result = clearance_index(&steel, &set, ClearanceOptions::default())?;
//! assert!(!result.clearable());
//! assert_eq!(result.dominant(1)[0].0, "Co60");
//! # Ok::<(), radiological_material_clearance_finder::Error>(())
//! ```
//!
//! Every regulation here uses the same arithmetic, a sum of activity-to-limit
//! ratios that must stay below one. What differs is the tables, and those are
//! data compiled into the crate: see [`limit_sets`] for what is available and
//! [`get_limit_set`] for the provenance of any one of them.
//!
//! A transport or transmutation code with its own material type hands over its
//! atom densities in atoms per barn-cm, and its volume if total activity is
//! wanted:
//!
//! ```
//! # use std::collections::HashMap;
//! use radiological_material_clearance_finder::Material;
//!
//! let atoms_per_barn_cm: HashMap<String, f64> =
//!     [("Fe56".to_string(), 0.0849), ("Co60".to_string(), 1e-10)].into();
//! let inventory = Material::from_atom_densities(atoms_per_barn_cm)?.with_volume(1000.0)?;
//! assert!((inventory.mass_density()? - 7.89).abs() < 0.01);
//! # Ok::<(), radiological_material_clearance_finder::Error>(())
//! ```

pub mod classify;
pub mod cooling;
pub mod decay;
pub mod error;
mod fmt;
pub mod index;
pub mod limits;
pub mod material;
pub mod nuclide;

pub use classify::{
    alpha_activity, beta_gamma_activity, nrc_waste_class, uk_waste_category, NrcWasteClass,
    UkCategory, UkWasteCategory,
};
pub use cooling::{index_series, time_to_clear};
pub use decay::{DecayData, AVOGADRO, BECQUEREL_PER_CURIE};
pub use error::{Error, Result};
pub use index::{
    clearable_routes, clearance_index, clearance_indices, ClearanceOptions, ClearanceResult,
    EQUILIBRIUM_TOLERANCE,
};
pub use limits::{
    get_limit_set, limit_sets, register_limit_set, DynamicRule, LimitSet, LimitSetData,
};
pub use material::{ActivityUnit, Material};

/// Sum of floats starting from +0.0.
///
/// `Iterator::sum` on floats starts from -0.0, so an empty sum, such as the
/// index of a material with nothing limited, would print as "-0".
pub(crate) fn total<I, T>(values: I) -> f64
where
    I: IntoIterator<Item = T>,
    T: std::borrow::Borrow<f64>,
{
    values.into_iter().fold(0.0, |acc, v| acc + *v.borrow())
}
