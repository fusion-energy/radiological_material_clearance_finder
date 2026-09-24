//! Python bindings for the `radiological-material-clearance-finder` crate.
//!
//! Everything here is conversion: Python mappings to Rust maps, Rust errors to
//! the Python exception classes, and Rust results to dicts and tuples. The
//! calculation lives in the core crate, which yamc and yani link directly.

use std::collections::BTreeMap;
use std::ffi::CString;
use std::path::PathBuf;
use std::sync::Arc;

use pyo3::create_exception;
use pyo3::exceptions::{PyKeyError, PyOSError, PyTypeError, PyUserWarning, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyFrozenSet, PyTuple};

use radiological_material_clearance_finder as rmcf;
use rmcf::{ActivityUnit, ClearanceOptions, DecayData, LimitSet, LimitSetData, Material};

create_exception!(
    radiological_material_clearance_finder._core,
    NuclideNameError,
    PyValueError,
    "Raised when a string cannot be read as a single nuclide."
);
create_exception!(
    radiological_material_clearance_finder._core,
    UnknownNuclideError,
    PyKeyError,
    "Raised when a nuclide is absent from the decay data tables.\n\n\
     The canonical nuclide name is available as the ``name`` attribute."
);
create_exception!(
    radiological_material_clearance_finder._core,
    InsufficientDataError,
    PyValueError,
    "Raised when a quantity cannot be derived from what the material was given."
);
create_exception!(
    radiological_material_clearance_finder._core,
    IngrowthError,
    PyValueError,
    "Raised when the index climbs back above the threshold after clearing.\n\n\
     An index that rises with cooling time means a daughter is growing in faster \
     than its parent decays, so a material can meet the limits at one time and fail \
     them later. Reporting only the first crossing would be actively misleading, so \
     this is raised instead."
);

/// Map a core error onto the Python exception a caller would catch.
fn to_py(error: rmcf::Error) -> PyErr {
    match error {
        rmcf::Error::NuclideName(m) => NuclideNameError::new_err(m),
        rmcf::Error::UnknownNuclide { name, message } => Python::attach(|py| {
            let err = UnknownNuclideError::new_err(message);
            // Best effort: the message already names the nuclide.
            let _ = err.value(py).setattr("name", name);
            err
        }),
        rmcf::Error::InsufficientData(m) => InsufficientDataError::new_err(m),
        rmcf::Error::UnknownLimitSet(m) => PyKeyError::new_err(m),
        rmcf::Error::Ingrowth(m) => IngrowthError::new_err(m),
        rmcf::Error::Io(m) => PyOSError::new_err(m),
        rmcf::Error::Invalid(m) => PyValueError::new_err(m),
    }
}

trait OrPy<T> {
    fn py_err(self) -> PyResult<T>;
}

impl<T> OrPy<T> for rmcf::Result<T> {
    fn py_err(self) -> PyResult<T> {
        self.map_err(to_py)
    }
}

fn warn(py: Python<'_>, message: &str) -> PyResult<()> {
    let message = CString::new(message).expect("warning messages hold no NUL");
    PyErr::warn(py, &py.get_type::<PyUserWarning>(), &message, 1)
}

fn units(text: &str) -> PyResult<ActivityUnit> {
    text.parse().py_err()
}

fn pairs_to_dict<'py, V>(py: Python<'py>, pairs: &[(String, V)]) -> PyResult<Bound<'py, PyDict>>
where
    V: Clone + IntoPyObject<'py>,
{
    let dict = PyDict::new(py);
    for (k, v) in pairs {
        dict.set_item(k, v.clone())?;
    }
    Ok(dict)
}

fn strings_to_tuple<'py>(py: Python<'py>, items: &[String]) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, items)
}

// ----------------------------------------------------------------------
// nuclide names
// ----------------------------------------------------------------------

/// A nuclide name argument, which has to be a string.
///
/// Anything else is a NuclideNameError rather than a TypeError, so callers
/// validating untrusted input have one exception to catch.
fn name_of(obj: &Bound<'_, PyAny>) -> PyResult<String> {
    match obj.extract::<String>() {
        Ok(name) => Ok(name),
        Err(_) => Err(NuclideNameError::new_err(format!(
            "nuclide name must be a string, got {}",
            obj.get_type().name()?
        ))),
    }
}

/// Split a nuclide name into element symbol, mass number and metastable level.
///
/// Args:
///     name: A nuclide name in any accepted spelling, such as ``Co60``,
///         ``Co-60``, ``Ag-108m`` or ``Ag108_m1``.
///
/// Returns:
///     ``(symbol, mass_number, metastable_level)``, where the level is 0 for a
///     ground state and 1 or 2 for a metastable state. The symbol is returned
///     in its canonical capitalisation regardless of the input's.
///
/// Raises:
///     NuclideNameError: If the name is not a single nuclide. A bare element
///         such as ``Fe`` raises, because a clearance limit applies to a
///         nuclide and silently guessing a mass number would be wrong.
#[pyfunction]
fn parse(name: &Bound<'_, PyAny>) -> PyResult<(&'static str, u32, u8)> {
    rmcf::nuclide::parse(&name_of(name)?).py_err()
}

/// Return a nuclide name in canonical form.
///
/// ``Co-60`` and ``co 60`` both become ``Co60``; ``Ag-108m`` becomes
/// ``Ag108_m1`` and ``Hf-178n`` becomes ``Hf178_m2``. A trailing ``+`` is
/// stripped, since it marks a secular equilibrium value rather than a different
/// nuclide, but ``sec`` is not accepted here because a ``sec`` value is a
/// separate table entry. Use `parse_regulatory` when reading a regulatory
/// table.
///
/// Args:
///     name: A nuclide name in any accepted spelling.
///
/// Returns:
///     The canonical name, such as ``Co60`` or ``Ag108_m1``.
///
/// Raises:
///     NuclideNameError: If the name is not a single nuclide.
#[pyfunction]
fn normalise(name: &Bound<'_, PyAny>) -> PyResult<String> {
    rmcf::nuclide::normalise(&name_of(name)?).py_err()
}

/// Read a nuclide name as spelled in a regulatory table.
///
/// Regulatory tables mark parent nuclides whose limit already accounts for
/// daughters in secular equilibrium, with ``+`` for the daughters tabulated
/// alongside and ``sec`` for the whole decay chain. The marker is not part of
/// the nuclide's identity but it does select a different limit value, so it is
/// returned separately rather than discarded.
///
/// Args:
///     label: A table cell such as ``Sr-90+``, ``U-238sec`` or ``Th-232 sec``.
///
/// Returns:
///     ``(canonical_name, marker)`` where marker is ``"+"``, ``"sec"`` or
///     ``None``.
///
/// Raises:
///     NuclideNameError: If the label is not a single nuclide.
#[pyfunction]
fn parse_regulatory(label: &Bound<'_, PyAny>) -> PyResult<(String, Option<&'static str>)> {
    let (name, marker) = rmcf::nuclide::parse_regulatory(&name_of(label)?).py_err()?;
    Ok((name, marker.map(|m| m.as_str())))
}

/// Return the element symbol of a nuclide, for example ``Co`` for ``Co60``.
#[pyfunction]
fn element(name: &Bound<'_, PyAny>) -> PyResult<&'static str> {
    rmcf::nuclide::element(&name_of(name)?).py_err()
}

/// Return the mass number of a nuclide, for example 60 for ``Co60``.
#[pyfunction]
fn mass_number(name: &Bound<'_, PyAny>) -> PyResult<u32> {
    rmcf::nuclide::mass_number(&name_of(name)?).py_err()
}

/// Return the metastable level of a nuclide, 0 for a ground state.
#[pyfunction]
fn metastable_state(name: &Bound<'_, PyAny>) -> PyResult<u8> {
    rmcf::nuclide::metastable_state(&name_of(name)?).py_err()
}

/// Return the proton number of a nuclide, for example 27 for ``Co60``.
#[pyfunction]
fn atomic_number(name: &Bound<'_, PyAny>) -> PyResult<u32> {
    rmcf::nuclide::atomic_number(&name_of(name)?).py_err()
}

/// Whether a value can be read as a single nuclide.
#[pyfunction]
fn is_valid(name: &Bound<'_, PyAny>) -> bool {
    name.extract::<String>()
        .is_ok_and(|n| rmcf::nuclide::is_valid(&n))
}

// ----------------------------------------------------------------------
// decay data
// ----------------------------------------------------------------------

/// Half-life and atomic mass lookups for a set of nuclides.
///
/// Every lookup accepts any spelling `normalise` does. A nuclide in neither
/// table raises `UnknownNuclideError` rather than being assumed stable.
///
/// Args:
///     half_lives: Half-lives in seconds keyed by nuclide name. A nuclide
///         present in ``atomic_masses`` but absent here is stable. Defaults to
///         the vendored ENDF/B-VIII.0 values.
///     atomic_masses: Relative atomic masses in u keyed by nuclide name.
///         Defaults to the vendored AME2020 values.
#[pyclass(
    name = "DecayData",
    module = "radiological_material_clearance_finder._core",
    frozen
)]
struct PyDecayData {
    inner: Arc<DecayData>,
}

#[pymethods]
impl PyDecayData {
    #[new]
    #[pyo3(signature = (half_lives=None, atomic_masses=None))]
    fn new(
        half_lives: Option<BTreeMap<String, f64>>,
        atomic_masses: Option<BTreeMap<String, f64>>,
    ) -> PyResult<Self> {
        let base = DecayData::shared_default();
        if half_lives.is_none() && atomic_masses.is_none() {
            return Ok(PyDecayData { inner: base });
        }
        let data = match (half_lives, atomic_masses) {
            (Some(h), Some(m)) => DecayData::new(h, m),
            (Some(h), None) => DecayData::new(h, base.atomic_masses().clone()),
            (None, Some(m)) => DecayData::new(base.half_lives().clone(), m),
            (None, None) => unreachable!(),
        };
        Ok(PyDecayData {
            inner: Arc::new(data.py_err()?),
        })
    }

    /// Take half-lives from an OpenMC depletion chain file.
    ///
    /// Parsed without OpenMC. Nuclides absent from the chain fall back to the
    /// vendored ENDF values, matching how OpenMC resolves a chain against its
    /// own data. A chain entry with no ``half_life`` marks the nuclide stable.
    ///
    /// Args:
    ///     path: Path to a depletion chain XML file.
    ///
    /// Returns:
    ///     A `DecayData` using the chain's half-lives.
    ///
    /// Raises:
    ///     ValueError: If the file is not a readable chain, or a half-life is
    ///         not a positive finite number.
    #[staticmethod]
    fn from_chain_xml(path: PathBuf) -> PyResult<Self> {
        Ok(PyDecayData {
            inner: Arc::new(DecayData::from_chain_xml(path).py_err()?),
        })
    }

    /// Return a copy with some half-lives replaced.
    ///
    /// Args:
    ///     half_lives: Half-lives in seconds, keyed by any accepted nuclide
    ///         spelling. A value of ``None`` marks a nuclide stable.
    fn with_overrides(&self, half_lives: BTreeMap<String, Option<f64>>) -> PyResult<Self> {
        Ok(PyDecayData {
            inner: Arc::new(self.inner.with_overrides(half_lives).py_err()?),
        })
    }

    /// Whether the tables cover this nuclide at all.
    fn knows(&self, name: &str) -> PyResult<bool> {
        self.inner.knows(name).py_err()
    }

    /// Return the half-life in seconds, or ``None`` if the nuclide is stable.
    ///
    /// Raises:
    ///     UnknownNuclideError: If the nuclide is in neither table, so
    ///         "stable" would be a guess rather than a fact.
    fn half_life(&self, name: &str) -> PyResult<Option<f64>> {
        self.inner.half_life(name).py_err()
    }

    /// Return the decay constant in 1/s, zero for a stable nuclide.
    fn decay_constant(&self, name: &str) -> PyResult<f64> {
        self.inner.decay_constant(name).py_err()
    }

    /// Return the relative atomic mass in u.
    ///
    /// A metastable state with no tabulated mass takes its ground state's, which
    /// is accurate to better than 1e-5.
    ///
    /// Raises:
    ///     UnknownNuclideError: If neither the nuclide nor its ground state has
    ///         a tabulated mass, since a guess would corrupt the mass of the
    ///         whole material.
    fn atomic_mass(&self, name: &str) -> PyResult<f64> {
        self.inner.atomic_mass(name).py_err()
    }

    /// Whether the nuclide has a tabulated half-life.
    fn is_radioactive(&self, name: &str) -> PyResult<bool> {
        self.inner.is_radioactive(name).py_err()
    }

    /// Fraction of this nuclide's decays that proceed by alpha emission.
    ///
    /// Returns 0.0 for a nuclide with no alpha branch, which includes every
    /// metastable state. Bi-212, which branches 35.94 percent alpha, returns
    /// 0.3594.
    fn alpha_fraction(&self, name: &str) -> PyResult<f64> {
        self.inner.alpha_fraction(name).py_err()
    }

    fn __repr__(&self) -> String {
        format!(
            "DecayData({} half-lives, {} masses)",
            self.inner.half_lives().len(),
            self.inner.atomic_masses().len()
        )
    }
}

/// Return the shared default tables.
#[pyfunction]
fn default_decay_data() -> PyDecayData {
    PyDecayData {
        inner: DecayData::shared_default(),
    }
}

/// Half-life in seconds from the default tables, ``None`` if stable.
#[pyfunction]
fn half_life(name: &str) -> PyResult<Option<f64>> {
    rmcf::decay::half_life(name).py_err()
}

/// Decay constant in 1/s from the default tables, zero if stable.
#[pyfunction]
fn decay_constant(name: &str) -> PyResult<f64> {
    rmcf::decay::decay_constant(name).py_err()
}

/// Relative atomic mass in u from the default tables.
#[pyfunction]
fn atomic_mass(name: &str) -> PyResult<f64> {
    rmcf::decay::atomic_mass(name).py_err()
}

/// Whether the nuclide is radioactive according to the default tables.
#[pyfunction]
fn is_radioactive(name: &str) -> PyResult<bool> {
    rmcf::decay::is_radioactive(name).py_err()
}

/// Fraction of decays proceeding by alpha emission, from the default tables.
#[pyfunction]
fn alpha_fraction(name: &str) -> PyResult<f64> {
    rmcf::decay::alpha_fraction(name).py_err()
}

// ----------------------------------------------------------------------
// material
// ----------------------------------------------------------------------

/// A nuclide inventory to assess for clearance.
///
/// The primary input is a mapping of nuclide name to atom count, which is what
/// an activation or depletion calculation produces:
///
///     >>> mat = Material({"Fe56": 8.4e22, "Co60": 1.2e12})
///
/// Specific activity in Bq/g is scale invariant, so atom counts, atom
/// densities and atom fractions all give the same answer for the Bq/g limit
/// sets and no density is needed. Volumetric limit sets (Ci/m3) do need one,
/// supplied as ``density`` or derived from atom densities.
///
/// Warning:
///     The Bq/g denominator is the mass of **everything** in the mapping, so
///     stable isotopes must be included. Passing only the radioactive nuclides
///     of an activated steel gives a mass thousands of times too small and a
///     specific activity thousands of times too high. An inventory with more
///     than one nuclide and no stable ones warns for that reason.
///
/// Args:
///     atoms: Nuclide name to atom count. Names may be spelled in any form
///         `normalise` accepts, and duplicates after normalising are summed.
///     density: Mass density in g/cm3. Only needed for volumetric limit sets.
///     volume: Volume in cm3. Only needed for total activity in Bq or Ci.
///     name: A label carried through to results, for reporting.
///     decay_data: Half-life and mass tables. Defaults to the vendored
///         ENDF/B-VIII.0 and AME2020 values.
#[pyclass(
    name = "Material",
    module = "radiological_material_clearance_finder._core"
)]
struct PyMaterial {
    inner: Material,
}

impl PyMaterial {
    fn finish(
        py: Python<'_>,
        material: rmcf::Result<Material>,
        density: Option<f64>,
        volume: Option<f64>,
        name: String,
        decay_data: Option<PyRef<'_, PyDecayData>>,
    ) -> PyResult<Self> {
        let mut material = material.py_err()?.with_name(name);
        if let Some(data) = decay_data {
            material = material.with_decay_data(data.inner.clone());
        }
        if let Some(d) = density {
            material = material.with_density(d).py_err()?;
        }
        if let Some(v) = volume {
            material = material.with_volume(v).py_err()?;
        }
        if material.looks_truncated() {
            warn(
                py,
                "every nuclide in this material is radioactive, so the total mass used \
                 for Bq/g is only the radioactive mass. If this inventory was filtered to \
                 its radioactive nuclides, add the stable ones back or the specific \
                 activity will be far too high.",
            )?;
        }
        Ok(PyMaterial { inner: material })
    }
}

#[pymethods]
impl PyMaterial {
    #[new]
    #[pyo3(signature = (atoms=None, *, density=None, volume=None, name=String::new(), decay_data=None))]
    fn new(
        py: Python<'_>,
        atoms: Option<BTreeMap<String, f64>>,
        density: Option<f64>,
        volume: Option<f64>,
        name: String,
        decay_data: Option<PyRef<'_, PyDecayData>>,
    ) -> PyResult<Self> {
        let atoms = atoms.ok_or_else(|| {
            PyValueError::new_err(
                "a material needs atom amounts. Use from_specific_activities to build one \
                 from Bq/g values instead.",
            )
        })?;
        Self::finish(
            py,
            Material::from_atom_counts(atoms),
            density,
            volume,
            name,
            decay_data,
        )
    }

    /// Build from absolute atom counts. Identical to the constructor.
    ///
    /// Args:
    ///     atoms: Nuclide name to atom count.
    ///     density: Mass density in g/cm3, only needed for volumetric sets.
    ///     volume: Volume in cm3, only needed for total activity.
    ///     name: A label carried through to results.
    ///     decay_data: Half-life and mass tables to use instead of the defaults.
    #[staticmethod]
    #[pyo3(signature = (atoms, *, density=None, volume=None, name=String::new(), decay_data=None))]
    fn from_atom_counts(
        py: Python<'_>,
        atoms: BTreeMap<String, f64>,
        density: Option<f64>,
        volume: Option<f64>,
        name: String,
        decay_data: Option<PyRef<'_, PyDecayData>>,
    ) -> PyResult<Self> {
        Self::finish(
            py,
            Material::from_atom_counts(atoms),
            density,
            volume,
            name,
            decay_data,
        )
    }

    /// Build from atom densities in atoms per barn-cm, OpenMC's unit.
    ///
    /// The mass density follows from the atom densities and the atomic masses,
    /// so no ``density`` argument is accepted or needed: supplying one could
    /// only contradict the inventory.
    ///
    /// Args:
    ///     densities: Nuclide name to atoms/barn-cm.
    ///     volume: Volume in cm3, only needed for total activity.
    ///     density: Not accepted. Present only so that passing it raises rather
    ///         than being swallowed and ignored.
    ///     name: A label carried through to results.
    ///     decay_data: Half-life and mass tables to use instead of the defaults.
    ///
    /// Raises:
    ///     TypeError: If ``density`` is given.
    #[staticmethod]
    #[pyo3(signature = (densities, *, volume=None, density=None, name=String::new(), decay_data=None))]
    fn from_atom_densities(
        py: Python<'_>,
        densities: BTreeMap<String, f64>,
        volume: Option<f64>,
        density: Option<f64>,
        name: String,
        decay_data: Option<PyRef<'_, PyDecayData>>,
    ) -> PyResult<Self> {
        if density.is_some() {
            return Err(PyTypeError::new_err(
                "from_atom_densities derives the mass density from the atom densities and \
                 the atomic masses, so passing density= could only contradict the \
                 inventory. Drop it, or use from_specific_activities if the density is the \
                 value you trust.",
            ));
        }
        Self::finish(
            py,
            Material::from_atom_densities(densities),
            None,
            volume,
            name,
            decay_data,
        )
    }

    /// Build from a mass in grams per nuclide.
    ///
    /// Args:
    ///     masses: Nuclide name to mass in grams.
    ///     density: Mass density in g/cm3, only needed for volumetric sets.
    ///     volume: Volume in cm3. The masses already fix the total mass, so a
    ///         volume that disagrees with density times volume is an error.
    ///     name: A label carried through to results.
    ///     decay_data: Half-life and mass tables to use instead of the defaults.
    #[staticmethod]
    #[pyo3(signature = (masses, *, density=None, volume=None, name=String::new(), decay_data=None))]
    fn from_masses(
        py: Python<'_>,
        masses: BTreeMap<String, f64>,
        density: Option<f64>,
        volume: Option<f64>,
        name: String,
        decay_data: Option<PyRef<'_, PyDecayData>>,
    ) -> PyResult<Self> {
        Self::finish(
            py,
            Material::from_masses(masses),
            density,
            volume,
            name,
            decay_data,
        )
    }

    /// Build from mass fractions, which need not sum to one.
    ///
    /// Args:
    ///     fractions: Nuclide name to mass fraction or weight percent.
    ///     density: Mass density in g/cm3, only needed for volumetric sets.
    ///     volume: Volume in cm3, only needed for total activity.
    ///     name: A label carried through to results.
    ///     decay_data: Half-life and mass tables to use instead of the defaults.
    #[staticmethod]
    #[pyo3(signature = (fractions, *, density=None, volume=None, name=String::new(), decay_data=None))]
    fn from_mass_fractions(
        py: Python<'_>,
        fractions: BTreeMap<String, f64>,
        density: Option<f64>,
        volume: Option<f64>,
        name: String,
        decay_data: Option<PyRef<'_, PyDecayData>>,
    ) -> PyResult<Self> {
        Self::finish(
            py,
            Material::from_mass_fractions(fractions),
            density,
            volume,
            name,
            decay_data,
        )
    }

    /// Build from specific activities in Bq/g, as an assay reports them.
    ///
    /// No half-life or atomic mass data is used, since the specific activity is
    /// the direct input to the index. Volumetric limit sets then need an
    /// explicit ``density``, because activity alone does not imply a mass.
    ///
    /// Args:
    ///     activities: Nuclide name to specific activity in Bq/g.
    ///     density: Mass density in g/cm3, only needed for volumetric sets.
    ///     volume: Volume in cm3, only needed for total activity.
    ///     name: A label carried through to results.
    ///     decay_data: Half-life tables, used by the scope and dynamic rules.
    #[staticmethod]
    #[pyo3(signature = (activities, *, density=None, volume=None, name=String::new(), decay_data=None))]
    fn from_specific_activities(
        py: Python<'_>,
        activities: BTreeMap<String, f64>,
        density: Option<f64>,
        volume: Option<f64>,
        name: String,
        decay_data: Option<PyRef<'_, PyDecayData>>,
    ) -> PyResult<Self> {
        Self::finish(
            py,
            Material::from_specific_activities(activities),
            density,
            volume,
            name,
            decay_data,
        )
    }

    /// A label carried through to results, for reporting.
    #[getter]
    fn name(&self) -> String {
        self.inner.name.clone()
    }

    #[setter]
    fn set_name(&mut self, name: String) {
        self.inner.name = name;
    }

    /// Volume in cm3, or ``None``.
    #[getter]
    fn volume(&self) -> Option<f64> {
        self.inner.volume()
    }

    #[setter]
    fn set_volume(&mut self, volume: Option<f64>) -> PyResult<()> {
        self.inner.set_volume(volume).py_err()
    }

    /// The half-life and mass tables this material uses.
    #[getter]
    fn decay_data(&self) -> PyDecayData {
        PyDecayData {
            inner: self.inner.decay_data().clone(),
        }
    }

    /// The nuclides present, in canonical form, sorted.
    #[getter]
    fn nuclides<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.inner.nuclides())
    }

    /// Atom amounts, keyed by canonical name.
    ///
    /// Raises:
    ///     InsufficientDataError: If the material was built from specific
    ///         activities, which do not determine atom counts.
    #[getter]
    fn atoms(&self) -> PyResult<BTreeMap<String, f64>> {
        self.inner.atoms().py_err()
    }

    /// Total mass in grams.
    ///
    /// Raises:
    ///     InsufficientDataError: If the material carries only relative
    ///         amounts, such as atom densities or mass fractions, with no
    ///         volume to scale them by.
    ///     ValueError: If the atom amounts and density times volume disagree
    ///         by more than one percent.
    #[getter]
    fn mass(&self) -> PyResult<f64> {
        self.inner.mass().py_err()
    }

    /// Mass density in g/cm3.
    ///
    /// Raises:
    ///     InsufficientDataError: If no density was supplied and none can be
    ///         derived, naming what to pass.
    #[getter]
    fn density(&self) -> PyResult<f64> {
        self.inner.mass_density().py_err()
    }

    /// Specific activity in Bq/g.
    ///
    /// Scale invariant, so this works from atom counts, atom densities or atom
    /// fractions alike.
    ///
    /// Args:
    ///     by_nuclide: Return a dict keyed by nuclide rather than the total.
    ///
    /// Returns:
    ///     Bq/g as a float, or a dict of them.
    #[pyo3(signature = (by_nuclide=false))]
    fn specific_activity(&self, py: Python<'_>, by_nuclide: bool) -> PyResult<Py<PyAny>> {
        self.activity(py, "Bq/g", by_nuclide)
    }

    /// Activity in the requested units.
    ///
    /// Args:
    ///     units: One of `ACTIVITY_UNITS`.
    ///     by_nuclide: Return a dict keyed by nuclide rather than the total.
    ///
    /// Returns:
    ///     The activity as a float, or a dict of them.
    ///
    /// Raises:
    ///     ValueError: If the units are not recognised.
    ///     InsufficientDataError: If the units need a density or volume the
    ///         material does not have.
    #[pyo3(signature = (units="Bq/g", by_nuclide=false))]
    fn activity(&self, py: Python<'_>, units: &str, by_nuclide: bool) -> PyResult<Py<PyAny>> {
        let units = self::units(units)?;
        if by_nuclide {
            let per_nuclide = self.inner.activities(units).py_err()?;
            Ok(per_nuclide.into_pyobject(py)?.into_any().unbind())
        } else {
            Ok(self
                .inner
                .activity(units)
                .py_err()?
                .into_pyobject(py)?
                .into_any()
                .unbind())
        }
    }

    fn __repr__(&self) -> String {
        let count = self.inner.nuclides().count();
        if self.inner.name.is_empty() {
            format!("Material({count} nuclides)")
        } else {
            format!("Material('{}', {count} nuclides)", self.inner.name)
        }
    }
}

// ----------------------------------------------------------------------
// limit sets
// ----------------------------------------------------------------------

/// A named table of per-nuclide activity limits.
///
/// Every nuclide name is canonicalised on construction, so a hand-written
/// ``"Co-60"`` behaves exactly like a loaded ``"Co60"``, and values that cannot
/// mean anything (a zero, negative or infinite limit, unknown units, a
/// non-positive threshold) raise ValueError.
///
/// Args:
///     name: Stable identifier, such as ``UK_EPR16_out_of_scope``.
///     label: Human readable description for reports.
///     units: Activity units the limits are in: ``Bq/g``, ``Ci/m3``, or ``Bq``
///         for a total activity limit, which needs the material's mass.
///     limits: Limit per nuclide name. A whole chain "sec" value goes under a
///         ``_sec`` key, such as ``U238_sec``.
///     jurisdiction: Issuing authority, such as ``UK`` or ``Germany``.
///     default_limit: Limit applied to a nuclide absent from ``limits``, or
///         ``None`` where the regulation has no catch-all. Three UK sets define
///         one: 0.01 Bq/g for ``UK_EPR16_out_of_scope`` and
///         ``UK_IRR17_notification``, and 0.1 Bq/g for
///         ``UK_IRR17_registration``.
///     unlimited: Nuclides the source explicitly places no limit on. They are
///         covered, and contribute nothing, which is different from being
///         absent.
///     secular_equilibrium: Parent to the daughters whose contribution the
///         parent's limit already includes, the "+" rows.
///     secular_equilibrium_sec: Parent to the daughters covered by its whole
///         chain "sec" value, usually a much longer list with a much stricter
///         limit. UK EPR 2016 gives U-238 three progeny at 1 Bq/g under
///         "U-238+" and fourteen at 0.01 Bq/g under "U-238sec".
///     limits_secular_equilibrium: The limit for a parent listed twice, plain
///         and marked "+", when its daughters are actually present. StrlSchV
///         gives Th-232 as 10 Bq/g plain and 0.01 Bq/g marked.
///     metal_overrides: Limits replacing or adding to ``limits`` when the
///         material is activated metal.
///     limits_per_gram: Limits in nCi/g, converted using the material density.
///     limits_upper: Upper end of a limit given as a range in the source, kept
///         for reference. ``limits`` holds the conservative lower end.
///     dynamic_rule: Name of a rule that depends on the material's own
///         nuclides, one of `DYNAMIC_RULES`.
///     min_half_life_scope: Half-life in seconds below which, if *every*
///         radionuclide present falls under it, the material is outside the
///         regulation altogether. A whole-material test, not a per-nuclide
///         filter.
///     threshold: Index value at or above which the material fails, normally 1.
///     source: The regulation the set is taken from.
///     url: Where the source was retrieved from.
///     retrieved: When it was retrieved.
///     notes: Anything else worth knowing.
#[pyclass(
    name = "LimitSet",
    module = "radiological_material_clearance_finder._core",
    frozen
)]
struct PyLimitSet {
    inner: Arc<LimitSet>,
}

fn daughters_to_dict<'py>(
    py: Python<'py>,
    map: &BTreeMap<String, Vec<String>>,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for (parent, daughters) in map {
        dict.set_item(parent, strings_to_tuple(py, daughters)?)?;
    }
    Ok(dict)
}

#[pymethods]
impl PyLimitSet {
    #[new]
    #[pyo3(signature = (
        name, label, units, limits, jurisdiction=String::new(), default_limit=None,
        unlimited=Vec::new(), secular_equilibrium=BTreeMap::new(),
        secular_equilibrium_sec=BTreeMap::new(), limits_secular_equilibrium=BTreeMap::new(),
        metal_overrides=BTreeMap::new(), limits_per_gram=BTreeMap::new(),
        limits_upper=BTreeMap::new(), dynamic_rule=None, min_half_life_scope=None,
        threshold=1.0, source=String::new(), url=String::new(), retrieved=String::new(),
        notes=String::new(),
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        name: String,
        label: String,
        units: String,
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
        dynamic_rule: Option<String>,
        min_half_life_scope: Option<f64>,
        threshold: f64,
        source: String,
        url: String,
        retrieved: String,
        notes: String,
    ) -> PyResult<Self> {
        let data = LimitSetData {
            name,
            label,
            units,
            limits,
            jurisdiction,
            default_limit,
            unlimited,
            secular_equilibrium,
            secular_equilibrium_sec,
            limits_secular_equilibrium,
            metal_overrides,
            limits_per_gram,
            limits_upper,
            dynamic_rule,
            min_half_life_scope,
            threshold,
            source,
            url,
            retrieved,
            notes,
        };
        Ok(PyLimitSet {
            inner: Arc::new(LimitSet::new(data).py_err()?),
        })
    }

    /// Stable identifier, such as ``UK_EPR16_out_of_scope``.
    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }
    /// Human readable description for reports.
    #[getter]
    fn label(&self) -> &str {
        self.inner.label()
    }
    /// Activity units the limits are in.
    #[getter]
    fn units(&self) -> &'static str {
        self.inner.units().as_str()
    }
    /// Limit per canonical nuclide name, with whole chain values under ``_sec``.
    #[getter]
    fn limits(&self) -> BTreeMap<String, f64> {
        self.inner.limits().clone()
    }
    /// Issuing authority.
    #[getter]
    fn jurisdiction(&self) -> &str {
        self.inner.jurisdiction()
    }
    /// The catch-all limit for unlisted nuclides, or ``None``.
    #[getter]
    fn default_limit(&self) -> Option<f64> {
        self.inner.default_limit()
    }
    /// Nuclides the source explicitly places no limit on.
    #[getter]
    fn unlimited<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        strings_to_tuple(py, self.inner.unlimited())
    }
    /// The "+" rows: parent to the daughters its limit covers.
    #[getter]
    fn secular_equilibrium<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        daughters_to_dict(py, self.inner.secular_equilibrium())
    }
    /// The "sec" rows: parent to the daughters its whole chain value covers.
    #[getter]
    fn secular_equilibrium_sec<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        daughters_to_dict(py, self.inner.secular_equilibrium_sec())
    }
    /// The limit for a parent listed twice, applied when its daughters are present.
    #[getter]
    fn limits_secular_equilibrium(&self) -> BTreeMap<String, f64> {
        self.inner.limits_secular_equilibrium().clone()
    }
    /// Limits replacing or adding to ``limits`` for activated metal.
    #[getter]
    fn metal_overrides(&self) -> BTreeMap<String, f64> {
        self.inner.metal_overrides().clone()
    }
    /// Limits in nCi/g.
    #[getter]
    fn limits_per_gram(&self) -> BTreeMap<String, f64> {
        self.inner.limits_per_gram().clone()
    }
    /// Upper end of limits the source gives as a range.
    #[getter]
    fn limits_upper(&self) -> BTreeMap<String, f64> {
        self.inner.limits_upper().clone()
    }
    /// Name of the rule that depends on the material's own nuclides, or ``None``.
    #[getter]
    fn dynamic_rule(&self) -> Option<&'static str> {
        self.inner.dynamic_rule().map(|r| r.as_str())
    }
    /// The whole-material short half-life scope test, in seconds, or ``None``.
    #[getter]
    fn min_half_life_scope(&self) -> Option<f64> {
        self.inner.min_half_life_scope()
    }
    /// The index value at or above which a material fails.
    #[getter]
    fn threshold(&self) -> f64 {
        self.inner.threshold()
    }
    /// The regulation this set is taken from.
    #[getter]
    fn source(&self) -> &str {
        self.inner.source()
    }
    /// Where the source was retrieved from.
    #[getter]
    fn url(&self) -> &str {
        self.inner.url()
    }
    /// When the source was retrieved.
    #[getter]
    fn retrieved(&self) -> &str {
        self.inner.retrieved()
    }
    /// Anything else worth knowing about the set.
    #[getter]
    fn notes(&self) -> &str {
        self.inner.notes()
    }

    /// Daughters whose activity this set's limit for ``parent`` already covers.
    fn daughters_of<'py>(&self, py: Python<'py>, parent: &str) -> PyResult<Bound<'py, PyTuple>> {
        strings_to_tuple(py, self.inner.daughters_of(parent))
    }

    /// Every nuclide this set names, whether limited or explicitly unlimited.
    ///
    /// The ``_sec`` keys are whole chain variants of a nuclide already counted,
    /// not nuclides of their own, so they are not included.
    #[getter]
    fn covered_nuclides<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyFrozenSet>> {
        PyFrozenSet::new(py, self.inner.covered_nuclides())
    }

    /// The number of nuclides limited, not counting whole chain variants.
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other
            .extract::<PyRef<'_, PyLimitSet>>()
            .is_ok_and(|o| *o.inner == *self.inner)
    }

    /// Hash by name and units, which is what the registry identifies a set by.
    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.name().hash(&mut hasher);
        self.inner.units().hash(&mut hasher);
        hasher.finish()
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }
}

/// Resolve a limit set argument that may be a name or a `LimitSet`.
fn resolve_set(obj: &Bound<'_, PyAny>) -> PyResult<Arc<LimitSet>> {
    if let Ok(set) = obj.extract::<PyRef<'_, PyLimitSet>>() {
        return Ok(set.inner.clone());
    }
    let name: String = obj
        .extract()
        .map_err(|_| PyTypeError::new_err("limit_set must be a registered name or a LimitSet"))?;
    rmcf::get_limit_set(&name).py_err()
}

/// Names of the available limit sets, sorted.
///
/// Args:
///     jurisdiction: Restrict to one issuing authority, such as ``UK``.
#[pyfunction]
#[pyo3(signature = (jurisdiction=None))]
fn limit_sets<'py>(py: Python<'py>, jurisdiction: Option<&str>) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, rmcf::limit_sets(jurisdiction))
}

/// Look up a limit set by name.
///
/// Args:
///     name: A registered name, or an already built `LimitSet`, which is
///         returned as is.
///
/// Raises:
///     KeyError: If no such set is registered, listing the ones that are.
#[pyfunction]
fn get_limit_set(name: &Bound<'_, PyAny>) -> PyResult<PyLimitSet> {
    Ok(PyLimitSet {
        inner: resolve_set(name)?,
    })
}

/// Add a limit set at runtime, for a site-specific or draft table.
///
/// Args:
///     limit_set: The set to register. Replaces any set of the same name,
///         warning first if that name came from the shipped regulatory data,
///         since shadowing a published table by accident would be hard to spot
///         in a result that only records the name.
#[pyfunction]
fn register_limit_set(py: Python<'_>, limit_set: PyRef<'_, PyLimitSet>) -> PyResult<()> {
    let name = limit_set.inner.name().to_string();
    if rmcf::register_limit_set((*limit_set.inner).clone()) {
        warn(
            py,
            &format!(
                "replacing the limit set '{name}', which came from the shipped regulatory \
                 data. Results will report that name while using the new table."
            ),
        )?;
    }
    Ok(())
}

// ----------------------------------------------------------------------
// the clearance index
// ----------------------------------------------------------------------

/// The outcome of assessing one material against one limit set.
///
/// The per-nuclide dicts are ordered as a report reads them, largest first.
#[pyclass(
    name = "ClearanceResult",
    module = "radiological_material_clearance_finder._core",
    frozen
)]
struct PyClearanceResult {
    inner: rmcf::ClearanceResult,
}

#[pymethods]
impl PyClearanceResult {
    /// Name of the set assessed against.
    #[getter]
    fn limit_set(&self) -> &str {
        &self.inner.limit_set
    }
    /// The sum of activity-to-limit ratios.
    #[getter]
    fn index(&self) -> f64 {
        self.inner.index
    }
    /// The value the index must stay below, normally 1.
    #[getter]
    fn threshold(&self) -> f64 {
        self.inner.threshold
    }
    /// Activity units the comparison was made in.
    #[getter]
    fn units(&self) -> &'static str {
        self.inner.units.as_str()
    }
    /// Each nuclide's contribution to the index, largest first.
    #[getter]
    fn by_nuclide<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        pairs_to_dict(py, &self.inner.by_nuclide)
    }
    /// Each nuclide's activity in ``units``, largest first.
    #[getter]
    fn activities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        pairs_to_dict(py, &self.inner.activities)
    }
    /// The limit applied to each nuclide, after any metal, per-gram, dynamic or
    /// secular equilibrium adjustment.
    #[getter]
    fn limits_used<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        pairs_to_dict(py, &self.inner.limits_used)
    }
    /// Nuclides that took the set's catch-all limit because the table does not
    /// list them.
    #[getter]
    fn defaulted<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        strings_to_tuple(py, &self.inner.defaulted)
    }
    /// Nuclide to the reason its activity was left out, which is always that a
    /// parent's limit already accounts for it in full.
    #[getter]
    fn excluded<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        pairs_to_dict(py, &self.inner.excluded)
    }
    /// Nuclide to the activity a parent accounted for, where the parent could
    /// only support part of it. The remainder was assessed against the
    /// nuclide's own limit, so for these the ratio in `by_nuclide` is computed
    /// from ``activities[name] - credited[name]``. `assessed_activity` returns
    /// that residual directly.
    #[getter]
    fn credited<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        pairs_to_dict(py, &self.inner.credited)
    }
    /// Activity present with no limit and no catch-all, so absent from the
    /// index entirely. The number to check before trusting a comfortable index.
    #[getter]
    fn uncovered<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        pairs_to_dict(py, &self.inner.uncovered)
    }
    /// Nuclides the source explicitly places no limit on.
    #[getter]
    fn unlimited<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        strings_to_tuple(py, &self.inner.unlimited)
    }
    /// Whether the regulation excludes this material outright, which for the UK
    /// sets means every radionuclide present is shorter lived than 100 seconds.
    #[getter]
    fn out_of_scope(&self) -> bool {
        self.inner.out_of_scope
    }
    /// The material's label, carried through for reporting.
    #[getter]
    fn material_name(&self) -> &str {
        &self.inner.material_name
    }
    /// Whether the material meets this set's limits.
    #[getter]
    fn clearable(&self) -> bool {
        self.inner.clearable()
    }
    /// Total activity with no limit, in this result's units.
    #[getter]
    fn uncovered_activity(&self) -> f64 {
        self.inner.uncovered_activity()
    }
    /// Share of total activity that falls outside the sum, from 0 to 1.
    ///
    /// A large value means the index understates the inventory: activity is
    /// present that no limit in this set applies to. Zero means every
    /// radionuclide present was either limited or explicitly unlimited.
    #[getter]
    fn uncovered_fraction(&self) -> f64 {
        self.inner.uncovered_fraction()
    }

    /// The activity actually charged against the limit for one nuclide.
    ///
    /// The same as its total activity, except where a parent accounted for part
    /// of it, in which case the remainder is what the ratio was computed from.
    ///
    /// Args:
    ///     nuclide: Canonical nuclide name.
    fn assessed_activity(&self, nuclide: &str) -> f64 {
        self.inner.assessed_activity(nuclide)
    }

    /// The nuclides contributing most to the index.
    ///
    /// Args:
    ///     count: How many to return.
    #[pyo3(signature = (count=10))]
    fn dominant(&self, count: usize) -> Vec<(String, f64)> {
        self.inner.dominant(count).to_vec()
    }

    /// A JSON-serialisable copy of the result.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let r = &self.inner;
        let d = PyDict::new(py);
        d.set_item("limit_set", &r.limit_set)?;
        d.set_item("material", &r.material_name)?;
        d.set_item("index", r.index)?;
        d.set_item("threshold", r.threshold)?;
        d.set_item("units", r.units.as_str())?;
        d.set_item("clearable", r.clearable())?;
        d.set_item("out_of_scope", r.out_of_scope)?;
        d.set_item("by_nuclide", pairs_to_dict(py, &r.by_nuclide)?)?;
        d.set_item("activities", pairs_to_dict(py, &r.activities)?)?;
        d.set_item("limits_used", pairs_to_dict(py, &r.limits_used)?)?;
        d.set_item("defaulted", r.defaulted.clone())?;
        d.set_item("excluded", pairs_to_dict(py, &r.excluded)?)?;
        d.set_item("credited", pairs_to_dict(py, &r.credited)?)?;
        d.set_item("uncovered", pairs_to_dict(py, &r.uncovered)?)?;
        d.set_item("unlimited", r.unlimited.clone())?;
        d.set_item("uncovered_fraction", r.uncovered_fraction())?;
        Ok(d)
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        format!(
            "ClearanceResult('{}', index={}, clearable={})",
            self.inner.limit_set,
            self.inner.index,
            if self.inner.clearable() {
                "True"
            } else {
                "False"
            }
        )
    }
}

fn options(metal: bool, apply_default_limit: bool, exclude_daughters: bool) -> ClearanceOptions {
    ClearanceOptions {
        metal,
        apply_default_limit,
        exclude_daughters,
    }
}

/// Assess a material against one set of clearance limits.
///
/// Args:
///     material: The inventory to assess.
///     limit_set: A registered limit set name, or a `LimitSet`.
///     metal: Whether the material is activated metal, which changes some NRC
///         limits and adds others.
///     apply_default_limit: Whether to apply the set's catch-all limit to
///         nuclides the table does not list. The UK regulations define one, so
///         leaving this on is what the regulation says; turning it off shows
///         what the listed nuclides alone contribute.
///     exclude_daughters: Whether to leave out daughters whose parent's limit
///         already covers them, per the regulation's own table. Turning this
///         off double counts them, which is conservative but not what the
///         regulation intends.
///
/// Returns:
///     A `ClearanceResult` carrying the index and everything needed to judge
///     it, including any activity that fell outside the sum.
///
/// Raises:
///     KeyError: If the limit set name is not registered.
///     InsufficientDataError: If a volumetric set is used and the material has
///         no density.
#[pyfunction]
#[pyo3(signature = (material, limit_set, *, metal=false, apply_default_limit=true, exclude_daughters=true))]
fn clearance_index(
    material: PyRef<'_, PyMaterial>,
    limit_set: &Bound<'_, PyAny>,
    metal: bool,
    apply_default_limit: bool,
    exclude_daughters: bool,
) -> PyResult<PyClearanceResult> {
    let set = resolve_set(limit_set)?;
    let opts = options(metal, apply_default_limit, exclude_daughters);
    Ok(PyClearanceResult {
        inner: rmcf::clearance_index(&material.inner, &set, opts).py_err()?,
    })
}

/// Assess a material against many limit sets at once.
///
/// Sets whose units the material cannot supply, such as a volumetric set for a
/// material with no density, are skipped rather than raising, so one missing
/// density does not hide every result that does not need it.
///
/// Args:
///     material: The inventory to assess.
///     names: Limit set names, defaulting to every registered set.
///     metal: As for `clearance_index`.
///     apply_default_limit: As for `clearance_index`.
///     exclude_daughters: As for `clearance_index`.
///
/// Returns:
///     Results keyed by limit set name.
#[pyfunction]
#[pyo3(signature = (material, names=None, *, metal=false, apply_default_limit=true, exclude_daughters=true))]
fn clearance_indices<'py>(
    py: Python<'py>,
    material: PyRef<'_, PyMaterial>,
    names: Option<Vec<String>>,
    metal: bool,
    apply_default_limit: bool,
    exclude_daughters: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let names: Option<Vec<&str>> = names
        .as_ref()
        .map(|n| n.iter().map(String::as_str).collect());
    let opts = options(metal, apply_default_limit, exclude_daughters);
    let results = rmcf::clearance_indices(&material.inner, names.as_deref(), opts).py_err()?;
    let dict = PyDict::new(py);
    for result in results {
        dict.set_item(
            result.limit_set.clone(),
            PyClearanceResult { inner: result },
        )?;
    }
    Ok(dict)
}

/// The limit sets this material already meets, best margin first.
///
/// Args:
///     material: The inventory to assess.
///     names: Limit set names, defaulting to every registered set.
///     metal: As for `clearance_index`.
///     apply_default_limit: As for `clearance_index`.
///     exclude_daughters: As for `clearance_index`.
///
/// Returns:
///     Names of the sets whose index is below their threshold.
#[pyfunction]
#[pyo3(signature = (material, names=None, *, metal=false, apply_default_limit=true, exclude_daughters=true))]
fn clearable_routes(
    material: PyRef<'_, PyMaterial>,
    names: Option<Vec<String>>,
    metal: bool,
    apply_default_limit: bool,
    exclude_daughters: bool,
) -> PyResult<Vec<String>> {
    let names: Option<Vec<&str>> = names
        .as_ref()
        .map(|n| n.iter().map(String::as_str).collect());
    let opts = options(metal, apply_default_limit, exclude_daughters);
    rmcf::clearable_routes(&material.inner, names.as_deref(), opts).py_err()
}

// ----------------------------------------------------------------------
// classification
// ----------------------------------------------------------------------

/// Classify a material for near-surface disposal under 10 CFR 61.55.
///
/// Applies the sum of fractions rule to Table 1 and Table 2 and combines the
/// results with the rules in paragraphs 61.55(a)(3) to (a)(7).
///
/// Args:
///     material: The inventory to classify. Needs a density, since the NRC
///         limits are volumetric.
///     metal: Whether the material is activated metal, which the regulation
///         treats with its own rows.
///
/// Returns:
///     ``"Class A"``, ``"Class B"``, ``"Class C"`` or ``"GTCC"`` for greater
///     than Class C.
#[pyfunction]
#[pyo3(signature = (material, *, metal=false))]
fn nrc_waste_class(material: PyRef<'_, PyMaterial>, metal: bool) -> PyResult<&'static str> {
    Ok(rmcf::nrc_waste_class(&material.inner, metal)
        .py_err()?
        .as_str())
}

/// Activity from alpha emission, weighted by each nuclide's alpha branch.
///
/// Bi-212 branches 35.94 percent alpha, so it contributes that share of its
/// activity here and the rest to `beta_gamma_activity`.
///
/// Args:
///     material: The inventory.
///     units: Any unit `Material.activity` accepts.
#[pyfunction]
#[pyo3(signature = (material, units="Bq/g"))]
fn alpha_activity(material: PyRef<'_, PyMaterial>, units: &str) -> PyResult<f64> {
    rmcf::alpha_activity(&material.inner, self::units(units)?).py_err()
}

/// Activity from everything that is not alpha emission.
///
/// Args:
///     material: The inventory.
///     units: Any unit `Material.activity` accepts.
#[pyfunction]
#[pyo3(signature = (material, units="Bq/g"))]
fn beta_gamma_activity(material: PyRef<'_, PyMaterial>, units: &str) -> PyResult<f64> {
    rmcf::beta_gamma_activity(&material.inner, self::units(units)?).py_err()
}

/// The UK waste category of a material, with the numbers behind it.
#[pyclass(
    name = "UKWasteCategory",
    module = "radiological_material_clearance_finder._core",
    frozen
)]
struct PyUkWasteCategory {
    inner: rmcf::UkWasteCategory,
}

#[pymethods]
impl PyUkWasteCategory {
    /// ``"VLLW"``, ``"LLW"`` or ``"ILW"``.
    #[getter]
    fn category(&self) -> &'static str {
        self.inner.category.as_str()
    }
    /// Alpha activity in Bq/g.
    #[getter]
    fn alpha(&self) -> f64 {
        self.inner.alpha
    }
    /// Beta and gamma activity in Bq/g.
    #[getter]
    fn beta_gamma(&self) -> f64 {
        self.inner.beta_gamma
    }
    /// Tritium plus carbon-14 activity in Bq/g, which share their own VLLW
    /// allowance.
    #[getter]
    fn tritium_and_c14(&self) -> f64 {
        self.inner.tritium_and_c14
    }
    /// Total activity in Bq/g.
    #[getter]
    fn total(&self) -> f64 {
        self.inner.total
    }
    /// Why this category and not the one below it.
    #[getter]
    fn reason(&self) -> &str {
        &self.inner.reason
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        format!(
            "UKWasteCategory('{}', total={} Bq/g)",
            self.inner.category, self.inner.total
        )
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
/// Two limits of this classification are worth stating rather than hiding.
/// The low volume VLLW category ("dustbin disposal") is defined per 0.1 cubic
/// metre and per item, so it is a property of a consignment rather than of a
/// material and is not decided here. The boundary above which waste becomes
/// high level rather than intermediate is a thermal one, about 2 kW/m3, and
/// needs decay heat rather than activity, so this returns ILW for everything
/// above LLW.
///
/// Args:
///     material: The inventory to classify.
///
/// Returns:
///     A `UKWasteCategory` holding the category and the activities it was
///     decided on.
#[pyfunction]
fn uk_waste_category(material: PyRef<'_, PyMaterial>) -> PyResult<PyUkWasteCategory> {
    Ok(PyUkWasteCategory {
        inner: rmcf::uk_waste_category(&material.inner).py_err()?,
    })
}

// ----------------------------------------------------------------------
// cooling
// ----------------------------------------------------------------------

fn series_items<'py>(series: &Bound<'py, PyAny>) -> PyResult<Vec<(f64, PyRef<'py, PyMaterial>)>> {
    series
        .call_method0("items")?
        .try_iter()?
        .map(|item| item?.extract())
        .collect()
}

/// Evaluate the clearance index at each cooling time.
///
/// Args:
///     series: Cooling time in seconds to the material at that time.
///     limit_set: A registered limit set name, or a `LimitSet`.
///     metal: As for `clearance_index`.
///     apply_default_limit: As for `clearance_index`.
///     exclude_daughters: As for `clearance_index`.
///
/// Returns:
///     Cooling time to index, ordered by time.
#[pyfunction]
#[pyo3(signature = (series, limit_set, *, metal=false, apply_default_limit=true, exclude_daughters=true))]
fn index_series<'py>(
    py: Python<'py>,
    series: &Bound<'py, PyAny>,
    limit_set: &Bound<'py, PyAny>,
    metal: bool,
    apply_default_limit: bool,
    exclude_daughters: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let set = resolve_set(limit_set)?;
    let items = series_items(series)?;
    let opts = options(metal, apply_default_limit, exclude_daughters);
    let indexes =
        rmcf::index_series(items.iter().map(|(t, m)| (*t, &m.inner)), &set, opts).py_err()?;
    let dict = PyDict::new(py);
    for (time, index) in indexes {
        dict.set_item(time, index)?;
    }
    Ok(dict)
}

/// Find the cooling time at which a material first meets a set of limits.
///
/// The index falls roughly exponentially with cooling time, so this
/// interpolates logarithmically in the index between the two samples that
/// bracket the threshold. Interpolating linearly instead would place the
/// crossing systematically late, by a factor that grows with the spacing of the
/// samples.
///
/// Args:
///     series: Cooling time in seconds to the material at that time. Times need
///         not be sorted or evenly spaced.
///     limit_set: A registered limit set name, or a `LimitSet`.
///     allow_ingrowth: Return the first crossing even when the index later
///         climbs back above the threshold, rather than raising
///         `IngrowthError`.
///     metal: As for `clearance_index`.
///     apply_default_limit: As for `clearance_index`.
///     exclude_daughters: As for `clearance_index`.
///
/// Returns:
///     The cooling time in seconds at which the index first falls below the
///     threshold. If the earliest sample already meets it, that sample's time
///     is returned rather than zero, since the series says nothing about
///     anything earlier. ``None`` if the series never meets it.
///
/// Raises:
///     ValueError: If fewer than two times are given.
///     IngrowthError: If the index climbs back above the threshold at a later
///         time, so the material does not stay clear.
#[pyfunction]
#[pyo3(signature = (series, limit_set, *, allow_ingrowth=false, metal=false, apply_default_limit=true, exclude_daughters=true))]
fn time_to_clear(
    series: &Bound<'_, PyAny>,
    limit_set: &Bound<'_, PyAny>,
    allow_ingrowth: bool,
    metal: bool,
    apply_default_limit: bool,
    exclude_daughters: bool,
) -> PyResult<Option<f64>> {
    let set = resolve_set(limit_set)?;
    let items = series_items(series)?;
    let opts = options(metal, apply_default_limit, exclude_daughters);
    rmcf::time_to_clear(
        items.iter().map(|(t, m)| (*t, &m.inner)),
        &set,
        opts,
        allow_ingrowth,
    )
    .py_err()
}

// ----------------------------------------------------------------------
// module
// ----------------------------------------------------------------------

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    m.add("NuclideNameError", py.get_type::<NuclideNameError>())?;
    m.add("UnknownNuclideError", py.get_type::<UnknownNuclideError>())?;
    m.add(
        "InsufficientDataError",
        py.get_type::<InsufficientDataError>(),
    )?;
    m.add("IngrowthError", py.get_type::<IngrowthError>())?;
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_function(wrap_pyfunction!(normalise, m)?)?;
    m.add_function(wrap_pyfunction!(parse_regulatory, m)?)?;
    m.add_function(wrap_pyfunction!(element, m)?)?;
    m.add_function(wrap_pyfunction!(mass_number, m)?)?;
    m.add_function(wrap_pyfunction!(metastable_state, m)?)?;
    m.add_function(wrap_pyfunction!(atomic_number, m)?)?;
    m.add_function(wrap_pyfunction!(is_valid, m)?)?;

    m.add("AVOGADRO", rmcf::AVOGADRO)?;
    m.add("BECQUEREL_PER_CURIE", rmcf::BECQUEREL_PER_CURIE)?;
    m.add_class::<PyDecayData>()?;
    m.add_function(wrap_pyfunction!(default_decay_data, m)?)?;
    m.add_function(wrap_pyfunction!(half_life, m)?)?;
    m.add_function(wrap_pyfunction!(decay_constant, m)?)?;
    m.add_function(wrap_pyfunction!(atomic_mass, m)?)?;
    m.add_function(wrap_pyfunction!(is_radioactive, m)?)?;
    m.add_function(wrap_pyfunction!(alpha_fraction, m)?)?;

    let units: Vec<&str> = ActivityUnit::ALL.iter().map(|u| u.as_str()).collect();
    m.add("ACTIVITY_UNITS", PyTuple::new(py, units)?)?;
    m.add_class::<PyMaterial>()?;

    let rules: Vec<&str> = rmcf::DynamicRule::ALL.iter().map(|r| r.as_str()).collect();
    m.add("DYNAMIC_RULES", PyTuple::new(py, rules)?)?;
    m.add_class::<PyLimitSet>()?;
    m.add_function(wrap_pyfunction!(limit_sets, m)?)?;
    m.add_function(wrap_pyfunction!(get_limit_set, m)?)?;
    m.add_function(wrap_pyfunction!(register_limit_set, m)?)?;

    m.add("EQUILIBRIUM_TOLERANCE", rmcf::EQUILIBRIUM_TOLERANCE)?;
    m.add_class::<PyClearanceResult>()?;
    m.add_function(wrap_pyfunction!(clearance_index, m)?)?;
    m.add_function(wrap_pyfunction!(clearance_indices, m)?)?;
    m.add_function(wrap_pyfunction!(clearable_routes, m)?)?;

    m.add(
        "UK_LLW_ALPHA_BQ_PER_G",
        rmcf::classify::UK_LLW_ALPHA_BQ_PER_G,
    )?;
    m.add(
        "UK_LLW_BETA_GAMMA_BQ_PER_G",
        rmcf::classify::UK_LLW_BETA_GAMMA_BQ_PER_G,
    )?;
    m.add(
        "UK_HIGH_VOLUME_VLLW_BQ_PER_G",
        rmcf::classify::UK_HIGH_VOLUME_VLLW_BQ_PER_G,
    )?;
    m.add(
        "UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G",
        rmcf::classify::UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G,
    )?;
    m.add_class::<PyUkWasteCategory>()?;
    m.add_function(wrap_pyfunction!(nrc_waste_class, m)?)?;
    m.add_function(wrap_pyfunction!(uk_waste_category, m)?)?;
    m.add_function(wrap_pyfunction!(alpha_activity, m)?)?;
    m.add_function(wrap_pyfunction!(beta_gamma_activity, m)?)?;

    m.add_function(wrap_pyfunction!(index_series, m)?)?;
    m.add_function(wrap_pyfunction!(time_to_clear, m)?)?;
    Ok(())
}
