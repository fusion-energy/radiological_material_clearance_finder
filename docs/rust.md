# Using it from Rust

The calculation is a Rust crate, `radiological-material-clearance-finder`, and
the Python package is a thin pyo3 binding over it. The crate has no Python in
it, so a transport or transmutation code can link it into its own material type
and report clearance indexes without a round trip through Python.

```toml
[dependencies]
radiological-material-clearance-finder = "0.1"
```

The regulatory tables, and the ENDF/B-VIII.0 half-lives and AME2020 masses, are
compiled into the crate. There is nothing to download or point at.

## The same example

```rust
use radiological_material_clearance_finder::{
    clearance_index, get_limit_set, ClearanceOptions, Material,
};

let steel = Material::from_atom_counts([("Fe56", 8.4e22), ("Co60", 2.1e9), ("Cs137", 3.1e8)])?;
let set = get_limit_set("UK_EPR16_out_of_scope")?;
let result = clearance_index(&steel, &set, ClearanceOptions::default())?;

assert_eq!(result.index, 11.244691988476553);
assert!(!result.clearable());
println!("{result}");
```

It gives the same number, to the last bit, as the Python example in
[getting started](getting-started.md), because it is the same code.

## From a code's own material

yamc and yani, like OpenMC, hold a material's composition as atom densities in
atoms per barn-cm. Those go straight in, and the mass density follows from the
atomic masses, so the volumetric US sets work as well as the Bq/g ones:

```rust
let inventory = Material::from_atom_densities(material.get_atoms_per_barn_cm()?)?
    .with_volume(volume_cm3)?     // only needed for total activity in Bq or Ci
    .with_name(name);

// Every set the material can be assessed against, skipping those it cannot.
let results = clearance_indices(&inventory, None, ClearanceOptions::default())?;
```

A code whose transmutation used its own half-lives should assess with the same
ones, or the index is computed from a different decay constant than the
inventory was:

```rust
let data = Arc::new(DecayData::new(my_half_lives, my_atomic_masses)?);
let inventory = Material::from_atom_densities(densities)?.with_decay_data(data);
```

`DecayData::from_chain_xml` reads an OpenMC depletion chain for the same
purpose.

## What maps to what

| Python | Rust |
| --- | --- |
| `Material(atoms, density=, volume=)` | `Material::from_atom_counts(atoms)?.with_density(d)?.with_volume(v)?` |
| `Material.from_atom_densities`, `from_masses`, ... | `Material::from_atom_densities`, `from_masses`, ... |
| `material.activity("Ci/m3", by_nuclide=True)` | `material.activities(ActivityUnit::CiPerM3)?` |
| `clearance_index(m, "name", metal=True)` | `clearance_index(&m, &get_limit_set("name")?, ClearanceOptions { metal: true, ..Default::default() })?` |
| `LimitSet(name=..., limits=...)` | `LimitSet::new(LimitSetData { .. })?` |
| `register_limit_set(s)` (warns on shadowing) | `register_limit_set(s)` returns `true` on shadowing |
| `time_to_clear({t: m, ...}, set)` | `time_to_clear([(t, &m), ...], &set, options, allow_ingrowth)?` |
| `nrc_waste_class(m)` returns `"Class A"` | `nrc_waste_class(&m, metal)?` returns `NrcWasteClass::A` |
| exceptions | one `Error` enum, with a variant per exception class |

`ClearanceResult` has the same fields. The per-nuclide ones are `Vec<(String,
f64)>` in the same order the Python dicts have, largest first, and the struct
implements `serde::Serialize` for writing results out.

Full API documentation is on [docs.rs](https://docs.rs/radiological-material-clearance-finder).
