# radiological-material-clearance-finder

Clearance indexes for a nuclide inventory, against the UK, German, US, EU and
IAEA limit sets. Give it a map of nuclides and atom amounts and it tells you
whether the material can be cleared, by which route, and what is driving the
answer.

This is the Rust core. The same library is on PyPI as
[`radiological-material-clearance-finder`](https://pypi.org/project/radiological-material-clearance-finder/),
built from this crate with pyo3.

The 22 regulatory tables and the ENDF/B-VIII.0 half-lives and AME2020 masses
are compiled in, so there is nothing to download or configure. There is no
pyo3 dependency here, so a transport or transmutation code can link it into its
own material type without pulling Python in.

```toml
[dependencies]
radiological-material-clearance-finder = "0.1"
```

```rust
use radiological_material_clearance_finder::{
    clearance_index, get_limit_set, ClearanceOptions, Material,
};

let steel = Material::from_atom_counts([("Fe56", 8.4e22), ("Co60", 2.1e9), ("Cs137", 3.1e8)])?;
let result = clearance_index(&steel, &get_limit_set("UK_EPR16_out_of_scope")?, ClearanceOptions::default())?;

result.index;        // 11.244691988476553
result.clearable();  // false
result.uncovered;    // activity no limit in this set applies to
```

## From a code's own material

A material that knows its atom densities in atoms per barn-cm (OpenMC's, yamc's
and yani's unit) hands them over directly. The mass density follows from the
atomic masses, so volumetric (Ci/m3) sets work too, and a volume enables total
activity:

```rust,ignore
let inventory = Material::from_atom_densities(material.get_atoms_per_barn_cm()?)?
    .with_volume(volume_cm3)?
    .with_name(name);
let results = clearance_indices(&inventory, None, ClearanceOptions::default())?;
```

A code with its own decay data, so that the half-lives match the ones its
transmutation used, passes them in with `DecayData::new` and
`Material::with_decay_data`. An OpenMC chain file works too, through
`DecayData::from_chain_xml`.

See the [documentation](https://fusion-energy.github.io/radiological_material_clearance_finder/)
for what each limit set means and where every number comes from.

## Licence

MIT
