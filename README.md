# radiological material clearance finder

[![test](https://github.com/fusion-energy/radiological_material_clearance_finder/actions/workflows/test.yml/badge.svg)](https://github.com/fusion-energy/radiological_material_clearance_finder/actions/workflows/test.yml)
[![docs](https://github.com/fusion-energy/radiological_material_clearance_finder/actions/workflows/docs.yml/badge.svg)](https://github.com/fusion-energy/radiological_material_clearance_finder/actions/workflows/docs.yml)
[![PyPI](https://img.shields.io/pypi/v/radiological-material-clearance-finder.svg)](https://pypi.org/project/radiological-material-clearance-finder/)
[![python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://pypi.org/project/radiological-material-clearance-finder/)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Clearance indexes for a nuclide inventory, against the UK, German, US, EU and
IAEA limit sets. Give it a dictionary of isotopes and atom numbers and it tells
you whether the material can be cleared, by which route, and what is driving the
answer.

No runtime dependencies. OpenMC is not required.

### 📖 [Read the documentation](https://fusion-energy.github.io/radiological_material_clearance_finder/)

---

```bash
pip install radiological-material-clearance-finder
```

```python
from radiological_material_clearance_finder import Material, clearance_index

steel = Material({"Fe56": 8.4e22, "Co60": 2.1e9, "Cs137": 3.1e8})
result = clearance_index(steel, "UK_EPR16_out_of_scope")

result.index      # 11.244691988476553
result.clearable  # False
result.uncovered  # activity no limit in this set applies to
```

22 limit sets covering clearance, exemption and disposal classification, all
generated from the official regulatory sources and carrying their provenance.

| | |
| --- | --- |
| [Getting started](https://fusion-energy.github.io/radiological_material_clearance_finder/getting-started/) | Installing, the five ways to build a material, reading a result |
| [Limit sets](https://fusion-energy.github.io/radiological_material_clearance_finder/limit-sets/) | All 22 sets, what each one means, and how to inspect its source |
| [Regulatory notes](https://fusion-energy.github.io/radiological_material_clearance_finder/regulatory-notes/) | Secular equilibrium, catch-all limits, scope rules |
| [Using OpenMC](https://fusion-energy.github.io/radiological_material_clearance_finder/openmc/) | Materials and depletions, and cooling time to clearance |
| [Where the data comes from](https://fusion-energy.github.io/radiological_material_clearance_finder/data-provenance/) | Every source, and what has and has not been verified |
| [API reference](https://fusion-energy.github.io/radiological_material_clearance_finder/api/) | Rendered from the docstrings |

## Contributing

```bash
git clone https://github.com/fusion-energy/radiological_material_clearance_finder.git
cd radiological_material_clearance_finder
pip install -e ".[test]"
pytest
```

The regulatory tables are generated, not hand written. See
[where the data comes from](https://fusion-energy.github.io/radiological_material_clearance_finder/data-provenance/)
before editing anything under `src/radiological_material_clearance_finder/data/`.

## Licence

[MIT](LICENSE)
