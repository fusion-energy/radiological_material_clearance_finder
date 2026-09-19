# API reference

Everything in the public API, rendered from the docstrings in the source.

The whole of it is importable from the top-level package:

```python
from radiological_material_clearance_finder import Material, clearance_index
```

| Module | What is in it |
| --- | --- |
| [Material](material.md) | The nuclide inventory being assessed, and its five constructors. |
| [index](clearance-index.md) | The clearance index itself, and the result it returns. |
| [limits](limits.md) | Limit sets: the regulatory tables a material is measured against. |
| [classify](classify.md) | Waste classification, where the answer is a category rather than an index. |
| [cooling](cooling.md) | Finding the cooling time at which a material becomes clearable. |
| [decay](decay.md) | Half-lives, decay constants, atomic masses and alpha branching. |
| [nuclide](nuclide.md) | Nuclide name parsing, including the spellings the regulations use. |
| [openmc interop](openmc_interop.md) | Getting an inventory out of OpenMC without depending on it. |

The one exception is `openmc_interop`, which imports OpenMC inside its functions
and is therefore imported from its own module rather than from the top level:

```python
from radiological_material_clearance_finder.openmc_interop import from_openmc_material
```
