# Using it with OpenMC

OpenMC is not a dependency. It is imported inside the two interop functions and
nowhere else, so installing this package pulls in nothing and every limit set
works without it.

```python
from radiological_material_clearance_finder.openmc_interop import from_openmc_material

material = from_openmc_material(openmc_material)
```

That reads the atom densities in atoms per barn-cm, which fixes the mass density
as well, so the result supports both the Bq/g sets and the volumetric US sets
with nothing else supplied. The OpenMC material's name and volume are carried
across where it has them.

## From a depletion

`from_depletion_results` turns a depletion straight into the series that
[`time_to_clear`](api/cooling.md) expects:

```python
import openmc.deplete
from radiological_material_clearance_finder import time_to_clear
from radiological_material_clearance_finder.openmc_interop import from_depletion_results

results = openmc.deplete.Results("depletion_results.h5")
series = from_depletion_results(results, material_id=1)

seconds = time_to_clear(series, "StrlSchV_metal_recycling")
```

The keys are seconds since the start of the results, and the values are
`Material` objects. Nothing about the series is special, so a series assembled by
hand from any other depletion code works identically.

## Cooling time

The index falls roughly exponentially with cooling time, so `time_to_clear`
interpolates logarithmically in the index between the two samples that bracket
the threshold. Interpolating linearly instead places the crossing systematically
late, by a factor that grows with how far apart the samples are.

This package computes no decay of its own. That is deliberate: the series comes
from whatever produced it, so nothing here can disagree with the depletion the
numbers came from.

!!! warning "An index that rises again"

    A daughter can grow in faster than its parent decays, so a material can meet
    the limits at one time and fail them later. Reporting only the first crossing
    would be actively misleading, so `time_to_clear` raises `IngrowthError` when
    the index climbs back above the threshold at a later sample. Pass
    `allow_ingrowth=True` to take the first crossing anyway.

## What the agreement with OpenMC does and does not show

The US limit sets are checked against `openmc.Material.waste_disposal_rating` on
every CI run. The agreement covers the index, the per-nuclide breakdown and the
10 CFR 61.55 waste classification, and holds to floating point round-off:

```text
tests/test_agreement_openmc.py .......................  24 passed
```

Be clear about what that establishes. It shows that this package's path from an
inventory to a specific or volumetric activity, and its sum of fractions, match
an independent implementation exactly. It says nothing about whether the
underlying limit tables are right, because both implementations read the same
Fetter and NRC numbers. Table correctness is a separate question, addressed in
[where the data comes from](data-provenance.md), and it is checked against the
regulations rather than against OpenMC.

Two of the three high severity defects found in this package during review were
in the tables, not the arithmetic, and the arithmetic agreed with OpenMC exactly
throughout. Machine agreement and correctness are not the same thing.
