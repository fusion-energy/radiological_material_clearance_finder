# radiological material clearance finder

Clearance indexes for a nuclide inventory, against the UK, German, US, EU and
IAEA limit sets. Give it a dictionary of isotopes and atom numbers and it tells
you whether the material can be cleared, by which route, and what is driving the
answer.

No dependencies. OpenMC is not required, and is imported only if you ask for the
OpenMC adapter.

```python
from radiological_material_clearance_finder import Material, clearance_index

steel = Material({"Fe56": 8.4e22, "Co60": 2.1e9, "Cs137": 3.1e8})
result = clearance_index(steel, "UK_EPR16_out_of_scope")

result.index        # 11.244691988476553
result.clearable    # False
result.dominant(2)  # [('Co60', 11.215678908581781), ('Cs137', 0.029013079894772045)]
print(result)
```

```
UK_EPR16_out_of_scope: index 11.24 of 1  [NOT CLEARABLE]
  nuclide        activity        limit    share
  Co60              1.122          0.1   99.7%
  Cs137           0.02901            1    0.3%
```

## The index

Every regulation here uses the same arithmetic: divide each radionuclide's
activity by its tabulated limit and add up the ratios. Below one, the material
meets the limits. Germany calls it the Summenformel, the UK the summation rule,
the NRC the sum of fractions rule, and Fetter calls the result a waste disposal
rating. They are the same calculation over different tables.

What differs between regulations, and what the result therefore records, is
everything that is not a plain lookup:

```python
result.excluded    # daughters a parent's limit already accounts for
result.defaulted   # nuclides that took a catch-all limit
result.uncovered   # activity with no limit at all, absent from the index
result.uncovered_fraction   # how much of the inventory that is
result.out_of_scope         # the regulation excludes this material outright
```

`uncovered` matters. Under StrlSchV a nuclide with no tabulated limit simply
contributes nothing, so an index can look comfortable while a large activity
sits outside the sum. Checking `uncovered_fraction` is how you find out.

## Limit sets

```python
from radiological_material_clearance_finder import limit_sets, get_limit_set

limit_sets()                     # every set
limit_sets(jurisdiction="UK")    # just the UK ones
get_limit_set("StrlSchV_unrestricted").source
```

| Set | Basis |
| --- | --- |
| `UK_EPR16_out_of_scope` | EPR 2016 Sch 23 Table 2, the UK out of scope test for artificial radionuclides |
| `UK_EPR16_norm` | EPR 2016 Sch 23 Table 1, NORM industrial activities |
| `UK_EPR16_exempt_material` | EPR 2016 Sch 23 Table 5, keeping and using radioactive material |
| `UK_IRR17_notification` | IRR 2017 Sch 7 Part 1 col 2, exemption from notification |
| `UK_IRR17_registration` | IRR 2017 Sch 7 Part 1 col 4, exemption from registration |
| `UK_IRR17_natural` | IRR 2017 Sch 7 Part 2, unprocessed natural radionuclides |
| `StrlSchV_unrestricted` and 7 more pathways | StrlSchV 2018 Anlage 4 Tabelle 1 |
| `StrlSchV_exemption_activity` | StrlSchV Spalte 2, a total activity limit in Bq |
| `Fetter` | Fetter, Cheng and Mann (1990), class C disposal of activated metal |
| `NRC_long`, `NRC_short_A/B/C` | 10 CFR 61.55 Tables 1 and 2 |
| `EU_BSS_clearance` | 2013/59/Euratom Annex VII Table A Part 1 |
| `IAEA_GSR3_clearance` | IAEA GSR Part 3 Schedule I Table I.2 |

Every table is generated from the official source by a script in `tools/`, with
its URL and retrieval date recorded on the limit set. Nothing reaches the
network at runtime.

## Inputs

Specific activity in Bq/g is scale invariant, so atom counts, atom densities and
atom fractions all give the same answer and no density is needed. Volumetric
limit sets (the US ones, in Ci/m3) do need a density, which atom densities
supply on their own.

```python
Material({"Fe56": 8.4e22, "Co60": 2.1e9})            # atom counts
Material.from_atom_densities({"Fe56": 0.0849})       # atoms/barn-cm, OpenMC's unit
Material.from_specific_activities({"Co60": 1.1})     # Bq/g, as an assay reports
Material.from_masses({"Fe56": 1000.0})               # grams
Material.from_mass_fractions({"Fe56": 0.98})         # weight fractions
```

> **Include the stable isotopes.** The Bq/g denominator is the mass of
> everything you pass. An inventory filtered down to its radioactive nuclides
> gives a mass thousands of times too small and a specific activity thousands of
> times too high. The `Material` constructor warns when nothing you passed is
> stable.

## From OpenMC

```python
from radiological_material_clearance_finder.openmc_interop import from_openmc_material

material = from_openmc_material(openmc_material)
```

Ten lines, and OpenMC is imported inside the call. The rest of the package never
touches it.

## When does it clear?

Pass a series of materials at increasing cooling times, which is what a
depletion calculation already gives you, and the crossing time is interpolated
logarithmically in the index:

```python
from radiological_material_clearance_finder import time_to_clear

seconds = time_to_clear({0.0: mat_0, 1e7: mat_1, 1e8: mat_2}, "StrlSchV_metal_recycling")
```

This package computes no decay of its own, so nothing here can disagree with the
depletion code that produced the series. An index that rises with cooling time
raises `IngrowthError` rather than reporting a crossing the material later
reverses.

## Waste classification

```python
from radiological_material_clearance_finder import nrc_waste_class, uk_waste_category

nrc_waste_class(material, metal=True)   # 'Class A', 'Class B', 'Class C' or 'GTCC'
uk_waste_category(material)             # 'VLLW', 'LLW' or 'ILW', with the numbers
```

## Secular equilibrium

Where a limit already accounts for daughters, the regulations mark the parent
with `+` or `sec` and publish the daughter list themselves: StrlSchV Anlage 4
Tabelle 2, EPR 2016 Schedule 23 Table 3, IRR 2017 Schedule 7 Note 2. Those
tables are used directly, so no depletion chain is needed and the exclusion is
the regulation's own rather than an approximation of it.

The lists genuinely differ between regulations. EPR 2016 gives Zr-95 to Nb-95m
where IRR 2017 gives Zr-95 to Nb-95, so each limit set carries its own.

A daughter is only excluded when its parent is actually present, and where a
source lists a parent twice, plain and marked, the marked value applies only
when the daughters are really there. StrlSchV gives Th-232 as 10 Bq/g plain and
0.01 Bq/g marked, and which applies depends on the material.

## Regenerating the tables

```bash
python tools/build_uk_epr16.py       # legislation.gov.uk XML
python tools/build_uk_irr17.py
python tools/build_de_strlschv.py    # gesetze-im-internet.de
python tools/build_eu_bss.py         # EUR-Lex consolidated XHTML
python tools/build_iaea.py           # IAEA PDF, needs pdftotext
python tools/build_us.py             # extracted from OpenMC's waste.py
python tools/build_decay_data.py     # AME2020 and ENDF/B-VIII.0
python tools/build_decay_modes.py    # IAEA Livechart
python tools/verify_tables.py        # re-download everything and diff
```

Each source writes its numbers differently and none of them fail loudly when
misread. The UK XML marks exponents with `<Superior>`, so a flattened parse
turns 10² into 102; the German text uses `1 E-1` and decimal commas; EUR-Lex
uses comma decimals and non-breaking-space thousands separators. Every build
script asserts that each limit is one significant figure times a power of ten,
which is the shape these regulations use without exception, so a misread fails
rather than shipping.

## Tests

```bash
pytest                      # no OpenMC needed
pytest -m openmc            # cross-checks against OpenMC, skipped when absent
```

The US sets agree with `openmc.Material.waste_disposal_rating` to floating point
round-off. The German tables are diffed against the transcription in
openmc-dev/openmc#3898, and the IAEA PDF extraction against the EU XHTML parse,
which agree on all 257 nuclides in common.
