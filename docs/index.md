# radiological material clearance finder

Give this package a nuclide inventory and it tells you whether the material
meets the clearance, exemption or waste disposal limits of the UK, German, US,
EU and IAEA regulations, by which route, and what is driving the answer.

It exists because the arithmetic is trivial and the tables are not. An
activation calculation hands you a few hundred nuclides; the regulation hands
you a table with its own units, its own catch-all rules, its own list of
daughters a parent already accounts for, and its own idea of what is in scope at
all. This package implements the arithmetic once and treats all 22 tables
(4706 limit values) as data, generated from the official sources.

## Install

Not on PyPI yet, so install from a checkout:

```bash
pip install radiological-material-clearance-finder
```

There is nothing else to install: in a fresh virtual environment, `pip freeze`
afterwards lists exactly one package, this one.

## The shortest useful example

A gram-scale piece of activated steel, as atom counts. Stable Fe-56 is in the
inventory because the mass of everything you pass is the Bq/g denominator.

```python
from radiological_material_clearance_finder import Material, clearance_index

steel = Material({"Fe56": 8.4e22, "Co60": 2.1e9, "Cs137": 3.1e8})
result = clearance_index(steel, "UK_EPR16_out_of_scope")

print(result.index)
print(result.clearable)
print(result)
```

```
11.244691988476553
False
UK_EPR16_out_of_scope: index 11.24 of 1  [NOT CLEARABLE]
  nuclide        activity        limit    share
  Co60              1.122          0.1   99.7%
  Cs137           0.02901            1    0.3%
```

Eleven times over the limit, and Co-60 is 99.7% of the reason. To ask which of
the shipped routes this material would already pass:

```python
from radiological_material_clearance_finder import clearable_routes

print(clearable_routes(steel))
```

```
['UK_EPR16_norm', 'UK_IRR17_natural', 'StrlSchV_exemption_activity', 'UK_EPR16_exempt_material', 'UK_IRR17_registration', 'StrlSchV_incineration_100', 'StrlSchV_landfill_100', 'StrlSchV_landfill_1000', 'StrlSchV_incineration_1000']
```

Two of those are not the good news they look like. See
[what it could not account for](#the-result-reports-what-it-could-not-account-for)
below.

## The sum of fractions

Every regulation here uses the same calculation: divide each radionuclide's
activity by its tabulated limit and add up the ratios. Below one, the material
meets the limits. Germany calls it the Summenformel, the UK the summation rule,
the NRC the sum of fractions rule, and Fetter calls the result a waste disposal
rating. They are one calculation over different tables, so it is implemented
once.

Specific activity in Bq/g is scale invariant, which is why the example above
needed no density and no volume: atom counts, atom densities in atoms/barn-cm
and mass fractions all give the same index, to floating point round-off. The
volumetric sets (the US ones, in Ci/m3) do need a density, and atom densities
supply one on their own.

## What is different about this one

### Zero dependencies

`dependencies = []`, and importing the package pulls in nothing outside the
Python standard library and its own compiled extension. OpenMC is not required
and is imported only inside the
two functions of [`openmc_interop`](api/openmc_interop.md), if you call them.
That is deliberate: a clearance check should be runnable anywhere, including
somewhere that cannot install a neutronics stack.

### A Rust core that transport codes can link

The calculation is a Rust crate, and the Python package is a thin binding over
it. A transport or transmutation code, such as yamc or yani, can depend on the
crate and assess its own material objects without going through Python. See
[using it from Rust](rust.md).

### The tables are data, with provenance

Adding a jurisdiction means adding a JSON file and a build script, never a
branch in the calculation. Every set records where its numbers came from and
when they were fetched, and you can ask it at runtime:

```python
from radiological_material_clearance_finder import get_limit_set

limits = get_limit_set("StrlSchV_metal_recycling")
print(limits.label)
print(limits.source)
print(limits.url)
print(limits.retrieved)
print(limits.limits["Co60"], limits.units)
```

```
Metal scrap for recycling
Strahlenschutzverordnung (StrlSchV) 2018, Anlage 4 Tabelle 1, Spalte 14
https://www.gesetze-im-internet.de/strlschv_2018/anlage_4.html
2026-09-10
0.6 Bq/g
```

Nothing reaches the network at run time. The scripts in `tools/` regenerate the
tables from the official sources on demand, and each asserts that every limit is
one significant figure times a power of ten, which is the shape these
regulations use without exception, so a misparsed exponent fails the build
instead of shipping quietly.

### The result reports what it could not account for

A bare number cannot tell you that it left something out. `ClearanceResult`
can, and this is the part worth knowing before you trust an index. Recall
`UK_EPR16_norm` from the passing list above:

```python
norm = clearance_index(steel, "UK_EPR16_norm")
print(norm)
print(norm.uncovered_fraction)
```

```
UK_EPR16_norm: index 0 of 1  [CLEARABLE]
  2 nuclide(s) with no limit, 100.00% of total activity
1.0
```

An index of zero, formally clearable, and every becquerel in the material
outside the sum. That table is the NORM one and covers only natural-series
radionuclides, so neither Co-60 nor Cs-137 appears in it, and this set has no
catch-all limit, so the two of them contribute nothing at all.
`uncovered_fraction` of 1.0 is how you find that out. A result carries five
such fields:

| Field | What it records |
| --- | --- |
| `uncovered` | activity with no limit at all, absent from the index |
| `defaulted` | nuclides that fell back on the set's catch-all limit |
| `excluded` | daughters a parent's limit already accounts for in full |
| `credited` | activity a parent could only partly account for, the rest still in the sum |
| `out_of_scope` | the regulation excludes this material outright |

The last three exist because the regulations publish their own
parent-to-daughter equilibrium tables, and those tables genuinely disagree: EPR
2016 pairs Zr-95 with Nb-95m where IRR 2017 pairs it with Nb-95. Each limit set
therefore carries its own, so no decay chain is assumed on your behalf.

!!! warning "Include the stable isotopes"

    The Bq/g denominator is the mass of everything you pass. Drop the stable
    Fe-56 from the example above and the same two radionuclides give an index
    of `3.14e+14` instead of `11.24`, because the inventory now weighs
    2.8e-13 g rather than 7.8 g. An inventory filtered to its radioactive
    nuclides is the single easiest way to get a wrong answer here. `Material`
    warns when nothing you passed is stable.

## Verification, honestly

The US sets agree with `openmc.Material.waste_disposal_rating` to floating point
round-off, all 81 Fetter lower bounds match the 1990 paper, the IAEA PDF
extraction and the EU XHTML parse agree on all 257 values they have in common by
independent routes, and four regulatory tables re-derived from the live sources
match the committed data. The suite currently stands at 277 tests, all passing.

All of that is machine-to-machine agreement. No human has checked any of these
tables against the published regulations, and this package is not a substitute
for reading them or for your own regulatory advice.

## Where to go next

- [API reference](api/index.md), the package overview and every public symbol.
- [Material](api/material.md) for the five ways to build an inventory: atom
  counts, atom densities, specific activities, masses and mass fractions.
- [Limit sets](api/limits.md) for the 22 shipped sets, their provenance, and
  registering your own.
- [Cooling](api/cooling.md) for `time_to_clear`, which interpolates a crossing
  time from a series of cooled materials and refuses to report one that ingrowth
  later reverses.
- [Classification](api/classify.md) for NRC waste classes and UK waste
  categories.
- [OpenMC interop](api/openmc_interop.md) for starting from an
  `openmc.Material` or a depletion results file.
