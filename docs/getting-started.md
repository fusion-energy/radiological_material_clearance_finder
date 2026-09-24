# Getting started

This page takes you from an empty environment to a clearance verdict you can
defend: how to install the package, how to describe your inventory with each of
the five `Material` constructors, which units need extra information and which
do not, and how to read every field of the result rather than just the number.

Every Python snippet on this page was run and its output pasted verbatim, so
none of the numbers below are illustrative.

## Installation

```bash
pip install radiological-material-clearance-finder
```

There are no runtime dependencies, so that is the whole of it. Python 3.10 or
newer. Wheels are published for Linux, macOS and Windows; the calculation is a
compiled Rust extension, so installing from source needs a Rust toolchain.

=== "For development"

    ```bash
    git clone https://github.com/fusion-energy/radiological_material_clearance_finder.git
    cd radiological_material_clearance_finder
    python -m pip install maturin
    maturin develop --release --extras test   # builds the Rust extension
    pytest
    cargo test                                # the Rust crate on its own
    ```

    Rebuild with `maturin develop --release` after changing any Rust, and
    after regenerating any table under
    `crates/radiological-material-clearance-finder/data/`, since the tables are
    compiled in.

=== "From Rust"

    The same library is a Rust crate, with no Python in it, for codes that want
    to assess their own material objects directly.

    ```toml
    [dependencies]
    radiological-material-clearance-finder = "0.1"
    ```

    See [using it from Rust](rust.md).

=== "With OpenMC"

    OpenMC is not required, and is imported only inside the two interop
    functions. Install it alongside if you want to build materials from an
    OpenMC model or a depletion result.

    ```bash
    pip install radiological-material-clearance-finder openmc
    ```

    See [using it with OpenMC](openmc.md).

There are **no runtime dependencies**. Python 3.10 or newer is enough: the wheel
is built against the stable ABI, so one wheel per platform covers every CPython
from 3.10 on, and the continuous integration matrix covers 3.10, 3.12 and 3.13.
OpenMC is not required:
it is imported lazily, inside the OpenMC adapter, and only if you call it.

Check the install:

```python
import radiological_material_clearance_finder as rmcf

print(rmcf.__version__)
```

```text
0.1.0
```

Note that the distribution name uses hyphens and the import name uses
underscores:

| | |
| --- | --- |
| distribution | `radiological-material-clearance-finder` |
| import | `radiological_material_clearance_finder` |

## The shortest complete example

Describe the material, name a limit set, print the result.

```python
from radiological_material_clearance_finder import Material, clearance_index

steel = Material(
    {
        "Fe56": 8.4e22,   # stable bulk: this is the mass the Bq/g divides by
        "Cr52": 1.6e22,
        "Ni58": 7.0e21,
        "Co60": 2.1e9,
        "Cs137": 3.1e8,
    },
    name="activated steel",
)

result = clearance_index(steel, "UK_EPR16_out_of_scope")
print(result)
```

```text
activated steel  UK_EPR16_out_of_scope: index 8.902 of 1  [NOT CLEARABLE]
  nuclide        activity        limit    share
  Co60             0.8879          0.1   99.7%
  Cs137           0.02297            1    0.3%
```

The index is a sum of activity-to-limit ratios. Below the threshold, normally 1,
the material meets that set's limits. Here it is nearly nine times over, and
Co-60 is 99.7% of the reason.

## Describing the material

There are five constructors. They all end up in the same place, a per-nuclide
specific activity in Bq/g, so pick whichever matches the numbers you actually
have rather than converting by hand.

| Constructor | Input unit | Derives a density? |
| --- | --- | --- |
| `Material` / `Material.from_atom_counts` | atoms | only with a `volume` |
| `Material.from_atom_densities` | atoms/barn-cm | yes, always |
| `Material.from_specific_activities` | Bq/g | no, pass `density=` |
| `Material.from_masses` | grams | no, pass `density=` |
| `Material.from_mass_fractions` | mass fraction or weight percent | no, pass `density=` |

### from_atom_counts

Atom counts are what an activation or depletion calculation produces. This is
also what the bare constructor takes, so `Material(...)` and
`Material.from_atom_counts(...)` are the same call.

```python
from radiological_material_clearance_finder import Material

steel = Material.from_atom_counts(
    {
        "Fe56": 8.4e22,
        "Cr52": 1.6e22,
        "Ni58": 7.0e21,
        "Co60": 2.1e9,
        "Cs137": 3.1e8,
    },
    name="activated steel",
)

print(steel)
print(f"mass              {steel.mass:.4g} g")
print(f"specific activity {steel.specific_activity():.6g} Bq/g")
print(f"total activity    {steel.activity('Bq'):.6g} Bq")
```

```text
Material('activated steel', 5 nuclides)
mass              9.856 g
specific activity 0.910855 Bq/g
total activity    8.97695 Bq
```

Absolute atom counts fix the total mass on their own, which is why the total
activity in Bq is available here with no volume and no density. Grams per
nuclide, below, are the other input that does this; atom densities and mass
fractions are relative and do not.

### from_atom_densities

Atom densities in atoms per barn-cm are OpenMC's unit. Because an atom density
carries both a composition and a number of atoms per unit volume, the mass
density follows from the atomic masses and is derived for you.

```python
from radiological_material_clearance_finder import Material

steel = Material.from_atom_densities(
    {
        "Fe56": 0.084,
        "Cr52": 0.016,
        "Ni58": 0.0070,
        "Co60": 2.1e-15,
        "Cs137": 3.1e-16,
    },
    name="activated steel",
)

print(f"derived density   {steel.density:.6g} g/cm3")
print(f"specific activity {steel.specific_activity():.6g} Bq/g")
print(f"volumetric        {steel.activity('Ci/m3'):.6g} Ci/m3")
```

```text
derived density   9.85552 g/cm3
specific activity 0.910855 Bq/g
volumetric        0.00024262 Ci/m3
```

Passing a `density=` as well is refused rather than silently ignored or
silently preferred, because it could only contradict the inventory:

```python
from radiological_material_clearance_finder import Material

try:
    Material.from_atom_densities({"Fe56": 0.084}, density=7.9)
except TypeError as exc:
    print(f"TypeError: {exc}")
```

```text
TypeError: from_atom_densities derives the mass density from the atom densities and the atomic masses, so passing density= could only contradict the inventory. Drop it, or use from_specific_activities if the density is the value you trust.
```

An atom density fixes a density but not a size, so the total mass is still
unknown until you add a `volume` in cm3:

```python
from radiological_material_clearance_finder import Material

steel = Material.from_atom_densities({"Fe56": 0.084, "Co60": 2.1e-15})
try:
    steel.mass
except Exception as exc:
    print(f"{type(exc).__name__}: {exc}")

sized = Material.from_atom_densities({"Fe56": 0.084, "Co60": 2.1e-15}, volume=1000.0)
print(f"{sized.mass:.6g} g")
```

```text
InsufficientDataError: total mass is unknown. This material holds relative amounts, so pass volume=, or build it with from_atom_counts or from_masses.
7802.1 g
```

### from_specific_activities

This is the assay route: someone measured Bq/g per nuclide and that is all you
have. No half-life or atomic mass data is used, because the specific activity is
already the direct input to the index.

```python
from radiological_material_clearance_finder import Material, clearance_index

assay = Material.from_specific_activities(
    {"Co60": 0.84, "Cs137": 0.016, "H3": 12.0},
    density=7.9,
    name="assay report",
)

print(f"specific activity {assay.specific_activity():.6g} Bq/g")
print(f"volumetric        {assay.activity('Ci/m3'):.6g} Ci/m3")
print()
print(clearance_index(assay, "StrlSchV_metal_recycling"))
```

```text
specific activity 12.856 Bq/g
volumetric        0.00274493 Ci/m3

assay report  StrlSchV_metal_recycling: index 1.439 of 1  [NOT CLEARABLE]
  nuclide        activity        limit    share
  Co60               0.84          0.6   97.3%
  Cs137             0.016          0.6    1.9%
  H3                   12         1000    0.8%
```

Activities alone do not imply a mass, so a `density=` in g/cm3 is needed for the
volumetric limit sets. It is the one constructor where you should pass the
density you trust rather than let one be derived.

### from_masses

Grams per nuclide, for when the inventory is a bill of materials rather than an
atom count.

```python
from radiological_material_clearance_finder import Material, clearance_index

offcut = Material.from_masses(
    {
        "Fe56": 700.0,      # grams
        "Cr52": 180.0,
        "Ni58": 120.0,
        "Co60": 2.0e-11,
        "Cs137": 5.0e-12,
    },
    density=7.9,
    name="1 kg offcut",
)

print(f"mass              {offcut.mass:.6g} g")
print(f"specific activity {offcut.specific_activity():.6g} Bq/g")
print(f"total activity    {offcut.activity('Bq'):.6g} Bq")
print()
print(clearance_index(offcut, "StrlSchV_metal_recycling"))
```

```text
mass              1000 g
specific activity 0.853448 Bq/g
total activity    853.448 Bq

1 kg offcut  StrlSchV_metal_recycling: index 1.422 of 1  [NOT CLEARABLE]
  nuclide        activity        limit    share
  Co60             0.8374          0.6   98.1%
  Cs137           0.01606          0.6    1.9%
```

### from_mass_fractions

Mass fractions, or weight percentages, which need not sum to one or to a hundred.
Useful when you know a composition and a contamination level but not a size.

```python
from radiological_material_clearance_finder import Material, clearance_index

alloy = Material.from_mass_fractions(
    {
        "Fe56": 70.0,     # weight percent, not normalised
        "Cr52": 18.0,
        "Ni58": 12.0,
        "Co60": 2.0e-12,
        "Cs137": 5.0e-13,
    },
    density=7.9,
    name="steel by weight percent",
)

print(f"specific activity {alloy.specific_activity():.6g} Bq/g")
print(clearance_index(alloy, "StrlSchV_metal_recycling"))
```

```text
specific activity 0.853448 Bq/g
steel by weight percent  StrlSchV_metal_recycling: index 1.422 of 1  [NOT CLEARABLE]
  nuclide        activity        limit    share
  Co60             0.8374          0.6   98.1%
  Cs137           0.01606          0.6    1.9%
```

That is the same index as the 1 kg offcut above, to every digit shown, because
the two describe the same material at one tenth of the scale. The two floats
actually differ in their last bit, which is the point of the next section.

## Units, honestly

### Bq/g is scale invariant

The specific activity divides an activity by a mass, and both scale together.
Double every number in the inventory and the Bq/g is unchanged. So atom counts,
atom densities and mass fractions all give the same answer for the Bq/g limit
sets, and none of them needs a density or a volume.

```python
from radiological_material_clearance_finder import Material

counts = Material({"Fe56": 8.4e22, "Cr52": 1.6e22, "Ni58": 7.0e21,
                   "Co60": 2.1e9, "Cs137": 3.1e8})
densities = Material.from_atom_densities({"Fe56": 0.084, "Cr52": 0.016, "Ni58": 0.0070,
                                          "Co60": 2.1e-15, "Cs137": 3.1e-16})

print(counts.specific_activity())
print(densities.specific_activity())
print(f"relative difference "
      f"{abs(counts.specific_activity() / densities.specific_activity() - 1):.1e}")
```

```text
0.9108552696589205
0.9108552696589204
relative difference 2.2e-16
```

The two differ in the last bit of the double, which is floating point round-off
from multiplying and dividing by different powers of ten, not a physical
difference. Do not test these for exact equality.

### When a density is needed

The five US sets (`Fetter`, `NRC_long`, `NRC_short_A`, `NRC_short_B`,
`NRC_short_C`) tabulate limits in Ci/m3. Converting Bq/g to Ci/m3 needs a mass
density, and there is no default: guessing one would quietly scale the verdict.

```python
from radiological_material_clearance_finder import Material, clearance_index

no_density = Material({"Fe56": 8.4e22, "Co60": 2.1e9})

# A Bq/g set is fine without one.
print(clearance_index(no_density, "StrlSchV_metal_recycling").index)

# A Ci/m3 set is not.
try:
    clearance_index(no_density, "Fetter")
except Exception as exc:
    print(f"{type(exc).__name__}: {exc}")
```

```text
1.8692798180969805
InsufficientDataError: mass density is unknown, and volumetric limits (Ci/m3) need it. Pass density= in g/cm3, or build the material with from_atom_densities, which derives it.
```

There are three ways to supply it, and they are checked against each other where
possible:

```python
from radiological_material_clearance_finder import Material

# 1. Derived from atom densities, which already carry it.
print(f"derived  {Material.from_atom_densities({'Fe56': 0.084, 'Co60': 2.1e-15}).density:.6g}")

# 2. Inferred from absolute atom counts plus a volume in cm3.
print(f"inferred {Material({'Fe56': 8.4e22, 'Co60': 2.1e9}, volume=1.0).density:.6g}")

# 3. Given outright.
print(f"given    {Material({'Fe56': 8.4e22}, density=7.9).density:.6g}")

# Over-determined and inconsistent is an error, not an average.
try:
    Material({"Fe56": 8.4e22}, density=7.9, volume=1.0).mass
except ValueError as exc:
    print(f"ValueError: {exc}")
```

```text
derived  7.8021
inferred 7.8021
given    7.9
ValueError: this material is over-determined and inconsistent: the atom amounts weigh 7.8021 g, but density times volume gives 7.9 g. Drop whichever of the two is not meant to describe this material.
```

The consistency check has a 1% tolerance and only fires when both a density and
a volume are given alongside absolute atom counts. In the third case above the
given 7.9 g/cm3 is simply used, and nothing cross-checks it against the 7.8021
g/cm3 the atom counts imply, because without a volume the atom counts fix no
density of their own.

### When a volume is needed

One set, `StrlSchV_exemption_activity`, is a total activity limit in Bq rather
than a concentration. That needs an absolute mass, so a material holding only
relative amounts needs a `volume` in cm3 before it can be assessed.

`clearance_indices` and `clearable_routes` skip any set whose units the material
cannot supply, rather than raising, so one missing number does not hide every
result that does not need it:

```python
from radiological_material_clearance_finder import Material, clearance_indices

# Atom densities: a density but no size, so the Bq set cannot be evaluated.
no_volume = Material.from_atom_densities({"Fe56": 0.084, "Co60": 2.1e-15})
print(len(clearance_indices(no_volume)), "of 22 sets")
print("StrlSchV_exemption_activity" in clearance_indices(no_volume))

# Mass fractions with no density: the five volumetric sets go too.
relative = Material.from_mass_fractions({"Fe56": 0.98, "Co60": 2.4e-14})
print(len(clearance_indices(relative)), "of 22 sets")
print(sorted(set(clearance_indices(no_volume)) - set(clearance_indices(relative))))
```

```text
21 of 22 sets
False
16 of 22 sets
['Fetter', 'NRC_long', 'NRC_short_A', 'NRC_short_B', 'NRC_short_C']
```

!!! tip "Skipping is deliberate, so count the results"
    A set that was skipped is simply absent from the returned dictionary. If you
    expected a US route and it is not there, the material has no density.

## Include the stable isotopes

!!! danger "This is the single most likely way to get a wrong answer"
    The Bq/g denominator is the mass of **everything** you pass, stable
    isotopes included. An inventory filtered down to just its radioactive
    nuclides has a mass that is orders of magnitude too small, so every
    specific activity, and the whole index, comes out orders of magnitude too
    high. In the example below the error is a factor of 3.5 x
    10<sup>13</sup>.

An activated steel is overwhelmingly stable by mass. Filtering the inventory to
the nuclides that "matter" removes the denominator:

```python
import warnings

from radiological_material_clearance_finder import Material, clearance_index

full = Material({"Fe56": 8.4e22, "Cr52": 1.6e22, "Ni58": 7.0e21,
                 "Co60": 2.1e9, "Cs137": 3.1e8})

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    filtered = Material({"Co60": 2.1e9, "Cs137": 3.1e8})
    for w in caught:
        print(f"{w.category.__name__}: {w.message}")

print()
print(f"full     {full.specific_activity():.6g} Bq/g, index "
      f"{clearance_index(full, 'UK_EPR16_out_of_scope').index:.6g}")
print(f"filtered {filtered.specific_activity():.6g} Bq/g, index "
      f"{clearance_index(filtered, 'UK_EPR16_out_of_scope').index:.6g}")
print(f"too high by a factor of "
      f"{filtered.specific_activity() / full.specific_activity():.4g}")
```

```text
UserWarning: every nuclide in this material is radioactive, so the total mass used for Bq/g is only the radioactive mass. If this inventory was filtered to its radioactive nuclides, add the stable ones back or the specific activity will be far too high.

full     0.910855 Bq/g, index 8.90184
filtered 3.21211e+13 Bq/g, index 3.13921e+14
too high by a factor of 3.526e+13
```

A factor of 3.5 x 10<sup>13</sup>. The warning is an ordinary `UserWarning`, so
in a plain script it appears on stderr and execution continues; the
`catch_warnings` block above is only there to make the output deterministic for
this page.

The check is deliberately narrow, so do not rely on it as a safety net. Three
things must all hold for it to fire: at least two nuclides are present, **every
one** of them is radioactive, and the material was built from atom amounts
rather than from specific activities. Miss any one and you get silence:

```python
from radiological_material_clearance_finder import Material

# One nuclide only: no warning, because a pure source is a legitimate input.
print(Material({"Co60": 2.1e9}).specific_activity())

# Built from Bq/g: no warning, because no mass was ever inferred.
print(Material.from_specific_activities({"Co60": 0.84}).specific_activity())
```

```text
41869403517512.44
0.84
```

The first number, 4.19 x 10<sup>13</sup> Bq/g, is the specific activity of pure
Co-60. That is the right answer for a pure source and badly wrong for a piece of
steel, and only you know which you meant.

## Nuclide names

Names are normalised on the way in, and spellings that mean the same nuclide are
added together rather than colliding.

```python
from radiological_material_clearance_finder import Material, normalise

for raw in ("Co60", "Co-60", "co60", "Ag108m", "Ag-108m", "Ag108_m1", "H-3", "60Co"):
    try:
        print(f"{raw:<10} -> {normalise(raw)}")
    except Exception as exc:
        print(f"{raw:<10} -> {type(exc).__name__}: {exc}")

print(Material.from_specific_activities({"Co60": 0.03, "Co-60": 0.02})
      .specific_activity(by_nuclide=True))
```

```text
Co60       -> Co60
Co-60      -> Co60
co60       -> Co60
Ag108m     -> Ag108_m1
Ag-108m    -> Ag108_m1
Ag108_m1   -> Ag108_m1
H-3        -> H3
60Co       -> NuclideNameError: cannot read '60Co' as a nuclide name
{'Co60': 0.05}
```

The leading-mass-number form `60Co` is **not** accepted. Negative amounts, NaN
amounts and a non-positive density are all rejected at construction:

```python
from radiological_material_clearance_finder import Material

for bad in ({"Co60": -1.0}, {"Co60": float("nan")}):
    try:
        Material.from_specific_activities(bad)
    except ValueError as exc:
        print(f"ValueError: {exc}")

try:
    Material.from_specific_activities({"Co60": 1.0}, density=0.0)
except ValueError as exc:
    print(f"ValueError: {exc}")
```

```text
ValueError: negative amount -1.0 for Co60
ValueError: amount for Co60 is NaN
ValueError: density must be greater than zero, got 0.0. A non-positive density would make every volumetric activity zero and so report the material as clearable.
```

A density of zero would make every volumetric activity zero, and so every Ci/m3
index zero, which reports as clearable. That is the most dangerous way for this
to fail, so it is rejected rather than allowed through.

## Reading a ClearanceResult

`clearance_index` returns a `ClearanceResult`. The bare ratio is `result.index`,
but the point of the object is everything a bare ratio throws away.

### The verdict

```python
from radiological_material_clearance_finder import Material, clearance_index

steel = Material({"Fe56": 8.4e22, "Cr52": 1.6e22, "Ni58": 7.0e21,
                  "Co60": 2.1e9, "Cs137": 3.1e8}, name="activated steel")
r = clearance_index(steel, "UK_EPR16_out_of_scope")

print("limit_set   ", r.limit_set)
print("index       ", r.index)
print("threshold   ", r.threshold)
print("units       ", r.units)
print("clearable   ", r.clearable)
print("dominant(2) ", r.dominant(2))
```

```text
limit_set    UK_EPR16_out_of_scope
index        8.901839343555865
threshold    1.0
units        Bq/g
clearable    False
dominant(2)  [('Co60', 8.878871193218828), ('Cs137', 0.022968150337037578)]
```

`clearable` is `index < threshold` or the material being out of scope, so read
it rather than comparing the index to 1 yourself.

Three parallel dictionaries carry the working: `activities` is each nuclide's
activity in the set's units, `limits_used` is the limit that was actually
applied after any adjustment, and `by_nuclide` is the ratio of the two, largest
first.

```python
from radiological_material_clearance_finder import Material, clearance_index

steel = Material({"Fe56": 8.4e22, "Co60": 2.1e9, "Cs137": 3.1e8})
r = clearance_index(steel, "StrlSchV_metal_recycling")

for name, ratio in r.by_nuclide.items():
    print(f"{name:<8} {r.activities[name]:>10.4g} Bq/g / {r.limits_used[name]:>6.4g} "
          f"= {ratio:.4g}")
```

```text
Co60          1.122 Bq/g /    0.6 = 1.869
Cs137       0.02901 Bq/g /    0.6 = 0.04836
```

Note that `activities` holds every nuclide present, including the stable ones at
zero activity, whereas `by_nuclide` and `limits_used` hold only the nuclides
that ended up in the sum.

### excluded and credited: secular equilibrium

Each regulation publishes its own parent-to-daughter table, so no decay chain
is inferred here and the lists are the regulation's rather than an approximation
of them. Where a parent's tabulated limit already accounts for a daughter, the
daughter drops out of the sum and the reason is recorded in `excluded`:

```python
from radiological_material_clearance_finder import Material, clearance_index

mat = Material.from_specific_activities({"Sr90": 0.5, "Y90": 0.5}, name="Sr90 + Y90")
r = clearance_index(mat, "UK_EPR16_out_of_scope")

print("index    ", r.index)
print("excluded ", r.excluded)
```

```text
index     0.5
excluded  {'Y90': 'in secular equilibrium with Sr90, whose limit already accounts for it'}
```

The credit is capped at the parent's own activity, because secular equilibrium
means equal activities. A daughter present in excess of its parent got there by
some other route, and the excess stays in the sum. That partial case appears in
`credited` rather than `excluded`:

```python
from radiological_material_clearance_finder import Material, clearance_index

mat = Material.from_specific_activities(
    {
        "Co60": 0.05,
        "Sr90": 0.5,
        "Y90": 2.0,     # four times its parent's activity
        "Ba133": 0.3,   # not listed in the EPR 2016 table at all
    },
    name="mixed",
)
r = clearance_index(mat, "UK_EPR16_out_of_scope")

print("excluded   ", r.excluded)
print("credited   ", r.credited)
print("defaulted  ", r.defaulted)
print("limits_used", r.limits_used)
print("by_nuclide ", r.by_nuclide)
```

```text
excluded    {}
credited    {'Y90': 0.5}
defaulted   ('Ba133',)
limits_used {'Co60': 0.1, 'Sr90': 1.0, 'Y90': 100.0, 'Ba133': 0.01}
by_nuclide  {'Ba133': 30.0, 'Co60': 0.5, 'Sr90': 0.5, 'Y90': 0.015}
```

Y-90's ratio is 0.015, which is (2.0 - 0.5) / 100 and not 2.0 / 100: the 0.5
Bq/g the parent could support was credited, and only the remaining 1.5 Bq/g was
assessed against Y-90's own limit. `activities` still reports the full 2.0 Bq/g,
so the two dictionaries genuinely differ and `credited` is how you reconcile
them.

Be clear about which part of this is quoted and which part is interpreted. The
daughter lists are the regulation's own, published in tables such as EPR 2016
Schedule 23 Table 3 and StrlSchV Anlage 4 Tabelle 2, and they genuinely differ
between regulations, so each set carries its own. The two conditions on the
credit are this package's reading of what secular equilibrium means: the parent
must itself have a limit in that set, otherwise the daughter would be dropped on
the strength of a limit that does not exist, and the credit is capped at the
parent's activity, otherwise a trace of Sr-90 would delete an arbitrarily large
Y-90 activity from the sum.

The whole mechanism can be switched off with `exclude_daughters=False`, which
double counts the daughter. That is conservative but not what the regulation
intends:

```python
from radiological_material_clearance_finder import Material, clearance_index

mat = Material.from_specific_activities({"Sr90": 0.5, "Y90": 0.5})
for flag in (True, False):
    r = clearance_index(mat, "UK_EPR16_out_of_scope", exclude_daughters=flag)
    print(f"exclude_daughters={flag!s:<5} index {r.index:.4g}  by_nuclide {r.by_nuclide}")
```

```text
exclude_daughters=True  index 0.5  by_nuclide {'Sr90': 0.5}
exclude_daughters=False index 0.505  by_nuclide {'Sr90': 0.5, 'Y90': 0.005}
```

### defaulted: the catch-all limit

Three sets define a catch-all limit for nuclides their table does not list: 0.01
Bq/g for `UK_EPR16_out_of_scope` and `UK_IRR17_notification`, and 0.1 Bq/g for
`UK_IRR17_registration`. Any nuclide that took it is named in `defaulted`. The
German, US, EU and IAEA sets define no such limit, which leads directly to the
most important field on the result.

### uncovered and uncovered_fraction: activity outside the sum

!!! warning "A low index can mean a small numerator or a missing limit"
    In a set with no catch-all, a nuclide the table does not list contributes
    **nothing** to the index. The sum stays comfortable while the activity is
    still there. `uncovered` names that activity and `uncovered_fraction` says
    how much of the inventory it is. Check it before trusting a low index.

Here is the same inventory under two regulations. Ti-44 (half-life 60 years) and
Ho-166m (1200 years) are both long lived, and neither has a row in the German
metal recycling table, StrlSchV 2018 Anlage 4 Tabelle 1 Spalte 14:

```python
from radiological_material_clearance_finder import Material, clearance_index

mat = Material.from_specific_activities(
    {"Co60": 0.3, "Ti44": 4.0, "Ho166_m1": 1.5},
    name="Ti bearing scrap",
)

german = clearance_index(mat, "StrlSchV_metal_recycling")
print(german)
print(f"uncovered          {german.uncovered}")
print(f"uncovered_fraction {german.uncovered_fraction:.4g}")
print(f"clearable          {german.clearable}")
print()

uk = clearance_index(mat, "UK_EPR16_out_of_scope")
print(uk)
print(f"defaulted          {uk.defaulted}")
print(f"uncovered_fraction {uk.uncovered_fraction:.4g}")
```

```text
Ti bearing scrap  StrlSchV_metal_recycling: index 0.5 of 1  [CLEARABLE]
  nuclide        activity        limit    share
  Co60                0.3          0.6  100.0%
  2 nuclide(s) with no limit, 94.83% of total activity
uncovered          {'Ti44': 4.0, 'Ho166_m1': 1.5}
uncovered_fraction 0.9483
clearable          True

Ti bearing scrap  UK_EPR16_out_of_scope: index 553 of 1  [NOT CLEARABLE]
  nuclide        activity        limit    share
  Ti44                  4         0.01   72.3%
  Ho166_m1            1.5         0.01   27.1%
  Co60                0.3          0.1    0.5%
defaulted          ('Ho166_m1', 'Ti44')
uncovered_fraction 0
```

Index 0.5 and `clearable` `True` under the German set, with 94.8% of the
activity outside the sum. The identical material is 553 times over the UK limit,
purely because the UK regulation applies 0.01 Bq/g to anything it does not list
and the German regulation applies nothing at all. `print(result)` shows the
uncovered count as a footer line, and `uncovered_fraction` is the number to
assert on in a script.

Turning the catch-all off shows the same effect from the other direction: the
UK index collapses from 553 to 3, and the activity reappears in `uncovered`.

```python
from radiological_material_clearance_finder import Material, clearance_index

mat = Material.from_specific_activities({"Co60": 0.3, "Ti44": 4.0, "Ho166_m1": 1.5})

on = clearance_index(mat, "UK_EPR16_out_of_scope")
off = clearance_index(mat, "UK_EPR16_out_of_scope", apply_default_limit=False)

print(f"catch-all on : index {on.index:>8.4g}  defaulted {on.defaulted}  "
      f"uncovered {on.uncovered_fraction:.1%}")
print(f"catch-all off: index {off.index:>8.4g}  defaulted {off.defaulted}  "
      f"uncovered {off.uncovered_fraction:.1%}")
```

```text
catch-all on : index      553  defaulted ('Ho166_m1', 'Ti44')  uncovered 0.0%
catch-all off: index        3  defaulted ()  uncovered 94.8%
```

`uncovered_fraction` is a share of total activity, not of the index, and it uses
the result's own units. It is 0.0 when every radionuclide present was either
limited or explicitly unlimited.

### unlimited: no limit on purpose

`uncovered` is activity the source is silent about. `unlimited` is activity the
source deliberately places no limit on, which is a different statement and is
kept separate. The Fetter set does this for 20 nuclides, tritium among them:

```python
from radiological_material_clearance_finder import Material, clearance_index

mat = Material.from_specific_activities(
    {"Co60": 1.0e3, "H3": 5.0e4, "Ni63": 2.0e4},
    density=7.9,
    name="Fetter check",
)
r = clearance_index(mat, "Fetter")

print(r)
print("unlimited ", r.unlimited)
print("uncovered ", r.uncovered)
print("activities", {k: round(v, 4) for k, v in r.activities.items()})
```

```text
Fetter check  Fetter: index 6.101e-06 of 1  [CLEARABLE]
  nuclide        activity        limit    share
  Ni63               4.27        7e+05  100.0%
  Co60             0.2135        3e+08    0.0%
unlimited  ('H3',)
uncovered  {}
activities {'H3': 10.6757, 'Ni63': 4.2703, 'Co60': 0.2135}
```

Tritium's 10.68 Ci/m3 is in `activities` and in neither the sum nor `uncovered`,
so `uncovered_fraction` stays at zero even though most of the activity is not in
the index. That is correct here and worth understanding before you use
`uncovered_fraction` as a completeness check.

### out_of_scope: a test on the whole material

The UK sets carry one whole-material scope test. As the limit set's own notes
record it, EPR 2016 Part 2 paragraph 7 places a substance outside the regulation
when **none** of its radionuclides has a half-life exceeding 100 seconds. That
is a test on the material, not a per-nuclide filter, so a short-lived nuclide
cannot be dropped from the sum on its own. When it applies, the verdict
overrides the index:

```python
from radiological_material_clearance_finder import Material, clearance_index, half_life

print("N16", half_life("N16"), "s")
print("O19", half_life("O19"), "s")
print("Al28", half_life("Al28"), "s")
print()

coolant = Material.from_specific_activities({"N16": 1e6, "O19": 5e5}, name="coolant")
r = clearance_index(coolant, "UK_EPR16_out_of_scope")
print(r)
print("index       ", r.index)
print("out_of_scope", r.out_of_scope)
print("clearable   ", r.clearable)
```

```text
N16 7.13 s
O19 26.88 s
Al28 134.484 s

coolant  UK_EPR16_out_of_scope: index 1.5e+08 of 1  [OUT OF SCOPE (every radionuclide is short lived)]
  nuclide        activity        limit    share
  N16               1e+06         0.01   66.7%
  O19               5e+05         0.01   33.3%
index        150000000.0
out_of_scope True
clearable    True
```

An index of 1.5 x 10<sup>8</sup> and `clearable` `True`, because the regulation
does not reach this material at all. Adding one nuclide with a half-life over
100 seconds brings the whole material into scope:

```python
from radiological_material_clearance_finder import Material, clearance_index

mat = Material.from_specific_activities({"N16": 1e6, "O19": 5e5, "Al28": 1.0})
r = clearance_index(mat, "UK_EPR16_out_of_scope")
print("out_of_scope", r.out_of_scope)
print("clearable   ", r.clearable)
print("index       ", r.index)
```

```text
out_of_scope False
clearable    False
index        150000100.0
```

### to_dict, for logging and serialising

```python
import json

from radiological_material_clearance_finder import Material, clearance_index

mat = Material.from_specific_activities({"Co60": 0.05, "Cs137": 0.02}, name="offcut")
r = clearance_index(mat, "StrlSchV_metal_recycling")
print(json.dumps(r.to_dict(), indent=2))
```

```text
{
  "limit_set": "StrlSchV_metal_recycling",
  "material": "offcut",
  "index": 0.11666666666666667,
  "threshold": 1.0,
  "units": "Bq/g",
  "clearable": true,
  "out_of_scope": false,
  "by_nuclide": {
    "Co60": 0.08333333333333334,
    "Cs137": 0.03333333333333333
  },
  "activities": {
    "Co60": 0.05,
    "Cs137": 0.02
  },
  "limits_used": {
    "Co60": 0.6,
    "Cs137": 0.6
  },
  "defaulted": [],
  "excluded": {},
  "credited": {},
  "uncovered": {},
  "unlimited": [],
  "uncovered_fraction": 0.0
}
```

`to_dict` flattens the two properties, `clearable` and `uncovered_fraction`,
into the mapping so a stored record does not need the class to be interpreted.

## Comparing every route at once

`clearance_indices` assesses many sets and `clearable_routes` returns just the
ones that pass, best margin first. Print the uncovered fraction alongside the
index, or the narrow tables will flatter you:

```python
from radiological_material_clearance_finder import Material, clearance_indices

steel = Material.from_atom_densities(
    {"Fe56": 0.084, "Cr52": 0.016, "Ni58": 0.0070, "Co60": 2.1e-15, "Cs137": 3.1e-16},
    volume=1000.0,
    name="activated steel",
)

print(f"{'limit set':<30} {'index':>10}  {'uncovered':>9}  clearable")
for name, r in sorted(clearance_indices(steel).items(),
                      key=lambda kv: kv[1].index / kv[1].threshold):
    print(f"{name:<30} {r.index:>10.4g}  {r.uncovered_fraction:>8.1%}  {r.clearable}")
```

```text
limit set                           index  uncovered  clearable
NRC_long                                0    100.0%  True
UK_EPR16_norm                           0    100.0%  True
UK_IRR17_natural                        0    100.0%  True
Fetter                          1.231e-10      0.0%  True
NRC_short_C                      1.33e-09     97.5%  True
NRC_short_B                      1.39e-07     97.5%  True
NRC_short_A                     6.456e-06      0.0%  True
UK_EPR16_exempt_material          0.09109      0.0%  True
UK_IRR17_registration             0.09109      0.0%  True
StrlSchV_exemption_activity        0.1101      0.0%  True
StrlSchV_incineration_100          0.1291      0.0%  True
StrlSchV_landfill_100              0.1503      0.0%  True
StrlSchV_landfill_1000             0.4468      0.0%  True
StrlSchV_incineration_1000         0.4516      0.0%  True
StrlSchV_metal_recycling            1.518      0.0%  False
UK_EPR16_out_of_scope               8.902      0.0%  False
EU_BSS_clearance                    9.109      0.0%  False
IAEA_GSR3_clearance                 9.109      0.0%  False
StrlSchV_unrestricted               9.109      0.0%  False
UK_IRR17_notification               9.109      0.0%  False
StrlSchV_rubble                     9.923      0.0%  False
StrlSchV_soil                       29.98      0.0%  False
```

The top three rows are index 0 with 100% of the activity uncovered. `NRC_long`
tabulates exactly three nuclides, C-14, I-129 and Tc-99; `UK_EPR16_norm` and
`UK_IRR17_natural` cover only naturally occurring series members such as Ra-226
and the thorium and uranium chains. None of them lists Co-60 or Cs-137, so
neither activity has anywhere to go. Those three are not clearance routes for
this material; they are tables that have nothing to say about it. Read the
uncovered column, always.

## Where the numbers came from

Every limit set records its own provenance, so you never have to take a limit on
trust from this documentation:

```python
from radiological_material_clearance_finder import get_limit_set

s = get_limit_set("UK_EPR16_out_of_scope")
print("jurisdiction  ", s.jurisdiction)
print("units         ", s.units)
print("rows          ", len(s.limits))
print("default_limit ", s.default_limit)
print("source        ", s.source)
print("url           ", s.url)
print("retrieved     ", s.retrieved)
```

```text
jurisdiction   UK
units          Bq/g
rows          282
default_limit  0.01
source         Environmental Permitting (England and Wales) Regulations 2016, Schedule 23, Part 3 Table 2
url            https://www.legislation.gov.uk/ukdsi/2016/9780111150184/schedule/23
retrieved      2026-09-10
```

There is also a `notes` field, which is where a set explains its scope test, its
catch-all and any reason it is not expected to agree with a neighbouring set.

## What is verified and what is not

Being straight about the provenance of the numbers on this page:

- Every code block above was executed and its output pasted verbatim, on
  CPython 3.14.4 with version 0.1.0 of the package, and the full test suite of
  277 tests passes.
- The limit tables are generated from official sources by scripts in `tools/`,
  and each set records the URL and retrieval date it came from. Nothing reaches
  the network at run time.
- The US sets agree with `openmc.Material.waste_disposal_rating` to floating
  point round-off. The IAEA PDF extraction and the EU XHTML parse agree on all
  257 values they have in common, by independent paths. All 81 Fetter lower
  bounds match the 1990 paper. Four regulatory tables re-derived from the live
  sources match the committed data.
- **No human has checked any of it against the source regulations.** All of the
  agreement above is machine against machine. Treat this package as a
  calculator, not as an authority, and check anything that matters against the
  regulation itself.

## Next steps

- [Limit sets](api/limits.md) for the 22 sets, their jurisdictions, units and
  recorded sources.
- [Material](api/material.md) for every constructor argument, the activity unit
  conversions and the exact error conditions.
- [Clearance index](api/clearance-index.md) for the full `ClearanceResult`
  reference and the `metal`, `apply_default_limit` and `exclude_daughters`
  keywords.
- [Cooling](api/cooling.md) for `time_to_clear`, which interpolates a crossing
  time from a series of materials at increasing cooling times.
- [OpenMC interop](api/openmc_interop.md) for starting from an
  `openmc.Material` or a depletion result.
