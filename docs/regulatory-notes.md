# Regulatory notes

Every limit set in this package is evaluated with the same arithmetic: divide
each radionuclide's activity by its tabulated limit, sum the ratios, and compare
the total against one. Germany calls that sum the Summenformel, the UK the
summation rule, the NRC the sum of fractions rule, and Fetter calls the result a
waste disposal rating.

The arithmetic is not where wrong answers come from. They come from the rules
around it: which nuclides are in the sum at all, which limit a nuclide takes
when the table lists it twice, what happens to a nuclide the table never
mentions, and whether the regulation applies to the material in the first place.
This page documents those rules, set by set, with the real numbers.

Every example below was run against the installed package and the output is
copied verbatim.

## Secular equilibrium

### No depletion chain is needed

A regulation that limits a parent nuclide usually says its value already
accounts for the daughters that grow in with it. The German regulation marks
such a parent with a `+` and lists the daughters in Anlage 4 Tabelle 2; the UK
and the EU and the IAEA all publish equivalent parent-and-progeny tables.

Because the regulation publishes the list, this package never has to work out
which daughters a parent has. It reads the list out of the regulation and
applies it. Nothing at runtime consults a decay chain, branching ratios or
half-life ordering for this purpose. The only nuclear data shipped for the
calculation is half-life, atomic mass and alpha fraction, and half-life is used
solely for the NRC five year rule and the UK scope test described below.

!!! note "Why this matters for accuracy"
    A physically derived chain would be a different answer, not a better one.
    Where a regulation's list is incomplete or idiosyncratic, the regulation's
    list is still the one the material is judged against. Substituting a
    correct-looking chain would produce a number no regulator would recognise.

### The lists genuinely differ between regulations

This is the part that surprises people. Two regulations covering the same
nuclide can name different daughters for it. Zirconium-95 is the standard
example: EPR 2016 Schedule 23 Part 3 Table 3 gives it Nb-95m, while IRR 2017
Schedule 7 Part 1 gives it Nb-95.

```python
from radiological_material_clearance_finder import Material, clearance_index

crud = Material.from_specific_activities({"Zr95": 0.5, "Nb95": 0.5}, name="clad crud")

for name in ("UK_IRR17_notification", "UK_EPR16_out_of_scope"):
    r = clearance_index(crud, name)
    print(name)
    print(f"  index       {r.index:.4g}")
    print(f"  limits used {r.limits_used}")
    print(f"  excluded    {r.excluded}")
```

```text
UK_IRR17_notification
  index       0.5
  limits used {'Zr95': 1.0}
  excluded    {'Nb95': 'in secular equilibrium with Zr95, whose limit already accounts for it'}
UK_EPR16_out_of_scope
  index       5.5
  limits used {'Zr95': 0.1, 'Nb95': 1.0}
  excluded    {}
```

Under IRR17 the Nb-95 leaves the sum entirely, because Zr-95's limit is declared
to account for it. Under EPR16 it stays in the sum on its own row, because that
regulation's Table 3 names Nb-95m rather than Nb-95, and no Nb-95m is present.

!!! warning "The two indexes above are not a like for like comparison"
    These two sets differ in their Zr-95 limit as well as in their daughter
    list, 1.0 Bq/g against 0.1 Bq/g, so the gap between 0.5 and 5.5 is caused by
    both differences together. The daughter list alone determines the
    *membership* of the sum, which is what `excluded` reports.

The disagreement is not confined to Zr-95. Of the 51 parents that EPR16
Schedule 23 Part 3 Table 3 and IRR 2017 Schedule 7 Part 1 have in common, 12
carry different daughters: Am-242m, Bi-212, Pb-212, Ra-224, Ra-226, Th-228,
Th-229, Th-234, U-232, U-235, U-238 and Zr-95.

It also happens *within* a single regulation. EPR 2016 Part 6 paragraph 29
directs the markers in its own Table 5 to Part 6 Table 8, not to the Part 3
Table 3 used by the out of scope values. Table 3 has 57 parents and Table 8 has
31; of the 27 they share, 16 carry different daughters. Table 3 gives U-235
twelve daughters (Th-231, Pa-231, Ac-227, Th-227, Fr-223, Ra-223, Rn-219,
Po-215, Pb-211, Bi-211, Tl-207, Po-211) where Table 8 gives it only Th-231.

### Which set uses which daughter table

Each limit set carries its own table, so no set can borrow another's. The 22
shipped sets resolve to six distinct tables:

| Daughter table | Parents | Limit sets |
| --- | --- | --- |
| StrlSchV Anlage 4 Tabelle 2 | 187 | the nine `StrlSchV_*` sets |
| EPR 2016 Sch 23 Part 3 Table 3 | 57 | `UK_EPR16_out_of_scope`, `UK_EPR16_norm` |
| IRR 2017 Sch 7 Part 1 progeny table | 56 | `UK_IRR17_notification`, `UK_IRR17_registration` |
| EU 2013/59 Annex VII Table A progeny footnote | 36 | `EU_BSS_clearance` |
| IAEA GSR Part 3 Table I.2 progeny footnote | 36 | `IAEA_GSR3_clearance` |
| EPR 2016 Sch 23 Part 6 Table 8 | 31 | `UK_EPR16_exempt_material` |
| IRR 2017 Sch 7 Part 2 progeny table | 4 | `UK_IRR17_natural` |
| none | 0 | `Fetter`, `NRC_long`, `NRC_short_A`, `NRC_short_B`, `NRC_short_C` |

The EU and IAEA tables are separate sources that happen to be identical, which
is the expected result since the directive adopts the IAEA values. You can
regenerate this grouping at any time:

```python
from radiological_material_clearance_finder import get_limit_set, limit_sets

groups = {}
for name in limit_sets():
    eq = get_limit_set(name).secular_equilibrium
    key = tuple(sorted((p, eq[p]) for p in eq))
    groups.setdefault(key, []).append(name)

for key, names in sorted(groups.items(), key=lambda kv: -len(kv[0])):
    print(f"{len(key):3d} parents  {', '.join(names)}")
```

```text
187 parents  StrlSchV_exemption_activity, StrlSchV_incineration_100, StrlSchV_incineration_1000, StrlSchV_landfill_100, StrlSchV_landfill_1000, StrlSchV_metal_recycling, StrlSchV_rubble, StrlSchV_soil, StrlSchV_unrestricted
 57 parents  UK_EPR16_norm, UK_EPR16_out_of_scope
 56 parents  UK_IRR17_notification, UK_IRR17_registration
 36 parents  EU_BSS_clearance, IAEA_GSR3_clearance
 31 parents  UK_EPR16_exempt_material
  4 parents  UK_IRR17_natural
  0 parents  Fetter, NRC_long, NRC_short_A, NRC_short_B, NRC_short_C
```

The US sets have no daughter table because 10 CFR 61.55 and Fetter do not
publish one. Nothing is credited to a parent in those sets, and no
`exclude_daughters` setting changes that.

### Two conditions before a parent can account for a daughter

The regulation's claim that a parent's value "already takes into account the
daughter radionuclides present" rests on two things being true. Both are checked
rather than assumed.

**The parent must itself be limited by this set.** A parent with no row in the
column being used contributes nothing to the sum, so crediting its daughters
against it would remove them on the strength of a limit that does not exist,
pushing the index towards zero. Nine of the shipped sets contain at least one
parent that appears in the daughter table but has no limit of its own.

The StrlSchV soil column is the sharpest case: Fe-60 heads the Tabelle 2 entry
for Co-60, but the soil column gives Fe-60 no value.

```python
from radiological_material_clearance_finder import Material, clearance_index

# Fe-60 heads the equilibrium table entry for Co-60 in StrlSchV Anlage 4
# Tabelle 2, but the soil column gives Fe-60 no limit of its own.
soil = Material.from_specific_activities({"Fe60": 0.001, "Co60": 0.02}, name="spoil")

r = clearance_index(soil, "StrlSchV_soil")
print("index      ", f"{r.index:.4g}")
print("limits used", r.limits_used)
print("excluded   ", r.excluded)
print("uncovered  ", r.uncovered)
print("uncovered_fraction", f"{r.uncovered_fraction:.4%}")
```

```text
index       0.6667
limits used {'Co60': 0.03}
excluded    {}
uncovered   {'Fe60': 0.001}
uncovered_fraction 4.7619%
```

The Co-60 stays in the sum against its own 0.03 Bq/g limit. Had the credit been
granted, a 1 mBq/g trace of Fe-60 would have deleted the entire Co-60
contribution and the index would have read zero. Other parents in the same
column with the same problem include Np-237 (whose daughter Pa-233 has a limit
of 0.4 Bq/g), Si-32 (P-32 at 0.02 Bq/g), Hf-182 (Ta-182 at 0.06 Bq/g), Pb-202
(Tl-202 at 0.2 Bq/g) and Ac-227 (Ra-223 at 0.01 Bq/g).

**The credit is bounded by the parent's own activity.** Secular equilibrium
means the daughter's activity equals its parent's, so that is all the parent can
account for. A daughter present in excess of its parent got there by some other
route, and the excess stays in the sum.

### The `credited` field

When a parent can account for only part of a daughter's activity, the daughter
is not excluded. The supported part is subtracted, the remainder is assessed
against the daughter's own limit, and the subtracted amount is reported in
`credited`.

```python
from radiological_material_clearance_finder import Material, clearance_index

# Secular equilibrium would put Y-90 at the same activity as its Sr-90 parent.
# Here there is three times more, so two thirds of it arrived by another route.
mat = Material.from_specific_activities({"Sr90": 1.0, "Y90": 3.0}, name="separated yttrium")

r = clearance_index(mat, "EU_BSS_clearance")
print("index       ", f"{r.index:.6g}")
print("activities  ", r.activities)
print("credited    ", r.credited)
print("excluded    ", r.excluded)
print("by_nuclide  ", {k: round(v, 6) for k, v in r.by_nuclide.items()})
```

```text
index        1.002
activities   {'Y90': 3.0, 'Sr90': 1.0}
credited     {'Y90': 1.0}
excluded     {}
by_nuclide   {'Sr90': 1.0, 'Y90': 0.002}
```

One Bq/g of the Y-90 is credited to the Sr-90. The remaining 2 Bq/g is divided by
Y-90's own limit of 1000 Bq/g, contributing 0.002 and taking the index from 1.000
to 1.002, which is the difference between clearable and not. Without the bound, a
trace of Sr-90 would have deleted an arbitrarily large Y-90 activity from the sum.

So `excluded` and `credited` are mutually exclusive per nuclide: a daughter is
fully accounted for and drops out, or it is partly accounted for and stays with a
reduced activity.

### The `exclude_daughters` switch is not the conservative direction

`clearance_index(..., exclude_daughters=False)` turns the whole mechanism off.
It is documented as double counting, which sounds conservative, but for a parent
whose limit switches on the presence of its daughters it is not.

```python
from radiological_material_clearance_finder import Material, clearance_index

aged = Material.from_specific_activities(
    {"Th232": 1.0, "Ra228": 1.0, "Ac228": 1.0, "Th228": 1.0, "Ra224": 1.0,
     "Rn220": 1.0, "Po216": 1.0, "Pb212": 1.0, "Bi212": 1.0, "Tl208": 0.36},
    name="aged thorium",
)

for flag in (True, False):
    r = clearance_index(aged, "StrlSchV_unrestricted", exclude_daughters=flag)
    print(f"exclude_daughters={flag!s:<6} index {r.index:>8.4g}  "
          f"Th-232 limit {r.limits_used['Th232']:>6g}  {len(r.by_nuclide)} nuclide(s) in the sum")
```

```text
exclude_daughters=True   index      100  Th-232 limit   0.01  1 nuclide(s) in the sum
exclude_daughters=False  index     20.5  Th-232 limit     10  8 nuclide(s) in the sum
```

Turning the switch off makes the index five times *lower*, for two reasons that
compound. The marked Th-232 limit of 0.01 Bq/g is only selected on the same code
path that performs the exclusion, so switching the path off reverts Th-232 to its
plain 10 Bq/g row. And two of the ten nuclides, Po-216 and Tl-208, have no row of
their own in the German column, so with no parent to account for them they fall
out of the sum into `uncovered` instead.

!!! warning "Leave `exclude_daughters` at its default"
    The default of `True` is what the regulations say. Use `False` only to
    inspect what the individual rows contribute, and read `uncovered` and
    `limits_used` when you do, because the number it produces can be either
    higher or lower than the regulatory one.

## Nuclides listed twice, plain and marked

Some nuclides have two rows in the same column: one plain, and one marked to
indicate the daughters are present. Which value applies depends on the material,
not on the table, so both are kept. The plain value is the limit in
`LimitSet.limits`, and the marked value sits in
`LimitSet.limits_secular_equilibrium`, applied only when at least one of that
parent's tabulated daughters is actually in the inventory.

### Th-232 in StrlSchV, a factor of 1000

Anlage 4 Tabelle 1 Spalte 3 gives Th-232 as 10 Bq/g plain and 0.01 Bq/g marked.

```python
from radiological_material_clearance_finder import Material, clearance_index

bare = Material.from_specific_activities({"Th232": 1.0}, name="fresh thorium")
aged = Material.from_specific_activities(
    {"Th232": 1.0, "Ra228": 1.0, "Ac228": 1.0, "Th228": 1.0, "Ra224": 1.0,
     "Rn220": 1.0, "Po216": 1.0, "Pb212": 1.0, "Bi212": 1.0, "Tl208": 0.36},
    name="aged thorium",
)

for mat in (bare, aged):
    r = clearance_index(mat, "StrlSchV_unrestricted")
    print(f"{mat.name:16s} Th-232 limit {r.limits_used['Th232']:>8g} Bq/g   "
          f"index {r.index:>8.4g}   in sum {sorted(r.by_nuclide)}")
```

```text
fresh thorium    Th-232 limit       10 Bq/g   index      0.1   in sum ['Th232']
aged thorium     Th-232 limit     0.01 Bq/g   index      100   in sum ['Th232']
```

The same 1 Bq/g of Th-232 clears comfortably on its own and fails by a factor of
100 once its chain is present, and the whole ten nuclide chain collapses into a
single Th-232 row. Zr-97 and Pb-212 are also listed twice in StrlSchV; Zr-97
carries the same value either way, and in the exemption activity set Pb-212 is
10<sup>7</sup> Bq plain against 10<sup>5</sup> Bq marked.

### U-240 in IRR 2017, a factor of 10 000 the other way

Do not assume the marked value is the stricter one. In IRR 2017 Schedule 7
Part 1 column 2 the marked U-240 value is four orders of magnitude *more*
lenient than the plain one, 100 Bq/g against 0.01 Bq/g.

```python
from radiological_material_clearance_finder import Material, clearance_index

alone = Material.from_specific_activities({"U240": 1.0}, name="U-240 alone")
with_np = Material.from_specific_activities({"U240": 1.0, "Np240": 1.0}, name="U-240 + Np-240")

for mat in (alone, with_np):
    r = clearance_index(mat, "UK_IRR17_notification")
    print(f"{mat.name:16s} U-240 limit {r.limits_used['U240']:>8g} Bq/g   "
          f"index {r.index:>8.4g}   excluded {sorted(r.excluded)}")
```

```text
U-240 alone      U-240 limit     0.01 Bq/g   index      100   excluded []
U-240 + Np-240   U-240 limit      100 Bq/g   index     0.01   excluded ['Np240']
```

This is why the two rows are kept apart rather than one overwriting the other.
Letting the marked row win unconditionally would make U-240 on its own ten
thousand times too lenient. The direction reverses between regulations, too:
EPR 2016 Table 5 and IRR 2017 column 4 both give U-240 1000 Bq/g plain against
10 Bq/g marked, the stricter direction, while the IRR 2017 notification column
above runs the other way.

### Whole-chain rows are carried but not selected

A third kind of row exists. Some sources give a parent a value for the parent
taken *with its entire natural series* in secular equilibrium, rather than with
a listed set of daughters. Those are stored under a `_sec` suffix so they stay
distinct from a plain row.

!!! warning "Known limitation"
    Unlike the marked `+` rows, a `_sec` row is **not** selected automatically
    when the chain is present. It is promoted to the plain limit only when it is
    the only row published for that nuclide, which is what makes
    `UK_IRR17_natural` give natural uranium its published 1 Bq/g. Where a plain
    row also exists, the plain row is always used and the whole-chain value is
    inert.

`UK_EPR16_norm` is where this bites. Schedule 23 Part 3 Table 1 gives U-238 a
plain 5 Bq/g and a whole-chain 0.5 Bq/g, and the whole-chain value is never
applied:

```python
from radiological_material_clearance_finder import Material, clearance_index, get_limit_set

norm = get_limit_set("UK_EPR16_norm")
print("Table 1 U-238 plain row     ", norm.limits["U238"], "Bq/g")
print("Table 1 U-238 whole-chain row", norm.limits["U238_sec"], "Bq/g")

ore = Material.from_specific_activities(
    {n: 1.0 for n in ("U238", "Th234", "U234", "Th230", "Ra226", "Pb210", "Po210")},
    name="uranium ore",
)
r = clearance_index(ore, "UK_EPR16_norm")
print("limit applied to U-238       ", r.limits_used["U238"], "Bq/g")
```

```text
Table 1 U-238 plain row      5.0 Bq/g
Table 1 U-238 whole-chain row 0.5 Bq/g
limit applied to U-238        5.0 Bq/g
```

The same applies to Th-232 (5 against 0.5) and U-235 (5 against 1) in that set.
For a NORM material genuinely in secular equilibrium this is a factor of ten
non-conservative, and the whole-chain value has to be applied by hand. Read the
values off `LimitSet.limits` under the `_sec` keys.

## Catch-all limits, and the activity that leaves the sum

A nuclide with no row in the table is handled in one of two completely different
ways depending on the jurisdiction, and the difference is easy to miss because
the index looks equally plausible either way.

Three sets define a catch-all limit for unlisted nuclides:

| Set | Catch-all |
| --- | --- |
| `UK_EPR16_out_of_scope` | 0.01 Bq/g |
| `UK_IRR17_notification` | 0.01 Bq/g |
| `UK_IRR17_registration` | 0.1 Bq/g |

IRR 2017 qualifies its catch-all with "unless the Executive has approved some
other quantity for that radionuclide". Every other set, including the other UK
sets and all the German, US, EU and IAEA sets, has none. There, an unlisted
nuclide's activity is simply absent from the sum.

W-188 makes the three-way difference concrete. It has a row in StrlSchV Anlage 4
Tabelle 1, but none in the EU, IAEA or IRR 2017 tables.

```python
from radiological_material_clearance_finder import Material, clearance_index

# W-188 has a row in StrlSchV Anlage 4 Tabelle 1 but none in the EU, IAEA or
# UK IRR 2017 tables.
armour = Material.from_specific_activities(
    {"W181": 5.0, "W185": 20.0, "W188": 3.0, "Ta182": 0.05},
    name="tungsten armour",
)

for name in ("StrlSchV_unrestricted", "UK_IRR17_notification", "EU_BSS_clearance"):
    r = clearance_index(armour, name)
    print(f"{name:24s} index {r.index:>9.4g}  "
          f"defaulted {r.defaulted!s:<11} uncovered {r.uncovered!s:<18} "
          f"uncovered_fraction {r.uncovered_fraction:.2%}")
```

```text
StrlSchV_unrestricted    index      1.32  defaulted ()          uncovered {}                 uncovered_fraction 0.00%
UK_IRR17_notification    index       301  defaulted ('W188',)   uncovered {'W188': 0}      uncovered_fraction 0.00%
EU_BSS_clearance         index      1.02  defaulted ()          uncovered {'W188': 3.0}      uncovered_fraction 10.70%
```

Three treatments of one nuclide:

- **Germany** has a real row for it, 10 Bq/g, giving an index of 1.32.
- **IRR 2017** has no row, so W-188 takes the 0.01 Bq/g catch-all and dominates
  everything else, giving 301. It is listed in `defaulted`.
- **The EU** has no row and no catch-all, so W-188 contributes nothing at all.
  The index of 1.02 looks reassuringly similar to Germany's 1.32 while silently
  omitting 10.7% of the activity in the material.

!!! danger "Check `uncovered_fraction` before you trust an index"
    A comfortable index from a set with no catch-all can mean the material is
    clean, or it can mean the activity that would have failed it is not in the
    sum. Nothing in the index itself distinguishes the two.
    `ClearanceResult.uncovered`, `uncovered_activity` and `uncovered_fraction`
    exist to be read, not to be optional. See
    [the clearance index API](api/clearance-index.md) for the full set of
    fields.

Note the distinction between `uncovered` and `unlimited`. A nuclide in
`unlimited` is one the source *explicitly* places no limit on, which means it is
covered and contributes zero. Only the Fetter set uses this, for 20 nuclides. A
nuclide in `uncovered` is one the table never mentions, which is not the same
thing at all.

Passing `apply_default_limit=False` turns the catch-all off, moving those
nuclides from `defaulted` into `uncovered`. That shows what the listed rows
alone contribute; it is not what the regulation says.

## The EPR 2016 hundred second scope rule

EPR 2016 Schedule 23 Part 2 paragraph 7 places a substance outside the
regulation altogether when *none* of the radionuclides it contains has a
half-life exceeding 100 seconds. This is a test on the whole material, not a
per-nuclide filter, and it is the one rule here that can declare a material
acceptable while its index is enormous. Only `UK_EPR16_out_of_scope` carries it.

```python
from radiological_material_clearance_finder import Material, clearance_index, half_life

for name in ("N16", "O19", "F20", "Co60"):
    print(f"{name:6s} half-life {half_life(name):>12,.2f} s")
print()

# Activated coolant: every radionuclide present is shorter lived than 100 s.
coolant = Material.from_specific_activities(
    {"N16": 5.0e4, "O19": 2.0e3, "F20": 1.0e2}, name="coolant, at power"
)
# The same coolant carrying a trace of corrosion product.
with_crud = Material.from_specific_activities(
    {"N16": 5.0e4, "O19": 2.0e3, "F20": 1.0e2, "Co60": 1.0e-3},
    name="coolant + Co-60 trace",
)

for mat in (coolant, with_crud):
    r = clearance_index(mat, "UK_EPR16_out_of_scope")
    print(f"{mat.name:24s} out_of_scope={r.out_of_scope!s:<6} "
          f"index={r.index:>12.4g}  clearable={r.clearable}")
    print(f"{'':24s} still in the sum: {sorted(r.by_nuclide)}")
```

```text
N16    half-life         7.13 s
O19    half-life        26.88 s
F20    half-life        11.16 s
Co60   half-life 166,344,200.00 s

coolant, at power        out_of_scope=True   index=    5.21e+06  clearable=True
                         still in the sum: ['F20', 'N16', 'O19']
coolant + Co-60 trace    out_of_scope=False  index=    5.21e+06  clearable=False
                         still in the sum: ['Co60', 'F20', 'N16', 'O19']
```

Three things to take from this.

The index is *identical* in both cases, 5.21 million. Only `out_of_scope`
changes, and `clearable` follows it. One mBq/g of Co-60, contributing 0.01 to an
index of five million, flips the verdict from acceptable to not, because it is
the only nuclide present with a half-life over 100 seconds.

The short-lived nuclides stay in the sum in both cases. The rule does not remove
them individually. Once anything long-lived is present, the full N-16, O-19 and
F-20 activity is assessed against the 0.01 Bq/g catch-all as normal, which is
what produces the huge index.

A material with no radionuclides at all is not treated as out of scope. The rule
is a statement about how long-lived the radionuclides are, and with none present
there is nothing for it to be true of; returning `True` would force `clearable`
regardless of the index.

## The EU natural series expansion

EU 2013/59 Annex VII Table A Part 2 gives a single value for a whole natural
decay series, "for naturally occurring radionuclides in solid materials in
secular equilibrium with their progeny", rather than a value per nuclide. The
U-238 and Th-232 series each get 1 Bq/g, and K-40 is listed on its own at
10 Bq/g.

A value attached to a series head cannot be looked up by nuclide, so it is
expanded across the series members when the data is built. Membership is walked
from the IAEA Nuclear Data Section decay mode table by `tools/_series.py`,
following alpha, beta and electron capture branches, rather than being written
out by hand. Where a nuclide also has an explicit Part 1 row, the Part 1 row
wins, because an explicit row is more specific than a whole-series value.

The expansion accounts for exactly the difference between the EU and IAEA sets
in this package: 33 series members at 1 Bq/g plus K-40 at 10 Bq/g, 34 nuclides.
All 257 values the two sets have in common are identical.

```python
from radiological_material_clearance_finder import Material, clearance_index

# Natural uranium ore in secular equilibrium: every member of the U-238 series
# at the same activity concentration.
ore = Material.from_specific_activities(
    {n: 0.05 for n in (
        "U238", "Th234", "Pa234_m1", "U234", "Th230", "Ra226", "Rn222",
        "Po218", "Pb214", "Bi214", "Po214", "Pb210", "Bi210", "Po210",
    )},
    name="uranium ore",
)

for name in ("EU_BSS_clearance", "IAEA_GSR3_clearance"):
    r = clearance_index(ore, name)
    print(f"{name:22s} index {r.index:>7.4g}  clearable={r.clearable!s:<6} "
          f"nuclides in sum {len(r.by_nuclide):>2d}  uncovered {len(r.uncovered):>2d}  "
          f"uncovered_fraction {r.uncovered_fraction:.0%}")
```

```text
EU_BSS_clearance       index     0.7  clearable=True   nuclides in sum 14  uncovered 14  uncovered_fraction 0%
IAEA_GSR3_clearance    index       0  clearable=True   nuclides in sum  0  uncovered 14  uncovered_fraction 100%
```

!!! danger "The IAEA result is the failure mode this page is about"
    `IAEA_GSR3_clearance` returns an index of exactly zero and `clearable` is
    `True`, for a material that is entirely radioactive. The shipped IAEA set
    carries GSR Part 3 Schedule I **Table I.2** only, which covers artificial
    radionuclides. The natural series values in Table I.3 are not included, so
    not one of these 14 nuclides has a limit and all of the activity sits in
    `uncovered`. The 100% `uncovered_fraction` is the only signal that the zero
    is meaningless.

This is a limitation of the shipped data rather than of the regulation. For
natural material use `EU_BSS_clearance`, which carries the expanded Part 2
values, or the UK natural sets `UK_IRR17_natural` and `UK_EPR16_norm`.

## What is and is not verified

Being precise about provenance matters more here than anywhere else in the
package, so the state of verification is stated plainly.

**Machine-to-machine agreement is good.** The US sets agree with
`openmc.Material.waste_disposal_rating` to floating point round-off across 24
tests in CI. The IAEA PDF extraction and the EU XHTML parse agree on all 257
values they have in common, having been derived by completely independent paths.
All 81 Fetter lower bounds match the 1990 paper. The four regulatory tables
re-derived from live sources match the committed data.

**No human has checked any of it against the source regulations.** Every
agreement claimed above is one machine-produced table matching another. The
citations on this page come from the limit set metadata and the build scripts in
`tools/`, which were themselves written against the sources; they have not been
independently confirmed against the published legal text. Where this page cites
a Part, Table or paragraph, treat it as a pointer to check rather than as a
checked fact.

**These values are not regulatory advice.** The index this package computes is
the regulation's own arithmetic applied to the inventory you supply. Whether the
right limit set, the right column, the right waste route and the right
assumptions about your material were chosen is not something the code can
verify.

## See also

- [Material](api/material.md), and in particular why stable isotopes must be
  included in the inventory: the Bq/g denominator is the mass of everything
  passed.
- [The clearance index and `ClearanceResult`](api/clearance-index.md) for every
  field mentioned here.
- [Limit sets](api/limits.md) for the `LimitSet` fields
  (`secular_equilibrium`, `limits_secular_equilibrium`, `default_limit`,
  `min_half_life_scope`) and the provenance of each set.
- [Cooling](api/cooling.md) for `time_to_clear`, which interpolates
  log-linearly in the index and raises `IngrowthError` if the index climbs back
  above the threshold at a later cooling time.


## Two published values for one parent

A parent can carry two values with two different daughter lists, and picking the
wrong one is the largest single error this package has had.

EPR 2016 Schedule 23 lists uranium-238 twice in the out of scope column:

| Row | Limit | Daughters the value accounts for |
| --- | --- | --- |
| `U-238+` | 1 Bq/g | Th-234, Pa-234m, Pa-234 |
| `U-238sec` | 0.01 Bq/g | the full fourteen member chain, down to Po-210 |

The whole chain value is a hundred times stricter, because it is covering a
hundred times more. Applying the `+` value while excluding the `sec` list charges
the parent against a limit that accounts for only three of the fourteen
radionuclides removed from the sum.

The two lists are kept apart, and the value whose own list matches what is
present is the one applied:

```python
from radiological_material_clearance_finder import Material, clearance_index, get_limit_set

limit_set = get_limit_set("UK_EPR16_out_of_scope")

# Natural uranium in secular equilibrium: the whole chain is there.
activities = {"U238": 0.5}
activities.update({d: 0.5 for d in limit_set.secular_equilibrium_sec["U238"]})
result = clearance_index(Material.from_specific_activities(activities),
                         "UK_EPR16_out_of_scope")
print(result.index, result.clearable, result.limits_used["U238"])
```

```text
50.0 False 0.01
```

With only the short lived progeny present, the `+` value applies instead:

```python
result = clearance_index(
    Material.from_specific_activities({"U238": 0.5, "Th234": 0.5, "Pa234": 0.5}),
    "UK_EPR16_out_of_scope",
)
print(result.index, result.limits_used["U238"], sorted(result.excluded))
```

```text
0.5 1.0 ['Pa234', 'Th234']
```

`secular_equilibrium` holds the `+` list and `secular_equilibrium_sec` the whole
chain list, and the stricter value sits in `limits` under a `_sec` key.
