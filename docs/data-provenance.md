# Where the data comes from

Every limit in this package is generated from an official source by a script in
`tools/`, and each limit set records the source, the URL and the date it was
retrieved:

```python
from radiological_material_clearance_finder import get_limit_set

limit_set = get_limit_set("StrlSchV_unrestricted")
limit_set.source     # 'Strahlenschutzverordnung (StrlSchV) 2018, Anlage 4 Tabelle 1, Spalte 3'
limit_set.url        # 'https://www.gesetze-im-internet.de/strlschv_2018/anlage_4.html'
limit_set.retrieved  # ISO date
```

Nothing reaches the network at runtime. The tables are shipped as JSON inside the
package.

## The sources

| Sets | Source | How it is read |
| --- | --- | --- |
| `UK_EPR16_*` | Environmental Permitting (England and Wales) Regulations 2016, Schedule 23 | XML from legislation.gov.uk |
| `UK_IRR17_*` | Ionising Radiations Regulations 2017, Schedule 7 | XML from legislation.gov.uk |
| `StrlSchV_*` | Strahlenschutzverordnung 2018, Anlage 4 | HTML from gesetze-im-internet.de |
| `EU_BSS_clearance` | Council Directive 2013/59/Euratom, Annex VII Table A | consolidated XHTML from EUR-Lex |
| `IAEA_GSR3_clearance` | IAEA GSR Part 3, Schedule I Table I.2 | PDF, via `pdftotext -layout` |
| `Fetter`, `NRC_*` | Fetter et al. (1990); 10 CFR 61.55 | extracted from OpenMC's `waste.py` with `ast` |

Decay data is ENDF/B-VIII.0 half-lives, AME2020 atomic masses, and decay modes
from the IAEA Nuclear Data Section.

## Every source writes its numbers differently

This is the real hazard, and none of it fails loudly. A misread parses fine and
produces a table that is quietly wrong:

- The two UK XML sources mark exponents with `<Superior>`, so flattening tags to
  text turns 10<sup>2</sup> into the literal digits `102`.
- The German source instead writes `1 E-1`, with a space, and uses decimal commas
  in its half-life column. It is ISO-8859-1, and decoding it as UTF-8 fails on
  the first umlaut.
- EUR-Lex uses comma decimal separators and non-breaking spaces as thousands
  separators, and emits the parent and progeny footnote table inside the same
  HTML table as the data, so its rows look like limit rows.
- The IAEA footnote is set as two column pairs on shared lines, and a long
  progeny list wraps onto lines that already carry an unrelated row in the other
  column.

Every build script therefore asserts that each limit is one significant figure
times a power of ten, which is the shape all of these regulations use without
exception. That single assertion catches the whole class: `102` fails it
immediately, where `100` would have passed unnoticed.

## Regenerating

```bash
python tools/build_uk_epr16.py       # legislation.gov.uk XML
python tools/build_uk_irr17.py
python tools/build_de_strlschv.py    # gesetze-im-internet.de
python tools/build_eu_bss.py         # EUR-Lex consolidated XHTML
python tools/build_iaea.py           # IAEA PDF, needs pdftotext
python tools/build_us.py             # extracted from OpenMC's waste.py
python tools/build_decay_data.py     # AME2020 and ENDF/B-VIII.0
python tools/build_decay_modes.py    # IAEA Livechart
```

`tools/verify_tables.py` re-downloads every downloadable source and diffs it
against the shipped tables, so an amendment to a regulation surfaces
deliberately. It runs weekly rather than on every pull request, because it
depends on reaching several government websites and a check that goes red for
network reasons teaches everyone to ignore it.

The US tables have no downloadable source. They are checked instead against
frozen copies in `tests/data`, which run on every commit.

## What has and has not been verified

Verified, and re-checked on every CI run:

- The US sets agree with OpenMC's own implementation to floating point
  round-off, across the index, the per-nuclide breakdown and the waste class.
- The IAEA PDF extraction and the EU XHTML parse agree on all 257 values they
  have in common, by completely independent paths.
- All 81 Fetter lower bounds match the 1990 paper.
- All four downloadable tables, re-derived from the live sources, reproduce the
  shipped data exactly.

!!! warning "No human has checked these tables against the regulations"

    Every confirmation above is machine to machine. Adversarial review of this
    package found several defects where the arithmetic agreed with an
    independent implementation exactly while the table underneath was wrong, and
    every one of them made a material look more clearable than it is.

    Before relying on a number here for an actual release decision, check the
    limits that drive it against the published regulation. `result.limits_used`
    gives you exactly which values to check, and `get_limit_set(name).url` gives
    you where to check them.
