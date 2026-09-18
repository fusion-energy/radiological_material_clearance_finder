"""Build the UK EPR 2016 Schedule 23 limit sets from legislation.gov.uk.

    python tools/build_uk_epr16.py

Schedule 23 of the Environmental Permitting (England and Wales) Regulations 2016
defines when a substance counts as radioactive material or radioactive waste at
all, which is the UK's out of scope, or clearance, test. Part 3 Table 2 gives the
concentration below which an artificial radionuclide is out of scope, Table 1
does the same for NORM industrial activities, and Table 3 lists the daughters
whose contribution a marked parent's value already includes.

Two features of this regulation that the German and US tables do not share are
carried into the limit sets rather than lost:

* Table 2 ends with a catch-all, so an unlisted artificial radionuclide is not
  unregulated, it takes 0.01 Bq/g, and
* Part 2 paragraph 7 puts a substance outside the regulation entirely when none
  of its radionuclides has a half-life exceeding 100 seconds.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from _tables import fetch as _fetch  # noqa: E402
from _tables import check_regulatory, parse_value, rows, tabulars  # noqa: E402
from radiological_material_clearance_finder import nuclide as nuc  # noqa: E402

URL = "https://www.legislation.gov.uk/ukdsi/2016/9780111150184/schedule/23/data.xml"
DATA = Path(__file__).resolve().parents[1] / "src" / "radiological_material_clearance_finder" / "data" / "limits"

CITATION = "Environmental Permitting (England and Wales) Regulations 2016, Schedule 23"
#: Part 2 paragraph 7, in seconds.
SHORT_HALF_LIFE_SCOPE = 100.0


def fetch(url: str) -> str:
    """Fetch the source, retrying on the transient failures CI runners see."""
    return _fetch(url)


def nuclide_rows(table_body: str):
    """Yield ``(canonical_key, marker, cells)`` for rows that start with a nuclide."""
    for cells in rows(table_body):
        if not cells:
            continue
        try:
            name, marker = nuc.parse_regulatory(cells[0])
        except nuc.NuclideNameError:
            continue
        key = f"{name}_sec" if marker == "sec" else name
        yield key, marker, cells


def parse_limits(table_body: str, column: int, label: str) -> tuple[dict, dict]:
    """Read one numeric column of a table, keyed by nuclide.

    A few nuclides appear twice, once plain and once marked "+", with different
    values: U-240 in Table 5 is 10^3 Bq/g plain and 10 Bq/g marked. Both are
    kept. The plain value is the limit, and the marked value is returned
    separately so it can be applied only when that parent's daughters are
    actually present, which is what the marker means.

    Returns:
        ``(limits, equilibrium_limits)`` keyed by nuclide.
    """
    limits: dict[str, float] = {}
    markers: dict[str, str | None] = {}
    alternatives: dict[str, float] = {}

    for key, marker, cells in nuclide_rows(table_body):
        if len(cells) <= column:
            continue
        value = parse_value(cells[column])
        if value is None:
            continue
        if key not in limits:
            limits[key], markers[key] = value, marker
            continue
        if limits[key] == value:
            continue
        # Keep the plain row's value as the limit and the "+" row's value as the
        # secular equilibrium alternative, whichever order the two rows appear in.
        if marker == "+" and markers[key] is None:
            alternatives[key] = value
        elif markers[key] == "+" and marker is None:
            alternatives[key] = limits[key]
            limits[key], markers[key] = value, marker
        else:
            raise ValueError(
                f"{label}: {key} appears twice with values {limits[key]} and {value} "
                f"and markers {markers[key]!r}, {marker!r}, which is not a case this "
                f"build script knows how to resolve"
            )
    check_regulatory(limits, label)
    return limits, alternatives


def parse_daughters(table_body: str) -> tuple[dict, dict]:
    """Read a secular equilibrium table, keeping the "+" and "sec" lists apart.

    A parent can appear twice with two different lists. U-238 has three short
    lived progeny under "U-238+" and a fourteen member chain under "U-238sec",
    and each list belongs to its own limit value: 1 Bq/g for the "+" row and
    0.01 Bq/g for the "sec" row in the out of scope column. Merging them lets
    the fourteen member chain be excluded while the parent is charged against
    the value that only accounts for three of them, which understates the index
    by a factor of a hundred.

    Returns:
        ``(plus, sec)``, each mapping a parent to the daughters its own value
        accounts for.
    """
    plus: dict[str, list[str]] = {}
    sec: dict[str, list[str]] = {}
    for cells in rows(table_body):
        if len(cells) < 2:
            continue
        try:
            parent, marker = nuc.parse_regulatory(cells[0])
        except nuc.NuclideNameError:
            continue
        daughters = sec if marker == "sec" else plus
        # Two Table 8 rows qualify their list with prose, as in "Where Ra-224+
        # is referred to in Table 5: Rn-220, Po-216, ...". Splitting on commas
        # without stripping that prefix glues it to the first daughter, which
        # then fails to parse and is dropped, silently losing Rn-220 and Rn-222.
        progeny = cells[1].split(":", 1)[-1] if ":" in cells[1] else cells[1]
        names = []
        for piece in re.split(r"[,;]", progeny):
            piece = piece.strip()
            if not piece:
                continue
            try:
                child, _ = nuc.parse_regulatory(piece)
            except nuc.NuclideNameError:
                continue
            names.append(child)
        if names:
            daughters.setdefault(parent, [])
            for child in names:
                if child not in daughters[parent]:
                    daughters[parent].append(child)
    return plus, sec


def find_catch_all(table_body: str) -> float | None:
    """Read Table 2's 'any other radionuclide' value.

    The value sits in a spanned cell that the XML emits as its own row, just
    above the row carrying the description, so this looks in both places rather
    than assuming a fixed column.
    """
    table = rows(table_body)
    for position, cells in enumerate(table):
        if not cells or not cells[0].lower().startswith("any other"):
            continue
        for cell in cells[1:]:
            value = _quiet_value(cell)
            if value is not None:
                return value
        if position > 0 and len(table[position - 1]) == 1:
            return _quiet_value(table[position - 1][0])
    return None


def _quiet_value(cell: str) -> float | None:
    """Parse a value, treating prose as absent rather than as an error."""
    try:
        return parse_value(cell)
    except ValueError:
        return None


def check_markers_have_daughters(table_body: str, maps: tuple, label: str) -> None:
    """Every marked parent must have a row in the matching daughter table."""
    plus, sec = maps
    missing = []
    for _key, marker, cells in nuclide_rows(table_body):
        if marker is None:
            continue
        parent, _ = nuc.parse_regulatory(cells[0])
        if parent not in (sec if marker == "sec" else plus):
            missing.append(cells[0])
    if missing:
        raise ValueError(
            f"{label}: {len(missing)} nuclide(s) marked '+' or 'sec' have no entry "
            f"in the secular equilibrium table: {missing[:8]}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cached", type=Path, help="local copy of the XML instead of downloading")
    args = parser.parse_args()

    document = args.cached.read_text() if args.cached else fetch(URL)
    tables = {number: (title, body) for number, title, body in tabulars(document)}
    for required in ("Table 1", "Table 2", "Table 3"):
        if required not in tables:
            raise SystemExit(f"{required} not found in the source document")

    _, table1 = tables["Table 1"]
    _, table2 = tables["Table 2"]
    _, table3 = tables["Table 3"]

    daughters = parse_daughters(table3)
    plus_daughters, sec_daughters = daughters
    check_markers_have_daughters(table2, daughters, "Table 2")
    check_markers_have_daughters(table1, daughters, "Table 1")

    artificial, artificial_alt = parse_limits(table2, 1, "EPR16 Table 2")
    catch_all = find_catch_all(table2)
    if catch_all is None:
        raise SystemExit("Table 2 catch-all row not found, which would change every index")

    norm_solid, norm_alt = parse_limits(table1, 1, "EPR16 Table 1 solid")

    sets = [
        {
            "name": "UK_EPR16_out_of_scope",
            "label": "UK out of scope concentration, artificial radionuclides in solids",
            "jurisdiction": "UK",
            "units": "Bq/g",
            "limits": artificial,
            "limits_secular_equilibrium": artificial_alt,
            "default_limit": catch_all,
            "secular_equilibrium": plus_daughters,
            "secular_equilibrium_sec": sec_daughters,
            "min_half_life_scope": SHORT_HALF_LIFE_SCOPE,
            "threshold": 1.0,
            "source": f"{CITATION}, Part 3 Table 2",
            "url": "https://www.legislation.gov.uk/ukdsi/2016/9780111150184/schedule/23",
            "retrieved": date.today().isoformat(),
            "notes": (
                "Below these concentrations a substance is not radioactive material "
                "or radioactive waste, which is the UK equivalent of clearance. The "
                f"catch-all of {catch_all} Bq/g applies to any artificial radionuclide "
                "the table does not list, so unlisted nuclides are limited rather than "
                "unregulated. Part 2 paragraph 7 places a substance outside the "
                "regulation when none of its radionuclides has a half-life exceeding "
                "100 seconds, which is a test on the whole material. These values "
                "are a separate derivation from the exemption values of EU directive "
                "2013/59/Euratom Annex VII, based on a 10 microsievert per year dose "
                "criterion via Euratom RP 122 part 1, and are generally ten times "
                "stricter, so they are not expected to match EU_BSS_clearance or "
                "UK_IRR17_notification."
            ),
        },
        {
            "name": "UK_EPR16_norm",
            "label": "UK out of scope concentration, NORM industrial activities in solids",
            "jurisdiction": "UK",
            "units": "Bq/g",
            "limits": norm_solid,
            "limits_secular_equilibrium": norm_alt,
            "secular_equilibrium": plus_daughters,
            "secular_equilibrium_sec": sec_daughters,
            "threshold": 1.0,
            "source": f"{CITATION}, Part 3 Table 1",
            "url": "https://www.legislation.gov.uk/ukdsi/2016/9780111150184/schedule/23",
            "retrieved": date.today().isoformat(),
            "notes": (
                "Applies to naturally occurring radionuclides in a NORM industrial "
                "activity, not to artificial ones. Entries keyed with a _sec suffix "
                "are the whole decay chain in secular equilibrium; the plain key is "
                "the value covering only the daughters listed in Table 3."
            ),
        },
    ]

    if "Table 5" in tables:
        # Part 6 paragraph 29 sends the "+" and "sec" markers in Table 5 to
        # Table 8, not to the Part 3 Table 3 used above. The two tables really do
        # differ: Table 3 has 57 parents and Table 8 has 31, and sixteen parents
        # they share carry different daughters. Table 3 gives U-235 twelve
        # daughters where Table 8 gives it only Th-231, so using the wrong one
        # deletes activity from the sum that nothing accounts for.
        if "Table 8" not in tables:
            raise SystemExit(
                "Table 5 is present but Table 8 is not, so the Part 6 secular "
                "equilibrium map cannot be built. The source layout has changed."
            )
        _, table5 = tables["Table 5"]
        _, table8 = tables["Table 8"]
        part6_daughters = parse_daughters(table8)
        part6_plus, part6_sec = part6_daughters
        check_markers_have_daughters(table5, part6_daughters, "Table 5")
        material_conc, material_alt = parse_limits(table5, 2, "EPR16 Table 5 concentration")
        if material_conc:
            sets.append(
                {
                    "name": "UK_EPR16_exempt_material",
                    "label": "UK exemption concentration for keeping and using radioactive material",
                    "jurisdiction": "UK",
                    "units": "Bq/g",
                    "limits": material_conc,
                    "limits_secular_equilibrium": material_alt,
                    "secular_equilibrium": part6_plus,
                    "secular_equilibrium_sec": part6_sec,
                    "threshold": 1.0,
                    "source": f"{CITATION}, Part 6 Table 5",
                    "url": "https://www.legislation.gov.uk/ukdsi/2016/9780111150184/schedule/23",
                    "retrieved": date.today().isoformat(),
                    "notes": (
                        "Exemption from the need for a permit to keep and use "
                        "radioactive material. The regulation pairs this concentration "
                        "with a maximum total activity on the premises, which is a "
                        "quantity limit rather than a concentration and is not part of "
                        "this index. The secular equilibrium daughters come from "
                        "Part 6 Table 8, which paragraph 29 designates for this Part, "
                        "and not from the Part 3 Table 3 used by the out of scope sets."
                    ),
                }
            )

    payload = {"_generated_by": "tools/build_uk_epr16.py", "_url": URL, "sets": sets}
    (DATA / "uk_epr16.json").write_text(json.dumps(payload, indent=1) + "\n")

    for entry in sets:
        print(f"{entry['name']:32} {len(entry['limits']):4} limits")
    print(f"catch-all {catch_all} Bq/g, {len(plus_daughters)} '+' parents, "
          f"{len(sec_daughters)} 'sec' parents")
    print(f"wrote {DATA / 'uk_epr16.json'}")


if __name__ == "__main__":
    main()
