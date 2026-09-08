"""Build the EU and IAEA clearance limit sets from the consolidated BSS text.

    python tools/build_eu_bss.py

Annex VII Table A of Council Directive 2013/59/Euratom gives the activity
concentrations for exemption and clearance. Part 1 covers artificial
radionuclides in kBq/kg, which is numerically Bq/g, and Part 2 covers natural
ones. These are the values IAEA GSR Part 3 Schedule I Table I.2 carries, and
that the UK IRR 2017 and EPR 2016 tables transcribe, so this set doubles as a
cross-check on the others.

The EUR-Lex XHTML has two traps. The parent and progeny footnote table is
emitted inside the *same* ``<table>`` element as Table A Part 1, so its rows
look like data rows, and it is separated here by testing whether the second cell
is a number or a list of nuclides. Values use a comma decimal separator and a
non-breaking space for thousands.

The ``TXT/XML`` endpoint returns a CELLAR metadata notice rather than the legal
text, so the XHTML is the machine-readable source.
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
from _tables import check_regulatory, clean, parse_value  # noqa: E402
from radiological_material_clearance_finder import nuclide as nuc  # noqa: E402

URL = "https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:02013L0059-20140117"
DATA = Path(__file__).resolve().parents[1] / "src" / "radiological_material_clearance_finder" / "data" / "limits"

_TABLE = re.compile(r"<table[^>]*>.*?</table>", re.S)
_ROW = re.compile(r"<tr[^>]*>(.*?)</tr>", re.S)
_CELL = re.compile(r"<t[dh][^>]*>(.*?)</t[dh]>", re.S)


def fetch(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=180) as handle:
        return handle.read().decode("utf-8", errors="replace")


def table_rows(block: str) -> list[list[str]]:
    """Return the cleaned cells of each row in an HTML table block."""
    return [[clean(cell) for cell in _CELL.findall(row)] for row in _ROW.findall(block)]


def find_table_a_part1(document: str) -> list[list[str]]:
    """The table holding Annex VII Table A Part 1, found by its header."""
    for block in _TABLE.findall(document):
        rows = table_rows(block)
        if not rows or len(rows[0]) < 2:
            continue
        header = " ".join(rows[0]).lower()
        if "radionuclide" in header and "activity concentration" in header:
            return rows
    raise SystemExit("Annex VII Table A Part 1 not found in the source document")


def split_limits_and_progeny(rows: list[list[str]]) -> tuple[dict, dict]:
    """Separate the limit rows from the parent and progeny rows.

    Both live in one HTML table. A limit row's second cell is a number; a
    progeny row's second cell is a list of nuclide names. Testing the cell
    rather than the row position keeps this working if the layout shifts.
    """
    limits: dict[str, float] = {}
    progeny: dict[str, list[str]] = {}

    for cells in rows:
        if len(cells) < 2:
            continue
        # Footnote markers render as an empty "()" after the name.
        name_text = re.sub(r"\(\s*\)", "", cells[0]).strip()
        try:
            name, _marker = nuc.parse_regulatory(name_text)
        except nuc.NuclideNameError:
            continue

        try:
            value = parse_value(cells[1], decimal_comma=True)
        except ValueError:
            value = None

        if value is not None:
            if name in limits and limits[name] != value:
                raise ValueError(f"{name} appears twice with different values")
            limits[name] = value
            continue

        children = []
        for piece in re.split(r"[,;]", cells[1]):
            piece = piece.strip()
            if not piece:
                continue
            try:
                child, _ = nuc.parse_regulatory(piece)
            except nuc.NuclideNameError:
                children = []
                break
            if child not in children:
                children.append(child)
        if children:
            progeny[name] = children

    check_regulatory(limits, "EU BSS Annex VII Table A Part 1")
    return limits, progeny


NOTES = (
    "Values are given in the directive as kBq/kg, which is numerically identical "
    "to Bq/g, and are the same values as IAEA GSR Part 3 (2014) Schedule I Table "
    "I.2, verified identical for all 257 nuclides in common. UK IRR 2017 Schedule "
    "7 Part 1 transcribes them and agrees for 255 of the 257, differing only for "
    "Na-24 (0.1 against 1) and Pt-197 (10 against 1000). UK EPR 2016 Schedule 23 "
    "Table 2 is a separate derivation from a 10 microsievert per year dose "
    "criterion via Euratom RP 122, is generally ten times stricter, and is not "
    "expected to match. IAEA RS-G-1.7, which carried an earlier form of these "
    "values, has been superseded by GSG-17 and GSG-18."
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cached", type=Path)
    args = parser.parse_args()

    document = args.cached.read_text(encoding="utf-8", errors="replace") if args.cached else fetch(URL)
    rows = find_table_a_part1(document)
    limits, progeny = split_limits_and_progeny(rows)

    if len(limits) < 250:
        raise SystemExit(
            f"only {len(limits)} limits parsed from Table A Part 1, which is far "
            f"short of the ~294 rows the table holds"
        )

    sets = [
        {
            "name": "EU_BSS_clearance",
            "label": "EU Basic Safety Standards exemption and clearance concentration",
            "jurisdiction": "EU",
            "units": "Bq/g",
            "limits": limits,
            "secular_equilibrium": progeny,
            "threshold": 1.0,
            "source": "Council Directive 2013/59/Euratom, Annex VII Table A Part 1",
            "url": URL,
            "retrieved": date.today().isoformat(),
            "notes": NOTES,
        }
    ]
    payload = {"_generated_by": "tools/build_eu_bss.py", "_url": URL, "sets": sets}
    (DATA / "eu_bss.json").write_text(json.dumps(payload, indent=1) + "\n")
    print(f"EU_BSS_clearance  {len(limits)} limits, {len(progeny)} progeny parents")
    print(f"wrote {DATA / 'eu_bss.json'}")


if __name__ == "__main__":
    main()
