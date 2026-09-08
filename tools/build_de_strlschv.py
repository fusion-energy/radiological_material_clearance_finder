"""Build the German StrlSchV 2018 clearance limit sets from the official text.

    python tools/build_de_strlschv.py [--verify-against PATH]

Anlage 4 of the Strahlenschutzverordnung 2018 carries the clearance values
(Freigabewerte) for each release pathway in Tabelle 1, and the parent to
daughter list behind the "+" marker in Tabelle 2. Parsing the regulation
directly means the daughter handling is the regulation's own, with no depletion
chain and no half-life heuristic involved.

Tabelle 1 numbers the columns it calls Spalten, and a cell's index is its Spalte
minus one. The pathway to Spalte mapping below is the regulation's own numbering.

``--verify-against`` diffs the result against the transcription in
openmc-dev/openmc#3898, which is an independent reading of the same table.
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import sys
import urllib.request
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from _tables import check_regulatory, clean, parse_value, rows  # noqa: E402
from radiological_material_clearance_finder import nuclide as nuc  # noqa: E402

URL = "https://www.gesetze-im-internet.de/strlschv_2018/anlage_4.html"
DATA = Path(__file__).resolve().parents[1] / "src" / "radiological_material_clearance_finder" / "data" / "limits"
CITATION = "Strahlenschutzverordnung (StrlSchV) 2018, Anlage 4 Tabelle 1"

#: Pathway name to (Spalte, English label, the regulation's own column heading).
PATHWAYS = {
    "StrlSchV_unrestricted": (
        3,
        "Unrestricted clearance of solid and liquid substances",
        "Freigrenze, uneingeschraenkte Freigabe von festen u. fluessigen Stoffen in Bq/g",
    ),
    "StrlSchV_rubble": (
        6,
        "Building rubble, more than 1000 t per year",
        "Bauschutt von mehr als 1000 Mg/a in Bq/g",
    ),
    "StrlSchV_soil": (7, "Soil surfaces", "Bodenflaechen in Bq/g"),
    "StrlSchV_landfill_100": (
        8,
        "Landfill, site taking up to 100 t per year",
        "festen Stoffen bis zu 100 Mg/a zur Beseitigung auf Deponien in Bq/g",
    ),
    "StrlSchV_incineration_100": (
        9,
        "Incineration, plant taking up to 100 t per year",
        "Stoffen bis zu 100 Mg/a zur Beseitigung in Verbrennungsanlagen in Bq/g",
    ),
    "StrlSchV_landfill_1000": (
        10,
        "Landfill, site taking up to 1000 t per year",
        "festen Stoffen bis zu 1000 Mg/a zur Beseitigung auf Deponien in Bq/g",
    ),
    "StrlSchV_incineration_1000": (
        11,
        "Incineration, plant taking up to 1000 t per year",
        "Stoffen bis zu 1000 Mg/a zur Beseitigung in Verbrennungsanlagen in Bq/g",
    ),
    "StrlSchV_metal_recycling": (
        14,
        "Metal scrap for recycling",
        "Metallschrott zum Recycling in Bq/g",
    ),
}

#: Spalte 2, an exemption on total activity rather than concentration.
EXEMPTION_SPALTE = 2

PATHWAY_NOTE = (
    "Mg/a is Megagramm im Kalenderjahr, tonnes per calendar year, and refers to "
    "the annual mass of cleared material the receiving facility accepts, not to "
    "the mass of any one item. A lower annual throughput permits a higher "
    "specific activity, because the collective dose it can produce is smaller. "
    "Anlage 4 Tabelle 1 also tabulates a high activity source threshold in TBq "
    "and two surface contamination columns in Bq/cm2, which this package does "
    "not model because it works from a bulk inventory rather than a surface."
)


#: gesetze-im-internet.de serves ISO-8859-1, and decoding it as UTF-8 fails on
#: the first umlaut, which appears in the column headings this script matches on.
ENCODING = "iso-8859-1"


def fetch(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=120) as handle:
        return handle.read().decode(ENCODING)


def find_tables(document: str) -> tuple[list[list[str]], list[list[str]]]:
    """Return Tabelle 1 and Tabelle 2 as row lists, identified by their headers."""
    tables = [rows(block) for block in re.findall(r"<table.*?</table>", document, re.S)]
    main = daughters = None
    for table in tables:
        header = " ".join(table[0]) if table else ""
        if table and table[0] and table[0][0].strip() == "Radionuklid":
            main = table
        elif "Mutternuklid" in header and "Tochternuklide" in header:
            daughters = table
    if main is None or daughters is None:
        raise SystemExit(
            "could not find Tabelle 1 and Tabelle 2 by their headers in the source"
        )
    return main, daughters


def check_column_numbering(table: list[list[str]]) -> None:
    """Confirm the Spalte numbering row is where the column mapping assumes."""
    for row in table[:6]:
        if row[:4] == ["1", "2", "3", "4"]:
            return
    raise SystemExit(
        "the Spalte numbering row was not found, so the column mapping cannot be "
        "trusted. The source layout has changed."
    )


def parse_limits(table: list[list[str]], spalte: int, label: str) -> tuple[dict, dict]:
    """Read one Spalte of Tabelle 1, keyed by nuclide.

    Three nuclides are listed twice, once plain and once marked "+": Th-232 is
    10 Bq/g plain and 0.01 Bq/g marked in Spalte 3, while Zr-97 and Pb-212 carry
    the same value either way. Both are kept. The plain value is the limit, and
    the marked value is returned separately so it can be applied only when the
    parent's daughters are actually present, which is what the marker means.

    Returns:
        ``(limits, equilibrium_limits)`` keyed by nuclide.
    """
    column = spalte - 1
    limits: dict[str, float] = {}
    markers: dict[str, str | None] = {}
    alternatives: dict[str, float] = {}

    for row in table:
        if not row or len(row) <= column:
            continue
        try:
            name, marker = nuc.parse_regulatory(row[0])
        except nuc.NuclideNameError:
            continue
        value = parse_value(row[column])
        if value is None:
            continue
        if name not in limits:
            limits[name], markers[name] = value, marker
            continue
        if limits[name] == value:
            continue
        if marker == "+" and markers[name] is None:
            alternatives[name] = value
        elif markers[name] == "+" and marker is None:
            alternatives[name] = limits[name]
            limits[name], markers[name] = value, marker
        else:
            raise ValueError(
                f"{label}: {name} appears twice with values {limits[name]} and "
                f"{value}, markers {markers[name]!r} and {marker!r}"
            )
    check_regulatory(limits, label)
    return limits, alternatives


def parse_daughters(table: list[list[str]]) -> dict[str, list[str]]:
    """Read Tabelle 2, Mutternuklid to Tochternuklide."""
    out: dict[str, list[str]] = {}
    for row in table:
        if len(row) < 2:
            continue
        try:
            parent, _ = nuc.parse_regulatory(row[0])
        except nuc.NuclideNameError:
            continue
        names = []
        for piece in re.split(r"[,;]", row[1]):
            piece = piece.strip()
            if not piece:
                continue
            try:
                child, _ = nuc.parse_regulatory(piece)
            except nuc.NuclideNameError:
                continue
            if child not in names:
                names.append(child)
        if names:
            out[parent] = names
    return out


def check_markers_have_daughters(table: list[list[str]], daughters: dict) -> None:
    missing = []
    for row in table:
        if not row:
            continue
        try:
            name, marker = nuc.parse_regulatory(row[0])
        except nuc.NuclideNameError:
            continue
        if marker == "+" and name not in daughters:
            missing.append(row[0])
    if missing:
        raise ValueError(
            f"{len(missing)} nuclide(s) marked '+' in Tabelle 1 have no row in "
            f"Tabelle 2: {missing[:8]}"
        )


def verify_against_pr(path: Path, built: dict[str, dict[str, float]]) -> None:
    """Diff the parsed tables against the transcription in openmc PR #3898."""
    text = path.read_text()
    print("\nverifying against", path)
    for name, limits in built.items():
        start = text.find(f"== '{name}'")
        if start < 0:
            continue
        end = text.find("units = 'Bq/g'", start)
        block = text[start:end]
        literal = block[block.find("{") : block.rfind("}") + 1]
        theirs = {nuc.normalise(k): v for k, v in ast.literal_eval(literal).items()}
        only_ours = sorted(set(limits) - set(theirs))
        only_theirs = sorted(set(theirs) - set(limits))
        differing = sorted(
            k for k in set(limits) & set(theirs)
            if abs(limits[k] - theirs[k]) > 1e-12 * max(abs(limits[k]), abs(theirs[k]))
        )
        status = "identical" if not (only_ours or only_theirs or differing) else "DIFFERS"
        print(f"  {name:28} ours {len(limits):4}  PR {len(theirs):4}  {status}")
        if only_ours:
            print(f"      only in ours ({len(only_ours)}): {only_ours[:10]}")
        if only_theirs:
            print(f"      only in PR   ({len(only_theirs)}): {only_theirs[:10]}")
        if differing:
            print(f"      differing values: "
                  f"{[(k, limits[k], theirs[k]) for k in differing[:6]]}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cached", type=Path)
    parser.add_argument("--verify-against", type=Path, help="_waste_pr3898.py to diff against")
    args = parser.parse_args()

    document = args.cached.read_text(encoding=ENCODING) if args.cached else fetch(URL)
    main_table, daughter_table = find_tables(document)
    check_column_numbering(main_table)

    daughters = parse_daughters(daughter_table)
    check_markers_have_daughters(main_table, daughters)

    built: dict[str, dict[str, float]] = {}
    sets = []
    for name, (spalte, label, german) in PATHWAYS.items():
        limits, alternatives = parse_limits(main_table, spalte, name)
        built[name] = limits
        sets.append(
            {
                "name": name,
                "label": label,
                "jurisdiction": "Germany",
                "units": "Bq/g",
                "limits": limits,
                "limits_secular_equilibrium": alternatives,
                "secular_equilibrium": daughters,
                "threshold": 1.0,
                "source": f"{CITATION}, Spalte {spalte}",
                "url": URL,
                "retrieved": date.today().isoformat(),
                "notes": f"German column heading: {german}. {PATHWAY_NOTE}",
            }
        )
        print(f"{name:28} Spalte {spalte:2}  {len(limits):4} limits")

    exemption, exemption_alt = parse_limits(main_table, EXEMPTION_SPALTE, "StrlSchV_exemption_activity")
    sets.append(
        {
            "name": "StrlSchV_exemption_activity",
            "label": "German exemption limit on total activity (Freigrenze in Bq)",
            "jurisdiction": "Germany",
            "units": "Bq",
            "limits": exemption,
            "limits_secular_equilibrium": exemption_alt,
            "secular_equilibrium": daughters,
            "threshold": 1.0,
            "source": f"{CITATION}, Spalte {EXEMPTION_SPALTE}",
            "url": URL,
            "retrieved": date.today().isoformat(),
            "notes": (
                "A limit on total activity rather than concentration, so the "
                "material needs an absolute mass or volume for this to be "
                "evaluated. The regulation applies the exemption when either this "
                "or the Spalte 3 concentration limit is met."
            ),
        }
    )
    print(f"{'StrlSchV_exemption_activity':28} Spalte {EXEMPTION_SPALTE:2}  {len(exemption):4} limits")
    print(f"Tabelle 2: {len(daughters)} parent nuclides")

    payload = {"_generated_by": "tools/build_de_strlschv.py", "_url": URL, "sets": sets}
    (DATA / "de_strlschv.json").write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")
    print(f"wrote {DATA / 'de_strlschv.json'}")

    if args.verify_against:
        verify_against_pr(args.verify_against, built)


if __name__ == "__main__":
    main()
