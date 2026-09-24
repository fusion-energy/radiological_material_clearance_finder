"""Build the UK IRR 2017 Schedule 7 limit sets from legislation.gov.uk.

    python tools/build_uk_irr17.py

Schedule 7 of the Ionising Radiations Regulations 2017 gives the activity
concentrations below which work with a radionuclide needs no notification or
registration. Part 1 covers artificial radionuclides and natural ones processed
for their radioactive properties, Part 2 covers unprocessed natural ones.

Three quirks of this source are handled explicitly rather than skipped, because
skipping any of them loses real data:

* footnote markers are glued to nuclide names, so Part 2's potassium row reads
  ``K-40`` followed by a superscript 1,
* some nuclides appear only in a chemical form, and tritium is one of them: the
  table's only hydrogen row is ``H-3 (tritiated compounds)``, so discarding
  qualified rows would drop tritium from a fusion inventory entirely, and
* the progeny lists here differ from the ones in EPR 2016 Schedule 23, giving
  Zr-95 to Nb-95 where EPR16 gives Zr-95 to Nb-95m, so the two must not share a
  table.
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
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from _tables import fetch as _fetch  # noqa: E402
from _tables import check_regulatory, parse_value, rows, tabulars  # noqa: E402
from radiological_material_clearance_finder import nuclide as nuc  # noqa: E402

URL = "https://www.legislation.gov.uk/uksi/2017/1075/schedule/7/data.xml"
DATA = Path(__file__).resolve().parents[1] / "crates" / "radiological-material-clearance-finder" / "data" / "limits"
CITATION = "Ionising Radiations Regulations 2017, Schedule 7"

# A trailing superscript on a name is a footnote reference, not part of the name.
_FOOTNOTE = re.compile(r"\^\d+\s*$")
_QUALIFIER = re.compile(r"^(?P<base>[^()]+?)\s*\((?P<qualifier>[^)]*)\)\s*$")


def fetch(url: str) -> str:
    """Fetch the source, retrying on the transient failures CI runners see."""
    return _fetch(url)


def read_name(cell: str) -> tuple[str, str | None, str | None]:
    """Read a nuclide name cell.

    Returns:
        ``(canonical_key, marker, qualifier)``, or raises if the cell is not a
        nuclide. A qualifier is a chemical or source form such as ``monoxide``.
    """
    text = _FOOTNOTE.sub("", cell).strip()
    if "/" in text:
        # Neutron sources such as Pu-239/Be-91 are not a single nuclide.
        raise nuc.NuclideNameError(f"{cell!r} names a source, not a nuclide")
    qualifier = None
    match = _QUALIFIER.match(text)
    if match:
        text = match.group("base")
        qualifier = match.group("qualifier")
    name, marker = nuc.parse_regulatory(text)
    key = f"{name}_sec" if marker == "sec" else name
    return key, marker, qualifier


def parse_limits(table: list[list[str]], column: int, label: str) -> tuple[dict, dict, dict]:
    """Read one numeric column, preferring unqualified rows over qualified ones.

    A nuclide listed both plain and marked "+" gets two different values, and
    they are kept apart rather than one overwriting the other. U-240 is the case
    that bites: the notification column gives 0.01 Bq/g plain and 100 Bq/g
    marked, so letting the marked row win makes the limit ten thousand times too
    lenient for U-240 on its own. Per the regulation's own note the marked value
    applies only when the progeny are there, which is what
    ``limits_secular_equilibrium`` expresses.

    Returns:
        ``(limits, qualifiers, equilibrium)``, where qualifiers records the
        chemical form for any nuclide whose only row carried one, and
        equilibrium holds the marked value for a nuclide that also has a plain
        row.
    """
    limits: dict[str, float] = {}
    qualifiers: dict[str, str] = {}
    markers: dict[str, str | None] = {}
    equilibrium: dict[str, float] = {}

    for cells in table:
        if not cells or len(cells) <= column:
            continue
        try:
            key, marker, qualifier = read_name(cells[0])
        except nuc.NuclideNameError:
            continue
        value = parse_value(cells[column])
        if value is None:
            continue

        if key not in limits:
            limits[key], markers[key] = value, marker
            if qualifier:
                qualifiers[key] = qualifier
            continue
        # An unqualified row always wins over a qualified one.
        if qualifier is None and key in qualifiers:
            limits[key], markers[key] = value, marker
            qualifiers.pop(key)
            continue
        if qualifier is not None or limits[key] == value:
            continue
        if marker == "+" and markers[key] is None:
            equilibrium[key] = value
        elif markers[key] == "+" and marker is None:
            equilibrium[key] = limits[key]
            limits[key], markers[key] = value, marker
        else:
            raise ValueError(
                f"{label}: {key} appears twice with values {limits[key]} and "
                f"{value}, markers {markers[key]!r} and {marker!r}"
            )
    check_regulatory(limits, label)
    return limits, qualifiers, equilibrium


def parse_progeny(table: list[list[str]]) -> dict[str, list[str]]:
    """Read a parent and progeny table into parent to daughter names."""
    progeny: dict[str, list[str]] = {}
    for cells in table:
        if len(cells) < 2:
            continue
        try:
            parent, _, _ = read_name(cells[0])
        except nuc.NuclideNameError:
            continue
        names = []
        for piece in re.split(r"[,;]", cells[1]):
            piece = piece.strip()
            if not piece:
                continue
            try:
                child, _, _ = read_name(piece)
            except nuc.NuclideNameError:
                continue
            if child not in names:
                names.append(child)
        if names:
            progeny[parent] = names
    return progeny


def find_catch_all(table: list[list[str]], column: int) -> float | None:
    """Read the 'other radionuclides not listed above' row.

    Its values sit in the following row, since the label spans the width of the
    table.
    """
    for position, cells in enumerate(table):
        if not cells or not cells[0].lower().startswith("other radionuclides"):
            continue
        for candidate in (cells[1:], table[position + 1] if position + 1 < len(table) else []):
            if len(candidate) >= column:
                # The label row has no leading value cell, so the following row
                # is offset by one relative to the header columns.
                value = _quiet_value(candidate[column - 1])
                if value is not None:
                    return value
    return None


def _quiet_value(cell: str) -> float | None:
    try:
        return parse_value(cell)
    except ValueError:
        return None


def check_markers_have_progeny(table: list[list[str]], progeny: dict, label: str) -> None:
    """Every "+" parent must appear in the progeny table.

    Only "+" is checked. A "sec" entry means the whole decay chain is in secular
    equilibrium, which the regulation does not enumerate because it is the
    entire natural series, so having no progeny row is correct for those.
    """
    missing = []
    for cells in table:
        if not cells:
            continue
        try:
            key, marker, _ = read_name(cells[0])
        except nuc.NuclideNameError:
            continue
        if marker != "+":
            continue
        parent = key
        if parent not in progeny:
            missing.append(cells[0])
    if missing:
        raise ValueError(f"{label}: marked nuclides with no progeny row: {missing[:8]}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cached", type=Path)
    args = parser.parse_args()

    document = args.cached.read_text() if args.cached else fetch(URL)
    tables = [rows(body) for _number, _title, body in tabulars(document)]
    if len(tables) < 4:
        raise SystemExit(f"expected at least 4 tables in Schedule 7, found {len(tables)}")

    part1, part1_progeny, part2, part2_progeny = tables[0], tables[1], tables[2], tables[3]
    progeny1 = parse_progeny(part1_progeny)
    progeny2 = parse_progeny(part2_progeny)
    check_markers_have_progeny(part1, progeny1, "IRR17 Part 1")
    check_markers_have_progeny(part2, progeny2, "IRR17 Part 2")

    notification, qualifiers, notification_eq = parse_limits(
        part1, 1, "IRR17 Part 1 notification"
    )
    registration, _, registration_eq = parse_limits(
        part1, 3, "IRR17 Part 1 registration"
    )
    natural, _, natural_eq = parse_limits(part2, 1, "IRR17 Part 2 natural")

    catch_all_notification = find_catch_all(part1, 1)
    catch_all_registration = find_catch_all(part1, 3)
    if catch_all_notification is None:
        raise SystemExit("IRR17 Part 1 catch-all not found, which would change every index")

    qualifier_note = (
        " Nuclides whose only row in the regulation carries a chemical form take "
        "that row: " + ", ".join(f"{k} ({v})" for k, v in sorted(qualifiers.items())) + "."
        if qualifiers
        else ""
    )
    common = {
        "jurisdiction": "UK",
        "units": "Bq/g",
        "url": "https://www.legislation.gov.uk/uksi/2017/1075/schedule/7",
        "retrieved": date.today().isoformat(),
        "threshold": 1.0,
    }
    sets = [
        {
            "name": "UK_IRR17_notification",
            "label": "UK exemption from notification, artificial radionuclides",
            "limits": notification,
            "limits_secular_equilibrium": notification_eq,
            "default_limit": catch_all_notification,
            "secular_equilibrium": progeny1,
            "source": f"{CITATION}, Part 1 column 2",
            "notes": (
                "Concentration below which work with the radionuclide needs no "
                "notification, and no registration for amounts over 1000 kg. The "
                f"catch-all of {catch_all_notification} Bq/g covers radionuclides the "
                "table does not list, and the regulation qualifies it with 'unless the "
                "Executive has approved some other quantity for that radionuclide'."
                + qualifier_note
            ),
            **common,
        },
        {
            "name": "UK_IRR17_registration",
            "label": "UK exemption from registration, amounts up to 1000 kg",
            "limits": registration,
            "limits_secular_equilibrium": registration_eq,
            "default_limit": catch_all_registration,
            "secular_equilibrium": progeny1,
            "source": f"{CITATION}, Part 1 column 4",
            "notes": (
                "Concentration below which registration is not needed for amounts of "
                "radioactive material up to 1000 kg." + qualifier_note
            ),
            **common,
        },
        {
            "name": "UK_IRR17_natural",
            "label": "UK exemption from notification, unprocessed natural radionuclides",
            "limits": natural,
            "limits_secular_equilibrium": natural_eq,
            "secular_equilibrium": progeny2,
            "source": f"{CITATION}, Part 2 column 2",
            "notes": (
                "Applies to naturally occurring radionuclides that are not processed "
                "for their radioactive, fissile or fertile properties, whether or not "
                "they are in secular equilibrium with their progeny."
            ),
            **common,
        },
    ]

    payload = {"_generated_by": "tools/build_uk_irr17.py", "_url": URL, "sets": sets}
    (DATA / "uk_irr17.json").write_text(json.dumps(payload, indent=1) + "\n")
    for entry in sets:
        print(f"{entry['name']:28} {len(entry['limits']):4} limits, "
              f"default {entry.get('default_limit')}")
    print(f"progeny: Part 1 {len(progeny1)} parents, Part 2 {len(progeny2)} parents")
    print(f"qualified-only nuclides: {qualifiers}")
    print(f"wrote {DATA / 'uk_irr17.json'}")


if __name__ == "__main__":
    main()
