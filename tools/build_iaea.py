"""Build the IAEA clearance limit set from GSR Part 3 Schedule I Table I.2.

    python tools/build_iaea.py

Table I.2 of IAEA General Safety Requirements Part 3 (2014) gives activity
concentrations in Bq/g for clearance of solid material. This is the table that
national regulations reference; RS-G-1.7, which carried an earlier form of these
values, has been superseded by GSG-17 and GSG-18, so it is not used here.

The IAEA publishes only PDFs, so this extracts the table with ``pdftotext
-layout``. The layout is regular: two Radionuclide and value column pairs per
page. Footnote markers are glued to the nuclide name, giving ``Sr-90a`` and
``Zn-69ma``, and thousands are separated by a space.

Extracting from a PDF is less certain than parsing XHTML, so the result is
checked against the EU BSS set, which carries the same values through a
machine-readable source.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
import urllib.request
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from _tables import check_regulatory  # noqa: E402
from radiological_material_clearance_finder import nuclide as nuc  # noqa: E402

URL = "https://www-pub.iaea.org/MTCD/Publications/PDF/Pub1578_web-57265295.pdf"
DATA = Path(__file__).resolve().parents[1] / "src" / "radiological_material_clearance_finder" / "data" / "limits"

TABLE_HEADING = "TABLE I.2. LEVELS FOR EXEMPTION OF BULK AMOUNTS OF"
NEXT_HEADING = "TABLE I.3."
PROGENY_HEADING = "Parent radionuclides, and their progeny"

#: The progeny footnote is laid out as two independent column pairs. This is the
#: character offset where the right pair starts.
_COLUMN_SPLIT = 51

_NUCLIDE_AT_START = re.compile(r"^([A-Z][a-z]?-\d{1,3}(?:[mn]\d?)?)\s+(.*)$")

# A nuclide followed by its value. The trailing [a-d] is a footnote marker, not
# part of the name, and is never a metastable state, which is m or n.
_PAIR = re.compile(
    r"([A-Z][a-z]?-\d{1,3}(?:[mn]\d?)?)([a-d])?\s+(\d{1,3}(?:\s\d{3})*(?:\.\d+)?)"
)


def fetch(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=300) as handle:
        destination.write_bytes(handle.read())


def to_text(pdf: Path) -> str:
    try:
        result = subprocess.run(
            ["pdftotext", "-layout", str(pdf), "-"],
            check=True, capture_output=True, text=True,
        )
    except FileNotFoundError:
        raise SystemExit(
            "pdftotext is not installed. It comes with poppler-utils and is "
            "needed because the IAEA publishes this table only as a PDF."
        ) from None
    return result.stdout


def parse_table(text: str) -> dict[str, float]:
    """Read Table I.2 out of the extracted text."""
    start = text.find(TABLE_HEADING)
    if start < 0:
        raise SystemExit("Table I.2 heading not found in the extracted text")
    end = text.find(NEXT_HEADING, start)
    section = text[start : end if end > 0 else len(text)]

    limits: dict[str, float] = {}
    for line in section.splitlines():
        for name_text, _footnote, value_text in _PAIR.findall(line):
            try:
                name, _marker = nuc.parse_regulatory(name_text)
            except nuc.NuclideNameError:
                continue
            value = float(value_text.replace(" ", ""))
            if name in limits and limits[name] != value:
                raise ValueError(
                    f"{name} read twice with different values, {limits[name]} and "
                    f"{value}, so the column layout was misread"
                )
            limits[name] = value
    check_regulatory(limits, "IAEA GSR Part 3 Table I.2")
    return limits


def parse_progeny(text: str) -> dict[str, list[str]]:
    """Read the parent and progeny footnote that follows Table I.2.

    The footnote is set as two independent column pairs on the same physical
    lines, and a long progeny list wraps onto following lines that already carry
    unrelated rows in the other column. U-232's chain wraps across three lines
    that share them with Ru-103 and Ru-105. So the running parent is tracked per
    column; tracking it globally attaches the continuation nuclides to whichever
    row happened to be printed alongside.
    """
    start = text.find(PROGENY_HEADING)
    if start < 0:
        raise SystemExit("the parent and progeny footnote was not found")
    end = text.find(NEXT_HEADING, start)
    section = text[start : end if end > 0 else len(text)]

    progeny: dict[str, list[str]] = {}
    current = [None, None]
    for line in section.splitlines()[1:]:
        for column, fragment in enumerate((line[:_COLUMN_SPLIT], line[_COLUMN_SPLIT:])):
            stripped = fragment.strip()
            if not stripped:
                continue
            match = _NUCLIDE_AT_START.match(stripped)
            if match and _looks_like_nuclide_list(match.group(2)):
                try:
                    parent, _ = nuc.parse_regulatory(match.group(1))
                except nuc.NuclideNameError:
                    continue
                current[column] = parent
                progeny.setdefault(parent, [])
                _extend(progeny[parent], match.group(2))
            elif current[column] and _looks_like_nuclide_list(stripped):
                _extend(progeny[current[column]], stripped)
            else:
                # Prose, a page number or a heading. Whatever the running parent
                # was, its list has ended.
                current[column] = None
    return {k: v for k, v in progeny.items() if v}


def _looks_like_nuclide_list(text: str) -> bool:
    """Whether a fragment is nothing but nuclide names and separators."""
    text = text.strip()
    if not text:
        return False
    pieces = [p.strip() for p in text.split(",") if p.strip()]
    if not pieces:
        return False
    for piece in pieces:
        try:
            nuc.parse_regulatory(piece)
        except nuc.NuclideNameError:
            return False
    return True


def _extend(target: list[str], fragment: str) -> None:
    for piece in fragment.split(","):
        piece = piece.strip()
        if not piece:
            continue
        try:
            child, _ = nuc.parse_regulatory(piece)
        except nuc.NuclideNameError:
            continue
        if child not in target:
            target.append(child)


def cross_check(limits: dict[str, float]) -> str:
    """Compare against the EU BSS set, which holds the same values."""
    from radiological_material_clearance_finder.limits import get_limit_set

    try:
        european = get_limit_set("EU_BSS_clearance").limits
    except KeyError:
        return "EU_BSS_clearance is not built, so no cross-check was possible."
    shared = set(limits) & set(european)
    differing = {
        name: (limits[name], european[name])
        for name in shared
        if abs(limits[name] - european[name]) > 1e-12 * max(limits[name], european[name])
    }
    print(f"cross-check against EU_BSS_clearance: {len(shared)} shared, "
          f"{len(differing)} differing")
    for name, (ours, theirs) in list(differing.items())[:10]:
        print(f"    {name}: IAEA {ours}, EU {theirs}")
    if differing:
        return (
            f"Differs from EU_BSS_clearance for {len(differing)} of {len(shared)} "
            f"shared nuclides."
        )
    return (
        f"Verified identical to EU_BSS_clearance for all {len(shared)} nuclides in "
        f"common, which is the expected result since the EU directive adopts these "
        f"values."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cached", type=Path, help="local copy of the PDF")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as workdir:
        pdf = args.cached
        if pdf is None:
            pdf = Path(workdir) / "gsr3.pdf"
            fetch(URL, pdf)
        extracted = to_text(pdf)
        limits = parse_table(extracted)
        progeny = parse_progeny(extracted)

    if len(limits) < 200:
        raise SystemExit(f"only {len(limits)} nuclides parsed, expected around 257")

    agreement = cross_check(limits)
    sets = [
        {
            "name": "IAEA_GSR3_clearance",
            "label": "IAEA clearance of solid material, artificial radionuclides",
            "jurisdiction": "IAEA",
            "units": "Bq/g",
            "limits": limits,
            "secular_equilibrium": progeny,
            "threshold": 1.0,
            "source": "IAEA GSR Part 3 (2014), Schedule I, Table I.2",
            "url": URL,
            "retrieved": date.today().isoformat(),
            "notes": (
                "Levels for exemption of bulk amounts of solid material without "
                "further consideration, and for clearance of solid material without "
                "further consideration. Extracted from the IAEA PDF, which is the "
                "only form published. IAEA RS-G-1.7 carried an earlier form of these "
                "values and has been superseded by GSG-17 and GSG-18. " + agreement
            ),
        }
    ]
    payload = {"_generated_by": "tools/build_iaea.py", "_url": URL, "sets": sets}
    (DATA / "iaea.json").write_text(json.dumps(payload, indent=1) + "\n")
    print(f"IAEA_GSR3_clearance  {len(limits)} limits, {len(progeny)} progeny parents")
    for parent in ("U232", "Zr97", "Cs137", "Pu244"):
        print(f"    {parent}: {progeny.get(parent)}")
    print(f"wrote {DATA / 'iaea.json'}")


if __name__ == "__main__":
    main()
