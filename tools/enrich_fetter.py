"""Add the ranges and unlimited entries from Fetter's own Table 2.

    python tools/enrich_fetter.py [--pdf PATH]

OpenMC's Fetter dict keeps only the conservative lower bound where the paper
gives a range, and drops the rows marked TMSA (theoretical maximum specific
activity, meaning effectively unlimited) entirely. Dropping them is the
misleading part: H-3, Kr-85, Zr-93 and others then look like nuclides the paper
says nothing about, when it says they need no limit at all.

This reads Table 2 out of the paper and adds two things to the existing limit
set without changing any value OpenMC agrees with: ``limits_upper`` for the top
of each range, and ``unlimited`` for the TMSA rows.

Source: Fetter, Cheng and Mann, "Long-term radioactive waste from fusion
reactors: Part II", Fusion Engineering and Design 13(2) 239-246 (1990),
doi:10.1016/0920-3796(90)90104-E.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from radiological_material_clearance_finder import nuclide as nuc  # noqa: E402

PDF_URL = (
    "https://fetter.it-prod-webhosting.aws.umd.edu/sites/default/files/fetter/files/"
    "1990-FED-RadWaste.pdf"
)
DATA = Path(__file__).resolve().parents[1] / "src" / "radiological_material_clearance_finder" / "data" / "limits"

_NUMBER = r"\d+\.?\d*\.?E[+-]\d+"
_ROW = re.compile(
    rf"^\s*([A-Z][a-z]?-\d{{1,3}}m?)\s+\S+\s*\S*\s+"
    rf"(TMSA|{_NUMBER}(?:\s*-\s*{_NUMBER})?)",
)


def fetch(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=180) as handle:
        destination.write_bytes(handle.read())


def to_text(pdf: Path) -> str:
    try:
        result = subprocess.run(
            ["pdftotext", "-layout", str(pdf), "-"],
            check=True, capture_output=True, text=True,
        )
    except FileNotFoundError:
        raise SystemExit("pdftotext is not installed (poppler-utils)") from None
    return result.stdout


def parse_table2(text: str) -> tuple[dict[str, float], dict[str, float], list[str]]:
    """Return ``(lower, upper, unlimited)`` from Table 2."""
    lower: dict[str, float] = {}
    upper: dict[str, float] = {}
    unlimited: list[str] = []

    for line in text.splitlines():
        match = _ROW.match(line)
        if not match:
            continue
        try:
            name, _marker = nuc.parse_regulatory(match.group(1))
        except nuc.NuclideNameError:
            continue
        value = match.group(2)
        if value == "TMSA":
            if name not in unlimited:
                unlimited.append(name)
            continue
        parts = [p.strip() for p in value.split("-") if p.strip()]
        # "6.E+02" or "6.E+02 - 6.E+03"; the exponent sign is inside the token,
        # so splitting on "-" is safe only after E+/E- have been consumed.
        parts = re.findall(_NUMBER, value)
        numbers = [float(p.replace(".E", "E")) for p in parts]
        lower[name] = numbers[0]
        if len(numbers) > 1:
            upper[name] = numbers[-1]
    return lower, upper, unlimited


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as workdir:
        pdf = args.pdf
        if pdf is None:
            pdf = Path(workdir) / "fetter.pdf"
            fetch(PDF_URL, pdf)
        lower, upper, unlimited = parse_table2(to_text(pdf))

    payload = json.loads((DATA / "us.json").read_text())
    entry = next(s for s in payload["sets"] if s["name"] == "Fetter")
    existing = {nuc.normalise(k): v for k, v in entry["limits"].items()}

    disagreeing = {
        name: (existing[name], lower[name])
        for name in set(existing) & set(lower)
        if abs(existing[name] - lower[name]) > 1e-9 * max(existing[name], lower[name])
    }
    print(f"parsed {len(lower)} limits, {len(upper)} ranges, {len(unlimited)} TMSA rows")
    print(f"OpenMC has {len(existing)}; {len(set(existing) & set(lower))} in common, "
          f"{len(disagreeing)} disagreeing")
    for name, (theirs, ours) in list(disagreeing.items())[:10]:
        print(f"    {name}: OpenMC {theirs}, paper {ours}")
    if disagreeing:
        raise SystemExit(
            "the paper and OpenMC disagree on a lower bound, so the parse is not "
            "trustworthy and nothing was written"
        )

    entry["limits_upper"] = {k: upper[k] for k in sorted(upper)}
    entry["unlimited"] = sorted(unlimited)
    entry["notes"] += (
        f" The paper gives {len(upper)} of these limits as a range, whose upper "
        f"bound is in limits_upper, and marks {len(unlimited)} nuclides TMSA "
        f"(theoretical maximum specific activity), meaning no limit applies. Those "
        f"are listed in unlimited so they are reported as covered rather than as "
        f"activity the index ignores."
    )
    (DATA / "us.json").write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {DATA / 'us.json'}")
    print(f"unlimited: {sorted(unlimited)}")


if __name__ == "__main__":
    main()
