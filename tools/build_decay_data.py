"""Regenerate the vendored half-life and atomic mass tables.

    python tools/build_decay_data.py

Atomic masses come from AME2020, parsed from ``mass_1.mas20.txt``. The copy in
an OpenMC checkout is byte-identical to the IAEA AMDC original (md5
6d28b75833cf53c7cc230223f63da6f6), and the file is downloaded from AMDC when no
local copy is given.

Half-lives are the ENDF/B-VIII.0 decay sublibrary values. Deriving them from the
ENDF files directly means downloading the whole sublibrary, so this reads
OpenMC's ``half_life.json``, which is that derivation and is MIT licensed. The
provenance is recorded in the output so the indirection is visible rather than
implied.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from radiological_material_clearance_finder import nuclide as nuc  # noqa: E402

AME2020_URL = "https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt"
AME2020_MD5 = "6d28b75833cf53c7cc230223f63da6f6"
DATA = Path(__file__).resolve().parents[1] / "src" / "radiological_material_clearance_finder" / "data"

# The file documents its own Fortran format on header line 21. These are the
# zero-based slices that format implies.
_SLICE_Z = slice(9, 14)
_SLICE_A = slice(14, 19)
_SLICE_SYMBOL = slice(20, 23)
_SLICE_MASS_INT = slice(106, 109)
_SLICE_MASS_FRAC = slice(110, 123)
_HEADER_LINES = 36


def parse_ame2020(text: str) -> dict[str, float]:
    """Return {canonical nuclide name: relative atomic mass in u}.

    Estimated (non-experimental) masses are included. AME2020 flags them by
    writing ``#`` in place of the decimal point, which raises in ``float()``,
    so the point is restored before conversion.
    """
    masses: dict[str, float] = {}
    for line in text.splitlines()[_HEADER_LINES:]:
        if not line.strip():
            continue
        symbol = line[_SLICE_SYMBOL].strip()
        if symbol == "n":
            # The free neutron is not a nuclide anyone clears.
            continue
        mass_number = int(line[_SLICE_A])
        integer_part = line[_SLICE_MASS_INT].strip()
        fractional = line[_SLICE_MASS_FRAC].replace("#", ".").strip()
        mass_u = float(integer_part) + 1e-6 * float(fractional)
        name = nuc.normalise(f"{symbol}{mass_number}")
        masses[name] = mass_u
    return masses


def add_metastable_masses(masses: dict[str, float], half_lives: dict[str, float]) -> None:
    """Give each metastable state the mass of its ground state.

    AME2020 tabulates ground states. An isomer's excitation energy is at most a
    few MeV, which is under 1e-5 of the mass of even a light nuclide, so using
    the ground state mass is far below the precision that matters for a mass
    denominator. Without this, every metastable nuclide would be missing a mass
    and would be silently absent from the total mass of a material.
    """
    for name in half_lives:
        if "_m" not in name:
            continue
        ground = name.split("_m")[0]
        if ground in masses and name not in masses:
            masses[name] = masses[ground]


def load_half_lives(path: Path) -> dict[str, float]:
    """Return {canonical nuclide name: half-life in seconds} from OpenMC's JSON."""
    raw = json.loads(path.read_text())
    out: dict[str, float] = {}
    for key, seconds in raw.items():
        if seconds is None:
            continue
        try:
            name = nuc.normalise(key)
        except nuc.NuclideNameError:
            continue
        out[name] = float(seconds)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ame", type=Path, help="local mass_1.mas20.txt, else downloaded")
    parser.add_argument(
        "--half-life-json",
        type=Path,
        default=Path.home() / "openmc" / "openmc" / "data" / "half_life.json",
        help="OpenMC's half_life.json (ENDF/B-VIII.0 derived)",
    )
    args = parser.parse_args()

    if args.ame is not None:
        ame_text = args.ame.read_text()
        ame_origin = str(args.ame)
    else:
        request = urllib.request.Request(AME2020_URL, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(request, timeout=120) as handle:
            ame_text = handle.read().decode("utf-8", errors="replace")
        ame_origin = AME2020_URL

    masses = parse_ame2020(ame_text)

    if not args.half_life_json.exists():
        raise SystemExit(
            f"half-life source not found at {args.half_life_json}. Pass "
            f"--half-life-json pointing at an OpenMC checkout's "
            f"openmc/data/half_life.json"
        )
    half_lives = load_half_lives(args.half_life_json)
    add_metastable_masses(masses, half_lives)

    _write(
        DATA / "atomic_mass.json",
        {
            "_source": "AME2020 atomic mass evaluation",
            "_url": AME2020_URL,
            "_origin": ame_origin,
            "_reference": "Wang et al., Chinese Physics C 45, 030003 (2021), doi:10.1088/1674-1137/abddaf",
            "_note": "Relative atomic masses in u. Metastable states take their ground state mass, which is accurate to better than 1e-5.",
            "values": dict(sorted(masses.items())),
        },
    )
    _write(
        DATA / "half_life.json",
        {
            "_source": "ENDF/B-VIII.0 decay sublibrary, via OpenMC openmc/data/half_life.json (MIT)",
            "_url": "https://www.nndc.bnl.gov/endf-b8.0/download.html",
            "_origin": str(args.half_life_json),
            "_note": "Half-lives in seconds. Nuclides absent from this table are treated as stable.",
            "values": dict(sorted(half_lives.items())),
        },
    )
    print(f"atomic_mass.json: {len(masses)} nuclides")
    print(f"half_life.json:   {len(half_lives)} nuclides")


def _write(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=1, sort_keys=False) + "\n")


if __name__ == "__main__":
    main()
