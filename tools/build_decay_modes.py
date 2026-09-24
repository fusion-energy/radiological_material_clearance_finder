"""Build the alpha decay fraction table from the IAEA Livechart API.

    python tools/build_decay_modes.py

The UK waste categories limit alpha activity and beta/gamma activity separately,
so splitting a material's activity needs to know what fraction of each nuclide's
decays are alpha. This fetches the IAEA Nuclear Data Section ground state table,
which gives up to three decay modes per nuclide with their branching
percentages, and keeps the alpha branch.

Only ground states are published by this endpoint. Metastable states decay by
isomeric transition or beta and essentially never by alpha, so a nuclide absent
from this table is treated as a beta or gamma emitter, which is recorded in the
output rather than left implicit.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import urllib.request
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from radiological_material_clearance_finder import nuclide as nuc  # noqa: E402

URL = "https://nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=all"
DATA = Path(__file__).resolve().parents[1] / "crates" / "radiological-material-clearance-finder" / "data"


def fetch(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=180) as handle:
        return handle.read().decode("utf-8", errors="replace")


def parse(text: str) -> dict[str, float]:
    """Return {nuclide: fraction of decays that are alpha}, omitting zeroes."""
    fractions: dict[str, float] = {}
    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        symbol = (row.get("symbol") or "").strip()
        if not symbol or symbol == "n":
            continue
        try:
            mass = int(row["z"]) + int(row["n"])
            name = nuc.normalise(f"{symbol}{mass}")
        except (ValueError, KeyError, nuc.NuclideNameError):
            continue
        # An alpha branch with no measured intensity gets only the intensity the
        # measured branches leave unaccounted for. Giving it a flat 100 percent
        # made nuclides that are almost entirely beta emitters look like pure
        # alpha emitters, which would put their activity in the wrong half of
        # the UK alpha and beta or gamma split.
        alpha = 0.0
        measured = 0.0
        unmeasured_alpha = False
        unmeasured_other = False
        for index in (1, 2, 3):
            mode = (row.get(f"decay_{index}") or "").strip().upper()
            if not mode:
                continue
            raw = (row.get(f"decay_{index}_%") or "").strip()
            try:
                value = float(raw)
            except ValueError:
                if mode == "A":
                    unmeasured_alpha = True
                else:
                    unmeasured_other = True
                continue
            measured += value
            if mode == "A":
                alpha += value
        # The remainder only belongs to alpha when alpha is the only branch
        # whose intensity is unknown. With another unmeasured branch competing
        # for it, handing the lot to alpha would overstate the alpha fraction
        # and put activity in the wrong half of the UK alpha and beta or gamma
        # split.
        if unmeasured_alpha and not unmeasured_other:
            alpha += max(0.0, 100.0 - measured)
        if alpha > 0.0:
            fractions[name] = min(alpha, 100.0) / 100.0
    return fractions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cached", type=Path)
    args = parser.parse_args()

    text = args.cached.read_text(errors="replace") if args.cached else fetch(URL)
    fractions = parse(text)
    if len(fractions) < 100:
        raise SystemExit(f"only {len(fractions)} alpha emitters found, expected several hundred")

    payload = {
        "_source": "IAEA Nuclear Data Section, Livechart ground state table",
        "_url": URL,
        "_retrieved": date.today().isoformat(),
        "_note": (
            "Fraction of decays proceeding by alpha emission, for nuclides with a "
            "non-zero alpha branch. Any nuclide absent from this table has no alpha "
            "branch and is counted as a beta or gamma emitter. Only ground states "
            "are published by this endpoint; metastable states decay by isomeric "
            "transition or beta and are therefore treated as beta or gamma."
        ),
        "values": dict(sorted(fractions.items())),
    }
    (DATA / "alpha_fraction.json").write_text(json.dumps(payload, indent=1) + "\n")
    print(f"alpha_fraction.json: {len(fractions)} nuclides with an alpha branch")
    for name in ("Am241", "Pu239", "U238", "Ra226", "Po210", "Bi212", "Sr90", "Cs137"):
        print(f"  {name:8} {fractions.get(name, 0.0):.4f}")


if __name__ == "__main__":
    main()
