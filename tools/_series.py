"""Derive the members of a natural decay series from IAEA decay mode data.

EU 2013/59 Annex VII Table A Part 2 and IAEA GSR Part 3 Table I.3 give one
value for a whole natural decay series rather than per nuclide, so applying
those values needs to know which nuclides are in the series. Rather than
hardcode a list, this walks the decay modes published by the IAEA Nuclear Data
Section, the same table ``build_decay_modes.py`` reads, so the membership is
derived from data with provenance and is regenerated with everything else.

Only alpha, beta and electron capture branches are followed, which is all a
natural series contains. Mass number never increases along any of them, so the
walk is bounded by construction.
"""
from __future__ import annotations

import csv
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from radiological_material_clearance_finder import nuclide as nuc  # noqa: E402

#: How each decay mode changes (proton number, neutron number).
_TRANSITIONS = {
    "A": (-2, -2),
    "B-": (1, -1),
    "B+": (-1, 1),
    "EC": (-1, 1),
    "2B-": (2, -2),
    "EC+B+": (-1, 1),
}

#: Branches rarer than this are not part of the series in any practical sense.
_MIN_BRANCH_PERCENT = 0.01


def _load(text: str) -> dict[tuple[int, int], dict]:
    """Index the Livechart ground state table by (proton number, neutron number)."""
    table: dict[tuple[int, int], dict] = {}
    for row in csv.DictReader(io.StringIO(text)):
        try:
            key = (int(row["z"]), int(row["n"]))
        except (KeyError, ValueError):
            continue
        table[key] = row
    return table


def series_members(text: str, head: str) -> list[str]:
    """Return every radionuclide reachable by decay from ``head``.

    Args:
        text: The Livechart ground state CSV.
        head: The series head, for example ``U238``.

    Returns:
        Canonical nuclide names, including the head, in sorted order. Stable end
        products are not included, since a limit on activity concentration
        cannot apply to them.
    """
    table = _load(text)
    by_symbol = {
        (row["symbol"].strip(), z + n): (z, n)
        for (z, n), row in table.items()
        for z, n in [(z, n)]
        if row.get("symbol")
    }
    symbol, mass, _ = nuc.parse(head)
    start = by_symbol.get((symbol, mass))
    if start is None:
        raise ValueError(f"{head} is not in the decay mode table")

    members: set[str] = set()
    stack = [start]
    seen: set[tuple[int, int]] = set()
    while stack:
        key = stack.pop()
        if key in seen:
            continue
        seen.add(key)
        row = table.get(key)
        if row is None:
            continue

        half_life = (row.get("half_life_sec") or "").strip()
        try:
            radioactive = float(half_life) > 0.0
        except ValueError:
            radioactive = False
        if not radioactive:
            # A stable end product carries no activity, so no limit applies.
            continue

        name = nuc.normalise(f"{row['symbol'].strip()}{key[0] + key[1]}")
        members.add(name)

        for index in (1, 2, 3):
            mode = (row.get(f"decay_{index}") or "").strip().upper()
            shift = _TRANSITIONS.get(mode)
            if shift is None:
                continue
            raw = (row.get(f"decay_{index}_%") or "").strip()
            try:
                if float(raw) < _MIN_BRANCH_PERCENT:
                    continue
            except ValueError:
                # A branch known to exist with no measured intensity still counts.
                pass
            stack.append((key[0] + shift[0], key[1] + shift[1]))

    # The Livechart ground state table carries no isomers, but a series can run
    # through one: Th-234 decays overwhelmingly to Pa-234m rather than to the
    # Pa-234 ground state, and the regulations list Pa-234m as a progeny in
    # their own tables. So isomers of members are included.
    #
    # This is deliberately over-inclusive. The source records which nuclide a
    # decay produces but not which level of it, so there is no way from this
    # data to tell that Pb-210 populates the Bi-210 ground state rather than
    # Bi-210m. Isomers such as Bi-210m, Bi-212m and Po-212m are therefore
    # included without positive evidence that the chain feeds them.
    #
    # The direction of that error is the safe one. Including an isomer applies
    # the series limit to it; excluding it would leave it with no limit at all
    # in a set that has no catch-all, so it would carry activity that the index
    # ignores. Over-inclusion is stricter, under-inclusion is a silent gap.
    from radiological_material_clearance_finder.decay import default_decay_data

    data = default_decay_data()
    for name in sorted(members):
        for level in (1, 2):
            isomer = f"{name}_m{level}"
            if data.knows(isomer) and data.half_life(isomer) is not None:
                members.add(isomer)
    return sorted(members)


