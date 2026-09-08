"""Nuclide name parsing, validation and normalisation.

The canonical form used throughout this package is the one OpenMC and GND use,
``Co60`` and ``Ag108_m1``, because the most common source of an inventory is an
OpenMC material and matching it avoids a translation layer at the boundary.

Regulatory tables spell nuclides differently again (``Co-60``, ``Ag-108m``,
``U-238sec``, ``Sr-90+``), so :func:`parse_regulatory` handles those and
separates the secular-equilibrium marker from the nuclide identity.
"""
from __future__ import annotations

import re

__all__ = [
    "normalise",
    "parse",
    "parse_regulatory",
    "element",
    "mass_number",
    "metastable_state",
    "atomic_number",
    "is_valid",
    "NuclideNameError",
]

#: Element symbols indexed by atomic number, with index 0 unused.
_SYMBOLS = (
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
    "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
    "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl",
    "Mc", "Lv", "Ts", "Og",
)

_Z_BY_SYMBOL = {sym: z for z, sym in enumerate(_SYMBOLS) if sym}
#: Case-insensitive lookup, so "co60", "CO-60" and "Co60" all resolve.
_Z_BY_LOWER = {sym.lower(): z for sym, z in _Z_BY_SYMBOL.items()}

# Element symbol, optional separator, mass number, optional metastable marker.
# The marker covers "m", "m1", "m2", "_m1", "_m2", the "n" that the German
# regulation uses for the second metastable state ("Hf-178n" is Hf178_m2), and
# the space-separated form EUR-Lex uses ("Zn-69 m").
_PATTERN = re.compile(
    r"""^\s*
    (?P<symbol>[A-Za-z]{1,2})
    [\s\-]?
    (?P<mass>\d{1,3})
    (?:[\s_]?(?P<meta>[mn])(?P<level>\d)?)?
    \s*$""",
    re.VERBOSE,
)

# Trailing regulatory markers: "+" means the tabulated daughters are covered,
# "sec" means the whole decay chain is. Either may be preceded by a space.
_MARKER = re.compile(r"[\s]*(?P<marker>\+|sec)\s*$", re.IGNORECASE)


class NuclideNameError(ValueError):
    """Raised when a string cannot be read as a single nuclide."""


def parse(name: str) -> tuple[str, int, int]:
    """Split a nuclide name into element symbol, mass number and metastable level.

    Args:
        name: A nuclide name in any accepted spelling, such as ``Co60``,
            ``Co-60``, ``Ag-108m`` or ``Ag108_m1``.

    Returns:
        ``(symbol, mass_number, metastable_level)``, where the level is 0 for a
        ground state and 1 or 2 for a metastable state. The symbol is returned
        in its canonical capitalisation regardless of the input's.

    Raises:
        NuclideNameError: If the name is not a single nuclide. A bare element
            such as ``Fe`` raises, because a clearance limit applies to a
            nuclide and silently guessing a mass number would be wrong.
    """
    if not isinstance(name, str):
        raise NuclideNameError(f"nuclide name must be a string, got {type(name).__name__}")

    match = _PATTERN.match(name)
    if match is None:
        if re.fullmatch(r"\s*[A-Za-z]{1,2}\s*", name or ""):
            raise NuclideNameError(
                f"{name!r} is an element, not a nuclide. Clearance limits are "
                f"per nuclide, so give a mass number, for example "
                f"{name.strip().capitalize()}56."
            )
        raise NuclideNameError(f"cannot read {name!r} as a nuclide name")

    symbol_in = match.group("symbol")
    z = _Z_BY_LOWER.get(symbol_in.lower())
    if z is None:
        raise NuclideNameError(f"unknown element symbol {symbol_in!r} in {name!r}")

    mass = int(match.group("mass"))
    # The parser validates syntax, not the drip lines, but these two bounds
    # catch the common typo class without pretending to know nuclear stability.
    # The heaviest nuclide in AME2020 is Og295.
    if mass > 300:
        raise NuclideNameError(
            f"{name!r} has mass number {mass}, above the heaviest known nuclide"
        )
    if mass < z:
        raise NuclideNameError(
            f"{name!r} has mass number {mass} below its proton number {z}, "
            f"which is not a real nuclide"
        )

    level = 0
    if match.group("meta"):
        if match.group("level"):
            level = int(match.group("level"))
        else:
            # "n" is the second metastable state, "m" the first.
            level = 2 if match.group("meta").lower() == "n" else 1
    return _SYMBOLS[z], mass, level


def normalise(name: str) -> str:
    """Return a nuclide name in canonical form.

    ``Co-60`` and ``co 60`` both become ``Co60``; ``Ag-108m`` becomes
    ``Ag108_m1`` and ``Hf-178n`` becomes ``Hf178_m2``. A trailing ``+`` is stripped, since it marks a secular
    equilibrium value rather than a different nuclide, but ``sec`` is not
    accepted here because a ``sec`` value is a separate table entry. Use
    :func:`parse_regulatory` when reading a regulatory table.

    Args:
        name: A nuclide name in any accepted spelling.

    Returns:
        The canonical name, such as ``Co60`` or ``Ag108_m1``.

    Raises:
        NuclideNameError: If the name is not a single nuclide.
    """
    if isinstance(name, str):
        name = re.sub(r"\s*\+\s*$", "", name)
    symbol, mass, level = parse(name)
    return f"{symbol}{mass}" + (f"_m{level}" if level else "")


def parse_regulatory(label: str) -> tuple[str, str | None]:
    """Read a nuclide name as spelled in a regulatory table.

    Regulatory tables mark parent nuclides whose limit already accounts for
    daughters in secular equilibrium, with ``+`` for the daughters tabulated
    alongside and ``sec`` for the whole decay chain. The marker is not part of
    the nuclide's identity but it does select a different limit value, so it is
    returned separately rather than discarded.

    Args:
        label: A table cell such as ``Sr-90+``, ``U-238sec`` or ``Th-232 sec``.

    Returns:
        ``(canonical_name, marker)`` where marker is ``"+"``, ``"sec"`` or
        ``None``.

    Raises:
        NuclideNameError: If the label is not a single nuclide.
    """
    marker = None
    text = label
    match = _MARKER.search(text)
    if match is not None:
        marker = match.group("marker").lower()
        text = text[: match.start()]
    symbol, mass, level = parse(text)
    return f"{symbol}{mass}" + (f"_m{level}" if level else ""), marker


def element(name: str) -> str:
    """Return the element symbol of a nuclide, for example ``Co`` for ``Co60``."""
    return parse(name)[0]


def mass_number(name: str) -> int:
    """Return the mass number of a nuclide, for example 60 for ``Co60``."""
    return parse(name)[1]


def metastable_state(name: str) -> int:
    """Return the metastable level of a nuclide, 0 for a ground state."""
    return parse(name)[2]


def atomic_number(name: str) -> int:
    """Return the proton number of a nuclide, for example 27 for ``Co60``."""
    return _Z_BY_SYMBOL[parse(name)[0]]


def is_valid(name: str) -> bool:
    """Whether a string can be read as a single nuclide."""
    try:
        parse(name)
    except NuclideNameError:
        return False
    return True
