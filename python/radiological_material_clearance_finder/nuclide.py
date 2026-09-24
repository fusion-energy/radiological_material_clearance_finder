"""Nuclide name parsing, validation and normalisation.

The canonical form used throughout this package is the one OpenMC and GND use,
``Co60`` and ``Ag108_m1``, because the most common source of an inventory is an
OpenMC material and matching it avoids a translation layer at the boundary.

Regulatory tables spell nuclides differently again (``Co-60``, ``Ag-108m``,
``U-238sec``, ``Sr-90+``), so `parse_regulatory` handles those and
separates the secular-equilibrium marker from the nuclide identity.
"""
from ._core import (
    NuclideNameError,
    atomic_number,
    element,
    is_valid,
    mass_number,
    metastable_state,
    normalise,
    parse,
    parse_regulatory,
)

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
