"""The clearance index itself: a sum of activity-to-limit ratios.

Every regulation here uses the same arithmetic. Each radionuclide's activity is
divided by its tabulated limit and the ratios are summed; a total below one
means the material meets the limits. The German regulation calls it the
Summenformel, the UK calls it the summation rule, the NRC calls it the sum of
fractions rule and Fetter calls the result a waste disposal rating. This module
implements it once.

What differs between regulations, and what `ClearanceResult` therefore
records, is what happens to a nuclide that is *not* simply looked up: one whose
parent already accounts for it, one the table does not list, and one the
regulation places outside its scope entirely.
"""
from ._core import (
    EQUILIBRIUM_TOLERANCE,
    ClearanceResult,
    clearable_routes,
    clearance_index,
    clearance_indices,
)

__all__ = [
    "clearance_index",
    "clearance_indices",
    "clearable_routes",
    "ClearanceResult",
    "EQUILIBRIUM_TOLERANCE",
]
