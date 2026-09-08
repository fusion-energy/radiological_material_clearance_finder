"""Clearance indexes for a radiological material.

Give it a nuclide inventory and it tells you whether the material meets the
clearance, exemption or disposal limits of the UK, German, US, EU and IAEA
regulations, and by how much:

    >>> from radiological_material_clearance_finder import Material, clearance_index
    >>> steel = Material({"Fe56": 8.4e22, "Co60": 2.1e9, "Cs137": 3.1e8})
    >>> result = clearance_index(steel, "UK_EPR16_out_of_scope")
    >>> result.clearable
    False
    >>> result.dominant(2)
    [('Co60', ...), ('Cs137', ...)]

Every regulation here uses the same arithmetic, a sum of activity-to-limit
ratios that must stay below one. What differs is the tables, and those are data:
see :func:`limit_sets` for what is available and
:func:`~radiological_material_clearance_finder.limits.get_limit_set` for the
provenance of any one of them.

OpenMC is not a dependency. To start from an OpenMC material, use
:func:`~radiological_material_clearance_finder.openmc_interop.from_openmc_material`,
which imports OpenMC only when called.
"""
from __future__ import annotations

from .classify import (
    UKWasteCategory,
    alpha_activity,
    beta_gamma_activity,
    nrc_waste_class,
    uk_waste_category,
)
from .cooling import IngrowthError, index_series, time_to_clear
from .decay import (
    DecayData,
    UnknownNuclideError,
    alpha_fraction,
    atomic_mass,
    decay_constant,
    half_life,
)
from .index import ClearanceResult, clearable_routes, clearance_index, clearance_indices
from .limits import LimitSet, get_limit_set, limit_sets, register_limit_set
from .material import ACTIVITY_UNITS, InsufficientDataError, Material
from .nuclide import NuclideNameError, normalise

__version__ = "0.1.0"

__all__ = [
    "Material",
    "ACTIVITY_UNITS",
    "clearance_index",
    "clearance_indices",
    "clearable_routes",
    "ClearanceResult",
    "limit_sets",
    "get_limit_set",
    "register_limit_set",
    "LimitSet",
    "nrc_waste_class",
    "uk_waste_category",
    "UKWasteCategory",
    "alpha_activity",
    "beta_gamma_activity",
    "time_to_clear",
    "index_series",
    "IngrowthError",
    "half_life",
    "decay_constant",
    "atomic_mass",
    "alpha_fraction",
    "DecayData",
    "normalise",
    "NuclideNameError",
    "UnknownNuclideError",
    "InsufficientDataError",
    "__version__",
]
