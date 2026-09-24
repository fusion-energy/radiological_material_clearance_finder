"""The material a clearance index is computed for.

A material here is just a nuclide inventory plus enough information to put it on
a per-gram or per-volume basis. It deliberately knows nothing about geometry,
temperature or cross sections, so it can be built from a dictionary typed by
hand as readily as from a depletion result.
"""
from ._core import ACTIVITY_UNITS, InsufficientDataError, Material

__all__ = ["Material", "ACTIVITY_UNITS", "InsufficientDataError"]
