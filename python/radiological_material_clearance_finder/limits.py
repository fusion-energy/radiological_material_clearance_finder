"""Limit sets: the regulatory tables a material is measured against.

Every limit set is data, compiled into the extension from the JSON in
``crates/radiological-material-clearance-finder/data/limits``. Adding a
jurisdiction means adding a file and a build script, never a branch in the
calculation. Only two things in these regulations cannot be expressed as a
table, and both are named explicitly rather than worked around:

* limits given per gram of material (the NRC's nCi/g transuranic entries), which
  need the material's density before they can be compared against Ci/m3, and
* limits that depend on the material's own nuclides (the NRC Class A rule that
  anything with a half-life under five years takes a 700 Ci/m3 limit), named in
  `DYNAMIC_RULES`.
"""
from ._core import DYNAMIC_RULES, LimitSet, get_limit_set, limit_sets, register_limit_set

__all__ = ["LimitSet", "limit_sets", "get_limit_set", "register_limit_set", "DYNAMIC_RULES"]
