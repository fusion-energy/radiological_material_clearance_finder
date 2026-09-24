"""Waste classification, where the answer is a category rather than an index.

The US and UK both sort waste into named classes rather than reporting a single
ratio. The NRC's classes come out of the same sum-of-fractions machinery as a
clearance index, applied to several tables and combined by rules in 10 CFR
61.55. The UK's categories instead compare gross alpha and beta or gamma
activity against concentration thresholds, so they need no per-nuclide table at
all and are computed directly.

The UK thresholds are stated per tonne in the 2007 policy and held here in
Bq/g: `UK_LLW_ALPHA_BQ_PER_G` (4 GBq/te), `UK_LLW_BETA_GAMMA_BQ_PER_G`
(12 GBq/te), `UK_HIGH_VOLUME_VLLW_BQ_PER_G` (4 MBq/te) and
`UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G` (40 MBq/te, a shared allowance
for tritium and carbon-14).
"""
from ._core import (
    UK_HIGH_VOLUME_VLLW_BQ_PER_G,
    UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G,
    UK_LLW_ALPHA_BQ_PER_G,
    UK_LLW_BETA_GAMMA_BQ_PER_G,
    UKWasteCategory,
    alpha_activity,
    beta_gamma_activity,
    nrc_waste_class,
    uk_waste_category,
)

__all__ = [
    "nrc_waste_class",
    "uk_waste_category",
    "UKWasteCategory",
    "alpha_activity",
    "beta_gamma_activity",
    "UK_LLW_ALPHA_BQ_PER_G",
    "UK_LLW_BETA_GAMMA_BQ_PER_G",
    "UK_HIGH_VOLUME_VLLW_BQ_PER_G",
    "UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G",
]
