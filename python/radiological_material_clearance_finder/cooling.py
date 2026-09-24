"""Finding when a material becomes clearable.

This package computes no decay of its own. The caller supplies a series of
materials at increasing cooling times, which is what a depletion calculation
already produces, and this finds where the clearance index crosses the
threshold.
"""
from ._core import IngrowthError, index_series, time_to_clear

__all__ = ["time_to_clear", "index_series", "IngrowthError"]
