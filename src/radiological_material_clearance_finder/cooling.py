"""Finding when a material becomes clearable.

This package computes no decay of its own. The caller supplies a series of
materials at increasing cooling times, which is what a depletion calculation
already produces, and this finds where the clearance index crosses the
threshold.
"""
from __future__ import annotations

import math
from typing import Mapping

from .index import clearance_index
from .limits import get_limit_set
from .material import Material

__all__ = ["time_to_clear", "index_series", "IngrowthError"]


class IngrowthError(ValueError):
    """Raised when the index climbs back above the threshold after clearing.

    An index that rises with cooling time means a daughter is growing in faster
    than its parent decays, so a material can meet the limits at one time and
    fail them later. Reporting only the first crossing would be actively
    misleading, so this is raised instead.
    """


def index_series(
    series: Mapping[float, Material],
    limit_set,
    **kwargs,
) -> dict[float, float]:
    """Evaluate the clearance index at each cooling time.

    Args:
        series: Cooling time in seconds to the material at that time.
        limit_set: A registered limit set name, or a
            :class:`~radiological_material_clearance_finder.limits.LimitSet`.
        **kwargs: Passed through to
            :func:`~radiological_material_clearance_finder.index.clearance_index`.

    Returns:
        Cooling time to index, ordered by time.
    """
    return {
        time: clearance_index(material, limit_set, **kwargs).index
        for time, material in sorted(series.items())
    }


def time_to_clear(
    series: Mapping[float, Material],
    limit_set,
    *,
    allow_ingrowth: bool = False,
    **kwargs,
) -> float | None:
    """Find the cooling time at which a material first meets a set of limits.

    The index falls roughly exponentially with cooling time, so this
    interpolates logarithmically in the index between the two samples that
    bracket the threshold. Interpolating linearly instead would place the
    crossing systematically late, by a factor that grows with the spacing of the
    samples.

    Args:
        series: Cooling time in seconds to the material at that time. Times need
            not be sorted or evenly spaced.
        limit_set: A registered limit set name, or a
            :class:`~radiological_material_clearance_finder.limits.LimitSet`.
        allow_ingrowth: Return the first crossing even when the index later
            climbs back above the threshold, rather than raising
            :class:`IngrowthError`.
        **kwargs: Passed through to
            :func:`~radiological_material_clearance_finder.index.clearance_index`.

    Returns:
        The cooling time in seconds at which the index first falls below the
        threshold, ``0.0`` if the first sample already meets it, or ``None`` if
        the series never does.

    Raises:
        ValueError: If fewer than two times are given.
        IngrowthError: If the index climbs back above the threshold at a later
            time, so the material does not stay clear.
    """
    threshold = get_limit_set(limit_set).threshold
    indexes = index_series(series, limit_set, **kwargs)
    if len(indexes) < 2:
        raise ValueError("a cooling series needs at least two times to interpolate between")

    times = list(indexes)
    values = [indexes[t] for t in times]

    crossing = None
    if values[0] < threshold:
        crossing = times[0]
    else:
        for position in range(1, len(times)):
            if values[position] >= threshold:
                continue
            crossing = _interpolate(
                times[position - 1],
                values[position - 1],
                times[position],
                values[position],
                threshold,
            )
            break

    if crossing is None:
        return None

    if not allow_ingrowth:
        rebound = [
            (time, value)
            for time, value in zip(times, values)
            if time > crossing and value >= threshold
        ]
        if rebound:
            time, value = rebound[0]
            raise IngrowthError(
                f"the index first falls below {threshold:g} at {crossing:.4g} s but "
                f"climbs back to {value:.4g} at {time:.4g} s, so the material does "
                f"not stay clear. Pass allow_ingrowth=True to take the first "
                f"crossing anyway."
            )
    return crossing


def _interpolate(
    time_a: float, index_a: float, time_b: float, index_b: float, threshold: float
) -> float:
    """Interpolate the crossing time, logarithmically in the index."""
    if index_a <= 0.0 or index_b <= 0.0 or index_a == index_b:
        # Fall back to linear when a log is not defined or the index is flat.
        if index_a == index_b:
            return time_b
        fraction = (index_a - threshold) / (index_a - index_b)
    else:
        fraction = (math.log(index_a) - math.log(threshold)) / (
            math.log(index_a) - math.log(index_b)
        )
    fraction = min(max(fraction, 0.0), 1.0)
    return time_a + fraction * (time_b - time_a)
