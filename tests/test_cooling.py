"""Finding the cooling time at which a material clears."""
import math

import pytest

from radiological_material_clearance_finder import (
    IngrowthError,
    LimitSet,
    Material,
    index_series,
    time_to_clear,
)
from radiological_material_clearance_finder.decay import half_life

COBALT = LimitSet(name="COBALT", label="toy", units="Bq/g", limits={"Co60": 1.0})


def decayed(initial, seconds):
    """A material holding Co-60 decayed for a time, as a depletion would give."""
    factor = math.exp(-math.log(2) * seconds / half_life("Co60"))
    return Material.from_specific_activities({"Co60": initial * factor})


def test_index_series_evaluates_each_time():
    series = {0.0: decayed(10.0, 0.0), 1e8: decayed(10.0, 1e8)}
    indexes = index_series(series, COBALT)
    assert indexes[0.0] == pytest.approx(10.0)
    assert indexes[1e8] < indexes[0.0]


def test_the_crossing_time_matches_the_analytic_answer():
    """Ten times the limit takes log2(10) half-lives to reach it."""
    expected = half_life("Co60") * math.log2(10.0)
    times = [0.0, 5e8, 1e9, 2e9]
    series = {t: decayed(10.0, t) for t in times}
    found = time_to_clear(series, COBALT)
    assert found == pytest.approx(expected, rel=1e-9)


def test_log_interpolation_beats_linear_on_coarse_samples():
    """The index falls exponentially, so linear interpolation lands late."""
    expected = half_life("Co60") * math.log2(10.0)
    series = {t: decayed(10.0, t) for t in (0.0, 4e9)}
    found = time_to_clear(series, COBALT)
    assert found == pytest.approx(expected, rel=1e-9)

    index_a, index_b = 10.0, index_series(series, COBALT)[4e9]
    linear = 4e9 * (index_a - 1.0) / (index_a - index_b)
    assert linear > found


def test_an_already_clear_material_returns_the_first_time():
    series = {100.0: decayed(0.5, 0.0), 1e8: decayed(0.5, 1e8)}
    assert time_to_clear(series, COBALT) == 100.0


def test_a_series_that_never_clears_returns_none():
    series = {t: decayed(1e12, t) for t in (0.0, 1e6, 1e7)}
    assert time_to_clear(series, COBALT) is None


def test_ingrowth_after_clearing_raises():
    """Clearing at one time and failing later must not be reported as clear."""
    series = {
        0.0: Material.from_specific_activities({"Co60": 5.0}),
        1e8: Material.from_specific_activities({"Co60": 0.5}),
        # A daughter grows in and pushes the index back up.
        2e8: Material.from_specific_activities({"Co60": 3.0}),
        3e8: Material.from_specific_activities({"Co60": 0.2}),
    }
    with pytest.raises(IngrowthError, match="does not stay clear"):
        time_to_clear(series, COBALT)


def test_ingrowth_can_be_accepted_explicitly():
    series = {
        0.0: Material.from_specific_activities({"Co60": 5.0}),
        1e8: Material.from_specific_activities({"Co60": 0.5}),
        2e8: Material.from_specific_activities({"Co60": 3.0}),
    }
    found = time_to_clear(series, COBALT, allow_ingrowth=True)
    assert found == pytest.approx(_expected_crossing(0.0, 5.0, 1e8, 0.5))


def test_a_monotonic_series_does_not_raise():
    series = {t: decayed(10.0, t) for t in (0.0, 5e8, 1e9, 2e9)}
    assert time_to_clear(series, COBALT) is not None


def _expected_crossing(time_a, index_a, time_b, index_b, threshold=1.0):
    fraction = (math.log(index_a) - math.log(threshold)) / (
        math.log(index_a) - math.log(index_b)
    )
    return time_a + fraction * (time_b - time_a)


def test_a_series_needs_at_least_two_points():
    with pytest.raises(ValueError, match="at least two"):
        time_to_clear({0.0: decayed(10.0, 0.0)}, COBALT)


def test_the_series_need_not_be_given_in_order():
    """The docstring promises this, so it is pinned."""
    times = [0.0, 5e8, 1e9, 2e9]
    ascending = {t: decayed(10.0, t) for t in times}
    shuffled = {t: decayed(10.0, t) for t in [1e9, 0.0, 2e9, 5e8]}
    assert list(index_series(shuffled, COBALT)) == sorted(times)
    assert time_to_clear(shuffled, COBALT) == pytest.approx(
        time_to_clear(ascending, COBALT)
    )


def test_an_already_clear_series_returns_its_first_sample_not_zero():
    """The series says nothing before its earliest sample, so it cannot claim zero."""
    series = {100.0: decayed(0.5, 0.0), 1e8: decayed(0.5, 1e8)}
    assert time_to_clear(series, COBALT) == 100.0
    later = {5e7: decayed(0.5, 0.0), 1e8: decayed(0.5, 1e8)}
    assert time_to_clear(later, COBALT) == 5e7
