"""Decay and mass data, checked against values that are known independently."""
import math

import pytest

from radiological_material_clearance_finder import decay


def test_half_lives_match_published_values():
    year = 365.25 * 86400
    assert decay.half_life("Co60") / year == pytest.approx(5.27, rel=1e-3)
    assert decay.half_life("Cs137") / year == pytest.approx(30.08, rel=1e-3)
    assert decay.half_life("H3") / year == pytest.approx(12.32, rel=1e-3)
    assert decay.half_life("Sr90") / year == pytest.approx(28.9, rel=1e-2)


def test_stable_nuclides_have_no_half_life():
    assert decay.half_life("Fe56") is None
    assert decay.decay_constant("Fe56") == 0.0
    assert not decay.is_radioactive("Fe56")


def test_decay_constant_follows_from_the_half_life():
    assert decay.decay_constant("Co60") == pytest.approx(
        math.log(2) / decay.half_life("Co60")
    )


def test_atomic_masses_match_ame2020():
    assert decay.atomic_mass("H1") == pytest.approx(1.007825031898, rel=1e-11)
    assert decay.atomic_mass("Co60") == pytest.approx(59.933815536, rel=1e-11)
    assert decay.atomic_mass("U238") == pytest.approx(238.050786936, rel=1e-11)


def test_metastable_states_take_the_ground_state_mass():
    assert decay.atomic_mass("Ag108_m1") == pytest.approx(decay.atomic_mass("Ag108"))


def test_specific_activity_of_pure_cobalt_60():
    """The textbook value is 1131 Ci/g."""
    per_gram = decay.decay_constant("Co60") * decay.AVOGADRO / decay.atomic_mass("Co60")
    assert per_gram / decay.BECQUEREL_PER_CURIE == pytest.approx(1131, rel=1e-3)


def test_alpha_fractions():
    assert decay.alpha_fraction("Am241") == pytest.approx(1.0)
    assert decay.alpha_fraction("Pu239") == pytest.approx(1.0)
    # Bi-212 branches 35.94 percent alpha and the rest beta.
    assert decay.alpha_fraction("Bi212") == pytest.approx(0.3594, abs=1e-4)
    assert decay.alpha_fraction("Cs137") == 0.0
    assert decay.alpha_fraction("Ag108_m1") == 0.0


def test_unknown_nuclide_raises_rather_than_being_assumed_stable():
    with pytest.raises(decay.UnknownNuclideError):
        decay.half_life("Og296")


def test_overrides_replace_half_lives():
    data = decay.DecayData().with_overrides({"Co-60": 1.0})
    assert data.half_life("Co60") == 1.0
    assert data.half_life("Cs137") == decay.half_life("Cs137")
    assert data.alpha_fraction("Am241") == pytest.approx(1.0)


def test_override_to_none_marks_a_nuclide_stable():
    data = decay.DecayData().with_overrides({"Co60": None})
    assert data.half_life("Co60") is None
