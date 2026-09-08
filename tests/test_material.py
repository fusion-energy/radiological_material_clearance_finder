"""Material construction and the unit handling behind a specific activity."""
import pytest

from radiological_material_clearance_finder import decay
from radiological_material_clearance_finder.material import (
    InsufficientDataError,
    Material,
)


def test_specific_activity_is_scale_invariant():
    """Bq/g does not care whether the amounts are counts, densities or fractions."""
    small = Material({"Fe56": 1e22, "Co60": 1e12}).specific_activity()
    large = Material({"Fe56": 1e28, "Co60": 1e18}).specific_activity()
    assert small == pytest.approx(large, rel=1e-12)


def test_pure_cobalt_60_has_the_textbook_specific_activity():
    material = Material.from_masses({"Co60": 1.0})
    curies_per_gram = material.specific_activity() / decay.BECQUEREL_PER_CURIE
    assert curies_per_gram == pytest.approx(1131, rel=1e-3)


def test_atom_densities_determine_the_mass_density():
    """0.0849 atoms/barn-cm of Fe-56 is iron at about 7.9 g/cm3."""
    material = Material.from_atom_densities({"Fe56": 0.0849})
    assert material.density == pytest.approx(7.89, rel=1e-2)


def test_stable_nuclides_dilute_the_specific_activity():
    active_only = Material.from_masses({"Co60": 1.0}).specific_activity()
    diluted = Material.from_masses({"Co60": 1.0, "Fe56": 999.0}).specific_activity()
    assert diluted == pytest.approx(active_only / 1000.0, rel=1e-9)


def test_an_all_radioactive_inventory_warns_about_a_truncated_mass():
    with pytest.warns(UserWarning, match="every nuclide in this material is radioactive"):
        Material({"Co60": 1e12, "Cs137": 1e12})


def test_a_single_pure_source_does_not_warn():
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        Material({"Co60": 1e12})


def test_a_material_with_stable_nuclides_does_not_warn():
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        Material({"Fe56": 1e22, "Co60": 1e12})


def test_volumetric_activity_needs_a_density_and_says_so():
    material = Material({"Fe56": 1e22, "Co60": 1e12})
    with pytest.raises(InsufficientDataError, match="density"):
        material.activity("Ci/m3")


def test_specific_activities_can_be_given_directly():
    material = Material.from_specific_activities({"Co60": 5.0, "Cs137": 2.0})
    assert material.specific_activity() == pytest.approx(7.0)
    assert material.specific_activity(by_nuclide=True)["Co60"] == 5.0


def test_activity_units_convert_consistently():
    material = Material.from_atom_densities({"Fe56": 0.0849, "Co60": 1e-10}, volume=1000.0)
    per_gram = material.activity("Bq/g")
    assert material.activity("Bq/kg") == pytest.approx(per_gram * 1000.0)
    assert material.activity("Bq/cm3") == pytest.approx(per_gram * material.density)
    assert material.activity("Bq/m3") == pytest.approx(per_gram * material.density * 1e6)
    assert material.activity("Ci/m3") == pytest.approx(
        material.activity("Bq/m3") / decay.BECQUEREL_PER_CURIE
    )
    assert material.activity("Bq") == pytest.approx(per_gram * material.mass)


def test_unknown_units_are_rejected():
    with pytest.raises(ValueError, match="unknown activity units"):
        Material({"Co60": 1e12}).activity("banana")


def test_names_are_normalised_and_duplicates_combined():
    material = Material({"Co-60": 1e12, "co60": 1e12})
    assert material.nuclides == ("Co60",)
    assert material.atoms["Co60"] == 2e12


def test_negative_amounts_are_rejected():
    with pytest.raises(ValueError, match="negative"):
        Material({"Co60": -1.0})


def test_mass_fractions_and_masses_agree():
    by_mass = Material.from_masses({"Fe56": 980.0, "Co60": 20.0}).specific_activity()
    by_fraction = Material.from_mass_fractions({"Fe56": 0.98, "Co60": 0.02}).specific_activity()
    assert by_mass == pytest.approx(by_fraction, rel=1e-12)
