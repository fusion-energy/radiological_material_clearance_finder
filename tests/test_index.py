"""The clearance index, and the regulatory behaviour around it."""
import pytest

from radiological_material_clearance_finder import (
    LimitSet,
    Material,
    clearable_routes,
    clearance_index,
    clearance_indices,
    get_limit_set,
    register_limit_set,
)


@pytest.fixture
def toy_set():
    """A small limit set, so expectations can be computed by hand."""
    return LimitSet(
        name="TOY",
        label="toy",
        units="Bq/g",
        limits={"Co60": 1.0, "Cs137": 10.0},
        secular_equilibrium={"Sr90": ("Y90",)},
    )


def test_a_nuclide_at_its_limit_gives_an_index_of_one(toy_set):
    material = Material.from_specific_activities({"Co60": 1.0})
    assert clearance_index(material, toy_set).index == pytest.approx(1.0)


def test_ratios_sum(toy_set):
    material = Material.from_specific_activities({"Co60": 0.5, "Cs137": 2.0})
    result = clearance_index(material, toy_set)
    assert result.index == pytest.approx(0.5 + 0.2)
    assert result.by_nuclide == {"Co60": pytest.approx(0.5), "Cs137": pytest.approx(0.2)}


def test_clearable_is_decided_against_the_threshold(toy_set):
    assert clearance_index(
        Material.from_specific_activities({"Co60": 0.99}), toy_set
    ).clearable
    assert not clearance_index(
        Material.from_specific_activities({"Co60": 1.01}), toy_set
    ).clearable


def test_nuclides_with_no_limit_are_reported_not_hidden(toy_set):
    material = Material.from_specific_activities({"Co60": 1.0, "Fe55": 1e6})
    result = clearance_index(material, toy_set)
    assert result.index == pytest.approx(1.0)
    assert result.uncovered == {"Fe55": 1e6}
    assert result.uncovered_fraction > 0.999


def test_a_catch_all_limit_covers_unlisted_nuclides():
    limits = LimitSet(
        name="TOY_DEFAULT", label="toy", units="Bq/g",
        limits={"Co60": 1.0}, default_limit=0.01,
    )
    material = Material.from_specific_activities({"Co60": 1.0, "Fe55": 0.02})
    result = clearance_index(material, limits)
    assert result.index == pytest.approx(1.0 + 2.0)
    assert result.defaulted == ("Fe55",)
    assert result.uncovered == {}


def test_the_catch_all_can_be_switched_off():
    limits = LimitSet(
        name="TOY_DEFAULT2", label="toy", units="Bq/g",
        limits={"Co60": 1.0}, default_limit=0.01,
    )
    material = Material.from_specific_activities({"Co60": 1.0, "Fe55": 0.02})
    result = clearance_index(material, limits, apply_default_limit=False)
    assert result.index == pytest.approx(1.0)
    assert result.uncovered == {"Fe55": 0.02}


def test_daughters_are_excluded_only_when_the_parent_is_present(toy_set):
    with_parent = Material.from_specific_activities({"Sr90": 1.0, "Y90": 1.0})
    without = Material.from_specific_activities({"Y90": 1.0})
    assert "Y90" in clearance_index(with_parent, toy_set).excluded
    assert clearance_index(without, toy_set).excluded == {}


def test_daughter_exclusion_can_be_switched_off():
    limits = LimitSet(
        name="TOY_SE", label="toy", units="Bq/g",
        limits={"Sr90": 1.0, "Y90": 1.0}, secular_equilibrium={"Sr90": ("Y90",)},
    )
    material = Material.from_specific_activities({"Sr90": 1.0, "Y90": 1.0})
    assert clearance_index(material, limits).index == pytest.approx(1.0)
    assert clearance_index(material, limits, exclude_daughters=False).index == pytest.approx(2.0)


def test_the_marked_limit_applies_only_when_daughters_are_present():
    """StrlSchV gives Th-232 as 10 Bq/g plain and 0.01 Bq/g with its chain."""
    limits = LimitSet(
        name="TOY_TH", label="toy", units="Bq/g",
        limits={"Th232": 10.0, "Ra228": 1.0},
        limits_secular_equilibrium={"Th232": 0.01},
        secular_equilibrium={"Th232": ("Ra228",)},
    )
    alone = Material.from_specific_activities({"Th232": 1.0})
    assert clearance_index(alone, limits).index == pytest.approx(0.1)

    with_chain = Material.from_specific_activities({"Th232": 1.0, "Ra228": 1.0})
    result = clearance_index(with_chain, limits)
    assert result.limits_used["Th232"] == 0.01
    assert result.index == pytest.approx(100.0)
    assert "Ra228" in result.excluded


def test_the_short_half_life_scope_rule_is_a_whole_material_test():
    """EPR 2016 puts a substance out of scope only if every radionuclide is short lived."""
    limits = LimitSet(
        name="TOY_SCOPE", label="toy", units="Bq/g",
        limits={"N16": 1.0, "Co60": 1.0}, min_half_life_scope=100.0,
    )
    short_only = Material.from_specific_activities({"N16": 1e6})
    assert clearance_index(short_only, limits).out_of_scope
    assert clearance_index(short_only, limits).clearable

    mixed = Material.from_specific_activities({"N16": 1e6, "Co60": 1e6})
    result = clearance_index(mixed, limits)
    assert not result.out_of_scope
    # The short lived nuclide still counts towards the sum.
    assert result.index == pytest.approx(2e6)


def test_explicitly_unlimited_nuclides_are_covered_not_uncovered():
    limits = LimitSet(
        name="TOY_UNLIM", label="toy", units="Bq/g",
        limits={"Co60": 1.0}, unlimited=("H3",),
    )
    material = Material.from_specific_activities({"Co60": 1.0, "H3": 1e9})
    result = clearance_index(material, limits)
    assert result.index == pytest.approx(1.0)
    assert result.unlimited == ("H3",)
    assert result.uncovered == {}


def test_unknown_limit_set_lists_the_known_ones():
    with pytest.raises(KeyError, match="Available"):
        get_limit_set("NOT_A_REAL_SET")


def test_result_serialises(toy_set):
    material = Material.from_specific_activities({"Co60": 2.0})
    payload = clearance_index(material, toy_set).to_dict()
    assert payload["index"] == pytest.approx(2.0)
    assert payload["clearable"] is False
    import json

    json.dumps(payload)


def test_clearance_indices_skips_sets_the_material_cannot_supply():
    """A material with no density must not hide the Bq/g results."""
    material = Material({"Fe56": 1e22, "Co60": 1e10})
    results = clearance_indices(material)
    assert "StrlSchV_unrestricted" in results
    assert "Fetter" not in results


def test_clearable_routes_are_ordered_by_margin():
    material = Material({"Fe56": 1e22, "Co60": 1e6})
    routes = clearable_routes(material)
    assert routes
    assert "StrlSchV_unrestricted" in routes


def test_register_limit_set_makes_it_available():
    register_limit_set(
        LimitSet(name="SITE_SPECIFIC", label="site", units="Bq/g", limits={"Co60": 1.0})
    )
    assert get_limit_set("SITE_SPECIFIC").limits["Co60"] == 1.0
