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
        # Sr90 needs a limit of its own: a parent with no limit cannot account
        # for its daughters, which is a separate case tested below.
        limits={"Co60": 1.0, "Cs137": 10.0, "Sr90": 1.0, "Y90": 1.0},
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


def test_a_parent_with_no_limit_cannot_account_for_its_daughters():
    """A parent that is not limited here contributes nothing, so it can credit nothing.

    Crediting against a limit that does not exist removes the daughter from the
    sum and drives the index towards zero, reporting clearable. StrlSchV_soil
    lists Pa-233 but not its parent Np-237, so this is real shipped data.
    """
    material = Material.from_specific_activities({"Np237": 1e4, "Pa233": 1e4})
    result = clearance_index(material, "StrlSchV_soil")
    assert result.excluded == {}
    assert result.index == pytest.approx(
        clearance_index(material, "StrlSchV_soil", exclude_daughters=False).index
    )
    assert not result.clearable
    assert "Np237" in result.uncovered


def test_a_parent_only_accounts_for_the_daughter_activity_it_supports():
    """Secular equilibrium means equal activities, so that is all a parent covers."""
    limits = LimitSet(
        name="TOY_CREDIT", label="toy", units="Bq/g",
        limits={"Sr90": 1.0, "Y90": 1.0}, secular_equilibrium={"Sr90": ("Y90",)},
    )
    # In equilibrium the daughter is fully accounted for.
    equilibrium = Material.from_specific_activities({"Sr90": 100.0, "Y90": 100.0})
    result = clearance_index(equilibrium, limits)
    assert result.excluded == {"Y90": result.excluded["Y90"]}
    assert result.index == pytest.approx(100.0)

    # A trace of parent cannot account for a large daughter.
    lopsided = Material.from_specific_activities({"Sr90": 1.0, "Y90": 1000.0})
    result = clearance_index(lopsided, limits)
    assert result.excluded == {}
    assert result.credited == {"Y90": pytest.approx(1.0)}
    # 1.0 of Sr90, plus the 999.0 of Y90 its parent cannot support.
    assert result.index == pytest.approx(1.0 + 999.0)


def test_a_sec_only_limit_is_applied_rather_than_left_unreachable():
    """UK_IRR17_natural publishes U-238 only as a whole-chain "sec" row."""
    result = clearance_index(
        Material.from_specific_activities({"U238": 1e6}), "UK_IRR17_natural"
    )
    assert result.limits_used["U238"] == 1.0
    assert result.uncovered == {}
    assert result.index == pytest.approx(1e6)


def test_the_scope_rule_does_not_fire_when_there_are_no_radionuclides():
    """It is a rule about half-lives, so with no radionuclides it cannot apply.

    Returning out of scope here would force clearable regardless of the index,
    which is reachable when activities are supplied directly for a nuclide the
    decay tables call stable.
    """
    limits = LimitSet(
        name="TOY_SCOPE2", label="toy", units="Bq/g",
        limits={"Fe56": 1.0}, min_half_life_scope=100.0,
    )
    result = clearance_index(Material.from_specific_activities({"Fe56": 1e6}), limits)
    assert not result.out_of_scope
    assert not result.clearable
    assert result.index == pytest.approx(1e6)


def test_metal_overrides_replace_and_add_limits():
    """The NRC lists activated metal as its own rows, so an override may add one."""
    limits = LimitSet(
        name="TOY_METAL", label="toy", units="Bq/g",
        limits={"C14": 8.0},
        metal_overrides={"C14": 80.0, "Ni59": 220.0},
    )
    material = Material.from_specific_activities({"C14": 8.0, "Ni59": 220.0})

    plain = clearance_index(material, limits, metal=False)
    assert plain.index == pytest.approx(1.0)
    assert plain.uncovered == {"Ni59": 220.0}, "Ni59 has no non-metal limit"

    activated = clearance_index(material, limits, metal=True)
    assert activated.index == pytest.approx(0.1 + 1.0)
    assert activated.uncovered == {}


def test_the_nrc_short_lived_rule_limits_by_half_life():
    """Table 2 gives everything under five years a 700 Ci/m3 Class A limit."""
    limits = LimitSet(
        name="TOY_DYNAMIC", label="toy", units="Bq/g",
        limits={}, dynamic_rule="nrc_short_lived_class_a",
    )
    # Mn-54 is 312 days, Cs-137 is 30 years.
    result = clearance_index(
        Material.from_specific_activities({"Mn54": 700.0, "Cs137": 1e6}), limits
    )
    assert result.limits_used == {"Mn54": 700.0}
    assert result.index == pytest.approx(1.0)
    assert result.uncovered == {"Cs137": 1e6}


def test_the_threshold_boundary_is_exclusive():
    """An index of exactly the threshold does not clear."""
    limits = LimitSet(name="TOY_EDGE", label="toy", units="Bq/g", limits={"Co60": 1.0})
    assert not clearance_index(
        Material.from_specific_activities({"Co60": 1.0}), limits
    ).clearable
    assert clearance_index(
        Material.from_specific_activities({"Co60": 0.999}), limits
    ).clearable


def test_clearable_routes_are_ordered_by_margin():
    material = Material({"Fe56": 1e22, "Co60": 1e6})
    routes = clearable_routes(material)
    assert routes
    margins = [
        clearance_index(material, name).index / get_limit_set(name).threshold
        for name in routes
    ]
    assert margins == sorted(margins), "routes are not ordered by margin"


def test_a_limit_set_built_in_python_canonicalises_its_keys():
    """Otherwise a hand-written "Co-60" would silently match nothing."""
    limits = LimitSet(
        name="TOY_CANON", label="toy", units="Bq/g",
        limits={"Co-60": 1.0}, secular_equilibrium={"Sr-90": ("Y-90",)},
    )
    assert limits.limits == {"Co60": 1.0}
    assert limits.secular_equilibrium == {"Sr90": ("Y90",)}
    result = clearance_index(Material.from_specific_activities({"Co60": 2.0}), limits)
    assert result.index == pytest.approx(2.0)


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"units": "sieverts"}, "not one of"),
        ({"threshold": 0.0}, "threshold must be positive"),
        ({"limits": {"Co60": 0.0}}, "must be positive"),
        ({"default_limit": -1.0}, "default_limit must be positive"),
    ],
)
def test_limit_sets_reject_values_that_cannot_mean_anything(kwargs, match):
    base = {"name": "TOY_BAD", "label": "toy", "units": "Bq/g", "limits": {"Co60": 1.0}}
    with pytest.raises(ValueError, match=match):
        LimitSet(**{**base, **kwargs})


def test_a_lenient_equilibrium_limit_needs_real_equilibrium():
    """IRR 2017 gives U-240 0.01 Bq/g plain and 100 marked, so the marked row is
    a ten thousand fold relief. Buying that with a trace of daughter twelve
    orders of magnitude below equilibrium would be absurd, so it is withheld
    until the daughters are actually there.
    """
    alone = clearance_index(
        Material.from_specific_activities({"U240": 1.0}), "UK_IRR17_notification"
    )
    assert alone.index == pytest.approx(100.0)
    assert not alone.clearable

    trace = clearance_index(
        Material.from_specific_activities({"U240": 1.0, "Np240_m1": 1e-12}),
        "UK_IRR17_notification",
    )
    assert trace.index == pytest.approx(100.0), "a trace must not unlock the lenient row"
    assert not trace.clearable

    equilibrium = clearance_index(
        Material.from_specific_activities(
            {"U240": 1.0, "Np240_m1": 1.0, "Np240": 1.0}
        ),
        "UK_IRR17_notification",
    )
    assert equilibrium.limits_used["U240"] == 100.0
    assert equilibrium.clearable


def test_a_stricter_equilibrium_limit_applies_as_soon_as_daughters_appear():
    """The conservative direction needs no equilibrium test.

    StrlSchV gives Th-232 10 Bq/g plain and 0.01 marked.
    """
    result = clearance_index(
        Material.from_specific_activities({"Th232": 1.0, "Ra228": 1.0}),
        "StrlSchV_unrestricted",
    )
    assert result.limits_used["Th232"] == 0.01


def test_an_unlimited_parent_cannot_account_for_its_daughters():
    """No limit on the parent says nothing about the daughters.

    Fetter marks some nuclides as having no limit at all. Crediting daughters
    against a limit that does not exist deletes their activity outright.
    """
    limits = LimitSet(
        name="TOY_UNLIM_PARENT", label="toy", units="Bq/g",
        limits={"Ba137_m1": 1.0}, unlimited=("Cs137",),
        secular_equilibrium={"Cs137": ("Ba137_m1",)},
    )
    result = clearance_index(
        Material.from_specific_activities({"Cs137": 1e6, "Ba137_m1": 1e6}), limits
    )
    assert result.excluded == {}
    assert result.index == pytest.approx(1e6)
    assert not result.clearable


def test_a_secular_equilibrium_cycle_does_not_delete_both_nuclides():
    """A nuclide already accounted for cannot also account for others."""
    limits = LimitSet(
        name="TOY_CYCLE", label="toy", units="Bq/g",
        limits={"Cs137": 1.0, "Ba137_m1": 1.0},
        secular_equilibrium={"Cs137": ("Ba137_m1",), "Ba137_m1": ("Cs137",)},
    )
    result = clearance_index(
        Material.from_specific_activities({"Cs137": 1e6, "Ba137_m1": 1e6}), limits
    )
    assert len(result.excluded) == 1, "exactly one side of the pair may be credited"
    assert result.index == pytest.approx(1e6)
    assert not result.clearable


def test_a_hand_built_set_with_only_a_whole_chain_row_still_applies_it():
    """Promotion has to happen however the set was built, not only from JSON."""
    limits = LimitSet(name="TOY_SEC_ONLY", label="toy", units="Bq/g",
                      limits={"U238_sec": 1.0})
    result = clearance_index(Material.from_specific_activities({"U238": 5.0}), limits)
    assert result.index == pytest.approx(5.0)
    assert result.uncovered == {}


def test_two_keys_normalising_to_one_nuclide_is_an_error():
    """Schedule 7 publishes U-240 and U-240+ ten thousand apart.

    Whichever landed last in the mapping would silently win.
    """
    with pytest.raises(ValueError, match="both normalise to"):
        LimitSet(name="TOY_COLLIDE", label="toy", units="Bq/g",
                 limits={"U240": 1000.0, "U-240": 0.1})
    # The same value twice is not a conflict.
    LimitSet(name="TOY_SAME", label="toy", units="Bq/g",
             limits={"U240": 1000.0, "U-240": 1000.0})


@pytest.mark.parametrize("bad", [float("inf"), float("nan")])
def test_non_finite_limits_are_rejected(bad):
    """An infinite limit divides every activity to nothing, reading as no limit."""
    with pytest.raises(ValueError, match="positive and finite"):
        LimitSet(name="TOY_INF", label="toy", units="Bq/g", limits={"Co60": bad})
