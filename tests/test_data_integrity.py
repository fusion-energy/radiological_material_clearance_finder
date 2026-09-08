"""Checks on the shipped regulatory tables.

These guard the failure mode that matters most: a table that parsed without
error but read the numbers wrongly. None of the sources raise when misread, so
the shape of the values has to be asserted instead.
"""
import math

import pytest

from radiological_material_clearance_finder import decay, nuclide
from radiological_material_clearance_finder.limits import get_limit_set, limit_sets

ALL_SETS = limit_sets()


def test_the_expected_sets_are_registered():
    for name in (
        "UK_EPR16_out_of_scope",
        "UK_IRR17_notification",
        "StrlSchV_unrestricted",
        "StrlSchV_metal_recycling",
        "Fetter",
        "NRC_long",
        "EU_BSS_clearance",
        "IAEA_GSR3_clearance",
    ):
        assert name in ALL_SETS


#: The sets parsed out of published markup, where a misread is silent. The US
#: sets are extracted from Python literals in OpenMC's source, so they carry no
#: markup to flatten, and 10 CFR 61.55 genuinely uses values like 3.5 and 4600.
MARKUP_SETS = tuple(
    name for name in ALL_SETS
    if get_limit_set(name).jurisdiction in {"UK", "Germany", "EU", "IAEA"}
)


@pytest.mark.parametrize("name", MARKUP_SETS)
def test_limits_are_one_significant_figure_times_a_power_of_ten(name):
    """The shape every one of these regulations uses, without exception.

    A value like 102 means superscript markup was flattened, turning 10^2 into
    the digits "102", which is the quiet failure this catches.
    """
    limit_set = get_limit_set(name)
    for nuclide_name, value in limit_set.limits.items():
        assert value > 0, f"{name}: {nuclide_name} has a non-positive limit"
        exponent = math.floor(math.log10(value))
        mantissa = value / 10.0**exponent
        assert abs(mantissa - round(mantissa)) < 1e-9 and 1 <= round(mantissa) <= 9, (
            f"{name}: {nuclide_name} = {value} is not one significant figure "
            f"times a power of ten"
        )


@pytest.mark.parametrize("name", ALL_SETS)
def test_every_key_is_a_valid_nuclide(name):
    limit_set = get_limit_set(name)
    for key in limit_set.limits:
        base = key[:-4] if key.endswith("_sec") else key
        assert nuclide.normalise(base) == base, f"{name}: {key} is not canonical"


@pytest.mark.parametrize("name", ALL_SETS)
def test_every_key_has_decay_data(name):
    """A limit for a nuclide we cannot compute an activity for is unusable."""
    data = decay.default_decay_data()
    for key in get_limit_set(name).limits:
        base = key[:-4] if key.endswith("_sec") else key
        assert data.knows(base), f"{name}: no decay data for {base}"


@pytest.mark.parametrize("name", ALL_SETS)
def test_secular_equilibrium_entries_are_valid_nuclides(name):
    limit_set = get_limit_set(name)
    for parent, daughters in limit_set.secular_equilibrium.items():
        assert nuclide.normalise(parent) == parent
        assert daughters, f"{name}: {parent} has an empty daughter list"
        for daughter in daughters:
            assert nuclide.normalise(daughter) == daughter


@pytest.mark.parametrize("name", ALL_SETS)
def test_provenance_is_recorded(name):
    limit_set = get_limit_set(name)
    assert limit_set.source, f"{name} has no source"
    assert limit_set.url, f"{name} has no url"
    assert limit_set.retrieved, f"{name} has no retrieval date"
    assert limit_set.units in {"Bq/g", "Ci/m3", "Bq"}


def test_table_sizes_are_what_the_sources_hold():
    """A short table means rows were skipped, which no other check would show."""
    expected = {
        "UK_EPR16_out_of_scope": 282,
        "UK_IRR17_notification": 300,
        "UK_IRR17_registration": 300,
        "StrlSchV_unrestricted": 763,
        "StrlSchV_metal_recycling": 283,
        "StrlSchV_soil": 113,
        "EU_BSS_clearance": 257,
        "IAEA_GSR3_clearance": 257,
        "Fetter": 81,
    }
    for name, count in expected.items():
        assert len(get_limit_set(name).limits) == count, (
            f"{name} has {len(get_limit_set(name).limits)} limits, expected {count}"
        )


def test_iaea_and_eu_agree_everywhere_they_overlap():
    """Two independent sources and extraction methods, one PDF and one XHTML."""
    iaea = get_limit_set("IAEA_GSR3_clearance").limits
    european = get_limit_set("EU_BSS_clearance").limits
    shared = set(iaea) & set(european)
    assert len(shared) == 257
    for name in shared:
        assert iaea[name] == pytest.approx(european[name], rel=1e-12), name


def test_uk_irr17_matches_the_eu_directive_it_transcribes():
    """IRR 2017 Schedule 7 carries the EU BSS Annex VII values.

    Two nuclides genuinely differ in the published UK text, and pinning them
    here means a future change to either table shows up as a test failure
    rather than as a silent drift.
    """
    uk = get_limit_set("UK_IRR17_notification").limits
    european = get_limit_set("EU_BSS_clearance").limits
    shared = set(uk) & set(european)
    differing = {n for n in shared if uk[n] != european[n]}
    assert differing == {"Na24", "Pt197"}


def test_uk_epr16_has_its_catch_all_and_scope_rule():
    limit_set = get_limit_set("UK_EPR16_out_of_scope")
    assert limit_set.default_limit == 0.01
    assert limit_set.min_half_life_scope == 100.0


def test_known_clearance_values():
    """Spot values read directly from the published tables."""
    assert get_limit_set("StrlSchV_unrestricted").limits["Co60"] == 0.1
    assert get_limit_set("StrlSchV_unrestricted").limits["H3"] == 100.0
    assert get_limit_set("StrlSchV_metal_recycling").limits["Co60"] == 0.6
    assert get_limit_set("StrlSchV_rubble").limits["Co60"] == 0.09
    assert get_limit_set("UK_EPR16_out_of_scope").limits["Co60"] == 0.1
    assert get_limit_set("UK_IRR17_notification").limits["H3"] == 100.0
    assert get_limit_set("EU_BSS_clearance").limits["Co60"] == 0.1
    assert get_limit_set("NRC_long").limits["C14"] == 8.0
    assert get_limit_set("NRC_long").metal_overrides["Ni59"] == 220.0
    assert get_limit_set("Fetter").limits["C14"] == 600.0


def test_progeny_lists_differ_between_regulations():
    """They are not interchangeable, which is why each set carries its own."""
    epr = get_limit_set("UK_EPR16_out_of_scope").secular_equilibrium
    irr = get_limit_set("UK_IRR17_notification").secular_equilibrium
    assert epr["Zr95"] == ("Nb95_m1",)
    assert irr["Zr95"] == ("Nb95",)
