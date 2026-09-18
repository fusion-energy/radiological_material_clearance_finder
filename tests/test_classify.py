"""Waste classification against the US and UK category schemes."""
import pytest

from radiological_material_clearance_finder import (
    Material,
    alpha_activity,
    beta_gamma_activity,
    nrc_waste_class,
    uk_waste_category,
)
from radiological_material_clearance_finder.classify import (
    UK_HIGH_VOLUME_VLLW_BQ_PER_G,
    UK_LLW_ALPHA_BQ_PER_G,
    UK_LLW_BETA_GAMMA_BQ_PER_G,
)


def steel(activities, density=7.87):
    return Material.from_specific_activities(activities, density=density)


def test_a_clean_material_is_class_a():
    assert nrc_waste_class(steel({"Fe55": 1.0})) == "Class A"


#: Bq/g per Ci/m3 for a material of this density, so the NRC's volumetric
#: limits can be aimed at directly.
def bq_per_g(curies_per_cubic_metre, density=7.87):
    return curies_per_cubic_metre * 3.7e10 / (density * 1e6)


def test_class_rises_with_activity():
    """Cs-137 has Class A, B and C limits of 1, 44 and 4600 Ci/m3."""
    classes = [
        nrc_waste_class(steel({"Cs137": bq_per_g(target)}))
        for target in (0.5, 10.0, 1000.0, 1e5)
    ]
    assert classes == ["Class A", "Class B", "Class C", "GTCC"]


def test_the_class_boundaries_sit_at_the_tabulated_limits():
    """Cs-137 crosses out of Class A at exactly 1 Ci/m3."""
    assert nrc_waste_class(steel({"Cs137": bq_per_g(0.99)})) == "Class A"
    assert nrc_waste_class(steel({"Cs137": bq_per_g(1.01)})) == "Class B"
    assert nrc_waste_class(steel({"Cs137": bq_per_g(43.0)})) == "Class B"
    assert nrc_waste_class(steel({"Cs137": bq_per_g(45.0)})) == "Class C"


def test_alpha_and_beta_gamma_split_by_decay_mode():
    material = steel({"Am241": 10.0, "Cs137": 90.0})
    assert alpha_activity(material) == pytest.approx(10.0)
    assert beta_gamma_activity(material) == pytest.approx(90.0)


def test_a_branching_nuclide_contributes_to_both():
    """Bi-212 branches 35.94 percent alpha."""
    material = steel({"Bi212": 100.0})
    assert alpha_activity(material) == pytest.approx(35.94, abs=1e-2)
    assert beta_gamma_activity(material) == pytest.approx(64.06, abs=1e-2)
    assert alpha_activity(material) + beta_gamma_activity(material) == pytest.approx(100.0)


def test_uk_vllw_boundary():
    below = steel({"Cs137": UK_HIGH_VOLUME_VLLW_BQ_PER_G * 0.9})
    above = steel({"Cs137": UK_HIGH_VOLUME_VLLW_BQ_PER_G * 1.1})
    assert uk_waste_category(below).category == "VLLW"
    assert uk_waste_category(above).category == "LLW"


def test_uk_tritium_has_its_own_vllw_allowance():
    """Tritium is allowed to 40 MBq/te where everything else stops at 4."""
    assert uk_waste_category(steel({"H3": 30.0})).category == "VLLW"
    assert uk_waste_category(steel({"H3": 50.0})).category == "LLW"
    # The same activity of anything else is LLW.
    assert uk_waste_category(steel({"Cs137": 30.0})).category == "LLW"


def test_uk_llw_upper_bounds():
    assert uk_waste_category(
        steel({"Cs137": UK_LLW_BETA_GAMMA_BQ_PER_G * 0.9})
    ).category == "LLW"
    assert uk_waste_category(
        steel({"Cs137": UK_LLW_BETA_GAMMA_BQ_PER_G * 1.1})
    ).category == "ILW"
    assert uk_waste_category(
        steel({"Am241": UK_LLW_ALPHA_BQ_PER_G * 1.1})
    ).category == "ILW"


def test_the_alpha_limit_is_stricter_than_the_beta_gamma_one():
    """4 GBq/te alpha against 12 GBq/te beta and gamma."""
    activity = 6.0e3
    assert uk_waste_category(steel({"Am241": activity})).category == "ILW"
    assert uk_waste_category(steel({"Cs137": activity})).category == "LLW"


def test_the_category_explains_itself():
    result = uk_waste_category(steel({"Am241": 1e5}))
    assert result.category == "ILW"
    assert "alpha" in result.reason
    assert "4 GBq/te" in result.reason


# The 61.55 combination rules, which decide the class when Table 1 and Table 2
# nuclides are both present. C-14 is a Table 1 nuclide with an 8 Ci/m3 limit and
# Cs-137 is a Table 2 one, so the two ratios can be aimed independently.
def test_table1_only_uses_rule_a3():
    """61.55(a)(3): Class A below 0.1, Class C below 1.0, GTCC at or above."""
    assert nrc_waste_class(steel({"C14": bq_per_g(0.4)})) == "Class A"
    assert nrc_waste_class(steel({"C14": bq_per_g(4.0)})) == "Class C"
    assert nrc_waste_class(steel({"C14": bq_per_g(9.0)})) == "GTCC"


def test_both_tables_use_rule_a5():
    """61.55(a)(5): Table 1 below 0.1 defers to Table 2, otherwise it dominates."""
    # Table 1 ratio 0.05, so the Table 2 column decides.
    assert nrc_waste_class(
        steel({"C14": bq_per_g(0.4), "Cs137": bq_per_g(10.0)})
    ) == "Class B"
    # Table 1 ratio 0.5, which caps the answer at Class C.
    assert nrc_waste_class(
        steel({"C14": bq_per_g(4.0), "Cs137": bq_per_g(10.0)})
    ) == "Class C"
    # Table 1 ratio 0.5 with Table 2 above its Class C column gives GTCC.
    assert nrc_waste_class(
        steel({"C14": bq_per_g(4.0), "Cs137": bq_per_g(1e5)})
    ) == "GTCC"
    # Table 1 at or above 1.0 is GTCC whatever Table 2 says.
    assert nrc_waste_class(
        steel({"C14": bq_per_g(9.0), "Cs137": bq_per_g(0.1)})
    ) == "GTCC"


def test_neither_table_uses_rule_a6():
    """61.55(a)(6): nothing from either table present means Class A."""
    assert nrc_waste_class(steel({"Fe55": 1e6})) == "Class A"


def test_carbon_14_shares_the_tritium_vllw_allowance():
    """The 2007 policy pairs them, so C-14 gets the 40 MBq/te allowance too.

    Counting it in the general 4 MBq/te term instead classifies a material a
    whole category too high.
    """
    # 30 Bq/g of C-14 alone: within the shared 40 Bq/g allowance.
    assert uk_waste_category(steel({"C14": 30.0})).category == "VLLW"
    # Anything without that allowance is LLW at the same activity.
    assert uk_waste_category(steel({"Cs137": 30.0})).category == "LLW"
    # Tritium and C-14 share one allowance, so 25 each exceeds it together.
    assert uk_waste_category(steel({"H3": 25.0, "C14": 25.0})).category == "LLW"
    assert uk_waste_category(steel({"H3": 15.0, "C14": 15.0})).category == "VLLW"


def test_the_shared_allowance_is_reported():
    result = uk_waste_category(steel({"H3": 10.0, "C14": 20.0}))
    assert result.tritium_and_c14 == pytest.approx(30.0)
    assert "carbon-14" in result.reason
