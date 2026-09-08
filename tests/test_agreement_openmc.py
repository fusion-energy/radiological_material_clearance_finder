"""Arithmetic agreement with OpenMC, for the limit sets it also implements.

This validates the path from an inventory to a specific or volumetric activity,
which the table agreement tests cannot reach. Marked ``openmc`` and skipped when
OpenMC is absent.

    pytest -m openmc
"""
import warnings

import pytest

from radiological_material_clearance_finder import Material, clearance_index, nrc_waste_class

openmc = pytest.importorskip("openmc", reason="OpenMC is not installed")
pytestmark = pytest.mark.openmc

#: An activated steel, as atom densities in atoms per barn-cm.
INVENTORY = {
    "Fe54": 0.0049, "Fe56": 0.0771, "Fe57": 0.00178, "Fe58": 0.000237,
    "Cr52": 0.0157, "Ni58": 0.0079, "Ni59": 1.2e-6, "Ni63": 3.4e-7,
    "Co60": 2.1e-9, "C14": 5.5e-11, "Tc99": 8.0e-12, "I129": 1.1e-13,
    "Nb94": 4.4e-13, "H3": 2.2e-9, "Cs137": 3.1e-11, "Sr90": 1.7e-11,
    "Pu241": 2.0e-16, "Am241": 5.0e-17, "Mn54": 9.0e-10, "Cm242": 1e-18,
}
US_SETS = ["Fetter", "NRC_long", "NRC_short_A", "NRC_short_B", "NRC_short_C"]


@pytest.fixture
def pair():
    """The same inventory as an OpenMC material and as one of ours."""
    used_before = set(openmc.Material.used_ids)
    try:
        their_material = openmc.Material()
        for name, density in INVENTORY.items():
            their_material.add_nuclide(name, density)
        their_material.set_density("atom/b-cm", sum(INVENTORY.values()))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ours = Material.from_atom_densities(INVENTORY, name="steel")
        yield their_material, ours
    finally:
        # Creating a material consumes an auto-increment id, and depletion
        # results are indexed by material id, so the id state is restored.
        openmc.Material.used_ids.clear()
        openmc.Material.used_ids.update(used_before)


def test_mass_density_agrees(pair):
    theirs, ours = pair
    assert ours.density == pytest.approx(theirs.get_mass_density(), rel=1e-12)


@pytest.mark.parametrize("units", ["Bq/g", "Bq/cm3", "Ci/m3"])
def test_activity_agrees(pair, units):
    theirs, ours = pair
    assert ours.activity(units) == pytest.approx(theirs.get_activity(units), rel=1e-12)


@pytest.mark.parametrize("name", US_SETS)
@pytest.mark.parametrize("metal", [False, True])
def test_us_indexes_agree(pair, name, metal):
    theirs, ours = pair
    expected = theirs.waste_disposal_rating(limits=name, metal=metal)
    assert clearance_index(ours, name, metal=metal).index == pytest.approx(
        expected, rel=1e-12
    )


@pytest.mark.parametrize("name", US_SETS)
def test_per_nuclide_breakdown_agrees(pair, name):
    theirs, ours = pair
    expected = theirs.waste_disposal_rating(limits=name, by_nuclide=True)
    got = clearance_index(ours, name).by_nuclide
    assert set(got) == set(expected)
    for nuclide_name, value in expected.items():
        assert got[nuclide_name] == pytest.approx(value, rel=1e-12), nuclide_name


@pytest.mark.parametrize("metal", [False, True])
def test_waste_classification_agrees(pair, metal):
    theirs, ours = pair
    if not hasattr(theirs, "waste_classification"):
        pytest.skip("this OpenMC has no waste_classification")
    assert nrc_waste_class(ours, metal=metal) == theirs.waste_classification(metal=metal)


def test_half_lives_agree_with_openmc():
    """Both are ENDF/B-VIII.0, so any difference would show up in every index."""
    from radiological_material_clearance_finder import decay

    checked = 0
    for name in ("Co60", "Cs137", "H3", "Sr90", "Ni63", "C14", "Am241", "Pu241"):
        assert decay.half_life(name) == pytest.approx(
            openmc.data.half_life(name), rel=1e-12
        ), name
        checked += 1
    assert checked == 8


def test_atomic_masses_agree_with_openmc():
    from radiological_material_clearance_finder import decay

    for name in ("H1", "Fe56", "Co60", "U238", "Cs137"):
        assert decay.atomic_mass(name) == pytest.approx(
            openmc.data.atomic_mass(name), rel=1e-12
        ), name


def test_from_openmc_material_round_trips(pair):
    from radiological_material_clearance_finder.openmc_interop import from_openmc_material

    theirs, ours = pair
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        converted = from_openmc_material(theirs)
    assert converted.density == pytest.approx(ours.density, rel=1e-12)
    assert converted.activity("Bq/g") == pytest.approx(ours.activity("Bq/g"), rel=1e-12)
