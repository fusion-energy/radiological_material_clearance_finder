"""Nuclide name parsing, including the spellings the regulations use."""
import pytest

from radiological_material_clearance_finder import nuclide as nuc


@pytest.mark.parametrize(
    "given, expected",
    [
        ("Co60", "Co60"),
        ("Co-60", "Co60"),
        ("co 60", "Co60"),
        ("CO-60", "Co60"),
        ("Ag-108m", "Ag108_m1"),
        ("Ag108m", "Ag108_m1"),
        ("Ag108_m1", "Ag108_m1"),
        ("Ag108m2", "Ag108_m2"),
        # "n" is the second metastable state, as StrlSchV writes it.
        ("Hf-178n", "Hf178_m2"),
        ("Ir-194n", "Ir194_m2"),
        # EUR-Lex separates the metastable marker with a space.
        ("Zn-69 m", "Zn69_m1"),
        ("Am-242 m", "Am242_m1"),
        # A trailing "+" marks a secular equilibrium value, not another nuclide.
        ("Sr-90+", "Sr90"),
        ("H-3", "H3"),
        ("Og295", "Og295"),
    ],
)
def test_normalise(given, expected):
    assert nuc.normalise(given) == expected


@pytest.mark.parametrize(
    "symbol", ["Rn222", "Zn69", "Sn113", "In111", "Mn54", "N14", "Ne20"]
)
def test_element_symbols_ending_in_n_are_not_metastable(symbol):
    """The "n" metastable marker must not eat an element symbol."""
    assert nuc.normalise(symbol) == symbol
    assert nuc.metastable_state(symbol) == 0


@pytest.mark.parametrize(
    "label, expected",
    [
        ("Sr-90+", ("Sr90", "+")),
        ("U-238sec", ("U238", "sec")),
        ("Th-232 sec", ("Th232", "sec")),
        ("Co-60", ("Co60", None)),
        ("Es-254m+", ("Es254_m1", "+")),
    ],
)
def test_parse_regulatory_separates_the_marker(label, expected):
    assert nuc.parse_regulatory(label) == expected


def test_bare_element_is_rejected_with_a_useful_message():
    with pytest.raises(nuc.NuclideNameError, match="element, not a nuclide"):
        nuc.normalise("Fe")


@pytest.mark.parametrize("bad", ["", "Xx60", "H-500", None, "not a nuclide", 42])
def test_invalid_names_raise(bad):
    with pytest.raises(nuc.NuclideNameError):
        nuc.parse(bad)


def test_mass_number_below_proton_number_is_rejected():
    with pytest.raises(nuc.NuclideNameError, match="below its proton number"):
        nuc.parse("U50")


def test_accessors():
    assert nuc.element("Ag108_m1") == "Ag"
    assert nuc.mass_number("Ag108_m1") == 108
    assert nuc.metastable_state("Ag108_m1") == 1
    assert nuc.atomic_number("Co60") == 27
    assert nuc.is_valid("Co60") and not nuc.is_valid("Fe")
