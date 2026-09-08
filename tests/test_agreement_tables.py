"""Table agreement with OpenMC, without importing OpenMC.

The tables are the biggest risk in this package: a single wrong number is
invisible in every other test. OpenMC's ``waste.py`` and the vendored copy of
openmc-dev/openmc#3898 are independent transcriptions of the same regulations,
so their dict literals are recovered here with :mod:`ast` and compared. Reading
the source as text rather than importing it means this runs whether or not
OpenMC is installed, and executes nothing from those files.
"""
import sys
from pathlib import Path

import pytest

from radiological_material_clearance_finder import nuclide
from radiological_material_clearance_finder.limits import get_limit_set

# Reuse the extractor the build script already uses, rather than writing a
# second one that could be wrong in a different way.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from build_us import _branches, _extract  # noqa: E402

OPENMC_WASTE = Path.home() / "openmc" / "openmc" / "waste.py"
PR_WASTE = (
    Path.home()
    / "Downloads"
    / "neutronics_analysis-main"
    / "openmc_strlschv"
    / "_waste_pr3898.py"
)


def literal_limits(path: Path, set_name: str) -> tuple[dict, dict, dict]:
    """Recover the limit dicts for one limit set from a source file.

    Returns:
        ``(base, metal_overrides, per_gram)``, each keyed by canonical name.
    """
    branches = _branches(path.read_text())
    if set_name not in branches:
        pytest.skip(f"{set_name} is not in {path.name}")
    base, metal, per_gram = _extract(branches[set_name])
    canonical = lambda mapping: {nuclide.normalise(k): v for k, v in mapping.items()}
    return canonical(base), canonical(metal), canonical(per_gram)


@pytest.mark.skipif(not OPENMC_WASTE.exists(), reason="no OpenMC source checkout")
@pytest.mark.parametrize(
    "name", ["Fetter", "NRC_long", "NRC_short_A", "NRC_short_B", "NRC_short_C"]
)
def test_us_tables_match_openmc(name):
    """These are extracted from OpenMC's source, so they must match exactly."""
    base, metal, per_gram = literal_limits(OPENMC_WASTE, name)
    limit_set = get_limit_set(name)
    assert limit_set.limits == pytest.approx(base)
    assert limit_set.metal_overrides == pytest.approx(metal)
    assert limit_set.limits_per_gram == pytest.approx(per_gram)


@pytest.mark.skipif(not PR_WASTE.exists(), reason="no vendored PR #3898 copy")
@pytest.mark.parametrize(
    "name",
    [
        "StrlSchV_unrestricted",
        "StrlSchV_rubble",
        "StrlSchV_soil",
        "StrlSchV_landfill_100",
        "StrlSchV_landfill_1000",
        "StrlSchV_incineration_100",
        "StrlSchV_incineration_1000",
        "StrlSchV_metal_recycling",
    ],
)
def test_german_tables_agree_with_the_pr_transcription(name):
    """Every value shared with PR #3898 must agree.

    Our tables are a superset: they are parsed from the full Anlage 4 Tabelle 1
    while the PR transcribed part of it. Two documented cases are excluded.

    Th-232 is listed twice in the regulation, plain and marked "+", and this
    package keeps both and picks between them from the material, while the PR
    kept only the plain value. Six nuclides, listed below, have no landfill
    value in the regulation at all: their Spalte 8 and 10 cells are empty and
    only the incineration columns are filled, which was checked against the
    published table cell by cell.
    """
    misplaced_by_the_pr = {"Au195", "Ho166_m1", "Lu177_m1", "Pm146", "Tb158", "Yb169"}
    theirs, _, _ = literal_limits(PR_WASTE, name)
    ours = get_limit_set(name).limits

    for key, value in theirs.items():
        if key == "Th232" or (key in misplaced_by_the_pr and "landfill" in name):
            continue
        assert key in ours, f"{name}: {key} is in the PR but missing from our table"
        assert ours[key] == pytest.approx(value, rel=1e-12), f"{name}: {key}"


@pytest.mark.skipif(not PR_WASTE.exists(), reason="no vendored PR #3898 copy")
def test_our_german_tables_are_a_superset_of_the_pr():
    """Parsing the official source picks up rows the transcription missed."""
    theirs, _, _ = literal_limits(PR_WASTE, "StrlSchV_unrestricted")
    ours = get_limit_set("StrlSchV_unrestricted").limits
    assert len(ours) > len(theirs)
    assert set(theirs) - set(ours) == set()
