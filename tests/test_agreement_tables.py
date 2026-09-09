"""Table agreement against frozen copies of two independent transcriptions.

The tables are the biggest risk in this package: one wrong number is invisible
in every other test. OpenMC's ``waste.py`` and openmc-dev/openmc#3898 are
independent transcriptions of the same regulations, so their dict literals were
extracted with :mod:`ast` and frozen into ``tests/data``.

Comparing against the frozen copies rather than against a checkout somewhere on
the machine is what makes these tests run at all. Keyed to a developer's home
directory they skipped everywhere else, including CI, which is the same as not
having them. ``tools/verify_tables.py`` still re-derives everything from the
live sources, and that is where drift in the upstream text is meant to surface.
"""
import json
from pathlib import Path

import pytest

from radiological_material_clearance_finder.limits import get_limit_set

DATA = Path(__file__).parent / "data"
OPENMC_TABLES = json.loads((DATA / "openmc_waste_tables.json").read_text())["sets"]
PR3898_TABLES = json.loads((DATA / "pr3898_strlschv_tables.json").read_text())["sets"]


@pytest.mark.parametrize("name", sorted(OPENMC_TABLES))
def test_us_tables_match_openmc(name):
    """Extracted from OpenMC's own source, so every field must match exactly."""
    expected = OPENMC_TABLES[name]
    limit_set = get_limit_set(name)
    assert limit_set.limits == pytest.approx(expected["limits"])
    assert limit_set.metal_overrides == pytest.approx(expected["metal_overrides"])
    assert limit_set.limits_per_gram == pytest.approx(expected["limits_per_gram"])


@pytest.mark.parametrize("name", sorted(PR3898_TABLES))
def test_german_tables_agree_with_the_pr_transcription(name):
    """Every value shared with PR #3898 must agree.

    Our tables are a superset, parsed from the whole of Anlage 4 Tabelle 1 while
    the PR transcribed part of it. Two documented differences are excluded.

    Th-232 appears twice in the regulation, plain and marked "+", and this
    package keeps both and chooses between them from the material, while the PR
    kept only the plain value.

    Six nuclides have no landfill value in the regulation at all. Their Spalte 8
    and Spalte 10 cells are empty and only the incineration columns are filled,
    checked cell by cell against the published table, so the PR placing them in
    the landfill sets is an error in the transcription rather than in the parse.
    """
    misplaced_by_the_pr = {"Au195", "Ho166_m1", "Lu177_m1", "Pm146", "Tb158", "Yb169"}
    theirs = PR3898_TABLES[name]
    ours = get_limit_set(name).limits

    for key, value in theirs.items():
        if key == "Th232" or (key in misplaced_by_the_pr and "landfill" in name):
            continue
        assert key in ours, f"{name}: {key} is in the PR but missing from our table"
        assert ours[key] == pytest.approx(value, rel=1e-12), f"{name}: {key}"


def test_our_german_tables_are_a_superset_of_the_pr():
    """Parsing the official source picks up rows the transcription missed."""
    theirs = PR3898_TABLES["StrlSchV_unrestricted"]
    ours = get_limit_set("StrlSchV_unrestricted").limits
    assert len(ours) > len(theirs)
    assert set(theirs) - set(ours) == set()


def test_the_frozen_fixtures_are_not_empty():
    """A fixture that failed to load would make every test above vacuous."""
    assert len(OPENMC_TABLES) == 5
    assert len(PR3898_TABLES) == 8
    assert sum(len(v) for v in PR3898_TABLES.values()) > 1500
    assert len(OPENMC_TABLES["Fetter"]["limits"]) == 81
