"""Table agreement against frozen copies of two independent transcriptions.

The tables are the biggest risk in this package: one wrong number is invisible
in every other test. OpenMC's ``waste.py`` and openmc-dev/openmc#3898 are
independent transcriptions of the same regulations, so their dict literals were
extracted with :mod:`ast` and frozen into ``tests/data``.

Comparing against the frozen copies rather than against a checkout somewhere on
the machine is what makes these tests run at all. Keyed to a developer's home
directory they skipped everywhere else, including CI, which is the same as not
having them.

Be clear about what this does and does not establish. The fixtures were produced
by the same extractor the build script uses, so these tests cannot catch a bug in
the extraction itself. What they do catch is the shipped table drifting away from
what was extracted, which is not hypothetical: a stray write once left NRC_long
ten times too lenient and this is what found it. Independent confirmation of the
values themselves comes from
`test_every_limit_set_has_at_least_one_pinned_value` in test_data_integrity.py,
whose expected numbers are read by hand from 10 CFR 61.55 and the published
tables rather than from any code. ``tools/verify_tables.py`` still re-derives everything from the
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

    Th-232 appears twice in the regulation, plain and marked "+". This package
    keeps both, the plain value in `limits` and the marked one in
    `limits_secular_equilibrium`, while the PR kept only one. Rather than skip
    the comparison, the PR's value is matched against whichever of the two it
    corresponds to, so all eight sets are still compared.

    Six nuclides have no landfill value in the regulation at all. Their Spalte 8
    and Spalte 10 cells are empty and only the incineration columns are filled,
    checked cell by cell against the published table, so the PR placing them in
    the landfill sets is an error in the transcription rather than in the parse.
    """
    misplaced_by_the_pr = {"Au195", "Ho166_m1", "Lu177_m1", "Pm146", "Tb158", "Yb169"}
    theirs = PR3898_TABLES[name]
    ours = get_limit_set(name).limits

    limit_set = get_limit_set(name)
    for key, value in theirs.items():
        if key in misplaced_by_the_pr and "landfill" in name:
            continue
        assert key in ours, f"{name}: {key} is in the PR but missing from our table"
        if key == "Th232":
            # One of our two values must be theirs, rather than neither.
            variant = limit_set.limits_secular_equilibrium.get("Th232")
            assert value == pytest.approx(ours[key], rel=1e-12) or (
                variant is not None and value == pytest.approx(variant, rel=1e-12)
            ), f"{name}: Th232 is {value} in the PR, ours are {ours[key]} and {variant}"
            continue
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
