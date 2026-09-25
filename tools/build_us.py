"""Build the US limit sets from OpenMC's waste.py.

    python tools/build_us.py [--waste-py PATH]

The Fetter and 10 CFR 61.55 tables are extracted from the source of
``openmc/waste.py`` with :mod:`ast`, not by importing OpenMC, so this runs
without OpenMC installed and cannot execute anything from that file. Extracting
rather than retyping is what makes the numbers agree with OpenMC by
construction.

The NRC values were verified against the eCFR text of 10 CFR 61.55 as of
2026-09-01, which is unchanged since 66 FR 55792 (2 November 2001).
"""
from __future__ import annotations

import argparse
import ast
import json
import sys
from datetime import date
from pathlib import Path

from _tables import portable_origin

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

DATA = Path(__file__).resolve().parents[1] / "crates" / "radiological-material-clearance-finder" / "data" / "limits"
DEFAULT_WASTE_PY = Path.home() / "openmc" / "openmc" / "waste.py"


def _branches(source: str) -> dict[str, list[ast.stmt]]:
    """Return the body of each ``limits == '<name>'`` branch, keyed by name."""
    tree = ast.parse(source)
    function = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "_waste_disposal_rating"
    )
    found: dict[str, list[ast.stmt]] = {}
    stack = list(function.body)
    while stack:
        node = stack.pop(0)
        if not isinstance(node, ast.If):
            continue
        test = node.test
        if (
            isinstance(test, ast.Compare)
            and isinstance(test.left, ast.Name)
            and test.left.id == "limits"
            and isinstance(test.comparators[0], ast.Constant)
        ):
            found[test.comparators[0].value] = node.body
        stack.extend(node.orelse)
    return found


def _extract(body: list[ast.stmt]) -> tuple[dict, dict, dict]:
    """Pull base limits, metal overrides and per-gram limits out of a branch."""
    base: dict[str, float] = {}
    metal: dict[str, float] = {}
    per_gram: dict[str, float] = {}

    for node in body:
        # limits = {...}, possibly with "X if metal else Y" values
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Dict):
            for key, value in zip(node.value.keys, node.value.values):
                name = key.value
                if isinstance(value, ast.IfExp):
                    # 35.0 if metal else 3.5
                    metal[name] = ast.literal_eval(value.body)
                    base[name] = ast.literal_eval(value.orelse)
                else:
                    base[name] = ast.literal_eval(value)

        # if metal: limits['C14'] = 80.0
        elif isinstance(node, ast.If) and isinstance(node.test, ast.Name) and node.test.id == "metal":
            for inner in node.body:
                if isinstance(inner, ast.Assign) and isinstance(inner.targets[0], ast.Subscript):
                    metal[inner.targets[0].slice.value] = ast.literal_eval(inner.value)

        # limits.update({'Pu241': 3500.0 * factor})  -> nCi/g entries
        elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            call = node.value
            if getattr(call.func, "attr", None) == "update" and isinstance(call.args[0], ast.Dict):
                for key, value in zip(call.args[0].keys, call.args[0].values):
                    if isinstance(value, ast.BinOp) and isinstance(value.op, ast.Mult):
                        per_gram[key.value] = ast.literal_eval(value.left)
                    else:
                        base[key.value] = ast.literal_eval(value)
    return base, metal, per_gram


FETTER_NOTES = (
    "Table 2 of Fetter, Cheng and Mann, Long-term radioactive waste from fusion "
    "reactors: Part II, Fusion Engineering and Design 13(2) 239-246 (1990), "
    "doi:10.1016/0920-3796(90)90104-E. Specific activity limits for class C "
    "disposal of activated metal, for nuclides with half-lives between 5 y and "
    "1e12 y. Where the paper gives a range, the conservative lower bound is used "
    "here and the upper bound is recorded in limits_upper."
)

NRC_NOTES = (
    "10 CFR 61.55 Tables 1 and 2, verified against the eCFR text dated "
    "2026-09-01. Values in Ci/m3 except the transuranic entries of Table 1, "
    "which the regulation gives in nCi/g and which are held in limits_per_gram. "
    "The activated metal values are separate rows in the regulation, so metal "
    "overrides both replace limits (C-14, Ni-63) and add nuclides that have no "
    "non-metal limit at all (Ni-59, Nb-94)."
)

LABELS = {
    "Fetter": ("Fetter class C disposal limits", FETTER_NOTES),
    "NRC_long": ("10 CFR 61.55 Table 1, long lived radionuclides", NRC_NOTES),
    "NRC_short_A": ("10 CFR 61.55 Table 2 column 1, Class A", NRC_NOTES),
    "NRC_short_B": ("10 CFR 61.55 Table 2 column 2, Class B", NRC_NOTES),
    "NRC_short_C": ("10 CFR 61.55 Table 2 column 3, Class C", NRC_NOTES),
}

URLS = {
    "Fetter": "https://doi.org/10.1016/0920-3796(90)90104-E",
}
NRC_URL = "https://www.ecfr.gov/current/title-10/chapter-I/part-61/subpart-C/section-61.55"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--waste-py", type=Path, default=DEFAULT_WASTE_PY)
    args = parser.parse_args()

    if not args.waste_py.exists():
        raise SystemExit(f"not found: {args.waste_py}")
    branches = _branches(args.waste_py.read_text())

    sets = []
    for name in ("Fetter", "NRC_long", "NRC_short_A", "NRC_short_B", "NRC_short_C"):
        if name not in branches:
            raise SystemExit(f"{name} branch not found in {args.waste_py}")
        base, metal, per_gram = _extract(branches[name])
        label, notes = LABELS[name]
        entry = {
            "name": name,
            "label": label,
            "jurisdiction": "US",
            "units": "Ci/m3",
            "limits": base,
            "threshold": 1.0,
            "source": "Fetter et al. (1990)" if name == "Fetter" else "10 CFR 61.55",
            "url": URLS.get(name, NRC_URL),
            "retrieved": date.today().isoformat(),
            "notes": notes,
        }
        if metal:
            entry["metal_overrides"] = metal
        if per_gram:
            entry["limits_per_gram"] = per_gram
        if name == "NRC_short_A":
            entry["dynamic_rule"] = "nrc_short_lived_class_a"
        sets.append(entry)
        print(f"{name:14} {len(base):3} limits, {len(metal)} metal, {len(per_gram)} nCi/g")

    payload = {
        "_generated_by": "tools/build_us.py",
        "_extracted_from": portable_origin(args.waste_py),
        "sets": sets,
    }
    (DATA / "us.json").write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {DATA / 'us.json'}")


if __name__ == "__main__":
    main()
