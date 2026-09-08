"""Re-download every official source and diff it against the shipped tables.

    python tools/verify_tables.py

Regulations are amended and published texts get corrected. This re-runs each
build script against a freshly downloaded source and reports any difference, so
drift is caught deliberately rather than discovered in a result.

The shipped tables are restored afterwards whatever happens, so this never
leaves the working tree modified. Exits non-zero if anything differs.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src" / "radiological_material_clearance_finder" / "data"
LIMITS = DATA / "limits"

BUILDERS = [
    ("uk_epr16.json", "build_uk_epr16.py"),
    ("uk_irr17.json", "build_uk_irr17.py"),
    ("de_strlschv.json", "build_de_strlschv.py"),
    ("eu_bss.json", "build_eu_bss.py"),
    ("iaea.json", "build_iaea.py"),
]


def limits_of(path: Path) -> dict:
    """The limit values only, ignoring the retrieval date and prose."""
    payload = json.loads(path.read_text())
    return {
        entry["name"]: {
            "limits": entry.get("limits", {}),
            "default_limit": entry.get("default_limit"),
            "secular_equilibrium": entry.get("secular_equilibrium", {}),
            "limits_secular_equilibrium": entry.get("limits_secular_equilibrium", {}),
        }
        for entry in payload["sets"]
    }


def compare(name: str, before: dict, after: dict) -> list[str]:
    problems = []
    for set_name in sorted(set(before) | set(after)):
        if set_name not in after:
            problems.append(f"{name}: {set_name} disappeared from the source")
            continue
        if set_name not in before:
            problems.append(f"{name}: {set_name} is new in the source")
            continue
        old, new = before[set_name], after[set_name]
        for field in ("limits", "default_limit", "secular_equilibrium",
                      "limits_secular_equilibrium"):
            if old[field] == new[field]:
                continue
            if field in ("default_limit",):
                problems.append(
                    f"{set_name}: {field} changed from {old[field]} to {new[field]}"
                )
                continue
            added = sorted(set(new[field]) - set(old[field]))
            removed = sorted(set(old[field]) - set(new[field]))
            changed = sorted(
                key for key in set(old[field]) & set(new[field])
                if old[field][key] != new[field][key]
            )
            if added:
                problems.append(f"{set_name}: {field} gained {len(added)}: {added[:8]}")
            if removed:
                problems.append(f"{set_name}: {field} lost {len(removed)}: {removed[:8]}")
            if changed:
                detail = [(k, old[field][k], new[field][k]) for k in changed[:6]]
                problems.append(f"{set_name}: {field} changed {len(changed)}: {detail}")
    return problems


def main() -> int:
    with tempfile.TemporaryDirectory() as backup_dir:
        backup = Path(backup_dir) / "limits"
        shutil.copytree(LIMITS, backup)
        problems: list[str] = []
        try:
            for filename, script in BUILDERS:
                print(f"--- {script}")
                before = limits_of(backup / filename)
                result = subprocess.run(
                    [sys.executable, str(ROOT / "tools" / script)],
                    capture_output=True, text=True,
                )
                if result.returncode != 0:
                    problems.append(
                        f"{script} failed: {result.stderr.strip().splitlines()[-1:]}"
                    )
                    continue
                problems.extend(compare(filename, before, limits_of(LIMITS / filename)))
        finally:
            shutil.rmtree(LIMITS)
            shutil.copytree(backup, LIMITS)

    if problems:
        print(f"\n{len(problems)} difference(s) between the sources and the shipped tables:")
        for problem in problems:
            print(f"  {problem}")
        print("\nThe shipped tables have been restored. Re-run the build scripts to adopt "
              "the change once it has been reviewed.")
        return 1
    print("\nEvery table matches its official source.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
