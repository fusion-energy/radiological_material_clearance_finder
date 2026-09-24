"""Stamp a release version into the workspace, from a GitHub release tag.

The crate and the wheel share one version, [workspace.package] version in the
root Cargo.toml, which cargo and maturin both read. The repository carries
0.0.0 there, and publish.yml runs this before every build so the published
artifacts carry the release tag's version instead. There is nothing to bump by
hand: tag the release v0.2.0 and that is the version.

    python tools/set_version.py v0.2.0          # rewrite Cargo.toml
    python tools/set_version.py --check v0.2.0  # only validate and print

Prints ``cargo=<version>`` and ``python=<version>`` lines, suitable for
``$GITHUB_OUTPUT``. The two differ only for pre-releases: cargo spells one
0.2.0-rc.1 and the wheel, following PEP 440 as maturin converts it, 0.2.0rc1.
"""
import argparse
import re
import sys
from pathlib import Path

CARGO_TOML = Path(__file__).resolve().parents[1] / "Cargo.toml"

# Plain releases and the three pre-release kinds that have a PEP 440 spelling.
# Anything else would give a crate version PyPI cannot represent, or the other
# way round, so it is refused rather than guessed at.
_TAG = re.compile(r"^v?(?P<base>\d+\.\d+\.\d+)(?:-(?P<kind>alpha|beta|rc)\.(?P<number>\d+))?$")
_PEP440 = {"alpha": "a", "beta": "b", "rc": "rc"}

# The version line inside [workspace.package], and no other section's.
_WORKSPACE_VERSION = re.compile(
    r'(^\[workspace\.package\]\n(?:(?!^\[).*\n)*?^version = ")[^"]*(")', re.MULTILINE
)


def versions(tag: str) -> tuple[str, str]:
    """The cargo and Python versions a release tag stands for."""
    match = _TAG.match(tag.strip())
    if match is None:
        raise ValueError(
            f"release tag {tag!r} is not a version. Use v1.2.3, or v1.2.3-rc.1, "
            f"-alpha.1 or -beta.1 for a pre-release."
        )
    base, kind, number = match.group("base", "kind", "number")
    if kind is None:
        return base, base
    return f"{base}-{kind}.{number}", f"{base}{_PEP440[kind]}{number}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("tag", help="release tag, such as v0.2.0")
    parser.add_argument("--check", action="store_true", help="validate only, change nothing")
    args = parser.parse_args()

    try:
        cargo, python = versions(args.tag)
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    if not args.check:
        text = CARGO_TOML.read_text()
        stamped, count = _WORKSPACE_VERSION.subn(rf"\g<1>{cargo}\g<2>", text)
        if count != 1:
            print(f"error: no [workspace.package] version found in {CARGO_TOML}", file=sys.stderr)
            return 1
        CARGO_TOML.write_text(stamped)

    print(f"cargo={cargo}")
    print(f"python={python}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
