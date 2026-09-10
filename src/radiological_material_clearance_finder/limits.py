"""Limit sets: the regulatory tables a material is measured against.

Every limit set is data, loaded from JSON in ``data/limits``. Adding a
jurisdiction means adding a file and a build script, never an ``elif`` in the
calculation. Only two things in these regulations cannot be expressed as a
table, and both are named explicitly rather than worked around:

* limits given per gram of material (the NRC's nCi/g transuranic entries), which
  need the material's density before they can be compared against Ci/m3, and
* limits that depend on the material's own nuclides (the NRC Class A rule that
  anything with a half-life under five years takes a 700 Ci/m3 limit), held in
  `DYNAMIC_RULES`.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping

from . import nuclide as _nuclide

__all__ = ["LimitSet", "limit_sets", "get_limit_set", "register_limit_set", "DYNAMIC_RULES"]

#: Units a limit set may be expressed in.
ACTIVITY_UNITS = frozenset({"Bq/g", "Ci/m3", "Bq"})

_LIMITS_DIR = Path(__file__).parent / "data" / "limits"

#: Five years in seconds, the NRC's short-lived and long-lived boundary.
_FIVE_YEARS = 5.0 * 365.25 * 86400.0


@dataclass(frozen=True)
class LimitSet:
    """A named table of per-nuclide activity limits.

    Attributes:
        name: Stable identifier, such as ``UK_EPR16_out_of_scope``.
        label: Human readable description for reports.
        jurisdiction: Issuing authority, such as ``UK`` or ``Germany``.
        units: Activity units the limits are in: ``Bq/g``, ``Ci/m3``, or ``Bq``
            for a total activity limit, which needs the material's mass.
        limits: Limit per canonical nuclide name.
        default_limit: Limit applied to a nuclide absent from ``limits``, or
            ``None`` where the regulation has no catch-all. Three UK sets define
            one: 0.01 Bq/g for ``UK_EPR16_out_of_scope`` and
            ``UK_IRR17_notification``, and 0.1 Bq/g for
            ``UK_IRR17_registration``. The other UK sets, and every German, US,
            EU and IAEA set, have none.
        unlimited: Nuclides the source explicitly places no limit on. They are
            covered, and contribute nothing, which is different from being
            absent.
        secular_equilibrium: Parent to the daughters whose contribution the
            parent's limit already includes, straight from the regulation's own
            table. These are the "+" rows.
        secular_equilibrium_sec: Parent to the daughters covered by its whole
            chain "sec" value, which is a different and usually much longer list
            carrying a different and much stricter limit. UK EPR 2016 gives
            U-238 three progeny at 1 Bq/g under "U-238+" and fourteen at 0.01
            Bq/g under "U-238sec". The stricter value is stored in ``limits``
            under a ``_sec`` key.
        limits_secular_equilibrium: The limit to use for a parent that the
            source lists twice, once plain and once marked "+", when its
            daughters are actually present. StrlSchV gives Th-232 as 10 Bq/g
            plain and 0.01 Bq/g marked, and which one applies depends on the
            material, not on the table.
        metal_overrides: Limits replacing or adding to ``limits`` when the
            material is activated metal.
        limits_per_gram: Limits in nCi/g, converted using the material density.
        limits_upper: Upper end of a limit given as a range in the source, kept
            for reference. ``limits`` holds the conservative lower end.
        dynamic_rule: Key into `DYNAMIC_RULES` for a rule that depends on
            the material's own nuclides.
        min_half_life_scope: Half-life in seconds below which, if *every*
            radionuclide present falls under it, the material is outside the
            regulation altogether. A whole-material test, not a per-nuclide
            filter.
        threshold: Index value at or above which the material fails, normally 1.
        source, url, retrieved, notes: Provenance.
    """

    name: str
    label: str
    units: str
    limits: dict[str, float]
    jurisdiction: str = ""
    default_limit: float | None = None
    unlimited: tuple[str, ...] = ()
    secular_equilibrium: dict[str, tuple[str, ...]] = field(default_factory=dict)
    secular_equilibrium_sec: dict[str, tuple[str, ...]] = field(default_factory=dict)
    limits_secular_equilibrium: dict[str, float] = field(default_factory=dict)
    metal_overrides: dict[str, float] = field(default_factory=dict)
    limits_per_gram: dict[str, float] = field(default_factory=dict)
    limits_upper: dict[str, float] = field(default_factory=dict)
    dynamic_rule: str | None = None
    min_half_life_scope: float | None = None
    threshold: float = 1.0
    source: str = ""
    url: str = ""
    retrieved: str = ""
    notes: str = ""

    def __post_init__(self):
        """Canonicalise and validate, so a hand-built set behaves like a loaded one.

        A set built in Python went through none of the normalisation the JSON
        loader applies, so ``{"Co-60": 0.1}`` silently matched nothing. Doing the
        work here means both paths cannot diverge.

        Limits are also required to be positive. A limit of zero means "no
        activity of this is permitted", the strictest possible value, but the
        sum of fractions cannot express that, and the evaluation used to skip it
        as though the nuclide had no limit at all, which is the opposite.
        Rejecting it here keeps that contradiction out of the data.
        """
        for field_name in (
            "limits",
            "metal_overrides",
            "limits_per_gram",
            "limits_upper",
            "limits_secular_equilibrium",
        ):
            mapping = getattr(self, field_name)
            object.__setattr__(self, field_name, _canonical_limits(mapping))

        object.__setattr__(
            self, "unlimited", tuple(_nuclide.normalise(n) for n in self.unlimited)
        )
        for field_name in ("secular_equilibrium", "secular_equilibrium_sec"):
            object.__setattr__(
                self,
                field_name,
                {
                    _nuclide.normalise(parent): tuple(
                        _nuclide.normalise(d) for d in daughters
                    )
                    for parent, daughters in getattr(self, field_name).items()
                },
            )

        if self.units not in ACTIVITY_UNITS:
            raise ValueError(
                f"{self.name}: units {self.units!r} is not one of "
                f"{', '.join(sorted(ACTIVITY_UNITS))}"
            )
        if not self.threshold > 0.0:
            raise ValueError(f"{self.name}: threshold must be positive, got {self.threshold}")
        if self.default_limit is not None and not self.default_limit > 0.0:
            raise ValueError(
                f"{self.name}: default_limit must be positive, got {self.default_limit}"
            )
        for field_name in ("limits", "metal_overrides", "limits_per_gram",
                           "limits_secular_equilibrium"):
            for nuclide_name, value in getattr(self, field_name).items():
                if not value > 0.0:
                    raise ValueError(
                        f"{self.name}: {field_name}[{nuclide_name}] is {value}, but a "
                        f"limit must be positive. A limit of zero cannot be expressed "
                        f"as a ratio and would be read as no limit at all."
                    )

    def daughters_of(self, parent: str) -> tuple[str, ...]:
        """Daughters whose activity this set's limit for ``parent`` already covers."""
        return self.secular_equilibrium.get(parent, ())

    @property
    def covered_nuclides(self) -> frozenset[str]:
        """Every nuclide this set names, whether limited or explicitly unlimited."""
        return frozenset(self.limits) | frozenset(self.unlimited) | frozenset(self.limits_per_gram)

    def __len__(self) -> int:
        return len(self.limits)

    def __repr__(self) -> str:
        return f"LimitSet({self.name!r}, {len(self.limits)} nuclides, {self.units})"


# ----------------------------------------------------------------------
# rules that depend on the material rather than only on the table
# ----------------------------------------------------------------------
def _nrc_short_lived_class_a(material, limit_set) -> dict[str, float]:
    """10 CFR 61.55 Table 2 row 1: all nuclides with a half-life under 5 years.

    The Class A column limits their total to 700 Ci/m3, so every such nuclide
    that does not have its own row takes that limit.
    """
    extra: dict[str, float] = {}
    for name in material.nuclides:
        half_life = material.decay_data.half_life(name)
        if half_life is not None and half_life < _FIVE_YEARS:
            extra[name] = 700.0
    return extra


#: Rules keyed by the name a limit set's ``dynamic_rule`` field refers to.
DYNAMIC_RULES: dict[str, Callable] = {
    "nrc_short_lived_class_a": _nrc_short_lived_class_a,
}


# ----------------------------------------------------------------------
# registry
# ----------------------------------------------------------------------
_REGISTRY: dict[str, LimitSet] = {}
_LOADED = False


def _canonical_limits(raw: Mapping[str, float]) -> dict[str, float]:
    out: dict[str, float] = {}
    for name, value in raw.items():
        # "sec" entries are stored under an explicit key so that the whole-chain
        # value and the tabulated-daughters value stay distinguishable.
        if name.endswith("_sec"):
            out[name] = float(value)
            continue
        out[_nuclide.normalise(name)] = float(value)
    return out


def _promote_sec_only_limits(limits: dict[str, float]) -> dict[str, float]:
    """Make a whole-chain value reachable when it is the only one published.

    A "sec" row is the value for a parent taken with its whole decay chain in
    secular equilibrium. It is stored under a "_sec" key so it stays distinct
    from a plain row, but that key can never match a nuclide in a material. For
    a nuclide whose only row is the "sec" one, keeping it there means the set
    applies no limit at all: UK_IRR17_natural would give unprocessed natural
    uranium a limit of nothing rather than the published 1 Bq/g. Where a plain
    row exists as well, it stays the limit and the "sec" value remains the
    variant selected when the chain is actually present.
    """
    promoted = dict(limits)
    for key, value in limits.items():
        if not key.endswith("_sec"):
            continue
        nuclide_name = key[: -len("_sec")]
        if nuclide_name not in promoted:
            promoted[nuclide_name] = value
    return promoted


def _from_dict(payload: Mapping) -> LimitSet:
    """Build a limit set from its JSON form. LimitSet.__post_init__ does the rest."""
    data = dict(payload)
    data["limits"] = _promote_sec_only_limits(_canonical_limits(data.get("limits", {})))
    known = {f for f in LimitSet.__dataclass_fields__}
    return LimitSet(**{k: v for k, v in data.items() if k in known})


def _load_all() -> None:
    global _LOADED
    if _LOADED:
        return
    for path in sorted(_LIMITS_DIR.glob("*.json")):
        payload = json.loads(path.read_text())
        for entry in payload.get("sets", []):
            limit_set = _from_dict(entry)
            _REGISTRY[limit_set.name] = limit_set
    _LOADED = True


def limit_sets(jurisdiction: str | None = None) -> tuple[str, ...]:
    """Names of the available limit sets, sorted.

    Args:
        jurisdiction: Restrict to one issuing authority, such as ``UK``.
    """
    _load_all()
    names = sorted(_REGISTRY)
    if jurisdiction is not None:
        names = [n for n in names if _REGISTRY[n].jurisdiction.lower() == jurisdiction.lower()]
    return tuple(names)


def get_limit_set(name: str | LimitSet) -> LimitSet:
    """Look up a limit set by name.

    Args:
        name: A registered name, or an already built `LimitSet`.

    Raises:
        KeyError: If no such set is registered, listing the ones that are.
    """
    if isinstance(name, LimitSet):
        return name
    _load_all()
    try:
        return _REGISTRY[name]
    except KeyError:
        raise KeyError(
            f"unknown limit set {name!r}. Available: {', '.join(sorted(_REGISTRY))}"
        ) from None


def register_limit_set(limit_set: LimitSet) -> None:
    """Add a limit set at runtime, for a site-specific or draft table.

    Args:
        limit_set: The set to register. Replaces any set of the same name,
            warning first if that name came from the shipped regulatory data,
            since shadowing a published table by accident would be hard to spot
            in a result that only records the name.
    """
    _load_all()
    if limit_set.name in _REGISTRY:
        warnings.warn(
            f"replacing the registered limit set {limit_set.name!r}, which came "
            f"from the shipped regulatory data. Results will report that name "
            f"while using the new table.",
            stacklevel=2,
        )
    _REGISTRY[limit_set.name] = limit_set
