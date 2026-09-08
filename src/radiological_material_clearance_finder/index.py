"""The clearance index itself: a sum of activity-to-limit ratios.

Every regulation here uses the same arithmetic. Each radionuclide's activity is
divided by its tabulated limit and the ratios are summed; a total below one
means the material meets the limits. The German regulation calls it the
Summenformel, the UK calls it the summation rule, the NRC calls it the sum of
fractions rule and Fetter calls the result a waste disposal rating. This module
implements it once.

What differs between regulations, and what :class:`ClearanceResult` therefore
records, is what happens to a nuclide that is *not* simply looked up: one whose
parent already accounts for it, one the table does not list, and one the
regulation places outside its scope entirely.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .decay import BECQUEREL_PER_CURIE
from .limits import DYNAMIC_RULES, LimitSet, get_limit_set, limit_sets
from .material import Material

__all__ = ["clearance_index", "clearance_indices", "clearable_routes", "ClearanceResult"]


@dataclass(frozen=True)
class ClearanceResult:
    """The outcome of assessing one material against one limit set.

    Attributes:
        limit_set: Name of the set assessed against.
        index: The sum of activity-to-limit ratios.
        threshold: The value the index must stay below, normally 1.
        units: Activity units the comparison was made in.
        by_nuclide: Each nuclide's contribution to the index, largest first.
        activities: Each nuclide's activity in ``units``.
        limits_used: The limit applied to each nuclide, after any metal,
            per-gram or dynamic adjustment.
        defaulted: Nuclides that took the set's catch-all limit because the
            table does not list them.
        excluded: Nuclide to the reason its activity was left out, which is
            always that a parent's limit already accounts for it.
        uncovered: Activity present with no limit and no catch-all, so absent
            from the index entirely. The number to check before trusting a
            comfortable index.
        unlimited: Nuclides the source explicitly places no limit on.
        out_of_scope: Whether the regulation excludes this material outright,
            which for the UK sets means every radionuclide present is shorter
            lived than 100 seconds.
        material_name: The material's label, carried through for reporting.
    """

    limit_set: str
    index: float
    threshold: float
    units: str
    by_nuclide: dict[str, float] = field(default_factory=dict)
    activities: dict[str, float] = field(default_factory=dict)
    limits_used: dict[str, float] = field(default_factory=dict)
    defaulted: tuple[str, ...] = ()
    excluded: dict[str, str] = field(default_factory=dict)
    uncovered: dict[str, float] = field(default_factory=dict)
    unlimited: tuple[str, ...] = ()
    out_of_scope: bool = False
    material_name: str = ""

    @property
    def clearable(self) -> bool:
        """Whether the material meets this set's limits."""
        return self.out_of_scope or self.index < self.threshold

    @property
    def uncovered_activity(self) -> float:
        """Total activity with no limit, in this result's units."""
        return sum(self.uncovered.values())

    @property
    def uncovered_fraction(self) -> float:
        """Share of total activity that falls outside the sum, from 0 to 1.

        A large value means the index understates the inventory: activity is
        present that no limit in this set applies to. Zero means every
        radionuclide present was either limited or explicitly unlimited.
        """
        total = sum(self.activities.values())
        return self.uncovered_activity / total if total > 0.0 else 0.0

    def dominant(self, count: int = 10) -> list[tuple[str, float]]:
        """The nuclides contributing most to the index.

        Args:
            count: How many to return.
        """
        return list(self.by_nuclide.items())[:count]

    def to_dict(self) -> dict:
        """A JSON-serialisable copy of the result."""
        return {
            "limit_set": self.limit_set,
            "material": self.material_name,
            "index": self.index,
            "threshold": self.threshold,
            "units": self.units,
            "clearable": self.clearable,
            "out_of_scope": self.out_of_scope,
            "by_nuclide": dict(self.by_nuclide),
            "activities": dict(self.activities),
            "limits_used": dict(self.limits_used),
            "defaulted": list(self.defaulted),
            "excluded": dict(self.excluded),
            "uncovered": dict(self.uncovered),
            "unlimited": list(self.unlimited),
            "uncovered_fraction": self.uncovered_fraction,
        }

    def __str__(self) -> str:
        header = f"{self.limit_set}: index {self.index:.4g} of {self.threshold:g}"
        if self.material_name:
            header = f"{self.material_name}  {header}"
        verdict = "CLEARABLE" if self.clearable else "NOT CLEARABLE"
        if self.out_of_scope:
            verdict = "OUT OF SCOPE (every radionuclide is short lived)"
        lines = [f"{header}  [{verdict}]"]
        if self.by_nuclide:
            lines.append(f"  {'nuclide':<10} {'activity':>12} {'limit':>12} {'share':>8}")
            for name, value in self.dominant(10):
                share = value / self.index if self.index else 0.0
                lines.append(
                    f"  {name:<10} {self.activities.get(name, 0.0):>12.4g} "
                    f"{self.limits_used.get(name, float('nan')):>12.4g} {share:>7.1%}"
                )
        if self.uncovered:
            lines.append(
                f"  {len(self.uncovered)} nuclide(s) with no limit, "
                f"{self.uncovered_fraction:.2%} of total activity"
            )
        return "\n".join(lines)


def _effective_limits(material: Material, limit_set: LimitSet, metal: bool) -> dict[str, float]:
    """Build the limit actually applied to each nuclide for this material."""
    limits = dict(limit_set.limits)

    if limit_set.dynamic_rule:
        rule = DYNAMIC_RULES[limit_set.dynamic_rule]
        for name, value in rule(material, limit_set).items():
            limits.setdefault(name, value)

    if metal:
        # Row based, not a modifier: the NRC lists "Ni-59 in activated metal" as
        # its own row, so an override may add a nuclide as well as replace one.
        limits.update(limit_set.metal_overrides)

    if limit_set.limits_per_gram:
        # nCi/g against a Ci/m3 comparison: 1 nCi/g is 1e-9 Ci/g, and
        # multiplying by g/cm3 then by 1e6 cm3/m3 gives Ci/m3.
        factor = 1e-9 * material.density * 1e6
        for name, value in limit_set.limits_per_gram.items():
            limits[name] = value * factor

    return limits


def _out_of_scope(material: Material, limit_set: LimitSet) -> bool:
    """Apply a whole-material scope test such as the UK 100 second rule.

    The regulation says a substance is not radioactive material where *none* of
    the radionuclides it contains has a half-life exceeding the threshold. That
    is a test on the material, so a short-lived nuclide cannot be dropped from
    the sum individually.
    """
    if limit_set.min_half_life_scope is None:
        return False
    half_lives = [
        material.decay_data.half_life(name)
        for name in material.nuclides
        if material.decay_data.is_radioactive(name)
    ]
    if not half_lives:
        return True
    return all(value <= limit_set.min_half_life_scope for value in half_lives)


def _resolve_equilibrium(material: Material, limit_set: LimitSet) -> tuple[dict, dict]:
    """Work out which daughters to exclude, and which parent limits that changes.

    A daughter is only excluded when its parent is actually present. Excluding
    it regardless would drop activity that arrived by some other route and is
    not in equilibrium with anything.

    Where the source lists a parent twice, plain and marked "+", the marked
    value applies only when the daughters really are there. StrlSchV gives
    Th-232 as 10 Bq/g plain and 0.01 Bq/g marked, so a thorium inventory with
    its chain present is held to a limit a thousand times stricter than freshly
    separated Th-232, which is what the regulation intends.

    Returns:
        ``(limit_overrides, excluded)``.
    """
    present = set(material.nuclides)
    excluded: dict[str, str] = {}
    overrides: dict[str, float] = {}
    for parent, daughters in limit_set.secular_equilibrium.items():
        if parent not in present:
            continue
        found = [d for d in daughters if d in present and d != parent]
        if not found:
            continue
        if parent in limit_set.limits_secular_equilibrium:
            overrides[parent] = limit_set.limits_secular_equilibrium[parent]
        for daughter in found:
            excluded[daughter] = (
                f"in secular equilibrium with {parent}, whose limit already "
                f"accounts for it"
            )
    return overrides, excluded


def clearance_index(
    material: Material,
    limit_set: str | LimitSet,
    *,
    metal: bool = False,
    apply_default_limit: bool = True,
    exclude_daughters: bool = True,
) -> ClearanceResult:
    """Assess a material against one set of clearance limits.

    Args:
        material: The inventory to assess.
        limit_set: A registered limit set name, or a :class:`LimitSet`.
        metal: Whether the material is activated metal, which changes some NRC
            limits and adds others.
        apply_default_limit: Whether to apply the set's catch-all limit to
            nuclides the table does not list. The UK regulations define one, so
            leaving this on is what the regulation says; turning it off shows
            what the listed nuclides alone contribute.
        exclude_daughters: Whether to leave out daughters whose parent's limit
            already covers them, per the regulation's own table. Turning this
            off double counts them, which is conservative but not what the
            regulation intends.

    Returns:
        A :class:`ClearanceResult` carrying the index and everything needed to
        judge it, including any activity that fell outside the sum.

    Raises:
        KeyError: If the limit set name is not registered.
        InsufficientDataError: If a volumetric set is used and the material has
            no density.
    """
    limit_set = get_limit_set(limit_set)
    activities = material.activity(units=limit_set.units, by_nuclide=True)
    limits = _effective_limits(material, limit_set, metal)
    if exclude_daughters:
        overrides, excluded = _resolve_equilibrium(material, limit_set)
        limits.update(overrides)
    else:
        excluded = {}
    unlimited_names = set(limit_set.unlimited)

    ratios: dict[str, float] = {}
    used: dict[str, float] = {}
    uncovered: dict[str, float] = {}
    defaulted: list[str] = []
    unlimited_present: list[str] = []

    for name, activity in activities.items():
        if activity <= 0.0:
            continue
        if name in excluded:
            continue
        if name in unlimited_names:
            unlimited_present.append(name)
            continue
        limit = limits.get(name)
        if limit is None:
            if apply_default_limit and limit_set.default_limit is not None:
                limit = limit_set.default_limit
                defaulted.append(name)
            else:
                uncovered[name] = activity
                continue
        if limit <= 0.0:
            uncovered[name] = activity
            continue
        ratios[name] = activity / limit
        used[name] = limit

    ordered = dict(sorted(ratios.items(), key=lambda item: -item[1]))
    return ClearanceResult(
        limit_set=limit_set.name,
        index=sum(ratios.values()),
        threshold=limit_set.threshold,
        units=limit_set.units,
        by_nuclide=ordered,
        activities=dict(sorted(activities.items(), key=lambda item: -item[1])),
        limits_used=used,
        defaulted=tuple(sorted(defaulted)),
        excluded=excluded,
        uncovered=dict(sorted(uncovered.items(), key=lambda item: -item[1])),
        unlimited=tuple(sorted(unlimited_present)),
        out_of_scope=_out_of_scope(material, limit_set),
        material_name=material.name,
    )


def clearance_indices(
    material: Material,
    names: "list[str] | tuple[str, ...] | None" = None,
    **kwargs,
) -> dict[str, ClearanceResult]:
    """Assess a material against many limit sets at once.

    Sets whose units the material cannot supply, such as a volumetric set for a
    material with no density, are skipped rather than raising, so one missing
    density does not hide every result that does not need it.

    Args:
        material: The inventory to assess.
        names: Limit set names, defaulting to every registered set.
        **kwargs: Passed through to :func:`clearance_index`.

    Returns:
        Results keyed by limit set name.
    """
    from .material import InsufficientDataError

    results: dict[str, ClearanceResult] = {}
    for name in names if names is not None else limit_sets():
        try:
            results[name] = clearance_index(material, name, **kwargs)
        except InsufficientDataError:
            continue
    return results


def clearable_routes(
    material: Material,
    names: "list[str] | tuple[str, ...] | None" = None,
    **kwargs,
) -> list[str]:
    """The limit sets this material already meets, best margin first.

    Args:
        material: The inventory to assess.
        names: Limit set names, defaulting to every registered set.
        **kwargs: Passed through to :func:`clearance_index`.

    Returns:
        Names of the sets whose index is below their threshold.
    """
    results = clearance_indices(material, names, **kwargs)
    passing = [r for r in results.values() if r.clearable]
    passing.sort(key=lambda r: r.index / r.threshold if r.threshold else r.index)
    return [r.limit_set for r in passing]
