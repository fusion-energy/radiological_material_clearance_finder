"""Waste classification, where the answer is a category rather than an index.

The US and UK both sort waste into named classes rather than reporting a single
ratio. The NRC's classes come out of the same sum-of-fractions machinery as a
clearance index, applied to several tables and combined by rules in 10 CFR
61.55. The UK's categories instead compare gross alpha and beta or gamma
activity against concentration thresholds, so they need no per-nuclide table at
all and are computed directly.
"""
from __future__ import annotations

from dataclasses import dataclass

from .index import clearance_index
from .material import Material

__all__ = [
    "nrc_waste_class",
    "uk_waste_category",
    "UKWasteCategory",
    "alpha_activity",
    "beta_gamma_activity",
    "UK_LLW_ALPHA_BQ_PER_G",
    "UK_LLW_BETA_GAMMA_BQ_PER_G",
    "UK_HIGH_VOLUME_VLLW_BQ_PER_G",
    "UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G",
]

# The 2007 policy states these per tonne. One GBq/tonne is 1e9 Bq per 1e6 g,
# which is 1e3 Bq/g, and one MBq/tonne is 1 Bq/g.
#: LLW upper bound on alpha activity, from 4 GBq/te.
UK_LLW_ALPHA_BQ_PER_G = 4.0e3
#: LLW upper bound on beta and gamma activity, from 12 GBq/te.
UK_LLW_BETA_GAMMA_BQ_PER_G = 12.0e3
#: High volume VLLW upper bound on total activity, from 4 MBq/te.
UK_HIGH_VOLUME_VLLW_BQ_PER_G = 4.0
#: High volume VLLW upper bound on tritium and carbon-14 together, from
#: 40 MBq/te. The 2007 policy gives them a shared allowance in both the low
#: volume and the high volume categories, so they are summed rather than tritium
#: being singled out.
UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G = 40.0


def nrc_waste_class(material: Material, *, metal: bool = False) -> str:
    """Classify a material for near-surface disposal under 10 CFR 61.55.

    Applies the sum of fractions rule to Table 1 and Table 2 and combines the
    results with the rules in paragraphs 61.55(a)(3) to (a)(7).

    Args:
        material: The inventory to classify. Needs a density, since the NRC
            limits are volumetric.
        metal: Whether the material is activated metal, which the regulation
            treats with its own rows.

    Returns:
        ``"Class A"``, ``"Class B"``, ``"Class C"`` or ``"GTCC"`` for greater
        than Class C.
    """
    table1 = clearance_index(material, "NRC_long", metal=metal).index
    table2 = [
        clearance_index(material, name, metal=metal).index
        for name in ("NRC_short_A", "NRC_short_B", "NRC_short_C")
    ]

    def by_table2(column_a: float, column_b: float, column_c: float) -> str:
        if column_a < 1.0:
            return "Class A"
        if column_b < 1.0:
            return "Class B"
        if column_c < 1.0:
            return "Class C"
        return "GTCC"

    table1_present = table1 > 0.0
    table2_present = any(value > 0.0 for value in table2)

    if table1_present and table2_present:
        # 61.55(a)(5)
        if table1 < 0.1:
            return by_table2(*table2)
        if table1 < 1.0:
            return "Class C" if table2[2] < 1.0 else "GTCC"
        return "GTCC"
    if table1_present:
        # 61.55(a)(3)
        if table1 < 0.1:
            return "Class A"
        return "Class C" if table1 < 1.0 else "GTCC"
    if table2_present:
        # 61.55(a)(4)
        return by_table2(*table2)
    # 61.55(a)(6)
    return "Class A"


def alpha_activity(material: Material, units: str = "Bq/g") -> float:
    """Activity from alpha emission, weighted by each nuclide's alpha branch.

    Bi-212 branches 35.94 percent alpha, so it contributes that share of its
    activity here and the rest to :func:`beta_gamma_activity`.

    Args:
        material: The inventory.
        units: Any unit :meth:`Material.activity` accepts.
    """
    activities = material.activity(units=units, by_nuclide=True)
    return sum(
        value * material.decay_data.alpha_fraction(name)
        for name, value in activities.items()
    )


def beta_gamma_activity(material: Material, units: str = "Bq/g") -> float:
    """Activity from everything that is not alpha emission.

    Args:
        material: The inventory.
        units: Any unit :meth:`Material.activity` accepts.
    """
    activities = material.activity(units=units, by_nuclide=True)
    return sum(
        value * (1.0 - material.decay_data.alpha_fraction(name))
        for name, value in activities.items()
    )


@dataclass(frozen=True)
class UKWasteCategory:
    """The UK waste category of a material, with the numbers behind it.

    Attributes:
        category: ``"VLLW"``, ``"LLW"`` or ``"ILW"``.
        alpha: Alpha activity in Bq/g.
        beta_gamma: Beta and gamma activity in Bq/g.
        tritium_and_c14: Tritium plus carbon-14 activity in Bq/g, which share
            their own VLLW allowance.
        total: Total activity in Bq/g.
        reason: Why this category and not the one below it.
    """

    category: str
    alpha: float
    beta_gamma: float
    tritium_and_c14: float
    total: float
    reason: str

    def __str__(self) -> str:
        return f"{self.category}: {self.reason}"


def uk_waste_category(material: Material) -> UKWasteCategory:
    """Classify a material against the UK solid low level waste categories.

    Thresholds are from *Policy for the Long Term Management of Solid Low Level
    Radioactive Waste in the United Kingdom* (Defra and the Devolved
    Administrations, 2007), which defines LLW as not exceeding 4 GBq/te alpha or
    12 GBq/te beta and gamma, and high volume VLLW as not exceeding 4 MBq/te
    total activity, with tritium allowed to 40 MBq/te.

    Two limits of this classification are worth stating rather than hiding.
    The low volume VLLW category ("dustbin disposal") is defined per 0.1 cubic
    metre and per item, so it is a property of a consignment rather than of a
    material and is not decided here. The boundary above which waste becomes
    high level rather than intermediate is a thermal one, about 2 kW/m3, and
    needs decay heat rather than activity, so this returns ILW for everything
    above LLW.

    Args:
        material: The inventory to classify.

    Returns:
        A :class:`UKWasteCategory` holding the category and the activities it
        was decided on.
    """
    per_nuclide = material.activity(units="Bq/g", by_nuclide=True)
    alpha = alpha_activity(material)
    beta_gamma = beta_gamma_activity(material)
    tritium_and_c14 = per_nuclide.get("H3", 0.0) + per_nuclide.get("C14", 0.0)
    total = sum(per_nuclide.values())
    remainder = total - tritium_and_c14

    if alpha > UK_LLW_ALPHA_BQ_PER_G:
        return UKWasteCategory(
            "ILW", alpha, beta_gamma, tritium_and_c14, total,
            f"alpha activity {alpha:.3g} Bq/g exceeds the LLW limit of "
            f"{UK_LLW_ALPHA_BQ_PER_G:g} Bq/g (4 GBq/te)",
        )
    if beta_gamma > UK_LLW_BETA_GAMMA_BQ_PER_G:
        return UKWasteCategory(
            "ILW", alpha, beta_gamma, tritium_and_c14, total,
            f"beta and gamma activity {beta_gamma:.3g} Bq/g exceeds the LLW limit "
            f"of {UK_LLW_BETA_GAMMA_BQ_PER_G:g} Bq/g (12 GBq/te)",
        )
    if remainder > UK_HIGH_VOLUME_VLLW_BQ_PER_G:
        return UKWasteCategory(
            "LLW", alpha, beta_gamma, tritium_and_c14, total,
            f"total activity {remainder:.3g} Bq/g excluding tritium and carbon-14 "
            f"exceeds the high volume VLLW limit of "
            f"{UK_HIGH_VOLUME_VLLW_BQ_PER_G:g} Bq/g (4 MBq/te), and it is within "
            f"both LLW limits",
        )
    if tritium_and_c14 > UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G:
        return UKWasteCategory(
            "LLW", alpha, beta_gamma, tritium_and_c14, total,
            f"tritium and carbon-14 activity {tritium_and_c14:.3g} Bq/g exceeds "
            f"their shared high volume VLLW allowance of "
            f"{UK_HIGH_VOLUME_VLLW_TRITIUM_AND_C14_BQ_PER_G:g} Bq/g (40 MBq/te)",
        )
    return UKWasteCategory(
        "VLLW", alpha, beta_gamma, tritium_and_c14, total,
        f"total activity {remainder:.3g} Bq/g excluding tritium and carbon-14 is "
        f"within the high volume VLLW limit of {UK_HIGH_VOLUME_VLLW_BQ_PER_G:g} "
        f"Bq/g (4 MBq/te)",
    )
