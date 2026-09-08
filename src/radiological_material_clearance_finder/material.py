"""The material a clearance index is computed for.

A material here is just a nuclide inventory plus enough information to put it on
a per-gram or per-volume basis. It deliberately knows nothing about geometry,
temperature or cross sections, so it can be built from a dictionary typed by
hand as readily as from a depletion result.
"""
from __future__ import annotations

import warnings
from typing import Mapping

from . import nuclide as _nuclide
from .decay import AVOGADRO, BECQUEREL_PER_CURIE, DecayData, default_decay_data

__all__ = ["Material", "ACTIVITY_UNITS"]

#: Activity units understood by :meth:`Material.activity`.
ACTIVITY_UNITS = ("Bq", "Bq/g", "Bq/kg", "Bq/cm3", "Bq/m3", "Ci", "Ci/m3")

#: 1 atom/barn-cm is this many atoms per cubic centimetre.
_ATOMS_PER_BARN_CM = 1e24


class InsufficientDataError(ValueError):
    """Raised when a quantity cannot be derived from what the material was given."""


class Material:
    """A nuclide inventory to assess for clearance.

    The primary input is a mapping of nuclide name to atom count, which is what
    an activation or depletion calculation produces:

        >>> mat = Material({"Fe56": 8.4e22, "Co60": 1.2e12})

    Specific activity in Bq/g is scale invariant, so atom counts, atom
    densities and atom fractions all give the same answer for the Bq/g limit
    sets and no density is needed. Volumetric limit sets (Ci/m3) do need one,
    supplied as ``density`` or derived from atom densities.

    .. warning::
        The Bq/g denominator is the mass of **everything** in the mapping, so
        stable isotopes must be included. Passing only the radioactive nuclides
        of an activated steel gives a mass thousands of times too small and a
        specific activity thousands of times too high.

    Args:
        atoms: Nuclide name to atom count. Names may be spelled in any form
            :func:`~radiological_material_clearance_finder.nuclide.normalise`
            accepts.
        density: Mass density in g/cm3. Only needed for volumetric limit sets.
        volume: Volume in cm3. Only needed for total activity in Bq or Ci.
        name: A label carried through to results, for reporting.
        decay_data: Half-life and mass tables. Defaults to the vendored
            ENDF/B-VIII.0 and AME2020 values.
    """

    def __init__(
        self,
        atoms: Mapping[str, float] | None = None,
        *,
        density: float | None = None,
        volume: float | None = None,
        name: str = "",
        decay_data: DecayData | None = None,
        _specific_activities: Mapping[str, float] | None = None,
        _atoms_absolute: bool = True,
    ):
        self.decay_data = decay_data if decay_data is not None else default_decay_data()
        self.name = name
        self.volume = float(volume) if volume is not None else None
        self._density = float(density) if density is not None else None
        self._atoms_absolute = _atoms_absolute

        if atoms is None and _specific_activities is None:
            raise ValueError("a material needs either atom amounts or specific activities")

        self._atoms = self._canonicalise(atoms) if atoms is not None else None
        self._given_activities = (
            self._canonicalise(_specific_activities) if _specific_activities is not None else None
        )
        self._warn_if_inventory_looks_truncated()

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    @classmethod
    def from_atom_counts(
        cls,
        atoms: Mapping[str, float],
        *,
        density: float | None = None,
        volume: float | None = None,
        **kwargs,
    ) -> "Material":
        """Build from absolute atom counts. Identical to the constructor."""
        return cls(atoms, density=density, volume=volume, **kwargs)

    @classmethod
    def from_atom_densities(
        cls,
        densities: Mapping[str, float],
        *,
        volume: float | None = None,
        **kwargs,
    ) -> "Material":
        """Build from atom densities in atoms per barn-cm, OpenMC's unit.

        The mass density follows from the atom densities and the atomic masses,
        so no ``density`` argument is accepted or needed: supplying one could
        only contradict the inventory.

        Args:
            densities: Nuclide name to atoms/barn-cm.
            volume: Volume in cm3, only needed for total activity.
        """
        material = cls(densities, volume=volume, _atoms_absolute=False, **kwargs)
        material._density = sum(
            amount * _ATOMS_PER_BARN_CM * material.decay_data.atomic_mass(nuc) / AVOGADRO
            for nuc, amount in material._atoms.items()
        )
        return material

    @classmethod
    def from_masses(
        cls, masses: Mapping[str, float], *, density: float | None = None, **kwargs
    ) -> "Material":
        """Build from a mass in grams per nuclide.

        Args:
            masses: Nuclide name to mass in grams.
            density: Mass density in g/cm3, only needed for volumetric sets.
        """
        decay_data = kwargs.get("decay_data") or default_decay_data()
        atoms = {
            name: grams * AVOGADRO / decay_data.atomic_mass(name)
            for name, grams in masses.items()
        }
        return cls(atoms, density=density, **kwargs)

    @classmethod
    def from_mass_fractions(
        cls, fractions: Mapping[str, float], *, density: float | None = None, **kwargs
    ) -> "Material":
        """Build from mass fractions, which need not sum to one.

        Args:
            fractions: Nuclide name to mass fraction or weight percent.
            density: Mass density in g/cm3, only needed for volumetric sets.
        """
        decay_data = kwargs.get("decay_data") or default_decay_data()
        atoms = {
            name: fraction * AVOGADRO / decay_data.atomic_mass(name)
            for name, fraction in fractions.items()
        }
        material = cls(atoms, density=density, _atoms_absolute=False, **kwargs)
        return material

    @classmethod
    def from_specific_activities(
        cls,
        activities: Mapping[str, float],
        *,
        density: float | None = None,
        **kwargs,
    ) -> "Material":
        """Build from specific activities in Bq/g, as an assay reports them.

        No half-life or atomic mass data is used, since the specific activity is
        the direct input to the index. Volumetric limit sets then need an
        explicit ``density``, because activity alone does not imply a mass.

        Args:
            activities: Nuclide name to specific activity in Bq/g.
            density: Mass density in g/cm3, only needed for volumetric sets.
        """
        return cls(None, density=density, _specific_activities=activities, **kwargs)

    def _canonicalise(self, mapping: Mapping[str, float]) -> dict[str, float]:
        out: dict[str, float] = {}
        for raw_name, value in mapping.items():
            name = _nuclide.normalise(raw_name)
            if value < 0.0:
                raise ValueError(f"negative amount {value!r} for {name}")
            out[name] = out.get(name, 0.0) + float(value)
        return out

    def _warn_if_inventory_looks_truncated(self) -> None:
        """Warn when the inventory holds no stable nuclides.

        An activated material is overwhelmingly stable by mass. An inventory
        that is entirely radioactive is usually one filtered down to its
        radioactive part, which makes every specific activity too high by the
        ratio of true mass to retained mass. A pure source is legitimate, so
        this warns rather than raises.
        """
        names = self.nuclides
        if len(names) < 2 or self._given_activities is not None:
            return
        try:
            if any(not self.decay_data.is_radioactive(name) for name in names):
                return
        except Exception:
            return
        warnings.warn(
            "every nuclide in this material is radioactive, so the total mass "
            "used for Bq/g is only the radioactive mass. If this inventory was "
            "filtered to its radioactive nuclides, add the stable ones back or "
            "the specific activity will be far too high.",
            stacklevel=3,
        )

    # ------------------------------------------------------------------
    # inventory
    # ------------------------------------------------------------------
    @property
    def nuclides(self) -> tuple[str, ...]:
        """The nuclides present, in canonical form, sorted."""
        source = self._atoms if self._atoms is not None else self._given_activities
        return tuple(sorted(source))

    @property
    def atoms(self) -> dict[str, float]:
        """Atom amounts as supplied, keyed by canonical name."""
        if self._atoms is None:
            raise InsufficientDataError(
                "this material was built from specific activities, which do not "
                "determine atom counts"
            )
        return dict(self._atoms)

    @property
    def mass(self) -> float:
        """Total mass in grams.

        Raises:
            InsufficientDataError: If the material carries only relative
                amounts, such as atom densities or mass fractions, with no
                volume to scale them by.
        """
        if self._atoms is not None and self._atoms_absolute:
            return self._relative_mass()
        if self._density is not None and self.volume is not None:
            return self._density * self.volume
        raise InsufficientDataError(
            "total mass is unknown. This material holds relative amounts, so "
            "pass volume=, or build it with from_atom_counts or from_masses."
        )

    def _relative_mass(self) -> float:
        """Mass in grams of the amounts as supplied, whatever their scale."""
        return sum(
            amount * self.decay_data.atomic_mass(name) / AVOGADRO
            for name, amount in self._atoms.items()
        )

    @property
    def density(self) -> float:
        """Mass density in g/cm3.

        Raises:
            InsufficientDataError: If no density was supplied and none can be
                derived, naming what to pass.
        """
        if self._density is not None:
            return self._density
        if self._atoms is not None and self._atoms_absolute and self.volume:
            return self._relative_mass() / self.volume
        raise InsufficientDataError(
            "mass density is unknown, and volumetric limits (Ci/m3) need it. "
            "Pass density= in g/cm3, or build the material with "
            "from_atom_densities, which derives it."
        )

    # ------------------------------------------------------------------
    # activity
    # ------------------------------------------------------------------
    def specific_activity(self, by_nuclide: bool = False):
        """Specific activity in Bq/g.

        Scale invariant, so this works from atom counts, atom densities or atom
        fractions alike.

        Args:
            by_nuclide: Return a dict keyed by nuclide rather than the total.

        Returns:
            Bq/g as a float, or a dict of them.
        """
        if self._given_activities is not None:
            result = dict(self._given_activities)
        else:
            mass = self._relative_mass()
            if mass <= 0.0:
                raise InsufficientDataError("material has zero mass")
            result = {
                name: self.decay_data.decay_constant(name) * amount / mass
                for name, amount in self._atoms.items()
            }
        return result if by_nuclide else sum(result.values())

    def activity(self, units: str = "Bq/g", by_nuclide: bool = False):
        """Activity in the requested units.

        Args:
            units: One of :data:`ACTIVITY_UNITS`.
            by_nuclide: Return a dict keyed by nuclide rather than the total.

        Returns:
            The activity as a float, or a dict of them.

        Raises:
            ValueError: If the units are not recognised.
            InsufficientDataError: If the units need a density or volume the
                material does not have.
        """
        if units not in ACTIVITY_UNITS:
            raise ValueError(
                f"unknown activity units {units!r}, expected one of {', '.join(ACTIVITY_UNITS)}"
            )
        per_gram = self.specific_activity(by_nuclide=True)

        if units == "Bq/g":
            factor = 1.0
        elif units == "Bq/kg":
            factor = 1000.0
        elif units == "Bq/cm3":
            factor = self.density
        elif units == "Bq/m3":
            factor = self.density * 1e6
        elif units == "Ci/m3":
            factor = self.density * 1e6 / BECQUEREL_PER_CURIE
        elif units == "Bq":
            factor = self.mass
        else:  # Ci
            factor = self.mass / BECQUEREL_PER_CURIE

        result = {name: value * factor for name, value in per_gram.items()}
        return result if by_nuclide else sum(result.values())

    def __repr__(self) -> str:
        label = f"{self.name!r}, " if self.name else ""
        return f"Material({label}{len(self.nuclides)} nuclides)"
