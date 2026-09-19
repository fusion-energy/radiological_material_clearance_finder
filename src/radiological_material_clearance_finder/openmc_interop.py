"""Getting an inventory out of OpenMC without depending on it.

OpenMC is imported only inside these functions and only when they are called, so
installing this package pulls in nothing, and everything except these two
functions works with no OpenMC present.
"""
from __future__ import annotations

from .material import Material

__all__ = ["from_openmc_material", "from_depletion_results"]


def from_openmc_material(openmc_material, *, decay_data=None) -> Material:
    """Convert an ``openmc.Material`` into a `Material`.

    Reads the atom densities in atoms per barn-cm, which fixes the mass density
    too, so the result supports both the Bq/g and the Ci/m3 limit sets.

    Args:
        openmc_material: An ``openmc.Material``. Its nuclides must be expanded,
            which ``get_nuclide_atom_densities`` does.
        decay_data: Half-life and mass tables to use instead of the defaults.

    Returns:
        The equivalent `Material`, carrying the OpenMC material's name
        and volume where it has them.
    """
    densities = openmc_material.get_nuclide_atom_densities()
    return Material.from_atom_densities(
        {name: float(value) for name, value in densities.items()},
        volume=getattr(openmc_material, "volume", None),
        name=getattr(openmc_material, "name", "") or "",
        decay_data=decay_data,
    )


def from_depletion_results(results, material_id, *, decay_data=None) -> dict[float, Material]:
    """Build a cooling time series from an OpenMC depletion results file.

    The result is in the form
    `time_to_clear`
    expects, so a depletion can be taken straight through to a clearance date.

    Args:
        results: An ``openmc.deplete.Results`` object.
        material_id: The depleted material's id, as a string or int.
        decay_data: Half-life and mass tables to use instead of the defaults.

    Returns:
        Time in seconds since the start of the results to the material at that
        time.
    """
    times = results.get_times(time_units="s")
    wanted = str(material_id)
    series: dict[float, Material] = {}
    for step, time in enumerate(times):
        found = None
        for candidate in results.export_to_materials(step):
            if str(candidate.id) == wanted:
                found = candidate
                break
        if found is None:
            raise KeyError(
                f"material {material_id!r} is not in the depletion results at step {step}"
            )
        series[float(time)] = from_openmc_material(found, decay_data=decay_data)
    return series
