"""Assess an activated steel against every clearance route.

    python examples/activated_steel.py

Uses a made-up inventory so it runs without any nuclear data. Replace it with a
depletion result and the rest works unchanged.
"""
from radiological_material_clearance_finder import (
    Material,
    clearance_index,
    clearance_indices,
    nrc_waste_class,
    time_to_clear,
    uk_waste_category,
)
from radiological_material_clearance_finder.decay import half_life

YEAR = 365.25 * 86400

#: Atom densities in atoms per barn-cm, as OpenMC reports them. The stable
#: isotopes matter: they are the mass the specific activity is divided by.
STEEL = {
    "Fe54": 4.90e-3, "Fe56": 7.71e-2, "Fe57": 1.78e-3, "Fe58": 2.37e-4,
    "Cr52": 1.57e-2, "Ni58": 7.90e-3, "Mn55": 1.40e-3, "Mo98": 3.00e-4,
    "Co59": 8.00e-5,
    # Activation products.
    "Co60": 2.1e-9, "Fe55": 8.0e-8, "Ni63": 3.4e-9, "Mn54": 9.0e-10,
    "H3": 2.2e-10, "Nb94": 4.4e-14, "Tc99": 8.0e-14, "C14": 5.5e-12,
}


def decayed(inventory, seconds):
    """Decay the activation products, leaving the stable bulk alone.

    A stand-in for a depletion result. It ignores ingrowth, which is why this
    package does not do decay itself.
    """
    out = {}
    for name, amount in inventory.items():
        life = half_life(name)
        out[name] = amount if life is None else amount * 0.5 ** (seconds / life)
    return out


def main() -> None:
    steel = Material.from_atom_densities(STEEL, name="activated steel")
    print(f"{steel!r}  density {steel.density:.3f} g/cm3  "
          f"{steel.activity('Bq/g'):.4g} Bq/g\n")

    print("Clearance index by route, at shutdown")
    print(f"  {'route':32} {'index':>12}  verdict")
    for name, result in sorted(
        clearance_indices(steel).items(), key=lambda item: item[1].index
    ):
        verdict = "clearable" if result.clearable else "not clearable"
        print(f"  {name:32} {result.index:>12.4g}  {verdict}")

    print(f"\nWaste classification")
    print(f"  NRC 10 CFR 61.55: {nrc_waste_class(steel, metal=True)}")
    print(f"  UK category:      {uk_waste_category(steel)}")

    print("\nWhat drives the German unrestricted clearance index")
    result = clearance_index(steel, "StrlSchV_unrestricted")
    print(result)
    if result.uncovered:
        print(f"  note: {result.uncovered_fraction:.2%} of the activity has no limit "
              f"in this set")

    print("\nCooling time to clear")
    times = [0.0, 1 * YEAR, 5 * YEAR, 10 * YEAR, 30 * YEAR, 100 * YEAR, 300 * YEAR]
    series = {t: Material.from_atom_densities(decayed(STEEL, t)) for t in times}
    for route in ("StrlSchV_metal_recycling", "StrlSchV_unrestricted",
                  "UK_EPR16_out_of_scope"):
        when = time_to_clear(series, route)
        if when is None:
            print(f"  {route:28} not within {times[-1] / YEAR:.0f} years")
        else:
            print(f"  {route:28} {when / YEAR:.1f} years")


if __name__ == "__main__":
    main()
