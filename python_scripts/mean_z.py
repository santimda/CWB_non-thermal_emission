#!/usr/bin/env python3
"""Estimate composition factors used by the non-thermal emission model.

The X/Y/Z calculation treats all metals as one representative species. 
The Wilms et al. (2000) case uses the ISM tabulated particle fractions for 
specific elements. 

The epsilon factors follow the expressions given in Padovani et al. (2018), Appendix A; 
the paper also reports the values for the Wilms et al. (2000) abundances."""

import argparse
import math
from dataclasses import dataclass


PADOVANI_WILMS = {"br": 2.24, "pp": 2.17}


@dataclass(frozen=True)
class Composition:
    name: str
    X: float
    Y: float
    Z: float
    metal_mass: float
    metal_charge: float
    note: str


def composition_from_elements(name, elements, metal_mass, metal_charge, note):
    """Convert abundances relative to He into bulk mass fractions."""
    total_mass = sum(A * abundance for _, A, _, abundance in elements)
    X = next(A * abundance for element, A, _, abundance in elements if element == "H") / total_mass
    Y = next(A * abundance for element, A, _, abundance in elements if element == "He") / total_mass
    return Composition(name, X, Y, 1.0 - X - Y, metal_mass, metal_charge, note)


# Zhekov et al. (2025), abundances by number relative to He.
WN_ELEMENTS = (
    ("H", 1, 1, 0.067),
    ("He", 4, 2, 1.0),
    ("C", 12, 6, 1.28e-4),
    ("N", 14, 7, 1.29e-3),
    ("O", 16, 8, 2.92e-4),
    ("Ne", 20, 10, 6.57e-4),
    ("Mg", 24, 12, 2.19e-4),
    ("Si", 28, 14, 2.16e-4),
    ("S", 32, 16, 5.11e-5),
    ("Ar", 40, 18, 2.0e-5),
    ("Ca", 40, 20, 2.0e-5),
    ("Fe", 56, 26, 1.28e-4),
)

WC_ELEMENTS = (
    ("H", 1, 1, 0.0),
    ("He", 4, 2, 1.0),
    ("C", 12, 6, 0.4),
    ("N", 14, 7, 0.0),
    ("O", 16, 8, 0.194),
    ("Ne", 20, 10, 1.86e-2),
    ("Mg", 24, 12, 2.72e-3),
    ("Si", 28, 14, 6.84e-4),
    ("S", 32, 16, 1.52e-4),
    ("Ar", 40, 18, 2.0e-5),
    ("Ca", 40, 20, 2.0e-5),
    ("Fe", 56, 26, 3.82e-4),
)


COMPOSITIONS = {
	# Solar composition values from Lodders+2025: https://ui.adsabs.harvard.edu/abs/2025SSRv..221...23L/abstract
	# 
    "solar": Composition(
        "solar X/Y/Z",
        0.739,
        0.245,
        0.016,
        17.0,
        17.0 / 2.02,
        "representative average mass and charge for metals",
    ),
    "WN": composition_from_elements(
        "WN star with Zhekov et al. (2025) abundances",
        WN_ELEMENTS,
        19.5,
        10.5,
        "rounded effective metal mass and charge",
    ),
    "WC": composition_from_elements(
        "WC star with Zhekov et al. (2025) abundances",
        WC_ELEMENTS,
        14.0,
        7.0,
        "rounded effective metal mass and charge",
    ),
    "WN_approx": Composition(
        "approximate WN star X/Y/Z",
        0.016,
        0.969,
        0.015,
        19.5,
        10.5,
        "Zhekov et al. (2025)-derived bulk fractions with rounded effective metal values",
    ),
    "WC_approx": Composition(
        "approximate WC star X/Y/Z",
        0.001,
        0.323,
        0.676,
        14.0,
        7.0,
        "Zhekov et al. (2025)-derived bulk fractions with rounded effective metal values",
    ),
}


# Wilms et al. (2000) particle fractions taken from Table A.1 in Padovani+ 2018.
WILMS_ELEMENTS = (
    ("H2", 2, 2, 0.835),
    ("He", 4, 2, 0.163),
    ("C", 12, 6, 4.01e-4),
    ("N", 14, 7, 1.27e-4),
    ("O", 16, 8, 8.19e-4),
    ("Ne", 20, 10, 1.46e-4),
    ("Na", 23, 11, 2.41e-6),
    ("Mg", 24, 12, 4.19e-5),
    ("Al", 27, 13, 3.57e-6),
    ("Si", 28, 14, 3.11e-5),
    ("P", 31, 15, 4.39e-7),
    ("S", 32, 16, 2.05e-5),
    ("Cl", 35, 17, 2.21e-7),
    ("Ar", 40, 18, 4.29e-6),
    ("Ca", 40, 20, 2.64e-6),
    ("Ti", 48, 22, 1.08e-7),
    ("Cr", 52, 24, 5.41e-7),
    ("Mn", 55, 25, 3.66e-7),
    ("Fe", 56, 26, 4.49e-5),
    ("Co", 59, 27, 1.39e-7),
    ("Ni", 59, 28, 1.87e-6),
)


def validate_xyz(X: float, Y: float, Z: float) -> None:
    if min(X, Y, Z) < 0.0:
        raise ValueError("mass fractions must be non-negative")
    if not math.isclose(X + Y + Z, 1.0, abs_tol=1.0e-12):
        raise ValueError(f"X + Y + Z must equal one, got {X + Y + Z:.8f}")


def xyz_factors(composition: Composition, full_ionization: bool = False) -> dict[str, float]:
    validate_xyz(composition.X, composition.Y, composition.Z)
    X, Y, Z = composition.X, composition.Y, composition.Z
    A = composition.metal_mass
    charge = composition.metal_charge

    if full_ionization:
        # For the WCR, H is atomic.  These expressions are the same
        # approximation used in system_parameters.f90 (A ~= 2 Z for metals).
        ion_per_mass = X + Y / 4.0 + Z / A
        mu_t = 4.0 / (3.0 + 5.0 * X - Z)
        mu_e = 2.0 / (1.0 + X)
        mu_i = 1.0 / (1.0 / mu_t - 1.0 / mu_e)
        zq_rms = math.sqrt(
            (X + (Y / 4.0) * 2.0**2 + (Z / A) * charge**2) / ion_per_mass
        )
        f_h = X / ion_per_mass
        f_he_ion = (Y / 4.0) / ion_per_mass
        f_metals_ion = (Z / A) / ion_per_mass
        return {
            "mu_t": mu_t,
            "mu_e": mu_e,
            "mu_i": mu_i,
            "zq_rms": zq_rms,
            "z_br": f_h + 3.0 * f_he_ion + f_metals_ion * charge * (charge + 1.0) / 2.0,
            "z_pp": f_h + f_he_ion * 4.0**0.79 + f_metals_ion * A**0.79,
        }

    # Padovani et al. use particle fractions for a molecular target.
    number_per_mass = X / 2.0 + Y / 4.0 + Z / A
    f_h2 = (X / 2.0) / number_per_mass
    f_he = (Y / 4.0) / number_per_mass
    f_metals = (Z / A) / number_per_mass
    return {
        "z_br": 2.0 * f_h2 + 3.0 * f_he + f_metals * charge * (charge + 1.0) / 2.0,
        "z_pp": 2.0 * f_h2 + f_he * 4.0**0.79 + f_metals * A**0.79,
    }


def wilms_factors(full_ionization: bool = False) -> dict[str, float]:
    total = sum(item[3] for item in WILMS_ELEMENTS)
    abundances = [(name, A, charge, fraction / total) for name, A, charge, fraction in WILMS_ELEMENTS]
    if full_ionization:
        # Split H2 into two atomic H ions for the WCR quantities, while
        # retaining H2 in the Padovani factors above.
        # Use the original table directly for mass, ion, and electron sums.
        mass = 2.0 * abundances[0][3] + sum(fraction * A for _, A, _, fraction in abundances[1:])
        ions = 2.0 * abundances[0][3] + sum(fraction for _, _, _, fraction in abundances[1:])
        electrons = 2.0 * abundances[0][3] + sum(fraction * charge for _, _, charge, fraction in abundances[1:])
        zq_rms = math.sqrt(
            (2.0 * abundances[0][3] + sum(fraction * charge**2 for _, _, charge, fraction in abundances[1:]))
            / ions
        )
        return {
            "mu_t": mass / (ions + electrons),
            "mu_e": mass / electrons,
            "mu_i": mass / ions,
            "zq_rms": zq_rms,
            "z_br": (2.0 * abundances[0][3] + sum(fraction * charge * (charge + 1.0) / 2.0 for _, _, charge, fraction in abundances[1:])) / ions,
            "z_pp": (2.0 * abundances[0][3] + sum(fraction * A**0.79 for _, A, _, fraction in abundances[1:])) / ions,
        }

    return {
        "z_br": 2.0 * abundances[0][3]
        + sum(fraction * charge * (charge + 1.0) / 2.0 for _, _, charge, fraction in abundances[1:]),
        "z_pp": 2.0 * abundances[0][3]
        + sum(fraction * A**0.79 for _, A, _, fraction in abundances[1:]),
    }


def element_factors(elements) -> dict[str, float]:
    """Calculate fully ionized factors directly from number abundances."""
    total = sum(abundance for _, _, _, abundance in elements)
    abundances = [
        (name, A, charge, abundance / total)
        for name, A, charge, abundance in elements
    ]
    return {
        "z_br": sum(fraction * charge * (charge + 1.0) / 2.0 for _, _, charge, fraction in abundances),
        "z_pp": sum(fraction * A**0.79 for _, A, _, fraction in abundances),
    }


def print_report(name: str, values: dict[str, float], reference: bool = False, full_ionization: bool = False) -> None:
    print(f"\n{name}")
    keys = ("z_br", "z_pp")
    if full_ionization:
        keys = ("mu_t", "mu_i", "mu_e", "zq_rms") + keys
    for key in keys:
        print(f"  {key:12s} = {values[key]:.3g}")
    if reference and not full_ionization:
        for key in ("br", "pp"):
            calculated = values[f"z_{key}"]
            expected = PADOVANI_WILMS[key]
            print(
                f"  Padovani z_{key:2s}: "
                f"{expected:.3g} "
                f"(ratio calculated/reference = {calculated / expected:.2f})"
            )


def print_element_report(name: str, elements) -> None:
    values = element_factors(elements)
    print(f"\n{name}, element-by-element, fully ionized")
    print(f"  {'z_br':12s} = {values['z_br']:.3g}")
    print(f"  {'z_pp':12s} = {values['z_pp']:.3g}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--composition",
        choices=("solar", "WN", "WC", "WN_approx", "WC_approx", "wilms", "all"),
        default="all",
        help="composition to report (default: all)",
    )
    parser.add_argument(
        "--full-ionization",
        action="store_true",
        help="also report fully ionized WCR quantities",
    )
    args = parser.parse_args()

    selected = ("solar", "WN", "WC", "WN_approx", "WC_approx") if args.composition == "all" else (args.composition,)
    for composition_name in selected:
        composition = COMPOSITIONS[composition_name]
        if args.full_ionization:
            print(f"\nComposition: {composition.name} ({composition.note})")
            print_report(
                "X/Y/Z approximation (fully ionized WCR)",
                xyz_factors(composition, full_ionization=True),
                full_ionization=True,
            )
        elif composition_name == "WN":
            print_element_report("WN star with Zhekov et al. (2025) abundances", WN_ELEMENTS)
        elif composition_name == "WC":
            print_element_report("WC star with Zhekov et al. (2025) abundances", WC_ELEMENTS)

    if args.composition in ("wilms", "all"):
        print_report(
            "Wilms et al. (2000)"
            + (" (fully ionized WCR)" if args.full_ionization else ""),
            wilms_factors(args.full_ionization),
            reference=True,
            full_ionization=args.full_ionization,
        )


if __name__ == "__main__":
    main()
