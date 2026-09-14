#!/usr/bin/env python3
"""Plot Fortran Gaunt factors and the analytical approximation."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from python_scripts.plot_utils import nice_fonts  

RESULTS = Path(__file__).resolve().parent / "results"


def analytical_gaunt(frequency, temperature, charge):
    return 9.77 * (1.0 + 0.13 * np.log10(temperature**1.5 / charge / frequency))


def main():
    values = np.loadtxt(RESULTS / "gaunt_fortran.dat", comments="#")
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)

    for axis, temperature, charge, label in zip(
        axes,
        (0.3 * 86.08e3, 0.3 * 60.14e3),
        (1.0, 1.005),
        (r"$T_{w1}$, $Z_1$", r"$T_{w2}$, $Z_2$"),
    ):
        wind = np.isclose(values[:, 0], temperature)
        frequency = values[wind, 2]
        tabulated = values[wind, 3]
        axis.semilogx(frequency, tabulated, label="van Hoof+2014")
        axis.semilogx(
            frequency,
            analytical_gaunt(frequency, temperature, charge),
            "--",
            label="analytical approx.",
        )
        axis.set_title(
            rf"{label}: $T={temperature:.3g}\,\mathrm{{K}}$, $Z={charge:g}$"
        )
        axis.set_xlabel(r"Frequency [Hz]")
        axis.grid(True, which="both", alpha=0.25)

    axes[0].set_ylabel(r"Free-free Gaunt factor $g_{\rm ff}$")
    axes[0].legend()
    figure.tight_layout()
    figure.savefig(RESULTS / "gaunt_factors.pdf")


if __name__ == "__main__":
    main()
