#!/usr/bin/env python3
"""Sample and plot combined gamma-gamma transmission from production tables."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from python_scripts.plot_utils import nice_fonts  # noqa: F401  (applies project plot settings)


OUTPUT = ROOT / "output"
RESULTS = Path(__file__).resolve().parent / "results"

M_E_GG = 600
M_ANG_GG = 400
E_MIN = 50.0e6
E_MAX = 50.0e12
ENERGIES = (20.0e9, 100.0e9, 1.0e12)
VIEW_ANGLES_DEG = (1.0, 60.0)

AU = 1.495978707e13
R_SUN = 6.957e10

# Current reference geometry from system_parameters.f90 and orbit.f90.
DISTANCE_KPC = 2.4
PROJECTED_SEPARATION_MAS = 47.0
REFERENCE_ANGLE_DEG = 85.0
RSTAR_AU = 6.0 * R_SUN / AU
SEPARATION_AU = PROJECTED_SEPARATION_MAS * DISTANCE_KPC / np.sin(
    np.deg2rad(REFERENCE_ANGLE_DEG)
)

NX = 201
NY = 201
X_RANGE = (-0.49 * SEPARATION_AU, 1.5 * SEPARATION_AU)
Y_RANGE = (-0.99 * SEPARATION_AU, 1.0 * SEPARATION_AU)


def read_table(path):
    table = np.loadtxt(path)
    expected_shape = (M_E_GG, M_ANG_GG)
    if table.shape != expected_shape:
        raise ValueError(f"{path} does not have shape {expected_shape}")
    if not np.isfinite(table).all() or (table < 0.0).any():
        raise ValueError(f"{path} contains invalid opacity values")
    return table


def table_indices(energy, angle):
    energy_int = (E_MAX / E_MIN) ** (1.0 / (M_E_GG - 1))
    energy_index = int(np.log(energy / E_MIN) / np.log(energy_int))
    angle_int = np.pi / M_ANG_GG
    angle_min = angle_int / 2.0
    angle_index = (M_ANG_GG * np.abs(angle - angle_min) / np.pi).astype(int)
    return energy_index, np.clip(angle_index, 0, M_ANG_GG - 1)


def transmission_map(tau1, tau2, energy, view_angle_deg):
    x = np.linspace(*X_RANGE, NX)
    y = np.linspace(*Y_RANGE, NY)
    xx, yy = np.meshgrid(x, y)
    observer_angle = np.deg2rad(view_angle_deg)
    observer = np.array([np.cos(observer_angle), np.sin(observer_angle)])

    dx1, dy1 = xx, yy
    dx2, dy2 = xx - SEPARATION_AU, yy
    d1 = np.hypot(dx1, dy1)
    d2 = np.hypot(dx2, dy2)
    angle1 = np.arccos(np.clip((dx1 * observer[0] + dy1 * observer[1]) / d1, -1.0, 1.0))
    angle2 = np.arccos(np.clip((dx2 * observer[0] + dy2 * observer[1]) / d2, -1.0, 1.0))

    energy_index, angle_index1 = table_indices(energy, angle1)
    _, angle_index2 = table_indices(energy, angle2)
    tau = (RSTAR_AU / d1) * tau1[energy_index, angle_index1]
    tau += (RSTAR_AU / d2) * tau2[energy_index, angle_index2]
    transmission = np.exp(-tau)
    transmission[(d1 < RSTAR_AU) | (d2 < RSTAR_AU)] = np.nan
    return xx, yy, transmission


def save_map(path, xx, yy, transmission):
    values = np.column_stack((xx.ravel() / SEPARATION_AU,
                              yy.ravel() / SEPARATION_AU,
                              transmission.ravel()))
    np.savetxt(path, values, header="x_over_D y_over_D combined_transmission")


def energy_label(energy):
    if energy >= 1.0e12:
        return rf"{energy / 1.0e12:g}\,\mathrm{{TeV}}"
    return rf"{energy / 1.0e9:g}\,\mathrm{{GeV}}"


def main():
    tau1 = read_table(OUTPUT / "opac_gg1.dat")
    tau2 = read_table(OUTPUT / "opac_gg2.dat")
    RESULTS.mkdir(parents=True, exist_ok=True)

    maps = {}
    for energy in ENERGIES:
        for view_angle in VIEW_ANGLES_DEG:
            result = transmission_map(tau1, tau2, energy, view_angle)
            maps[(energy, view_angle)] = result
            label = f"{energy / 1.0e9:g}GeV_view{view_angle:g}deg"
            save_map(RESULTS / f"transmission_{label}.dat", *result)

    fig, axes = plt.subplots(
        len(ENERGIES), len(VIEW_ANGLES_DEG),
        figsize=(10, 13), sharex=True, sharey=True,
    )
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.05, top=0.98,
                        wspace=0.03, hspace=0.0)
    contour_levels = (0.1, 0.7, 0.9)
    contour_labels = {0.1: "90\%", 0.7: "30\%", 0.9: "10\%"}

    for row, energy in enumerate(ENERGIES):
        for column, view_angle in enumerate(VIEW_ANGLES_DEG):
            xx, yy, transmission = maps[(energy, view_angle)]
            axis = axes[row, column]
            image = axis.pcolormesh(
                xx / SEPARATION_AU, yy / SEPARATION_AU, transmission,
                shading="auto", vmin=0.0, vmax=1.0, cmap="magma",
            )
            contours = axis.contour(
                xx / SEPARATION_AU, yy / SEPARATION_AU, transmission,
                levels=contour_levels, colors="white", linewidths=1.0,
            )
            axis.clabel(contours, fmt=contour_labels, inline=True, fontsize=11)

            # The physical stellar radius is too small to see at this scale;
            # use 0.02 D circles consistently as visual markers.
            axis.add_patch(Circle((0.0, 0.0), 0.02, facecolor="cyan", edgecolor="darkturquoise", lw=1.0))
            axis.add_patch(Circle((1.0, 0.0), 0.02, facecolor="cyan", edgecolor="darkturquoise", lw=1.0))
            axis.text(0.0, 0.0, "1", ha="left", va="bottom", fontsize=14, color="darkturquoise")
            axis.text(1.0, 0.0, "2", ha="left", va="bottom", fontsize=14, color="darkturquoise")
            axis.set_aspect("equal")
            axis.text(
                0.03, 0.96,
                rf"$E_{{\rm ph}} = {energy_label(energy)}$, "
                rf"$\psi = {view_angle:g}^\circ$",
                transform=axis.transAxes, ha="left", va="top", color="black",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.65, pad=2.0),
            )
            if row == len(ENERGIES) - 1:
                axis.set_xlabel(r"$x / D$")
            else:
                axis.tick_params(labelbottom=False)
            if column == 0:
                axis.set_ylabel(r"$y / D$")

    colorbar_axis = fig.add_axes([0.90, 0.065, 0.02, 0.9])
    fig.colorbar(image, cax=colorbar_axis, label="Net transmission")
    fig.savefig(RESULTS / "gamma-gamma_transmission.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
