#!/usr/bin/env python3
"""Plot free-free transmission maps written by the Fortran validation driver."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from python_scripts.plot_utils import nice_fonts  


FREQUENCIES = (0.3e9, 1.0e9, 5.0e9)
VIEW_ANGLES_DEG = (75.0, 105.0)
NX = 100
NY = 100
X_RANGE = (-1.49, 2.5)
Y_RANGE = (-1.99, 2.0)


def read_map(path):
    values = np.loadtxt(path)
    expected_shape = (NX * NY, 2 + len(FREQUENCIES))
    if values.shape != expected_shape:
        raise ValueError(f"{path} has shape {values.shape}, expected {expected_shape}")
    if not np.isfinite(values).all():
        raise ValueError(f"{path} contains non-finite values")
    x = values[:, 0].reshape(NY, NX)
    y = values[:, 1].reshape(NY, NX)
    transmission = values[:, 2:].T.reshape(len(FREQUENCIES), NY, NX)
    return x, y, transmission


def read_geometry(output):
    orbit = np.atleast_2d(np.loadtxt(output / "orbit.dat"))
    separation = orbit[0, 2]
    cd = np.atleast_2d(np.loadtxt(output / "CD.dat"))
    cd_approx = np.atleast_2d(np.loadtxt(output / "CD_approx.dat"))

    cd_x = cd[:, 0] / separation
    cd_y = cd[:, 1] / separation
    approx_positive = cd_approx[cd_approx[:, 1] >= 0.0]
    approx_x = approx_positive[:, 0] / separation
    approx_y = approx_positive[:, 1] / separation

    def mirrored_branches(x_values, y_values):
        x_values = np.concatenate((x_values, [np.nan], x_values))
        y_values = np.concatenate((y_values, [np.nan], -y_values))
        return x_values, y_values

    cd_x, cd_y = mirrored_branches(cd_x, cd_y)
    approx_x, approx_y = mirrored_branches(approx_x, approx_y)
    coefficients = np.atleast_1d(np.loadtxt(output / "a_CD.dat"))
    return separation, cd_x, cd_y, approx_x, approx_y, coefficients


def lower_hyperbola_mask(x, y, separation, coefficients, view_angle):
    a_cd, b_cd = coefficients
    x_v = a_cd / separation
    b_over_d = b_cd / separation
    tangent = abs(np.tan(np.deg2rad(view_angle)))
    if view_angle < 90.0:
        x_line = x_v + y / tangent
    else:
        x_line = x_v - y / tangent
    x_hyperbola = x_v * np.sqrt(1.0 + (y / b_over_d)**2)
    return (y < 0.0) & (x_line < x) & (x < x_hyperbola)


def frequency_label(frequency):
    return rf"{frequency / 1.0e9:g}\,\mathrm{{GHz}}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()

    output = args.work_dir / "output"
    results = Path(__file__).resolve().parent / "results"
    results.mkdir(parents=True, exist_ok=True)

    maps = {
        view_angle: read_map(output / f"ffa_map_psi{view_angle:g}.dat")
        for view_angle in VIEW_ANGLES_DEG
    }
    separation, cd_x, cd_y, approx_x, approx_y, coefficients = read_geometry(output)

    fig, axes = plt.subplots(
        len(FREQUENCIES), len(VIEW_ANGLES_DEG),
        figsize=(10, 13), sharex=True, sharey=True,
    )
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.05, top=0.98,
                        wspace=0.03, hspace=0.0)
    contour_levels = (0.1, 0.7, 0.9)
    contour_labels = {0.1: "10\%", 0.7: "70\%", 0.9: "90\%"}

    for row, frequency in enumerate(FREQUENCIES):
        for column, view_angle in enumerate(VIEW_ANGLES_DEG):
            xx, yy, transmissions = maps[view_angle]
            axis = axes[row, column]
            transmission = transmissions[row]
            mask = lower_hyperbola_mask(xx, yy, separation, coefficients, view_angle)
            transmission = transmission.copy()
            transmission[mask] = np.nan
            image = axis.pcolormesh(
                xx, yy, transmission, shading="auto", vmin=0.0, vmax=1.0, cmap="magma",
            )
            contours = axis.contour(
                xx, yy, transmission,
                levels=contour_levels, colors="white", linewidths=1.0,
            )
            axis.clabel(contours, fmt=contour_labels, inline=True, fontsize=11)
            axis.plot(cd_x, cd_y, color="darkgreen", linewidth=2.0, zorder=4)
            axis.plot(approx_x, approx_y, color="limegreen", linewidth=2.5,
                      linestyle="--", zorder=4)

            # The physical stellar radius is too small to see at this scale;
            # use 0.02 D circles consistently as visual markers.
            axis.add_patch(Circle((0.0, 0.0), 0.02, facecolor="cyan",
                                  edgecolor="darkturquoise", lw=1.0, zorder=5))
            axis.add_patch(Circle((1.0, 0.0), 0.02, facecolor="cyan",
                                  edgecolor="darkturquoise", lw=1.0, zorder=5))
            axis.text(0.0, 0.0, "1", ha="left", va="bottom", fontsize=14,
                      color="darkturquoise", zorder=6)
            axis.text(1.0, 0.0, "2", ha="left", va="bottom", fontsize=14,
                      color="darkturquoise", zorder=6)
            axis.set_aspect("equal")
            axis.set_xlim(*X_RANGE)
            axis.set_ylim(*Y_RANGE)
            axis.text(
                0.03, 0.96,
                rf"$\nu = {frequency_label(frequency)}$, "
                rf"$\psi = {view_angle:g}^\circ$",
                transform=axis.transAxes, ha="left", va="top", color="black",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.65, pad=2.0),
            )
            if row == len(FREQUENCIES) - 1:
                axis.set_xlabel(r"$x / D$")
            else:
                axis.tick_params(labelbottom=False)
            if column == 0:
                axis.set_ylabel(r"$y / D$")

    colorbar_axis = fig.add_axes([0.90, 0.065, 0.02, 0.9])
    fig.colorbar(image, cax=colorbar_axis, label="Net transmission")
    fig.savefig(results / "free-free_transmission.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
