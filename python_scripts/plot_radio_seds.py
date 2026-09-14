#!/usr/bin/env python3
"""Plot the Apep radio SED Mdot scan and report observed spectral indices."""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from python_scripts.plot_utils import lighten_color


DATA = Path(__file__).resolve().parent / "data_radio_apep.dat"
RESULTS = Path(__file__).resolve().parent / "results"

MJDY = 1.0e-26
KPC = 3.0856e21
DISTANCE = 2.4 * KPC
H = 6.6260755e-27
EV = 1.602176e-12
DISTANCE_FACTOR = 4.0 * np.pi * DISTANCE**2

CALIBRATION_ERRORS = {
    "ATCA": 0.10,
    "RACS": 0.10,
    "uGMRT": 0.10,
    "ALMA-B3": 0.07,
    "ALMA-B6": 0.10,
    "ALMA-B7": 0.10,
    "ALMA-ACA-B6": 0.10,
    "ALMA-ACA-B7": 0.10,
    "MOST": 0.10,
}


def read_sed(path):
    """Read an emitted SED and return frequency in GHz and flux in mJy."""
    sed = np.loadtxt(path)
    energy = sed[:, 0]
    luminosity_absorbed = sed[:, 2]
    frequency = energy * EV / H * 1.0e-9
    flux = luminosity_absorbed * H / DISTANCE_FACTOR / MJDY
    return frequency, flux


def power_law(frequency, normalization, alpha, reference_frequency=90.0):
    return normalization * (frequency / reference_frequency) ** alpha


def fit_power_law(frequency, flux, errors=None, reference_frequency=90.0, initial=None):
    """Fit a power law using the same weighted ``curve_fit`` method as the notebook."""
    if initial is None:
        initial = (float(np.median(flux)), 0.0)
    return curve_fit(
        lambda nu, normalization, alpha: power_law(
            nu, normalization, alpha, reference_frequency
        ),
        frequency,
        flux,
        p0=initial,
        sigma=errors,
    )


def print_observed_spectral_indices(data):
    """Print the fitted spectral index for each ALMA band used in the paper."""
    bands = (
        ("ALMA-B3", 86.0, (-0.1, 10.0)),
        ("ALMA-ACA-B6", 225.0, (0.5, 7.0)),
        ("ALMA-ACA-B7", 320.0, (0.5, 7.0)),
    )
    for instrument, reference_frequency, initial in bands:
        band = data[data["instrument"] == instrument]
        parameters, covariance = fit_power_law(
            band["nu"].to_numpy(),
            band["flux"].to_numpy(),
            errors=band["flux_err"].to_numpy(),
            reference_frequency=reference_frequency,
            initial=initial,
        )
        uncertainty = np.sqrt(np.diag(covariance))
        print(f"{instrument}: alpha = {parameters[1]:.2f} +/- {uncertainty[1]:.2f}")


def load_observations():
    data = pd.read_csv(DATA, sep=r"\s+", comment="#")
    data["f_sys"] = data["instrument"].map(CALIBRATION_ERRORS)
    data["flux_err"] = np.sqrt(
        data["flux_err_sta"] ** 2 + (data["f_sys"] * data["flux"]) ** 2
    )
    return data


def load_model_seds():
    syn_paths = sorted(RESULTS.glob("rad_syn_WCR_M*.dat"))
    ff_paths = sorted(RESULTS.glob("windemiff_M*.dat"))
    if len(syn_paths) != 6 or len(ff_paths) != 6:
        raise RuntimeError("Expected six synchrotron and six free-free Mdot SEDs")

    mdots = [float(path.stem.split("_M")[-1]) for path in syn_paths]
    seds_syn = [read_sed(path) for path in syn_paths]
    seds_ff = []
    for path in ff_paths:
        sed = np.loadtxt(path)
        sed[:, 0] *= 1.0e-9
        seds_ff.append(sed)

    seds_total = []
    for sed_syn, sed_ff in zip(seds_syn, seds_ff):
        parameters, _ = fit_power_law(
            sed_ff[:, 0], sed_ff[:, 1], reference_frequency=1.0, initial=(0.1, 0.6)
        )
        ff_at_syn_frequency = power_law(
            sed_syn[0], *parameters, reference_frequency=1.0
        )
        seds_total.append(sed_syn[1] + ff_at_syn_frequency)

    return mdots, seds_syn, seds_ff, seds_total


def make_figure(data, mdots, seds_syn, seds_ff, seds_total):
    xmin, xmax = 4.0e-1, 9.0e2
    ymin, ymax = 3.0, 3.0e2

    figure, axis = plt.subplots(1)
    axis.set_xlabel(r"$\nu$ [GHz]")
    axis.set_ylabel(r"$S_\nu$ [mJy]")
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)

    count = len(mdots)
    colors = plt.cm.copper_r(np.linspace(0.05, 0.8, count))
    color_map = matplotlib.colormaps["copper_r"].resampled(count)
    norm = matplotlib.colors.BoundaryNorm(np.arange(count + 1) + 0.5, count)
    scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=color_map)

    for index in range(count):
        axis.plot(seds_syn[index][0], seds_total[index], color=colors[index], linewidth=2)
        axis.plot(seds_ff[index][:, 0], seds_ff[index][:, 1], ":", color=colors[index], alpha=0.7)
    axis.plot(seds_syn[-1][0], seds_syn[-1][1], "--", color=colors[-1], alpha=0.6)

    styles = {
        "ATCA": dict(marker="v", color="blue", whitening=0.5, alpha=0.9, ms=5, label="ATCA"),
        "RACS": dict(marker="v", color="magenta", whitening=0.5, alpha=0.9, ms=5, label="RACS"),
        "uGMRT": dict(marker="^", color="g", whitening=0.5, alpha=0.9, ms=5, label="uGMRT"),
        "ALMA-B3": dict(marker="D", color="red", whitening=0.5, alpha=0.9, ms=5, label="ALMA"),
        "ALMA-ACA-B6": dict(marker="D", color="red", whitening=0.5, alpha=0.9, ms=5, label=None),
        "ALMA-ACA-B7": dict(marker="D", color="red", whitening=0.5, alpha=0.9, ms=5, label=None),
        "ALMA-B6": dict(marker="D", color="red", whitening=0.5, alpha=0.9, ms=4, label=None),
        "ALMA-B7": dict(marker="D", color="red", whitening=0.5, alpha=0.9, ms=4, label=None),
        "MOST": dict(marker="o", color="y", whitening=0.5, alpha=0.9, ms=5, label="MOST"),
    }
    for instrument in np.unique(data["instrument"]):
        band = data[data["instrument"] == instrument]
        style = styles[instrument]
        axis.errorbar(
            x=band["nu"],
            y=band["flux"],
            xerr=band["bw"] / 2,
            yerr=band["flux_err"],
            fmt=style["marker"],
            ms=style["ms"],
            color=style["color"],
            ecolor=style["color"],
            markerfacecolor=(*lighten_color(style["color"], style["whitening"]), style["alpha"]),
            markeredgecolor=style["color"],
            alpha=style["alpha"],
            capsize=1.8,
            linestyle="none",
            label=style["label"],
        )

    axis.legend(loc="best")
    axis.loglog()
    axis.xaxis.labelpad = -0.1
    colorbar = figure.colorbar(
        scalar_map,
        ticks=np.arange(1, count + 1),
        label=r"$\dot{M}_\mathrm{WN}$ [$10^{-5}$ M$_\odot$ yr$^{-1}$]",
        pad=0.01,
        ax=axis,
    )
    colorbar.ax.set_yticklabels(mdots)
    axis.tick_params(axis="both", which="both", right=True, left=True, top=True, direction="in")
    figure.savefig(RESULTS / "SEDs_Mdot.pdf", bbox_inches="tight", pad_inches=0)
    plt.close(figure)


def main():
    data = load_observations()
    print_observed_spectral_indices(data)
    mdots, seds_syn, seds_ff, seds_total = load_model_seds()
    make_figure(data, mdots, seds_syn, seds_ff, seds_total)
    print(f"Wrote {RESULTS / 'SEDs_Mdot.pdf'}")


if __name__ == "__main__":
    main()
