#!/usr/bin/env python3
"""Plot the available absorbed, beam-convolved synthetic radio maps."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse

from python_scripts.plot_utils import nice_fonts  # noqa: F401


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output"
RESULTS = Path(__file__).resolve().parent / "results"

SIGMA_X = 5.6
SIGMA_Y = 11.3
D = 47.0

def map_paths():
    """Return available convolved maps ordered by their filename frequency."""
    paths = list(OUTPUT.glob("WCR_conv_*Gabs.dat"))
    if not paths:
        raise FileNotFoundError("No absorbed convolved maps found in output/")

    def frequency(path):
        label = path.stem.removeprefix("WCR_conv_").removesuffix("Gabs")
        return int(label)

    return sorted(((frequency(path), path) for path in paths), key=lambda item: item[0])


def read_map(frequency, path):
    """Read one convolved absorbed map."""
    data = np.loadtxt(path)
    ny = np.count_nonzero(data[:, 0] == data[0, 0])
    nx = len(data) // ny

    x = data[::ny, 0]
    y = data[:ny, 1]
    absorbed = data[:, 2].reshape(nx, ny).T

    # Keep the orientation used by the former map notebook.
    return frequency, y, -x, absorbed


def make_plot(maps):
    palette = LinearSegmentedColormap.from_list(
        "map_gray", [(1.0, 1.0, 1.0), (0.15, 0.15, 0.15)]
    )
    figure, axes = plt.subplots(
        1,
        len(maps),
        figsize=(5 * len(maps), 4.8),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)
    images = []
    contour_ticks = []

    for axis, (frequency, y, x, value) in zip(axes, maps):
        maximum = np.nanmax(value)
        fill_levels = np.linspace(0.0, maximum, 101)
        contour_levels = maximum * np.array([0.03, 0.1, 0.3, 0.5, 0.7, 0.9])
        contour_ticks.append(contour_levels)
        image = axis.contourf(
            x, y, value, levels=fill_levels, cmap=palette, extend="max"
        )
        images.append(image)
        axis.contour(x, y, value, levels=contour_levels, colors="black", linewidths=0.8)
        axis.set_title(rf"$\nu={frequency}\,\mathrm{{GHz}}$")
        axis.set_xlim(-90, 30)
        axis.set_ylim(-70, 70)
        axis.add_patch(
            Ellipse(
                (-78, -58),
                2.0 * SIGMA_X,
                2.0 * SIGMA_Y,
                color="grey",
                alpha=0.8,
                zorder=4,
            )
        )
        axis.scatter(
            [-D, 0.0],
            [0.0, 0.0],
            s=180,
            marker="*",
            color="DeepSkyBlue",
            zorder=5,
        )
        axis.annotate("WC", (-D + 3.0, -8.0), color="DeepSkyBlue", zorder=6)
        axis.annotate("WN", (3.0, -8.0), color="DeepSkyBlue", zorder=6)
        axis.set_aspect("equal")

    axes[0].set_ylabel("$y$ [mas]")
    axes[1].set_xlabel("$x$ [mas]")
    for axis, image, ticks in zip(axes, images, contour_ticks):
        colorbar = figure.colorbar(image, ax=axis, label=r"$S_\nu$ [mJy beam$^{-1}$]")
        colorbar.set_ticks(ticks)
        colorbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])
    figure.suptitle("Synthetic radio maps")
    figure.savefig(RESULTS / "maps_absorbed.pdf", bbox_inches="tight")
    plt.close(figure)


def main():
    RESULTS.mkdir(exist_ok=True)
    maps = [read_map(frequency, path) for frequency, path in map_paths()]
    make_plot(maps)


if __name__ == "__main__":
    main()
