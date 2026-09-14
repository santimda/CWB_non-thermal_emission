"""Shared plotting style and small helpers for project figures."""

import matplotlib
import numpy as np


# Set ``text.usetex`` to False if LaTeX is not installed locally.
nice_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 18,
    "font.size": 18,
    "axes.linewidth": 1.2,
    "axes.titlesize": 16,
    "legend.fontsize": 12,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "xtick.major.size": 5.5,
    "xtick.minor.size": 3,
    "ytick.major.size": 5.5,
    "ytick.minor.size": 3,
    "xtick.major.width": 1.3,
    "xtick.minor.width": 1.3,
    "ytick.major.width": 1.3,
    "ytick.minor.width": 1.3,
}

matplotlib.rcParams.update(nice_fonts)


def lighten_color(color, amount=0.4):
    """Mix an RGB color with white; ``amount=0`` leaves it unchanged."""
    color_array = np.array(matplotlib.colors.to_rgb(color))
    return tuple(1 - (1 - color_array) * (1 - amount))
