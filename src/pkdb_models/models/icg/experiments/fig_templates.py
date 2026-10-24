"""Templates of the matplotlib figures of the indocyanine green experiments."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment


def scan_color(length: int, reverse: bool = False) -> ListedColormap:
    """Get a colormap from blue to red (red to blue for `reverse`).

    Args:
        length: number of colors.
        reverse: the colormap from red to blue.
    """
    red_steps = np.linspace(0, 1, length)
    blue_steps = np.linspace(1, 0, length)
    if reverse:
        red_steps = np.linspace(1, 0, length)
        blue_steps = np.linspace(0, 1, length)
    colors = [(red_steps[k], 0, blue_steps[k]) for k in range(length)]
    return ListedColormap(colors)


def base_figure() -> tuple[Figure, Axes]:
    """Get a single matplotlib figure with the fonts of the experiments."""
    f, ax = plt.subplots(figsize=(6, 6))
    f.subplots_adjust(wspace=0.3, hspace=0.3)
    ax.tick_params(axis="x", labelsize=ICGSimulationExperiment.tick_font_size)
    ax.tick_params(axis="y", labelsize=ICGSimulationExperiment.tick_font_size)
    return f, ax
