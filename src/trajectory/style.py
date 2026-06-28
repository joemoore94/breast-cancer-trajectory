"""Shared figure style for all breast-cancer-trajectory plots."""

from __future__ import annotations

import matplotlib as mpl
import seaborn as sns


def apply_style(scatter: bool = False) -> None:
    """Set a consistent matplotlib/seaborn style for all figures.

    Parameters
    ----------
    scatter:
        Pass ``True`` for scatter/spatial figures (UMAP, diffusion maps) where
        a grid background is distracting.
    """
    style = "white" if scatter else "whitegrid"
    sns.set_theme(style=style, context="talk", font_scale=1.6)
    mpl.rcParams.update({
        "axes.titlesize": 26,
        "axes.labelsize": 24,
        "axes.titlepad": 14,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "legend.fontsize": 20,
        "legend.title_fontsize": 20,
        "figure.titlesize": 28,
        "font.size": 20,
    })
