"""
Custom matplotlib colormaps.

Importing :mod:`ggcmpy` registers these colormaps with matplotlib, so they can
be used by name, e.g. ``plt.imshow(data, cmap="openggcm")``.
"""

from __future__ import annotations

import matplotlib as mpl  # type: ignore[import-not-found]
import matplotlib.colors as mcolors  # type: ignore[import-not-found]

OPENGGCM_CMAP_NAME = "openggcm"

_openggcm_colors = [
    (0.0, "#ffffff"),  # 0 / 4
    (0.25, "#87ceeb"),  # 1 / 4
    (0.5, "#90ff90"),  # 2 / 4
    (0.75, "#daa520"),  # 3 / 4
    (1.0, "#8b0000"),  # 4 / 4
]


def _register() -> None:
    """Register all colormaps defined here with matplotlib (idempotent)."""
    if OPENGGCM_CMAP_NAME in mpl.colormaps:
        return
    cmap = mcolors.LinearSegmentedColormap.from_list(
        OPENGGCM_CMAP_NAME, _openggcm_colors
    )
    mpl.colormaps.register(cmap, name=OPENGGCM_CMAP_NAME)


_register()
