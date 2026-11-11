from __future__ import annotations

from pathlib import Path

from . import colormaps, iono
from ._version import version as __version__

sample_dir = Path(__file__).parent / "sample"

__all__ = [
    "__version__",
    "colormaps",
    "iono",
    "sample_dir",
]
