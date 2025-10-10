from __future__ import annotations

from pathlib import Path

from . import iono
from ._version import version as __version__

sample_dir = Path(__file__).parent / "sample"

__all__ = [
    "__version__",
    "iono",
    "sample_dir",
]
