from __future__ import annotations

import os

from ._version import version as __version__

sample_dir = os.path.join(os.path.dirname(__file__), "sample")

__all__ = [
    "__version__",
    "sample_dir",
]
