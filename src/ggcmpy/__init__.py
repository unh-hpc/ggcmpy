import os

__version__ = "0.0.1"
sample_dir = os.path.join(os.path.dirname(__file__), "sample")

__all__ = [
    "__version__",
    "sample_dir",
]
