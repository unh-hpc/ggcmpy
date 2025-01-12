from __future__ import annotations

import importlib.metadata

import ggcmpy as m


def test_version():
    assert importlib.metadata.version("ggcmpy") == m.__version__
