from __future__ import annotations

import pathlib

import pytest

import ggcmpy
from ggcmpy.backends.jrrle import JrrleFile


def test_jrrle_file_non_existent():
    with pytest.raises(RuntimeError):
        JrrleFile("/non-existent")


def test_jrrle_file_open():
    JrrleFile(str(pathlib.Path(ggcmpy.sample_dir) / "coupling0001.iof.000030"))
