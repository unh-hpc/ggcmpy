from __future__ import annotations

import pathlib

import pytest
from conftest import sample_iof

import ggcmpy
from ggcmpy.backends.jrrle import JrrleFile


def test_jrrle_file_non_existent():
    with pytest.raises(RuntimeError):
        JrrleFile("/non-existent")


def test_jrrle_file_open():
    JrrleFile(str(pathlib.Path(ggcmpy.sample_dir) / "coupling0001.iof.000030"))


def test_jrrle_file_open_path():
    JrrleFile(pathlib.Path(ggcmpy.sample_dir) / "coupling0001.iof.000030")


def test_jrrle_file_iter():
    with JrrleFile(pathlib.Path(ggcmpy.sample_dir) / "coupling0001.iof.000030") as file:
        for (var_name, _), ref_name in zip(file, sample_iof["data_vars"], strict=False):
            assert var_name == ref_name
            file.advance_one_line()
