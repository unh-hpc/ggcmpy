from __future__ import annotations

import pathlib

import numpy as np
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
        for var_name, ref_name in zip(file.vars, sample_iof["data_vars"], strict=False):
            assert var_name == ref_name


def test_jrrle_file_inquire():
    vars: list[str] = sample_iof["data_vars"]  # type: ignore[assignment]
    with JrrleFile(pathlib.Path(ggcmpy.sample_dir) / "coupling0001.iof.000030") as file:
        for var_name in vars[:3]:
            file._inquire(var_name)
        assert vars[1] in file.fields_seen
        assert vars[3] not in file.fields_seen
        file._inquire(vars[3])
        assert vars[3] in file.fields_seen


def test_jrrle_file_read_field():
    vars: list[str] = sample_iof["data_vars"]  # type: ignore[assignment]
    with JrrleFile(pathlib.Path(ggcmpy.sample_dir) / "coupling0001.iof.000030") as file:
        meta, var0 = file.read_field(vars[0])
        assert np.isclose(np.sum(var0), 0.014970759)
        meta, var1 = file.read_field(vars[1])
        assert np.isclose(np.sum(var1), 36606.758)
        meta, var2 = file.read_field(vars[2])
        assert np.isclose(np.sum(var2), 7.3803894e-05)

        meta, var1 = file.read_field(vars[1])
        assert np.isclose(np.sum(var1), 36606.758)

        meta, var3 = file.read_field(vars[3])
        assert np.isclose(np.sum(var3), 38640.805)

        with pytest.raises(KeyError):
            meta, var = file.read_field("nowhere")
