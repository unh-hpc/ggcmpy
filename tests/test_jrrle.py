from __future__ import annotations

import hashlib
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


def test_jrrle_file_inquire():
    vars: list[str] = sample_iof["data_vars"]  # type: ignore[assignment]
    with JrrleFile(pathlib.Path(ggcmpy.sample_dir) / "coupling0001.iof.000030") as file:
        for var_name in vars[:3]:
            file.inquire(var_name)
        assert vars[1] in file.fields_seen
        assert vars[3] not in file.fields_seen
        file.inquire(vars[3])
        assert vars[3] in file.fields_seen


def test_jrrle_file_read_field():
    vars: list[str] = sample_iof["data_vars"]  # type: ignore[assignment]
    with JrrleFile(pathlib.Path(ggcmpy.sample_dir) / "coupling0001.iof.000030") as file:
        meta, var0 = file.read_field(vars[0], 2)
        assert hashlib.sha256(var0.tobytes()).hexdigest().startswith("6f8b81d9")
        meta, var1 = file.read_field(vars[1], 2)
        assert hashlib.sha256(var1.tobytes()).hexdigest().startswith("572fed56")
        meta, var2 = file.read_field(vars[2], 2)
        assert hashlib.sha256(var2.tobytes()).hexdigest().startswith("cd58c4af")

        meta, var1 = file.read_field(vars[1], 2)
        assert hashlib.sha256(var1.tobytes()).hexdigest().startswith("572fed56")

        meta, var3 = file.read_field(vars[3], 2)
        assert hashlib.sha256(var3.tobytes()).hexdigest().startswith("7960c540")
