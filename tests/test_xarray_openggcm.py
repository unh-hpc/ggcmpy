from __future__ import annotations

import logging

import adios2py
import numpy as np
import pytest
import xarray as xr
from conftest import sample_3df, sample_iof, sample_py
from xarray_adios2 import Adios2Store

import ggcmpy
from ggcmpy import openggcm

logger = logging.getLogger(__name__)


def test_open_dataset_iof():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/coupling0001.iof.000030")
    assert set(ds.coords.keys()) >= set(sample_iof["coords"])
    assert set(ds.keys()) == set(sample_iof["data_vars"])
    assert ds.time == sample_iof["time"][0]  # type: ignore[index]
    assert ds.pot.sizes == sample_iof["sizes"]


def test_open_dataset_3df():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.3df.001200")
    assert set(ds.coords.keys()) >= set(sample_3df["coords"])
    assert set(ds.keys()) == set(sample_3df["data_vars"])
    assert ds.time == sample_3df["time"]
    assert ds.rr.sizes == sample_3df["sizes"]


def test_open_dataset_py():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.py_0.001200")
    assert set(ds.coords.keys()) >= set(sample_py["coords"])
    assert set(ds.keys()) == set(sample_py["data_vars"])
    assert ds.time == sample_py["time"]
    assert ds.rr.sizes == sample_py["sizes"]


def test_read_iof_jrrle_mfdataset():
    filenames = [
        f"{ggcmpy.sample_dir}/coupling0001.iof.{step:06d}" for step in [30, 60]
    ]
    logger.debug("before open")
    ds = xr.open_mfdataset(filenames)
    logger.debug("after open")
    assert set(ds.coords.keys()) >= set(sample_iof["coords"])
    assert set(ds.keys()) == set(sample_iof["data_vars"])
    assert np.all(ds.time == sample_iof["time"])
    logger.debug("before var0")
    var0 = ds[sample_iof["data_vars"][0]]  # type: ignore[index]
    var0 = var0.isel(time=0)
    logger.debug("before to_numpy()")
    assert np.isclose(np.sum(var0), 0.014970759)


@pytest.fixture
def time_dataset():
    amie_encoding = {"units": "time_array", "dtype": "int32"}

    return xr.Dataset(
        {
            "time_scalar": xr.Variable(
                (), np.datetime64("2021-01-01", "ns"), encoding=amie_encoding
            )
        },
        coords={
            "time": xr.Variable(
                "time",
                [np.datetime64("2020-01-01", "ns"), np.datetime64("2020-01-02")],
                encoding=amie_encoding,
            )
        },
    )


def test_time_write_read_by_step(tmp_path, time_dataset):
    filename = tmp_path / "test_ta.bp"
    time_dataset.attrs["step_dimension"] = "time"
    time_dataset.dump_to_store(
        Adios2Store.open(filename, "w"), encoder=openggcm.encode_openggcm
    )
    ds = xr.open_dataset(filename, decode_times=openggcm.AmieTimeArrayCoder())
    # FIXME, should mark time_scalar as invariant, so it doesn't get read back twice
    ds["time_scalar"] = ds.time_scalar.isel(time=0)
    assert time_dataset.equals(ds)


def test_time_write_read_one_step(tmp_path, time_dataset):
    filename = tmp_path / "test_ta_one.bp"
    time_dataset.dump_to_store(
        Adios2Store.open(filename, "w"), encoder=openggcm.encode_openggcm
    )
    with adios2py.File(filename, "r") as file:
        for step in file.steps:
            ds_read = xr.open_dataset(
                Adios2Store(step), decode_times=openggcm.AmieTimeArrayCoder()
            )
            assert time_dataset.equals(ds_read)
