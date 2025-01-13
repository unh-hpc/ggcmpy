from __future__ import annotations

import numpy as np
import xarray as xr
from conftest import sample_3df, sample_iof, sample_py

import ggcmpy


def test_open_dataset_iof():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/coupling0001.iof.000030")
    assert set(ds.coords.keys()) == sample_iof["coords"]
    assert set(ds.keys()) == sample_iof["data_vars"]
    assert ds.time == sample_iof["time"][0]  # type: ignore[index]
    assert ds.pot.sizes == sample_iof["sizes"]


def test_open_dataset_3df():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.3df.001200")
    assert set(ds.coords.keys()) == sample_3df["coords"]
    assert set(ds.keys()) == sample_3df["data_vars"]
    assert ds.time == sample_3df["time"]
    assert ds.rr.sizes == sample_3df["sizes"]


def test_open_dataset_py():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.py_0.001200")
    assert set(ds.coords.keys()) == sample_py["coords"]
    assert set(ds.keys()) == sample_py["data_vars"]
    assert ds.time == sample_py["time"]
    assert ds.rr.sizes == sample_py["sizes"]


def test_read_iof_jrrle_mfdataset():
    filenames = [
        f"{ggcmpy.sample_dir}/coupling0001.iof.{step:06d}" for step in [30, 60]
    ]
    ds = xr.open_mfdataset(filenames)
    assert set(ds.coords.keys()) == sample_iof["coords"]
    assert set(ds.keys()) == sample_iof["data_vars"]
    assert np.all(ds.time == sample_iof["time"])
