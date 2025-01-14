from __future__ import annotations

import hashlib
import logging

import numpy as np
import xarray as xr
from conftest import sample_3df, sample_iof, sample_py

import ggcmpy

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
    assert hashlib.sha256(var0.to_numpy().tobytes()).hexdigest().startswith("6f8b81d9")
