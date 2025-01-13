from __future__ import annotations

import numpy as np
import xarray as xr

import ggcmpy

sample_iof = {
    "coords": {"lats", "longs", "time"},
    "data_vars": {
        "pot",
        "sigp",
        "sigh",
        "rrio",
        "ppio",
        "ttio",
        "prec_e_e0_1",
        "prec_e_fe_1",
        "prec_e_e0_2",
        "prec_e_fe_2",
        "fac_dyn",
        "fac_tot",
        "vdown",
        "ctaut",
        "delphi",
        "tau",
        "pacurr",
        "xjh",
        "delbp",
        "epio",
        "ctiot",
        "delbr",
        "ctaup",
        "etio",
        "delbt",
        "cpolt",
        "cpolp",
        "ctiop",
    },
    "time": np.asarray(
        ["1967-01-01T00:00:30.150", "1967-01-01T00:01:00.101"], dtype=np.datetime64
    ),
}

sample_3df = {
    "coords": {"x", "y", "z"},
    "data_vars": {
        "rr",
        "pp",
        "vx",
        "vy",
        "vz",
        "bx",
        "by",
        "bz",
        "xjx",
        "xjy",
        "xjz",
        "xtra1",
        "xtra2",
        "resis",
    },
    "time": {},
}


def test_open_dataset_iof():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/coupling0001.iof.000030")
    assert set(ds.coords.keys()) == sample_iof["coords"]
    assert set(ds.keys()) == sample_iof["data_vars"]
    assert ds.time == sample_iof["time"][0]  # type: ignore[index]


def test_open_dataset_3df():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.3df.001200")
    assert set(ds.coords.keys()) == sample_3df["coords"]
    assert set(ds.keys()) == sample_3df["data_vars"]
    # assert ds.time == sample_3df["time"]


def test_read_iof_jrrle_mfdataset():
    filenames = [
        f"{ggcmpy.sample_dir}/coupling0001.iof.{step:06d}" for step in [30, 60]
    ]
    ds = xr.open_mfdataset(filenames)
    assert set(ds.coords.keys()) == sample_iof["coords"]
    assert set(ds.keys()) == sample_iof["data_vars"]
    assert np.all(ds.time == sample_iof["time"])
