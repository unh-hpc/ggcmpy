from __future__ import annotations

import numpy as np
import xarray as xr

import ggcmpy

iof_coords = {"lats", "longs", "time"}
iof_vars = {
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
}
iof_times = np.asarray(
    ["1967-01-01T00:00:30.150", "1967-01-01T00:01:00.101"], dtype=np.datetime64
)


def test_open_dataset():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/coupling0001.iof.000030")
    assert set(ds.coords.keys()) == iof_coords
    assert set(ds.keys()) == iof_vars
    assert ds.time == iof_times[0]


def test_read_iof_jrrle_mfdataset():
    filenames = [
        f"{ggcmpy.sample_dir}/coupling0001.iof.{step:06d}" for step in [30, 60]
    ]
    ds = xr.open_mfdataset(filenames)
    assert set(ds.coords.keys()) == iof_coords
    assert set(ds.keys()) == iof_vars
    assert np.all(ds.time.to_numpy() == iof_times)
