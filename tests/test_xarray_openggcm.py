from __future__ import annotations

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


def test_open_dataset():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/coupling0001.iof.000030")
    assert set(ds.coords.keys()) == iof_coords
    assert set(ds.keys()) == iof_vars
