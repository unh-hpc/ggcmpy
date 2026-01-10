from __future__ import annotations

import xarray as xr

import ggcmpy
from ggcmpy import _jrrle


def test_load_fields():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.3df.001200")
    ds = ds.isel(time=0)
    # FIXME, don't have sample data with electric fields
    _jrrle.particle_tracing_f2py.load(ds.bx, ds.by, ds.bz, ds.xjx, ds.xjy, ds.xjz)
    idx = 5, 6, 7, 2
    val = _jrrle.particle_tracing_f2py.at(*idx)
    assert val == ds.bz[idx[0:3]]
