from __future__ import annotations

import xarray as xr

import ggcmpy
import ggcmpy.tracing
from ggcmpy import _jrrle  # type: ignore[attr-defined]


def load_sample_data() -> xr.Dataset:
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.3df.001200")
    # FIXME, don't have sample data with electric fields
    _jrrle.particle_tracing_f2py.load(
        ds.bx, ds.by, ds.bz, ds.xjx, ds.xjy, ds.xjz, ds.x, ds.y, ds.z
    )
    return ds


def test_load_at():
    ds = load_sample_data()
    idx = 5, 6, 7
    val = _jrrle.particle_tracing_f2py.at(*idx, 2)
    assert val == ds.bz[idx]


def test_interpolate():
    ds = load_sample_data()
    idx = 5, 6, 7
    # interpolate at the cell center, should be same as at()
    val = _jrrle.particle_tracing_f2py.interpolate(
        ds.x[idx[0]], ds.y[idx[1]], ds.z[idx[2]], 2
    )
    assert val == ds.bz[idx]

    # interpolate halfway between two grid points in x direction
    val = _jrrle.particle_tracing_f2py.interpolate(
        0.5 * (ds.x[idx[0]] + ds.x[idx[0] + 1]),
        ds.y[idx[1]],
        ds.z[idx[2]],
        0,
    )
    assert val == 0.5 * (
        ds.bx[idx[0], idx[1], idx[2]] + ds.bx[idx[0] + 1, idx[1], idx[2]]
    )


def test_load_fields():
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.3df.001200")
    ggcmpy.tracing.load_fields(ds)
    idx = 5, 6, 7
    assert ggcmpy.tracing.at(*idx, 2) == ds.bz[idx]
    assert (
        ggcmpy.tracing.interpolate(ds.x[idx[0]], ds.y[idx[1]], ds.z[idx[2]], 2)
        == ds.bz[idx]
    )
