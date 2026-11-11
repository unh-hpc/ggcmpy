from __future__ import annotations

import numpy as np
import xarray as xr

from ggcmpy.tracing import emfields

NX, NY, NZ = 6, 5, 4


def make_dataset() -> xr.Dataset:
    x = np.linspace(0.0, 1.0, NX)
    y = np.linspace(0.0, 1.0, NY)
    z = np.linspace(0.0, 1.0, NZ)
    rng = np.random.default_rng(0)
    data = {
        name: (("z", "y", "x"), rng.random((NZ, NY, NX)))
        for name in ["bx1", "by1", "bz1", "eflx", "efly", "eflz"]
    }
    return xr.Dataset(data, coords={"x": x, "y": y, "z": z})


def test_prepare_yee_dataset_coords_and_shapes():
    ds = emfields.prepare_yee_dataset(make_dataset())

    assert ds.x_nc.size == NX - 1
    assert ds.y_nc.size == NY - 1
    assert ds.z_nc.size == NZ - 1

    assert ds.bx1.dims == ("x_nc", "y", "z")
    assert ds.by1.dims == ("x", "y_nc", "z")
    assert ds.bz1.dims == ("x", "y", "z_nc")
    assert ds.eflx.dims == ("x", "y_nc", "z_nc")
    assert ds.efly.dims == ("x_nc", "y", "z_nc")
    assert ds.eflz.dims == ("x_nc", "y_nc", "z")

    assert ds.bx1.shape == (NX - 1, NY, NZ)
    assert ds.eflz.shape == (NX - 1, NY - 1, NZ)


def test_prepare_yee_dataset_midpoints():
    ds = make_dataset()
    out = emfields.prepare_yee_dataset(ds)

    for d in ["x", "y", "z"]:
        c = ds[d].to_numpy()
        np.testing.assert_allclose(out[f"{d}_nc"].to_numpy(), 0.5 * (c[1:] + c[:-1]))


def test_prepare_yee_dataset_works_with_interpolator():
    out = emfields.prepare_yee_dataset(make_dataset())

    fields = emfields.yee_cic_python(out)
    point = np.array([0.5, 0.5, 0.5])
    assert np.all(np.isfinite(fields.B(point)))
    assert np.all(np.isfinite(fields.E(point)))
