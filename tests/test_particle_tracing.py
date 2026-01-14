from __future__ import annotations

import numpy as np
import pytest
import xarray as xr
from scipy import constants  # type: ignore[import-untyped]

import ggcmpy
import ggcmpy.tracing
from ggcmpy import _jrrle  # type: ignore[attr-defined]


def load_sample_data() -> xr.Dataset:
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.3df.001200")
    ds["ex"] = xr.zeros_like(ds.bx)
    ds["ey"] = xr.zeros_like(ds.by)
    ds["ez"] = xr.zeros_like(ds.bz)
    # FIXME, don't have sample data with electric fields
    _jrrle.particle_tracing_f2py.load(
        ds.bx, ds.by, ds.bz, ds.ex, ds.ey, ds.ez, ds.x, ds.y, ds.z
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


@pytest.mark.parametrize(
    "FieldInterpolator",
    [
        ggcmpy.tracing.FieldInterpolator_python,
        ggcmpy.tracing.FieldInterpolator_f2py,
    ],
)
def test_FieldInterpolator(FieldInterpolator):
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.3df.001200")
    ds["ex"] = xr.zeros_like(ds.bx)
    ds["ey"] = xr.zeros_like(ds.by)
    ds["ez"] = xr.zeros_like(ds.bz)
    interpolator = FieldInterpolator(ds)
    idx = 5, 6, 7
    assert interpolator.B((ds.x[idx[0]], ds.y[idx[1]], ds.z[idx[2]]))[2] == ds.bz[idx]


@pytest.mark.parametrize(
    "Integrator",
    [ggcmpy.tracing.BorisIntegrator_python, ggcmpy.tracing.BorisIntegrator_f2py],
)
def test_BorisIntegrator(Integrator):
    """particle gyrating in a uniform magnetic field"""
    q = constants.e  # [C]
    m = constants.m_e  # [kg]
    B_0 = 1e-8  # [T]
    E_0 = 0.0  # [V/m]
    shape = (10, 5, 15)
    crd = [np.linspace(-1.0, 1.0, d) for d in shape]
    df = xr.Dataset(
        data_vars={
            "ex": (("x", "y", "z"), np.zeros(shape)),
            "ey": (("x", "y", "z"), np.zeros(shape)),
            "ez": (("x", "y", "z"), E_0 * np.ones(shape)),
            "bx": (("x", "y", "z"), np.zeros(shape)),
            "by": (("x", "y", "z"), np.zeros(shape)),
            "bz": (("x", "y", "z"), B_0 * np.ones(shape)),
        },
        coords={"x": ("x", crd[0]), "y": ("y", crd[1]), "z": ("z", crd[2])},
    )
    x0 = np.array([0.0, 0.0, 0.0])  # [m]
    v0 = np.array([0.0, 100.0, 0.0])  # [m/s]
    om_ce = q * B_0 / m  # [rad/s]
    r_ce = np.linalg.norm(v0) / om_ce  # [m]
    t_max = 2 * np.pi / om_ce  # one gyroperiod # [s]
    steps = 100
    dt = t_max / steps  # [s]

    boris = Integrator(df, q, m)
    df = boris.integrate(x0, v0, t_max, dt)

    assert len(df) == steps + 1

    assert np.allclose(df.vx, np.sin(om_ce * df.time) * v0[1], atol=1.0)
    assert np.allclose(df.vy, np.cos(om_ce * df.time) * v0[1], atol=1.0)
    assert np.allclose(df.vz, 0.0)

    assert np.allclose(df.x, r_ce * (1 - np.cos(om_ce * df.time)), atol=1e-3)
    assert np.allclose(df.y, r_ce * (np.sin(om_ce * df.time)), atol=1e-3)
    assert np.allclose(df.z, 0.0)
