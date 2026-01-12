from __future__ import annotations

import numpy as np
import xarray as xr
from scipy import constants  # type: ignore[import-untyped]

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
    """test the high-level f2py interfaces"""
    ds = xr.open_dataset(f"{ggcmpy.sample_dir}/sample_jrrle.3df.001200")
    ggcmpy.tracing.load_fields(ds)
    idx = 5, 6, 7
    assert ggcmpy.tracing.at(*idx, 2) == ds.bz[idx]
    assert (
        ggcmpy.tracing.interpolate(ds.x[idx[0]], ds.y[idx[1]], ds.z[idx[2]], 2)
        == ds.bz[idx]
    )


def test_BorisIntegrator_python():
    """particle gyrating in a uniform magnetic field"""
    q = constants.e  # [C]
    m = constants.m_e  # [kg]
    B_0 = 1e-8  # [T]
    E_0 = 0.0  # [V/m]
    x0 = np.array([0.0, 0.0, 0.0])  # [m]
    v0 = np.array([0.0, 100.0, 0.0])  # [m/s]
    om_ce = q * B_0 / m  # [rad/s]
    t_max = 2 * np.pi / om_ce  # one gyroperiod # [s]
    steps = 100
    dt = t_max / steps  # [s]

    def get_B(x):  # noqa: ARG001
        return np.array([0.0, 0.0, B_0])  # [T]

    def get_E(x):  # noqa: ARG001
        return np.array([0.0, 0.0, E_0])  # [V/m]

    boris = ggcmpy.tracing.BorisIntegrator_python(get_B, get_E, q, m)
    df = boris.integrate(x0, v0, t_max, dt)

    assert len(df) == steps + 1
    assert np.allclose(df.iloc[0][["x", "y", "z"]], x0)
    # after half a gyroperiod, should have moved from initial position
    assert not np.allclose(df.iloc[steps // 2][["x", "y", "z"]], x0, atol=1e-3)
    # after one gyroperiod, should return to near the initial position
    assert np.allclose(df.iloc[-1][["x", "y", "z"]], x0, atol=1e-3)


def test_BorisIntegrator_f2py():
    """particle gyrating in a uniform magnetic field using the f2py interface"""
    q = constants.e  # [C]
    m = constants.m_e  # [kg]
    B_0 = 1e-8  # [T]
    # E_0 = 0.0  # [V/m]
    x0 = np.array([0.0, 0.0, 0.0])  # [m]
    v0 = np.array([0.0, 100.0, 0.0])  # [m/s]
    om_ce = q * B_0 / m  # [rad/s]
    t_max = 2 * np.pi / om_ce  # one gyroperiod # [s]
    steps = 100
    dt = t_max / steps  # [s]

    boris = ggcmpy.tracing.BorisIntegrator_f2py(None, None, q, m)
    boris.integrate(x0, v0, t_max, dt)

    # assert len(df) == steps + 1
    # assert np.allclose(df.iloc[0][["x", "y", "z"]], x0)
    # # after half a gyroperiod, should have moved from initial position
    # assert not np.allclose(df.iloc[steps // 2][["x", "y", "z"]], x0, atol=1e-3)
    # # after one gyroperiod, should return to near the initial position
    # assert np.allclose(df.iloc[-1][["x", "y", "z"]], x0, atol=1e-3)
