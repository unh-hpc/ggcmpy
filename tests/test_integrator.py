from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ggcmpy import constants
from ggcmpy.tracing import emfields, integrator


@pytest.mark.parametrize(
    "integrator",
    [
        integrator.boris_python,
        integrator.boris_cxx,
    ],
)
def test_boris_integrator_uniform(integrator):
    """particle gyrating in a uniform magnetic field"""
    q = constants.e  # [C]
    m = constants.m_e  # [kg]
    B_0 = 1e-8  # [T]
    v_0 = 0.5 * constants.c
    fields = emfields.uniform_cxx(B_0=np.array([0.0, 0.0, B_0]))
    x0 = np.array([0.0, 0.0, 0.0])  # [m]
    v0 = np.array([0.0, v_0, 0.0])  # [m/s]
    gamma = 1.0 / np.sqrt(1 - (np.linalg.norm(v0) / constants.c) ** 2)
    u0 = gamma * v0 / constants.c
    om_ce = q * B_0 / (gamma * m)  # [rad/s]
    r_ce = m * np.linalg.norm(u0) * constants.c / (np.abs(q) * B_0)  # [m]
    t_max = 2 * np.pi / om_ce  # one gyroperiod # [s]
    steps = 100
    dt = t_max / steps  # [s]
    prts_df = pd.DataFrame(
        np.array([[0.0, *x0, *u0]]), columns=["time", "x", "y", "z", "ux", "uy", "uz"]
    )

    boris = integrator(fields, q, m)
    df = boris.integrate(prts_df, t_max, dt)

    assert len(df) >= steps
    assert len(df) <= steps + 2

    assert np.allclose(df.ux, np.sin(om_ce * df.time) * u0[1], atol=1e-2 * u0[1])
    assert np.allclose(df.uy, np.cos(om_ce * df.time) * u0[1], atol=1e-2 * u0[1])
    assert np.allclose(df.uz, 0.0)

    assert np.allclose(df.x, r_ce * (1 - np.cos(om_ce * df.time)), atol=1e-2 * r_ce)
    assert np.allclose(df.y, r_ce * (np.sin(om_ce * df.time)), atol=1e-2 * r_ce)
    assert np.allclose(df.z, 0.0)
