from __future__ import annotations

import matplotlib.pyplot as plt  # type: ignore[import-not-found]
import numpy as np
import pandas as pd
import pytest

import ggcmpy.tracing
from ggcmpy import constants
from ggcmpy.tracing import emfields, integrator

R_E = constants.radius_earth  # [m]


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


@pytest.mark.mpl_image_compare
def test_boris_integrator_dipole():
    """particle gyrating / bouncing in a dipole magnetic field"""

    def gyro_frequency(B: float, q: float, m: float, gamma: float) -> float:
        return np.abs(q) * B / (gamma * m)  # type: ignore[no-any-return]

    fields = emfields.dipole_cxx(m=constants.dipole_moment_earth)  # [A m^2]

    q = -constants.e
    m = constants.m_e
    x0 = np.array([5.0 * R_E, 0.0, 0.0])  # [m]
    B_0 = np.linalg.norm(fields.B(x0))
    E_kin = 1000.0 * 1e3 * constants.e  # 1000 keV in J
    gamma = 1.0 + E_kin / (m * constants.c**2)
    v_e = constants.c * np.sqrt(1.0 - 1.0 / gamma**2)

    v0 = np.array([0.0, v_e / np.sqrt(2.0), v_e / np.sqrt(2.0)])  # [m/s]
    u0 = gamma * v0 / constants.c

    om_ce = gyro_frequency(B_0, q, m, gamma)
    r_ce = m * np.linalg.norm(u0) * constants.c / (np.abs(q) * B_0)  # [m]

    print(f"B={B_0} [T] om_ce={om_ce:.2f} [1/s] r_ce={r_ce:.2f} [m]")

    t_ce = 2.0 * np.pi / om_ce  # [s]
    t_max = 100.0 * t_ce  # [s]

    prts = pd.DataFrame(
        np.array([[0.0, *x0, *u0]]), columns=["time", "x", "y", "z", "ux", "uy", "uz"]
    )

    boris = ggcmpy.tracing.integrator.boris_cxx(fields, q, m)
    df = boris.integrate(prts, t_max)

    B_final = np.linalg.norm(fields.B(df.loc[df.index[-1], ["x", "y", "z"]].to_numpy()))
    om_ce_final = gyro_frequency(B_final, q, m, gamma)
    t_ce_final = 2.0 * np.pi / om_ce_final

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    df[df.time < 5.0 * t_ce].plot(
        x="x", y="z", style=".-", ax=axs[0], title="First 5 t_ce"
    )
    df[df.time >= t_max - 5.0 * t_ce_final].plot(
        x="x", y="z", style=".-", ax=axs[1], title="Last 5 t_ce "
    )
    df.plot(x="x", y="z", style="-", ax=axs[2], title="All steps")
    fig.tight_layout()

    return fig
