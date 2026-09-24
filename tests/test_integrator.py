from __future__ import annotations

import matplotlib.pyplot as plt  # type: ignore[import-not-found]
import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest

import ggcmpy.tracing
from ggcmpy import constants
from ggcmpy.tracing import emfields, integrator

R_E = constants.radius_earth  # [m]


def make_particle(
    id: float, x0: npt.ArrayLike, v0: npt.ArrayLike
) -> tuple[float, np.ndarray, np.ndarray]:
    x0, v0 = np.asarray(x0), np.asarray(v0)
    gamma = 1.0 / np.sqrt(1 - (np.linalg.norm(v0) / constants.c) ** 2)
    u0 = gamma * v0 / constants.c
    return id, 0.0, *x0, *u0


def to_prts_df(particles: list[tuple[float, np.ndarray, np.ndarray]]) -> pd.DataFrame:
    return pd.DataFrame(
        np.array(particles), columns=["id", "time", "x", "y", "z", "ux", "uy", "uz"]
    )


def gyro_frequency(B: float, q: float, m: float, u: np.ndarray) -> float:
    gamma = np.sqrt(1.0 + np.linalg.norm(u) ** 2)
    return np.abs(q) * B / (gamma * m)  # type: ignore[no-any-return]


def gyro_radius(B: float, q: float, m: float, u: np.ndarray) -> float:
    return m * np.linalg.norm(u) * constants.c / (np.abs(q) * B)  # type: ignore[no-any-return]


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
    prts_df = to_prts_df([make_particle(0, x0, v0)])

    u0 = prts_df.loc[0, ["ux", "uy", "uz"]].to_numpy()
    om_ce = gyro_frequency(B_0, q, m, u0)
    r_ce = gyro_radius(B_0, q, m, u0)

    t_final = 2 * np.pi / om_ce  # one gyroperiod # [s]
    steps = 100

    boris = integrator(fields, q, m)
    df = boris.integrate(
        prts_df, t_final=t_final, dt_max_gyro=1.0 / steps, snapshot_interval_steps=1
    )

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

    fields = emfields.dipole_cxx(m=constants.dipole_moment_earth)  # [A m^2]

    q = -constants.e
    m = constants.m_e
    x0 = np.array([5.0 * R_E, 0.0, 0.0])  # [m]
    B_0 = np.linalg.norm(fields.B(x0))
    E_kin = 1000.0 * 1e3 * constants.e  # 1000 keV in J
    gamma = 1.0 + E_kin / (m * constants.c**2)
    v_e = constants.c * np.sqrt(1.0 - 1.0 / gamma**2)

    v0 = np.array([0.0, v_e / np.sqrt(2.0), v_e / np.sqrt(2.0)])  # [m/s]
    prts = to_prts_df([make_particle(0, x0, v0)])
    u0 = prts.loc[0, ["ux", "uy", "uz"]].to_numpy()
    om_ce = gyro_frequency(B_0, q, m, u0)
    r_ce = gyro_radius(B_0, q, m, u0)

    print(f"B={B_0} [T] om_ce={om_ce:.2f} [1/s] r_ce={r_ce:.2f} [m]")

    t_ce = 2.0 * np.pi / om_ce  # [s]
    t_final = 100.0 * t_ce  # [s]

    boris = ggcmpy.tracing.integrator.boris_cxx(fields, q, m)
    df = boris.integrate(prts, t_final=t_final, snapshot_interval_steps=1)

    B_final = np.linalg.norm(fields.B(df.loc[df.index[-1], ["x", "y", "z"]].to_numpy()))
    om_ce_final = gyro_frequency(B_final, q, m, u0)
    t_ce_final = 2.0 * np.pi / om_ce_final

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    df[df.time < 5.0 * t_ce].plot(
        x="x", y="z", style=".-", ax=axs[0], title="First 5 t_ce"
    )
    df[df.time >= t_final - 5.0 * t_ce_final].plot(
        x="x", y="z", style=".-", ax=axs[1], title="Last 5 t_ce "
    )
    df.plot(x="x", y="z", style="-", ax=axs[2], title="All steps")
    fig.tight_layout()

    return fig


@pytest.mark.parametrize(
    "integrator",
    [
        integrator.boris_python,
        integrator.boris_cxx,
    ],
)
def test_boris_integrator_snapshot(integrator):
    """
    Integrate particle gyrating / bouncing in a dipole magnetic field.

    Taking snapshots at different intervals should not change the result of the integration.
    """

    fields = emfields.dipole_cxx(m=constants.dipole_moment_earth)  # [A m^2]

    q = -constants.e
    m = constants.m_e
    x0 = np.array([5.0 * R_E, 0.0, 0.0])  # [m]
    B_0 = np.linalg.norm(fields.B(x0))
    E_kin = 1000.0 * 1e3 * constants.e  # 1000 keV in J
    gamma = 1.0 + E_kin / (m * constants.c**2)
    v_e = constants.c * np.sqrt(1.0 - 1.0 / gamma**2)

    v0 = np.array([0.0, v_e / np.sqrt(2.0), v_e / np.sqrt(2.0)])  # [m/s]
    prts = to_prts_df([make_particle(0, x0, v0)])

    u0 = prts.loc[0, ["ux", "uy", "uz"]].to_numpy()
    om_ce = gyro_frequency(B_0, q, m, u0)
    r_ce = gyro_radius(B_0, q, m, u0)

    print(f"B={B_0} [T] om_ce={om_ce:.2f} [1/s] r_ce={r_ce:.2f} [m]")

    t_ce = 2.0 * np.pi / om_ce  # [s]
    t_final = 1.0 * t_ce  # [s]

    boris = integrator(fields, q, m)
    df = boris.integrate(prts, t_final=t_final, snapshot_interval_steps=1)
    df2 = boris.integrate(prts, t_final=t_final, snapshot_interval_steps=10)
    df = df.iloc[::10]
    df2 = df2.iloc[
        : len(df)
    ]  # this one might have an additional final step, which we drop here

    assert np.allclose(df.to_numpy(), df2.to_numpy())


@pytest.mark.mpl_image_compare(filename="test_boris_integrator_multiple.png")
@pytest.mark.parametrize(
    "integrator",
    [
        integrator.boris_python,
        integrator.boris_cxx,
    ],
)
def test_boris_integrator_multiple(integrator):
    """multiple particle gyrating in a uniform magnetic field"""
    q = constants.e  # [C]
    m = constants.m_e  # [kg]
    B_0 = 1e-8  # [T]
    fields = emfields.uniform_cxx(B_0=np.array([0.0, 0.0, B_0]))

    prts_df = to_prts_df(
        [
            make_particle(id, [0.0, 0.0, 0.0], [0.0, v_0, 0.0])
            for id, v_0 in enumerate([0.3 * constants.c, 0.5 * constants.c])
        ]
    )

    om_ce = gyro_frequency(B_0, q, m, prts_df.loc[0, ["ux", "uy", "uz"]].to_numpy())

    t_final = 2 * np.pi / om_ce  # one gyroperiod # [s]
    dt_max_gyro = 1.0 / 100

    boris = integrator(fields, q, m)
    df = boris.integrate(
        prts_df, t_final=t_final, dt_max_gyro=dt_max_gyro, snapshot_interval_steps=1
    )

    fig, ax = plt.subplots()
    df.plot.scatter(
        x="x",
        y="y",
        c="id",
        style=".",
        ax=ax,
        cmap="viridis",
        title="Particle trajectories",
    )
    ax.set_aspect("equal")
    fig.tight_layout()

    return fig
