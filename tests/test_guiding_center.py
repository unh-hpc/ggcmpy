from __future__ import annotations

import numpy as np
import pandas as pd

from ggcmpy import constants
from ggcmpy.tracing import emfields, guiding_center, integrator


def test_guiding_center_uniform_gyration():
    """guiding center stays fixed while a particle gyrates in a uniform B field"""
    q = constants.e  # [C]
    m = constants.m_e  # [kg]
    B_0 = 1e-8  # [T]
    v_0 = 0.5 * constants.c
    steps = 100

    fields = emfields.uniform_cxx(B_0=np.array([0.0, 0.0, B_0]))
    u0 = np.array([0.0, v_0, 0.0]) / constants.c
    gamma = np.sqrt(1.0 + np.linalg.norm(u0) ** 2)
    om_ce = q * B_0 / (gamma * m)
    r_ce = np.linalg.norm(u0) * constants.c / abs(om_ce)

    prts_df = pd.DataFrame(
        [
            {
                "id": 0.0,
                "time": 0.0,
                "x": 0.0,
                "y": 0.0,
                "z": 0.0,
                "ux": u0[0],
                "uy": u0[1],
                "uz": u0[2],
            }
        ]
    )
    t_final = 2 * np.pi / abs(om_ce)
    boris = integrator.boris_cxx(fields, q, m)
    _, df = boris.integrate(
        prts_df, t_final=t_final, dt_max_gyro=1.0 / steps, snapshot_interval_steps=1
    )

    df_gc = guiding_center.to_guiding_center(df, q, m, fields)

    assert len(df_gc) == len(df)
    # the guiding center should not move while the particle gyrates
    assert np.std(df_gc.x_gc) < 1e-2 * r_ce
    assert np.std(df_gc.y_gc) < 1e-2 * r_ce
    assert np.allclose(df_gc.z_gc, 0.0)
    # parallel velocity is zero, and the magnetic moment is conserved
    assert np.allclose(df_gc.v_par, 0.0, atol=1e-2 * np.linalg.norm(u0))
    assert np.std(df_gc.mu) < 1e-2 * np.mean(df_gc.mu)
