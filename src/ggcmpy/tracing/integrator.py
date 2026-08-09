"""
integrator.py

Test Particle Integrators
"""

from __future__ import annotations

import pandas as pd
import xarray as xr

import ggcmpy
from ggcmpy import (
    constants,
)
from ggcmpy.tracing import boris_push_cxx, emfields


class boris_base:
    """
    Base class for Boris particle pusher.

    Methods:
        push(prts, t_max, dt_max, gyro_max):
            Pushes the particles in `prts` from their current state to a maximum time of `t_max`,
            using a maximum time step of `dt_max` and a maximum gyro period of `gyro_max`.
    """

    def __init__(
        self,
        fields: emfields.emfields,
        q=constants.e,
        m=constants.m_e,
    ):
        self._fields = fields
        self._q = q
        self._m = m


class boris_python(boris_base):
    """
    Boris particle pusher implemented in pure Python

    Methods:
        push(prts, t_max, dt_max, gyro_max):
            Pushes the particles in `prts` from their current state to a maximum time of `t_max`,
            using a maximum time step of `dt_max` and a maximum gyro period of `gyro_max`.
    """

    def __init__(
        self,
        fields: xr.Dataset | emfields.emfields,
        q=constants.e,
        m=constants.m_e,
    ):
        if isinstance(fields, xr.Dataset):
            fields = emfields.yee_cic_python(fields)
        super().__init__(fields, q, m)

    def integrate(
        self, prts_df: pd.DataFrame, t_max, dt_max=1.0, gyro_max=0.1
    ) -> pd.DataFrame:
        boris = ggcmpy.tracing.BorisIntegrator_python(self._fields, self._q, self._m)

        x0 = prts_df.iloc[0][["x", "y", "z"]].to_numpy()
        u0 = prts_df.iloc[0][["ux", "uy", "uz"]].to_numpy()
        return boris.integrate(x0, u0, t_max, dt_max, gyro_max)


class boris_cxx(boris_base):
    """
    Boris particle pusher implemented in C++.

    Methods:
        push(prts, t_max, dt_max, gyro_max):
            Pushes the particles in `prts` from their current state to a maximum time of `t_max`,
            using a maximum time step of `dt_max` and a maximum gyro period of `gyro_max`.
    """

    def __init__(
        self,
        fields: xr.Dataset | emfields.emfields,
        q=constants.e,
        m=constants.m_e,
    ):
        if isinstance(fields, xr.Dataset):
            fields = emfields.yee_cic_cxx(fields)
        super().__init__(fields, q, m)

    def integrate(self, prts_df, t_max, dt_max=1.0, gyro_max=0.1) -> pd.DataFrame:
        boris = boris_push_cxx(self._fields, self._q, self._m)

        snapshots = [prts_df]

        while prts_df.iloc[0].time < t_max:
            # hack to make the boris push do just one time step
            prts_df = boris.push(prts_df, prts_df.iloc[0].time + 1e-7, dt_max, gyro_max)
            snapshots.append(prts_df)

        snapshots = pd.concat(snapshots, ignore_index=True)
        assert isinstance(snapshots, pd.DataFrame)
        return snapshots
