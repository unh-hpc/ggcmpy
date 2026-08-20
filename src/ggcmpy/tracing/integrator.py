"""
integrator.py

Test Particle Integrators
"""

from __future__ import annotations

import pandas as pd
import xarray as xr

from scipy import constants  # type: ignore[import-untyped]
from ggcmpy.tracing import boris_push_cxx, boris_push_python, emfields

# pylint: disable=C0103


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
        boris_push_cls=None,
    ):
        self._fields = fields
        self._q = q
        self._m = m
        self._boris_push_cls = boris_push_cls

    def integrate(
        self, prts_df: pd.DataFrame, t_max, dt_max=1.0, gyro_max=0.1
    ) -> pd.DataFrame:
        boris = self._boris_push_cls(self._fields, self._q, self._m)

        snapshots = [prts_df]

        while prts_df.iloc[0].time < t_max:
            # hack to make the boris push do just one time step
            prts_df = boris.push(prts_df, prts_df.iloc[0].time + 1e-7, dt_max, gyro_max)
            snapshots.append(prts_df)

        return pd.concat(snapshots, ignore_index=True)


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

        super().__init__(fields, q, m, boris_push_cls=boris_push_python)


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

        super().__init__(fields, q, m, boris_push_cls=boris_push_cxx)
