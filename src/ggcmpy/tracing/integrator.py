"""
integrator.py

Test Particle Integrators
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

import ggcmpy
from ggcmpy import (
    constants,
)
from ggcmpy.tracing import boris_cxx as boris_cxx_impl
from ggcmpy.tracing import emfields


class boris_python:
    """
    Boris particle pusher implemented in pure Python

    Methods:
        push(prts, t_max, dt_max, gyro_max):
            Pushes the particles in `prts` from their current state to a maximum time of `t_max`,
            using a maximum time step of `dt_max` and a maximum gyro period of `gyro_max`.
    """

    def __init__(
        self,
        df: xr.Dataset | emfields.emfields,
        q=constants.e,
        m=constants.m_e,
    ):
        if isinstance(df, xr.Dataset):
            self._emfields: emfields.emfields = emfields.yee_cic_python(df)
        else:
            self._emfields = df

        self._q = q
        self._m = m

    def integrate(self, x0, u0, t_max, dt_max=1.0, gyro_max=0.1) -> pd.DataFrame:
        boris = ggcmpy.tracing.BorisIntegrator_python(self._emfields, self._q, self._m)

        return boris.integrate(x0, u0, t_max, dt_max, gyro_max)


class boris_cxx:
    """
    Boris particle pusher implemented in C++.

    Methods:
        push(prts, t_max, dt_max, gyro_max):
            Pushes the particles in `prts` from their current state to a maximum time of `t_max`,
            using a maximum time step of `dt_max` and a maximum gyro period of `gyro_max`.
    """

    def __init__(
        self,
        df: xr.Dataset | emfields.emfields,
        q=constants.e,
        m=constants.m_e,
    ):
        if isinstance(df, xr.Dataset):
            self._emfields: emfields.emfields = emfields.yee_cic_cxx(df)
        else:
            assert isinstance(df, (emfields.uniform_cxx, emfields.yee_cic_cxx))
            self._emfields = df

        self._q = q
        self._m = m

    def integrate(self, x0, u0, t_max, dt_max=1.0, gyro_max=0.1) -> pd.DataFrame:
        boris = boris_cxx_impl(self._emfields, self._q, self._m)

        prts_df = pd.DataFrame(
            np.array([[0.0, *x0, *u0]]),
            columns=["time", "x", "y", "z", "ux", "uy", "uz"],
        )
        snapshots = [prts_df]

        while prts_df.iloc[0].time < t_max:
            # hack to make the boris push do just one time step
            prts_df = boris.push(prts_df, prts_df.iloc[0].time + 1e-7, dt_max, gyro_max)
            snapshots.append(prts_df)

        return pd.concat(snapshots, ignore_index=True)
