"""
tracing.py

Particle tracing
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr
from scipy import constants  # type: ignore[import-untyped]

from ggcmpy import _jrrle  # type: ignore[attr-defined]


class FieldInterpolator_python:
    def __init__(self, ds: xr.Dataset) -> None:
        assert {"bx", "by", "bz", "ex", "ey", "ez"} <= ds.data_vars.keys()
        self._ds = ds

    def __call__(self, point: np.ndarray, m: int) -> np.ndarray:
        return (
            self._get_component(m).interp(x=point[0], y=point[1], z=point[2]).to_numpy()
        )

    def at(self, idx: np.ndarray, m: int) -> np.ndarray:
        return self._get_component(m)[idx].to_numpy()

    def _get_component(self, m: int) -> xr.DataArray:
        components = ["bx", "by", "bz", "ex", "ey", "ez"]
        return self._ds[components[m]]


class FieldInterpolator_f2py:
    def __init__(self, ds: xr.Dataset) -> None:
        _jrrle.particle_tracing_f2py.load(
            ds.bx, ds.by, ds.bz, ds.ex, ds.ey, ds.ez, ds.x, ds.y, ds.z
        )

    def __call__(self, point: np.ndarray, m: int) -> float:
        """Interpolate the field value at the given spatial coordinates.

        Parameters
        ----------
        x : float
            Spatial coordinate x.
        y : float
            Spatial coordinate y.
        z : float
            Spatial coordinate z.
        m : int
            Field component index (0: bx, 1: by, 2: bz, 3: ex, 4: ey, 5: ez).

        Returns
        -------
        float
            Interpolated field value at the specified coordinates and component.
        """
        val = _jrrle.particle_tracing_f2py.interpolate(*point, m)
        assert isinstance(val, float)
        return val

    def at(self, idx: np.ndarray, m: int) -> float:
        """Get the field value at the given grid index.

        Parameters
        ----------
        i : int
            Grid index i.
        j : int
            Grid index j.
        k : int
            Grid index k.
        m : int
            Field component index (0: bx, 1: by, 2: bz, 3: ex, 4: ey, 5: ez).

        Returns
        -------
        float
            Field value at the specified index and component.
        """
        val = _jrrle.particle_tracing_f2py.at(*idx, m)
        assert isinstance(val, float)
        return val


class BorisIntegrator_python:
    def __init__(self, ds, q=constants.e, m=constants.m_e) -> None:
        self.q = q
        self.m = m
        if isinstance(ds, xr.Dataset):
            self._interpolator = FieldInterpolator_python(ds)
        else:
            self._interpolator = ds  # expect a callable

    def integrate(self, x0, v0, t_max, dt) -> pd.DataFrame:
        t = 0.0
        x = x0.copy()
        v = v0.copy()
        qprime = 0.5 * dt * self.q / self.m
        times, positions, velocities = [], [], []
        while t < t_max:
            times.append(t)
            positions.append(x.copy())
            velocities.append(v.copy())
            B = self._get_B(x)
            E = self._get_E(x)
            x += 0.5 * dt * v
            v += qprime * E
            h = qprime * B
            s = 2 * h / (1 + np.abs(h) ** 2)
            v += np.cross(v + np.cross(v, h), s)
            v += qprime * E
            x += 0.5 * dt * v
            t += dt

        return pd.DataFrame(
            np.column_stack((times, positions, velocities)),
            columns=["time", "x", "y", "z", "vx", "vy", "vz"],
        )

    def _get_B(self, x: np.ndarray) -> np.ndarray:
        if isinstance(self._interpolator, FieldInterpolator_python):
            return np.array([self._interpolator(x, d) for d in range(3)])
        return self._interpolator[0](x)  # type: ignore[unreachable]

    def _get_E(self, x: np.ndarray) -> np.ndarray:
        if isinstance(self._interpolator, FieldInterpolator_python):
            return np.array([self._interpolator(x, d + 3) for d in range(3)])
        return self._interpolator[1](x)  # type: ignore[unreachable]


class BorisIntegrator_f2py:
    def __init__(self, df, q=constants.e, m=constants.m_e) -> None:
        _jrrle.particle_tracing_f2py.boris_init(q, m)
        self._interpolator = FieldInterpolator_f2py(df)

    def integrate(self, x0, v0, t_max, dt) -> pd.DataFrame:
        n_steps = int(t_max / dt) + 2  # add some extra space for round-off issues
        data = np.zeros((7, n_steps), dtype=np.float32, order="F")
        n_out = _jrrle.particle_tracing_f2py.boris_integrate(x0, v0, t_max, dt, data)
        return pd.DataFrame(
            data.T[:n_out], columns=["time", "x", "y", "z", "vx", "vy", "vz"]
        )


BorisIntegrator = BorisIntegrator_python
