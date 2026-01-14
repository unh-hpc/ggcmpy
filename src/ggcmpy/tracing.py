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

    def B(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [
                self._ds[fld].interp(x=point[0], y=point[1], z=point[2]).to_numpy()
                for fld in ["bx", "by", "bz"]
            ]
        )

    def E(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [
                self._ds[fld].interp(x=point[0], y=point[1], z=point[2]).to_numpy()
                for fld in ["ex", "ey", "ez"]
            ]
        )


class FieldInterpolator_f2py:
    def __init__(self, ds: xr.Dataset) -> None:
        _jrrle.particle_tracing_f2py.load(
            ds.bx, ds.by, ds.bz, ds.ex, ds.ey, ds.ez, ds.x, ds.y, ds.z
        )

    def B(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [_jrrle.particle_tracing_f2py.interpolate(*point, d) for d in range(3)]
        )

    def E(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [_jrrle.particle_tracing_f2py.interpolate(*point, d + 3) for d in range(3)]
        )


class UniformField:
    def __init__(
        self,
        B_0: np.ndarray | None = None,
        E_0: np.ndarray | None = None,
    ) -> None:
        self.B_0 = B_0 if B_0 is not None else np.array([0.0, 0.0, 0.0])
        self.E_0 = E_0 if E_0 is not None else np.array([0.0, 0.0, 0.0])

    def B(self, x: np.ndarray) -> np.ndarray:  # noqa: ARG002
        return self.B_0

    def E(self, x: np.ndarray) -> np.ndarray:  # noqa: ARG002
        return self.E_0


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
            return self._interpolator.B(x)
        return self._interpolator[0](x)  # type: ignore[unreachable]

    def _get_E(self, x: np.ndarray) -> np.ndarray:
        if isinstance(self._interpolator, FieldInterpolator_python):
            return self._interpolator.E(x)
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
