"""
tracing.py

Particle tracing
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr
from scipy import constants  # type: ignore[import-untyped]

from ggcmpy import _jrrle  # type: ignore[attr-defined]


def make_vector_field(
    grid: Sequence[tuple[str, tuple[str, ...]]],
    coords: dict[str, np.ndarray],
    vector_field: Callable[[np.ndarray], np.ndarray],
) -> dict[str, Any]:
    flds = {}
    for d, (fld_name, dims) in enumerate(grid):
        crds = [coords[dim] for dim in dims]
        fld = np.empty(tuple(len(c) for c in crds))
        for i in range(fld.shape[0]):
            for j in range(fld.shape[1]):
                for k in range(fld.shape[2]):
                    val = vector_field(np.array([crds[0][i], crds[1][j], crds[2][k]]))
                    fld[i, j, k] = val[d]
        flds[fld_name] = (dims, fld)

    return flds


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


class FieldInterpolatorYee_python:
    def __init__(self, ds: xr.Dataset) -> None:
        assert {"bx1", "by1", "bz1", "ex1", "ey1", "ez1"} <= ds.data_vars.keys()
        self._ds = ds

    def B(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [self._interpolate(self._ds[fld], point) for fld in ["bx1", "by1", "bz1"]]
        )

    def E(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [self._interpolate(self._ds[fld], point) for fld in ["ex1", "ey1", "ez1"]]
        )

    def _interpolate(self, da: xr.DataArray, point: np.ndarray) -> float:
        val = da.interp(dict(zip(da.dims, point, strict=True))).to_numpy()
        return float(val)


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


class FieldInterpolatorYee_f2py:
    def __init__(self, ds: xr.Dataset) -> None:
        _jrrle.particle_tracing_f2py.load_yee(
            ds.bx1,
            ds.by1,
            ds.bz1,
            ds.ex1,
            ds.ey1,
            ds.ez1,
            ds.x,
            ds.y,
            ds.z,
            ds.x_nc,
            ds.y_nc,
            ds.z_nc,
        )

    def B(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [_jrrle.particle_tracing_f2py.interpolate_yee(*point, d) for d in range(3)]
        )

    def E(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [
                _jrrle.particle_tracing_f2py.interpolate_yee(*point, d + 3)
                for d in range(3)
            ]
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


class DipoleField:
    def __init__(self, m):
        self.m = m

    def B(self, r):
        rhat = r / np.linalg.norm(r)
        return (
            constants.mu_0
            / (4 * np.pi)
            * (3 * np.dot(self.m, rhat) * rhat - self.m)
            / np.linalg.norm(r) ** 3
        )

    def E(self, r):  # noqa: ARG002
        return np.array([0.0, 0.0, 0.0])


class BorisIntegrator_python:
    def __init__(self, ds, q=constants.e, m=constants.m_e) -> None:
        self.q = q
        self.m = m
        if isinstance(ds, xr.Dataset):
            self._interpolator = FieldInterpolator_python(ds)
        else:
            self._interpolator = ds  # assume it's already an interpolator

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
            B = self._interpolator.B(x)
            E = self._interpolator.E(x)
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


class BorisIntegrator_f2py:
    def __init__(self, df, q=constants.e, m=constants.m_e) -> None:
        _jrrle.particle_tracing_f2py.boris_init(q, m)
        if isinstance(df, xr.Dataset):
            self._interpolator = FieldInterpolator_f2py(df)
        else:
            assert isinstance(df, FieldInterpolator_f2py)
            self._interpolator = df

    def integrate(self, x0, v0, t_max, dt) -> pd.DataFrame:
        n_steps = int(t_max / dt) + 2  # add some extra space for round-off issues
        data = np.zeros((7, n_steps), dtype=np.float32, order="F")
        n_out = _jrrle.particle_tracing_f2py.boris_integrate(x0, v0, t_max, dt, data)
        return pd.DataFrame(
            data.T[:n_out], columns=["time", "x", "y", "z", "vx", "vy", "vz"]
        )


BorisIntegrator = BorisIntegrator_python
