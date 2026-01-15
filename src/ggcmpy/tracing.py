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

# pylint: disable=C0103,I1101


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
    """
    FieldInterpolator_python provides interpolation of electromagnetic field components
    from an xarray.Dataset at arbitrary 3D points.

    Methods:
        B(point: np.ndarray) -> np.ndarray:
            Interpolates and returns the magnetic field vector [bx, by, bz] at the given 3D point.

        E(point: np.ndarray) -> np.ndarray:
            Interpolates and returns the electric field vector [ex, ey, ez] at the given 3D point.

    Args:
        ds (xr.Dataset): An xarray dataset containing the required field components as data variables.
    """

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
    """
    FieldInterpolatorYee_python provides interpolation of electromagnetic field components
    (B and E fields) from an xarray.Dataset on a Yee grid.

    Methods:
        B(point: np.ndarray) -> np.ndarray:
            Interpolates and returns the magnetic field vector (B) at the specified point.

        E(point: np.ndarray) -> np.ndarray:
            Interpolates and returns the electric field vector (E) at the specified point.

    Args:
        ds (xr.Dataset): An xarray dataset containing the required field components.
    """

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
    """
    FieldInterpolator_f2py provides an interface to interpolate electromagnetic field values
    from a given xarray.Dataset using a Fortran backend via f2py.

    Methods:
        __init__(ds: xr.Dataset)
            Initializes the interpolator by loading field data (bx, by, bz, ex, ey, ez, x, y, z)
            from the provided xarray.Dataset into the Fortran backend.
        B(point: np.ndarray) -> np.ndarray
            Interpolates and returns the magnetic field vector (B) at the specified spatial point.
        E(point: np.ndarray) -> np.ndarray
            Interpolates and returns the electric field vector (E) at the specified spatial point.

    Args:
        ds (xr.Dataset): An xarray dataset containing the required field components.
    """

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
    """
    FieldInterpolatorYee_f2py provides an interface for interpolating electromagnetic field components
    (B and E fields) at arbitrary points using Yee grid data loaded from an xarray.Dataset.

    Methods:
        __init__(ds: xr.Dataset)
            Initializes the interpolator by loading Yee grid field data from the provided xarray.Dataset.
        B(point: np.ndarray) -> np.ndarray
            Interpolates and returns the magnetic field vector (B) at the specified spatial point.
        E(point: np.ndarray) -> np.ndarray
            Interpolates and returns the electric field vector (E) at the specified spatial point.

    Args:
        ds (xr.Dataset): An xarray dataset containing the required Yee grid field components.
    """

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
    """
    A class representing a uniform electromagnetic field.

    Methods:
        B(x: np.ndarray) -> np.ndarray:
            Returns the uniform magnetic field vector, independent of position x.
        E(x: np.ndarray) -> np.ndarray:
            Returns the uniform electric field vector, independent of position x.
    """

    def __init__(
        self,
        B_0: np.ndarray | None = None,
        E_0: np.ndarray | None = None,
    ) -> None:
        self.B_0 = B_0 if B_0 is not None else np.array([0.0, 0.0, 0.0])
        self.E_0 = E_0 if E_0 is not None else np.array([0.0, 0.0, 0.0])

    def B(self, x: np.ndarray) -> np.ndarray:  # noqa: ARG002 pylint: disable=unused-argument
        return self.B_0

    def E(self, x: np.ndarray) -> np.ndarray:  # noqa: ARG002 pylint: disable=unused-argument
        return self.E_0


class DipoleField:
    """
    Represents a magnetic dipole field.

    Methods:
        B(r):
            Calculate the magnetic field vector at position r due to the dipole.

        E(r):
            Return the electric field vector at position r (always zero for static dipole).
    """

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

    def E(self, r):  # noqa: ARG002 pylint: disable=unused-argument
        return np.array([0.0, 0.0, 0.0])


class BorisIntegrator_python:
    """
    BorisIntegrator_python implements the Boris algorithm for integrating the motion of charged particles in electromagnetic fields.
    This class supports both Yee and non-Yee field interpolators, automatically selecting the appropriate interpolator based on the input dataset.

    Args:
        ds (xr.Dataset or FieldInterpolator_python or FieldInterpolatorYee_python):
            The dataset containing electromagnetic field data, or a pre-initialized field interpolator.
        q (float, optional):
            Particle charge in Coulombs. Defaults to the elementary charge (constants.e).
        m (float, optional):
            Particle mass in kilograms. Defaults to the electron mass (constants.m_e).

    Attributes:
        q (float): Particle charge.
        m (float): Particle mass.

    Methods:
        integrate(x0, v0, t_max, dt) -> pd.DataFrame:
            Integrates the particle trajectory using the Boris algorithm.
    """

    def __init__(self, ds, q=constants.e, m=constants.m_e) -> None:
        self.q = q
        self.m = m
        self._interpolator: FieldInterpolator_python | FieldInterpolatorYee_python
        if isinstance(ds, xr.Dataset):
            if {"bx1", "by1", "bz1", "ex1", "ey1", "ez1"} <= ds.data_vars.keys():
                self._interpolator = FieldInterpolatorYee_python(ds)
            else:
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
            x += 0.5 * dt * v
            B = self._interpolator.B(x)
            E = self._interpolator.E(x)
            v += qprime * E
            h = qprime * B
            s = 2 * h / (1 + np.linalg.norm(h) ** 2)
            v += np.cross(v + np.cross(v, h), s)
            v += qprime * E
            x += 0.5 * dt * v
            t += dt

        return pd.DataFrame(
            np.column_stack((times, positions, velocities)),
            columns=["time", "x", "y", "z", "vx", "vy", "vz"],
        )


class BorisIntegrator_f2py:
    """
    BorisIntegrator_f2py provides an interface for integrating charged particle trajectories
    using the Boris algorithm, with field interpolation via f2py-wrapped Fortran routines.

    Args:
        ds (xr.Dataset or FieldInterpolator_python or FieldInterpolatorYee_python):
            The dataset containing electromagnetic field data, or a pre-initialized field interpolator.
        q (float, optional):
            Particle charge in Coulombs. Defaults to the elementary charge (constants.e).
        m (float, optional):
            Particle mass in kilograms. Defaults to the electron mass (constants.m_e).

    Attributes:
        q (float): Particle charge.
        m (float): Particle mass.

    Methods:
        integrate(x0, v0, t_max, dt) -> pd.DataFrame:
            Integrates the particle trajectory using the Boris algorithm.
    """

    def __init__(self, df, q=constants.e, m=constants.m_e) -> None:
        _jrrle.particle_tracing_f2py.boris_init(q, m)
        if isinstance(df, xr.Dataset):
            self._interpolator = FieldInterpolatorYee_f2py(df)
        else:
            assert isinstance(df, FieldInterpolatorYee_f2py)
            self._interpolator = df

    def integrate(self, x0, v0, t_max, dt) -> pd.DataFrame:
        n_steps = int(t_max / dt) + 2  # add some extra space for round-off issues
        data = np.zeros((7, n_steps), dtype=np.float32, order="F")
        n_out = _jrrle.particle_tracing_f2py.boris_integrate(x0, v0, t_max, dt, data)
        return pd.DataFrame(
            data.T[:n_out], columns=["time", "x", "y", "z", "vx", "vy", "vz"]
        )


BorisIntegrator = BorisIntegrator_python
