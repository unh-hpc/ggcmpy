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

from ggcmpy import (  # type: ignore[attr-defined]
    _jrrle,
    _openggcm,
)

from . import emfields

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
        # fld = xr.apply_ufunc(lambda x, y, z: vector_field(np.array([x, y, z]))[d], *crds, vectorize=True)
        flds[fld_name] = (dims, fld)

    return flds


def discretize_emfields_cc(coords: dict[str, np.ndarray], fields: Any) -> xr.Dataset:
    b_grid = [("bx", ("x", "y", "z")), ("by", ("x", "y", "z")), ("bz", ("x", "y", "z"))]
    e_grid = [("ex", ("x", "y", "z")), ("ey", ("x", "y", "z")), ("ez", ("x", "y", "z"))]

    return xr.Dataset(
        make_vector_field(b_grid, coords, fields.B)
        | make_vector_field(e_grid, coords, fields.E),
        coords=coords,
    )


def discretize_emfields_yee(coords: dict[str, np.ndarray], fields: Any) -> xr.Dataset:
    b1_grid = [
        ("bx1", ("x_nc", "y", "z")),
        ("by1", ("x", "y_nc", "z")),
        ("bz1", ("x", "y", "z_nc")),
    ]
    e1_grid = [
        ("eflx", ("x", "y_nc", "z_nc")),
        ("efly", ("x_nc", "y", "z_nc")),
        ("eflz", ("x_nc", "y_nc", "z")),
    ]

    return xr.Dataset(
        make_vector_field(b1_grid, coords, fields.B)
        | make_vector_field(e1_grid, coords, fields.E),
        coords=coords,
    )


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
            ds.eflx,
            ds.efly,
            ds.eflz,
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


class BorisIntegrator_python:
    """
    BorisIntegrator_python implements the Boris algorithm for integrating the motion of charged particles in electromagnetic fields.
    This class supports both Yee and non-Yee field interpolators, automatically selecting the appropriate interpolator based on the input dataset.

    Args:
        ds (xr.Dataset or emfields.interpolator_python or emfields.interpolator_yee_python):
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
        self._interpolator: (
            emfields.uniform | emfields.interpolator_python | emfields.yee_cic_python
        )
        if isinstance(ds, xr.Dataset):
            if {"bx1", "by1", "bz1", "eflx", "efly", "eflz"} <= ds.data_vars.keys():
                self._interpolator = emfields.yee_cic_python(ds)
            else:
                self._interpolator = emfields.interpolator_python(ds)
        else:
            self._interpolator = ds  # assume it's already an emfields

    def integrate(self, x0, u0, t_max, dt_max=1.0, gyro_max=0.1) -> pd.DataFrame:
        t = 0.0
        x = x0.copy()
        u = u0.copy()
        qprime = 0.5 * self.q / self.m

        B = self._interpolator.B(x)
        times, positions, momenta = [], [], []
        while t < t_max:
            times.append(t)
            positions.append(x.copy())
            momenta.append(u.copy())

            om_c = np.linalg.norm(
                2.0 * qprime * B
            )  # gyro frequency (based on previous B)
            dt = min(dt_max, gyro_max * 2.0 * np.pi / om_c)

            x = self.push_x(x, u, 0.5 * dt)
            B = self._interpolator.B(x)
            E = self._interpolator.E(x)
            u = self.push_u(u, E, B, qprime * dt)
            x = self.push_x(x, u, 0.5 * dt)
            t += dt

        return pd.DataFrame(
            np.column_stack((times, positions, momenta)),
            columns=["time", "x", "y", "z", "ux", "uy", "uz"],
        )

    @staticmethod
    def push_x(x: np.ndarray, u: np.ndarray, dt: float):
        gamma = np.sqrt(1 + np.linalg.norm(u) ** 2)
        return x + dt * u * constants.c / gamma

    @staticmethod
    def push_u(u: np.ndarray, E: np.ndarray, B: np.ndarray, dq: float):
        um = u + dq * E / constants.c

        # Rotation due to magnetic field
        root = dq / np.sqrt(1.0 + np.linalg.norm(um) ** 2)
        tau = root * B
        tau_norm = 1.0 / (1.0 + np.linalg.norm(tau) ** 2)
        up = np.array(
            [
                (
                    (1.0 + (tau[0]) ** 2 - (tau[1]) ** 2 - (tau[2]) ** 2) * um[0]
                    + (2.0 * tau[0] * tau[1] + 2.0 * tau[2]) * um[1]
                    + (2.0 * tau[0] * tau[2] - 2.0 * tau[1]) * um[2]
                )
                * tau_norm,
                (
                    (2.0 * tau[0] * tau[1] - 2.0 * tau[2]) * um[0]
                    + (1.0 - (tau[0]) ** 2 + (tau[1]) ** 2 - (tau[2]) ** 2) * um[1]
                    + (2.0 * tau[1] * tau[2] + 2.0 * tau[0]) * um[2]
                )
                * tau_norm,
                (
                    (2.0 * tau[0] * tau[2] + 2.0 * tau[1]) * um[0]
                    + (2.0 * tau[1] * tau[2] - 2.0 * tau[0]) * um[1]
                    + (1.0 - (tau[0]) ** 2 - (tau[1]) ** 2 + (tau[2]) ** 2) * um[2]
                )
                * tau_norm,
            ]
        )

        return up + dq * E / constants.c


class BorisIntegrator_f2py:
    """
    BorisIntegrator_f2py provides an interface for integrating charged particle trajectories
    using the Boris algorithm, with field interpolation via f2py-wrapped Fortran routines.

    Args:
        ds (xr.Dataset or emfields.interpolator_python or emfields.interpolator_yee_python):
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

    def integrate(self, x0, u0, t_max, dt_max=1.0, gyro_max=0.1) -> pd.DataFrame:
        n_steps = int(t_max / dt_max) + 2  # add some extra space for round-off issues
        data = np.zeros((7, n_steps), dtype=np.float32, order="F")
        n_out = _jrrle.particle_tracing_f2py.boris_integrate(
            x0, u0, t_max, dt_max, gyro_max, data
        )
        return pd.DataFrame(
            data.T[:n_out], columns=["time", "x", "y", "z", "ux", "uy", "uz"]
        )


class particles_cxx(_openggcm.tracing.particles):  # type: ignore[misc]
    def to_dataframe(self) -> pd.DataFrame:
        t, r, u = self.to_tuple()
        return pd.DataFrame(
            np.column_stack((t, r, u)),
            columns=("time", "x", "y", "z", "ux", "uy", "uz"),
        )


class BorisIntegrator_cxx:
    def __init__(self, df, q=constants.e, m=constants.m_e):
        from . import emfields

        if isinstance(df, xr.Dataset):
            self._emfields = emfields.yee_cic_cxx(df)
        else:
            assert isinstance(df, (emfields.uniform_cxx, emfields.yee_cic_cxx))
            self._emfields = df

        self._q = q
        self._m = m

    def integrate(self, x0, u0, t_max, dt_max=1.0, gyro_max=0.1) -> pd.DataFrame:
        boris = _openggcm.tracing.boris(self._emfields)

        prt = _openggcm.tracing.particle(x=x0, u=u0, q=self._q, m=self._m)

        particles = particles_cxx()
        boris.push(prt, t_max, dt_max, gyro_max, particles)

        return particles.to_dataframe()


BorisIntegrator = BorisIntegrator_python
