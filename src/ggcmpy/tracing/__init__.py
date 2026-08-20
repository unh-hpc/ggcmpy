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


class boris_push_python:
    """
    A class implementing the Boris particle pusher algorithm in pure Python.
    """

    def __init__(self, fields: emfields.emfields, q=constants.e, m=constants.m_e):
        self._fields = fields
        self._q = q
        self._m = m

    def push(
        self, prts_df: pd.DataFrame, t_max: float, dt_max: float, gyro_max: float
    ) -> pd.DataFrame:
        qprime = 0.5 * self._q / self._m
        B = self._fields.B(prts_df.loc[0, ["x", "y", "z"]])  # type: ignore[arg-type]
        om_c = 2.0 * np.abs(qprime) * np.linalg.norm(B)
        dt = min(dt_max, gyro_max * 2.0 * np.pi / om_c)

        while prts_df.loc[0, "time"] < t_max:  # type: ignore[operator]
            prts_df.loc[0, ["x", "y", "z"]] = self.push_x(
                prts_df.loc[0, ["x", "y", "z"]].to_numpy(),
                prts_df.loc[0, ["ux", "uy", "uz"]].to_numpy(),
                0.5 * dt,
            )
            B = self._fields.B(prts_df.loc[0, ["x", "y", "z"]].to_numpy())
            E = self._fields.E(prts_df.loc[0, ["x", "y", "z"]].to_numpy())
            prts_df.loc[0, ["ux", "uy", "uz"]] = self.push_u(
                prts_df.iloc[0][["ux", "uy", "uz"]].to_numpy(), E, B, qprime * dt
            )
            prts_df.loc[0, ["x", "y", "z"]] = self.push_x(
                prts_df.loc[0, ["x", "y", "z"]].to_numpy(),
                prts_df.loc[0, ["ux", "uy", "uz"]].to_numpy(),
                0.5 * dt,
            )
            prts_df.loc[0, "time"] += dt

        return prts_df

    @staticmethod
    def push_x(x: np.ndarray, u: np.ndarray, dt: float) -> np.ndarray:
        inv_gamma = 1.0 / np.sqrt(1 + np.linalg.norm(u) ** 2)
        return x + dt * u * inv_gamma * constants.c  # type: ignore[no-any-return]

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


class BorisIntegratorBase:
    """
    Base class for Boris integrators.
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

    def integrate(self, x0, u0, t_max, dt_max=1.0, gyro_max=0.1) -> pd.DataFrame:
        boris_push = self._boris_push_cls(self._fields, self._q, self._m)

        prts_df = pd.DataFrame(
            np.array([[0.0, *x0, *u0]]),
            columns=["time", "x", "y", "z", "ux", "uy", "uz"],
        )
        snapshots = [prts_df.copy()]

        while prts_df.loc[0, "time"] < t_max:
            prts_df = boris_push.push(
                prts_df,
                prts_df.loc[0, "time"] + 1e-7,  # type: ignore[operator]
                dt_max,
                gyro_max,
            )
            snapshots.append(prts_df.copy())

        return pd.concat(snapshots, ignore_index=True)  # type: ignore[no-any-return]


class BorisIntegrator_python(BorisIntegratorBase):
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
        if isinstance(ds, xr.Dataset):
            if {"bx1", "by1", "bz1", "eflx", "efly", "eflz"} <= ds.data_vars.keys():
                fields: emfields.emfields = emfields.yee_cic_python(ds)
            else:
                fields = emfields.interpolator_python(ds)
        else:
            fields = ds  # assume it's already an emfields.emfields

        super().__init__(fields, q, m, boris_push_cls=boris_push_python)


class BorisIntegrator_f2py(BorisIntegratorBase):
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
            fields = FieldInterpolatorYee_f2py(df)
        else:
            assert isinstance(df, FieldInterpolatorYee_f2py)
            fields = df

        super().__init__(fields, q, m)

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
    """Wrapper class for the C++ particles class, providing a convenient interface for particle data management."""

    def __new__(cls, df: pd.DataFrame) -> particles_cxx:
        t = df["time"].to_numpy()
        r = df[["x", "y", "z"]].to_numpy()
        u = df[["ux", "uy", "uz"]].to_numpy()
        return super().__new__(cls, t, r, u)  # type: ignore[no-any-return] # pylint: disable=E1121

    def __init__(self, df: pd.DataFrame) -> None:
        pass

    def to_dataframe(self) -> pd.DataFrame:
        t, r, u = self.to_tuple()
        return pd.DataFrame(
            np.column_stack((t, r, u)),
            columns=("time", "x", "y", "z", "ux", "uy", "uz"),
        )


class boris_push_cxx(_openggcm.tracing.boris):  # type: ignore[misc]
    """Wrapper class for the C++ boris class, providing a convenient interface for particle integration."""

    def push(
        self, prts_df: pd.DataFrame, t_max: float, dt_max: float, gyro_max: float
    ) -> pd.DataFrame:
        prts = particles_cxx(prts_df)
        super().push(prts, t_max, dt_max, gyro_max)
        return prts.to_dataframe()


class BorisIntegrator_cxx(BorisIntegratorBase):
    """
    BorisIntegrator_cxx provides an interface for integrating charged particle trajectories
    using the Boris algorithm, with field interpolation via C++ routines.

    Args:
        df (xr.Dataset or emfields.interpolator_cxx or emfields.interpolator_yee_cxx):
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

    def __init__(self, df, q=constants.e, m=constants.m_e):
        if isinstance(df, xr.Dataset):
            fields = emfields.yee_cic_cxx(df)
        else:
            assert isinstance(
                df, (emfields.uniform_cxx, emfields.dipole_cxx, emfields.yee_cic_cxx)
            )
            fields = df

        super().__init__(fields, q, m, boris_push_cls=boris_push_cxx)


BorisIntegrator = BorisIntegrator_python
