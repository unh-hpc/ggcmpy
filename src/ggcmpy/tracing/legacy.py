"""
legacy.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from ggcmpy import (  # type: ignore[attr-defined]
    _jrrle,
    constants,
)
from ggcmpy.tracing import emfields, integrator

# pylint: disable=C0103,I1101


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


class BorisIntegratorBase:
    """
    Base class for Boris integrators.
    """

    def __init__(
        self,
        fields: emfields.emfields,
        q=constants.e,
        m=constants.m_e,
        integrator_boris_cls=None,
    ):
        self._fields = fields
        self._q = q
        self._m = m
        self._integrator_boris_cls = integrator_boris_cls

    def integrate(self, x0, u0, t_final, dt_max=1.0, dt_max_gyro=0.1) -> pd.DataFrame:
        integrator_boris: integrator.boris_base = self._integrator_boris_cls(
            self._fields, self._q, self._m
        )

        prts_df = pd.DataFrame(
            np.array([[0, 0.0, *x0, *u0]]),
            columns=["id", "time", "x", "y", "z", "ux", "uy", "uz"],
        )

        return integrator_boris.integrate(
            prts_df,
            t_final=t_final,
            snapshot_interval_steps=1,
            dt_max=dt_max,
            dt_max_gyro=dt_max_gyro,
        )


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
        integrate(x0, v0, t_final, dt) -> pd.DataFrame:
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

        super().__init__(fields, q, m, integrator_boris_cls=integrator.boris_python)


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
        integrate(x0, v0, t_final, dt) -> pd.DataFrame:
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

        super().__init__(fields, q, m, integrator_boris_cls=integrator.boris_cxx)


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
        integrate(x0, v0, t_final, dt) -> pd.DataFrame:
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

    def integrate(self, x0, u0, t_final, dt_max=1.0, dt_max_gyro=0.1) -> pd.DataFrame:
        n_steps = int(t_final / dt_max) + 2  # add some extra space for round-off issues
        data = np.zeros((7, n_steps), dtype=np.float32, order="F")
        n_out = _jrrle.particle_tracing_f2py.boris_integrate(
            x0, u0, t_final, dt_max, dt_max_gyro, data
        )
        return pd.DataFrame(
            data.T[:n_out], columns=["time", "x", "y", "z", "ux", "uy", "uz"]
        )
