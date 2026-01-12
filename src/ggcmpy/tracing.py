"""
tracing.py

Particle tracing
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import constants  # type: ignore[import-untyped]

from ggcmpy import _jrrle  # type: ignore[attr-defined]


def load_fields(ds) -> None:
    """Load the field data into the Fortran backend for particle tracing.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing the necessary field variables:
        - bx, by, bz: Magnetic field components
        - ex, ey, ez: Electric field components
        - x, y, z: Grid coordinates
    """
    # FIXME, need to actually use electric fields when available
    _jrrle.particle_tracing_f2py.load(
        ds.bx, ds.by, ds.bz, ds.ex, ds.ey, ds.ez, ds.x, ds.y, ds.z
    )


def at(i: int, j: int, k: int, m: int) -> float:
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
    val = _jrrle.particle_tracing_f2py.at(i, j, k, m)
    assert isinstance(val, float)
    return val


def interpolate(x: float, y: float, z: float, m: int) -> float:
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
    val = _jrrle.particle_tracing_f2py.interpolate(x, y, z, m)
    assert isinstance(val, float)
    return val


class BorisIntegrator_python:
    def __init__(self, get_B, get_E, q=constants.e, m=constants.m_e) -> None:
        self.q = q
        self.m = m
        self.get_B = get_B
        self.get_E = get_E

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
            B = self.get_B(x)
            E = self.get_E(x)
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
    def __init__(self, get_B, get_E, q=constants.e, m=constants.m_e) -> None:  # noqa: ARG002
        _jrrle.particle_tracing_f2py.boris_init(q, m)

    def integrate(self, x0, v0, t_max, dt) -> pd.DataFrame:
        n_steps = int(t_max / dt) + 2  # add some extra space for round-off issues
        data = np.zeros((7, n_steps), dtype=np.float32, order="F")
        n_out = _jrrle.particle_tracing_f2py.boris_integrate(x0, v0, t_max, dt, data)
        return pd.DataFrame(
            data.T[:n_out], columns=["time", "x", "y", "z", "vx", "vy", "vz"]
        )


BorisIntegrator = BorisIntegrator_python
