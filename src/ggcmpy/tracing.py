"""
tracing.py

Particle tracing
"""

from __future__ import annotations

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
        ds.bx, ds.by, ds.bz, ds.xjx, ds.xjy, ds.xjz, ds.x, ds.y, ds.z
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
