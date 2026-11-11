"""
__init__.py

Particle tracing
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

import numpy as np
import xarray as xr

from . import emfields as emfields  # pylint: disable=useless-import-alias
from . import guiding_center as guiding_center  # pylint: disable=useless-import-alias
from . import integrator as integrator  # pylint: disable=useless-import-alias
from . import legacy as legacy  # pylint: disable=useless-import-alias

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
