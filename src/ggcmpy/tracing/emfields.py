from __future__ import annotations

import numpy as np
import scipy.constants  # type: ignore[import-untyped]
import xarray as xr
from numpy.typing import ArrayLike

# pylint: disable=C0103


class uniform_python:
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
        B_0: ArrayLike | None = None,
        E_0: ArrayLike | None = None,
    ) -> None:
        self.B_0 = np.asarray(B_0) if B_0 is not None else np.array([0.0, 0.0, 0.0])
        self.E_0 = np.asarray(E_0) if E_0 is not None else np.array([0.0, 0.0, 0.0])

    def B(self, r: ArrayLike) -> np.ndarray:  # noqa: ARG002 pylint: disable=unused-argument
        return self.B_0

    def E(self, r: ArrayLike) -> np.ndarray:  # noqa: ARG002 pylint: disable=unused-argument
        return self.E_0


class dipole_python:
    """
    Represents a magnetic dipole field.

    Methods:
        B(r):
            Calculate the magnetic field vector at position r due to the dipole.

        E(r):
            Return the electric field vector at position r (always zero for static dipole).
    """

    def __init__(self, m: ArrayLike):
        self.m = m

    def B(self, r: ArrayLike):
        r = np.asarray(r)
        rhat = r / np.linalg.norm(r)
        return (
            scipy.constants.mu_0
            / (4 * np.pi)
            * (3 * np.dot(self.m, rhat) * rhat - self.m)
            / np.linalg.norm(r) ** 3
        )

    def E(self, r: ArrayLike):  # noqa: ARG002 pylint: disable=unused-argument
        return np.array([0.0, 0.0, 0.0])


class interpolator_python:
    """
    A class for interpolating electromagnetic field components from a xarray.Dataset.
    """

    def __init__(self, ds: xr.Dataset) -> None:
        assert {"bx", "by", "bz", "ex", "ey", "ez"} <= ds.data_vars.keys()
        self._ds = ds

    def B(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [self._interpolate(self._ds[fld], point) for fld in ["bx", "by", "bz"]]
        )

    def E(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [self._interpolate(self._ds[fld], point) for fld in ["ex", "ey", "ez"]]
        )

    def _interpolate(self, da: xr.DataArray, point: np.ndarray) -> float:
        val = da.interp(dict(zip(da.dims, point, strict=True))).to_numpy()
        return float(val)


class interpolator_yee_python:
    """
    A class for interpolating electromagnetic field components from a xarray.Dataset on a Yee grid.
    """

    def __init__(self, ds: xr.Dataset) -> None:
        assert {"bx1", "by1", "bz1", "eflx", "efly", "eflz"} <= ds.data_vars.keys()
        self._ds = ds

    def B(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [self._interpolate(self._ds[fld], point) for fld in ["bx1", "by1", "bz1"]]
        )

    def E(self, point: np.ndarray) -> np.ndarray:
        return np.array(
            [
                self._interpolate(self._ds[fld], point)
                for fld in ["eflx", "efly", "eflz"]
            ]
        )

    def _interpolate(self, da: xr.DataArray, point: np.ndarray) -> float:
        val = da.interp(dict(zip(da.dims, point, strict=True))).to_numpy()
        return float(val)


# Default implementations

uniform = uniform_python
dipole = dipole_python
