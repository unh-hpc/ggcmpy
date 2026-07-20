from __future__ import annotations

import numpy as np
import scipy.constants  # type: ignore[import-untyped]


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
        B_0: np.ndarray | None = None,
        E_0: np.ndarray | None = None,
    ) -> None:
        self.B_0 = B_0 if B_0 is not None else np.array([0.0, 0.0, 0.0])
        self.E_0 = E_0 if E_0 is not None else np.array([0.0, 0.0, 0.0])

    def B(self, r: np.ndarray) -> np.ndarray:  # noqa: ARG002 pylint: disable=unused-argument
        return self.B_0

    def E(self, r: np.ndarray) -> np.ndarray:  # noqa: ARG002 pylint: disable=unused-argument
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

    def __init__(self, m):
        self.m = m

    def B(self, r):
        rhat = r / np.linalg.norm(r)
        return (
            scipy.constants.mu_0
            / (4 * np.pi)
            * (3 * np.dot(self.m, rhat) * rhat - self.m)
            / np.linalg.norm(r) ** 3
        )

    def E(self, r):  # noqa: ARG002 pylint: disable=unused-argument
        return np.array([0.0, 0.0, 0.0])


# Default implementations

uniform = uniform_python
dipole = dipole_python
