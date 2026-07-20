from __future__ import annotations

import numpy as np


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

    def B(self, x: np.ndarray) -> np.ndarray:  # noqa: ARG002 pylint: disable=unused-argument
        return self.B_0

    def E(self, x: np.ndarray) -> np.ndarray:  # noqa: ARG002 pylint: disable=unused-argument
        return self.E_0


uniform = uniform_python
