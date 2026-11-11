"""
guiding_center.py

Calculation of guiding center positions along particle trajectories.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ggcmpy import constants
from ggcmpy.tracing import emfields

# pylint: disable=C0103


class calc_guiding_center:
    """
    Calculate the guiding center position for a given trajectory.
    """

    def __init__(
        self,
        fields: emfields.emfields,
        q: float,
        m: float,
        scale_B: float = 1.0,
        scale_E: float = 1.0,
        scale_r: float = 1.0,
    ):
        self.fields = fields
        self.q = q
        self.m = m
        self.scale_B = scale_B
        self.scale_E = scale_E
        self.scale_r = scale_r
        self.q_over_m = q / m

    def __call__(self, *args: Any) -> Any:
        """
        Calculate the guiding center position for a given row of trajectory data.
        """
        if len(args) == 1:
            row = args[0]
            assert isinstance(row, pd.Series), "Expected a pandas Series"
            x = np.array([row.x, row.y, row.z])
            u = np.array([row.ux, row.uy, row.uz])

            x_gc, y_gc, z_gc, v_par, mu, E_kin = self._to_guiding_center(x, u)
            return pd.Series(
                [row.id, row.time, x_gc, y_gc, z_gc, v_par, mu, E_kin],
                index=["id", "time", "x_gc", "y_gc", "z_gc", "v_par", "mu", "E_kin"],
            )

        assert len(args) == 2, "Expected either a single row or two arguments (x, u)"
        x, u = args
        return self._to_guiding_center(x, u)

    def _to_guiding_center(
        self, x: np.ndarray, u: np.ndarray
    ) -> tuple[float, float, float, float, float, float]:
        """
        Calculate the guiding center position for a given position and velocity.

        Returns (x_gc, y_gc, z_gc, u_par, mu, E_kin), with E_kin in eV.
        """
        B = self.scale_B * self.fields.B(x)
        B_mag = np.linalg.norm(B)
        b_unit = B / B_mag
        x_gc = x + constants.c * np.cross(u, b_unit) / (
            self.q_over_m * B_mag * self.scale_r
        )

        # get B, E at guiding center position
        B = self.scale_B * self.fields.B(x_gc)
        B_mag = np.linalg.norm(B)
        b_unit = B / B_mag
        u_par = np.dot(u, b_unit)
        u_perp = u - u_par * b_unit
        gamma = np.sqrt(1.0 + np.linalg.norm(u) ** 2)
        E_kin = (gamma - 1.0) * self.m * constants.c**2
        mu = 0.5 * self.m * constants.c**2 * np.linalg.norm(u_perp) ** 2 / B_mag
        return x_gc[0], x_gc[1], x_gc[2], u_par, mu, E_kin / constants.e


def to_guiding_center(
    df: pd.DataFrame,
    q: float,
    m: float,
    fields: emfields.emfields,
    scale_B: float = 1.0,
    scale_E: float = 1.0,
    scale_r: float = 1.0,
) -> pd.DataFrame:
    """
    Calculate the guiding center position for each row of a trajectory DataFrame.
    """
    calc_gc = calc_guiding_center(
        fields, q, m, scale_B=scale_B, scale_E=scale_E, scale_r=scale_r
    )
    return pd.DataFrame(df.apply(calc_gc, axis=1, result_type="expand"))
