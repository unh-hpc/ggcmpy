"""
boris_push.py

Boris particle pusher implementations.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ggcmpy import _openggcm, constants  # type: ignore[attr-defined]
from ggcmpy.tracing import emfields

# pylint: disable=C0103,I1101


class boris_push_python:
    """
    A class implementing the Boris particle pusher algorithm in pure Python.
    """

    def __init__(self, fields: emfields.emfields, q=constants.e, m=constants.m_e):
        self._fields = fields
        self._q = q
        self._m = m

    def push(
        self,
        prts_df: pd.DataFrame,
        t_final: float | None = None,
        max_steps: int | None = None,
        dt_max: float | None = None,
        dt_max_gyro: float = 0.1,
    ) -> pd.DataFrame:
        """
        Pushes the particles in `prts_df` from their current state.

        Args:
            prts_df (pd.DataFrame): DataFrame containing particle states with columns ['time', 'x', 'y', 'z', 'ux', 'uy', 'uz'].
            t_final (float | None): Final time to push the particles to.
            max_steps (int | None): Maximum number of steps to take.
            dt_max (float | None): Maximum time step for the integration.
            dt_max_gyro (float): Maximum time step as fraction of the gyroperiod.
        """
        prts_df = prts_df.copy()  # don't modify input

        qprime = 0.5 * self._q / self._m

        assert t_final is not None or max_steps is not None

        for n in range(len(prts_df)):
            step = 0
            while True:
                if t_final is not None and prts_df.loc[n, "time"] >= t_final:  # type: ignore[operator]
                    break

                if max_steps is not None and step >= max_steps:
                    break

                B = self._fields.B(prts_df.loc[n, ["x", "y", "z"]].to_numpy())
                u = prts_df.loc[n, ["ux", "uy", "uz"]].to_numpy()
                gamma = np.sqrt(1 + np.linalg.norm(u) ** 2)
                om_c = 2.0 * np.abs(qprime) * np.linalg.norm(B) / gamma
                dt = dt_max_gyro * 2.0 * np.pi / om_c
                if dt_max is not None:
                    dt = min(dt_max, dt)

                prts_df.loc[n, ["x", "y", "z"]] = self.push_x(
                    prts_df.loc[n, ["x", "y", "z"]].to_numpy(),
                    prts_df.loc[n, ["ux", "uy", "uz"]].to_numpy(),
                    0.5 * dt,
                )
                B = np.asarray(
                    self._fields.B(prts_df.loc[n, ["x", "y", "z"]].to_numpy())
                )
                E = np.asarray(
                    self._fields.E(prts_df.loc[n, ["x", "y", "z"]].to_numpy())
                )
                prts_df.loc[n, ["ux", "uy", "uz"]] = self.push_u(
                    prts_df.iloc[n][["ux", "uy", "uz"]].to_numpy(), E, B, qprime * dt
                )
                prts_df.loc[n, ["x", "y", "z"]] = self.push_x(
                    prts_df.loc[n, ["x", "y", "z"]].to_numpy(),
                    prts_df.loc[n, ["ux", "uy", "uz"]].to_numpy(),
                    0.5 * dt,
                )
                prts_df.loc[n, "time"] += dt
                step += 1

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


class particles_cxx(_openggcm.tracing.particles):  # type: ignore[misc]
    """Wrapper class for the C++ particles class, providing a convenient interface for particle data management."""

    def __new__(cls, df: pd.DataFrame) -> particles_cxx:
        return super().__new__(  # type: ignore[no-any-return] # pylint: disable=E1121
            cls,
            df["id"].to_numpy(),
            df["time"].to_numpy(),
            df[["x", "y", "z"]].to_numpy(),
            df[["ux", "uy", "uz"]].to_numpy(),
        )

    def __init__(self, df: pd.DataFrame) -> None:
        pass

    def to_dataframe(self) -> pd.DataFrame:
        _id, t, r, u = self.to_tuple()
        df = pd.DataFrame(
            np.column_stack((t, r, u)),
            columns=("time", "x", "y", "z", "ux", "uy", "uz"),
        )
        df["id"] = _id
        return df


class boris_push_cxx(_openggcm.tracing.boris):  # type: ignore[misc]
    """Wrapper class for the C++ boris class, providing a convenient interface for particle integration."""

    def push(
        self,
        prts_df: pd.DataFrame,
        t_final: float | None = None,
        max_steps: int | None = None,
        dt_max: float | None = None,
        dt_max_gyro: float = 0.1,
    ) -> pd.DataFrame:
        """
        Pushes the particles in `prts_df` from their current state.

        Args:
            prts_df (pd.DataFrame): DataFrame containing particle states with columns ['time', 'x', 'y', 'z', 'ux', 'uy', 'uz'].
            t_final (float | None): Final time to push the particles to.
            max_steps (int | None): Maximum number of steps to take.
            dt_max (float | None): Maximum time step for the integration.
            dt_max_gyro (float): Maximum time step as fraction of the gyroperiod.
        """
        prts = particles_cxx(prts_df)
        super().push(
            prts,
            t_final=t_final,
            max_steps=max_steps,
            dt_max=dt_max,
            dt_max_gyro=dt_max_gyro,
        )
        return prts.to_dataframe()
