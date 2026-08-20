from __future__ import annotations

import numpy as np
from scipy.constants import (  # type: ignore[import-untyped]
    c,
    e,
    epsilon_0,
    m_e,
    m_p,
    mu_0,
)

__all__ = [
    "c",
    "dipole_moment_earth",
    "e",
    "epsilon_0",
    "m_e",
    "m_p",
    "mu_0",
    "radius_earth",
]

# pylint: disable=C0103

radius_earth = 6.371e6

dipole_moment_earth = np.array([0.0, 0.0, -7.8e22])
