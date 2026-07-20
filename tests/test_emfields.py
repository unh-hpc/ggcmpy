from __future__ import annotations

import numpy as np
import scipy.constants  # type: ignore[import-untyped]

import ggcmpy.constants
from ggcmpy import tracing

R_E = ggcmpy.constants.radius_earth
m_E = ggcmpy.constants.dipole_moment_earth


def test_emfields_uniform():
    emfields = tracing.emfields.uniform(E_0=[1.0, 2.0, 3.0], B_0=[4.0, 5.0, 6.0])
    pos = [6.0, 7.0, 8.0]
    assert np.allclose(emfields.E(pos), [1.0, 2.0, 3.0])
    assert np.allclose(emfields.B(pos), [4.0, 5.0, 6.0])


def test_emfields_dipole():
    emfields = tracing.emfields.dipole(m=m_E)
    r = np.array([R_E, 0.0, 0.0])
    r_hat = r / np.linalg.norm(r)
    B_expected = (
        scipy.constants.mu_0
        / (4.0 * np.pi)
        * (3 * r_hat * np.dot(r_hat, m_E) - m_E)
        / (np.linalg.norm(r) ** 3)
    )
    assert np.allclose(emfields.E(r), [0.0, 0.0, 0.0])
    assert np.allclose(emfields.B(r), B_expected)
