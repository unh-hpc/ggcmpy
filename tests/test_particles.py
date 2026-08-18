from __future__ import annotations

import numpy as np
from ggcmpy._openggcm.tracing import particle  # type: ignore[import-not-found]


def test_particle_ctor():
    prt = particle(
        x=np.array([0.0, 0.0, 0.0]), u=np.array([1.0, 0.0, 0.0]), q=-1.0, m=1.0
    )
    assert np.allclose(prt.x, np.array([0.0, 0.0, 0.0]))
    assert np.allclose(prt.u, np.array([1.0, 0.0, 0.0]))
    prt.x = [1.0, 2.0, 3.0]
    assert np.allclose(prt.x, np.array([1.0, 2.0, 3.0]))
    assert np.isclose(prt.q, -1.0)
    assert np.isclose(prt.m, 1.0)
    assert np.isclose(prt.gamma, np.sqrt(2.0))
