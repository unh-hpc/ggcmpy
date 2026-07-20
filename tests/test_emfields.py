from __future__ import annotations

import numpy as np

from ggcmpy import tracing


def test_emfields_uniform():
    emfields = tracing.emfields.uniform(E_0=[1.0, 2.0, 3.0], B_0=[4.0, 5.0, 6.0])
    pos = [6.0, 7.0, 8.0]
    assert np.allclose(emfields.E(pos), [1.0, 2.0, 3.0])
    assert np.allclose(emfields.B(pos), [4.0, 5.0, 6.0])
