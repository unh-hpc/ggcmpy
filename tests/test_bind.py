from __future__ import annotations

import numpy as np

from ggcmpy import _openggcm  # type: ignore[attr-defined]


def test_xt_fixed_from_python():
    a = [1.0, 2.0, 3.0]
    b = _openggcm.xt_fixed_from_python(a)
    assert np.array_equal(b, a)
