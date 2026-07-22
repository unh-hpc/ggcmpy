from __future__ import annotations

import numpy as np

from ggcmpy import _openggcm  # type: ignore[attr-defined]


def test_xt_fixed_from_python():
    a = [1.0, 2.0, 3.0]
    b = _openggcm.xt_fixed_from_python(a)
    assert np.array_equal(b, a)


def test_test_ndarray():
    import sys

    a = np.array([1.0, 2.0, 3.0])
    assert sys.getrefcount(a) == 2
    b = _openggcm.test_ndarray(a)
    assert sys.getrefcount(a) == 4
    assert b[1] == 2.0
    a[1] = 5.0
    assert b[1] == 5.0
    del b
    assert sys.getrefcount(a) == 2


def test_test_xtadapt():
    import sys

    a = np.array([1.0, 2.0, 3.0])
    b = _openggcm.test_xtadapt(a)
    assert b[1] == 2.0
    a[1] = 5.0
    assert b[1] == 5.0
    del b
    assert sys.getrefcount(a) == 2
