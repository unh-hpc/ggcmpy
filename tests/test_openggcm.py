from __future__ import annotations

from ggcmpy import openggcm
import numpy as np

sample_time_array = np.array([[2010, 1, 1, 13, 0, 0, 100], [2010, 1, 1, 13, 1, 0, 100]])
sample_datetime64 = [
    np.datetime64("2010-01-01T13:00:00.100000000"),
    np.datetime64("2010-01-01T13:01:00.100000000"),
]


def test_to_dt64():
    assert np.all(openggcm._time_array_to_dt64(sample_time_array) == sample_datetime64)


def test_to_time_array():
    assert np.all(
        openggcm._dt64_to_time_array(sample_datetime64, dtype=np.int32)
        == sample_time_array
    )
