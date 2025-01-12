from __future__ import annotations

from ggcmpy import openggcm
import xarray as xr
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


def test_decode_openggcm_variable():
    var = xr.Variable(("time", "time_array"), sample_time_array)
    var_dt64 = xr.Variable("time", sample_datetime64)

    decoded_var = openggcm._decode_openggcm_variable(var, "name")
    assert decoded_var.equals(var)

    var.attrs["units"] = "time_array"
    decoded_var = openggcm._decode_openggcm_variable(var, "name")
    assert decoded_var.equals(var_dt64)
    assert var.attrs == {"units": "time_array"}  # original var not changed
    assert decoded_var.encoding == {"units": "time_array", "dtype": var.dtype}
