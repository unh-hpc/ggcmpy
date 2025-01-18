from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from ggcmpy import openggcm

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


@pytest.mark.parametrize(
    ("dims", "data_time_array", "data_datetime64"),
    [
        (("time", "time_array"), sample_time_array, sample_datetime64),
        (("time_array",), sample_time_array[0], sample_datetime64[0]),
    ],
)
def test_AmieTimeArrayCoder(dims, data_time_array, data_datetime64):
    var = xr.Variable(dims, data_time_array)
    var_dt64 = xr.Variable(dims[:-1], data_datetime64)

    coder = openggcm.AmieTimeArrayCoder()
    decoded_var = coder.decode(var, "name")
    assert decoded_var.equals(var)  # type: ignore[no-untyped-call]

    var.attrs["units"] = "time_array"
    decoded_var = coder.decode(var, "name")
    assert decoded_var.equals(var_dt64)  # type: ignore[no-untyped-call]
    assert var.attrs == {"units": "time_array"}  # original var not changed
    assert decoded_var.encoding == {"units": "time_array", "dtype": var.dtype}


@pytest.mark.parametrize(
    ("dims", "data_time_array", "data_datetime64"),
    [
        (("time", "time_array"), sample_time_array, sample_datetime64),  # time array
        (("time_array",), sample_time_array[0], sample_datetime64[0]),  # scalar time
    ],
)
def test_encode_openggcm_variable(dims, data_time_array, data_datetime64):
    var = xr.Variable(dims, data_time_array)
    var_dt64 = xr.Variable(dims[:-1], data_datetime64)

    var.attrs["units"] = "time_array"
    coder = openggcm.AmieTimeArrayCoder()
    decoded_var = coder.decode(var, "name")
    assert decoded_var.equals(var_dt64)  # type: ignore[no-untyped-call]
    encoded_var = coder.encode(decoded_var)
    assert encoded_var.equals(var)  # type: ignore[no-untyped-call]
    assert encoded_var.encoding == {}


def test_coords():
    ds = xr.Dataset(
        {"longs": np.linspace(-180, 180, 61)}, {"lats": np.linspace(90, -90, 181)}
    )
    assert np.allclose(ds.ggcm.coords["mlts"], np.linspace(0, 24, 61))
    assert np.allclose(ds.ggcm.coords["colats"], np.linspace(0, 180, 181))
    assert np.allclose(ds.ggcm.mlts, np.linspace(0, 24, 61))
    assert np.allclose(ds.ggcm.colats, np.linspace(0, 180, 181))
