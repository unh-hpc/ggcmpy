from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

import ggcmpy
from ggcmpy import openggcm

sample_time_array = np.array([[2010, 1, 1, 13, 0, 0, 100], [2010, 1, 1, 13, 1, 0, 100]])
sample_datetime64 = [
    np.datetime64("2010-01-01T13:00:00.100000000"),
    np.datetime64("2010-01-01T13:01:00.100000000"),
]


def test_to_dt64():
    assert np.all(
        openggcm.timearray_to_dt64(sample_time_array[0]) == sample_datetime64[0]
    )


def test_to_timearray():
    assert np.all(
        openggcm._dt64_to_timearray(sample_datetime64, dtype=np.int32)
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
    assert decoded_var.attrs == {}
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


def test_iopar():
    iof = xr.open_dataset(f"{ggcmpy.sample_dir}/coupling0001.iof.000030")
    iof = iof.isel(time=0)
    iox = ggcmpy.iono.iopar(iof)  # type: ignore[no-untyped-call]
    xr.testing.assert_allclose(iox["delbt"], iof["delbt"])


def test_epoch1966():
    testdate = np.datetime64("1967-02-03T04:05:06.100")
    dsecs = ggcmpy.openggcm.epoch1966(testdate)

    import datetime

    date0 = np.datetime64("1966-01-01T00:00:00.000").astype(datetime.datetime)
    assert np.isclose(dsecs, (testdate - date0).total_seconds())


def test_cotr_gse_to_mhd():
    testdate = np.datetime64("1967-02-03T04:05:06.100")

    r1 = np.array([1.0, 2.0, 3.0])
    r2 = ggcmpy.openggcm.cotr(testdate, "gse", "mhd", r1)
    assert np.allclose(r2, [-1.0, -2.0, 3.0])


def test_cotr_gse_to_gsm():
    # special time indicating no dipole tilt
    testdate = np.datetime64("1967-01-01T00:00:00.000")

    r1 = np.array([1.0, 2.0, 3.0])
    r2 = ggcmpy.openggcm.cotr(testdate, "gse", "gsm", r1)
    assert np.allclose(r2, [1.0, 2.0, 3.0])


def test_cpcp():
    iof = xr.open_dataset(ggcmpy.sample_dir / "sample_jrrle.iof.010491")
    cpcp = iof.ggcm.cpcp()
    assert np.isclose(cpcp, 85588.0)
    assert cpcp.attrs["long_name"] == "Cross Polar Cap Potential"
    assert cpcp.attrs["units"] == "V"


def test_cpcp_time_series():
    files = sorted((ggcmpy.sample_dir / "cir07_19970227_liang_norcm").glob("*.iof.*"))
    files = files[:10]  # limit to first 10 files for test speed
    iof = xr.open_mfdataset(files)
    cpcp = iof.ggcm.cpcp()
    assert cpcp.sizes == {"time": 10}


def test_lat_lon_to_cart():
    lat = np.array([0.0, 0.0, 90.0, -90.0])
    lon = np.array([0.0, 90.0, 0.0, 0.0])
    x, y, z = ggcmpy.openggcm._lat_lon_to_cart(lat, lon)
    assert np.allclose(x, [1.0, 0.0, 0.0, 0.0])
    assert np.allclose(y, [0.0, 1.0, 0.0, 0.0])
    assert np.allclose(z, [0.0, 0.0, 1.0, -1.0])


def test_cart_to_lat_lon():
    x = np.array([1.0, 0.0, 0.0, 0.0])
    y = np.array([0.0, 1.0, 0.0, 0.0])
    z = np.array([0.0, 0.0, 1.0, -1.0])
    lat, lon = ggcmpy.openggcm._cart_to_lat_lon(np.array([x, y, z]))
    assert np.allclose(lat, [0.0, 0.0, 90.0, -90.0])
    assert np.allclose(lon, [0.0, 90.0, 0.0, 0.0])


def test_cotr_geo_sm_lat_lon():
    lat_geo, lon_geo = 80.0, 10.0
    time = np.datetime64("2010-01-01T00:00:00")
    lat_sm, lon_sm = ggcmpy.openggcm._cotr_geo_sm_lat_lon(time, lat_geo, lon_geo)
    assert np.isclose(lat_sm, 76.7981)
    assert np.isclose(lon_sm, -126.36199)
    time += np.timedelta64(1, "h")  # add 1 hour
    lat_sm, lon_sm = ggcmpy.openggcm._cotr_geo_sm_lat_lon(time, lat_geo, lon_geo)
    assert np.isclose(lat_sm, 76.7981)
    assert np.isclose(lon_sm, -126.3620 + 1.0 / 24.0 * 360.0, atol=1.0)
