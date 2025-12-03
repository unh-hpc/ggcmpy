from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import ggcmpy
from ggcmpy.timeseries import read_ggcm_solarwind_file


def test_read_ggcm_solarwind_file():
    rr = read_ggcm_solarwind_file(
        ggcmpy.sample_dir / "cir07_19970227_liang_norcm/tmp.minvar/om1997.rr"
    )
    assert set(rr.columns) == {"om1997.rr"}
    assert pd.api.types.is_datetime64_ns_dtype(rr.index)
    assert np.isclose(rr["om1997.rr"].mean(), 4.8191107)


def test_read_ggcm_solarwind_directory():
    bfield = ggcmpy.timeseries.read_ggcm_solarwind_directory(
        ggcmpy.sample_dir / "cir07_19970227_liang_norcm/tmp.minvar", glob="om1997.b*gse"
    )
    assert set(bfield.columns) == {"om1997.bxgse", "om1997.bygse", "om1997.bzgse"}
    assert pd.api.types.is_datetime64_ns_dtype(bfield.index)
    assert np.isclose(bfield["om1997.bzgse"].mean(), -1.87348468)


def test_store_to_pyspedas_pandas():
    pyspedas = pytest.importorskip("pyspedas")

    bfield = ggcmpy.timeseries.read_ggcm_solarwind_directory(
        ggcmpy.sample_dir / "cir07_19970227_liang_norcm/tmp.minvar", glob="om1997.b*gse"
    )
    ggcmpy.timeseries.store_to_pyspedas(bfield)
    assert "om1997.bxgse" in pyspedas.tplot_names()


def test_store_to_pyspedas_xarray():
    pyspedas = pytest.importorskip("pyspedas")

    rng = np.random.default_rng()

    data = xr.DataArray(
        data=rng.random(100),
        dims=["time"],
        coords={"time": pd.date_range("2023-01-01", periods=100, freq="min")},
        attrs={"long_name": "Cross Polar Cap Potential", "name": "CPCP", "units": "V"},
        name="cpcp",
    )
    ggcmpy.timeseries.store_to_pyspedas(data)
    assert "cpcp" in pyspedas.tplot_names()

    dset = xr.Dataset(data_vars={"var1": data, "var2": data * 2})
    ggcmpy.timeseries.store_to_pyspedas(dset)
    assert "var1" in pyspedas.tplot_names()
    assert "var2" in pyspedas.tplot_names()
