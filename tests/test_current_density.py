from __future__ import annotations

import numpy as np
import xarray as xr

import ggcmpy.openggcm  # noqa: F401  # registers the .ggcm accessor


def test_ggcm_current_density():
    rng = np.random.default_rng(0)
    shape = (3, 4)
    ds = xr.Dataset(
        {
            "xjx": (("lats", "longs"), rng.standard_normal(shape)),
            "xjy": (("lats", "longs"), rng.standard_normal(shape)),
            "xjz": (("lats", "longs"), rng.standard_normal(shape)),
        }
    )

    xj = ds.ggcm.current_density()

    expected = np.sqrt(ds.xjx**2 + ds.xjy**2 + ds.xjz**2)
    np.testing.assert_allclose(xj.to_numpy(), expected.to_numpy())
    assert xj.name == "xj"
    assert xj.attrs["long_name"] == "Current density magnitude"
    assert "xj" not in ds  # the dataset itself is not modified
