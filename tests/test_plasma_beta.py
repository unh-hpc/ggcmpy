from __future__ import annotations

import numpy as np
import xarray as xr

import ggcmpy.openggcm  # noqa: F401  # registers the .ggcm accessor
from ggcmpy import constants


def test_ggcm_plasma_beta_units():
    # pp in pPa, bx/by/bz in nT: |B| = 5 nT, p = 10 pPa
    ds = xr.Dataset(
        {
            "pp": (("lats",), np.array([10.0])),
            "bx": (("lats",), np.array([3.0])),
            "by": (("lats",), np.array([4.0])),
            "bz": (("lats",), np.array([0.0])),
        }
    )

    beta = ds.ggcm.plasma_beta()

    # beta = 2 mu_0 p / B^2 with p = 10e-12 Pa and B = 5e-9 T
    expected = 2 * constants.mu_0 * 10e-12 / (5e-9) ** 2
    np.testing.assert_allclose(beta.to_numpy(), [expected])
    assert beta.name == "beta"
    assert beta.attrs["long_name"] == "Plasma beta"


def test_ggcm_plasma_beta_2d_and_3d():
    rng = np.random.default_rng(1)
    for dims in [("x", "y"), ("x", "y", "z")]:
        shape = (3,) * len(dims)
        pp = rng.uniform(0.1, 10.0, shape)
        bx, by, bz = (rng.uniform(-5.0, 5.0, shape) for _ in range(3))
        ds = xr.Dataset(
            {
                "pp": (dims, pp),
                "bx": (dims, bx),
                "by": (dims, by),
                "bz": (dims, bz),
            }
        )

        beta = ds.ggcm.plasma_beta()

        assert beta.dims == dims
        assert beta.shape == shape
        expected = 2 * constants.mu_0 * pp * 1e-12 / ((bx**2 + by**2 + bz**2) * 1e-18)
        np.testing.assert_allclose(beta.to_numpy(), expected)
