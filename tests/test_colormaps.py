from __future__ import annotations

import matplotlib as mpl  # type: ignore[import-not-found]
import matplotlib.colors as mcolors  # type: ignore[import-not-found]
import numpy as np

import ggcmpy  # noqa: F401  # registers the colormaps


def test_openggcm_registered():
    assert "openggcm" in mpl.colormaps
    cmap = mpl.colormaps["openggcm"]
    assert isinstance(cmap, mcolors.Colormap)


def test_openggcm_endpoints_and_stops():
    cmap = mpl.colormaps["openggcm"]
    np.testing.assert_allclose(cmap(0.0)[:3], mcolors.to_rgb("#ffffff"), atol=1e-2)
    np.testing.assert_allclose(cmap(1.0)[:3], mcolors.to_rgb("#8b0000"), atol=1e-2)
    np.testing.assert_allclose(cmap(0.5)[:3], mcolors.to_rgb("#90ff90"), atol=1e-2)


def test_openggcm_usable_by_name():
    import matplotlib.pyplot as plt  # type: ignore[import-not-found]  # pylint: disable=import-outside-toplevel

    fig, ax = plt.subplots()
    im = ax.imshow(np.arange(4).reshape(2, 2), cmap="openggcm")
    assert im.get_cmap().name == "openggcm"
    plt.close(fig)
