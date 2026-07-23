from __future__ import annotations

import numpy as np
import pytest
import scipy.constants  # type: ignore[import-untyped]
import xarray as xr
from numpy.typing import ArrayLike

import ggcmpy.constants
import ggcmpy.tracing

R_E = ggcmpy.constants.radius_earth
m_E = ggcmpy.constants.dipole_moment_earth


b_grid = [("bx", ("x", "y", "z")), ("by", ("x", "y", "z")), ("bz", ("x", "y", "z"))]
e_grid = [("ex", ("x", "y", "z")), ("ey", ("x", "y", "z")), ("ez", ("x", "y", "z"))]

b1_grid = [
    ("bx1", ("x_nc", "y", "z")),
    ("by1", ("x", "y_nc", "z")),
    ("bz1", ("x", "y", "z_nc")),
]
e1_grid = [
    ("eflx", ("x", "y_nc", "z_nc")),
    ("efly", ("x_nc", "y", "z_nc")),
    ("eflz", ("x_nc", "y_nc", "z")),
]


def make_coords() -> dict[str, np.ndarray]:
    coords = {f"{dir}_nc": np.linspace(-1.0, 1.0, 11) for dir in ["x", "y", "z"]}
    coords |= {
        dir: 0.5 * (coords[f"{dir}_nc"][:-1] + coords[f"{dir}_nc"][1:])
        for dir in ["x", "y", "z"]
    }
    return coords


class emfields_test:
    def B(self, r: ArrayLike) -> np.ndarray:
        r = np.asarray(r)
        return np.array(
            [
                2 * r[0] + 3 * r[1] + 5 * r[2],
                3 * r[0] + 5 * r[1] + 7 * r[2],
                5 * r[0] + 7 * r[1] + 11 * r[2],
            ]
        )

    def E(self, r: ArrayLike) -> np.ndarray:
        r = np.asarray(r)
        return np.array(
            [
                7 * r[0] + 11 * r[1] + 13 * r[2],
                11 * r[0] + 13 * r[1] + 17 * r[2],
                13 * r[0] + 17 * r[1] + 19 * r[2],
            ]
        )


@pytest.mark.parametrize(
    "emfields_uniform",
    [
        ggcmpy.tracing.emfields.uniform_python,
        ggcmpy.tracing.emfields.uniform_cxx,
    ],
)
def test_emfields_uniform(emfields_uniform):
    emfields = emfields_uniform(E_0=[1.0, 2.0, 3.0], B_0=[4.0, 5.0, 6.0])
    pos = [6.0, 7.0, 8.0]
    assert np.allclose(emfields.E(pos), [1.0, 2.0, 3.0])
    assert np.allclose(emfields.B(pos), [4.0, 5.0, 6.0])


@pytest.mark.parametrize(
    "emfields_dipole",
    [
        ggcmpy.tracing.emfields.dipole_python,
        ggcmpy.tracing.emfields.dipole_cxx,
    ],
)
def test_emfields_dipole(emfields_dipole):
    emfields = emfields_dipole(m=m_E)
    r = np.array([R_E, 0.0, 0.0])
    r_hat = r / np.linalg.norm(r)
    B_expected = (
        scipy.constants.mu_0
        / (4.0 * np.pi)
        * (3 * r_hat * np.dot(r_hat, m_E) - m_E)
        / (np.linalg.norm(r) ** 3)
    )
    assert np.allclose(emfields.E(r), [0.0, 0.0, 0.0])
    assert np.allclose(emfields.B(r), B_expected)


@pytest.mark.parametrize(
    "interpolator",
    [
        ggcmpy.tracing.emfields.interpolator_python,
        ggcmpy.tracing.FieldInterpolator_f2py,
    ],
)
def test_emfields_interpolator(interpolator):
    emfields = emfields_test()

    coords = make_coords()
    emfields_cc = xr.Dataset(
        ggcmpy.tracing.make_vector_field(b_grid, coords, emfields.B)
        | ggcmpy.tracing.make_vector_field(e_grid, coords, emfields.E),
        coords=coords,
    )

    emfields_ip = interpolator(emfields_cc)
    point = (0.1, 0.25, 0.3)
    # since the original field is linear, the interpolation should be exact
    assert np.allclose(emfields_ip.B(point), emfields.B(point))
    assert np.allclose(emfields_ip.E(point), emfields.E(point))


@pytest.mark.parametrize(
    "interpolator",
    [
        ggcmpy.tracing.emfields.interpolator_yee_python,
        ggcmpy.tracing.FieldInterpolatorYee_f2py,
        ggcmpy.tracing.emfields.yee_cic_cxx,
    ],
)
def test_emfields_yee(interpolator):
    emfields = emfields_test()

    coords = make_coords()
    emfields_yee = xr.Dataset(
        ggcmpy.tracing.make_vector_field(b1_grid, coords, emfields.B)
        | ggcmpy.tracing.make_vector_field(e1_grid, coords, emfields.E),
        coords=coords,
    )

    emfields_ip = interpolator(emfields_yee)
    point = np.array([0.1, 0.25, 0.3])
    # since the original field is linear, the interpolation should be exact
    assert np.allclose(emfields_ip.B(point), emfields.B(point), atol=1e-7)
    assert np.allclose(emfields_ip.E(point), emfields.E(point), atol=1e-7)
