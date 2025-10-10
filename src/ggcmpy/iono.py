"""
iono.py

Ionosphere diagnostics
"""

from __future__ import annotations

import xarray as xr

from ggcmpy import _jrrle  # type: ignore[attr-defined]


def gradpt(pot):
    """Compute the gradient of potential field on a spherical surface
    (at ionospheric height) using the Fortran routine.

    Parameters
    ----------
    pot : xarray.DataArray
        2D array of the potential field on a spherical surface.
        Coordinates must be 'longs' and 'lats'.

    Returns
    -------
    xarray.Dataset
        Dataset containing the two components of the gradient field:
        'et' (theta-component) and 'ep' (phi-component).
    """
    pot = pot.transpose("longs", "lats")  # ensure correct order for Fortran

    et, ep = _jrrle.ground_perturbation_m.gradpt(pot.values)
    return xr.Dataset(
        {"et": (("longs", "lats"), et), "ep": (("longs", "lats"), ep)},
        coords=pot.coords,
    )


def iopar(ds):
    """Postprocess ionospheric quantities using the Fortran routine.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing the necessary input variables:
        - fac_tot: FAC
        - pot: Potential
        - sigp: Pedersen conductivity
        - sigh: Hall conductivity

    Returns
    -------
    xarray.Dataset
        Dataset containing the computed ionospheric parameters:
        - ep: Electric field phi-component
        - et: Electric field theta-component
        - ctiot: Total ionospheric current theta-component
        - ctiop: Total ionospheric current phi-component
        - tau:
        - ctaup: current phi-component
        - ctaut: current theta-component
        - cpolp: current phi-component
        - cpolt: current theta-component
        - delbr: ground magnetic perturbation radial component
        - delbp: ground magnetic perturbation phi-component
        - delbt: ground magnetic perturbation theta-component
        - xjh: Joule heating rate"""
    if "time" in ds.sizes:
        msg = "Dataset must be 2D without time dimension. Use .isel(time=0) to select a time slice."
        raise ValueError(msg)
    ds = ds.transpose("longs", "lats")  # ensure correct order for Fortran

    ep, et, ctiot, ctiop, tau, ctaup, ctaut, cpolp, cpolt, delbr, delbp, delbt, xjh = (
        _jrrle.ground_perturbation_m.iopar(ds.fac_tot, ds.pot, ds.sigp, ds.sigh)
    )
    fields = {
        "ep": ep,
        "et": et,
        "ctiot": ctiot,
        "ctiop": ctiop,
        "tau": tau,
        "ctaup": ctaup,
        "ctaut": ctaut,
        "cpolp": cpolp,
        "cpolt": cpolt,
        "delbr": delbr,
        "delbp": delbp,
        "delbt": delbt,
        "xjh": xjh,
    }
    return xr.Dataset(
        {name: (("longs", "lats"), data) for name, data in fields.items()},
        coords=ds.coords,
    )
