from __future__ import annotations

import numpy as np
import scipy.interpolate  # type: ignore[import-untyped]
import xarray as xr

from ggcmpy.openggcm import _cotr_geo_sm_lat_lon, _cotr_sm_geo_lat_lon


def transform_geo_to_sm(da: xr.DataArray) -> xr.DataArray:
    """
    Transform GEO data to SM coordinates using Inverse Lookup Bilinear Interpolation.
    """

    lat_name = "lat" if "lat" in da.coords else "lats"
    lon_name = "lon" if "lon" in da.coords else "longs"

    # If in SM:
    if lat_name == "lats":
        return da

    plot_data = da.squeeze().transpose(lat_name, lon_name)
    lat_vals = plot_data.coords[lat_name].values
    lon_vals = plot_data.coords[lon_name].values
    z_vals = plot_data.values

    time_val = np.datetime64(np.atleast_1d(plot_data.coords["time"].values)[0])

    # 1. Prepare the source grid (GEO).
    if lat_vals[0] > lat_vals[-1]:
        lat_vals = lat_vals[::-1]
        z_vals = z_vals[::-1, :]

    if np.isclose(lon_vals[0] + 360, lon_vals[-1]):
        lon_ext = np.concatenate([lon_vals[:-1] - 360, lon_vals, lon_vals[1:] + 360])
        z_ext = np.concatenate([z_vals[:, :-1], z_vals, z_vals[:, 1:]], axis=1)
    else:
        lon_ext = np.concatenate([lon_vals - 360, lon_vals, lon_vals + 360])
        z_ext = np.concatenate([z_vals, z_vals, z_vals], axis=1)

    interp = scipy.interpolate.RegularGridInterpolator(
        (lat_vals, lon_ext),
        z_ext,
        method="linear",
        bounds_error=False,
        fill_value=None,
    )

    # 2. Create the target grid (SM).
    sm_lats = np.linspace(-90, 90, 181)
    sm_lons = np.linspace(0, 360, 361)
    sm_lon_mesh, sm_lat_mesh = np.meshgrid(sm_lons, sm_lats)

    # 3. Inverse lookup
    vec_cotr = np.vectorize(
        lambda y, x: _cotr_sm_geo_lat_lon(time_val, float(y), float(x))
    )
    geo_lat_target, geo_lon_target = vec_cotr(sm_lat_mesh, sm_lon_mesh)

    # 4. Bilinear interpolation
    query_points = np.stack([geo_lat_target, geo_lon_target], axis=-1)

    sm_z = interp(query_points)
    sm_z = np.clip(sm_z, np.nanmin(z_vals), np.nanmax(z_vals))

    return xr.DataArray(
        data=sm_z,
        coords={"mlat": sm_lats, "mlon": sm_lons},
        dims=["mlat", "mlon"],
        name=da.name,
        attrs=da.attrs,
    ).assign_coords(time=time_val)


def transform_sm_to_geo(da: xr.DataArray) -> xr.DataArray:
    """
    Transform SM data to GEO coordinates using Inverse Lookup Bilinear Interpolation.
    """

    lat_name = "lat" if "lat" in da.coords else "lats"
    lon_name = "lon" if "lon" in da.coords else "longs"

    # If in GEO:
    if lat_name == "lat":
        return da

    plot_data = da.squeeze().transpose(lat_name, lon_name)
    lat_vals = plot_data.coords[lat_name].values
    lon_vals = plot_data.coords[lon_name].values
    z_vals = plot_data.values

    time_val = np.datetime64(np.atleast_1d(plot_data.coords["time"].values)[0])

    # 1. Prepare the source grid (SM).
    if lat_vals[0] > lat_vals[-1]:
        lat_vals = lat_vals[::-1]
        z_vals = z_vals[::-1, :]

    if np.isclose(lon_vals[0] + 360, lon_vals[-1]):
        lon_ext = np.concatenate([lon_vals[:-1] - 360, lon_vals, lon_vals[1:] + 360])
        z_ext = np.concatenate([z_vals[:, :-1], z_vals, z_vals[:, 1:]], axis=1)
    else:
        lon_ext = np.concatenate([lon_vals - 360, lon_vals, lon_vals + 360])
        z_ext = np.concatenate([z_vals, z_vals, z_vals], axis=1)

    interp = scipy.interpolate.RegularGridInterpolator(
        (lat_vals, lon_ext),
        z_ext,
        method="linear",
        bounds_error=False,
        fill_value=None,
    )

    # 2. Create the target grid (GEO).
    geo_lats = np.linspace(-90, 90, 181)
    geo_lons = np.linspace(0, 360, 361)
    geo_lon_mesh, geo_lat_mesh = np.meshgrid(geo_lons, geo_lats)

    # 3. Inverse lookup
    vec_cotr = np.vectorize(
        lambda y, x: _cotr_geo_sm_lat_lon(time_val, float(y), float(x))
    )
    sm_lat_target, sm_lon_target = vec_cotr(geo_lat_mesh, geo_lon_mesh)

    # 4. Bilinear interpolation
    query_points = np.stack([sm_lat_target, sm_lon_target], axis=-1)

    geo_z = interp(query_points)
    geo_z = np.clip(geo_z, np.nanmin(z_vals), np.nanmax(z_vals))

    return xr.DataArray(
        data=geo_z,
        coords={"lat": geo_lats, "lon": geo_lons},
        dims=["lat", "lon"],
        name=da.name,
        attrs=da.attrs,
    ).assign_coords(time=time_val)
