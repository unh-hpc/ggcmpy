# pylint: disable=import-outside-toplevel, cyclic-import
from __future__ import annotations

import argparse
import sys
import warnings
from typing import Any

import cartopy.feature as cfeature  # pylint: disable=import-error
import matplotlib.pyplot as plt  # type: ignore[import-not-found]
import numpy as np
import xarray as xr

# Define longitude choices for plot tick labels.
grids_theta_mlt = (
    "12",
    "14",
    "16",
    "18",
    "20",
    "22",
    "0",
    "2",
    "4",
    "6",
    "8",
    "10",
)

grids_theta_deg = (
    "0",
    "30",
    "60",
    "90",
    "120",
    "150",
    "180",
    "210",
    "240",
    "270",
    "300",
    "330",
)


class InvalidLatitudesException(Exception):
    """Exception is raised when invalid latitude boundaries are provided."""


def lats_invalid() -> None:
    """Print an error message and exit the program if latitude arguments are invalid."""
    msg = (
        "Invalid latitudes! Please enter one of the following:\n"
        "1) '-n' or '-s' as an option, OR\n"
        "2) 0 <= lats_min < lats_max <= 90 for the northern hemisphere, OR\n"
        "3) -90 <= lats_min < lats_max <= 0 for the southern hemisphere."
    )
    sys.exit(msg)


def get_plot_params(
    lats_max: int, lats_min: int, spacing: int
) -> tuple[range, tuple[str, ...]]:
    """Generate the radial limits and labels for the polar plot grid."""
    if lats_min >= 0:
        # Northern hemisphere
        range_r = range(int(90 - lats_max), int(90 - lats_min), spacing)
        grids_r = tuple(f"{90 - r}" for r in range_r)
    elif lats_min < 0:
        # Southern hemisphere
        range_r = range(int(90 + lats_min), int(90 + lats_max), spacing)
        grids_r = tuple(f"{r - 90}" for r in range_r)
    else:
        raise InvalidLatitudesException

    return range_r, grids_r


def draw_coastlines_polar(ax: Any, lats_min: int, time: np.datetime64 | None) -> None:
    """Draw geographic coastlines transformed into Solar Magnetic (SM) coordinates."""
    from .openggcm import _cotr_geo_sm_lat_lon

    if time is None:
        warnings.warn(
            "Coastlines require a time variable for SM transformation. Skipping.",
            stacklevel=2,
        )
        return

    feature = cfeature.COASTLINE.with_scale("110m")

    # Vectorize the transformation function.
    vec_cotr = np.vectorize(
        lambda lat, lon: _cotr_geo_sm_lat_lon(time, float(lat), float(lon))
    )

    for geom in feature.geometries():
        lines = geom.geoms if geom.geom_type == "MultiLineString" else [geom]
        for line in lines:
            coords = np.asarray(line.coords)
            plot_lats, plot_lons = vec_cotr(coords[:, 1], coords[:, 0])
            theta = np.deg2rad(plot_lons)
            r = 90.0 - plot_lats if lats_min >= 0 else 90.0 + plot_lats
            ax.plot(theta, r, color="black", linewidth=0.4)


def draw_magnetometers(
    ax: Any,
    lats_min: int,
    time: np.datetime64 | None,
    highlight: str | None = None,
    network: str = "AL",
) -> None:
    """Plot ground magnetometer stations in SM coordinates."""
    from .openggcm import (
        CANOPUS_MAGNETOMETERS,
        MAGNETOMETERS,
        _cotr_geo_sm_lat_lon,
    )

    if time is None:
        warnings.warn(
            "Magnetometers require a time variable for SM transformation. Skipping.",
            stacklevel=2,
        )
        return

    stations_dict = CANOPUS_MAGNETOMETERS if network.upper() == "CL" else MAGNETOMETERS

    for name, station in stations_dict.items():
        stn_lat = float(station["lat"])  # type: ignore[arg-type]
        stn_lon = float(station["lon"])  # type: ignore[arg-type]

        plot_lat, plot_lon = _cotr_geo_sm_lat_lon(time, stn_lat, stn_lon)
        theta = np.deg2rad(plot_lon)
        r = 90.0 - plot_lat if lats_min >= 0 else 90.0 + plot_lat

        if name == highlight:
            ax.scatter(theta, r, s=40, color="green", zorder=6)
        else:
            ax.scatter(theta, r, s=20, color="black", zorder=5)

        ax.annotate(
            name,
            xy=(theta, r),
            xytext=(2.0, 2.0),
            textcoords="offset points",
            fontsize=7,
        )


def plot_from_dataarray(
    da: xr.DataArray,
    lats_max: int,
    lats_min: int,
    spacing: int,
    mlt: bool,
    levels: Any | None = None,
    cmap: str = "bwr",
    extend: str = "both",
    coastlines: bool = False,
    stations: bool = False,
    timestamp: bool = False,
    highlight_station: str | None = None,
    network: str = "AL",
    **kwargs: Any,
) -> None:
    """Core function to format and plot 2D polar data in SM coordinates"""

    # Safely extract formatting metadata from the Xarray dataset.
    da_name = str(da.name) if da.name is not None else "Variable"
    long_name = da.attrs.get("long_name", da_name)
    units = da.attrs.get("units", "")
    plot_title = f"{long_name} [{units}]" if units else long_name

    # Extract single temporal value for SM transformations (coastlines/stations).
    time_val: np.datetime64 | None = None
    if "time" in da.coords:
        time_val = np.datetime64(np.atleast_1d(da.coords["time"].values)[0])

    # Guarantee the DataArray is in SM coordinates before slicing.
    from .coord_transform import transform_geo_to_sm

    da = transform_geo_to_sm(da)

    # Determine the coordinate names.
    lat_name = "mlat" if "mlat" in da.coords else "lats"
    lon_name = "mlon" if "mlon" in da.coords else "longs"

    # Ensure the latitude slice direction is valid.
    lat_vals = da.coords[lat_name].values
    is_descending = lat_vals[0] > lat_vals[-1]

    if is_descending:
        da_sliced = da.sel({lat_name: slice(lats_max, lats_min)})
    else:
        da_sliced = da.sel({lat_name: slice(lats_min, lats_max)})

    # Ensure consistent matrix orientation: (Latitude, Longitude).
    plot_data = da_sliced.squeeze().transpose(lat_name, lon_name)
    z_vals = plot_data.values
    lat_vals = plot_data.coords[lat_name].values
    lon_vals = plot_data.coords[lon_name].values

    # Prepare the grid.
    r_1d = 90.0 - lat_vals if lats_min >= 0 else 90.0 + lat_vals
    theta_1d = np.deg2rad(lon_vals)
    theta_grid, r_grid = np.meshgrid(theta_1d, r_1d)
    plot_z = z_vals

    # Initialize the plot.
    fig, ax = plt.subplots(
        subplot_kw={"projection": "polar", "theta_offset": np.pi / 2}
    )

    if levels is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            abs_max = np.nanmax(np.abs(plot_z))

        if np.isnan(abs_max) or abs_max == 0:
            abs_max = 1.0

        levels = np.linspace(-abs_max, abs_max, 51)

    if coastlines:
        draw_coastlines_polar(ax, lats_min, time=time_val)

    if stations:
        draw_magnetometers(
            ax,
            lats_min=lats_min,
            time=time_val,
            highlight=highlight_station,
            network=network,
        )

    # Hardcode the timestamp printing.
    if timestamp and time_val is not None:
        ax.text(np.pi, 48, str(time_val), ha="center", va="center", fontsize=10)

    # Apply axis labels and grids.
    range_r, grids_r = get_plot_params(lats_max, lats_min, spacing)
    plt.thetagrids(np.arange(0, 360, 30), grids_theta_mlt if mlt else grids_theta_deg)
    plt.rgrids(range_r, grids_r)

    ax.set_title(plot_title, pad=8, fontsize=13)
    ax.set_axisbelow(False)

    # Plot the final transformed data.
    mesh = ax.contourf(
        theta_grid,
        r_grid,
        plot_z,
        cmap=cmap,
        levels=levels,
        extend=extend,
        **kwargs,
    )

    plot_rmax = 90 - lats_min if lats_min >= 0 else 90 + lats_max
    ax.set_rmax(plot_rmax)
    ax.set_ylim(0, plot_rmax)

    fig.colorbar(mesh, pad=0.08, shrink=0.85)


def get_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Plot OpenGGCM polar data.")
    group = parser.add_mutually_exclusive_group()

    group.add_argument("-n", "--north", action="store_true")
    group.add_argument("-s", "--south", action="store_true")

    parser.add_argument("file")
    parser.add_argument("var")
    parser.add_argument("lats_max", type=int, nargs="?")
    parser.add_argument("lats_min", type=int, nargs="?")
    parser.add_argument("spacing", type=int, nargs="?", default=10)
    parser.add_argument("mlt", nargs="?", default="true")
    parser.add_argument("--coastlines", action="store_true")
    parser.add_argument("--stations", action="store_true")
    parser.add_argument(
        "--network",
        type=str,
        choices=["AL", "CL", "al", "cl"],
        default="AL",
        help="Magnetometer network to plot (AL or CL)",
    )
    parser.add_argument(
        "--highlight",
        type=str,
        default=None,
        help="Station code to highlight (e.g., FCC)",
    )
    parser.add_argument("--timestamp", action="store_true")

    args = parser.parse_args()

    # Reject overlapping constraints.
    if (args.north or args.south) and (args.lats_max or args.lats_min):
        lats_invalid()

    # Set default plot boundaries if none are provided.
    if args.south:
        args.lats_max = -50 if args.lats_max is None else args.lats_max
        args.lats_min = -90 if args.lats_min is None else args.lats_min
    else:
        args.lats_max = 90 if args.lats_max is None else args.lats_max
        args.lats_min = 50 if args.lats_min is None else args.lats_min

    return args


def main() -> None:
    args = get_args()
    try:
        with xr.open_dataset(args.file) as ds:
            ds[args.var].ggcm.plot(
                lats_max=args.lats_max,
                lats_min=args.lats_min,
                spacing=args.spacing,
                mlt=(args.mlt.lower() == "true"),
                coastlines=args.coastlines,
                stations=args.stations,
                network=args.network,
                highlight_station=args.highlight,
                timestamp=args.timestamp,
            )
            plt.show()
    except InvalidLatitudesException:
        lats_invalid()


if __name__ == "__main__":
    main()
