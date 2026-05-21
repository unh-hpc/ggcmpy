# pylint: disable=import-outside-toplevel, cyclic-import
from __future__ import annotations

import argparse
import sys
from typing import Any

import cartopy.feature as cfeature  # pylint: disable=import-error
import matplotlib.pyplot as plt  # type: ignore[import-not-found]
import numpy as np
import xarray as xr

# Define longitude choices.
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
    """Handle latitude exceptions."""


# Exit the program if given invalid latitudes.
def lats_invalid() -> None:
    msg = (
        "Invalid latitudes! Please enter one of the following:\n"
        "1) '-n' or '-s' as an option, OR\n"
        "2) 0 <= lats_min < lats_max <= 90 for the northern hemisphere, OR\n"
        "3) -90 <= lats_min < lats_max <= 0 for the southern hemisphere."
    )
    sys.exit(msg)


# Return a tuple of the plot parameters.
def get_plot_params(
    lats_max: int, lats_min: int, spacing: int, da: xr.DataArray
) -> tuple[range, tuple[str, ...], xr.DataArray]:
    if 0 <= lats_min < lats_max <= 90:
        range_r = range(int(90 - lats_max), int(lats_max - lats_min), spacing)
        grids_r = tuple(
            f"{lat}" for lat in range(int(lats_max), int(lats_min), -spacing)
        )
        coord_ns = da.ggcm.coords["colats"]
    elif -90 <= lats_min < lats_max <= 0:
        range_r = range(int(lats_min), int(lats_max), spacing)
        grids_r = tuple(
            f"{lat}" for lat in range(int(lats_min), int(lats_max), spacing)
        )
        coord_ns = da.coords["lats"]
    else:
        raise InvalidLatitudesException

    return range_r, grids_r, coord_ns


def draw_coastlines_polar(ax: Any, lats_min: int, time: np.datetime64) -> None:
    from .openggcm import _cotr_geo_sm_lat_lon

    feature = cfeature.COASTLINE.with_scale("110m")

    for geom in feature.geometries():
        lines = geom.geoms if geom.geom_type == "MultiLineString" else [geom]
        for line in lines:
            coords = np.asarray(line.coords)
            lon = coords[:, 0]
            lat = coords[:, 1]

            # Transform coastline geographic coordinates to Solar Magnetic (SM).
            mlats = []
            mlons = []
            for geo_lat, geo_lon in zip(lat, lon, strict=True):
                mlat, mlon = _cotr_geo_sm_lat_lon(time, float(geo_lat), float(geo_lon))
                mlats.append(mlat)
                mlons.append(mlon)

            mlats_array = np.array(mlats)
            mlons_array = np.array(mlons)

            # Mask based on the new magnetic latitudes.
            mask = mlats_array >= lats_min
            if not np.any(mask):
                continue

            # Plot using the transformed magnetic coordinates.
            theta = np.deg2rad(mlons_array[mask])
            r = np.clip(90 - mlats_array[mask], 0, 90 - lats_min)

            ax.plot(theta, r, color="black", linewidth=0.4)


def draw_magnetometers(
    ax: Any,
    time: np.datetime64,
    highlight: str | None = None,
    network: str = "AL",
) -> None:
    from .openggcm import (
        CANOPUS_MAGNETOMETERS,
        MAGNETOMETERS,
        _cotr_geo_sm_lat_lon,
    )

    stations_dict = CANOPUS_MAGNETOMETERS if network.upper() == "CL" else MAGNETOMETERS

    for name, station in stations_dict.items():
        # Transform geographic coordinates to Solar Magnetic (SM) for plotting.
        mlat, mlon = _cotr_geo_sm_lat_lon(
            time,
            float(station["lat"]),  # type: ignore[arg-type]
            float(station["lon"]),  # type: ignore[arg-type]
        )

        # Use the transformed mlon and mlat.
        theta = np.deg2rad(mlon)
        r = 90 - mlat

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
    fig, ax = plt.subplots(
        subplot_kw={"projection": "polar", "theta_offset": np.pi / 2}
    )

    plot_title = f"{da.attrs['long_name']} [{da.attrs['units']}]"
    da_sliced = da.sel(lats=slice(lats_max, lats_min))

    if levels is None:
        abs_max = np.abs(da_sliced.values).max()
        levels = np.linspace(-abs_max, abs_max, 51)

    if coastlines:
        time_val = da_sliced.coords["time"].values
        draw_coastlines_polar(ax, lats_min, time=time_val)

    if stations:
        time_val = da_sliced.coords["time"].values
        draw_magnetometers(
            ax, time=time_val, highlight=highlight_station, network=network
        )

    if timestamp:
        current_time = da_sliced.time.values
        ax.text(np.pi, 48, current_time, ha="center", va="center", fontsize=10)

    range_theta = range(0, 360, 30)
    range_r, grids_r, coord_ns = get_plot_params(lats_max, lats_min, spacing, da_sliced)

    plt.thetagrids(range_theta, grids_theta_mlt if mlt else grids_theta_deg)
    plt.rgrids(range_r, grids_r)

    ax.set_title(plot_title, pad=8, fontsize=13)
    ax.set_axisbelow(False)

    lon = da_sliced.coords["longs"]

    mesh = ax.contourf(
        np.deg2rad(lon),
        coord_ns,
        da_sliced.T,
        cmap=cmap,
        levels=levels,
        extend=extend,
        **kwargs,
    )

    fig.colorbar(mesh, pad=0.08, shrink=0.85)


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot OpenGGCM polar data")
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

    if (args.north or args.south) and (args.lats_max or args.lats_min):
        lats_invalid()

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
    except InvalidLatitudesException:
        lats_invalid()


if __name__ == "__main__":
    main()
