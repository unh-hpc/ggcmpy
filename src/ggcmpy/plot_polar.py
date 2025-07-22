from __future__ import annotations

import argparse
import sys
from typing import Any

import matplotlib.pyplot as plt  # type: ignore[import-not-found]
import numpy as np
import xarray as xr

# Define longitude choices.
grids_theta_mlt = (
    "12",
    "",
    "14",
    "",
    "16",
    "",
    "18",
    "",
    "20",
    "",
    "22",
    "",
    "0",
    "",
    "2",
    "",
    "4",
    "",
    "6",
    "",
    "8",
    "",
    "10",
    "",
)
grids_theta_deg = (
    "0",
    "",
    "30",
    "",
    "60",
    "",
    "90",
    "",
    "120",
    "",
    "150",
    "",
    "180",
    "",
    "210",
    "",
    "240",
    "",
    "270",
    "",
    "300",
    "",
    "330",
    "",
)

# Define title choices.
title_dict = {
    "fac_tot": "Total Field-Aligned Current [A/m^2]",
    "pot": "Potential [V]",
}


class InvalidLatitudesException(Exception):
    """Custom exception for invalid latitudes."""


# Exit the program if given invalid latitudes.
def lats_invalid() -> None:
    error_message = (
        "Invalid latitudes! Please enter one of the following:\n"
        "1) '-n' or '-s' as an option, OR\n"
        "2) 0 <= lats_min < lats_max <= 90 for the northern hemisphere, OR\n"
        "3) -90 <= lats_min < lats_max <= 0 for the southern hemisphere."
    )
    sys.exit(error_message)


# Return a tuple of the plot parameters.
def get_plot_params(
    lats_max: int, lats_min: int, spacing: int, da: xr.DataArray
) -> tuple[range, tuple[str, ...], xr.DataArray]:
    if 0 <= lats_min < lats_max <= 90:
        range_r = range(int(90 - lats_max), int(lats_max - lats_min), spacing)
        grids_r = tuple(
            f"{lat}" for lat in range(int(lats_max), int(lats_min), -spacing)
        )
        coord_ns = da.coords["colats"]
    elif -90 <= lats_min < lats_max <= 0:
        range_r = range(int(lats_min), int(lats_max), spacing)
        grids_r = tuple(
            f"{lat}" for lat in range(int(lats_min), int(lats_max), spacing)
        )
        coord_ns = da.coords["lats"]
    else:
        error_message = (
            "Invalid latitudes! Please enter one of the following:\n"
            "1) 0 <= lats_min < lats_max <= 90 for the northern hemisphere, OR\n"
            "2) -90 <= lats_min < lats_max <= 0 for the southern hemisphere."
        )
        raise InvalidLatitudesException(error_message)
    return range_r, grids_r, coord_ns


# Render Matplotlib.
def render_plot(
    da: xr.DataArray,
    plot_title: str,
    lats_max: int,
    lats_min: int,
    spacing: int,
    mlt: bool,
    levels=None,
    cmap="bwr",
    extend="both",
) -> None:
    da_sliced = da.sel(lats=slice(lats_max, lats_min))
    if levels is None:
        abs_max = np.abs(da_sliced.values).max()
        levels = np.linspace(-abs_max, abs_max, 51)

    lon = da_sliced.coords["longs"]

    fig, ax = plt.subplots(
        subplot_kw={"projection": "polar", "theta_offset": np.pi / 2}
    )
    range_theta = range(0, 360, 15)
    plt.thetagrids(range_theta, grids_theta_mlt if mlt else grids_theta_deg)

    range_r, grids_r, coord_ns = get_plot_params(lats_max, lats_min, spacing, da_sliced)
    plt.rgrids(range_r, grids_r)
    ax.set_title(plot_title)

    mesh = ax.contourf(
        np.deg2rad(lon),
        coord_ns,
        da_sliced.T,
        cmap=cmap,
        levels=levels,
        extend=extend,
    )
    fig.colorbar(mesh)
    plt.show()


# Plot the data.
def plot_from_file(
    file: str,
    var: str,
    lats_max: int,
    lats_min: int,
    spacing: int,
    mlt: bool,
    **kwargs: Any,
) -> None:
    with xr.open_dataset(file) as ds:
        ds.coords["colats"] = 90 - ds.coords["lats"]
        da_variable = ds[var]
        plot_title = title_dict.get(var, "")

        render_plot(
            da=da_variable,
            plot_title=plot_title,
            lats_max=lats_max,
            lats_min=lats_min,
            spacing=spacing,
            mlt=mlt,
            **kwargs,
        )
        return


# Plot the data using the Xarray accessor.
def plot_from_dataarray(
    da: xr.DataArray,
    lats_max: int,
    lats_min: int,
    spacing: int,
    mlt: bool,
    **kwargs: Any,
) -> None:
    da = da.copy(deep=True)
    da.coords["colats"] = 90 - da.coords["lats"]
    name_as_key = da.name
    plot_title = title_dict.get(name_as_key, "") if isinstance(name_as_key, str) else ""

    render_plot(
        da=da,
        plot_title=plot_title,
        lats_max=lats_max,
        lats_min=lats_min,
        spacing=spacing,
        mlt=mlt,
        **kwargs,
    )


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("This module plots data from an OpenGGCM output file.")
    )

    # Define mutually exclusive options.
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-n", "--north", action="store_true")
    group.add_argument("-s", "--south", action="store_true")

    # Define positional arguments.
    parser.add_argument("file", help="data file")
    parser.add_argument("var", help="variable to be plotted")
    parser.add_argument("lats_max", type=int, nargs="?", help="maximum latitude")
    parser.add_argument("lats_min", type=int, nargs="?", help="minimum latitude")
    parser.add_argument(
        "spacing",
        type=int,
        nargs="?",
        default=10,
        help="spacing between grids",
    )
    parser.add_argument(
        "mlt",
        nargs="?",
        default="true",
        help="magnetic local time coordinates on/off",
    )

    args = parser.parse_args()

    # Validate the arguments.
    if (args.north or args.south) and (args.lats_max or args.lats_min):
        lats_invalid()

    # Define conditionally default latitudes.
    if args.north:
        if args.lats_max is None:
            args.lats_max = 90
        if args.lats_min is None:
            args.lats_min = 50
    if args.south:
        if args.lats_max is None:
            args.lats_max = -50
        if args.lats_min is None:
            args.lats_min = -90

    return args


def main() -> None:
    args = get_args()
    try:
        plot_from_file(
            args.file,
            args.var,
            args.lats_max,
            args.lats_min,
            args.spacing,
            args.mlt.lower() == "true",
        )
    except InvalidLatitudesException:
        lats_invalid()


if __name__ == "__main__":
    main()
