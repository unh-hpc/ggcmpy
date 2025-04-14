"""
Choong Min Um <choongmin.um@unh.edu>

This module plots data from an OpenGGCM output file.
"""

from __future__ import annotations

import argparse
import sys

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


# Plot the data.
def plot(
    file: str, var: str, lats_max: int, lats_min: int, spacing: int, mlt: bool
) -> None:
    with xr.open_dataset(file) as ds:
        ds.coords["colats"] = 90 - ds.coords["lats"]
        da = ds[var].sel(lats=slice(int(lats_max), int(lats_min)))
        lon = da.coords["longs"]
        fig, ax = plt.subplots(
            subplot_kw={"projection": "polar", "theta_offset": np.pi / 2}
        )
        range_theta = range(0, 360, 15)
        plt.thetagrids(range_theta, grids_theta_mlt if mlt else grids_theta_deg)
        range_r, grids_r, coord_ns = get_plot_params(lats_max, lats_min, spacing, da)
        plt.rgrids(range_r, grids_r)
        ax.set_title(title_dict[var] or "")
        mesh = ax.contourf(
            np.deg2rad(lon),
            coord_ns,
            da.T,
            cmap="bwr",
            levels=np.linspace(-1e-6, 1e-6, 51),
            extend="both",
        )
        fig.colorbar(mesh)
        plt.show()
        return


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
        plot(
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
