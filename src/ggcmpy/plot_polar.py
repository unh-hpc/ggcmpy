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
        coord_ns = da.ggcm.coords["colats"]
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


def plot_from_dataarray(
    da: xr.DataArray,
    lats_max: int,
    lats_min: int,
    spacing: int,
    mlt: bool,
    levels=None,
    cmap="bwr",
    extend="both",
    **kwargs,
) -> None:
    da_sliced = da.sel(lats=slice(lats_max, lats_min))
    if levels is None:
        abs_max = np.abs(da_sliced.values).max()
        levels = np.linspace(-abs_max, abs_max, 51)

    plot_title = f"{da.attrs['long_name']} [{da.attrs['units']}]"
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
        **kwargs,
    )
    fig.colorbar(mesh)
    plt.show()


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
    if args.south:
        if args.lats_max is None:
            args.lats_max = -50
        if args.lats_min is None:
            args.lats_min = -90
    else:  # Default to northern hemisphere.
        if args.lats_max is None:
            args.lats_max = 90
        if args.lats_min is None:
            args.lats_min = 50

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
            )
    except InvalidLatitudesException:
        lats_invalid()


if __name__ == "__main__":
    main()
