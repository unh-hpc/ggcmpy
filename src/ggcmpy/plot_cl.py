from __future__ import annotations

import os
import re
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt  # type: ignore[import-not-found]
import numpy as np
import pandas as pd
import xarray as xr
from spacepy import coordinates as coord  # type: ignore[import-not-found]
from spacepy.time import Ticktock  # type: ignore[import-not-found]

dir_cur = Path.cwd()
dir_par = dir_cur.parent
dir_par_name = dir_cur.parent.name

dir_input = Path("../inp")
dir_cl_observed = Path()
dir_model = Path("../target")

file_bz = dir_input / "wi.bzgse"
file_rr = dir_input / "wi.rr"
file_cl = dir_cl_observed / "CN_K0_MARI_2629956.txt"

# Geographic coordinates (lat, lon) for magnetometer stations
magnetometers_geo = [
    ("BACK", 57.707, 265.794-360.0),
    ("CONT", 65.754, 248.750-360.0),
    ("DAWS", 64.048, 220.890-360.0),
    ("ESKI", 61.106, 265.950-360.0),
    ("FCHU", 58.763, 265.920-360.0),
    ("FSIM", 61.756, 238.770-360.0),
    ("FSMI", 60.017, 248.050-360.0),
    ("GILL", 56.376, 265.360-360.0),
    ("ISLL", 53.856, 265.340-360.0),
    ("MCMU", 56.657, 248.790-360.0),
    ("PINA", 50.199, 263.960-360.0),
    ("RABB", 58.222, 256.320-360.0),
    ("RANK", 62.824, 267.890-360.0),
    ("TALO", 69.540, 266.450-360.0),
]

# Unit conversion factor for the model output
UNIT_CONVERSION = [
    1,
    1e9,
]


def load_data_wind(filepath: Any, value_col: Any):
    """Loads and parses time-series data from WIND files"""
    try:
        df_wind = pd.read_csv(
            filepath,
            sep=r"\s+",
            header=None,
            names=[
                "year",
                "month",
                "day",
                "hour",
                "minute",
                "second_float",
                value_col,
            ],
        )
        time = pd.to_datetime(df_wind[["year", "month", "day", "hour", "minute"]])
        df_wind["time"] = time + pd.to_timedelta(df_wind["second_float"], unit="s")
        return df_wind.set_index("time")[[value_col]]
    except OSError as e:
        print(f"Error loading {filepath}: {e}")  # noqa: T201
        return pd.DataFrame(columns=[value_col])


def load_data_cl_observed(filepath: Any):
    """Loads, cleans, and interpolates CL data"""
    try:
        with Path.open(filepath, encoding="utf-8") as f:
            for i, line in enumerate(f):
                if "dd-mm-yyyy" in line:
                    skip_rows = i + 1
                    break
            else:
                err_msg = "Header line 'dd-mm-yyyy' not found."
                raise ValueError(err_msg)

        df_obs = pd.read_csv(
            filepath,
            sep=r"\s+",
            header=None,
            names=["date", "time", "cl"],
            skiprows=skip_rows,
            comment="#",
        )

        df_obs["datetime"] = pd.to_datetime(
            df_obs["date"] + " " + df_obs["time"],
            format="%d-%m-%Y %H:%M:%S.%f",
        )
        df_obs = df_obs.set_index("datetime").sort_index()

        # Clean and process the data using a chained, readable method.
        df_obs["cl"] = df_obs["cl"].replace(-1.0000000000000001e31, np.nan)
        full_index = pd.date_range(df_obs.index.min(), df_obs.index.max(), freq="1min")
        df_obs = df_obs.reindex(full_index)

        df_obs["cl_interp"] = df_obs["cl"].interpolate(limit=5, limit_direction="both")
        mask = df_obs["cl_interp"].notna() & (
            df_obs["cl_interp"].diff().abs().astype(float) < 500.0
        )
        df_obs["cl_clean"] = df_obs["cl_interp"].where(mask)
        return df_obs[["cl_clean"]]
    except OSError as e:
        print(f"Error loading {filepath}: {e}")  # noqa: T201
        return pd.DataFrame(columns=["cl_clean"])


def transform_coords(timestamp: Any, magnetometers: Any):
    """Converts from Geographic (GEO) to Solar Magnetic (SM) coords"""
    tick = Ticktock([timestamp], "UTC")
    magnetometers_sm = []

    for station, lat_geo, lon_geo in magnetometers:
        try:
            c = coord.Coords([[1.0, lat_geo, lon_geo]], "GEO", "sph")
            c.ticks = tick
            lat_sm, lon_sm = c.convert("SM", "sph").data[0][1:3]
            magnetometers_sm.append((station, lat_sm, lon_sm))

        except (AttributeError, ValueError, TypeError) as err:
            print(err)  # noqa: T201
    return magnetometers_sm


def _get_cl_model(timestamp: Any, ds: Any, magnetometers: Any):
    """Calculates the model CL value for a single timestamp"""
    delbt = ds["delbt"]
    values = []

    for _station, lat_sm, lon_sm in transform_coords(timestamp, magnetometers):
        try:
            # Find the interpolated model data at the station's location.
            val = delbt.interp(lats=lat_sm, longs=lon_sm).item()
            values.append(val)
        except (KeyError, TypeError, ValueError):
            values.append(np.nan)
    return np.nanmin(values) if values else np.nan


def _get_cpcp_model(ds: Any):
    """Calculates the model CPCP value for a single timestamp"""
    cpcp = (
        ds.sel(lats=slice(90.0, 0.0)).pot.max()
        - ds.sel(lats=slice(90.0, 0.0)).pot.min()
    )
    return cpcp if cpcp else np.nan


def load_data_cl_model(dir_model: Any, start_time_ggcm: Any, factor: Any):
    """Loads and processes a time series of model-generated CL data"""
    pattern = re.compile(r"\.iof\.(\d{6})$")
    data_points = []

    for fname in sorted(os.listdir(dir_model)):
        match = pattern.search(fname)
        if not match:
            continue

        try:
            seconds = int(match.group(1))
            timestamp = start_time_ggcm + timedelta(seconds=seconds)
            ds = xr.open_dataset(dir_model / fname)

            if "delbt" in ds:
                cl_model = _get_cl_model(timestamp, ds, magnetometers_geo)
                if not np.isnan(cl_model):
                    # Apply the unit conversion.
                    converted_val = cl_model * factor
                    data_points.append((timestamp, converted_val))
        except OSError as e:
            print(f"Skipping file {fname}: {e}")  # noqa: T201

    if not data_points:
        return pd.DataFrame(columns=["cl_model"])

    # Create a dataframe from the collected data points.
    return pd.DataFrame(data_points, columns=["datetime", "cl_model"]).set_index(
        "datetime"
    )


def load_data_cpcp_model(dir_model: Any, start_time_ggcm: Any):
    """Loads and processes a time series of model-generated CPCP data"""
    pattern = re.compile(r"\.iof\.(\d{6})$")
    data_points = []

    for fname in sorted(os.listdir(dir_model)):
        match = pattern.search(fname)
        if not match:
            continue

        try:
            seconds = int(match.group(1))
            timestamp = start_time_ggcm + timedelta(seconds=seconds)
            ds = xr.open_dataset(dir_model / fname)

            if "pot" in ds:
                cpcp_model = _get_cpcp_model(ds)
                if not np.isnan(cpcp_model):
                    # Convert from V to kV.
                    converted_val = cpcp_model / 1e3
                    data_points.append((timestamp, converted_val))
        except OSError as e:
            print(f"Skipping file {fname}: {e}")  # noqa: T201

    if not data_points:
        return pd.DataFrame(columns=["cpcp_model"])

    # Create a dataframe from the collected data points.
    return pd.DataFrame(data_points, columns=["datetime", "cpcp_model"]).set_index("datetime")


def plot(factor: Any):
    """Loads all data and generates a plot comparing observations and model"""
    start_time_ggcm = datetime(1997, 1, 9, 22, 0)
    start_time_plot = datetime(1997, 1, 10)

    datasets = {
        "bz": load_data_wind(file_bz, "bz"),
        "n": load_data_wind(file_rr, "n"),
        "cl": load_data_cl_observed(file_cl),
        "cl_model": load_data_cl_model(dir_model, start_time_ggcm, factor),
        "cpcp_model": load_data_cpcp_model(dir_model, start_time_ggcm),
    }

    # Calculate hours from the start time for all dataframes.
    for df in datasets.values():
        if not df.empty:
            df["hours"] = (df.index - start_time_plot).total_seconds() / 3600

    fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True)
    fig.suptitle(
        f"cl_index_{start_time_plot.date()!s}_{dir_par_name}_{factor:.0e}.png",
        fontsize=16,
    )

    # Plot Panels 1 & 2 (WIND Data).
    plot_configs = [
        {
            "ax": axes[0],
            "df": datasets["bz"],
            "col": "bz",
            "label": "WIND-Bz [nT]",
            "ylim": (-20, 20),
        },
        {
            "ax": axes[1],
            "df": datasets["n"],
            "col": "n",
            "label": "WIND-N [cm$^{-3}$]",
            "ylim": (0, 20),
        },
    ]

    for config in plot_configs:
        if not config["df"].empty:
            config["ax"].plot(
                config["df"]["hours"],
                config["df"][config["col"]],
                color="black",
                lw=1,
            )
            config["ax"].set_ylabel(config["label"])
            config["ax"].set_ylim(config["ylim"])
            config["ax"].grid(True, linestyle=":", alpha=0.7)
            config["ax"].tick_params(direction="in", top=True, right=True)
    axes[0].axhline(0, color="gray", linestyle="--", linewidth=0.8)

    # Plot Panel 3 (Observed CL & Model CL).
    ax3 = axes[2]
    if not datasets["cl"].empty:
        ax3.plot(
            datasets["cl"]["hours"],
            datasets["cl"]["cl_clean"],
            label="Observed CL",
            color="black",
            lw=1.5,
        )
    if not datasets["cl_model"].empty:
        ax3.plot(
            datasets["cl_model"]["hours"],
            datasets["cl_model"]["cl_model"],
            label="Model CL",
            color="red",
            lw=1,
        )
    ax3.axhline(0, color="gray", linestyle="--", linewidth=0.8)

    # Dynamically set Y-axis limits based on available data.
    # combined_cl = pd.concat(
    #     [datasets["cl"]["cl_clean"], datasets["cl_model"]["cl_model"]]
    # ).dropna()
    # if not combined_cl.empty:
    #     ax3.set_ylim(combined_cl.min() - 50, combined_cl.max() + 50)
    # else:
    #    ax3.set_ylim(-2000, 0)

    # Manually set Y-axis limits.
    ax3.set_ylim(-2000, 0)

    ax3.set_ylabel("CL [nT]")
    ax3.set_xlim(0, 14)
    ax3.set_xticks(range(15))
    ax3.grid(True, linestyle=":", alpha=0.7)
    ax3.tick_params(direction="in", top=True, right=True)
    if not datasets["cl"].empty or not datasets["cl_model"].empty:
        ax3.legend()

    # Plot Panel 4 (Model CPCP).
    ax4 = axes[3]
    if not datasets["cpcp_model"].empty:
        ax4.plot(
            datasets["cpcp_model"]["hours"],
            datasets["cpcp_model"]["cpcp_model"],
            label="Model CPCP",
            color="black",
            lw=1,
        )
    ax4.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax4.set_ylim(0, 700)
    ax4.set_ylabel("CPCP [kV]")
    ax4.set_xlabel(f"Time (hours from {start_time_plot.strftime('%Y-%m-%d %H:%M')} UT)")
    ax4.set_xlim(0, 14)
    ax4.set_xticks(range(15))
    ax4.grid(True, linestyle=":", alpha=0.7)
    ax4.tick_params(direction="in", top=True, right=True)
    if not datasets["cpcp_model"].empty:
        ax4.legend()

    filename = f"cl_{start_time_plot.date()!s}_{dir_par_name}_{factor:.0e}.png"
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(filename)
    # plt.show()


if __name__ == "__main__":
    for num in UNIT_CONVERSION:
        plot(num)
