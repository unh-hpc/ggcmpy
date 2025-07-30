import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
import os
import re
from spacepy import coordinates as coord
from spacepy.time import Ticktock

dir_cur = os.getcwd()
dir_par = os.path.dirname(dir_cur)
dir_par_name = os.path.basename(dir_par)

dir_input = "../inp/"
dir_cl_observed = "./"
dir_cl_model = "../target/"

file_bz = os.path.join(dir_input, "wi.bzgse")
file_rr = os.path.join(dir_input, "wi.rr")
file_cl = os.path.join(dir_cl_observed, "CN_K0_MARI_2629956.txt")

# Geographic coordinates (lat, lon) for magnetometer stations
magnetometers_cl = [
    ("BACK", 64.31, -96.02),
    ("CONT", 65.75, -111.25),
    ("DAWS", 64.05, -139.40),
    ("ESKI", 61.09, -94.06),
    ("FCHU", 58.76, -94.09),
    ("MCMU", 56.65, -111.40),
    ("FSIM", 61.76, -121.26),
    ("FSMI", 60.02, -111.93),
    ("GILL", 56.38, -94.65),
    ("ISLL", 53.88, -94.69),
    ("PINA", 50.21, -96.04),
    ("RABB", 58.22, -103.68),
    ("RANK", 62.82, -92.11),
]
# Unit conversion factor for the model output
UNIT_CONVERSION = [1, 1e9, 1e10]


def load_data_wind(filepath, value_col):
    """Loads and parses time-series data from WIND files"""
    try:
        df = pd.read_csv(
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
        time = pd.to_datetime(df[["year", "month", "day", "hour", "minute"]])
        df["time"] = time + pd.to_timedelta(df["second_float"], unit="s")
        return df.set_index("time")[[value_col]]
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return pd.DataFrame(columns=[value_col])


def load_data_cl_observed(filepath):
    """Loads, cleans, and interpolates CL data"""
    try:
        with open(filepath, "r") as f:
            for i, line in enumerate(f):
                if "dd-mm-yyyy" in line:
                    skip_rows = i + 1
                    break
            else:
                raise ValueError("Header line 'dd-mm-yyyy' not found.")

        df = pd.read_csv(
            filepath,
            sep=r"\s+",
            header=None,
            names=["date", "time", "cl"],
            skiprows=skip_rows,
            comment="#",
        )

        df["datetime"] = pd.to_datetime(
            df["date"] + " " + df["time"], format="%d-%m-%Y %H:%M:%S.%f"
        )
        df = df.set_index("datetime").sort_index()

        # Clean and process the data using a chained, readable method.
        df["cl"] = df["cl"].replace(-1.0000000000000001e31, np.nan)
        full_index = pd.date_range(df.index.min(), df.index.max(), freq="1min")
        df = df.reindex(full_index)

        df["cl_interp"] = df["cl"].interpolate(limit=5, limit_direction="both")
        mask = df["cl_interp"].notna() & (df["cl_interp"].diff().abs() < 500)
        df["cl_clean"] = df["cl_interp"].where(mask)
        return df[["cl_clean"]]
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return pd.DataFrame(columns=["cl_clean"])


def _get_cl_model(timestamp, ds, magnetometers):
    """Calculate the model CL value for a single timestamp."""
    delbt = ds["delbt"]
    tick = Ticktock([timestamp], "UTC")
    values = []

    for _, geo_lat, geo_lon in magnetometers:
        try:
            # Convert from Geographic (GEO) to Solar Magnetic (SM) coords.
            c = coord.Coords([[geo_lat, geo_lon, 110.0]], "GEO", "sph")
            c.ticks = tick
            sm_lat, sm_lon = c.convert("SM", "sph").data[0][:2]
            # Find the interpolated model data at the station's location.
            val = delbt.interp(longs=sm_lon, lats=sm_lat).item()
            values.append(val)
        except Exception:
            values.append(np.nan)
    return np.nanmin(values) if values else np.nan


def load_data_cl_model(dir_model, start_time, factor):
    """Loads and processes a time series of model-generated CL data"""
    pattern = re.compile(r"\.iof\.(\d{6})$")
    data_points = []

    for fname in sorted(os.listdir(dir_model)):
        match = pattern.search(fname)
        if not match:
            continue

        try:
            seconds = int(match.group(1))
            timestamp = start_time + timedelta(seconds=seconds)
            ds = xr.open_dataset(os.path.join(dir_model, fname))

            if "delbt" in ds:
                cl_model = _get_cl_model(timestamp, ds, magnetometers_cl)
                if not np.isnan(cl_model):
                    # Apply unit conversion and store the valid data point.
                    converted_val = cl_model * factor
                    data_points.append((timestamp, converted_val))
        except Exception as e:
            print(f"Skipping file {fname}: {e}")

    if not data_points:
        return pd.DataFrame(columns=["cl_model"])

    # Create DataFrame from the collected data points.
    return pd.DataFrame(
        data_points, columns=["datetime", "cl_model"]
    ).set_index("datetime")


def plot(factor):
    """Loads all data and generates a plot comparing observations and model."""
    start_time = datetime(1997, 1, 10)

    datasets = {
        "bz": load_data_wind(file_bz, "bz"),
        "n": load_data_wind(file_rr, "n"),
        "cl": load_data_cl_observed(file_cl),
        "cl_model": load_data_cl_model(dir_cl_model, start_time, factor),
    }

    # Calculate hours from the start time for all dataframes.
    for df in datasets.values():
        if not df.empty:
            df["hours"] = (df.index - start_time).total_seconds() / 3600

    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    fig.suptitle(
        "Geomagnetic Data for January 10, 1997",
        fontsize=16,
    )

    # Plotting Panels 1 & 2 (WIND Data)
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

    # Plotting Panel 3 (Observed vs. Model CL)
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
    combined_cl = pd.concat(
        [datasets["cl"]["cl_clean"], datasets["cl_model"]["cl_model"]]
    ).dropna()
    if not combined_cl.empty:
        ax3.set_ylim(combined_cl.min() - 50, combined_cl.max() + 50)
    else:
        ax3.set_ylim(-2000, 0)
    
    # Manually set Y-axis limits.
    # ax3.set_ylim(-2000, 0)

    ax3.set_ylabel("CL [nT]")
    ax3.set_xlabel(
        f"Time (hours from {start_time.strftime('%Y-%m-%d %H:%M')} UT)"
    )
    ax3.set_xlim(0, 14)
    ax3.set_xticks(range(0, 15))
    ax3.grid(True, linestyle=":", alpha=0.7)
    ax3.tick_params(direction="in", top=True, right=True)
    if not datasets["cl"].empty or not datasets["cl_model"].empty:
        ax3.legend()

    filename = f"cl_{str(start_time.date())}_{dir_par_name}_{factor:.0e}.png"
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(filename)
    # plt.show()


if __name__ == "__main__":
    for factor in UNIT_CONVERSION:
        plot(factor)
