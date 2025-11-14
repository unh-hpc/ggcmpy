from __future__ import annotations

import pathlib
from warnings import warn

import pandas as pd

SATELLITES = {
    "om1997": "OMNI 1997",
    "wi": "Wind",
}

QUANTITIES: dict[str, dict[str, str | float]] = {
    "vxgse": {
        "name": "Vx (GSE)",
        "long_name": "Vx in GSE Coordinates",
        "units": "km/s",
    },
    "vygse": {
        "name": "Vy (GSE)",
        "long_name": "Vy in GSE Coordinates",
        "units": "km/s",
    },
    "vzgse": {
        "name": "Vz (GSE)",
        "long_name": "Vz in GSE Coordinates",
        "units": "km/s",
    },
    "bxgse": {"name": "Bx (GSE)", "long_name": "Bx in GSE Coordinates", "units": "nT"},
    "bygse": {"name": "By (GSE)", "long_name": "By in GSE Coordinates", "units": "nT"},
    "bzgse": {"name": "Bz (GSE)", "long_name": "Bz in GSE Coordinates", "units": "nT"},
    "rr": {"name": "N", "long_name": "Number Density", "units": "cm^-3"},
    "pp": {"name": "P", "long_name": "Pressure", "units": "nPa"},
    "tt": {"name": "T", "long_name": "Temperature", "factor": 1e-3, "units": "keV"},
}


def read_ggcm_solarwind_file(
    filename: pathlib.Path, varname: str | None = None
) -> pd.DataFrame:
    """Reads a single field time series (e.g., "wi.bzgse") file into a pandas DataFrame.

    Parameters
    ----------
    filename : pathlib.Path
        Path to the time series file.
    varname : str | None, optional
        Name of the variable. If None, the filename is used.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the time series data.
    """
    if varname is None:
        varname = filename.name
    data = pd.read_csv(
        filename,
        sep=r"\s+",
        names=["year", "month", "day", "hour", "minute", "second", varname],
    )
    if data["year"].max() < 100:
        data["year"] += 1900
    data["datetime"] = pd.to_datetime(
        data[["year", "month", "day", "hour", "minute", "second"]]
    )
    data = data.set_index("datetime")

    data = data[[varname]]  # keep only the variable column
    if data.index.has_duplicates:
        warn(
            f"duplicate times found in {filename}, dropping duplicates",
            stacklevel=2,
        )
        # Remove duplicate times, keeping the first occurrence
        data = data[~data.index.duplicated(keep="first")]
    sat, sfx = varname.split(".")  # e.g., om1997.vxgse
    sat = SATELLITES.get(sat, sat)
    if sfx in QUANTITIES:
        quantity = QUANTITIES[sfx]
        if "factor" in quantity:
            data[varname] *= quantity["factor"]
        data[varname].attrs["long_name"] = f"{sat} {quantity['long_name']}"
        data[varname].attrs["name"] = quantity["name"]
        data[varname].attrs["units"] = quantity["units"]

    return data


def read_ggcm_solarwind_directory(directory: pathlib.Path, glob: str = "*"):
    """Reads all single field time series files in a GGCM solar wind directory into a pandas DataFrame.

    Parameters
    ----------
    directory : pathlib.Path
        Path to the directory containing time series files.
    glob : str, optional
        Glob pattern to match files, by default "*".

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the concatenated time series data from all files.
    """

    dfs, attrs = [], {}
    for file in directory.glob(glob):
        df_read = read_ggcm_solarwind_file(file)
        varname = df_read.columns[0]
        dfs.append(df_read)
        attrs[varname] = df_read[varname].attrs

    df_combined = pd.concat(dfs, axis=1)
    for varname in df_combined.columns:
        df_combined[varname].attrs = attrs[varname]

    return df_combined


def store_to_pyspedas(df: pd.DataFrame):
    """Stores a pandas DataFrame into pyspedas tplot variable.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to be stored.
    """

    try:
        # pylint: disable=C0415
        import pyspedas  # type: ignore[import-not-found]
    except ImportError as err:
        msg = "pyspedas is required for this function. Please install it first."
        raise ImportError(msg) from err

    for varname in df.columns:
        pyspedas.store_data(varname, data={"x": df.index, "y": df[varname]})
        da = pyspedas.get_data(varname, xarray=True)
        attrs = df[varname].attrs
        if "long_name" in attrs:
            da.attrs["plot_options"]["yaxis_opt"]["axis_label"] = attrs["name"]
        if "name" in attrs:
            da.attrs["plot_options"]["yaxis_opt"]["legend_names"] = [attrs["name"]]
        if "units" in attrs:
            da.attrs["plot_options"]["yaxis_opt"]["axis_subtitle"] = attrs["units"]
            # pyspedas.set_units(varname, attrs["units"]) # FIXME, doesn't seem to work

        if varname.endswith(".pp"):
            pyspedas.options(varname, "ylog", True)
