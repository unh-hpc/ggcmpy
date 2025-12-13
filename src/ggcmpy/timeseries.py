from __future__ import annotations

import pathlib
from collections.abc import Hashable
from typing import Any
from warnings import warn

import pandas as pd
import xarray as xr
from numpy.typing import ArrayLike

SATELLITES = {
    "om1997": "OMNI 1997",
    "wi": "Wind",
}

QUANTITIES: dict[str, dict[Hashable, str | float]] = {
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
            assert isinstance(quantity["factor"], float)
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


def write_ggcm_solarwind_file(filename: pathlib.Path, field: xr.DataArray):
    with filename.open("w") as f:
        for v in field:
            st = v.time.dt.strftime("%Y %m %d %H %M %S.%f").item()
            f.write("%s %f\n" % (st, v))


def write_ggcm_solarwind_files(sw_data: xr.Dataset, opt: Any):
    for v in _GGCM_SOLARWIND_VARIABLES:
        if v not in sw_data:
            continue

        filename = pathlib.Path(f"{opt.sat}.{v}")
        if opt.debug:
            print("Writing:%s" % filename)  # noqa: T201
        write_ggcm_solarwind_file(filename, sw_data[v])


def store_to_pyspedas(data: pd.DataFrame | xr.DataArray | xr.Dataset):
    """Stores a pandas DataFrame or xarray DataArray into pyspedas tplot variable.

    Parameters
    ----------
    data : pd.DataFrame | xr.DataArray
        DataFrame or DataArray to be stored.
    """

    if isinstance(data, pd.DataFrame):
        for varname in data.columns:
            _store_to_pyspedas(varname, data.index, data[varname], data[varname].attrs)
    elif isinstance(data, xr.Dataset):
        for key in data.data_vars:
            var = data[key]
            _store_to_pyspedas(key, var.coords["time"], var, var.attrs)
    elif isinstance(data, xr.DataArray):
        _store_to_pyspedas(data.name, data.coords["time"], data, data.attrs)
    else:
        msg = "Data must be a pandas DataFrame or xarray DataArray"  # type: ignore[unreachable]
        raise TypeError(msg)


def _store_to_pyspedas(
    varname: Hashable, x: ArrayLike, y: ArrayLike, attrs: dict[Hashable, Any]
):
    try:
        # pylint: disable=C0415
        import pyspedas  # type: ignore[import-not-found]
    except ImportError as err:
        msg = "pyspedas is required for this function. Please install it first."
        raise ImportError(msg) from err

    pyspedas.store_data(varname, data={"x": x, "y": y})

    ytitle = attrs.get("long_name", attrs.get("name"))
    if ytitle is not None:
        pyspedas.options(varname, "ytitle", ytitle)

    if "units" in attrs:
        pyspedas.options(varname, "ysubtitle", f"[{attrs['units']}]")

    if "name" in attrs:
        pyspedas.options(varname, "legend_names", [attrs["name"]])


_GGCM_SOLARWIND_VARIABLES = [
    "xgse",
    "ygse",
    "zgse",
    "bxgse",
    "bygse",
    "bzgse",
    "vxgse",
    "vygse",
    "vzgse",
    "pp",
    "rr",
    "np",
    "temp",
    "vth",
    "tkev",
    "tev",
    "btot",
    "vtot",
]
