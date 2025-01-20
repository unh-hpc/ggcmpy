from __future__ import annotations

import datetime as dt
import os
from collections.abc import Hashable, Iterable, Mapping, Sequence
from itertools import islice
from pathlib import Path
from typing import Any
from warnings import warn

import numpy as np
import pandas as pd
import xarray as xr
from numpy.typing import ArrayLike, DTypeLike, NDArray
from typing_extensions import override
from xarray.coding.times import CFDatetimeCoder


def read_grid2(filename: os.PathLike[Any] | str) -> dict[str, NDArray[Any]]:
    # load the cell centered grid
    with Path(filename).open(encoding="utf-8") as fin:
        nx = int(next(fin).split()[0])
        gx = list(islice(fin, 0, nx, 1))

        ny = int(next(fin).split()[0])
        gy = list(islice(fin, 0, ny, 1))

        nz = int(next(fin).split()[0])
        gz = list(islice(fin, 0, nz, 1))

    return {
        "x": np.array(gx, dtype="f4"),
        "y": np.array(gy, dtype="f4"),
        "z": np.array(gz, dtype="f4"),
    }


def parse_timestring(timestr: str) -> dict[str, Any]:
    if not timestr.startswith("time="):
        msg = f"Time string '{timestr}' is malformed"
        raise ValueError(msg)
    time_comps = timestr.removeprefix("time=").split()

    return {
        "elapsed_time": float(time_comps[0]),
        "time": np.datetime64(
            dt.datetime.strptime(time_comps[2], "%Y:%m:%d:%H:%M:%S.%f")
        ),
    }


def decode_openggcm(ds: xr.Dataset) -> xr.Dataset:
    warn(
        "decode_openggcm() is deprecated. Pass 'decode_times=openggcm.AmieTimeArrayDecoder()' "
        "to open_dataset() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    coder = AmieTimeArrayCoder()
    for name, var in ds.variables.items():
        ds[name] = coder.decode(var, name)

    return ds


def encode_openggcm(
    variables: Mapping[str, xr.Variable], attributes: Mapping[str, Any]
) -> tuple[Mapping[str, xr.Variable], Mapping[str, Any]]:
    coder = AmieTimeArrayCoder()
    new_variables = {name: coder.encode(var, name) for name, var in variables.items()}

    return new_variables, attributes


class AmieTimeArrayCoder(CFDatetimeCoder):
    """
    A custom coder for encoding and decoding time arrays in xarray Variables.
    This class extends the CFDatetimeCoder to handle variables with a custom
    "time_array" encoding, and otherwise falls back to the usual CF coding/encoding.

    Methods
    -------
    encode(variable: xr.Variable, name: Hashable | None = None) -> xr.Variable
        Encodes an xarray Variable with a "time_array" unit into a custom format.
    decode(variable: xr.Variable, name: Hashable | None = None) -> xr.Variable
        Decodes an xarray Variable with a "time_array" unit from a custom format.
    """

    def __init__(self, use_cftime: bool | None = None):
        super().__init__(use_cftime=use_cftime)

    @override
    def encode(
        self, variable: xr.Variable, name: Hashable | None = None
    ) -> xr.Variable:
        if variable.encoding.get("units") == "time_array":
            attrs = variable.attrs.copy()
            attrs["units"] = "time_array"
            attrs.pop("_FillValue", None)

            return xr.Variable(
                dims=(*variable.dims, "time_array"),
                data=_dt64_to_time_array(
                    variable,
                    variable.encoding.get("dtype", "int32"),
                ),
                attrs=attrs,
            )

        return super().encode(variable, name)

    @override
    def decode(
        self, variable: xr.Variable, name: Hashable | None = None
    ) -> xr.Variable:
        if variable.attrs.get("units") == "time_array":
            times: Any = variable.to_numpy().tolist()
            if variable.ndim == 1:
                times = _time_array_to_dt64([times])[0]
            else:
                times = _time_array_to_dt64(times)

            dims = (dim for dim in variable.dims if dim != "time_array")
            attrs = variable.attrs.copy()
            encoding = {"units": attrs.pop("units"), "dtype": variable.dtype}
            return xr.Variable(dims=dims, data=times, attrs=attrs, encoding=encoding)

        return super().decode(variable, name)


def _time_array_to_dt64(times: Iterable[Sequence[int]]) -> Sequence[np.datetime64]:
    return [
        np.datetime64(
            dt.datetime(
                year=time[0],
                month=time[1],
                day=time[2],
                hour=time[3],
                minute=time[4],
                second=time[5],
                microsecond=time[6] * 1000,
            ),
            "ns",
        )
        for time in times
    ]


def _dt64_to_time_array(times: ArrayLike, dtype: DTypeLike) -> ArrayLike:
    dt_times = pd.to_datetime(np.asarray(times))
    return np.array(
        [
            dt_times.year,
            dt_times.month,
            dt_times.day,
            dt_times.hour,
            dt_times.minute,
            dt_times.second,
            dt_times.microsecond // 1000,
        ],
        dtype=dtype,
    ).T


@xr.register_dataarray_accessor("ggcm")  # type: ignore[no-untyped-call]
@xr.register_dataset_accessor("ggcm")  # type: ignore[no-untyped-call]
class OpenGGCMAccessor:
    """
    Xarray accessor to add OpenGGCM-specific features

    As of now, that is ds.ggcm.coords which provides mlts and colats.
    """

    def __init__(self, xarray_obj: Any):
        self._obj = xarray_obj
        self._coords: xr.Coordinates = self._obj.coords

        if "colats" not in self._coords and "lats" in self._coords:
            self._coords = self._coords.assign(colats=90 - self._coords["lats"])

        if "mlts" not in self._coords and "longs" in self._coords:
            self._coords = self._coords.assign(
                mlts=(self._coords["longs"] + 180) * 24 / 360,
            )

    def __repr__(self) -> str:
        return "OpenGGCM accessor\n" + repr(self._coords)

    @property
    def coords(self) -> xr.Coordinates:
        return self._coords

    def __getattr__(self, name: str) -> Any:
        if name in self._coords:
            return self._coords[name]

        msg = f"Attribute {name} not found"
        raise AttributeError(msg)
