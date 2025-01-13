from __future__ import annotations

import os
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from typing_extensions import override
from xarray.backends import BackendEntrypoint
from xarray.backends.common import AbstractDataStore
from xarray.core.datatree import DataTree
from xarray.core.types import ReadBuffer

from . import openggcm
from .backends import jrrle
from .jrrle_store import JrrleStore


class JrrleEntrypoint(BackendEntrypoint):
    """Entrypoint that lets xarray recognize and read OpenGGCM jrrle (custom binary) output."""

    # url = "https://link_to/your_backend/documentation"  # FIXME

    def open_dataset(
        self,
        filename_or_obj,
        *,
        drop_variables=None,
        # other backend specific keyword arguments
        # `chunks` and `cache` DO NOT go here, they are handled by xarray
    ):
        return jrrle_open_dataset(filename_or_obj, drop_variables=drop_variables)

    open_dataset_parameters = ("filename_or_obj", "drop_variables")

    def guess_can_open(self, filename_or_obj):
        if not isinstance(filename_or_obj, str | os.PathLike):
            return False

        try:
            jrrle.parse_filename(filename_or_obj)
        except ValueError:
            return False

        return True

    @override
    def open_datatree(
        self,
        filename_or_obj: str | os.PathLike[Any] | ReadBuffer[Any] | AbstractDataStore,
        **kwargs: Any,
    ) -> DataTree:
        raise NotImplementedError()


def jrrle_open_dataset(
    filename_or_obj: str,
    *,
    drop_variables: Iterable[str] | None = None,  # pylint: disable=W0613  # noqa: ARG001
):
    meta = jrrle.parse_filename(filename_or_obj)

    coords: dict[str, Any]
    if meta["type"] in {"2df", "3df"}:
        grid2_filename = Path(meta["dirname"] / f"{meta['run']}.grid2")
        coords = openggcm.read_grid2(grid2_filename)
        shape = coords["x"].shape[0], coords["y"].shape[0], coords["z"].shape[0]
    else:
        shape = None

    if meta["type"] == "2df":
        if meta["plane"] == "x":
            data_dims = ["y", "z"]
            coords["x"] = [meta["plane_location"]]
        elif meta["plane"] == "y":
            data_dims = ["x", "z"]
            coords["y"] = [meta["plane_location"]]
        elif meta["plane"] == "z":
            data_dims = ["x", "y"]
            coords["z"] = [meta["plane_location"]]
    elif meta["type"] == "3df":
        data_dims = ["x", "y", "z"]
    elif meta["type"] == "iof":
        data_dims = ["longs", "lats"]
    else:
        msg = f"unknown type {type}"
        raise RuntimeError(msg)

    time: None | str = None

    store = JrrleStore.open(filename_or_obj)
    with store.acquire() as f:
        f.inquire_all_fields()

        flds = f.fields_seen
        variables = {}
        for fld in flds:
            ndim = len(flds[fld]["shape"])
            fld_info, arr = f.read_field(fld, ndim)
            if shape is None:
                shape = fld_info["shape"]
            parsed = openggcm.parse_timestring(fld_info["timestr"])
            # FIXME, should all be variables
            data_attrs = {
                "inttime": fld_info["inttime"],
                "elapsed_time": parsed["elapsed_time"],
            }

            if time is not None:
                assert time == parsed["time"], "inconsistent time info in jrrle file"
            time = parsed["time"]

            variables[fld] = xr.DataArray(data=arr, dims=data_dims, attrs=data_attrs)  # pylint: disable=E0606

    assert time is not None
    assert shape
    if meta["type"] == "iof":
        coords = {
            "lats": np.linspace(90.0, -90.0, shape[1]),
            "longs": np.linspace(-180.0, 180.0, shape[0]),
        }

    coords["time"] = [np.datetime64(time, "ns")]

    attrs = {"run": meta["run"], "shape": shape}

    return xr.Dataset(variables, coords=coords, attrs=attrs)
    #    ds.set_close(my_close_method)
