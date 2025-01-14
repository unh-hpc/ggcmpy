from __future__ import annotations

import os
import pathlib
from typing import Any, Protocol

import numpy as np
from xarray.backends import CachingFileManager, FileManager
from xarray.backends.common import AbstractDataStore
from xarray.backends.locks import SerializableLock, ensure_lock
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset

from ggcmpy import openggcm

from .backends import jrrle

# not sure this is needed
JRRLE_LOCK = SerializableLock()


class Lock(Protocol):
    """Provides duck typing for xarray locks, which do not inherit from a common base class."""

    def acquire(self, blocking: bool = True) -> bool: ...
    def release(self) -> None: ...
    def __enter__(self) -> None: ...
    def __exit__(self, *args: Any) -> None: ...
    def locked(self) -> bool: ...


class JrrleStore(AbstractDataStore):
    def __init__(
        self,
        manager: FileManager,
        filename: str | os.PathLike[Any],
        mode: str | None = None,
        lock: Lock = JRRLE_LOCK,
        autoclose: bool = False,
    ):
        assert isinstance(manager, FileManager)
        self._manager = manager
        self._mode = mode
        self.lock = ensure_lock(lock)  # type: ignore[no-untyped-call]
        self.autoclose = autoclose
        self._filename = filename

    @classmethod
    def open(
        cls,
        filename: str | os.PathLike[Any],
        mode: str = "r",
        lock: Lock | None = None,
        autoclose: bool = False,
    ) -> JrrleStore:
        if lock is None:
            if mode == "r":
                lock = JRRLE_LOCK
            else:
                raise NotImplementedError()

        assert isinstance(filename, str | os.PathLike)

        manager = CachingFileManager(jrrle.JrrleFile, filename, mode=mode)
        return cls(manager, filename, mode=mode, lock=lock, autoclose=autoclose)

    def acquire(self, needs_lock: bool = True) -> jrrle.JrrleFile:
        with self._manager.acquire_context(needs_lock) as file:  # type: ignore[no-untyped-call]
            ds = file
        assert isinstance(ds, jrrle.JrrleFile)
        return ds

    @property
    def ds(self) -> jrrle.JrrleFile:
        return self.acquire()

    def _init_coords_dims(self, meta: dict[str, Any]) -> None:
        self._coords = dict[str, Any]()
        if meta["type"] in {"2df", "3df"}:
            grid2_filename = pathlib.Path(meta["dirname"] / f"{meta['run']}.grid2")
            self._coords = openggcm.read_grid2(grid2_filename)

        if meta["type"] == "2df":
            if meta["plane"] == "x":
                self._data_dims = ["y", "z"]
                self._coords["x"] = [meta["plane_location"]]
            elif meta["plane"] == "y":
                self._data_dims = ["x", "z"]
                self._coords["y"] = [meta["plane_location"]]
            elif meta["plane"] == "z":
                self._data_dims = ["x", "y"]
                self._coords["z"] = [meta["plane_location"]]
        elif meta["type"] == "3df":
            self._data_dims = ["x", "y", "z"]
        elif meta["type"] == "iof":
            self._data_dims = ["longs", "lats"]
        else:
            msg = f"unknown type {meta['type']}"
            raise RuntimeError(msg)

    def open_dataset(self) -> Dataset:
        meta = jrrle.parse_filename(self._filename)
        self._init_coords_dims(meta)

        shape: tuple[int, ...] | None = None
        time: str | None = None
        inttime: int | None = None
        elapsed_time: float | None = None

        with self.acquire() as f:
            variables = {}
            for fld, fld_info in f.vars.items():
                _, arr = f.read_field(fld)

                if shape is not None:
                    assert shape == arr.shape, "inconsistent shapes in jrrle file"
                if time is not None:
                    assert (
                        time == fld_info["time"]
                    ), "inconsistent time info in jrrle file"

                shape = arr.shape
                time = fld_info["time"]
                inttime = fld_info["inttime"]
                elapsed_time = fld_info["elapsed_time"]

                variables[fld] = DataArray(data=arr, dims=self._data_dims)

        assert time is not None
        assert shape is not None
        if meta["type"] == "iof":
            self._coords = {
                "lats": np.linspace(90.0, -90.0, shape[1]),
                "longs": np.linspace(-180.0, 180.0, shape[0]),
            }

        coords = self._coords
        coords["time"] = [np.datetime64(time, "ns")]
        coords["inttime"] = inttime
        coords["elapsed_time"] = elapsed_time

        attrs = {"run": meta["run"]}

        return Dataset(variables, coords=coords, attrs=attrs)
