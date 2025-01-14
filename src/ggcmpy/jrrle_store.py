from __future__ import annotations

import os
import pathlib
from collections.abc import Mapping
from typing import Any, Protocol

import numpy as np
from xarray.backends import CachingFileManager, FileManager
from xarray.backends.common import AbstractDataStore
from xarray.backends.locks import SerializableLock, ensure_lock
from xarray.core.dataset import Dataset
from xarray.core.variable import Variable

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

        self._meta = jrrle.parse_filename(self._filename)

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

    def dims(self) -> tuple[str, ...]:
        meta = self._meta
        if meta["type"] == "3df":
            return ("x", "y", "z")
        if meta["type"] == "2df":
            return tuple(dim for dim in ("x", "y", "z") if dim != meta["plane"])
        assert meta["type"] == "iof"
        return ("longs", "lats")

    def coords(self, shape: tuple[int, ...]) -> dict[str, Any]:
        meta = self._meta
        if meta["type"] in {"2df", "3df"}:
            grid2_filename = pathlib.Path(meta["dirname"] / f"{meta['run']}.grid2")
            coords: dict[str, Any] = openggcm.read_grid2(grid2_filename)
            if meta["type"] == "2df":
                coords[meta["plane"]] = [meta["plane_location"]]
            return coords

        assert meta["type"] == "iof"
        return {
            "lats": np.linspace(90.0, -90.0, shape[1]),
            "longs": np.linspace(-180.0, 180.0, shape[0]),
        }

    def open_store_variable(
        self,
        name: str,
        fld_info: Mapping[str, Any],  # noqa: ARG002
        dims: tuple[str, ...],
    ) -> Variable:
        _, arr = self.ds.read_field(name)

        return Variable(dims=dims, data=arr)

    def open_dataset(self) -> Dataset:
        shape: tuple[int, ...] | None = None
        time: str | None = None
        inttime: int | None = None
        elapsed_time: float | None = None

        dims = self.dims()

        with self.acquire() as f:
            variables = dict[str, Any]()
            for fld, fld_info in f.vars.items():
                if shape is not None:
                    assert (
                        shape == fld_info["shape"]
                    ), "inconsistent shapes in jrrle file"
                if time is not None:
                    assert (
                        time == fld_info["time"]
                    ), "inconsistent time info in jrrle file"

                shape = fld_info["shape"]
                time = fld_info["time"]
                inttime = fld_info["inttime"]
                elapsed_time = fld_info["elapsed_time"]
                variables[fld] = self.open_store_variable(fld, fld_info, dims)

        assert shape is not None
        coords = self.coords(shape)

        coords["time"] = [np.datetime64(time, "ns")]
        variables["inttime"] = inttime
        variables["elapsed_time"] = elapsed_time
        variables.update(coords)

        attrs = {"run": self._meta["run"]}

        return Dataset(variables, attrs=attrs)
