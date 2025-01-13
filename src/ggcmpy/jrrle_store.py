from __future__ import annotations

import os
from typing import Any, Protocol

from xarray.backends import CachingFileManager, FileManager
from xarray.backends.common import AbstractDataStore
from xarray.backends.locks import SerializableLock, ensure_lock

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
        mode: str | None = None,
        lock: Lock = JRRLE_LOCK,
        autoclose: bool = False,
    ):
        assert isinstance(manager, FileManager)
        self._manager = manager
        self._mode = mode
        self.lock = ensure_lock(lock)  # type: ignore[no-untyped-call]
        self.autoclose = autoclose

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
        return cls(manager, mode=mode, lock=lock, autoclose=autoclose)

    def acquire(self, needs_lock: bool = True) -> jrrle.JrrleFile:
        with self._manager.acquire_context(needs_lock) as file:  # type: ignore[no-untyped-call]
            ds = file
        assert isinstance(ds, jrrle.JrrleFile)
        return ds

    @property
    def ds(self) -> jrrle.JrrleFile:
        return self.acquire()
