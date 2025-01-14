from __future__ import annotations

import os
from collections.abc import Iterable
from typing import Any

from typing_extensions import override
from xarray.backends import BackendEntrypoint
from xarray.backends.common import AbstractDataStore
from xarray.backends.store import StoreBackendEntrypoint
from xarray.core.dataset import Dataset
from xarray.core.datatree import DataTree
from xarray.core.types import ReadBuffer

from .backends import jrrle
from .jrrle_store import JrrleStore


class JrrleEntrypoint(BackendEntrypoint):
    """Entrypoint that lets xarray recognize and read OpenGGCM jrrle (custom binary) output."""

    # url = "https://link_to/your_backend/documentation"  # FIXME

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
    def open_dataset(
        self,
        filename_or_obj: str | os.PathLike[Any] | ReadBuffer[Any] | AbstractDataStore,
        *,
        mask_and_scale: bool = True,
        decode_times: bool = True,
        concat_characters: bool = True,
        decode_coords: bool = True,
        drop_variables: str | Iterable[str] | None = None,
        use_cftime: bool | None = None,
        decode_timedelta: bool | None = None,
    ) -> Dataset:
        if isinstance(filename_or_obj, str | os.PathLike):
            store = JrrleStore.open(filename_or_obj)
        else:
            msg = f"unknown {filename_or_obj=}"
            raise TypeError(msg)

        store_entrypoint = StoreBackendEntrypoint()

        return store_entrypoint.open_dataset(
            store,
            mask_and_scale=mask_and_scale,
            decode_times=decode_times,
            concat_characters=concat_characters,
            decode_coords=decode_coords,
            drop_variables=drop_variables,
            use_cftime=use_cftime,
            decode_timedelta=decode_timedelta,
        )

    @override
    def open_datatree(
        self,
        filename_or_obj: str | os.PathLike[Any] | ReadBuffer[Any] | AbstractDataStore,
        **kwargs: Any,
    ) -> DataTree:
        raise NotImplementedError()
