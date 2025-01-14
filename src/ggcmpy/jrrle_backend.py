from __future__ import annotations

import os
from collections.abc import Iterable
from typing import Any

from typing_extensions import override
from xarray.backends import BackendEntrypoint
from xarray.backends.common import AbstractDataStore
from xarray.core.datatree import DataTree
from xarray.core.types import ReadBuffer

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
    store = JrrleStore.open(filename_or_obj)

    meta = jrrle.parse_filename(filename_or_obj)

    return store.open_dataset(meta)
    #    ds.set_close(my_close_method)
