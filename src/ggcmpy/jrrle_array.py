from __future__ import annotations

import logging
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray
from xarray.backends.common import BackendArray
from xarray.core import indexing

if TYPE_CHECKING:
    from .jrrle_store import JrrleStore

logger = logging.getLogger(__name__)


class JrrleArray(BackendArray):
    """Lazy evaluation of a variable stored in an adios2 file.

    This also takes care of slicing out the specific component of the data stored as 4-d array.
    """

    def __init__(
        self, variable_name: str, datastore: JrrleStore, fld_info: Mapping[str, Any]
    ) -> None:
        self.variable_name = variable_name
        self.datastore = datastore
        self.shape = fld_info["shape"]
        self.dtype = np.dtype(np.float32)

    def get_array(self, needs_lock: bool = True) -> NDArray[Any]:
        _, arr = self.datastore.acquire(needs_lock).read_field(self.variable_name)
        return arr

    def __getitem__(self, key: indexing.ExplicitIndexer) -> NDArray[Any]:
        return indexing.explicit_indexing_adapter(  # type: ignore[no-any-return]
            key, self.shape, indexing.IndexingSupport.BASIC, self._getitem
        )

    def _getitem(self, key) -> NDArray[Any]:
        with self.datastore.lock:
            return self.get_array(needs_lock=False)[key]  # type: ignore[no-any-return]
