from __future__ import annotations

from typing import TYPE_CHECKING, Any

from numpy.typing import NDArray
from xarray.backends.common import BackendArray
from xarray.core import indexing

if TYPE_CHECKING:
    from .jrrle_store import JrrleStore


class JrrleArray(BackendArray):
    """Lazy evaluation of a variable stored in an adios2 file.

    This also takes care of slicing out the specific component of the data stored as 4-d array.
    """

    def __init__(
        self,
        variable_name: str,
        datastore: JrrleStore,
    ) -> None:
        self.variable_name = variable_name
        self.datastore = datastore
        array = self.get_array()
        self.shape = array.shape
        self.dtype = array.dtype

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
