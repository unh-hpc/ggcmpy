from __future__ import annotations

import asyncio
import threading
from typing import Any, cast

import numpy as np
from xarray.core import indexing

from ggcmpy.jrrle_array import JrrleArray


def test_async_getitem_runs_in_worker_thread(monkeypatch):
    class DummyStore:
        def __init__(self) -> None:
            self.lock = threading.Lock()

    array = JrrleArray("rr", cast(Any, DummyStore()), {"shape": (1,)})
    key = indexing.BasicIndexer((slice(None),))
    main_thread_id = threading.get_ident()
    worker_thread_id = None

    def fake_getitem(_, indexer):
        nonlocal worker_thread_id
        worker_thread_id = threading.get_ident()
        assert indexer == key
        return np.array([1.0], dtype=np.float32)

    monkeypatch.setattr(JrrleArray, "__getitem__", fake_getitem)

    result = asyncio.run(array.async_getitem(key))

    np.testing.assert_array_equal(result, np.array([1.0], dtype=np.float32))
    assert worker_thread_id is not None
    assert worker_thread_id != main_thread_id
