from __future__ import annotations

import asyncio
import importlib.util
import threading
from pathlib import Path

import numpy as np
from xarray.core import indexing


def test_async_getitem_runs_in_worker_thread(monkeypatch):
    module_path = Path(__file__).resolve().parents[1] / "src" / "ggcmpy" / "jrrle_array.py"
    spec = importlib.util.spec_from_file_location("test_jrrle_array_module", module_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    class DummyStore:
        def __init__(self) -> None:
            self.lock = threading.Lock()

    array = module.JrrleArray("rr", DummyStore(), {"shape": (1,)})
    key = indexing.BasicIndexer((slice(None),))
    main_thread_id = threading.get_ident()
    worker_thread_id = None

    def fake_getitem(_, indexer):
        nonlocal worker_thread_id
        worker_thread_id = threading.get_ident()
        assert indexer == key
        return np.array([1.0], dtype=np.float32)

    monkeypatch.setattr(module.JrrleArray, "__getitem__", fake_getitem)

    result = asyncio.run(array.async_getitem(key))

    np.testing.assert_array_equal(result, np.array([1.0], dtype=np.float32))
    assert worker_thread_id is not None
    assert worker_thread_id != main_thread_id
