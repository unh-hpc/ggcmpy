# mostly taken from viscid... thanks Kris and Matt
# mostly stolen from pyggcm... thanks Matt
from __future__ import annotations

import os
from threading import Lock

from typing_extensions import Any, Self

from ggcmpy import _jrrle  # type: ignore[attr-defined]

# this lock is to prevent multiple threads from grabbing the same
# fortran file unit since checking for available units and opening
# the file may not be atomic. this should not be an issue for multiple
# processes
fortfile_open_lock = Lock()


class FortranFile:
    """
    small wrapper to allow the manipulation of fortran files from python
    """

    _unit: None | int = None

    def __init__(self, name: str | os.PathLike[Any], debug: int = 0):
        self.filename = os.fspath(name)
        self.debug = debug
        self._open()

    # Make sure we close it when we're done
    def __del__(self) -> None:
        self.close()

    def _open(self) -> None:
        if self.isopen:
            msg = f"Fortran file '{self.filename}' already open"
            raise RuntimeError(msg)

        with fortfile_open_lock:
            unit: int = _jrrle.fopen(self.filename, uu=-1, debug=self.debug)
            if unit < 0:
                msg = f"Fortran open error ({unit}) on '{self.filename}'"
                raise RuntimeError(msg)
            self._unit = unit

    def close(self) -> None:
        if self.isopen:
            _jrrle.fclose(self._unit, debug=self.debug)
            self._unit = None

    def seek(self, offset: int, whence: int = 0) -> int:
        status = _jrrle.seek(self.unit, offset, whence)
        if status != 0:
            msg = f"status != 0: {status}"
            raise AssertionError(msg)
        return status  # type: ignore[no-any-return]

    def tell(self) -> int:
        pos = _jrrle.tell(self.unit)
        assert pos >= 0
        return pos  # type: ignore[no-any-return]

    @property
    def isopen(self) -> bool:
        if self._unit is None:
            return False

        if bool(_jrrle.fisopen(self._unit)):
            return True

        msg = "File has a valid unit, but fortran says it's closed?"
        raise RuntimeError(msg)

    @property
    def unit(self) -> int:
        assert self.isopen
        assert self._unit is not None  # implied by the above, but again for typing
        return self._unit

    def rewind(self) -> None:
        _jrrle.frewind(self.unit, debug=self.debug)

    def advance_one_line(self) -> int:
        return _jrrle.fadvance_one_line(self.unit, debug=self.debug)  # type: ignore[no-any-return]

    def backspace(self) -> None:
        _jrrle.fbackspace(self.unit, debug=self.debug)

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.close()
