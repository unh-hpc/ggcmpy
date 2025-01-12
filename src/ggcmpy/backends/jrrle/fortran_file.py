# mostly taken from viscid... thanks Kris and Matt
# mostly stolen from pyggcm... thanks Matt

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

    _unit = -1

    filename = None
    debug = None

    def __init__(self, name: str, debug: int = 0):
        self.filename = name
        self.debug = debug

    # Make sure we close it when we're done
    def __del__(self) -> None:
        self.close()

    def open(self) -> None:
        if self.isopen:
            msg = f"Fortran file '{self.filename}' already open"
            raise RuntimeError(msg)

        with fortfile_open_lock:
            self._unit = _jrrle.fopen(self.filename, uu=-1, debug=self.debug)

        if self._unit < 0:
            msg = f"Fortran open error ({self._unit}) on '{self.filename}'"
            raise RuntimeError(msg)

    def close(self) -> None:
        if self.isopen:
            _jrrle.fclose(self._unit, debug=self.debug)
            self._unit = -1

    def seek(self, offset: int, whence: int = 0) -> int:
        assert self.isopen
        status = _jrrle.seek(self._unit, offset, whence)
        if status != 0:
            msg = f"status != 0: {status}"
            raise AssertionError(msg)
        return status  # type: ignore[no-any-return]

    def tell(self) -> int:
        assert self.isopen
        pos = _jrrle.tell(self._unit)
        assert pos >= 0
        return pos  # type: ignore[no-any-return]

    @property
    def isopen(self) -> bool:
        if self._unit > 0:
            if bool(_jrrle.fisopen(self._unit)):
                return True
            msg = "File has a valid unit, but fortran says " "it's closed?"
            raise RuntimeError(msg)
        return False

    @property
    def unit(self) -> int:
        return self._unit

    def rewind(self) -> None:
        _jrrle.frewind(self._unit, debug=self.debug)

    def advance_one_line(self) -> int:
        return _jrrle.fadvance_one_line(self._unit, debug=self.debug)  # type: ignore[no-any-return]

    def backspace(self) -> None:
        _jrrle.fbackspace(self._unit, debug=self.debug)

    def __enter__(self) -> Self:
        if not self.isopen:
            self.open()
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        if self.isopen:
            self.close()
