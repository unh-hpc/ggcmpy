# mostly taken from viscid... thanks Kris
from __future__ import annotations

import contextlib
import os
from collections import OrderedDict
from typing import Any

import numpy as np

from ggcmpy import _jrrle  # type: ignore[attr-defined]

from .fortran_file import FortranFile

read_ascii = False


def _jrrle_inquire_next(
    file: FortranFile,
) -> tuple[str, dict[str, Any]]:
    b_varname = np.array(" " * 80, dtype="S80")
    b_tstring = np.array(" " * 80, dtype="S80")
    found_field, ndim, nx, ny, nz, it = _jrrle.inquire_next(
        file.unit, b_varname, b_tstring
    )
    if not found_field:
        raise StopIteration
    shape = (nx, ny, nz)[:ndim]
    varname = str(np.char.decode(b_varname)).strip()
    timestr = str(np.char.decode(b_tstring)).strip()
    return varname, {
        "shape": shape,
        "inttime": it,
        "timestr": timestr,
        "file_position": file.tell(),
    }


class JrrleFile(FortranFile):
    """Interface for actually opening / reading a jrrle file"""

    def __init__(self, filename: str | os.PathLike[Any], mode: str = "r"):
        assert mode == "r"
        self._read_func = [
            _jrrle.read_jrrle1d,
            _jrrle.read_jrrle2d,
            _jrrle.read_jrrle3d,
        ]

        self.fields_seen: OrderedDict[str, Any] = OrderedDict()
        self.seen_all_fields = False
        super().__init__(filename)

    def read_field(self, fld_name, ndim) -> tuple[Any, Any]:
        """Read a field given a seekable location

        Parameters:
            loc(int): position in file we can seek to
            ndim(int): dimensionality of field

        Returns:
            tuple (field name, dict of meta data, array)
        """
        meta = self.inquire(fld_name)
        arr = np.empty(meta["shape"], dtype="float32", order="F")
        success = self._read_func[ndim - 1](self.unit, arr, fld_name, read_ascii)
        if not success:
            msg = "read_func failed"
            raise RuntimeError(msg)
        return meta, arr

    def inquire_all_fields(self) -> None:
        if self.seen_all_fields:
            return

        self.rewind()
        with contextlib.suppress(StopIteration):
            while True:
                self.inquire_next()
                self.advance_one_line()

    def inquire(self, fld_name: str) -> Any:
        if fld_name in self.fields_seen:
            meta = self.fields_seen[fld_name]
            # For some mysterious reason, the rewind is necessary -- without it,
            # reading after the seek fails. Mysteriously, a self.tell() also
            # sometimes makes things work.
            self.rewind()
            self.seek(meta["file_position"])
            return meta

        if self.fields_seen:
            # seek to past last previously read field
            last_added = next(reversed(self.fields_seen))
            self.seek(self.fields_seen[last_added]["file_position"])
            self.advance_one_line()

        try:
            while True:
                found_fld_name, meta = self.inquire_next()
                if found_fld_name == fld_name:
                    return meta
                self.advance_one_line()
        except StopIteration as err:
            msg = f"file '{self.filename}' has no field '{fld_name}'"
            raise KeyError(msg) from err

    def inquire_next(self) -> tuple[str, Any]:
        """Collect the meta-data from the next field in the file

        Returns:
            tuple (field name, dict of meta data) both
            of which will be None if there are no more Fields

        Note:
            After this operation is done, the file-pointer will be reset
            to the position it was before the inquiry.
        """

        try:
            varname, meta = _jrrle_inquire_next(self)
        except StopIteration:
            self.seen_all_fields = True
            raise

        if varname in self.fields_seen:
            assert meta == self.fields_seen[varname]
        else:
            self.fields_seen[varname] = meta

        return varname, meta
