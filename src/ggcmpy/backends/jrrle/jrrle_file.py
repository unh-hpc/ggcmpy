# mostly taken from viscid... thanks Kris
from __future__ import annotations

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
    return varname, {"shape": shape, "inttime": it, "timestr": timestr}


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
        self._read_func[ndim - 1](self.unit, arr, fld_name, read_ascii)
        return meta, arr

    def inquire_all_fields(self) -> None:
        if self.seen_all_fields:
            return

        self.rewind()
        for _ in self:
            self.advance_one_line()

    def __iter__(self):
        return self

    def __next__(self):
        return self.inquire_next()

    def inquire(self, fld_name: str) -> Any:
        try:
            meta = self.fields_seen[fld_name]
            self.rewind()  # FIXME
            self.seek(meta["file_position"])
            return meta
        except KeyError as key_error:
            try:
                last_added = next(reversed(self.fields_seen))
                self.seek(self.fields_seen[last_added]["file_position"])
                self.advance_one_line()
            except StopIteration:
                pass  # we haven't read any fields yet, that's ok

            for found_fld_name, meta in self:
                if found_fld_name == fld_name:
                    return meta
                self.advance_one_line()

            msg = f"file '{self.filename}' has no field '{fld_name}'"
            raise KeyError(msg) from key_error

    def inquire_next(self) -> tuple[str | None, Any]:
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

        meta["file_position"] = self.tell()

        if varname in self.fields_seen:
            assert meta == self.fields_seen[varname]
        else:
            self.fields_seen[varname] = meta

        return varname, meta
