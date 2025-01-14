# mostly taken from viscid... thanks Kris
from __future__ import annotations

import logging
import os
from collections import OrderedDict
from collections.abc import Mapping
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ggcmpy import _jrrle, openggcm  # type: ignore[attr-defined]

from .fortran_file import FortranFile

logger = logging.getLogger(__name__)

read_ascii = False

_jrrle_read_func = (
    _jrrle.read_jrrle1d,
    _jrrle.read_jrrle2d,
    _jrrle.read_jrrle3d,
)


def _jrrle_read_field(file, fld_name: str, meta: dict[str, Any]) -> NDArray[Any]:
    shape = meta["shape"]
    ndim = len(shape)
    arr = np.empty(shape, dtype="float32", order="F")
    success = _jrrle_read_func[ndim - 1](file.unit, arr, fld_name, read_ascii)
    if not success:
        msg = "read_func failed"
        raise RuntimeError(msg)
    return arr


def _jrrle_inquire_next(
    file: FortranFile,
) -> tuple[str, dict[str, Any]] | None:
    b_varname = np.array(" " * 80, dtype="S80")
    b_tstring = np.array(" " * 80, dtype="S80")
    found_field, ndim, nx, ny, nz, it = _jrrle.inquire_next(
        file.unit, b_varname, b_tstring
    )
    if not found_field:
        return None

    shape = (nx, ny, nz)[:ndim]
    varname = str(np.char.decode(b_varname)).strip()
    timestr = str(np.char.decode(b_tstring)).strip()
    parsed = openggcm.parse_timestring(timestr)

    return varname, {
        "shape": shape,
        "inttime": it,
        "time": parsed["time"],
        "elapsed_time": parsed["elapsed_time"],
        "file_position": file.tell(),
    }


class JrrleFile(FortranFile):
    """Interface for actually opening / reading a jrrle file"""

    def __init__(self, filename: str | os.PathLike[Any], mode: str = "r"):
        assert mode == "r"

        self.fields_seen: OrderedDict[str, Any] = OrderedDict()
        self.seen_all_fields = False
        super().__init__(filename)

    def read_field(self, fld_name) -> tuple[Any, NDArray[Any]]:
        """Read a field"""
        logger.debug("read_field(%s)" % fld_name)
        meta = self._inquire(fld_name)
        return meta, _jrrle_read_field(self, fld_name, meta)

    @property
    def vars(self) -> Mapping[str, Mapping[str, Any]]:
        self._inquire_all_fields()
        return self.fields_seen

    def _inquire_all_fields(self) -> None:
        if self.seen_all_fields:
            return

        logger.debug("_inquire_all_fields() (%s)" % self.filename)
        self.rewind()
        while True:
            rv = self._inquire_next()
            if not rv:
                break
            self.advance_one_line()

    def _inquire(self, fld_name: str) -> Any:
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

        while True:
            rv = self._inquire_next()
            if not rv:
                break
            found_fld_name, meta = rv
            if found_fld_name == fld_name:
                return meta
            self.advance_one_line()

        msg = f"file '{self.filename}' has no field '{fld_name}'"
        raise KeyError(msg)

    def _inquire_next(self) -> tuple[str, Any] | None:
        """Collect the meta-data from the next field in the file

        Returns:
            tuple (field name, dict of meta data) both
            of which will be None if there are no more Fields

        Note:
            After this operation is done, the file-pointer will be reset
            to the position it was before the inquiry.
        """

        rv = _jrrle_inquire_next(self)
        if not rv:
            # end of file
            self.seen_all_fields = True
            return None

        varname, meta = rv
        if varname in self.fields_seen:
            assert meta == self.fields_seen[varname]
        else:
            self.fields_seen[varname] = meta

        return varname, meta
