# mostly taken from viscid... thanks Kris

from collections import OrderedDict
import numpy as np

from .fortran_file import FortranFile
from ggcmpy import _jrrle  # type: ignore[attr-defined]

read_ascii = False


class JrrleFile(FortranFile):
    """Interface for actually opening / reading a jrrle file"""

    fields_seen = None
    seen_all_fields = None

    def __init__(self, filename):
        self._read_func = [
            _jrrle.read_jrrle1d,
            _jrrle.read_jrrle2d,
            _jrrle.read_jrrle3d,
        ]

        self.fields_seen = OrderedDict()
        self.seen_all_fields = False
        super(JrrleFile, self).__init__(filename)

    def read_field(self, fld_name, ndim):
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

    def inquire_all_fields(self, reinquire=False):
        if reinquire:
            self.seen_all_fields = False
            self.fields_seen = OrderedDict()

        if self.seen_all_fields:
            return

        self.rewind()
        while not self.seen_all_fields:
            self.inquire_next()
            # last_seen, meta = self.inquire_next()
            # if meta is not None:
            #     print(last_seen, "lives at", meta["file_position"])
            self.advance_one_line()

    def inquire(self, fld_name):
        try:
            meta = self.fields_seen[fld_name]
            self.rewind()  # FIXME
            self.seek(meta["file_position"])
            return meta
        except KeyError:
            try:
                last_added = next(reversed(self.fields_seen))
                self.seek(self.fields_seen[last_added]["file_position"])
                self.advance_one_line()
            except StopIteration:
                pass  # we haven't read any fields yet, that's ok

            while not self.seen_all_fields:
                found_fld_name, meta = self.inquire_next()
                if found_fld_name == fld_name:
                    return meta
                self.advance_one_line()

            raise KeyError(
                "file '{0}' has no field '{1}'" "".format(self.filename, fld_name)
            )

    def inquire_next(self):
        """Collect the meta-data from the next field in the file

        Returns:
            tuple (field name, dict of meta data) both
            of which will be None if there are no more Fields

        Note:
            After this operation is done, the file-pointer will be reset
            to the position it was before the inquiry.
        """
        if not self.isopen:
            raise RuntimeError("file is not open")

        varname = np.array(" " * 80, dtype="S80")
        tstring = np.array(" " * 80, dtype="S80")
        found_field, ndim, nx, ny, nz, it = _jrrle.inquire_next(
            self._unit, varname, tstring
        )
        varname = str(np.char.decode(varname)).strip()
        tstring = str(np.char.decode(tstring)).strip()

        if not found_field:
            self.seen_all_fields = True
            return None, None

        if varname in self.fields_seen:
            meta = self.fields_seen[varname]
        else:
            meta = dict(
                timestr=tstring,
                inttime=it,
                ndim=ndim,
                shape=tuple(x for x in (nx, ny, nz) if x > 0),
                file_position=self.tell(),
            )
            self.fields_seen[varname] = meta

        return varname, meta
