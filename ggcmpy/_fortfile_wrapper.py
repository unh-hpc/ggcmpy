# mostly stolen from pyggcm... thanks Matt

import ctypes as ct
from collections import OrderedDict
from numpy.ctypeslib import ndpointer
import numpy as np
from threading import Lock

_libggcm = ct.CDLL('_libggcm.cpython-311-darwin.so')

_libggcm.fopen.argtypes = [
    ct.c_char_p,
    ct.c_int]
_libggcm.fopen.restype = ct.c_int

_libggcm.fclose.argtypes = [
    ct.c_int,
    ct.c_int]

_libggcm.fisopen.argtypes = [
    ct.c_int,
    ct.c_int]
_libggcm.fisopen.restype = ct.c_int

_libggcm.tell.argtypes = [
    ct.c_int]
_libggcm.tell.restype = ct.c_int

_libggcm.fadvance_one_line.argtypes = [
    ct.c_int, ct.c_int]
_libggcm.fadvance_one_line.restype = ct.c_int

_libggcm.inquire_next.argtypes = [
    ct.c_int,                 # unit
    ct.POINTER(ct.c_int),     # found_field
    ct.POINTER(ct.c_int),     # ndim
    ct.POINTER(ct.c_int),     # nx
    ct.POINTER(ct.c_int),     # ny
    ct.POINTER(ct.c_int),     # nz
    ct.POINTER(ct.c_int),     # it
    ct.POINTER(ct.c_char_p),  # varname
    ct.POINTER(ct.c_char_p),  # tstring
]

# this lock is to prevent multiple threads from grabbing the same
# fortran file unit since checking for available units and opening
# the file may not be atomic. this should not be an issue for multiple
# processes
fortfile_open_lock = Lock()


class FortranFile(object):
    """
    small wrapper to allow the manipulation of fortran files from python
    """
    _unit = -1

    filename = None
    debug = None

    def __init__(self, name, debug=0):
        self.filename = name
        self.debug = debug

    # Make sure we close it when we're done
    def __del__(self):
        self.close()

    def open(self):
        if self.isopen:
            raise RuntimeError("Fortran file '{0}' already open"
                               "".format(self.filename))

        with fortfile_open_lock:
            self._unit = _libggcm.fopen(
                self.filename.encode('utf-8'), self.debug)

        if self._unit < 0:
            raise RuntimeError("Fortran open error ({0}) on '{1}'"
                               "".format(self._unit, self.filename))

    def close(self):
        if self.isopen:
            _libggcm.fclose(self._unit, self.debug)
            self._unit = -1

    def seek(self, offset, whence=0):
        assert self.isopen
        status = _jrrle.seek(self._unit, offset, whence)
        if status != 0:
            raise AssertionError("status != 0: {0}".format(status))
        return status

    def tell(self):
        assert self.isopen
        print(f"unit {self._unit}")
        pos = _libggcm.tell(self._unit)
        print(f"pos {pos}")
        assert pos >= 0
        return pos

    @property
    def isopen(self):
        if self._unit > 0:
            if bool(_libggcm.fisopen(self._unit, self.debug)):
                return True
            else:
                raise RuntimeError("File has a valid unit, but fortran says "
                                   "it's closed?")
        return False

    @property
    def unit(self):
        return self._unit

    def rewind(self):
        _libggcm.frewind(self._unit, self.debug)

    def advance_one_line(self):
        return _libggcm.fadvance_one_line(self._unit, self.debug)

    def backspace(self):
        _jrrle.fbackspace(self._unit, debug=self.debug)

    def __enter__(self):
        if not self.isopen:
            self.open()
        return self

    def __exit__(self, exc_type, value, traceback):
        if self.isopen:
            self.close()


class JrrleFileWrapper(FortranFile):
    """Interface for actually opening / reading a jrrle file"""
    fields_seen = None
    seen_all_fields = None

    def __init__(self, filename):
        # self._read_func = [_jrrle.read_jrrle1d, _jrrle.read_jrrle2d,
        #                    _jrrle.read_jrrle3d]

        self.fields_seen = OrderedDict()
        self.seen_all_fields = False
        super(JrrleFileWrapper, self).__init__(filename)

    def inquire_all_fields(self, reinquire=False):
        if reinquire:
            self.seen_all_fields = False
            self.fields_seen = OrderedDict()

        if self.seen_all_fields:
            return

        self.rewind()
        while not self.seen_all_fields:
            last_seen, meta = self.inquire_next()
            if meta is not None:
                print(last_seen, "lives at", meta["file_position"])
            self.advance_one_line()

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

        c_found_field = ct.c_int()
        c_ndim = ct.c_int()
        c_nx = ct.c_int()
        c_ny = ct.c_int()
        c_nz = ct.c_int()
        c_it = ct.c_int()
        c_varname = ct.c_char_p()
        c_tstring = ct.c_char_p()
        # found_field, ndim, nx, ny, nz, it = _jrrle.inquire_next(self._unit,
        #                                                         varname,
        #                                                         tstring)
        _libggcm.inquire_next(self._unit, c_found_field,
                              c_ndim, c_nx, c_ny, c_nz, c_it, c_varname, c_tstring)
        found_field = bool(c_found_field.value)
        ndim, nx, ny, nz, it = [x.value for x in [
            c_ndim, c_nx, c_ny, c_nz, c_it]]

        if not found_field:
            self.seen_all_fields = True
            return None, None

        varname = c_varname.value.decode('utf-8')
        tstring = c_tstring.value.decode('utf-8')
        print(f"varname {varname} tstring {tstring}")

        if varname in self.fields_seen:
            meta = self.fields_seen[varname]
        else:
            dims = tuple(x for x in (nx, ny, nz) if x > 0)

            meta = dict(timestr=tstring,
                        inttime=it,
                        ndim=ndim,
                        dims=dims,
                        file_position=self.tell())
            self.fields_seen[varname] = meta

        return varname, meta
