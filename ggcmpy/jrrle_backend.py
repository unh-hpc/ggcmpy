import xarray as xr
from xarray.backends import BackendEntrypoint
import numpy as np
from itertools import islice
import os
import re

from ggcmpy.ggcm_jrrle import JrrleFileWrapper


class JrrleEntrypoint(BackendEntrypoint):
    def open_dataset(
        self,
        filename_or_obj,
        *,
        drop_variables=None,
        # other backend specific keyword arguments
        # `chunks` and `cache` DO NOT go here, they are handled by xarray
    ):
        return jrrle_open_dataset(filename_or_obj, drop_variables=drop_variables)

    open_dataset_parameters = ["filename_or_obj", "drop_variables"]

    def guess_can_open(self, filename_or_obj):
        try:
            meta = _jrrle_parse_filename(filename_or_obj)
        except:
            return False

        return True

    description = "Use OpenGGCM jrrle files in Xarray"

    url = "https://link_to/your_backend/documentation"  # FIXME


def _jrrle_parse_filename(filename):
    dirname = os.path.dirname(filename)
    run, ifx, step = os.path.basename(filename).split(".")
    meta = dict(dirname=dirname, run=run, ifx=ifx, step=int(step))
    if ifx.startswith("px_") or ifx.startswith("py_") or ifx.startswith("pz_"):
        meta["type"] = "2df"
        meta["plane"] = ifx[1]
        # given as integer in tenths of RE
        meta["plane_location"] = float(ifx[3:]) / 10
    elif ifx == "3df":
        meta["type"] = "3df"
    else:
        raise ValueError(f"Parse Error: unknown ifx {ifx}")
    return meta


def _jrrle_read_grid2(filename):
    meta = _jrrle_parse_filename(filename)
    filename = os.path.join(meta['dirname'], f"{meta['run']}.grid2")
    # load the cell centered grid
    with open(filename, 'r') as fin:
        nx = int(next(fin).split()[0])
        gx = list(islice(fin, 0, nx, 1))

        ny = int(next(fin).split()[0])
        gy = list(islice(fin, 0, ny, 1))

        nz = int(next(fin).split()[0])
        gz = list(islice(fin, 0, nz, 1))

    gx = np.array(gx, dtype='f4')
    gy = np.array(gy, dtype='f4')
    gz = np.array(gz, dtype='f4')

    return {"x": gx, "y": gy, "z": gz}


def as_isotime(time):
    """Try to convert times in string format to ISO 8601

    Raises:
        TypeError: Elements are not strings
        ValueError: numpy.datetime64(time) fails
    """
    if isinstance(time, (list, tuple, np.ndarray)):
        scalar = False
    else:
        scalar = True
        time = [time]

    ret = [None] * len(time)
    for i, t in enumerate(time):
        t = t.strip().upper().lstrip('UT')
        if re.match(r"^[0-9]{2}([0-9]{2}:){3,5}[0-9]{1,2}(\.[0-9]*)?$", t):
            # Handle YYYY:MM:DD:hh:mm:ss.ms -> YYYY-MM-DDThh:mm:ss.ms
            #        YYYY:MM:DD:hh:mm:s.ms  -> YYYY-MM-DDThh:mm:s.ms
            #        YYYY:MM:DD:hh:mm:ss    -> YYYY-MM-DDThh:mm:ss
            #        YYYY:MM:DD:hh:mm       -> YYYY-MM-DDThh:mm
            #        YYYY:MM:DD:hh          -> YYYY-MM-DDThh
            # -- all this _tsp nonsense is to take care of s.ms; annoying
            _tsp = t.replace('.', ':').split(':')
            _tsp[0] = _tsp[0].zfill(4)
            _tsp[1:6] = [_s.zfill(2) for _s in _tsp[1:6]]
            t = ":".join(_tsp[:6])
            if len(_tsp) > 6:
                t += "." + _tsp[6]
            # --
            ret[i] = t[:10].replace(':', '-') + 'T' + t[11:]

        elif re.match(r"^[0-9]{2}([0-9]{2}:){2}[0-9]{2}$", t):
            # Handle YYYY:MM:DD -> YYYY-MM-DD
            ret[i] = t.replace(':', '-')
        else:
            ret[i] = t

        try:
            np.datetime64(ret[i])
        except ValueError:
            raise

    if scalar:
        return ret[0]
    else:
        if isinstance(time, np.ndarray):
            return np.array(time, dtype=time.dtype)
        else:
            return ret


def _openggcm_parse_timestring(timestr):
    prefix = 'time='
    timestr.strip()
    if not timestr.startswith(prefix):
        raise ValueError("Time string '{0}' is malformed".format(timestr))
    timestr = timestr[len(prefix):].split()
    t = float(timestr[0])
    t = np.timedelta64(int(1000*t), "ms")
    uttime = np.datetime64(as_isotime(timestr[2]))
    return t, uttime


def jrrle_open_dataset(filename_or_obj, *, drop_variables=None):
    meta = _jrrle_parse_filename(filename_or_obj)

    coords = _jrrle_read_grid2(filename_or_obj)
    dims = coords["x"].shape[0], coords["y"].shape[0], coords["z"].shape[0]
    attrs = dict(run=meta["run"], dims=dims)

    if meta["type"] == "2df":
        if meta["plane"] == "x":
            data_dims = ["y", "z"]
            coords["x"] = [meta["plane_location"]]
        elif meta["plane"] == "y":
            data_dims = ["x", "z"]
            coords["y"] = [meta["plane_location"]]
        elif meta["plane"] == "z":
            data_dims = ["x", "y"]
            coords["z"] = [meta["plane_location"]]
    elif meta["type"] == "3df":
        data_dims = ["x", "y", "z"]

    file_wrapper = JrrleFileWrapper(filename_or_obj)
    file_wrapper.open()
    file_wrapper.inquire_all_fields()

    with file_wrapper as f:
        flds = f.fields_seen
        vars = {}
        for fld in flds.keys():
            ndim = flds[fld]["ndim"]
            fld_info, arr = f.read_field(fld, ndim)
            time, uttime = _openggcm_parse_timestring(fld_info["timestr"])
            data_attrs = dict(
                inttime=fld_info["inttime"], time=time, uttime=uttime)

            vars[fld] = xr.DataArray(
                data=arr,
                dims=data_dims,
                attrs=data_attrs)
        # vars, attrs, coords = my_decode_variables(
        #     vars, attrs, decode_times, decode_timedelta, decode_coords
        # )  #  see also conventions.decode_cf_variables

    ds = xr.Dataset(vars, coords=coords, attrs=attrs)
#    ds.set_close(my_close_method)

    return ds
