import xarray as xr
from xarray.backends import BackendEntrypoint
import numpy as np
import os

from .backends import jrrle
from . import openggcm


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
            jrrle.parse_filename(filename_or_obj)
        except TypeError:
            return False

        return True

    description = "Use OpenGGCM jrrle files in Xarray"

    url = "https://link_to/your_backend/documentation"  # FIXME


def jrrle_open_dataset(filename_or_obj, *, drop_variables=None):
    meta = jrrle.parse_filename(filename_or_obj)

    if meta["type"] in {"2df", "3df"}:
        grid2_filename = os.path.join(meta["dirname"], f"{meta['run']}.grid2")
        coords = openggcm.read_grid2(grid2_filename)
        shape = coords["x"].shape[0], coords["y"].shape[0], coords["z"].shape[0]
    else:
        shape = None

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
    elif meta["type"] == "iof":
        data_dims = ["longs", "lats"]

    file_wrapper = jrrle.JrrleFile(filename_or_obj)
    file_wrapper.open()
    file_wrapper.inquire_all_fields()

    time = None
    with file_wrapper as f:
        flds = f.fields_seen
        vars = {}
        for fld in flds.keys():
            ndim = flds[fld]["ndim"]
            fld_info, arr = f.read_field(fld, ndim)
            if shape is None:
                shape = fld_info["dims"]
            data_attrs = dict(inttime=fld_info["inttime"])
            data_attrs.update(openggcm.parse_timestring(fld_info["timestr"]))

            if time is not None:
                assert (
                    time == data_attrs["time"]
                ), "inconsistent time info in jrrle file"
            time = data_attrs["time"]

            vars[fld] = xr.DataArray(
                data=arr, dims=data_dims, attrs=data_attrs)
        # vars, attrs, coords = my_decode_variables(
        #     vars, attrs, decode_times, decode_timedelta, decode_coords
        # )  #  see also conventions.decode_cf_variables

    if meta["type"] == "iof":
        coords = dict(
            lats=("lats", np.linspace(90.0, -90.0, shape[1])),
            colats=("lats", np.linspace(0.0, 180.0, shape[1])),
            longs=("longs", np.linspace(-180.0, 180.0, shape[0])),
            mlts=("longs", np.linspace(0.0, 24.0, shape[0])),
        )

    attrs = dict(run=meta["run"], shape=shape, time=time)

    ds = xr.Dataset(vars, coords=coords, attrs=attrs)
    #    ds.set_close(my_close_method)

    return ds
