import xarray as xr
from xarray.backends import BackendEntrypoint
import numpy as np

import os

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
            fname = os.path.basename(filename_or_obj)
        except TypeError:
            return False

        try:
            run, ifx, step = fname.split('.')
        except:
            return False

        return ifx in {"3df"}

    description = "Use OpenGGCM jrrle files in Xarray"

    url = "https://link_to/your_backend/documentation"  # FIXME


def jrrle_open_dataset(filename_or_obj, *, drop_variables=None):

    file_wrapper = JrrleFileWrapper(filename_or_obj)
    file_wrapper.open()
    file_wrapper.inquire_all_fields()
    file_wrapper.rewind()  # FIXME why?

    with file_wrapper as f:
        flds = f.fields_seen
        vars = {}
        for fld in flds.keys():
            ndim = flds[fld]["ndim"]
            meta, arr = f.read_field(fld, ndim)
            # dims = meta["dims"]
            # arr = np.array(arr.flatten(order='F').reshape(dims),  # [::-1]),
            #                order='F')

            vars[fld] = ("x", "y", "z"), arr
        # vars[fld] = xr.DataArray(
        #     data=arr,
        #     dims=["x", "y", "z"])
        # vars, attrs, coords = my_decode_variables(
        #     vars, attrs, decode_times, decode_timedelta, decode_coords
        # )  #  see also conventions.decode_cf_variables

    ds = xr.Dataset(vars)  # , attrs=attrs, coords=coords)
#    ds.set_close(my_close_method)

    return ds
