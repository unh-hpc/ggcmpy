import os
from typing import Any


def parse_filename(filename: str | os.PathLike[Any]) -> Any:
    dirname = os.path.dirname(filename)
    run, ifx, step = os.path.basename(filename).split(".")
    meta = {"dirname": dirname, "run": run, "step": int(step)}
    if ifx.startswith(("px_", "py_", "pz_")):
        meta["type"] = "2df"
        meta["plane"] = ifx[1]
        # given as integer in tenths of RE
        meta["plane_location"] = float(ifx[3:]) / 10
    else:
        meta["type"] = ifx
    return meta
