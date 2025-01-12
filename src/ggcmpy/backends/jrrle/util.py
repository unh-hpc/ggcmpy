from __future__ import annotations

import os
from pathlib import Path
from typing import Any


def parse_filename(filename: str | os.PathLike[Any]) -> Any:
    filename = Path(filename)
    dirname = filename.parent
    run, ifx, step = filename.name.split(".")
    meta = {"dirname": dirname, "run": run, "step": int(step)}
    if ifx.startswith(("px_", "py_", "pz_")):
        meta["type"] = "2df"
        meta["plane"] = ifx[1]
        # given as integer in tenths of RE
        meta["plane_location"] = float(ifx[3:]) / 10
    else:
        meta["type"] = ifx
    return meta
