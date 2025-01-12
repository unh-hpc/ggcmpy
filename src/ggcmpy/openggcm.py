from __future__ import annotations

import datetime as dt
import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, DTypeLike
from itertools import islice

from collections.abc import Sequence, Iterable


def read_grid2(filename):
    # load the cell centered grid
    with open(filename, "r") as fin:
        nx = int(next(fin).split()[0])
        gx = list(islice(fin, 0, nx, 1))

        ny = int(next(fin).split()[0])
        gy = list(islice(fin, 0, ny, 1))

        nz = int(next(fin).split()[0])
        gz = list(islice(fin, 0, nz, 1))

    gx = np.array(gx, dtype="f4")
    gy = np.array(gy, dtype="f4")
    gz = np.array(gz, dtype="f4")

    return {"x": gx, "y": gy, "z": gz}


def parse_timestring(timestr):
    if not timestr.startswith("time="):
        raise ValueError("Time string '{0}' is malformed".format(timestr))
    timestr = timestr[len("time=") :].split()

    return dict(
        elapsed_time=float(timestr[0]),
        time=np.datetime64(dt.datetime.strptime(timestr[2], "%Y:%m:%d:%H:%M:%S.%f")),
    )


def _time_array_to_dt64(times: Iterable[Sequence[int]]) -> Sequence[np.datetime64]:
    return [
        np.datetime64(
            dt.datetime(
                year=time[0],
                month=time[1],
                day=time[2],
                hour=time[3],
                minute=time[4],
                second=time[5],
                microsecond=time[6] * 1000,
            ),
            "ns",
        )
        for time in times
    ]


def _dt64_to_time_array(times: ArrayLike, dtype: DTypeLike) -> ArrayLike:
    dt_times = pd.to_datetime(times)  # type: ignore[arg-type]
    return np.array(
        [
            dt_times.year,
            dt_times.month,
            dt_times.day,
            dt_times.hour,
            dt_times.minute,
            dt_times.second,
            dt_times.microsecond // 1000,
        ],
        dtype=dtype,
    ).T
