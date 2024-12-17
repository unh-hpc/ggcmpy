import datetime as dt
import numpy as np
from itertools import islice


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
