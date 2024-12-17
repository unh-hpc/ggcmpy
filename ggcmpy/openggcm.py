import re
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


def _as_isotime(time):
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
        t = t.strip().upper().lstrip("UT")
        if re.match(r"^[0-9]{2}([0-9]{2}:){3,5}[0-9]{1,2}(\.[0-9]*)?$", t):
            # Handle YYYY:MM:DD:hh:mm:ss.ms -> YYYY-MM-DDThh:mm:ss.ms
            #        YYYY:MM:DD:hh:mm:s.ms  -> YYYY-MM-DDThh:mm:s.ms
            #        YYYY:MM:DD:hh:mm:ss    -> YYYY-MM-DDThh:mm:ss
            #        YYYY:MM:DD:hh:mm       -> YYYY-MM-DDThh:mm
            #        YYYY:MM:DD:hh          -> YYYY-MM-DDThh
            # -- all this _tsp nonsense is to take care of s.ms; annoying
            _tsp = t.replace(".", ":").split(":")
            _tsp[0] = _tsp[0].zfill(4)
            _tsp[1:6] = [_s.zfill(2) for _s in _tsp[1:6]]
            t = ":".join(_tsp[:6])
            if len(_tsp) > 6:
                t += "." + _tsp[6]
            # --
            ret[i] = t[:10].replace(":", "-") + "T" + t[11:]

        elif re.match(r"^[0-9]{2}([0-9]{2}:){2}[0-9]{2}$", t):
            # Handle YYYY:MM:DD -> YYYY-MM-DD
            ret[i] = t.replace(":", "-")
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


def parse_timestring(timestr):
    prefix = "time="
    timestr.strip()
    if not timestr.startswith(prefix):
        raise ValueError("Time string '{0}' is malformed".format(timestr))
    timestr = timestr[len(prefix) :].split()
    t = float(timestr[0])
    t = np.timedelta64(int(1000 * t), "ms")
    uttime = np.datetime64(_as_isotime(timestr[2]))
    return t, uttime
