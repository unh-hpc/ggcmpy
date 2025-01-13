from __future__ import annotations

import numpy as np

sample_iof = {
    "coords": {"lats", "longs", "time"},
    "data_vars": [
        "prec_e_fe_1",
        "prec_e_e0_1",
        "prec_e_fe_2",
        "prec_e_e0_2",
        "pacurr",
        "fac_tot",
        "pot",
        "delphi",
        "ppio",
        "rrio",
        "ttio",
        "vdown",
        "sigh",
        "sigp",
        "fac_dyn",
        "ctiot",
        "ctiop",
        "tau",
        "ctaup",
        "ctaut",
        "cpolp",
        "cpolt",
        "delbr",
        "delbp",
        "delbt",
        "epio",
        "etio",
        "xjh",
    ],
    "time": np.asarray(
        ["1967-01-01T00:00:30.150", "1967-01-01T00:01:00.101"], dtype=np.datetime64
    ),
    "sizes": {"lats": 181, "longs": 61},
}

sample_3df = {
    "coords": {"x", "y", "z", "time"},
    "data_vars": [
        "rr",
        "pp",
        "vx",
        "vy",
        "vz",
        "bx",
        "by",
        "bz",
        "xjx",
        "xjy",
        "xjz",
        "xtra1",
        "xtra2",
        "resis",
    ],
    "time": np.asarray(["1967-01-01T00:20:00.033000"], dtype=np.datetime64),
    "sizes": {"x": 64, "y": 32, "z": 32},
}

sample_py = {
    "coords": {"x", "y", "z", "time"},
    "data_vars": [
        "rr",
        "pp",
        "vx",
        "vy",
        "vz",
        "bx",
        "by",
        "bz",
        "xjx",
        "xjy",
        "xjz",
        "xtra1",
        "xtra2",
        "resis",
    ],
    "time": np.asarray(["1967-01-01T00:20:00.033000"], dtype=np.datetime64),
    "sizes": {"x": 64, "z": 32},
}
