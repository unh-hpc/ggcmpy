from __future__ import annotations

import numpy as np
import pandas as pd

import ggcmpy
from ggcmpy.timeseries import read_ggcm_solarwind_file


def test_read_ggcm_solarwind_file():
    rr = read_ggcm_solarwind_file(
        ggcmpy.sample_dir / "cir07_19970227_liang_norcm/tmp.minvar/om1997.rr"
    )
    assert set(rr.columns) == {"om1997.rr"}
    assert pd.api.types.is_datetime64_ns_dtype(rr.index)
    assert np.isclose(rr["om1997.rr"].mean(), 4.8191107)
