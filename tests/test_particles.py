from __future__ import annotations

import pandas as pd

from ggcmpy.tracing.boris_push import particles_cxx_from_pandas, particles_cxx_to_pandas


def test_particles_cxx():
    prts_df = pd.DataFrame(
        {
            "time": [0.0, 1.0, 2.0],
            "x": [0.0, 1.0, 2.0],
            "y": [0.0, 1.0, 2.0],
            "z": [0.0, 1.0, 2.0],
            "ux": [1.0, 0.0, 0.0],
            "uy": [0.0, 1.0, 0.0],
            "uz": [0.0, 0.0, 1.0],
            "id": [0, 1, 2],
        }
    )
    prts = particles_cxx_from_pandas(prts_df)
    assert len(prts_df) == 3
    assert repr(prts) == "particles(size=3)"
    prts_df2 = particles_cxx_to_pandas(prts)
    assert prts_df.equals(prts_df2)
