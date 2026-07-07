import numpy as np
import pandas as pd
import netCDF4
from pipelines.channel_attrs.assemble import assemble


def test_assemble_matches_global_schema(tmp_path):
    global_comids = np.array([10, 20, 30, 40], dtype="int64")
    frames = {
        "corridor_impervious": pd.DataFrame({"COMID": [10, 30], "corridor_impervious": [0.9, 0.1]}),
        "bfi": pd.DataFrame({"COMID": [10, 20, 30], "bfi": [0.5, 0.6, 0.7]}),
    }
    out = tmp_path / "test.nc"
    assemble(global_comids, frames, out)
    ds = netCDF4.Dataset(out)
    assert len(ds.dimensions["COMID"]) == 4
    assert str(ds["COMID"].dtype) == "int64"
    assert str(ds["corridor_impervious"].dtype) == "float64"
    v = ds["corridor_impervious"][:]
    assert v[0] == 0.9 and np.isnan(v[1]) and v[2] == 0.1 and np.isnan(v[3])
