import numpy as np
import pandas as pd
import pytest
from pipelines.channel_attrs.transfer import weighted_transfer


def test_weighted_transfer_length_weights():
    # COMID 1 maps to foreign reaches A (3000 m) and B (1000 m).
    xwalk = pd.DataFrame({
        "COMID": [1, 1, 2],
        "foreign_id": ["A", "B", "C"],
        "part_len": [3000.0, 1000.0, 500.0],
    })
    attrs = pd.DataFrame({"foreign_id": ["A", "B", "C"], "width": [100.0, 20.0, 7.0]})
    out = weighted_transfer(xwalk, attrs, value_col="width")
    # (3000*100 + 1000*20) / 4000 = 80
    assert out.loc[out.COMID == 1, "width"].item() == 80.0
    assert out.loc[out.COMID == 2, "width"].item() == 7.0


def test_weighted_transfer_drops_unmatched_foreign():
    xwalk = pd.DataFrame({"COMID": [1], "foreign_id": ["Z"], "part_len": [100.0]})
    attrs = pd.DataFrame({"foreign_id": ["A"], "width": [5.0]})
    out = weighted_transfer(xwalk, attrs, value_col="width")
    assert out.loc[out.COMID == 1, "width"].isna().all()
