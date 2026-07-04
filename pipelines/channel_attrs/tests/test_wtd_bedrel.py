"""Tests for wtd_bedrel.bed_relative — sign convention + NaN propagation."""
import numpy as np
import pandas as pd

from pipelines.channel_attrs.wtd_bedrel import bed_relative


def test_bed_relative_sign_convention():
    """Spec-mandated sign test: positive = water table BELOW bed = losing-possible."""
    wtd = pd.DataFrame({"COMID": [1, 2], "wtd_channel_zs": [10.0, 0.5]})
    bank = pd.DataFrame({"COMID": [1, 2], "bankfull_depth": [2.0, 2.0]})
    out = bed_relative(wtd, bank)
    # 10 - 2 = 8 (water table well below bed → losing-possible)
    assert out.loc[out["COMID"] == 1, "channel_wtd_bed_rel"].item() == 8.0
    # 0.5 - 2 = -1.5 (water table above bed → gaining)
    assert out.loc[out["COMID"] == 2, "channel_wtd_bed_rel"].item() == -1.5


def test_bed_relative_nan_wtd():
    """NaN WTD propagates to NaN bed-relative head."""
    wtd = pd.DataFrame({"COMID": [1], "wtd_channel_zs": [np.nan]})
    bank = pd.DataFrame({"COMID": [1], "bankfull_depth": [2.0]})
    out = bed_relative(wtd, bank)
    assert np.isnan(out.loc[out["COMID"] == 1, "channel_wtd_bed_rel"].item())


def test_bed_relative_nan_bankfull():
    """NaN bankfull_depth propagates to NaN bed-relative head."""
    wtd = pd.DataFrame({"COMID": [1], "wtd_channel_zs": [5.0]})
    bank = pd.DataFrame({"COMID": [1], "bankfull_depth": [np.nan]})
    out = bed_relative(wtd, bank)
    assert np.isnan(out.loc[out["COMID"] == 1, "channel_wtd_bed_rel"].item())


def test_bed_relative_unmatched_comid_is_nan():
    """COMIDs in wtd but not bankfull → NaN (left join behaviour)."""
    wtd = pd.DataFrame({"COMID": [1, 99], "wtd_channel_zs": [5.0, 3.0]})
    bank = pd.DataFrame({"COMID": [1], "bankfull_depth": [2.0]})
    out = bed_relative(wtd, bank)
    assert out.loc[out["COMID"] == 1, "channel_wtd_bed_rel"].item() == 3.0
    assert np.isnan(out.loc[out["COMID"] == 99, "channel_wtd_bed_rel"].item())
