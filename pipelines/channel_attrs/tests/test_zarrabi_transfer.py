"""Tests for zarrabi_transfer.load_zarrabi."""
import textwrap
import numpy as np
import pandas as pd
import pytest

from pipelines.channel_attrs.zarrabi_transfer import load_zarrabi


def _write_zarrabi_csv(tmp_path, rows):
    """Write a minimal Bankfull_Meanflow_CONUS.txt stub."""
    header = ",COMID,REACHCODE,TotDASqKM,StreamOrde,bnk_depth,bnk_width,mf_depth,mf_width\n"
    content = header + "".join(
        f"{i},{r['COMID']},{r.get('REACHCODE','A')},{r.get('TotDASqKM',1.0)},"
        f"{r.get('StreamOrde',3)},{r['bnk_depth']},{r['bnk_width']},"
        f"{r.get('mf_depth',0.5)},{r.get('mf_width',5.0)}\n"
        for i, r in enumerate(rows)
    )
    p = tmp_path / "Bankfull_Meanflow_CONUS.txt"
    p.write_text(content)
    return p


def test_load_zarrabi_renames_columns(tmp_path):
    p = _write_zarrabi_csv(tmp_path, [
        {"COMID": 1001, "bnk_depth": 1.5, "bnk_width": 20.0},
        {"COMID": 1002, "bnk_depth": 2.0, "bnk_width": 30.0},
    ])
    df = load_zarrabi(p)
    assert "foreign_id" in df.columns
    assert "bankfull_depth" in df.columns
    assert "bankfull_width" in df.columns
    assert "COMID" not in df.columns
    assert "bnk_depth" not in df.columns
    assert "bnk_width" not in df.columns


def test_load_zarrabi_foreign_id_dtype(tmp_path):
    p = _write_zarrabi_csv(tmp_path, [{"COMID": 1001, "bnk_depth": 1.2, "bnk_width": 15.0}])
    df = load_zarrabi(p)
    assert df["foreign_id"].dtype == np.dtype("int32")


def test_load_zarrabi_values_preserved(tmp_path):
    p = _write_zarrabi_csv(tmp_path, [
        {"COMID": 9999, "bnk_depth": 3.14, "bnk_width": 42.0},
    ])
    df = load_zarrabi(p)
    assert df["foreign_id"].iloc[0] == 9999
    assert df["bankfull_depth"].iloc[0] == pytest.approx(3.14)
    assert df["bankfull_width"].iloc[0] == pytest.approx(42.0)
