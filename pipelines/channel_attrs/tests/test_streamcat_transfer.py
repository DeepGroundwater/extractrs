"""Unit test for streamcat_transfer.load_streamcat JSON parsing helper.

The weighted_transfer core is already tested in test_transfer.py.
This test exercises only the JSON parsing (structure, column renaming, dtype).
"""
import json
import tempfile
from pathlib import Path

import pandas as pd

from pipelines.channel_attrs.streamcat_transfer import load_streamcat


def _write_region_json(tmp_dir: Path, region: str, records: list[dict]) -> None:
    p = tmp_dir / f"pctimp2019_{region}.json"
    p.write_text(json.dumps({"items": records}))


def test_load_streamcat_parses_items_and_renames():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        _write_region_json(tmp_path, "Region01", [
            {"comid": 100, "pctimp2019catrp100": 5.0, "pctimp2019cat": 3.0},
            {"comid": 200, "pctimp2019catrp100": None, "pctimp2019cat": 1.0},
        ])
        _write_region_json(tmp_path, "Region02", [
            {"comid": 300, "pctimp2019catrp100": 72.5, "pctimp2019cat": 60.0},
        ])
        df = load_streamcat(tmp_path)

    assert set(df.columns) == {"foreign_id", "corridor_impervious"}
    assert len(df) == 3
    assert df.loc[df["foreign_id"] == 100, "corridor_impervious"].item() == 5.0
    assert pd.isna(df.loc[df["foreign_id"] == 200, "corridor_impervious"].item())
    assert df.loc[df["foreign_id"] == 300, "corridor_impervious"].item() == 72.5


def test_load_streamcat_foreign_id_dtype():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        _write_region_json(tmp_path, "Region01", [{"comid": 999, "pctimp2019catrp100": 10.0}])
        df = load_streamcat(tmp_path)

    # Must be int32 to match crosswalk foreign_id dtype.
    assert df["foreign_id"].dtype == "int32"


def test_load_streamcat_multi_region_concatenates():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        for i, region in enumerate(["Region01", "Region02", "Region03N"]):
            _write_region_json(tmp_path, region, [
                {"comid": (i + 1) * 1000 + j, "pctimp2019catrp100": float(j)}
                for j in range(5)
            ])
        df = load_streamcat(tmp_path)

    assert len(df) == 15
