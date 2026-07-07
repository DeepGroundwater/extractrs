"""Tests for sword_width.load_sword_reaches and build_crosswalk."""
import numpy as np
import pandas as pd
import geopandas as gpd
import pytest
from shapely.geometry import Polygon

from pipelines.channel_attrs import sword_width


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

class _FakeVar:
    """Minimal netCDF4 Variable stub: slice returns self, .data returns array."""
    def __init__(self, data):
        self._arr = np.ma.array(data)

    def __getitem__(self, _):
        return self

    @property
    def data(self):
        return self._arr.data if isinstance(self._arr, np.ma.MaskedArray) else self._arr


def _fake_dataset(reach_ids, xs, ys, widths, lengths, lakeflags):
    """Build a stub that mimics netCDF4.Dataset(...).groups['reaches'][var][:]."""
    arrays = {
        "reach_id": _FakeVar(np.array(reach_ids, dtype=np.int64)),
        "x": _FakeVar(np.array(xs, dtype=float)),
        "y": _FakeVar(np.array(ys, dtype=float)),
        "width": _FakeVar(np.array(widths, dtype=float)),
        "reach_length": _FakeVar(np.array(lengths, dtype=float)),
        "lakeflag": _FakeVar(np.array(lakeflags, dtype=int)),
    }

    class _FakeGroup:
        def __getitem__(self, key):
            return arrays[key]

    class _FakeDS:
        groups = {"reaches": _FakeGroup()}
        def close(self):
            pass

    return _FakeDS()


@pytest.fixture()
def patch_sword(monkeypatch):
    """Provide a factory that replaces netCDF4.Dataset with the stub."""
    import netCDF4

    def _install(reach_ids, xs, ys, widths, lengths, lakeflags):
        ds = _fake_dataset(reach_ids, xs, ys, widths, lengths, lakeflags)
        monkeypatch.setattr(netCDF4, "Dataset", lambda _path: ds)

    return _install


# ---------------------------------------------------------------------------
# load_sword_reaches tests
# ---------------------------------------------------------------------------

def test_load_sword_reaches_keeps_conus_river_reaches(patch_sword):
    # reach 71xxx → CONUS pfaf-2 = 71; reach 11xxx → non-CONUS
    patch_sword(
        reach_ids=[71_000_000_001, 71_000_000_002, 11_000_000_003],
        xs=[-90.0, -91.0, -80.0],
        ys=[30.0, 31.0, 35.0],
        widths=[50.0, 100.0, 60.0],
        lengths=[1000.0, 2000.0, 1500.0],
        lakeflags=[0, 0, 0],
    )
    df = sword_width.load_sword_reaches()
    assert len(df) == 2
    assert set(df["foreign_id"].tolist()) == {71_000_000_001, 71_000_000_002}


def test_load_sword_reaches_excludes_nonpositive_width(patch_sword):
    # negative width is SWORD missing-data sentinel
    patch_sword(
        reach_ids=[71_000_000_001, 71_000_000_002, 71_000_000_003],
        xs=[-90.0, -91.0, -92.0],
        ys=[30.0, 31.0, 32.0],
        widths=[50.0, -1.0, 0.0],
        lengths=[1000.0, 500.0, 800.0],
        lakeflags=[0, 0, 0],
    )
    df = sword_width.load_sword_reaches()
    assert len(df) == 1
    assert df["foreign_id"].iloc[0] == 71_000_000_001


def test_load_sword_reaches_excludes_lakes(patch_sword):
    # lakeflag != 0 → lake / reservoir
    patch_sword(
        reach_ids=[71_000_000_001, 71_000_000_002],
        xs=[-90.0, -91.0],
        ys=[30.0, 31.0],
        widths=[50.0, 80.0],
        lengths=[1000.0, 2000.0],
        lakeflags=[0, 1],
    )
    df = sword_width.load_sword_reaches()
    assert len(df) == 1
    assert df["foreign_id"].iloc[0] == 71_000_000_001


# ---------------------------------------------------------------------------
# build_crosswalk tests
# ---------------------------------------------------------------------------

@pytest.fixture()
def two_catchments():
    """Two non-overlapping 1°×1° boxes as a synthetic MERIT_CAT GeoDataFrame."""
    return gpd.GeoDataFrame(
        {"COMID": [1001, 1002]},
        geometry=[
            Polygon([(-95, 30), (-94, 30), (-94, 31), (-95, 31)]),
            Polygon([(-93, 30), (-92, 30), (-92, 31), (-93, 31)]),
        ],
        crs="EPSG:4326",
    )


def test_build_crosswalk_assigns_comid_by_containment(monkeypatch, two_catchments):
    monkeypatch.setattr(gpd, "read_file", lambda _path, **kw: two_catchments)

    sword_df = pd.DataFrame({
        "foreign_id": [71_001, 71_002],
        "x": [-94.5, -92.5],   # inside box 1001 and 1002 respectively
        "y": [30.5, 30.5],
        "reach_length": [1000.0, 2000.0],
        "width": [50.0, 80.0],
    })
    xwalk = sword_width.build_crosswalk(sword_df)

    assert set(xwalk.columns) >= {"COMID", "foreign_id", "part_len"}
    row1 = xwalk[xwalk["foreign_id"] == 71_001].iloc[0]
    row2 = xwalk[xwalk["foreign_id"] == 71_002].iloc[0]
    assert int(row1["COMID"]) == 1001
    assert int(row2["COMID"]) == 1002
    assert row1["part_len"] == pytest.approx(1000.0)
    assert row2["part_len"] == pytest.approx(2000.0)


def test_build_crosswalk_comid_dtype_is_int(monkeypatch, two_catchments):
    monkeypatch.setattr(gpd, "read_file", lambda _path, **kw: two_catchments)

    sword_df = pd.DataFrame({
        "foreign_id": [71_001],
        "x": [-94.5],
        "y": [30.5],
        "reach_length": [500.0],
        "width": [40.0],
    })
    xwalk = sword_width.build_crosswalk(sword_df)
    assert xwalk["COMID"].dtype == int or np.issubdtype(xwalk["COMID"].dtype, np.integer)
