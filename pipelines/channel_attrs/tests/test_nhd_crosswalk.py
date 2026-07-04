import geopandas as gpd
import pytest
from shapely.geometry import LineString

from pipelines.channel_attrs.nhd_crosswalk import build_crosswalk


def test_crosswalk_length_weighted_matching():
    merit = gpd.GeoDataFrame(
        {"COMID": [10]},
        geometry=[LineString([(0, 0), (2000, 0)])],
        crs="EPSG:5070",
    )
    # NHD reach A runs along the first 1500 m (offset 50 m north);
    # NHD reach B is 5 km away (no match).
    nhd = gpd.GeoDataFrame(
        {"nhd_comid": [900, 901]},
        geometry=[
            LineString([(0, 50), (1500, 50)]),
            LineString([(0, 5000), (1500, 5000)]),
        ],
        crs="EPSG:5070",
    )
    xw = build_crosswalk(merit, nhd, buffer_m=300.0, top_k=40)
    assert set(xw["foreign_id"]) == {900}
    row = xw.iloc[0]
    assert row["COMID"] == 10
    assert abs(row["part_len"] - 1500.0) < 1.0  # clipped NHD length inside the buffer
    # quality: 1500 m of a 2000 m MERIT reach matched
    assert abs(row["match_frac"] - 0.75) < 0.01


def test_crosswalk_no_overlap_returns_empty():
    merit = gpd.GeoDataFrame(
        {"COMID": [10]},
        geometry=[LineString([(0, 0), (2000, 0)])],
        crs="EPSG:5070",
    )
    nhd = gpd.GeoDataFrame(
        {"nhd_comid": [900]},
        geometry=[LineString([(0, 10000), (2000, 10000)])],
        crs="EPSG:5070",
    )
    xw = build_crosswalk(merit, nhd, buffer_m=300.0, top_k=40)
    assert xw.empty


def test_crosswalk_multiple_merit_reaches():
    """Two MERIT reaches matched by different NHD reaches."""
    merit = gpd.GeoDataFrame(
        {"COMID": [10, 11]},
        geometry=[
            LineString([(0, 0), (1000, 0)]),
            LineString([(5000, 0), (6000, 0)]),
        ],
        crs="EPSG:5070",
    )
    nhd = gpd.GeoDataFrame(
        {"nhd_comid": [900, 901]},
        geometry=[
            LineString([(0, 50), (1000, 50)]),
            LineString([(5000, 50), (6000, 50)]),
        ],
        crs="EPSG:5070",
    )
    xw = build_crosswalk(merit, nhd, buffer_m=300.0, top_k=40)
    assert set(xw["COMID"]) == {10, 11}
    assert set(xw["foreign_id"]) == {900, 901}
    # Both should be fully matched (1000 m of 1000 m reach)
    for comid in [10, 11]:
        frac = xw.loc[xw["COMID"] == comid, "match_frac"].iloc[0]
        assert abs(frac - 1.0) < 0.01
