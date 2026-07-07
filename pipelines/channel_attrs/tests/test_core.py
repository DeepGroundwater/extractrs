"""Tests for the pure-function API: sample_along_lines and weighted_transfer."""
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import pytest
from shapely.geometry import LineString

import pipelines.channel_attrs as ca
from pipelines.channel_attrs import sample_along_lines, weighted_transfer


def test_public_api_exports_all_symbols():
    """Every name listed in __all__ must be importable and callable."""
    import inspect
    for name in ca.__all__:
        obj = getattr(ca, name, None)
        assert obj is not None, f"__all__ lists '{name}' but it is not on the package"
        assert callable(obj), f"'{name}' is not callable"


# ── sample_along_lines ────────────────────────────────────────────────────────

def _make_raster(ny, nx, dx=100.0, dy=100.0, value=1.0, crs="EPSG:5070"):
    """Uniform-value 2-D DataArray with cell-centre coordinates."""
    x = np.arange(nx) * dx + dx / 2
    y = np.arange(ny) * dy + dy / 2
    data = np.full((ny, nx), value, dtype=np.float64)
    da = xr.DataArray(data, dims=["y", "x"], coords={"y": y, "x": x})
    da = da.rio.write_crs(crs)
    return da


def _make_gdf(lines, comids, crs="EPSG:5070"):
    return gpd.GeoDataFrame({"COMID": comids}, geometry=lines, crs=crs)


def test_sample_along_lines_uniform_raster():
    # A uniform raster should return that value for every reach.
    da = _make_raster(10, 10, value=2.5)
    gdf = _make_gdf(
        [LineString([(50, 50), (950, 50)]), LineString([(50, 150), (50, 950)])],
        [1, 2],
    )
    result = sample_along_lines(da, gdf)
    assert list(result["COMID"]) == [1, 2]
    assert np.allclose(result["value_mean"], 2.5, atol=0.01)
    assert (result["n"] > 0).all()


def test_sample_along_lines_preserves_comid_order():
    da = _make_raster(10, 10, value=1.0)
    gdf = _make_gdf(
        [LineString([(50, 50), (100, 50)]), LineString([(200, 50), (300, 50)])],
        [99, 7],
    )
    result = sample_along_lines(da, gdf)
    assert list(result["COMID"]) == [99, 7]


def test_sample_along_lines_empty_gdf():
    da = _make_raster(5, 5)
    gdf = gpd.GeoDataFrame({"COMID": []}, geometry=[], crs="EPSG:5070")
    result = sample_along_lines(da, gdf)
    assert len(result) == 0


def test_sample_along_lines_nodata_excluded():
    # A raster whose cells are all NaN yields n=0 and NaN mean.
    data = np.full((5, 5), np.nan)
    da = xr.DataArray(data, dims=["y", "x"],
                      coords={"y": np.arange(5) * 100.0 + 50, "x": np.arange(5) * 100.0 + 50})
    da = da.rio.write_crs("EPSG:5070")
    gdf = _make_gdf([LineString([(50, 50), (450, 50)])], [1])
    result = sample_along_lines(da, gdf)
    assert result["n"].iloc[0] == 0
    assert np.isnan(result["value_mean"].iloc[0])


def test_sample_along_lines_reprojects_gdf():
    # GDF in 4326 must be reprojected to match a 5070 raster before sampling.
    da = _make_raster(20, 20, dx=50000.0, dy=50000.0, value=5.0)
    # Construct a line in 4326 that overlaps the raster's projected extent.
    from pyproj import Transformer
    t = Transformer.from_crs("EPSG:5070", "EPSG:4326", always_xy=True)
    lon1, lat1 = t.transform(50000, 50000)
    lon2, lat2 = t.transform(950000, 50000)
    gdf = _make_gdf([LineString([(lon1, lat1), (lon2, lat2)])], [1], crs="EPSG:4326")
    result = sample_along_lines(da, gdf)
    assert result["n"].iloc[0] > 0
    assert np.isclose(result["value_mean"].iloc[0], 5.0, atol=0.1)


def test_sample_along_lines_bed_depths_frac_below():
    # Constant raster value 3.0 m; bed depth 2.0 m → every vertex is below bed.
    da = _make_raster(10, 10, value=3.0)
    gdf = _make_gdf([LineString([(50, 50), (950, 50)])], [1])
    result = sample_along_lines(da, gdf, bed_depths={1: 2.0})
    assert "frac_below" in result.columns
    assert np.isclose(result["frac_below"].iloc[0], 1.0)


def test_sample_along_lines_bed_depths_above_bed():
    # Raster value 1.0 m < bed depth 5.0 m → no vertices below bed.
    da = _make_raster(10, 10, value=1.0)
    gdf = _make_gdf([LineString([(50, 50), (950, 50)])], [1])
    result = sample_along_lines(da, gdf, bed_depths={1: 5.0})
    assert np.isclose(result["frac_below"].iloc[0], 0.0)


# ── weighted_transfer ─────────────────────────────────────────────────────────

def _xwalk(*rows):
    """Build a crosswalk DataFrame from (COMID, foreign_id, part_len) tuples."""
    return pd.DataFrame(rows, columns=["COMID", "foreign_id", "part_len"])


def test_weighted_transfer_equal_lengths():
    # Two NHD reaches of equal length → unweighted average.
    xw = _xwalk((1, 10, 1000.0), (1, 20, 1000.0))
    attrs = pd.DataFrame({"foreign_id": [10, 20], "val": [2.0, 4.0]})
    out = weighted_transfer(xw, attrs, "val")
    row = out.loc[out["COMID"] == 1, "val"].iloc[0]
    assert np.isclose(row, 3.0)


def test_weighted_transfer_unequal_lengths():
    # One reach contributes 3× the length → weighted mean closer to its value.
    xw = _xwalk((1, 10, 750.0), (1, 20, 250.0))
    attrs = pd.DataFrame({"foreign_id": [10, 20], "val": [4.0, 0.0]})
    out = weighted_transfer(xw, attrs, "val")
    row = out.loc[out["COMID"] == 1, "val"].iloc[0]
    assert np.isclose(row, 3.0)  # 4*750 + 0*250 / 1000


def test_weighted_transfer_all_null_attrs():
    # No matched foreign IDs → output COMID row has NaN.
    xw = _xwalk((1, 99, 1000.0))
    attrs = pd.DataFrame({"foreign_id": [10], "val": [5.0]})
    out = weighted_transfer(xw, attrs, "val")
    assert len(out) == 1
    assert np.isnan(out.loc[out["COMID"] == 1, "val"].iloc[0])


def test_weighted_transfer_skips_zero_length():
    # A part_len=0 row must be excluded from the weighted average.
    xw = _xwalk((1, 10, 0.0), (1, 20, 1000.0))
    attrs = pd.DataFrame({"foreign_id": [10, 20], "val": [999.0, 2.0]})
    out = weighted_transfer(xw, attrs, "val")
    assert np.isclose(out.loc[out["COMID"] == 1, "val"].iloc[0], 2.0)


def test_weighted_transfer_multiple_merit_comids():
    xw = _xwalk((1, 10, 500.0), (2, 10, 500.0), (2, 20, 500.0))
    attrs = pd.DataFrame({"foreign_id": [10, 20], "val": [1.0, 3.0]})
    out = weighted_transfer(xw, attrs, "val")
    assert np.isclose(out.loc[out["COMID"] == 1, "val"].iloc[0], 1.0)
    assert np.isclose(out.loc[out["COMID"] == 2, "val"].iloc[0], 2.0)


def test_weighted_transfer_output_covers_all_xwalk_comids():
    # Even a COMID with no attribute match must appear in the output (NaN).
    xw = _xwalk((1, 10, 500.0), (2, 99, 500.0))
    attrs = pd.DataFrame({"foreign_id": [10], "val": [1.0]})
    out = weighted_transfer(xw, attrs, "val")
    assert set(out["COMID"]) == {1, 2}
    assert np.isnan(out.loc[out["COMID"] == 2, "val"].iloc[0])
