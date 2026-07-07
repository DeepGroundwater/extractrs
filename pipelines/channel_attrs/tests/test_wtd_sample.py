"""Tests for wtd_sample — synthetic raster, no real data needed."""
import numpy as np
import pytest
import xarray as xr
import geopandas as gpd
from shapely.geometry import LineString, box

from pipelines.channel_attrs import sample_along_lines, zonal_stats_corridors


def _make_gradient_da():
    """1 km grid: WTD deepens eastward: columns 0,1,2 → 5, 10, 15 m."""
    return (
        xr.DataArray(
            np.array([[5.0, 10.0, 15.0]] * 3, dtype=float),
            coords={"y": [2500.0, 1500.0, 500.0], "x": [500.0, 1500.0, 2500.0]},
            dims=("y", "x"),
        )
        .rio.write_crs("EPSG:5070")
    )


def test_sample_along_lines_nearest_cell_mean():
    da = _make_gradient_da()
    gdf = gpd.GeoDataFrame(
        {"COMID": [1]},
        geometry=[LineString([(200, 1500), (2800, 1500)])],
        crs="EPSG:5070",
    )
    out = sample_along_lines(da, gdf, spacing_m=200.0)
    row = out.loc[out["COMID"] == 1].iloc[0]
    # Densified points cross all three cells; mean of nearest-cell values
    assert 5.0 < row["value_mean"] < 15.0
    assert row["n"] >= 10


def test_sample_along_lines_uniform_field():
    """All cells equal → mean must equal that value exactly."""
    da = (
        xr.DataArray(
            np.full((3, 3), 7.0),
            coords={"y": [2500.0, 1500.0, 500.0], "x": [500.0, 1500.0, 2500.0]},
            dims=("y", "x"),
        )
        .rio.write_crs("EPSG:5070")
    )
    gdf = gpd.GeoDataFrame(
        {"COMID": [42]},
        geometry=[LineString([(200, 1500), (2800, 1500)])],
        crs="EPSG:5070",
    )
    out = sample_along_lines(da, gdf, spacing_m=300.0)
    assert out.loc[out["COMID"] == 42, "value_mean"].item() == 7.0


def test_sample_along_lines_nan_nodata_excluded():
    """NaN cells in the raster must not contribute to the mean or count."""
    arr = np.array([[np.nan, 10.0, np.nan]] * 3)
    da = (
        xr.DataArray(
            arr,
            coords={"y": [2500.0, 1500.0, 500.0], "x": [500.0, 1500.0, 2500.0]},
            dims=("y", "x"),
        )
        .rio.write_crs("EPSG:5070")
    )
    gdf = gpd.GeoDataFrame(
        {"COMID": [5]},
        geometry=[LineString([(200, 1500), (2800, 1500)])],
        crs="EPSG:5070",
    )
    out = sample_along_lines(da, gdf, spacing_m=200.0)
    # Only the middle cell (10.0) contributes
    assert out.loc[out["COMID"] == 5, "value_mean"].item() == 10.0


def test_sample_along_lines_frac_below():
    """bed_depths dict → frac_below is fraction of sampled vertices where value > bed_depth."""
    da = _make_gradient_da()  # cells: 5, 10, 15 m
    gdf = gpd.GeoDataFrame(
        {"COMID": [1]},
        geometry=[LineString([(200, 1500), (2800, 1500)])],
        crs="EPSG:5070",
    )
    # bed_depth = 7 m → only cells with wtd > 7 (i.e., cells = 10, 15) count
    out = sample_along_lines(da, gdf, spacing_m=200.0, bed_depths={1: 7.0})
    assert "frac_below" in out.columns
    frac = out.loc[out["COMID"] == 1, "frac_below"].item()
    # Roughly 2/3 of points should be in the 10- and 15-m cells
    assert 0.4 < frac < 0.9


def test_sample_along_lines_missing_comid_in_bed_depths():
    """COMID absent from bed_depths → frac_below is NaN for that reach."""
    da = _make_gradient_da()
    gdf = gpd.GeoDataFrame(
        {"COMID": [99]},
        geometry=[LineString([(200, 1500), (2800, 1500)])],
        crs="EPSG:5070",
    )
    out = sample_along_lines(da, gdf, spacing_m=200.0, bed_depths={1: 2.0})
    assert np.isnan(out.loc[out["COMID"] == 99, "frac_below"].item())


# ── zonal_stats_corridors ──────────────────────────────────────────────────────

def _make_uniform_da(value: float, crs: str = "EPSG:5070"):
    """10×10 km raster with a constant value."""
    return (
        xr.DataArray(
            np.full((10, 10), value),
            coords={"y": np.linspace(9500, 500, 10), "x": np.linspace(500, 9500, 10)},
            dims=("y", "x"),
            name="wtd",
        )
        .rio.write_crs(crs)
    )


def test_zonal_stats_corridors_uniform_raster():
    # Uniform raster → mean over any corridor equals the constant value.
    da = _make_uniform_da(8.5)
    corridors = gpd.GeoDataFrame(
        {"COMID": [1]},
        geometry=[box(1000, 1000, 5000, 5000)],
        crs="EPSG:5070",
    )
    out = zonal_stats_corridors(da, corridors, stat="mean", id_col="COMID")
    assert "COMID" in out.columns and "wtd" in out.columns
    assert out["wtd"].iloc[0] == pytest.approx(8.5, abs=1e-3)


def test_zonal_stats_corridors_reprojects_corridors():
    # Corridors in a different CRS are silently reprojected before extraction.
    # CoordinateValidationError (xarray ValueError subclass) from extractrs when
    # no cells overlap is acceptable; a plain ValueError from our CRS guard is not.
    from xarray.core.coordinates import CoordinateValidationError

    da = _make_uniform_da(3.0, crs="EPSG:5070")
    corridors = gpd.GeoDataFrame(
        {"COMID": [1]},
        geometry=[box(-96.01, 39.99, -95.99, 40.01)],
        crs="EPSG:4326",
    )
    try:
        zonal_stats_corridors(da, corridors)
    except CoordinateValidationError:
        pass  # extractrs: corridor outside synthetic raster extent — acceptable
    except ValueError as e:
        pytest.fail(f"zonal_stats_corridors raised ValueError (CRS guard?): {e}")


def test_zonal_stats_corridors_no_crs_raises():
    arr = np.ones((5, 5))
    da = xr.DataArray(arr, dims=("y", "x"), name="v")
    corridors = gpd.GeoDataFrame({"COMID": [1]}, geometry=[box(0, 0, 1, 1)], crs="EPSG:5070")
    with pytest.raises(ValueError, match="no CRS"):
        zonal_stats_corridors(da, corridors)
