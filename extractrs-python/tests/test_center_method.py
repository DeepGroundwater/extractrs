"""Tests for the method='center' coverage mode.

Validates that extractrs with method='center' produces results identical
to xarray-spatial's cell-center rasterization approach.
"""

import numpy as np
import pytest
import xarray as xr
import geopandas as gpd
from shapely.geometry import box
from shapely.affinity import rotate

import extractrs  # noqa: F401 — registers .extrs accessor
from extractrs._extractrs import build_cache, apply_stat


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_da(ny, nx):
    """Create a 2D DataArray with deterministic values."""
    rng = np.random.default_rng(42)
    y = np.arange(ny, dtype=np.float64) + 0.5
    x = np.arange(nx, dtype=np.float64) + 0.5
    data = rng.standard_normal((ny, nx)).astype(np.float64)
    return xr.DataArray(data, dims=["y", "x"], coords={"y": y, "x": x})


def make_gdf(polys, ids):
    return gpd.GeoDataFrame({"zone_id": ids}, geometry=polys)


def rasterize_for_xrspatial(gdf, da):
    """Rasterize zones using rasterio (same logic as our benchmark)."""
    from rasterio.features import rasterize as rio_rasterize
    from rasterio.transform import from_origin

    y = da["y"].values
    x = da["x"].values
    dy = abs(float(y[1] - y[0]))
    dx = abs(float(x[1] - x[0]))
    transform = from_origin(float(x.min()) - dx / 2, float(y.max()) + dy / 2, dx, dy)
    shapes = [(geom, zid) for geom, zid in zip(gdf.geometry, gdf["zone_id"])]
    zone_arr = rio_rasterize(
        shapes, out_shape=(len(y), len(x)), transform=transform,
        fill=0, dtype=np.int32, all_touched=False,
    )
    return xr.DataArray(zone_arr, dims=["y", "x"], coords={"y": y, "x": x})


# ---------------------------------------------------------------------------
# Unit tests: build_cache with method parameter
# ---------------------------------------------------------------------------

class TestBuildCacheMethod:
    """Test that build_cache accepts and dispatches on the method parameter."""

    def test_default_is_exact(self):
        """Calling without method should produce exact (fractional) coverage."""
        da = make_da(4, 4)
        # Triangle that partially covers cells
        poly = box(0.0, 0.0, 2.5, 2.5)
        gdf = make_gdf([poly], [1])
        wkb_list = [g.wkb for g in gdf.geometry]
        cache = build_cache(wkb_list, [1], 0.0, 0.0, 4.0, 4.0, 1.0, 1.0)
        result = np.asarray(apply_stat(cache, da.values, -9999.0, "count"))
        # exact: partial cells contribute fractional count
        # center: only whole cells count
        # 2.5 x 2.5 box → exact count should be 6.25 (area in cell units)
        assert result[0] == pytest.approx(6.25, abs=0.01)

    def test_center_method_binary(self):
        """method='center' should produce integer-valued counts (binary coverage)."""
        da = make_da(4, 4)
        poly = box(0.0, 0.0, 2.5, 2.5)
        gdf = make_gdf([poly], [1])
        wkb_list = [g.wkb for g in gdf.geometry]
        cache = build_cache(wkb_list, [1], 0.0, 0.0, 4.0, 4.0, 1.0, 1.0, "center")
        result = np.asarray(apply_stat(cache, da.values, -9999.0, "count"))
        # Center of cells (0.5,0.5), (0.5,1.5), (1.5,0.5), (1.5,1.5) are inside
        # (2.5,0.5), (2.5,1.5), (0.5,2.5), (1.5,2.5) have center on boundary → excluded by geo::Contains
        # So count should be an integer
        assert result[0] == round(result[0])

    def test_invalid_method_raises(self):
        """Unknown method should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown method"):
            build_cache([], [], 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, "invalid")


# ---------------------------------------------------------------------------
# Accessor tests
# ---------------------------------------------------------------------------

class TestAccessorMethod:
    """Test that the .extrs accessor passes method through correctly."""

    def test_dataarray_accessor_method_param(self):
        da = make_da(10, 10)
        gdf = make_gdf([box(1.0, 1.0, 5.0, 5.0)], [1])
        exact = da.extrs.zonal_stats(gdf, stat="mean", id_col="zone_id", method="exact")
        center = da.extrs.zonal_stats(gdf, stat="mean", id_col="zone_id", method="center")
        # Both should return a DataArray
        assert isinstance(exact, xr.DataArray)
        assert isinstance(center, xr.DataArray)

    def test_dataset_accessor_method_param(self):
        da = make_da(10, 10)
        ds = da.to_dataset(name="var")
        gdf = make_gdf([box(1.0, 1.0, 5.0, 5.0)], [1])
        exact = ds.extrs.zonal_stats(gdf, stat="mean", id_col="zone_id", method="exact")
        center = ds.extrs.zonal_stats(gdf, stat="mean", id_col="zone_id", method="center")
        assert isinstance(exact, xr.Dataset)
        assert isinstance(center, xr.Dataset)

    def test_grid_aligned_methods_agree(self):
        """For grid-aligned polygons, exact and center should produce the same mean."""
        da = make_da(10, 10)
        gdf = make_gdf([box(2.0, 2.0, 6.0, 6.0)], [1])
        exact = da.extrs.zonal_stats(gdf, stat="mean", id_col="zone_id", method="exact")
        center = da.extrs.zonal_stats(gdf, stat="mean", id_col="zone_id", method="center")
        np.testing.assert_allclose(exact.values, center.values, atol=1e-14)


# ---------------------------------------------------------------------------
# Validation against xarray-spatial
# ---------------------------------------------------------------------------

class TestMatchesXarraySpatial:
    """Validate that method='center' produces identical results to xarray-spatial."""

    @pytest.fixture
    def rotated_zones_data(self):
        """200x200 grid with 100 rotated rectangular zones."""
        da = make_da(200, 200)
        nzy, nzx, angle = 10, 10, 25
        zone_h, zone_w = 200 / nzy, 200 / nzx
        polys, ids = [], []
        zid = 1
        for iy in range(nzy):
            for ix in range(nzx):
                cx, cy = (ix + 0.5) * zone_w, (iy + 0.5) * zone_h
                rect = box(cx - zone_w * 0.4, cy - zone_h * 0.4,
                           cx + zone_w * 0.4, cy + zone_h * 0.4)
                polys.append(rotate(rect, angle, origin=(cx, cy)))
                ids.append(zid)
                zid += 1
        gdf = make_gdf(polys, ids)
        return da, gdf

    def test_mean_matches_xrspatial(self, rotated_zones_data):
        """extractrs center mean matches xarray-spatial mean to machine precision."""
        from xrspatial.zonal import stats as xrs_stats

        da, gdf = rotated_zones_data
        zones_da = rasterize_for_xrspatial(gdf, da)

        xrs_result = xrs_stats(zones=zones_da, values=da, stats_funcs=["mean"])
        xrs_means = xrs_result.set_index("zone")["mean"].sort_index()
        xrs_means = xrs_means[xrs_means.index != 0]

        ext_result = da.extrs.zonal_stats(gdf, stat="mean", id_col="zone_id", method="center")
        ext_means = ext_result.to_series().sort_index()

        common = xrs_means.index.intersection(ext_means.index)
        assert len(common) == 100

        diff = np.abs(xrs_means.loc[common].values - ext_means.loc[common].values)
        assert diff.max() < 1e-14, f"max diff = {diff.max()}"

    def test_sum_matches_xrspatial(self, rotated_zones_data):
        """extractrs center sum matches xarray-spatial sum to machine precision."""
        from xrspatial.zonal import stats as xrs_stats

        da, gdf = rotated_zones_data
        zones_da = rasterize_for_xrspatial(gdf, da)

        xrs_result = xrs_stats(zones=zones_da, values=da, stats_funcs=["sum"])
        xrs_sums = xrs_result.set_index("zone")["sum"].sort_index()
        xrs_sums = xrs_sums[xrs_sums.index != 0]

        ext_result = da.extrs.zonal_stats(gdf, stat="sum", id_col="zone_id", method="center")
        ext_sums = ext_result.to_series().sort_index()

        common = xrs_sums.index.intersection(ext_sums.index)
        diff = np.abs(xrs_sums.loc[common].values - ext_sums.loc[common].values)
        assert diff.max() < 1e-10, f"max diff = {diff.max()}"

    def test_count_matches_xrspatial(self, rotated_zones_data):
        """extractrs center count matches xarray-spatial count exactly."""
        from xrspatial.zonal import stats as xrs_stats

        da, gdf = rotated_zones_data
        zones_da = rasterize_for_xrspatial(gdf, da)

        xrs_result = xrs_stats(zones=zones_da, values=da, stats_funcs=["count"])
        xrs_counts = xrs_result.set_index("zone")["count"].sort_index()
        xrs_counts = xrs_counts[xrs_counts.index != 0]

        ext_result = da.extrs.zonal_stats(gdf, stat="count", id_col="zone_id", method="center")
        ext_counts = ext_result.to_series().sort_index()

        common = xrs_counts.index.intersection(ext_counts.index)
        diff = np.abs(xrs_counts.loc[common].values - ext_counts.loc[common].values)
        assert diff.max() == 0.0, f"count mismatch: max diff = {diff.max()}"

    def test_exact_differs_from_xrspatial(self, rotated_zones_data):
        """Verify that exact method produces different results from xarray-spatial
        on rotated polygons (confirming the methods are actually different)."""
        from xrspatial.zonal import stats as xrs_stats

        da, gdf = rotated_zones_data
        zones_da = rasterize_for_xrspatial(gdf, da)

        xrs_result = xrs_stats(zones=zones_da, values=da, stats_funcs=["mean"])
        xrs_means = xrs_result.set_index("zone")["mean"].sort_index()
        xrs_means = xrs_means[xrs_means.index != 0]

        ext_result = da.extrs.zonal_stats(gdf, stat="mean", id_col="zone_id", method="exact")
        ext_means = ext_result.to_series().sort_index()

        common = xrs_means.index.intersection(ext_means.index)
        diff = np.abs(xrs_means.loc[common].values - ext_means.loc[common].values)
        # Exact and xrspatial should NOT agree on rotated zones
        assert diff.max() > 1e-4, (
            f"exact and xrspatial unexpectedly agree (max diff = {diff.max()})"
        )
