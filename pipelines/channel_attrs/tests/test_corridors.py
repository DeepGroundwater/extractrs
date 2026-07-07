import numpy as np
import geopandas as gpd
from shapely.geometry import LineString
from pipelines.channel_attrs import paths
from pipelines.channel_attrs.corridors import (
    build_corridors,
    build_scaled_corridors,
    order_bankfull_width_m,
    scaled_half_width_m,
)


def test_build_corridors_buffers_and_preserves_ids():
    gdf = gpd.GeoDataFrame(
        {"COMID": [1, 2]},
        geometry=[LineString([(0, 0), (1000, 0)]), LineString([(1000, 0), (1000, 800)])],
        crs="EPSG:5070",
    )
    out = build_corridors(gdf, half_width_m=100.0)
    assert list(out["COMID"]) == [1, 2]
    assert out.crs.to_epsg() == 5070
    # 1000 m line buffered 100 m: area ~ 1000*200 + pi*100^2 (round caps)
    assert abs(out.geometry.iloc[0].area - (1000 * 200 + 3.14159 * 100**2)) < 500


def test_build_corridors_reprojects_from_4326():
    # Lon/lat input must be reprojected to the equal-area CRS before buffering,
    # otherwise buffer(100) would be 100 *degrees*. ~0.01 deg lon at 40N is
    # ~850 m on the ground, so the corridor area lands near 850*200 + a cap.
    gdf = gpd.GeoDataFrame(
        {"COMID": [1]},
        geometry=[LineString([(-95.0, 40.0), (-94.99, 40.0)])],
        crs="EPSG:4326",
    )
    out = build_corridors(gdf, half_width_m=100.0)
    assert out.crs.to_epsg() == 5070
    assert out.geometry.iloc[0].area > 100_000  # degrees-buffer would be ~1e-4


def test_order_fallback_width_matches_wrf_hydro_table():
    # Values from paths.WRF_HYDRO_BW_M (NCAR wrf_hydro_functions.py Mannings_Bw).
    from pipelines.channel_attrs import paths
    for order, expected in paths.WRF_HYDRO_BW_M.items():
        assert abs(order_bankfull_width_m(order) - expected) < 1e-9, (
            f"order {order}: got {order_bankfull_width_m(order)}, expected {expected}"
        )

def test_order_fallback_clamps_to_table_bounds():
    from pipelines.channel_attrs import paths
    # order < 1 → clamped to order 1
    assert order_bankfull_width_m(0) == paths.WRF_HYDRO_BW_M[1]
    # order > 10 → clamped to order 10
    assert order_bankfull_width_m(12) == paths.WRF_HYDRO_BW_M[10]


def test_scaled_half_width_hinge():
    # Small width -> pinned to the 10 m floor; large width scales.
    assert scaled_half_width_m(2.0) == paths.E_POS_M    # 1.5*2 = 3 < 10
    assert scaled_half_width_m(400.0) == 600.0           # 1.5*400 = 600 > 10


import pytest

@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_scaled_corridor_floors_small_streams(order):
    # Orders 1-4: 1.5 * Bw < E_POS_M (e.g. order 4: 1.5*5.3=7.95 < 10 m).
    gdf = gpd.GeoDataFrame(
        {"COMID": [1], "order": [order]},
        geometry=[LineString([(0, 0), (1000, 0)])],
        crs="EPSG:5070",
    )
    out = build_scaled_corridors(gdf)
    assert abs(out["half_width_m"].iloc[0] - paths.E_POS_M) < 1e-6
    assert out.crs.to_epsg() == 5070


@pytest.mark.parametrize("order", [5, 6, 7, 8, 9, 10])
def test_scaled_corridor_widens_large_rivers(order):
    # Orders 5-10: 1.5 * Bw > E_POS_M — half_width equals 1.5 * WRF-Hydro Bw.
    bw = paths.WRF_HYDRO_BW_M[order]
    gdf = gpd.GeoDataFrame(
        {"COMID": [1], "order": [order]},
        geometry=[LineString([(0, 0), (1000, 0)])],
        crs="EPSG:5070",
    )
    out = build_scaled_corridors(gdf)
    assert np.isclose(out["half_width_m"].iloc[0], paths.ALPHA_HALF * bw)


def test_scaled_corridor_prefers_observed_width():
    # Observed width overrides the small order fallback (order 2 -> ~4 m).
    gdf = gpd.GeoDataFrame(
        {"COMID": [1], "order": [2], "channel_width_obs": [400.0]},
        geometry=[LineString([(0, 0), (1000, 0)])],
        crs="EPSG:5070",
    )
    out = build_scaled_corridors(gdf)
    assert abs(out["half_width_m"].iloc[0] - 600.0) < 1e-6  # 1.5 * 400


def test_scaled_corridor_ignores_invalid_observed_width():
    # NaN/zero observed width falls back to the order estimate, not to 0/NaN.
    gdf = gpd.GeoDataFrame(
        {"COMID": [1, 2], "order": [9, 9], "channel_width_obs": [np.nan, 0.0]},
        geometry=[LineString([(0, 0), (1000, 0)]), LineString([(0, 0), (1000, 0)])],
        crs="EPSG:5070",
    )
    out = build_scaled_corridors(gdf)
    expected = scaled_half_width_m(order_bankfull_width_m(9))
    assert np.allclose(out["half_width_m"].to_numpy(), expected)


def test_width_cap_applied():
    # An estuary-scale width (> WIDTH_CAP_M) must be clipped to the cap.
    gdf = gpd.GeoDataFrame(
        {"COMID": [1], "order": [1], "channel_width_obs": [paths.WIDTH_CAP_M * 2]},
        geometry=[LineString([(0, 0), (1000, 0)])],
        crs="EPSG:5070",
    )
    out = build_scaled_corridors(gdf)
    expected_hw = scaled_half_width_m(paths.WIDTH_CAP_M)  # cap applied before hinge
    assert np.isclose(out["half_width_m"].iloc[0], expected_hw), (
        f"Expected {expected_hw}, got {out['half_width_m'].iloc[0]}"
    )


def test_merge_preserves_row_count():
    # Merging width columns onto a mini-riv frame must not change the row count
    # (left-join semantics: all MERIT reaches kept, NaN where no width data).
    import pandas as pd

    gdf = gpd.GeoDataFrame(
        {"COMID": [1, 2, 3], "order": [3, 5, 8]},
        geometry=[LineString([(0, 0), (1000, 0)])] * 3,
        crs="EPSG:5070",
    )
    cw = pd.DataFrame({"COMID": [2], "channel_width_obs": [300.0]})
    bf = pd.DataFrame({"COMID": [3], "bankfull_width": [150.0]})
    merged = gdf.merge(cw, on="COMID", how="left").merge(bf, on="COMID", how="left")
    assert len(merged) == 3, f"Row count changed: {len(merged)}"
    assert np.isnan(merged.loc[merged.COMID == 1, "channel_width_obs"].iloc[0])
    assert merged.loc[merged.COMID == 2, "channel_width_obs"].iloc[0] == 300.0
