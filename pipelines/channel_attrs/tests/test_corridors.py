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


def test_order_fallback_width_matches_spec_table():
    # w(omega) = 2 * 1.9^(omega-1): order 1 -> 2 m, order 6 -> ~50 m (spec §hinge).
    assert abs(order_bankfull_width_m(1) - 2.0) < 1e-9
    assert abs(order_bankfull_width_m(6) - 2.0 * 1.9**5) < 1e-9
    assert 45 < order_bankfull_width_m(6) < 55


def test_scaled_half_width_hinge():
    # Small width -> pinned to the 100 m positional-error floor; large width scales.
    assert scaled_half_width_m(2.0) == 100.0        # 1.5*2 = 3 << 100
    assert scaled_half_width_m(400.0) == 600.0      # 1.5*400 = 600 > 100


def test_scaled_corridor_floors_small_streams():
    # Order 3 (w ~7 m) is far below the hinge: half_width stays at the floor.
    gdf = gpd.GeoDataFrame(
        {"COMID": [1], "order": [3]},
        geometry=[LineString([(0, 0), (1000, 0)])],
        crs="EPSG:5070",
    )
    out = build_scaled_corridors(gdf)
    assert abs(out["half_width_m"].iloc[0] - paths.E_POS_M) < 1e-6
    assert out.crs.to_epsg() == 5070


def test_scaled_corridor_widens_large_rivers():
    # Order 9 clears the hinge: w = 2*1.9^8 ~ 340 m, half_width = 1.5*w ~ 510 m.
    gdf = gpd.GeoDataFrame(
        {"COMID": [1], "order": [9]},
        geometry=[LineString([(0, 0), (1000, 0)])],
        crs="EPSG:5070",
    )
    out = build_scaled_corridors(gdf)
    assert out["half_width_m"].iloc[0] > 400.0


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
