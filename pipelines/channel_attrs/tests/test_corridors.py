import geopandas as gpd
from shapely.geometry import LineString
from pipelines.channel_attrs.corridors import build_corridors


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
