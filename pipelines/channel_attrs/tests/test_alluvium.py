"""Tests for alluvium.overlay_chunk."""
import numpy as np
import pandas as pd
import geopandas as gpd
import pytest
from shapely.geometry import box

from pipelines.channel_attrs.alluvium import overlay_chunk


def _corridors(boxes, comids, crs="EPSG:5070"):
    return gpd.GeoDataFrame({"COMID": comids}, geometry=boxes, crs=crs)


def _alluvial(polys, crs="EPSG:5070"):
    gdf = gpd.GeoDataFrame({"is_alluvial": [1] * len(polys)}, geometry=polys, crs=crs)
    return gdf


def test_overlay_chunk_full_coverage():
    # Corridor entirely inside the alluvial polygon → fraction = 1.0.
    corridor = box(0, 0, 100, 100)        # 10 000 m²
    alluvium = box(-50, -50, 200, 200)    # wraps the corridor completely
    result = overlay_chunk(_corridors([corridor], [1]), _alluvial([alluvium]))
    assert len(result) == 1
    assert result["COMID"].iloc[0] == 1
    assert result["alluvium_fraction"].iloc[0] == pytest.approx(1.0, abs=1e-6)


def test_overlay_chunk_no_coverage():
    # Corridor entirely outside alluvial polygon → empty result (caller fills 0).
    corridor = box(0, 0, 100, 100)
    alluvium = box(500, 500, 600, 600)    # nowhere near the corridor
    result = overlay_chunk(_corridors([corridor], [1]), _alluvial([alluvium]))
    assert len(result) == 0


def test_overlay_chunk_partial_coverage():
    # Corridor is a 100×100 box; alluvial polygon covers exactly the left half.
    corridor = box(0, 0, 100, 100)        # 10 000 m²
    alluvium = box(-10, -10, 50, 110)     # left half + margin
    result = overlay_chunk(_corridors([corridor], [1]), _alluvial([alluvium]))
    assert len(result) == 1
    frac = result["alluvium_fraction"].iloc[0]
    assert 0.45 < frac < 0.55


def test_overlay_chunk_preserves_comid():
    # COMID values should pass through unchanged.
    corridors = _corridors(
        [box(0, 0, 100, 100), box(200, 200, 300, 300)],
        comids=[42, 99],
    )
    alluvium = _alluvial([box(-500, -500, 1000, 1000)])
    result = overlay_chunk(corridors, alluvium)
    assert set(result["COMID"].tolist()) == {42, 99}


def test_overlay_chunk_multiple_comids_independent():
    # Two corridors at different alluvium coverage levels; values are independent.
    corridors = _corridors(
        [box(0, 0, 100, 100), box(1000, 1000, 1100, 1100)],
        comids=[1, 2],
    )
    # Alluvium covers corridor 1 fully but not corridor 2 at all.
    alluvium = _alluvial([box(-10, -10, 110, 110)])
    result = overlay_chunk(corridors, alluvium)
    assert len(result) == 1   # only corridor 1 has alluvium
    assert result["COMID"].iloc[0] == 1
    assert result["alluvium_fraction"].iloc[0] == pytest.approx(1.0, abs=1e-4)
