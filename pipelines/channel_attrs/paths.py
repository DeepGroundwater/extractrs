"""Canonical paths for the channel-attributes pipeline (leakance gate, Phase A).

Note: MERIT_RIV is at ddr/data/merit/ (not ddr/data/ as the plan draft said);
corrected against the actual filesystem layout on 2026-07-04.
"""
from pathlib import Path

RAW = Path("/mnt/ssd1/data/channel_attrs/raw")
DERIVED = Path("/mnt/ssd1/data/channel_attrs/derived")

MERIT_RIV = Path(
    "/home/tbindas/projects/ddr/data/merit/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"
)
GLOBAL_NC = Path("/home/tbindas/projects/ddr/data/merit_global_attributes_v2.nc")
STATS_DIR = Path("/home/tbindas/projects/ddr/data/statistics")

OUT_NC = Path("/home/tbindas/projects/ddr/data/merit_channel_attributes_v1.nc")
OUT_STATS = Path("/home/tbindas/projects/ddr/data/statistics/merit_channel_attributes_v1.json")

# Equal-area CRS for all buffering/length math (CONUS Albers).
CRS_EQUAL_AREA = "EPSG:5070"

CORRIDOR_HALF_WIDTH_M = 100.0   # StreamCat precedent (Hill et al. 2016)
CORRIDOR_WIDE_M = 200.0         # flat-valley widening (Amatulli 2022 error tail)
CROSSWALK_BUFFER_M = 300.0      # NHD->MERIT matching envelope (MERIT lateral error 100-300 m)
CROSSWALK_TOP_K = 40            # mirror Wade et al. 2025 table shape
