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
MERIT_CAT = Path(
    "/home/tbindas/projects/ddr/data/merit/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"
)
GLOBAL_NC = Path("/home/tbindas/projects/ddr/data/merit_global_attributes_v2.nc")
STATS_DIR = Path("/home/tbindas/projects/ddr/data/statistics")

OUT_NC = Path("/home/tbindas/projects/ddr/data/merit_channel_attributes_v1.nc")
OUT_STATS = Path("/home/tbindas/projects/ddr/data/statistics/merit_channel_attributes_v1.json")

# Equal-area CRS for all buffering/length math (CONUS Albers).
CRS_EQUAL_AREA = "EPSG:5070"

# Corridor buffer scaling (spec: specs/2026-07-04-corridor-buffer-scaling.md).
# half_width = max(E_POS_M, ALPHA_HALF * bankfull_width). Below ~order 6 the
# positional-error floor dominates and the scaled set equals the fixed 100 m set.
E_POS_M = 100.0             # positional-error floor (Amatulli 2022; Hill et al. 2016)
ALPHA_HALF = 1.5           # half-width per unit bankfull width (banks + hyporheic margin)
ORDER_WIDTH_BASE_M = 2.0   # order-1 bankfull width, m (Downing et al. 2012)
ORDER_WIDTH_RATIO = 1.9    # bankfull-width growth per Strahler order (L&M 1953 + Horton)

CROSSWALK_BUFFER_M = 300.0      # NHD->MERIT matching envelope (MERIT lateral error 100-300 m)
CROSSWALK_TOP_K = 40            # mirror Wade et al. 2025 table shape

# Width cap for corridor scaling: widths above this are SWORD estuary/lake
# artefacts, not channel widths.  SWORD measured max ~61 km on pre-bugfix1
# fabrics; observed max on bugfix1 is 16.7 km.  Cap at 3 km: a >3 km
# "channel" is not a corridor-extraction target.
WIDTH_CAP_M = 3000.0
