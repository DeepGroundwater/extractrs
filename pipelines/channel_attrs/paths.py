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
# half_width = max(E_POS_M, ALPHA_HALF * bankfull_width).
E_POS_M = 10.0              # minimum corridor half-width (m); floor for headwater streams
ALPHA_HALF = 1.5            # half-width per unit bankfull width (banks + hyporheic margin)
# WRF-Hydro CONUS order-based bottom-channel-width table (m), Strahler orders 1–10.
# Source: NCAR/wrf_hydro_arcgis_preprocessor wrf_hydro_functions.py ``Mannings_Bw``
# (LR 7/01/2020, confirmed JMC 6/18/21). These are trapezoidal channel BOTTOM widths
# (narrower than bankfull surface width). Used as the order fallback in
# ``_resolve_width_m`` when observed/modelled width columns are absent.
# With the 1.5× hinge rule: floor dominates for orders 1–9 (1.5×26=39 m < 100 m);
# only order 10 clears it (1.5×110=165 m). Observed/modelled widths remain the
# primary corridor-scaling source.
WRF_HYDRO_BW_M = {
    1:   1.6,
    2:   2.4,
    3:   3.5,
    4:   5.3,
    5:   7.4,
    6:  11.0,
    7:  14.0,
    8:  16.0,
    9:  26.0,
    10: 110.0,
}

CROSSWALK_BUFFER_M = 300.0      # NHD->MERIT matching envelope (MERIT lateral error 100-300 m)
CROSSWALK_TOP_K = 40            # mirror Wade et al. 2025 table shape

# Width cap for corridor scaling: widths above this are SWORD estuary/lake
# artefacts, not channel widths.  SWORD measured max ~61 km on pre-bugfix1
# fabrics; observed max on bugfix1 is 16.7 km.  Cap at 3 km: a >3 km
# "channel" is not a corridor-extraction target.
WIDTH_CAP_M = 3000.0
