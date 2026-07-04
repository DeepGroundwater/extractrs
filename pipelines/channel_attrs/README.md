# Channel-Corridor Attribute Pipeline (Phase A — Leakance Gate)

Extracts per-COMID channel-corridor and groundwater attributes for the ddrs
leakance feasibility gate. Final output:
`/home/tbindas/projects/ddr/data/merit_channel_attributes_v1.nc` (same schema
as `merit_global_attributes_v2.nc`).

Spec: `ddrs docs/superpowers/specs/2026-07-04-leakance-gate-program-design.md`
Plan: `ddrs docs/superpowers/plans/2026-07-04-phase-a-channel-attributes.md`

---

## Environment

- Python 3.13 venv at `.venv-pipelines` (`uv venv --python 3.13 .venv-pipelines`)
- `extractrs` installed from PyPI (v0.1.3) — no maturin build needed for py3.13
- Install: `uv pip install --python .venv-pipelines/bin/python extractrs geopandas pyogrio rioxarray xarray netCDF4 pandas pyarrow requests tqdm pytest shapely`
- Run tests: `.venv-pipelines/bin/python -m pytest pipelines/channel_attrs/tests/ -v`

---

## Schema Reconnaissance (Task 0)

### Global attributes NC schema (`merit_global_attributes_v2.nc`)

```
dims:  {COMID: 2939404}
vars:  30 total
  COMID:   (COMID,) int64    — the COMID coordinate variable
  T_sand:  (COMID,) float64  — example attribute (all 29 attrs same shape/type)
  ...
  log10_uparea: (COMID,) float64
```

All attribute variables are **float64 on the single `(COMID,)` dimension**.
The `COMID` variable is `int64`. `merit_channel_attributes_v1.nc` must match
this schema exactly (no extra dims, no groups).

### Stats JSON convention (`ddr/data/statistics/merit_attribute_statistics_merit_global_attributes_v2.nc.json`)

Per-variable dict with exactly **6 keys**:
```json
{
  "T_sand": {"min": 0.0, "max": 99.0, "mean": 48.64, "std": 23.80, "p10": 19.0, "p90": 84.0},
  ...
}
```
Task 10 must emit the identical structure (`min`, `max`, `mean`, `std`, `p10`,
`p90`) over finite values only.

### MERIT pfaf-7 flowline shapefile fields

File: `/home/tbindas/projects/ddr/data/merit/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp`
(NOT at `ddr/data/` as the plan draft said — actual location is `ddr/data/merit/`)

```
fields: COMID, lengthkm, lengthdir, sinuosity, slope, uparea, order,
        strmDrop_t, slope_taud, NextDownID, maxup, up1, up2, up3, up4
dtypes: int64, float64, float64, float64, float64, float64, int64,
        float64, float64, int64, int64, int64, int64, int64, int64
CRS:    EPSG:4326
count:  346,327 features (CONUS pfaf-7 subset)
```

Key fields for downstream tasks: `COMID` (int64), `lengthkm` (float64),
`uparea` (float64), `order` (int64). Note: `MERIT COMID` and `NHDPlus COMID`
are **unrelated ID spaces** sharing the same name (spec §A0).

`paths.py` uses `CRS_EQUAL_AREA = "EPSG:5070"` (CONUS Albers) for all
buffering and length operations; all reads reproject from EPSG:4326.

---

## Download Manifest

### Summary table (as of 2026-07-04)

| Dataset | Target | URL | Size | Status |
|---|---|---|---|---|
| MERIT-SWORD translation tables | `raw/merit_sword/ms_translate.zip` | `https://zenodo.org/api/records/13152826/files/ms_translate.zip/content` | 17 MB | **COMPLETE** |
| SWORD v16 netCDF (all continents) | `raw/sword/SWORD_v16_netcdf.zip` | `https://zenodo.org/api/records/10013982/files/SWORD_v16_netcdf.zip/content` | 1.57 GB | **IN PROGRESS** |
| NHDPlusV2 National Seamless GDB | `raw/nhdplusv2/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07.7z` | `https://dmap-data-commons-ow.s3.amazonaws.com/NHDPlusV21/Data/NationalData/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07.7z` | ~7 GB | **IN PROGRESS** |
| StreamCat `pctimp2019` (21 regions) | `raw/streamcat/pctimp2019_Region*.json` | `https://api.epa.gov/StreamCat/streams/metrics?name=pctimp2019&region=<R>` | ~6.8 MB/region | **IN PROGRESS** |
| Zarrabi bankfull geometry | `raw/zarrabi/Bankfull_Meanflow_CONUS.txt` | `https://zenodo.org/api/records/13883263/files/Bankfull_Meanflow_CONUS.txt/content` | 207 MB | **COMPLETE** |
| Zell & Sanford 2020 WTD | `raw/zs_wtd/Output_CONUS_trans_dtw.zip` | `https://www.sciencebase.gov/catalog/file/get/631405c5d34e36012efa3190?f=__disk__3b%2Fa3%2F73%2F3ba373573a407d90329925e56a1c82e55daadbe1` | 918 MB | **IN PROGRESS** |
| Fan 2013 WTD North America | `raw/fan_wtd/NAMERICA_WTD_annualmean.nc` | `http://thredds-gfnl.usc.es/thredds/fileServer/GLOBALWTDFTP/annualmeans/NAMERICA_WTD_annualmean.nc` | 66 MB | **COMPLETE** |
| USGS BFI (NHD-indexed parquet+zip) | `raw/bfi/BFI_CONUS.zip` | `https://www.sciencebase.gov/catalog/file/get/5669a8e3e4b08895842a1d4f?f=__disk__8b%2Fd6%2F90%2F8bd6901b9d1bf7940a8219dc500f0b04cb381206` | 17 MB | **COMPLETE** |
| USGS Principal Aquifers | `raw/aquifers/aquifers_us.zip` | `https://www.sciencebase.gov/catalog/file/get/63140610d34e36012efa385d?f=__disk__38%2F75%2F7d%2F38757d8db3921426b941b8efe3599b5c03c56917` | 7.5 MB | **COMPLETE** |
| GFPLAIN250m (global TIFF rar) | `raw/gfplain/GFPLAIN250m_TIFF.rar` | `https://ndownloader.figshare.com/files/12186356` | 33 MB | **COMPLETE** |

### Dataset-specific notes

**MERIT-SWORD (Zenodo 13152826):**
`ms_translate.zip` (17 MB on disk, Zenodo lists 20.9 MB — difference is
rounding in Zenodo's approximation; zip validates with 127 entries including
`ms_translate/mb_to_sword/` and `ms_translate/sword_to_mb/` subdirectories
of per-pfaf-2-region NetCDF files). The Wade et al. 2025 translation table
maps MERIT-Basins COMID ↔ SWORD reach_id with `part_len` overlap lengths.
CONUS = pfaf-2 regions 71–78 (files inside the zip; confirm exact file names
at Task 2 time).

**SWORD (Zenodo 10013982):**
The plan refers to "SWORD v2" but the current Zenodo record is **SWORD v16**
(public version numbering restarted at the SWOT mission launch). `SWORD_v16_netcdf.zip`
contains all-continent netCDF files; the NA file will be `na_sword_v16.nc`.
Download in progress (~1.57 GB total).

**NHDPlusV2:**
Downloading the national seamless GDB (`NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07.7z`,
~7 GB) — per-VPU shapefiles would be smaller but the national GDB avoids
managing 21 separate downloads. Download in progress; expect several hours.
The `.7z` archive requires `p7zip` to extract.

**StreamCat:**
The EPA FTP (`gaftp.epa.gov/EPADataCommons/ORD/NHDPlusLandscapeAttributes/StreamCat/HydroRegions/`)
has been retired. The REST API moved to `https://api.epa.gov/StreamCat/streams`.
The metric name in the new API is **lowercase** (`pctimp2019`), returning
columns `pctimp2019cat`, `pctimp2019ws`, `pctimp2019catrp100`,
`pctimp2019wsrp100` per COMID. The column we want is `pctimp2019catrp100`
(= plan's `PctImp2019Rp100Cat`, 100m riparian buffer, local catchment).
Downloading all 21 CONUS HUC2 regions as JSON via nohup script.
Regions: Region01, Region02, Region03N, Region03S, Region03W, Region04–Region18.

**Zell & Sanford 2020 WTD:**
ScienceBase item 631405c5d34e36012efa3190 contains MODFLOW-6 output archives.
`Output_CONUS_trans_dtw.zip` (918 MB) is the depth-to-water product; the
`trans_dtw` name suggests it may contain both transmissivity and DTW rasters.
Inspect archive contents after download to identify the usable raster layer(s).
Decision recorded: **Zell & Sanford 2020 is primary WTD source**; Fan 2013
(`NAMERICA_WTD_annualmean.nc`) is the secondary/cross-check source (as in the
plan's spec §3 sensitivity analysis). If the Zell & Sanford archive lacks a
directly usable single-band DTW raster, note here and promote Fan 2013 to primary.

**BFI:**
Primary URL `http://water.usgs.gov/GIS/dsdl/bfi48grd.zip` is dead
(water.usgs.gov port 80 refused — server decommissioned). Found a better
alternative: `BFI_CONUS.zip` from the NHDPlus attributes ScienceBase item
(10.5066/P91LFFN1 analog, USGS:5669a8e3e4b08895842a1d4f) — this is the BFI
already indexed to NHD COMIDs, which means Task 9 can transfer directly to
MERIT via the Task-3 crosswalk without needing rioxarray zonal stats. Preferred
over the raw bfi48grd grid. Downloaded (17 MB). Also attempted `bfi_cat.parquet`
from ScienceBase manager URL but got redirect/error (4.2 KB).

**GFPLAIN250m:**
figshare article 6665165 has two files: `GFPLAIN250m TIFF.rar` (33 MB) and
`GFPLAIN250m ASCII.rar` (47 MB). These are **global** archives; the NA tile
is inside the rar as a continental raster (file name to confirm on extraction).
No standalone NA tile exists as a separate download — the rar contains all
continents. Requires `unrar` or `python-unrar` to extract. Downloaded TIFF rar.

**Fan 2013 THREDDS:**
Catalog at `http://thredds-gfnl.usc.es/thredds/catalog/GLOBALWTDFTP/annualmeans/catalog.html`.
Available files: AFRICA, EURASIA, NAMERICA, OCEANIA, SAMERICA `_WTD_annualmean.nc`.
fileServer URL: `http://thredds-gfnl.usc.es/thredds/fileServer/GLOBALWTDFTP/annualmeans/NAMERICA_WTD_annualmean.nc`
Downloaded (66 MB). Server is slow (WebFetch times out; wget -c works directly).

---

## Running the pipeline

```bash
# Activate environment
source .venv-pipelines/bin/activate

# Run tests
python -m pytest pipelines/channel_attrs/tests/ -v

# Task 1: Corridor geometries (once MERIT shapefile accessible)
python -m pipelines.channel_attrs.corridors

# (subsequent tasks depend on downloads completing)
```
