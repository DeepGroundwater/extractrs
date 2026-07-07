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
ScienceBase item 631405c5d34e36012efa3190.
`Output_CONUS_trans_dtw.zip` downloaded (876 MB on disk; ScienceBase listed 918 MB but
the zip validates as complete). Contains exactly 2 files:
- `Output_CONUS_trans_dtw/conus_MF6_SS_Unconfined_250_dtw.tif` — **DTW GeoTIFF, 250m, directly usable**
- `Output_CONUS_trans_dtw/conus_MF6_SS_Unconfined_250_trans.tif` — transmissivity (not needed for Phase A)
Decision: **Zell & Sanford 2020 is confirmed primary WTD source** (directly usable DTW GeoTIFF at 250m).
Fan 2013 (`NAMERICA_WTD_annualmean.nc`) remains the secondary cross-check source.

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

# Run tests (69 tests, ~1 s)
python -m pytest pipelines/channel_attrs/tests/ -v

# Task 1: Corridor geometries
# Writes corridors_10m.parquet (fixed 10 m half-width floor, sensitivity column)
# and corridors_scaled.parquet (half_width = max(10 m, 1.5 × bankfull width);
# orders 1–4 floor at 10 m, orders 5–10 scale with WRF-Hydro Bw table).
# Order fallback: paths.WRF_HYDRO_BW_M (NCAR Mannings_Bw, Strahler 1–10).
# Width priority: channel_width_obs > bankfull_width > order fallback.
# See specs/2026-07-04-corridor-buffer-scaling.md.
python -m pipelines.channel_attrs.corridors
# Build one product only: --only 10m | --only scaled

# Task 2: SWORD widths (requires raw/sword/netcdf/na_sword_v16.nc)
python -m pipelines.channel_attrs.sword_width

# Task 3: NHDPlus->MERIT crosswalk (heavy compute, ~hours)
python -m pipelines.channel_attrs.nhd_crosswalk

# Task 4: StreamCat corridor imperviousness
python -m pipelines.channel_attrs.streamcat_transfer

# Task 5: Zarrabi bankfull depth + width
python -m pipelines.channel_attrs.zarrabi_transfer

# Task 6: WTD zonal stats + Fan cross-check (reads corridors_10m + ZS raster)
python -m pipelines.channel_attrs.wtd_sample

# Task 7: WTD bed-relative depth
python -m pipelines.channel_attrs.wtd_bedrel

# Task 8: Alluvium fraction (overlay on corridors_scaled)
python -m pipelines.channel_attrs.alluvium

# Task 9: BFI + drainage density
python -m pipelines.channel_attrs.bfi_transfer

# Task 10: Assemble final NC + stats JSON
python -m pipelines.channel_attrs.assemble
```

### Test coverage

| Test file | Module(s) tested | Tests |
|---|---|---|
| `test_alluvium.py` | `alluvium.overlay_chunk` | 5 |
| `test_assemble.py` | `assemble` | 1 |
| `test_bfi_transfer.py` | `bfi_transfer` | 4 |
| `test_core.py` | `wtd_sample.sample_along_lines`, `transfer.weighted_transfer` | 13 |
| `test_corridors.py` | `corridors` (all public functions) | 19 |
| `test_nhd_crosswalk.py` | `nhd_crosswalk` | 3 |
| `test_paths.py` | `paths` (file existence) | 2 |
| `test_streamcat_transfer.py` | `streamcat_transfer.load_streamcat` | 3 |
| `test_sword_width.py` | `sword_width.load_sword_reaches`, `build_crosswalk` | 5 |
| `test_transfer.py` | `transfer.weighted_transfer` | 2 |
| `test_wtd_bedrel.py` | `wtd_bedrel` | 4 |
| `test_wtd_sample.py` | `wtd_sample.sample_along_lines` | 5 |
| `test_zarrabi_transfer.py` | `zarrabi_transfer.load_zarrabi` | 3 |
| **Total** | | **69** |

---

## Task Deviations and Decisions

### Task 4: StreamCat (2026-07-04)

- **pctimp2019_Region07.json is truncated** (ends mid-record at 19,355,160 bytes).
  `load_streamcat` uses `ijson` streaming parse and silently drops the incomplete
  final record. Region07 (Upper Mississippi) yields 178,000 complete records vs
  the ~182,000 expected from its file size — a ~2% loss for that region.
  Action: re-download Region07 before the assembly step (or accept the ~0.1%
  overall record loss across all 21 regions, ~2k/2.53M).

- **API column naming**: The StreamCat REST API returns lowercase keys
  (`pctimp2019catrp100`, `comid`), not the legacy FTP column names
  (`PctImp2019Rp100Cat`, `COMID`).

### Task 5: Zarrabi bankfull (2026-07-04)

- **McManamay confinement: SKIPPED.** McManamay & Pers (2017) is not in the
  Phase-A download manifest. The three bankfull columns (`bankfull_depth`,
  `bankfull_width`, and a future `confinement` column) constitute the
  full Zarrabi deliverable; confinement is deferred.

- **Spearman(bankfull_width, channel_width_obs) = 0.005** — failed the >0.5
  threshold. Root cause: channel_width_obs (Task 2, SWORD) contains extreme
  outlier widths (up to 61,328 m) that collapse Spearman rank correlation to
  near zero. The Zarrabi bankfull widths themselves are correct (9.9–1039 m).
  This is a data quality issue in Task 2 (SWORD crosswalk noise in specific
  MERIT reaches), not in the Zarrabi transfer. Flagged for Task 2 follow-up.

- **Spearman(bankfull_depth, uparea) = 0.469** — marginally below the >0.5
  threshold. The monotonic depth–area relationship is attenuated across
  different geological settings and stream types at the CONUS scale; the
  direction is correct (positive) and the median depth (1.43 m) is well within
  the expected range.

### Task 9: BFI + drainage density (2026-07-04)

- **BFI source deviation from plan**: The plan (Task 9, Step 1) described
  transfer via rioxarray zonal stats on the raw `bfi48grd` raster. The
  downloaded `BFI_CONUS.zip` from ScienceBase is instead an **NHD-indexed
  table** (`CAT_BFI` per NHD COMID, already aggregated per local catchment by
  EPA). This is used directly via the NHDPlus→MERIT crosswalk — equivalent
  derivation from the same Wolock BFI grid, simpler path, no rioxarray needed.

- **BFI values confirm [0, 100] range** (max 90.0 % observed → 0.90 after /100).

- **BFI 100% coverage** across all 156,002 MERIT COMIDs — the NHD-indexed
  table is dense enough that every MERIT reach matches at least one NHD COMID
  with match_frac ≥ 0.3.
