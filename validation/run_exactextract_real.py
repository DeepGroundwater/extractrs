"""Run exactextract on real Daymet tmin file with correct CRS handling."""

import time
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.windows import from_bounds
import exactextract
import os

NC_PATH = "/mnt/ssd1/data/daymet/daymet_v4_daily_na_tmin_2024.nc"
SHP_PATH = "/mnt/ssd1/data/merit/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"
OUT_DIR = "/tmp/extractrs_validation"
BBOX_WGS84 = (-84.5, 34.5, -84.0, 35.0)
TIME_IDX = 182  # July 1 (0-indexed)

# ── Step 1: Open raster and get its exact CRS ──
subdataset = f"NETCDF:{NC_PATH}:tmin"
rast = rasterio.open(subdataset)
raster_crs = rast.crs
print(f"Raster CRS from file: {raster_crs.to_proj4()}")
print(f"Raster: {rast.width}x{rast.height}, {rast.count} bands, nodata={rast.nodata}")

# ── Step 2: Read basins, set WGS84, reproject to raster's exact CRS ──
basins = gpd.read_file(SHP_PATH, bbox=BBOX_WGS84)
basins = basins.set_crs("EPSG:4326", allow_override=True)
print(f"Loaded {len(basins)} basins in WGS84")

# Reproject to raster's CRS (LCC) using the exact CRS object from the file
basins_proj = basins.to_crs(raster_crs)
print(f"Reprojected basins to: {basins_proj.crs.to_proj4()}")

# ── Step 3: Extract bbox in projected coords and read raster subset ──
total_bounds = basins_proj.total_bounds  # [xmin, ymin, xmax, ymax]
# Add 2km padding
pad = 2000
lcc_xmin = total_bounds[0] - pad
lcc_ymin = total_bounds[1] - pad
lcc_xmax = total_bounds[2] + pad
lcc_ymax = total_bounds[3] + pad

window = from_bounds(lcc_xmin, lcc_ymin, lcc_xmax, lcc_ymax, rast.transform)
band_idx = TIME_IDX + 1  # rasterio 1-indexed

t_read = time.time()
data = rast.read(band_idx, window=window)
win_transform = rast.window_transform(window)
read_time = time.time() - t_read
print(f"Read band {band_idx} window: {data.shape}, in {read_time:.3f}s")

# Write subset as GeoTIFF with raster's exact CRS
subset_tif = os.path.join(OUT_DIR, "daymet_tmin_subset.tif")
with rasterio.open(
    subset_tif, "w", driver="GTiff",
    height=data.shape[0], width=data.shape[1], count=1,
    dtype=data.dtype, crs=raster_crs,
    transform=win_transform, nodata=-9999.0,
) as dst:
    dst.write(data, 1)
print(f"Wrote subset: {subset_tif} ({data.shape[1]}x{data.shape[0]})")

rast.close()
rast_sub = rasterio.open(subset_tif)

# Verify CRS match
print(f"\nCRS check:")
print(f"  Raster CRS:  {rast_sub.crs.to_proj4()}")
print(f"  Basins CRS:  {basins_proj.crs.to_proj4()}")
print(f"  Match: {rast_sub.crs == basins_proj.crs}")

# ── Step 4: Run exactextract ──
import warnings
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    t0 = time.time()
    results = exactextract.exact_extract(
        rast_sub, basins_proj,
        ops=["mean", "count", "sum"],
        include_cols=["COMID"],
        output="pandas",
    )
    t1 = time.time()
    ee_time = t1 - t0
    if w:
        for warning in w:
            print(f"  WARNING: {warning.message}")
    else:
        print("  No CRS warnings!")

print(f"\nexactextract: {ee_time:.4f}s for {len(basins)} basins")
print(f"\nResults preview:")
print(results.head(10).to_string())

# Save
csv_path = os.path.join(OUT_DIR, "exactextract_results.csv")
results.to_csv(csv_path, index=False)

valid = results["mean"].dropna()
print(f"\nValid: {len(valid)}/{len(results)}")
if len(valid) > 0:
    print(f"  Mean tmin: {valid.mean():.4f}°C")
    print(f"  Range: [{valid.min():.4f}, {valid.max():.4f}]°C")

# ── Step 5: Benchmark ──
times = []
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for _ in range(20):
        t0 = time.time()
        exactextract.exact_extract(
            rast_sub, basins_proj,
            ops=["mean"], include_cols=["COMID"], output="pandas"
        )
        times.append(time.time() - t0)

print(f"\n=== exactextract Benchmark (20 runs, {len(basins)} basins) ===")
print(f"  Mean:   {np.mean(times)*1000:.2f}ms")
print(f"  Std:    {np.std(times)*1000:.2f}ms")
print(f"  Min:    {np.min(times)*1000:.2f}ms")

rast_sub.close()
