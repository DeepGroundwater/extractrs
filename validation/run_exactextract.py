"""Run exactextract on the test data and save results."""

import time
import geopandas as gpd
import rasterio
import exactextract
import numpy as np
import os

OUT_DIR = "/tmp/extractrs_validation"
TIF_PATH = os.path.join(OUT_DIR, "test_prcp.tif")
SHP_PATH = "/mnt/ssd1/data/merit/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"

BBOX_WGS84 = (-84.5, 34.5, -84.0, 35.0)
DAYMET_LCC = "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +datum=WGS84 +units=m"

# Read basins directly from shapefile with bbox filter
basins = gpd.read_file(SHP_PATH, bbox=BBOX_WGS84)
basins = basins.set_crs("EPSG:4326", allow_override=True)
basins = basins.head(100)
print(f"Loaded {len(basins)} basins, CRS={basins.crs}")

# Reproject basins to LCC to match raster
basins_lcc = basins.to_crs(DAYMET_LCC)
print(f"Reprojected to LCC")

# Open raster
rast = rasterio.open(TIF_PATH)
print(f"Raster: {rast.width}x{rast.height}, CRS={rast.crs}, nodata={rast.nodata}")

# Run exactextract with LCC basins — time it
t0 = time.time()
results = exactextract.exact_extract(
    rast,
    basins_lcc,
    ops=["mean", "count", "sum"],
    include_cols=["COMID"],
    output="pandas",
)
t1 = time.time()
exactextract_time = t1 - t0

print(f"\nexactextract completed in {exactextract_time:.4f}s for {len(basins)} basins")
print(f"\nResults preview:")
print(results.head(20).to_string())

# Save results
csv_path = os.path.join(OUT_DIR, "exactextract_results.csv")
results.to_csv(csv_path, index=False)
print(f"\nSaved to {csv_path}")

# Summary stats
valid = results["mean"].dropna()
print(f"\nSummary:")
print(f"  Valid basins: {len(valid)}/{len(results)}")
if len(valid) > 0:
    print(f"  Mean prcp: {valid.mean():.4f}")
    print(f"  Min: {valid.min():.4f}")
    print(f"  Max: {valid.max():.4f}")

# Benchmark: run 20 times
times = []
for _ in range(20):
    t0 = time.time()
    exactextract.exact_extract(rast, basins_lcc, ops=["mean"], include_cols=["COMID"], output="pandas")
    times.append(time.time() - t0)

print(f"\n=== exactextract Benchmark (20 runs) ===")
print(f"  Mean: {np.mean(times):.4f}s")
print(f"  Std:  {np.std(times):.4f}s")
print(f"  Min:  {np.min(times):.4f}s")
