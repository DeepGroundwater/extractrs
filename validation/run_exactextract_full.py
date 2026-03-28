"""Run exactextract on CONUS MERIT basins with subprocess isolation per chunk."""

import subprocess, sys, time, os, gc
import pandas as pd
import numpy as np

SHP_PATH = "/mnt/ssd1/data/merit/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"
OUT_DIR = "/tmp/extractrs_validation"
TIF_PATH = os.path.join(OUT_DIR, "daymet_tmin_full.tif")
CHUNK = 25000

# Known bad basin indices (cause segfault in exactextract)
BAD_INDICES = {60764}

os.makedirs(OUT_DIR, exist_ok=True)

# Ensure GeoTIFF exists
if not os.path.exists(TIF_PATH):
    import rasterio
    NC_PATH = "/mnt/ssd1/data/daymet/daymet_v4_daily_na_tmin_2024.nc"
    rast = rasterio.open(f"NETCDF:{NC_PATH}:tmin")
    data = rast.read(183)
    with rasterio.open(TIF_PATH, "w", driver="GTiff",
        height=data.shape[0], width=data.shape[1], count=1,
        dtype=data.dtype, crs=rast.crs,
        transform=rast.transform, nodata=-9999.0, compress="lzw") as dst:
        dst.write(data, 1)
    rast.close()

# Count basins
import geopandas as gpd, rasterio
basins = gpd.read_file(SHP_PATH, bbox=(-125,24,-66,50))
n_basins = len(basins)
del basins; gc.collect()
print(f"{n_basins} CONUS basins", flush=True)

# Build ranges that skip bad indices
ranges = []
for i in range(0, n_basins, CHUNK):
    hi = min(i + CHUNK, n_basins)
    # Check if any bad indices fall in this range
    bads = sorted(b for b in BAD_INDICES if i <= b < hi)
    if not bads:
        ranges.append((i, hi))
    else:
        # Split around bad basins
        prev = i
        for b in bads:
            if b > prev:
                ranges.append((prev, b))
            prev = b + 1
        if prev < hi:
            ranges.append((prev, hi))

CHUNK_SCRIPT = '''
import warnings, gc, time, sys
import geopandas as gpd, rasterio, exactextract

lo, hi = int(sys.argv[1]), int(sys.argv[2])
out = sys.argv[3]

rast = rasterio.open("{tif}")
b = gpd.read_file("{shp}", bbox=(-125,24,-66,50))
b = b.set_crs("EPSG:4326", allow_override=True).to_crs(rast.crs)
chunk = b.iloc[lo:hi].copy()
del b; gc.collect()

t0 = time.time()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    r = exactextract.exact_extract(rast, chunk,
        ops=["mean","count","sum"], include_cols=["COMID"], output="pandas")
dt = time.time() - t0
r.to_csv(out, index=False)
rast.close()
print(f"{{dt:.1f}}")
'''.format(tif=TIF_PATH, shp=SHP_PATH)

# Run chunks via subprocess
total_ee = 0
all_csvs = []
processed = 0

for lo, hi in ranges:
    csv_out = os.path.join(OUT_DIR, f"ee_chunk_{lo}.csv")
    r = subprocess.run(
        [sys.executable, "-c", CHUNK_SCRIPT, str(lo), str(hi), csv_out],
        capture_output=True, text=True, timeout=300
    )
    n = hi - lo
    if r.returncode == 0:
        dt = float(r.stdout.strip())
        total_ee += dt
        all_csvs.append(csv_out)
        processed += n
        print(f"  {lo:>6}-{hi:<6} ({n:>5}): {dt:.1f}s", flush=True)
    else:
        print(f"  {lo:>6}-{hi:<6} ({n:>5}): CRASHED", flush=True)

# Merge
dfs = [pd.read_csv(f) for f in all_csvs]
results = pd.concat(dfs, ignore_index=True)
results.to_csv(os.path.join(OUT_DIR, "exactextract_results.csv"), index=False)
for f in all_csvs:
    os.remove(f)

valid = results["mean"].dropna()
print(f"\nexactextract: {total_ee:.1f}s for {processed} basins ({total_ee/processed*1000:.3f}ms/basin)")
print(f"Valid: {len(valid)}/{processed}")
print(f"  Mean tmin: {valid.mean():.4f}, range [{valid.min():.4f}, {valid.max():.4f}]")
print(f"  Skipped: {len(BAD_INDICES)} bad basins")
