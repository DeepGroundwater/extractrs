"""Find the basin that causes exactextract to segfault via subprocess isolation."""
import subprocess, sys

def test_range(lo, hi):
    """Test basins [lo, hi) in a subprocess. Returns True if OK."""
    code = f"""
import warnings, gc
import geopandas as gpd, rasterio, exactextract
rast = rasterio.open('/tmp/extractrs_validation/daymet_tmin_full.tif')
b = gpd.read_file('/mnt/ssd1/data/merit/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp', bbox=(-125,24,-66,50))
b = b.set_crs('EPSG:4326', allow_override=True).to_crs(rast.crs)
chunk = b.iloc[{lo}:{hi}].copy()
del b; gc.collect()
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    exactextract.exact_extract(rast, chunk, ops=['mean'], include_cols=['COMID'], output='pandas')
rast.close()
print('OK')
"""
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    return r.returncode == 0

# Binary search
lo, hi = 50000, 75000
print(f"Searching [{lo}, {hi})...")
while hi - lo > 1:
    mid = (lo + hi) // 2
    ok = test_range(lo, mid)
    if ok:
        print(f"  [{lo},{mid}): OK", flush=True)
        lo = mid
    else:
        print(f"  [{lo},{mid}): CRASH", flush=True)
        hi = mid

print(f"\nBad basin at index {lo}")

# Get COMID
import geopandas as gpd, rasterio
rast = rasterio.open('/tmp/extractrs_validation/daymet_tmin_full.tif')
b = gpd.read_file('/mnt/ssd1/data/merit/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp', bbox=(-125,24,-66,50))
b = b.set_crs('EPSG:4326', allow_override=True).to_crs(rast.crs)
row = b.iloc[lo]
print(f"COMID: {row['COMID']}, geom_type: {row.geometry.geom_type}, bounds: {row.geometry.bounds}")
rast.close()
