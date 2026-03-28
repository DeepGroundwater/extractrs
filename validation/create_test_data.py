"""Create synthetic test raster + subset of MERIT basins for validation."""

import numpy as np
import geopandas as gpd
from pyproj import Transformer
import netCDF4 as nc
import os

SHAPEFILE = "/mnt/ssd1/data/merit/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"
OUT_DIR = "/tmp/extractrs_validation"
os.makedirs(OUT_DIR, exist_ok=True)

# Test region: small area in North Georgia
BBOX_WGS84 = (-84.5, 34.5, -84.0, 35.0)

# Daymet LCC projection
DAYMET_LCC = "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
transformer = Transformer.from_crs("EPSG:4326", DAYMET_LCC, always_xy=True)

# Transform bbox to LCC to set up grid
corners = [
    (BBOX_WGS84[0], BBOX_WGS84[1]),  # SW
    (BBOX_WGS84[2], BBOX_WGS84[1]),  # SE
    (BBOX_WGS84[2], BBOX_WGS84[3]),  # NE
    (BBOX_WGS84[0], BBOX_WGS84[3]),  # NW
]
lcc_corners = [transformer.transform(lon, lat) for lon, lat in corners]
lcc_xmin = min(c[0] for c in lcc_corners)
lcc_xmax = max(c[0] for c in lcc_corners)
lcc_ymin = min(c[1] for c in lcc_corners)
lcc_ymax = max(c[1] for c in lcc_corners)

# Create 1km grid
dx = dy = 1000.0
# Snap to grid
x0 = np.floor(lcc_xmin / dx) * dx - dx
x1 = np.ceil(lcc_xmax / dx) * dx + dx
y0 = np.floor(lcc_ymin / dy) * dy - dy
y1 = np.ceil(lcc_ymax / dy) * dy + dy

x_centers = np.arange(x0 + dx/2, x1, dx)
y_centers = np.arange(y1 - dy/2, y0, -dy)  # descending (north to south)

nx = len(x_centers)
ny = len(y_centers)

print(f"Grid: {nx} x {ny} cells")
print(f"LCC extent: x=[{x0:.0f}, {x1:.0f}] y=[{y0:.0f}, {y1:.0f}]")
print(f"Cell size: {dx}m")

# Create synthetic precipitation: gradient + noise
np.random.seed(42)
rng = np.random.default_rng(42)
# Base gradient: increases with latitude (y-index decreases = more precip)
base = np.linspace(5.0, 15.0, ny).reshape(-1, 1) * np.ones((1, nx))
noise = rng.normal(0, 2, (ny, nx))
prcp = np.maximum(0, base + noise).astype(np.float32)

# Write NetCDF
nc_path = os.path.join(OUT_DIR, "test_prcp.nc")
with nc.Dataset(nc_path, "w", format="NETCDF4") as ds:
    ds.createDimension("x", nx)
    ds.createDimension("y", ny)
    ds.createDimension("time", 1)

    x_var = ds.createVariable("x", "f8", ("x",))
    x_var[:] = x_centers
    x_var.units = "m"
    x_var.standard_name = "projection_x_coordinate"

    y_var = ds.createVariable("y", "f8", ("y",))
    y_var[:] = y_centers
    y_var.units = "m"
    y_var.standard_name = "projection_y_coordinate"

    t_var = ds.createVariable("time", "f4", ("time",))
    t_var[:] = [0]

    prcp_var = ds.createVariable("prcp", "f4", ("time", "y", "x"), fill_value=-9999.0)
    prcp_var[0, :, :] = prcp
    prcp_var.units = "mm/day"

    # Add CRS variable for rioxarray
    crs_var = ds.createVariable("lambert_conformal_conic", "i4")
    crs_var.grid_mapping_name = "lambert_conformal_conic"
    crs_var.latitude_of_projection_origin = 42.5
    crs_var.longitude_of_central_meridian = -100.0
    crs_var.standard_parallel = [25.0, 60.0]
    crs_var.false_easting = 0.0
    crs_var.false_northing = 0.0
    crs_var.crs_wkt = f'PROJCRS["Daymet LCC",BASEGEOGCRS["WGS 84",DATUM["World Geodetic System 1984",ELLIPSOID["WGS 84",6378137,298.257223563]],UNIT["degree",0.0174532925199433]],CONVERSION["Lambert Conformal Conic",METHOD["Lambert Conic Conformal (2SP)"],PARAMETER["Latitude of false origin",42.5,UNIT["degree",0.0174532925199433]],PARAMETER["Longitude of false origin",-100,UNIT["degree",0.0174532925199433]],PARAMETER["Latitude of 1st standard parallel",25,UNIT["degree",0.0174532925199433]],PARAMETER["Latitude of 2nd standard parallel",60,UNIT["degree",0.0174532925199433]],PARAMETER["Easting at false origin",0,UNIT["metre",1]],PARAMETER["Northing at false origin",0,UNIT["metre",1]]],CS[Cartesian,2],AXIS["easting",east,UNIT["metre",1]],AXIS["northing",north,UNIT["metre",1]]]'

    prcp_var.grid_mapping = "lambert_conformal_conic"

print(f"Wrote {nc_path}: {nx}x{ny} grid, prcp range [{prcp.min():.2f}, {prcp.max():.2f}]")

# Read basins subset
print(f"\nReading basins from {SHAPEFILE}...")
basins = gpd.read_file(SHAPEFILE, bbox=BBOX_WGS84)
print(f"  Found {len(basins)} basins in bbox")

# Limit to first 100 for testing
basins = basins.head(100)

# Save as GeoJSON for easier comparison
geojson_path = os.path.join(OUT_DIR, "test_basins.geojson")
basins.to_file(geojson_path, driver="GeoJSON")
print(f"  Saved {len(basins)} basins to {geojson_path}")

# Save COMID list
comid_path = os.path.join(OUT_DIR, "test_comids.txt")
basins["COMID"].to_csv(comid_path, index=False, header=False)
print(f"  COMIDs: {basins['COMID'].tolist()[:10]}...")

# Also write a GeoTIFF of the test data for exactextract
try:
    import rasterio
    from rasterio.transform import from_bounds
    from rasterio.crs import CRS

    tif_path = os.path.join(OUT_DIR, "test_prcp.tif")
    transform = from_bounds(x0, y0, x1, y1, nx, ny)
    crs = CRS.from_proj4(DAYMET_LCC)

    with rasterio.open(
        tif_path, "w", driver="GTiff",
        height=ny, width=nx, count=1,
        dtype="float32", crs=crs,
        transform=transform, nodata=-9999.0,
    ) as dst:
        dst.write(prcp, 1)
    print(f"\nWrote {tif_path}")
except Exception as e:
    print(f"Warning: could not write GeoTIFF: {e}")
    tif_path = None

print("\nTest data created successfully.")
