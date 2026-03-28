"""Diagnose the top mismatch basin to understand coverage fraction errors."""

import geopandas as gpd
import rasterio
import exactextract
import numpy as np
import os

OUT_DIR = "/tmp/extractrs_validation"
SHP_PATH = "/mnt/ssd1/data/merit/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"
TIF_PATH = os.path.join(OUT_DIR, "daymet_tmin_subset.tif")
BBOX_WGS84 = (-84.5, 34.5, -84.0, 35.0)
TARGET_COMID = 73013825

basins = gpd.read_file(SHP_PATH, bbox=BBOX_WGS84)
basins = basins.set_crs("EPSG:4326", allow_override=True)
rast = rasterio.open(TIF_PATH)
basins_proj = basins.to_crs(rast.crs)

target = basins_proj[basins_proj["COMID"] == TARGET_COMID]
geom = target.geometry.iloc[0]
print(f"COMID {TARGET_COMID}:")
print(f"  Geometry type: {geom.geom_type}")
if geom.geom_type == "MultiPolygon":
    for i, poly in enumerate(geom.geoms):
        print(f"  Part {i}: {len(poly.exterior.coords)} vertices, area={poly.area:.0f} m²")
    total_verts = sum(len(p.exterior.coords) for p in geom.geoms)
else:
    total_verts = len(geom.exterior.coords)
    print(f"  Vertices: {total_verts}")
print(f"  Total area: {geom.area:.0f} m² ({geom.area/1e6:.1f} km²)")
print(f"  Bounds (LCC): {target.total_bounds}")

results = exactextract.exact_extract(
    rast, target,
    ops=["mean", "count", "sum", "min", "max"],
    include_cols=["COMID"],
    output="pandas",
)
print(f"\n  exactextract:")
print(f"    mean:  {results['mean'].iloc[0]:.6f}")
print(f"    count: {results['count'].iloc[0]:.1f} effective cells")

bounds = target.total_bounds
n_cells_x = int(np.ceil((bounds[2] - bounds[0]) / 1000))
n_cells_y = int(np.ceil((bounds[3] - bounds[1]) / 1000))
print(f"\n  Bounding box spans: {n_cells_x} x {n_cells_y} = {n_cells_x * n_cells_y} cells")
print(f"  Coverage fraction density: {results['count'].iloc[0]:.1f}/{n_cells_x * n_cells_y} = {results['count'].iloc[0]/(n_cells_x*n_cells_y)*100:.0f}%")
print(f"  → Many boundary cells where coverage fraction accuracy matters")

rast.close()
