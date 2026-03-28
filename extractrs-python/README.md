# extractrs

Fast exact zonal statistics for xarray — powered by Rust.

## Install

```bash
pip install extractrs
```

## Usage

```python
import xarray as xr
import geopandas as gpd
import extractrs as extrs

ds = xr.open_dataset("raster.nc")
gdf = gpd.read_file("polygons.shp")

result = ds.extrs.zonal_stats(gdf, stat="mean", id_col="COMID")
result.to_zarr("output.zarr")
```
