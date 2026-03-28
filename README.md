# extractrs

Fast exact zonal statistics for xarray — backend powered by Rust.

## Install

```bash
pip install extractrs
pip install extractrs[rio]  # adds automatic CRS reprojection via rioxarray
```

## Quick start

```python
import xarray as xr
import geopandas as gpd
import extractrs

ds = xr.open_dataset("temperature.nc")
basins = gpd.read_file("basins.shp")

result = ds.extrs.zonal_stats(basins, stat="mean", id_col="COMID")
# result is an xarray Dataset with dims (time, COMID)

result.to_zarr("output.zarr")
```

## API

### `ds.extrs.zonal_stats(...)`

```python
ds.extrs.zonal_stats(
    gdf,                # GeoDataFrame of polygons
    stat="mean",        # statistic to compute
    id_col=None,        # column to use as zone IDs (default: integer index)
    var=None,           # single variable name to process
    vars=None,          # list of variable names to process
)
```

Returns an `xarray.Dataset` with spatial dimensions replaced by the zone ID dimension. When neither `var` nor `vars` is specified, all data variables are processed. The same interface is available on DataArrays via `da.extrs.zonal_stats(gdf, stat, id_col)`.

### Supported statistics

`mean` | `sum` | `count` | `min` | `max` | `variance` | `stdev`

`mean`, `sum`, `variance`, and `stdev` are coverage-weighted.

### Low-level API

See `extractrs.build_cache()` and `extractrs.apply_stat()` for direct cache control when processing many timesteps against the same geometries.
