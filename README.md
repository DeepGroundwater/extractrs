# extractrs

```
                          .                                    .                    
                        .o8                                  .o8                    
 .ooooo.  oooo    ooo .o888oo oooo d8b  .oooo.    .ooooo.  .o888oo oooo d8b  .oooo.o
d88' `88b  `88b..8P'    888   `888""8P `P  )88b  d88' `"Y8   888   `888""8P d88(  "8
888ooo888    Y888'      888    888      .oP"888  888         888    888     `"Y88b. 
888    .o  .o8"'88b     888 .  888     d8(  888  888   .o8   888 .  888     o.  )88b
`Y8bod8P' o88'   888o   "888" d888b    `Y888""8o `Y8bod8P'   "888" d888b    8""888P'
```


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
