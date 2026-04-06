# Getting Started

This guide walks through extractrs from installation to saving results.

## Installation

Install extractrs with pip or uv:

=== "pip"

    ```bash
    pip install extractrs
    ```

=== "uv"

    ```bash
    uv add extractrs
    ```

For CRS mismatch detection, install the `rio` extra (adds [rioxarray](https://corteva.github.io/rioxarray/)):

=== "pip"

    ```bash
    pip install extractrs[rio]
    ```

=== "uv"

    ```bash
    uv add extractrs[rio]
    ```

## Loading Data

extractrs works with xarray Datasets and DataArrays paired with geopandas GeoDataFrames. The raster and polygons must share the same coordinate reference system (CRS).

```python
import xarray as xr
import geopandas as gpd
import extractrs  # registers the .extrs xarray accessor

# Load a gridded dataset (e.g., Daymet, ERA5, PRISM)
ds = xr.open_dataset("daymet_tmin_2024.nc")

# Load polygon boundaries (e.g., catchments, counties, watersheds)
basins = gpd.read_file("merit_basins.shp")
```

extractrs works with any data xarray can open, including GRIB via the cfgrib engine.

### Data Requirements

**Raster** — The xarray object must have:

- Two spatial dimensions named `y`/`lat`/`latitude` and `x`/`lon`/`longitude` (case-insensitive). Other names like `south_north`, `rlat`, or `easting` are not recognized — rename them first: `ds = ds.rename({"south_north": "y", "west_east": "x"})`
- Regular spacing (uniform cell size). `dx` and `dy` can differ from each other, but each axis must be uniformly spaced. Curvilinear or Gaussian grids are not supported.
- Optionally a time dimension named `time` or `t` (case-insensitive). Other names like `valid_time` or `day` must be renamed: `ds = ds.rename({"valid_time": "time"})`

**Polygons** — The GeoDataFrame must contain:

- A geometry column with Polygon or MultiPolygon geometries
- Optionally an integer ID column (e.g., `COMID`, `basin_id`). The column must contain values castable to `int64`. String identifiers like HUC codes are not supported — map them to integers first.

## Computing Zonal Statistics

### DataArray (single variable)

```python
result = ds["tmin"].extrs.zonal_stats(basins, stat="mean", id_col="COMID")
```

Returns an `xarray.DataArray` with spatial dimensions replaced by the zone ID dimension.

### Dataset (multiple variables)

```python
result = ds.extrs.zonal_stats(
    basins,
    stat="mean",
    id_col="COMID",
)
```

This processes **all** data variables in `ds` and returns an `xarray.Dataset`. To select specific variables:

```python
# Single variable (var takes a string)
result = ds.extrs.zonal_stats(basins, stat="mean", id_col="COMID", var="tmin")

# Multiple variables (vars takes a list)
result = ds.extrs.zonal_stats(basins, stat="mean", id_col="COMID", vars=["tmin", "tmax"])
```

`var` accepts a single string. `vars` accepts a list of strings. Do not pass a list to `var` — use `vars` instead. Passing both `var` and `vars` is not supported.

### Output Shape

| Input shape | Output shape | Return type |
|-------------|--------------|-------------|
| `(y, x)` | `(COMID,)` | `DataArray` or `Dataset` |
| `(time, y, x)` | `(time, COMID)` | `DataArray` or `Dataset` |

When `id_col` is not specified, the output dimension is named `basin` and indexed from 0. Result values are `float64`.

Polygons that fall entirely outside the raster extent or have empty geometry are silently excluded from the output. The output may have fewer zones than the input GeoDataFrame.

!!! warning "Extra dimensions"

    extractrs expects 2D `(y, x)` or 3D `(time, y, x)` arrays. If your data has additional dimensions (e.g., `level`, `ensemble`), select a single slice first: `ds["tmin"].sel(level=1000)`.

## Coverage Methods

The `method` parameter controls how cell coverage is computed. The default is `"exact"`.

```python
# Exact fractional coverage (default)
result = ds["tmin"].extrs.zonal_stats(basins, stat="mean", method="exact")

# Cell-center point-in-polygon
result = ds["tmin"].extrs.zonal_stats(basins, stat="mean", method="center")
```

`"exact"` computes the true fractional area of each cell covered by the polygon. A cell half-covered by a polygon contributes 50% of its value to the statistic.

`"center"` tests whether each cell's center point falls inside the polygon. Coverage is binary: 1.0 (inside) or 0.0 (outside). This is faster but less accurate — a polygon smaller than a single grid cell may have zero cell centers inside it, producing no coverage and a NaN result.

## Handling NoData

extractrs automatically detects nodata values by checking, in order:

1. `encoding["_FillValue"]` — where xarray stores the original fill value after decoding
2. `rio.nodata` — if rioxarray is installed
3. `attrs["_FillValue"]` or `attrs["missing_value"]`

If none of these are found, extractrs falls back to `-9999.0` as a sentinel. Cells matching the nodata value are excluded from all statistics.

NaN values are **always** treated as missing regardless of the nodata setting. For most xarray-decoded NetCDF files, `_FillValue` has already been replaced with NaN, so the nodata detection is a secondary safeguard.

!!! warning "Nodata override"

    There is currently no parameter to override the nodata value. If your data uses a non-standard sentinel that extractrs does not detect, set `_FillValue` in your data's encoding before calling `zonal_stats`: `da.encoding["_FillValue"] = -32768.0`.

## Memory

extractrs loads the full array eagerly into memory. If you open datasets with Dask chunks (`xr.open_dataset(..., chunks={"time": 30})`), the entire array will be materialized when `zonal_stats` is called.

For large datasets, use the low-level API to process in batches:

```python
from extractrs import build_cache, apply_stat_batch
import numpy as np

# Build cache once (see API Reference for parameter details)
# The accessor does this internally — the low-level API exposes it directly
cache = build_cache(wkb_list, id_list, xmin, ymin, xmax, ymax, dx, dy)

# Process year by year
for year, group in ds["tmin"].groupby("time.year"):
    data = group.values.astype(np.float64)
    result = apply_stat_batch(cache, data, nodata, "mean")
```

## Saving Results

The output is a standard xarray object. Save it with any format xarray supports:

```python
# Zarr (recommended for large outputs)
result.to_zarr("output.zarr")

# NetCDF
result.to_netcdf("output.nc")

# CSV (for small outputs — use reset_index() for a flat table)
result.to_dataframe().reset_index().to_csv("output.csv")
```

## CRS Requirements

The raster and polygon geometries must be in the same CRS. extractrs does not reproject data.

If rioxarray is installed (`pip install extractrs[rio]`), extractrs will:

- **Raise a `ValueError`** if both the raster and GeoDataFrame have a CRS set and they differ
- **Warn** (`UserWarning`) if the GeoDataFrame has no CRS set
- **Warn** (`UserWarning`) if you use a geographic CRS (e.g., EPSG:4326) with statistics other than `weighted_mean` or `weighted_sum`, because cell areas vary by latitude in geographic projections

If rioxarray is **not** installed, no CRS checking is performed and mismatched CRS will produce silently wrong results.

### Geographic CRS and Area Weighting

With a **projected CRS** (e.g., Lambert Conformal Conic, UTM), grid cells have approximately equal area and coverage-weighted statistics are correct.

With a **geographic CRS** (e.g., EPSG:4326), cells near the equator cover more ground area than cells near the poles. An unweighted `stat="mean"` treats all cells equally by grid fraction, which over-weights equatorial cells relative to their true area. For geographic grids, either:

- Reproject your data to an equal-area projection before calling `zonal_stats`, or
- Use `stat="weighted_mean"` with cell-area weights via the low-level API (`apply_stat_weights`)
