# extractrs

Fast exact zonal statistics for xarray — powered by Rust.

extractrs computes coverage-weighted zonal statistics over raster grids using polygon geometries. It provides an alternative to exactextract for common zonal statistics, with a cache-then-apply architecture that eliminates redundant polygon-cell intersection work across timesteps.

## Install

```bash
pip install extractrs
```

For CRS mismatch detection, install the `rio` extra (adds [rioxarray](https://corteva.github.io/rioxarray/)):

```bash
pip install extractrs[rio]
```

## Quick Example

```python
import xarray as xr
import geopandas as gpd
import extractrs  # registers the .extrs xarray accessor

ds = xr.open_dataset("temperature.nc")
basins = gpd.read_file("basins.shp")

result = ds["tmin"].extrs.zonal_stats(basins, stat="mean", id_col="COMID")
```

`result` is an `xarray.DataArray` with the spatial dimensions (y, x) replaced by the zone ID dimension (`COMID`). If the input has a time dimension, the output shape is `(time, COMID)`.

`stat="mean"` computes a coverage-weighted mean — cells partially covered by a polygon contribute proportionally to the result. See [Coverage Methods](#coverage-methods) for details.

## How It Works

1. **Build a coverage cache** — For each polygon, extractrs computes the fractional overlap with every raster cell it touches. This is done once and reused across all timesteps.
2. **Apply a statistic** — For each timestep, the cached coverage fractions are combined with raster values to produce one result per zone.

The cache is the expensive step (~8s for 225K basins). Once built, applying statistics costs ~166ms per timestep. For a single-timestep workflow, total time is ~8.3s vs. ~138s for exactextract (~17x speedup). The speedup grows with more timesteps as the cache cost is amortized.

## Coverage Methods

extractrs supports two methods for computing cell coverage:

| Method | Description | Use case |
|--------|-------------|----------|
| `"exact"` (default) | Fractional pixel coverage via polygon-cell intersection | Production accuracy — uses the same intersection algorithm as exactextract |
| `"center"` | Binary point-in-polygon test on cell centers (1.0 or 0.0) | Fast approximation — comparable to cell-center-based tools like xarray-spatial |

With `"center"`, a polygon smaller than a single grid cell may have zero cell centers inside it, producing no coverage and a NaN result.

## Supported Statistics

| Statistic | Formula | Description |
|-----------|---------|-------------|
| `mean` | `sum(value * coverage) / sum(coverage)` | Coverage-weighted mean |
| `sum` | `sum(value * coverage)` | Coverage-weighted sum (a cell with value 10 and 50% coverage contributes 5.0) |
| `count` | `sum(coverage)` | Sum of coverage fractions for cells with valid data (returns a float, not an integer cell count) |
| `min` | `min(value)` | Minimum value of any touched cell |
| `max` | `max(value)` | Maximum value of any touched cell |
| `variance` | coverage-weighted | Population variance (denominator is total coverage weight, not N-1) |
| `stdev` | coverage-weighted | Population standard deviation |
| `weighted_mean` | `sum(value * coverage * weight) / sum(coverage * weight)` | Mean weighted by an external weight raster (**low-level API only** — see [API Reference](api-reference.md)) |
| `weighted_sum` | `sum(value * coverage * weight)` | Sum weighted by an external weight raster (**low-level API only**) |

## Requirements

- Python >= 3.9
- numpy, xarray, geopandas
- Raster must be on a **regular (equispaced) grid** — curvilinear, rotated, or Gaussian grids are not supported
- Spatial dimensions must be named `y`/`lat`/`latitude` and `x`/`lon`/`longitude` (case-insensitive)
- The `id_col` column must contain **integer values** — string identifiers (e.g., HUC codes) are not supported
- Geometries and raster must share the same CRS — see [CRS Requirements](getting-started.md#crs-requirements)
- extractrs loads data eagerly into memory — Dask-backed arrays will be fully materialized when `zonal_stats` is called
