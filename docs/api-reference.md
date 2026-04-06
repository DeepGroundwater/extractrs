# API Reference

## xarray Accessors

The primary interface. Registered automatically when you `import extractrs`.

Use the accessor for quick, one-off computations. For multi-variable or multi-file workflows where the grid and polygons stay the same, use the [Low-Level API](#low-level-api) to build the cache once and reuse it.

### Dataset Accessor

::: extractrs.ExtractrsDatasetAccessor
    options:
      show_source: false
      members:
        - zonal_stats

### DataArray Accessor

::: extractrs.ExtractrsDataArrayAccessor
    options:
      show_source: false
      members:
        - zonal_stats

## Low-Level API

For workflows that process many timesteps or multiple variables against the same geometries, use the low-level API to build the cache once and apply statistics repeatedly. This avoids redundant cache builds — the accessor rebuilds the cache on every call.

The functions come in four variants:

| | Single timestep (2D array) | Batch (3D array) |
|---|---|---|
| **Unweighted** | `apply_stat` | `apply_stat_batch` |
| **Weighted** | `apply_stat_weights` | `apply_stat_batch_weights` |

### Example

```python
import numpy as np
import geopandas as gpd
import xarray as xr
from extractrs import build_cache, apply_stat_batch

ds = xr.open_dataset("daymet_tmin_2024.nc")
basins = gpd.read_file("basins.shp")

# Prepare inputs for build_cache
wkb_list = [geom.wkb for geom in basins.geometry]
id_list = basins["COMID"].astype(np.int64).tolist()

# Compute grid extent from cell centers (add half-cell padding for cell edges)
x = ds["x"].values.astype(np.float64)
y = ds["y"].values.astype(np.float64)
dx = abs(float(x[1] - x[0]))
dy = abs(float(y[1] - y[0]))
xmin = float(x.min()) - dx / 2
xmax = float(x.max()) + dx / 2
ymin = float(y.min()) - dy / 2
ymax = float(y.max()) + dy / 2

# Build cache once
cache = build_cache(wkb_list, id_list, xmin, ymin, xmax, ymax, dx, dy)

# Apply to multiple variables using the same cache
tmin = ds["tmin"].values.astype(np.float64)  # shape: (time, y, x)
tmax = ds["tmax"].values.astype(np.float64)

nodata = -9999.0
tmin_means = apply_stat_batch(cache, tmin, nodata, "mean")  # shape: (time, n_basins)
tmax_means = apply_stat_batch(cache, tmax, nodata, "mean")
```

!!! note "Grid extent is cell edges, not cell centers"

    `xmin, ymin, xmax, ymax` passed to `build_cache` define the grid extent at **cell edges**. If your coordinates are cell centers (as is typical in xarray), add half a cell width: `xmin = x_centers.min() - dx / 2`. The xarray accessor does this automatically.

### build_cache

```python
build_cache(
    wkb_list: list[bytes],
    id_list: list[int],
    xmin: float,
    ymin: float,
    xmax: float,
    ymax: float,
    dx: float,
    dy: float,
    method: str = "exact",
) -> CoverageCache
```

Build a coverage cache from WKB geometries and grid parameters.

**Parameters:**

- **wkb_list** — List of WKB-encoded geometry bytes (one per polygon). Obtain from a GeoDataFrame with `[geom.wkb for geom in gdf.geometry]`. Supports Polygon and MultiPolygon types.
- **id_list** — List of integer IDs, same length as `wkb_list`.
- **xmin, ymin, xmax, ymax** — Grid extent at cell edges (not cell centers).
- **dx, dy** — Cell size in the x and y directions.
- **method** — `"exact"` (default) for fractional pixel coverage, `"center"` for cell-center point-in-polygon.

**Returns:** A `CoverageCache` object.

**Raises:** `ValueError` if `wkb_list` and `id_list` have different lengths, or if `method` is not recognized.

### apply_stat

```python
apply_stat(
    cache: CoverageCache,
    raster: numpy.ndarray,  # shape: (rows, cols), dtype: float64
    nodata: float,
    stat: str,
) -> numpy.ndarray  # shape: (n_basins,), dtype: float64
```

Apply a statistic to a single 2D raster using the coverage cache.

**Parameters:**

- **cache** — `CoverageCache` from `build_cache()`.
- **raster** — 2D numpy array of raster values. Must be `float64`.
- **nodata** — Nodata value to exclude from statistics. NaN values are always excluded regardless.
- **stat** — Statistic name (see [Supported Statistics](#supported-statistics)).

**Returns:** 1D numpy array of length `n_basins` with one result per zone.

### apply_stat_batch

```python
apply_stat_batch(
    cache: CoverageCache,
    raster: numpy.ndarray,  # shape: (time, rows, cols), dtype: float64
    nodata: float,
    stat: str,
) -> numpy.ndarray  # shape: (time, n_basins), dtype: float64
```

Apply a statistic to a 3D raster (time, rows, cols) using the coverage cache. Processes all timesteps in a single call.

**Parameters:**

- **cache** — `CoverageCache` from `build_cache()`.
- **raster** — 3D numpy array. Must be `float64`.
- **nodata** — Nodata value to exclude.
- **stat** — Statistic name.

**Returns:** 2D numpy array of shape `(time, n_basins)`.

### apply_stat_weights

```python
apply_stat_weights(
    cache: CoverageCache,
    raster: numpy.ndarray,  # shape: (rows, cols), dtype: float64
    weights: numpy.ndarray,  # shape: (rows, cols), dtype: float64
    nodata: float,
    weight_nodata: float,
    stat: str,
) -> numpy.ndarray  # shape: (n_basins,), dtype: float64
```

Apply a weighted statistic to a single 2D raster. Use this for `weighted_mean` or `weighted_sum` with an external weight raster (e.g., cell-area weights for geographic grids).

**Parameters:**

- **cache** — `CoverageCache` from `build_cache()`.
- **raster** — 2D numpy array of raster values. Must be `float64`.
- **weights** — 2D numpy array of weight values, same shape as `raster`. Must be `float64`.
- **nodata** — Nodata value for the raster.
- **weight_nodata** — Nodata value for the weights array.
- **stat** — Statistic name.

**Returns:** 1D numpy array of length `n_basins`.

### apply_stat_batch_weights

```python
apply_stat_batch_weights(
    cache: CoverageCache,
    raster: numpy.ndarray,  # shape: (time, rows, cols), dtype: float64
    weights: numpy.ndarray,  # shape: (rows, cols), dtype: float64
    nodata: float,
    weight_nodata: float,
    stat: str,
) -> numpy.ndarray  # shape: (time, n_basins), dtype: float64
```

Apply a weighted statistic to a 3D raster with 2D weights. The same weight array is applied to every timestep.

**Parameters:**

- **cache** — `CoverageCache` from `build_cache()`.
- **raster** — 3D numpy array. Must be `float64`.
- **weights** — 2D numpy array, same spatial shape as each raster slice. Must be `float64`.
- **nodata** — Nodata value for the raster.
- **weight_nodata** — Nodata value for the weights.
- **stat** — Statistic name.

**Returns:** 2D numpy array of shape `(time, n_basins)`.

### CoverageCache

An opaque object returned by `build_cache`. Pass it to the `apply_stat*` functions. It stores the precomputed coverage fractions for each polygon against the grid.

- Tied to a specific grid geometry — reuse it only with rasters on the same grid
- Cannot be serialized to disk or modified after creation
- Not thread-safe

**Properties:**

- **`n_basins`** (`int`) — Number of basins in the cache. May be less than the number of input geometries if some polygons fall outside the grid extent or have empty geometry.

**Methods:**

- **`ids()`** → `numpy.ndarray[int64]` — Basin IDs in cache order.

## Supported Statistics

| Stat | Formula | Notes |
|------|---------|-------|
| `mean` | `sum(x * c) / sum(c)` | Coverage-weighted mean |
| `sum` | `sum(x * c)` | Coverage-weighted sum |
| `count` | `sum(c)` | Sum of coverage fractions (returns a float, not an integer cell count) |
| `min` | `min(x)` where `c > 0` | Minimum value of any covered cell |
| `max` | `max(x)` where `c > 0` | Maximum value of any covered cell |
| `variance` | coverage-weighted | Population variance (divides by total weight, not N-1) |
| `stdev` | coverage-weighted | Population standard deviation |
| `weighted_mean` | `sum(x * c * w) / sum(c * w)` | Requires `apply_stat_weights` or `apply_stat_batch_weights` |
| `weighted_sum` | `sum(x * c * w)` | Requires `apply_stat_weights` or `apply_stat_batch_weights` |

Where `x` = cell value, `c` = coverage fraction, `w` = external weight.

Cells with nodata or NaN values are excluded from all statistics. A zone where all covered cells are nodata returns NaN.

`weighted_mean` and `weighted_sum` are only available through the low-level API (`apply_stat_weights`, `apply_stat_batch_weights`). Passing these stat names to the xarray accessor will run with all weights set to 1.0, producing results identical to `mean` and `sum`.
