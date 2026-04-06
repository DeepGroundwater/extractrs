# Benchmarks

!!! note

    All benchmark values are taken from the examples/ folder and the included Jupyter Notebooks. 

All benchmarks compare extractrs against [exactextract](https://github.com/isciences/exactextract), the reference implementation for exact zonal statistics.

## CONUS-Scale Benchmark

**Dataset:** Daymet daily minimum temperature (2024), 1km Lambert Conformal Conic grid, 7814 x 8075 cells (~63M cells).

**Polygons:** 225,064 MERIT Pfafstetter level-7 sub-basins covering the contiguous United States.

### Python Package

Median of 10 runs for extractrs (cache pre-built) vs. median of 3 runs for exactextract (fewer runs due to the longer runtime):

| Step | extractrs | exactextract |
|------|----------:|-------------:|
| Cache build (one-time) | 8,181 ms | — |
| Per-timestep stat | 166 ms | 138,392 ms |
| **Single-timestep total** | **8,347 ms** | **138,392 ms** |

**Single timestep: 17x speedup.** For multi-timestep workflows, the cache cost amortizes away.

**365 daily timesteps (1 year):** extractrs total = 8.2s + 365 x 0.17s = **69s**. exactextract total = 365 x 138s = **14 hours**. Speedup: **733x**.

### Rust Pipeline

| Step | Time |
|------|-----:|
| Basin read (shapefile, filtered to CONUS bounding box) | 2.4 s |
| Raster read (14 GB NetCDF) | 4.9 s |
| Coverage cache build (one-time) | 12.0 s |
| Per-day zonal stats | 0.125 s |
| exactextract per-day equivalent | 140.4 s |

The Python package per-timestep time (166ms) is higher than the Rust pipeline (125ms) due to numpy array marshaling overhead across the PyO3 boundary.

### Why the Cache Matters

exactextract recomputes polygon-cell intersections on every call. extractrs separates the intersection step (cache build) from the aggregation step (apply stat). For multi-timestep workflows — processing a year of daily data across thousands of basins — the intersection geometry never changes. Computing it once and reusing it eliminates redundant work.

## Accuracy

extractrs uses the same polygon-cell intersection algorithm as exactextract. On a validation subset of 82 CONUS basins (spanning diverse sizes from small headwater catchments to large river basins):

- **Pearson R = 1.0**
- **Maximum absolute difference = 0.0°C** (bit-identical results)

At full CONUS scale (225,064 basins), the Rust pipeline achieves:

- **R = 0.9999999999**

The deviation at scale comes from floating-point accumulation order differences between the two implementations, not algorithmic differences. Both implementations compute the same polygon-cell intersections — the only variation is the order in which partial-cell contributions are summed, which affects the least significant bits of the result.

## Environment

- **CPU:** Single core — no parallelism in the stat application step. The cache build is also single-threaded.
- **extractrs version:** 0.1.0
- **exactextract version:** 0.2.0 (PyPI)
- **Python:** 3.13
- **Dependencies:** numpy, xarray, geopandas
