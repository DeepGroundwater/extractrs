"""extractrs: Fast zonal statistics for xarray — powered by Rust.

Two coverage methods are available:

- ``method="exact"`` (default): Computes exact fractional pixel coverage
  at polygon boundaries. Matches exactextract to floating-point precision.

- ``method="center"``: Binary cell-center point-in-polygon test (1.0 if
  the cell center is inside the polygon, 0.0 otherwise). Produces results
  identical to xarray-spatial's zonal statistics.
"""

from extractrs._extractrs import CoverageCache, build_cache, apply_stat, apply_stat_batch
from extractrs.accessor import ExtractrsDatasetAccessor, ExtractrsDataArrayAccessor

__all__ = [
    "CoverageCache",
    "build_cache",
    "apply_stat",
    "apply_stat_batch",
]
