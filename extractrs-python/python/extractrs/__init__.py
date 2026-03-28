"""extractrs: Fast exact zonal statistics for xarray — powered by Rust."""

from extractrs._extractrs import CoverageCache, build_cache, apply_stat, apply_stat_batch
from extractrs.accessor import ExtractrsDatasetAccessor, ExtractrsDataArrayAccessor

__all__ = [
    "CoverageCache",
    "build_cache",
    "apply_stat",
    "apply_stat_batch",
]
