"""Combined benchmark: extractrs vs exactextract timing + 40-year projection."""

import pandas as pd
import subprocess
import time
import os

OUT_DIR = "/tmp/extractrs_validation"

# === exactextract timing (already measured) ===
# Mean: 0.0429s per run for 82 basins on 56x64 grid
ee_time = 0.0429

# === extractrs timing (from validate run) ===
# Cache: 0.0065s, Process: 0.000019s
rs_cache = 0.0065
rs_process = 0.000019

print("=" * 60)
print("BENCHMARK: extractrs vs exactextract")
print("=" * 60)
print(f"\nTest configuration:")
print(f"  Basins: 82")
print(f"  Grid: 56x64 cells (1km LCC)")
print(f"  Variable: prcp (synthetic)")

print(f"\n--- Per-invocation timing ---")
print(f"  exactextract:  {ee_time*1000:.1f}ms (includes coverage + stats)")
print(f"  extractrs:")
print(f"    cache build: {rs_cache*1000:.1f}ms (one-time)")
print(f"    process day: {rs_process*1000:.3f}ms (per day, using cached coverage)")
print(f"  Speedup (process only): {ee_time/rs_process:.0f}x")

print(f"\n--- 40-year projection (82 basins, same grid) ---")
days_40 = 365 * 40  # 14,600 days

# exactextract: must recompute coverage each time (no cache)
ee_total_40 = ee_time * days_40
print(f"\n  exactextract (no cache):")
print(f"    {ee_time:.4f}s × {days_40} days = {ee_total_40:.0f}s = {ee_total_40/60:.1f} min")

# extractrs: cache once, reuse
rs_total_40 = rs_cache + rs_process * days_40
print(f"\n  extractrs (cached coverage):")
print(f"    cache: {rs_cache:.4f}s (one-time)")
print(f"    process: {rs_process:.6f}s × {days_40} days = {rs_process*days_40:.3f}s")
print(f"    total: {rs_total_40:.3f}s = {rs_total_40/60:.4f} min")
print(f"    Speedup vs exactextract: {ee_total_40/rs_total_40:.0f}x")

# Scale to 346K basins (full pfaf_7)
print(f"\n--- 40-year projection (346K basins, continental scale) ---")
basin_scale = 346327 / 82
# Assume grid scales to ~7000x6000 for continental coverage
grid_scale = (7000 * 6000) / (56 * 64)

# exactextract scales linearly with basins × grid
ee_total_346k = ee_time * basin_scale * days_40
print(f"\n  exactextract (346K basins, no cache):")
print(f"    {ee_total_346k:.0f}s = {ee_total_346k/3600:.0f} hours = {ee_total_346k/86400:.1f} days")

# extractrs: cache scales with basins, process scales with basins × avg_pixels_per_basin
# avg pixels per basin ≈ 50 at 1km for MERIT catchments
rs_cache_346k = rs_cache * basin_scale  # ~27s for cache
rs_process_346k = rs_process * basin_scale  # per day
rs_total_346k = rs_cache_346k + rs_process_346k * days_40
print(f"\n  extractrs (346K basins, cached coverage):")
print(f"    cache: {rs_cache_346k:.0f}s ({rs_cache_346k/60:.1f} min, one-time)")
print(f"    process/day: {rs_process_346k:.3f}s")
print(f"    process total: {rs_process_346k*days_40:.0f}s = {rs_process_346k*days_40/60:.1f} min")
print(f"    total: {rs_total_346k:.0f}s = {rs_total_346k/3600:.1f} hours")
print(f"    Speedup vs exactextract: {ee_total_346k/rs_total_346k:.0f}x")

# With download overhead
download_per_day = 30  # seconds for NCSS continental subset
download_total = download_per_day * days_40
print(f"\n--- With network download ({download_per_day}s/day) ---")
print(f"  Download total: {download_total/3600:.0f} hours")
print(f"  extractrs + download: {(rs_total_346k + download_total)/3600:.1f} hours ({(rs_total_346k + download_total)/86400:.1f} days)")

# With pre-downloaded annual files
print(f"\n--- With pre-downloaded annual files ---")
print(f"  Storage: ~15 GB/year × 40 years = ~600 GB")
print(f"  extractrs compute only: {rs_total_346k/3600:.1f} hours")
print(f"  This is the recommended approach for production 40-year runs")
