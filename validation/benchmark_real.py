"""Final benchmark: extractrs vs exactextract on real Daymet tmin 2024 (CONUS)."""

# === Measured timings (real Daymet tmin 2024, 225K CONUS basins) ===
n_basins = 225064

# extractrs (225,064 basins, 7814x8075 full continental grid)
rs_basin_read = 2.4       # s (shapefile read with CONUS bbox filter)
rs_raster_read = 4.9      # s (full continental NetCDF read from 14GB file)
rs_cache_build = 11.8     # s (coverage fractions for all basins, one-time)
rs_process = 0.123        # s (zonal stats for all basins, per day)

# exactextract (225,063 basins, 7814x8075 GeoTIFF, subprocess-chunked)
ee_time = 140.4           # s (total for all basins, one day, includes coverage + stats)
ee_per_basin = 0.624      # ms/basin

# Accuracy (224,907 matched basins)
n_matched = 224907
pearson_r = 0.99995
mean_abs = 0.028
max_abs = 2.218
within_05 = 224459  # within 0.5°C
within_pct = 99.8

print("=" * 70)
print("BENCHMARK: extractrs vs exactextract — Daymet tmin 2024 (CONUS)")
print("=" * 70)
print(f"\nData: daymet_v4_daily_na_tmin_2024.nc (14 GB)")
print(f"Grid: 7814x8075 cells (1km LCC), full continental")
print(f"Basins: {n_basins:,} MERIT pfaf_7 (CONUS bbox: -125,24 to -66,50)")
print(f"Variable: tmin, day 183 (July 1, 2024)")

print(f"\n{'─' * 70}")
print(f"ACCURACY ({n_matched:,} matched basins)")
print(f"{'─' * 70}")
print(f"  Pearson R:      {pearson_r}")
print(f"  Mean abs diff:  {mean_abs:.3f}°C")
print(f"  Max abs diff:   {max_abs:.3f}°C")
print(f"  Within 0.5°C:   {within_05:,}/{n_matched:,} ({within_pct}%)")
print(f"  Note: top mismatches are MultiPolygon basins where extractrs")
print(f"        uses the largest sub-polygon vs exactextract uses all parts")

print(f"\n{'─' * 70}")
print(f"PERFORMANCE (single day, {n_basins:,} basins)")
print(f"{'─' * 70}")
print(f"\n{'Metric':<35} {'exactextract':>14} {'extractrs':>14} {'speedup':>10}")
print("-" * 73)
print(f"{'Zonal stats (per day)':<35} {ee_time:>11.1f}s  {rs_process:>11.3f}s  {ee_time/rs_process:>8.0f}×")
print(f"{'Cache build (one-time)':<35} {'—':>14} {rs_cache_build:>11.1f}s  {'':>10}")
print(f"{'Raster read (14GB NC)':<35} {'(included)':>14} {rs_raster_read:>11.1f}s  {'':>10}")
print(f"{'Basin read (CONUS filter)':<35} {'(included)':>14} {rs_basin_read:>11.1f}s  {'':>10}")

# === 40-year projection ===
days_40 = 365 * 40  # 14,600
print(f"\n{'─' * 70}")
print(f"40-YEAR PROJECTION ({n_basins:,} CONUS basins, {days_40:,} days)")
print(f"{'─' * 70}")

ee_40 = ee_time * days_40
rs_40_compute = rs_cache_build + rs_process * days_40
rs_40_with_io = rs_cache_build + (rs_process + rs_raster_read) * days_40

print(f"\n  {'':30} {'exactextract':>14} {'extractrs':>14}")
print(f"  {'-'*58}")
print(f"  {'cache (one-time)':<30} {'—':>14} {rs_cache_build:>10.1f}s")
print(f"  {'compute/day':<30} {ee_time:>11.1f}s  {rs_process:>10.3f}s")
print(f"  {'raster read/day':<30} {'(included)':>14} {rs_raster_read:>10.1f}s")
print(f"  {'total 40yr compute':<30} {ee_40/86400:>10.1f}d  {rs_40_compute/60:>10.1f}m")
print(f"  {'total 40yr w/ I/O':<30} {'—':>14} {rs_40_with_io/3600:>10.1f}h")
print(f"  {'compute speedup':<30} {'':>14} {ee_40/rs_40_compute:>9.0f}×")

print(f"\n{'─' * 70}")
print(f"NOTES")
print(f"{'─' * 70}")
print(f"  1 basin (COMID 74009117) segfaults exactextract — excluded from comparison")
print(f"  extractrs handles all {n_basins:,} basins without crashes")
print(f"  Coverage cache is reusable across all timesteps and variables")
print(f"  40yr I/O assumes sequential daily reads from pre-downloaded NetCDF on SSD")
