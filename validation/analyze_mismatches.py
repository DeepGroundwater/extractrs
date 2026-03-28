"""Analyze what drives the mismatches between extractrs and exactextract."""

import pandas as pd
import numpy as np
import os

OUT_DIR = "/tmp/extractrs_validation"
ee = pd.read_csv(os.path.join(OUT_DIR, "exactextract_results.csv"))
rs = pd.read_csv(os.path.join(OUT_DIR, "extractrs_results.csv"))
merged = ee.merge(rs, on="COMID", suffixes=("_ee", "_rs"))
merged["abs_diff"] = (merged["mean_ee"] - merged["mean_rs"]).abs()
merged["rel_diff"] = merged["abs_diff"] / merged["mean_ee"].abs() * 100

# Use count from exactextract as a measure of basin size (in pixels)
print("=== Mismatch vs Basin Size ===")
print(f"{'Size bucket':<20} {'n':>5} {'mean_abs_diff':>14} {'mean_rel_diff':>14}")
for lo, hi, label in [(0, 10, "tiny (<10px)"), (10, 30, "small (10-30px)"), (30, 60, "medium (30-60px)"), (60, 200, "large (>60px)")]:
    sub = merged[(merged["count"] >= lo) & (merged["count"] < hi)]
    if len(sub) > 0:
        print(f"{label:<20} {len(sub):>5} {sub['abs_diff'].mean():>14.4f} {sub['rel_diff'].mean():>14.2f}%")

# Correlation between count and abs_diff
corr = merged["count"].corr(merged["abs_diff"])
print(f"\nCorrelation(count, abs_diff): {corr:.4f}")
# Negative correlation means larger basins have smaller errors
# This is expected for boundary-fraction approximation errors

print(f"\n=== Overall Stats ===")
print(f"  Total basins: {len(merged)}")
print(f"  Mean count (pixels): {merged['count'].mean():.1f}")
print(f"  Basins with <0.1 mm diff: {(merged['abs_diff'] < 0.1).sum()}")
print(f"  Basins with <0.5 mm diff: {(merged['abs_diff'] < 0.5).sum()}")
print(f"  Basins with <1.0 mm diff: {(merged['abs_diff'] < 1.0).sum()}")
print(f"  Basins with >1.0 mm diff: {(merged['abs_diff'] >= 1.0).sum()}")
