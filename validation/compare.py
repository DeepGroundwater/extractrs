"""Compare extractrs and exactextract results."""

import pandas as pd
import numpy as np
import os

OUT_DIR = "/tmp/extractrs_validation"

# Load both results
ee = pd.read_csv(os.path.join(OUT_DIR, "exactextract_results.csv"))
rs = pd.read_csv(os.path.join(OUT_DIR, "extractrs_results.csv"))

print(f"exactextract: {len(ee)} basins")
print(f"extractrs:    {len(rs)} basins")

# Merge on COMID
merged = ee.merge(rs, on="COMID", suffixes=("_ee", "_rs"))
print(f"Matched: {len(merged)} basins")

# Compare means
merged["abs_diff"] = (merged["mean_ee"] - merged["mean_rs"]).abs()
merged["rel_diff"] = merged["abs_diff"] / merged["mean_ee"].abs() * 100

print(f"\n=== Comparison ===")
print(f"  Basins compared: {len(merged)}")
print(f"  Max abs diff:    {merged['abs_diff'].max():.6f} mm/day")
print(f"  Mean abs diff:   {merged['abs_diff'].mean():.6f} mm/day")
print(f"  Median abs diff: {merged['abs_diff'].median():.6f} mm/day")
print(f"  Max rel diff:    {merged['rel_diff'].max():.2f}%")
print(f"  Mean rel diff:   {merged['rel_diff'].mean():.2f}%")
print(f"  Median rel diff: {merged['rel_diff'].median():.2f}%")

# Show top mismatches
print(f"\nTop 10 mismatches:")
top = merged.nlargest(10, "abs_diff")[["COMID", "mean_ee", "mean_rs", "abs_diff", "rel_diff"]]
top.columns = ["COMID", "exactextract", "extractrs", "abs_diff", "rel_diff_%"]
print(top.to_string(index=False))

# Correlation
corr = merged["mean_ee"].corr(merged["mean_rs"])
print(f"\nPearson R: {corr:.8f}")

# Pass/fail
tol = 0.5  # mm/day
n_pass = (merged["abs_diff"] < tol).sum()
print(f"\nWithin {tol} mm/day tolerance: {n_pass}/{len(merged)} ({100*n_pass/len(merged):.1f}%)")

if merged["abs_diff"].max() < 1.0:
    print("\nVERDICT: PASS — results are consistent")
else:
    print(f"\nVERDICT: INVESTIGATE — max diff {merged['abs_diff'].max():.4f} exceeds 1.0 mm/day")
