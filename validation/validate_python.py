"""
Validate extractrs Python package against exactextract.

Runs both tools on the same real Daymet raster + MERIT basins,
compares results, and benchmarks performance.

Usage:
    python validation/validate_python.py
    python validation/validate_python.py --conus        # full CONUS (225K basins)
    python validation/validate_python.py --bbox -85,34,-84,35  # custom region
"""

import argparse
import os
import time
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from rasterio.windows import from_bounds

import exactextract
import extractrs as extrs

# ── Defaults ──
NC_PATH = "/mnt/ssd1/data/daymet/daymet_v4_daily_na_tmin_2024.nc"
SHP_PATH = "/mnt/ssd1/data/merit/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"
OUT_DIR = "/tmp/extractrs_validation"
TIME_IDX = 182  # July 1 (0-indexed)
VAR_NAME = "tmin"

# Known bad basins that segfault exactextract
BAD_COMIDS = {74009117}


def parse_args():
    parser = argparse.ArgumentParser(description="Validate extractrs vs exactextract")
    parser.add_argument("--conus", action="store_true", help="Full CONUS (-125,24,-66,50)")
    parser.add_argument("--bbox", type=str, default="-84.5,34.5,-84.0,35.0",
                        help="WGS84 bbox: west,south,east,north")
    parser.add_argument("--nc", type=str, default=NC_PATH, help="Daymet NetCDF path")
    parser.add_argument("--shp", type=str, default=SHP_PATH, help="MERIT shapefile path")
    parser.add_argument("--time-idx", type=int, default=TIME_IDX, help="Time index (0-based)")
    parser.add_argument("--var", type=str, default=VAR_NAME, help="Variable name")
    parser.add_argument("--n-bench", type=int, default=10, help="Benchmark iterations")
    return parser.parse_args()


def load_basins(shp_path, bbox_wgs84):
    """Load MERIT basins with WGS84 bbox filter."""
    basins = gpd.read_file(shp_path, bbox=bbox_wgs84)
    basins = basins.set_crs("EPSG:4326", allow_override=True)
    # Exclude known bad basins
    basins = basins[~basins["COMID"].isin(BAD_COMIDS)]
    return basins


def run_exactextract(rast_sub, basins_proj, stat="mean", weights=None):
    """Run exactextract on a rasterio dataset + projected basins."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kwargs = dict(
            include_cols=["COMID"],
            output="pandas",
        )
        if weights is not None:
            kwargs["weights"] = weights
        results = exactextract.exact_extract(
            rast_sub, basins_proj,
            ops=[stat],
            **kwargs,
        )
    return results


def run_extractrs(ds, basins_proj, id_col="COMID", stat="mean"):
    """Run extractrs on an xarray Dataset + GeoDataFrame (must be same CRS)."""
    result = ds.extrs.zonal_stats(basins_proj, stat=stat, id_col=id_col)
    return result


def prepare_exactextract_raster(nc_path, var_name, time_idx, basins_proj, out_dir):
    """Create a GeoTIFF subset for exactextract from the NetCDF."""
    subdataset = f"NETCDF:{nc_path}:{var_name}"
    rast = rasterio.open(subdataset)

    # Compute bbox from projected basins with padding
    total_bounds = basins_proj.total_bounds
    pad = 2000  # 2km
    lcc_bbox = (
        total_bounds[0] - pad, total_bounds[1] - pad,
        total_bounds[2] + pad, total_bounds[3] + pad,
    )

    # Snap window to integer pixel boundaries so the GeoTIFF grid
    # aligns exactly with the source NetCDF grid
    window = from_bounds(*lcc_bbox, rast.transform)
    window = window.round_offsets().round_lengths(op="ceil")
    band_idx = time_idx + 1  # rasterio is 1-indexed

    data = rast.read(band_idx, window=window)
    win_transform = rast.window_transform(window)

    tif_path = os.path.join(out_dir, "daymet_subset.tif")
    with rasterio.open(
        tif_path, "w", driver="GTiff",
        height=data.shape[0], width=data.shape[1], count=1,
        dtype=data.dtype, crs=rast.crs,
        transform=win_transform, nodata=-9999.0,
    ) as dst:
        dst.write(data, 1)

    raster_crs = rast.crs
    rast.close()
    return tif_path, raster_crs


def prepare_extractrs_dataset(nc_path, var_name, time_idx):
    """Open a single time slice as an xarray Dataset with rioxarray CRS."""
    ds = xr.open_dataset(nc_path)
    ds_slice = ds[[var_name]].isel(time=time_idx)
    # Attach CRS via rioxarray for auto-reprojection
    import rioxarray  # noqa: F401 — registers .rio accessor
    ds_slice = ds_slice.rio.write_crs(
        "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 "
        "+x_0=0 +y_0=0 +datum=WGS84 +units=m"
    )
    return ds_slice


def compare(ee_df, rs_ds, var_name, stat="mean"):
    """Compare exactextract (DataFrame) vs extractrs (xarray Dataset) results."""
    # Extract extractrs results as DataFrame
    rs_vals = rs_ds[var_name].values
    rs_ids = rs_ds["COMID"].values
    rs_df = pd.DataFrame({"COMID": rs_ids, f"{stat}_rs": rs_vals})

    # Rename exactextract column
    ee_df = ee_df.rename(columns={stat: f"{stat}_ee"})

    # Merge on COMID
    merged = ee_df[["COMID", f"{stat}_ee"]].merge(rs_df, on="COMID")

    # Drop rows where either is NaN
    merged = merged.dropna(subset=[f"{stat}_ee", f"{stat}_rs"])

    merged["abs_diff"] = (merged[f"{stat}_ee"] - merged[f"{stat}_rs"]).abs()
    merged["rel_diff"] = merged["abs_diff"] / merged[f"{stat}_ee"].abs().clip(lower=1e-10) * 100

    return merged


def print_report(merged, n_basins, ee_time, rs_time, rs_cache_time, n_bench,
                 stat="mean", unit="°C"):
    """Print the validation report."""
    ee_col = f"{stat}_ee"
    rs_col = f"{stat}_rs"

    print("\n" + "=" * 70)
    print(f"VALIDATION: extractrs vs exactextract [{stat}]")
    print("=" * 70)

    # Accuracy
    print(f"\n{'─' * 70}")
    print(f"ACCURACY ({len(merged):,} matched basins)")
    print(f"{'─' * 70}")

    corr = merged[ee_col].corr(merged[rs_col])
    print(f"  Pearson R:       {corr:.10f}")
    print(f"  Mean abs diff:   {merged['abs_diff'].mean():.6f}{unit}")
    print(f"  Max abs diff:    {merged['abs_diff'].max():.6f}{unit}")
    print(f"  Median abs diff: {merged['abs_diff'].median():.6f}{unit}")

    for tol in [0.001, 0.01, 0.1, 0.5]:
        n_pass = (merged["abs_diff"] < tol).sum()
        pct = 100 * n_pass / len(merged)
        print(f"  Within {tol}{unit}:  {n_pass:>8,}/{len(merged):,} ({pct:.2f}%)")

    # Top mismatches
    print(f"\n  Top 5 mismatches:")
    top = merged.nlargest(5, "abs_diff")[["COMID", ee_col, rs_col, "abs_diff"]]
    for _, row in top.iterrows():
        print(f"    COMID {int(row['COMID']):>10}: "
              f"ee={row[ee_col]:.4f}  rs={row[rs_col]:.4f}  "
              f"diff={row['abs_diff']:.6f}{unit}")

    # Performance
    print(f"\n{'─' * 70}")
    print(f"PERFORMANCE ({n_basins:,} basins, {n_bench} iterations)")
    print(f"{'─' * 70}")
    print(f"  exactextract:     {ee_time*1000:.1f}ms/run")
    print(f"  extractrs:")
    print(f"    cache build:    {rs_cache_time*1000:.1f}ms (one-time)")
    print(f"    zonal stats:    {rs_time*1000:.1f}ms/run")
    if ee_time > 0:
        print(f"    speedup:        {ee_time/rs_time:.0f}×")

    # Verdict
    print(f"\n{'─' * 70}")
    max_diff = merged["abs_diff"].max()
    if max_diff < 0.1:
        print(f"VERDICT: PASS — max diff {max_diff:.6f}{unit} < 0.1{unit}")
    elif max_diff < 0.5:
        print(f"VERDICT: ACCEPTABLE — max diff {max_diff:.4f}{unit} < 0.5{unit}")
    else:
        print(f"VERDICT: INVESTIGATE — max diff {max_diff:.4f}{unit} >= 0.5{unit}")
    print("=" * 70)

    return merged


def main():
    args = parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    if args.conus:
        bbox = (-125.0, 24.0, -66.0, 50.0)
    else:
        bbox = tuple(float(x) for x in args.bbox.split(","))

    print(f"Config:")
    print(f"  NetCDF:    {args.nc}")
    print(f"  Shapefile: {args.shp}")
    print(f"  Variable:  {args.var}, time_idx={args.time_idx}")
    print(f"  BBox:      {bbox}")
    print(f"  Bench:     {args.n_bench} iterations")

    # ── Load basins ──
    print(f"\nLoading basins...")
    t0 = time.time()
    basins = load_basins(args.shp, bbox)
    print(f"  {len(basins)} basins in {time.time()-t0:.1f}s")

    # ── Prepare extractrs inputs ──
    # Build the xarray dataset first so we can use its CRS for everything
    print(f"\nPreparing extractrs dataset...")
    ds = prepare_extractrs_dataset(args.nc, args.var, args.time_idx)
    raster_crs = ds.rio.crs
    basins_proj = basins.to_crs(raster_crs)

    # ── Prepare exactextract inputs ──
    print(f"\nPreparing exactextract raster...")
    tif_path, _ = prepare_exactextract_raster(
        args.nc, args.var, args.time_idx, basins_proj, OUT_DIR
    )
    rast_sub = rasterio.open(tif_path)
    print(f"  GeoTIFF: {rast_sub.width}x{rast_sub.height}")
    print(f"  Dataset: {dict(ds.dims)}")

    # ── Build shared extractrs cache ──
    print(f"\nBuilding extractrs cache...")
    t0 = time.time()
    cache = extrs.build_cache(
        [geom.wkb for geom in basins_proj.geometry],
        basins_proj["COMID"].astype(np.int64).tolist(),
        *_grid_params_from_ds(ds),
    )
    cache_time = time.time() - t0
    print(f"  Cache build: {cache_time*1000:.1f}ms ({cache.n_basins} basins)")

    data_2d = ds[args.var].values.astype(np.float64)
    nodata = float(ds[args.var].attrs.get("_FillValue", -9999.0))

    # ── Create weight raster (y-coordinate gradient) ──
    # Use row index as weight — provides non-trivial spatial variation so
    # weighted_mean differs from mean, while being trivially reproducible.
    y_vals = ds["y"].values.astype(np.float64)
    weight_2d = np.abs(y_vals[:, np.newaxis]) * np.ones((1, data_2d.shape[1]))
    # Write weight GeoTIFF for exactextract
    weight_tif = os.path.join(OUT_DIR, "weights.tif")
    tif_meta = rasterio.open(tif_path)
    with rasterio.open(
        weight_tif, "w", driver="GTiff",
        height=tif_meta.height, width=tif_meta.width, count=1,
        dtype="float64", crs=tif_meta.crs,
        transform=tif_meta.transform, nodata=-9999.0,
    ) as wdst:
        # Crop weight_2d to match the GeoTIFF subset extent
        from rasterio.windows import from_bounds
        win = from_bounds(
            *tif_meta.bounds, ds.rio.transform()
        ).round_offsets().round_lengths(op="ceil")
        r0, c0 = int(win.row_off), int(win.col_off)
        wdata = weight_2d[r0:r0 + tif_meta.height, c0:c0 + tif_meta.width]
        wdst.write(wdata, 1)
    weight_rast = rasterio.open(weight_tif)
    tif_meta.close()

    # ══════════════════════════════════════════════════════════════════════
    # PASS 1: mean (unweighted) — existing validation
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'═' * 70}")
    print("PASS 1: mean (unweighted)")
    print(f"{'═' * 70}")

    ee_times = []
    ee_result = None
    for i in range(args.n_bench):
        t0 = time.time()
        ee_result = run_exactextract(rast_sub, basins_proj, stat="mean")
        ee_times.append(time.time() - t0)
    ee_time = np.median(ee_times)
    print(f"  exactextract median: {ee_time*1000:.1f}ms")

    rs_times = []
    for i in range(args.n_bench):
        t0 = time.time()
        extrs.apply_stat(cache, data_2d, nodata, "mean")
        rs_times.append(time.time() - t0)
    rs_time = np.median(rs_times)
    print(f"  extractrs median:    {rs_time*1000:.3f}ms")

    rs_result = run_extractrs(ds, basins_proj, stat="mean")
    merged_mean = compare(ee_result, rs_result, args.var, stat="mean")
    print_report(merged_mean, len(basins), ee_time, rs_time, cache_time, args.n_bench,
                 stat="mean")

    # ══════════════════════════════════════════════════════════════════════
    # PASS 2: weighted_mean
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'═' * 70}")
    print("PASS 2: weighted_mean (y-coordinate weights)")
    print(f"{'═' * 70}")

    ee_w_times = []
    ee_w_result = None
    for i in range(args.n_bench):
        t0 = time.time()
        ee_w_result = run_exactextract(
            rast_sub, basins_proj, stat="weighted_mean", weights=weight_rast,
        )
        ee_w_times.append(time.time() - t0)
    ee_w_time = np.median(ee_w_times)
    print(f"  exactextract median: {ee_w_time*1000:.1f}ms")

    rs_w_times = []
    for i in range(args.n_bench):
        t0 = time.time()
        extrs.apply_stat_weights(cache, data_2d, weight_2d, nodata, np.nan, "weighted_mean")
        rs_w_times.append(time.time() - t0)
    rs_w_time = np.median(rs_w_times)
    print(f"  extractrs median:    {rs_w_time*1000:.3f}ms")

    # Build extractrs weighted result as xarray for compare()
    rs_w_vals = extrs.apply_stat_weights(
        cache, data_2d, weight_2d, nodata, np.nan, "weighted_mean",
    )
    id_list = basins_proj["COMID"].astype(np.int64).tolist()
    id_map = dict(zip(
        [e for e in cache.ids()], id_list,
    ))
    rs_w_ids = [id_map.get(i, i) for i in cache.ids()]
    rs_w_ds = xr.Dataset({
        args.var: xr.DataArray(
            np.asarray(rs_w_vals), dims=["COMID"], coords={"COMID": rs_w_ids},
        ),
    })

    merged_wmean = compare(ee_w_result, rs_w_ds, args.var, stat="weighted_mean")
    print_report(merged_wmean, len(basins), ee_w_time, rs_w_time, cache_time, args.n_bench,
                 stat="weighted_mean")

    weight_rast.close()

    # ── Save CSVs ──
    merged_mean.to_csv(os.path.join(OUT_DIR, "comparison_mean.csv"), index=False)
    merged_wmean.to_csv(os.path.join(OUT_DIR, "comparison_weighted_mean.csv"), index=False)
    print(f"\nResults saved to {OUT_DIR}/")

    rast_sub.close()


def _grid_params_from_ds(ds):
    """Extract grid params (xmin, ymin, xmax, ymax, dx, dy) from xarray Dataset."""
    # Find x and y coordinate names
    for name in ds.dims:
        low = name.lower()
        if low in ("y", "lat", "latitude"):
            y_name = name
        elif low in ("x", "lon", "longitude"):
            x_name = name

    x_vals = ds[x_name].values.astype(np.float64)
    y_vals = ds[y_name].values.astype(np.float64)

    dx = abs(float(x_vals[1] - x_vals[0]))
    dy = abs(float(y_vals[1] - y_vals[0]))

    xmin = float(x_vals.min()) - dx / 2
    xmax = float(x_vals.max()) + dx / 2
    ymin = float(y_vals.min()) - dy / 2
    ymax = float(y_vals.max()) + dy / 2

    return xmin, ymin, xmax, ymax, dx, dy


if __name__ == "__main__":
    main()
