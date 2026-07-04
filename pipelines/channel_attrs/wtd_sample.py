"""Nearest-channel-cell water-table sampling along MERIT flowlines.

Two WTD sources, one corridor sensitivity column:

- **Zell & Sanford 2020** (ZS, 250 m): CONUS transient groundwater model,
  depth-to-water raster.  Sign convention (verified empirically): positive =
  water table below land surface (95.2 % of CONUS cells ≥ 0; small negative
  values are artesian/gaining reaches — clipped to 0 in output).
  CRS: WGS-84 Albers Equal Area (params identical to EPSG:5070 but WGS-84
  datum; no measurable coordinate shift for CONUS).

- **Fan et al. 2013** (Fan, ~1 km): global annual-mean water-table depth.
  Sign convention: values are NEGATIVE metres (−10 = 10 m below surface).
  Convert: ``wtd_below_surface = −WTD``.

- **Corridor sensitivity** (ZS, extractrs): polygon-mean over the 100 m
  corridor from Task 1, cross-checked against the nearest-cell ZS column
  (spearman ρ > 0.8 expected; divergence → surface and stop).

Coarse WTD grids get NO fine buffer: a 250 m cell already spans MERIT's
positional error (spec §3 buffer strategy; RiverATLAS-style native-grid
association).  Lines are densified to ``spacing_m`` vertices and the nearest
cell value is taken per vertex.

Output: ``derived/wtd_channel.parquet``
  Columns: COMID, wtd_channel_zs, wtd_channel_fan, wtd_corridor100_zs,
           wtd_channel_n
  All values are positive below-surface metres after sign conversion.
"""
import time
from typing import Optional

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rioxarray as rxr  # noqa: F401 — registers .rio accessor

from . import paths

# ── ZS raster constants ────────────────────────────────────────────────────────
ZS_TIF = paths.RAW / "zs_wtd/Output_CONUS_trans_dtw/conus_MF6_SS_Unconfined_250_dtw.tif"
FAN_NC = paths.RAW / "fan_wtd/NAMERICA_WTD_annualmean.nc"

# Densification spacing (≤ ZS cell size so no cell is skipped)
SPACING_M = 200.0


# ── core functions ─────────────────────────────────────────────────────────────

def _densify(line, spacing_m: float):
    """Return evenly-spaced points along *line* at intervals of *spacing_m*."""
    n = max(2, int(line.length // spacing_m) + 1)
    return [line.interpolate(d) for d in np.linspace(0.0, line.length, n)]


def sample_along_lines(
    da: xr.DataArray,
    gdf: gpd.GeoDataFrame,
    spacing_m: float = 200.0,
    bed_depths: Optional[dict] = None,
) -> pd.DataFrame:
    """Sample *da* (2-D DataArray with x/y dims) along MERIT flowlines.

    Parameters
    ----------
    da
        2-D DataArray (dims ``y``, ``x``) with a valid CRS set via
        ``da.rio.write_crs()``.  NaN cells are excluded from mean and count.
    gdf
        GeoDataFrame with a ``COMID`` column; reprojected internally to *da*'s
        CRS before sampling.
    spacing_m
        Densification interval (metres in the equal-area CRS).
    bed_depths
        Optional ``{COMID: bankfull_depth_m}`` mapping.  When supplied, adds a
        ``frac_below`` column: fraction of sampled vertices where
        ``wtd_value > bankfull_depth`` (i.e., water table below the bed —
        losing-possible).  NaN for COMIDs absent from the mapping.

    Returns
    -------
    pd.DataFrame
        Columns: COMID, value_mean, n[, frac_below].

    Implementation note
    -------------------
    All reach vertices are batched into a single vectorised ``da.sel(x, y,
    method="nearest")`` call (one xarray lookup instead of one per reach),
    then grouped by COMID for per-reach aggregation.  This gives a ~10–50×
    speedup over per-reach looping at scale.
    """
    gdf = gdf.to_crs(da.rio.crs)
    want_frac = bed_depths is not None

    # Build flattened arrays of (comid, x, y) for all densified vertices.
    comid_arr: list = []
    x_arr: list = []
    y_arr: list = []
    for comid, line in zip(gdf["COMID"], gdf.geometry):
        pts = _densify(line, spacing_m)
        c = int(comid)
        for p in pts:
            comid_arr.append(c)
            x_arr.append(p.x)
            y_arr.append(p.y)

    if not comid_arr:
        cols = {"COMID": gdf["COMID"].astype(int), "value_mean": np.nan, "n": 0}
        if want_frac:
            cols["frac_below"] = np.nan
        return pd.DataFrame(cols)

    # Single vectorised nearest-cell lookup for all vertices.
    xs = xr.DataArray(x_arr, dims="pt")
    ys = xr.DataArray(y_arr, dims="pt")
    raw = da.sel(x=xs, y=ys, method="nearest").values.astype(float)

    # Per-COMID aggregation.
    comid_np = np.array(comid_arr, dtype=int)
    rows = []
    for comid, idx in _group_indices(comid_np):
        vals = raw[idx]
        finite = vals[np.isfinite(vals)]
        row: dict = {
            "COMID": comid,
            "value_mean": float(np.mean(finite)) if len(finite) else np.nan,
            "n": int(len(finite)),
        }
        if want_frac:
            bd = bed_depths.get(comid, np.nan)
            row["frac_below"] = (
                float(np.mean(finite > bd))
                if len(finite) and np.isfinite(bd)
                else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


def _group_indices(comid_arr: np.ndarray):
    """Yield (comid, indices) pairs for each unique COMID in order."""
    # Preserve the original encounter order (same as input GDF order).
    seen: dict = {}
    order: list = []
    for i, c in enumerate(comid_arr):
        if c not in seen:
            seen[c] = []
            order.append(c)
        seen[c].append(i)
    for c in order:
        yield c, np.array(seen[c], dtype=int)


# ── raster helpers ─────────────────────────────────────────────────────────────

def _open_zs() -> xr.DataArray:
    """Open ZS raster lazily, squeeze band, mask nodata → positive-down DTW."""
    da = rxr.open_rasterio(ZS_TIF, cache=False).squeeze("band")
    nd = da.rio.nodata
    return da.where(da != nd)


def _open_fan() -> xr.DataArray:
    """Open Fan 2013 NA annual-mean, apply sign convention, return EPSG:4326 DA.

    Fan stores WTD as negative metres (−10 = 10 m below surface).
    Convert: ``wtd_below_surface = −WTD``.
    """
    import xarray as xr

    ds = xr.open_dataset(FAN_NC)
    fan = ds["WTD"].isel(time=0)
    # Apply sign conversion to positive-below-surface.
    fan = -fan
    # Rename lat/lon → y/x so sample_along_lines can call da.sel(x=..., y=...).
    fan = fan.rename({"lat": "y", "lon": "x"})
    fan = fan.rio.write_crs("EPSG:4326")
    return fan


def _pfaf4(comids: pd.Series) -> pd.Series:
    """First 4 pfaf-level digits from 8-digit MERIT COMID (e.g. 71001234→7100)."""
    return (comids // 10000).astype(int)


def _tile_window(da: xr.DataArray, bounds, buf_m: float = 500.0):
    """Return a view of *da* clipped to *bounds* + buf_m on each side.

    Works for both ascending-y (Fan lat/lon raster) and descending-y (ZS
    projected raster) coordinate layouts.  xarray's ``sel(y=slice(a, b))``
    honours the direction of the coordinate array: descending y requires
    ``slice(high, low)``; ascending y requires ``slice(low, high)``.
    """
    xmin, ymin, xmax, ymax = bounds
    x_slice = slice(xmin - buf_m, xmax + buf_m)
    y_vals = da.y.values
    if len(y_vals) > 1 and float(y_vals[0]) > float(y_vals[-1]):
        # Decreasing y (top-to-bottom raster, e.g. ZS GeoTIFF)
        y_slice = slice(ymax + buf_m, ymin - buf_m)
    else:
        # Increasing y (lat-ordered, e.g. Fan netCDF)
        y_slice = slice(ymin - buf_m, ymax + buf_m)
    return da.sel(x=x_slice, y=y_slice)


# ── corridor-mean (extractrs) ──────────────────────────────────────────────────

def _corridor_mean_zs(zs_da: xr.DataArray, riv: gpd.GeoDataFrame) -> pd.DataFrame:
    """Compute extractrs corridor-mean of ZS raster over 100 m corridors.

    Processes by pfaf-4 tile to keep peak memory bounded.  Corridors are
    reprojected once to the ZS native CRS (WGS-84 Albers, datum shift <1 m).
    """
    import extractrs  # noqa: F401

    corr = gpd.read_parquet(paths.DERIVED / "corridors_100m.parquet")
    corr = corr.to_crs(zs_da.rio.crs)
    riv_native = riv.to_crs(zs_da.rio.crs)

    riv_native["pfaf4"] = _pfaf4(riv_native["COMID"])
    tiles = riv_native.groupby("pfaf4", sort=False)

    frames = []
    for tile_id, tile_riv in tiles:
        tile_comids = set(tile_riv["COMID"])
        tile_corr = corr[corr["COMID"].isin(tile_comids)]
        if tile_corr.empty:
            continue

        bounds = tile_riv.total_bounds
        tile_da = _tile_window(zs_da, bounds, buf_m=500.0)
        if tile_da.size == 0:
            continue

        ds = tile_da.to_dataset(name="wtd")
        try:
            result = ds.extrs.zonal_stats(tile_corr, stat="mean", id_col="COMID")
            df = result["wtd"].to_dataframe().reset_index()[["COMID", "wtd"]]
            frames.append(df)
        except Exception:
            # No corridor polygons overlap this tile window → skip silently.
            pass

    if not frames:
        return pd.DataFrame({"COMID": riv["COMID"].astype(int), "wtd": np.nan})
    return pd.concat(frames, ignore_index=True).rename(columns={"wtd": "wtd_corridor100_zs"})


# ── main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    t0 = time.time()
    print("Reading MERIT flowlines ...")
    riv = gpd.read_file(paths.MERIT_RIV, columns=["COMID", "lengthkm"])
    n_reaches = len(riv)
    print(f"  {n_reaches:,} reaches  (CRS: {riv.crs})")

    # Bankfull depths for Task-7 frac_below (keyed by int COMID).
    print("Reading bankfull depths ...")
    bnk = pd.read_parquet(paths.DERIVED / "bankfull.parquet", columns=["COMID", "bankfull_depth"])
    bed_depths: dict = dict(zip(bnk["COMID"].astype(int), bnk["bankfull_depth"]))
    print(f"  {len(bed_depths):,} reaches with bankfull depth")

    # ── Open rasters ──────────────────────────────────────────────────────────
    print("Opening ZS raster (lazy) ...")
    zs_da = _open_zs()
    print(f"  ZS shape={zs_da.shape}, CRS={zs_da.rio.crs.to_string()[:60]}")

    print("Opening Fan raster ...")
    fan_da = _open_fan()
    print(f"  Fan shape={fan_da.shape}, CRS={fan_da.rio.crs}")

    # ── Empirical ZS convention check ────────────────────────────────────────
    # Sample mid-CONUS window (x≈[−500k, 500k], y≈[1.5M, 2.5M] in ZS Albers)
    zs_sample = zs_da.sel(x=slice(-500_000, 500_000), y=slice(2_500_000, 1_500_000))
    zs_vals = zs_sample.values.ravel().astype(float)
    zs_vals = zs_vals[np.isfinite(zs_vals)]
    pct_nonneg = np.mean(zs_vals >= 0) * 100
    print(
        f"\nZS convention check:  min={zs_vals.min():.3f}  "
        f"p50={np.median(zs_vals):.3f}  max={zs_vals.max():.3f}  "
        f"pct_nonneg={pct_nonneg:.1f}%"
    )
    assert pct_nonneg >= 90.0, (
        f"ZS raster has {100-pct_nonneg:.1f}% negative values — "
        "convention unexpected; check the raster file."
    )
    print("  ZS sign convention: positive-down depth-to-water (as expected)")

    # Prepare per-tile processing
    riv["pfaf4"] = _pfaf4(riv["COMID"])
    tiles = riv.groupby("pfaf4", sort=False)
    n_tiles = riv["pfaf4"].nunique()
    print(f"\nProcessing {n_tiles} pfaf-4 tiles ...")

    zs_rows: list = []
    fan_rows: list = []

    for i, (tile_id, tile_riv) in enumerate(tiles, 1):
        if i % 5 == 0 or i == n_tiles:
            print(f"  tile {i}/{n_tiles}  ({tile_id})  ...")

        # ── ZS nearest-cell (with bed_depths for frac_below) ─────────────────
        tile_riv_zs = tile_riv.to_crs(zs_da.rio.crs)
        bounds_zs = tile_riv_zs.total_bounds
        tile_zs = _tile_window(zs_da, bounds_zs, buf_m=500.0)
        if tile_zs.size > 0:
            df_zs = sample_along_lines(tile_zs, tile_riv, spacing_m=SPACING_M, bed_depths=bed_depths)
            zs_rows.append(df_zs)

        # ── Fan nearest-cell ──────────────────────────────────────────────────
        tile_riv_fan = tile_riv  # already in EPSG:4326
        bounds_fan_4326 = tile_riv_fan.total_bounds  # (xmin, ymin, xmax, ymax) in 4326
        tile_fan = _tile_window(fan_da, bounds_fan_4326, buf_m=1.0)  # buf in degrees
        if tile_fan.size > 0:
            df_fan = sample_along_lines(tile_fan, tile_riv_fan, spacing_m=SPACING_M)
            fan_rows.append(df_fan)

    print("Concatenating tile results ...")
    zs_df = pd.concat(zs_rows, ignore_index=True)
    fan_df = pd.concat(fan_rows, ignore_index=True)

    # ── Apply sign conventions ────────────────────────────────────────────────
    # ZS: positive-down; clip artesian (negative) to 0
    zs_df["wtd_channel_zs"] = zs_df["value_mean"].clip(lower=0.0)
    # ZS frac_below: already computed per-vertex vs bankfull_depth
    if "frac_below" in zs_df.columns:
        zs_df = zs_df.rename(columns={"frac_below": "_frac_below_zs"})

    # Fan: already negated in _open_fan (values now ≥ 0); clip any numerical noise
    fan_df["wtd_channel_fan"] = fan_df["value_mean"].clip(lower=0.0)

    # ── Fan sign assertions ───────────────────────────────────────────────────
    fan_finite = fan_df["wtd_channel_fan"].dropna()
    assert (fan_finite >= 0).all(), "Fan WTD after conversion has negative values"
    fan_p50 = float(fan_finite.median())
    assert 0 <= fan_p50 <= 50, f"Fan WTD p50={fan_p50:.2f} m outside plausible 0–50 m band"
    print(f"Fan WTD post-conversion: p50={fan_p50:.2f} m  (assert ≥0 and 0–50 m: OK)")

    # ── Corridor-mean (extractrs) over ZS ────────────────────────────────────
    print("Computing corridor-mean (extractrs) over ZS raster ...")
    corr_df = _corridor_mean_zs(zs_da, riv)
    # Apply same ZS sign/clip: positive-down, clip artesian to 0
    if "wtd_corridor100_zs" in corr_df.columns:
        corr_df["wtd_corridor100_zs"] = corr_df["wtd_corridor100_zs"].clip(lower=0.0)

    # ── Merge ─────────────────────────────────────────────────────────────────
    out = (
        riv[["COMID"]].astype({"COMID": int})
        .merge(
            zs_df[["COMID", "wtd_channel_zs", "n"]].rename(columns={"n": "wtd_channel_n"}),
            on="COMID", how="left"
        )
        .merge(fan_df[["COMID", "wtd_channel_fan"]], on="COMID", how="left")
        .merge(corr_df[["COMID", "wtd_corridor100_zs"]], on="COMID", how="left")
    )

    # Persist frac_below alongside ZS for Task 7 (written to the same parquet
    # so wtd_bedrel.py can read it without re-running the raster sampling).
    if "_frac_below_zs" in zs_df.columns:
        out = out.merge(
            zs_df[["COMID", "_frac_below_zs"]].rename(columns={"_frac_below_zs": "frac_below_zs"}),
            on="COMID", how="left"
        )

    # ── Cross-checks ──────────────────────────────────────────────────────────
    from scipy.stats import spearmanr

    both = out.dropna(subset=["wtd_channel_zs", "wtd_channel_fan"])
    rho_zs_fan, _ = spearmanr(both["wtd_channel_zs"], both["wtd_channel_fan"])
    print(f"\nCross-check  spearman(ZS, Fan) = {rho_zs_fan:.3f}  (expect > 0.4)")
    if rho_zs_fan <= 0.4:
        print("  WARNING: spearman(ZS, Fan) below 0.4 — buffer strategy may matter")

    both2 = out.dropna(subset=["wtd_channel_zs", "wtd_corridor100_zs"])
    rho_near_corr, _ = spearmanr(both2["wtd_channel_zs"], both2["wtd_corridor100_zs"])
    print(
        f"Cross-check  spearman(nearest-cell ZS, corridor-mean ZS) = {rho_near_corr:.3f}"
        "  (expect > 0.8)"
    )
    if rho_near_corr <= 0.8:
        print(
            "  WARNING: spearman(nearest-cell, corridor-mean) below 0.8 — "
            "buffer strategy diverges more than literature suggests; REVIEW before proceeding."
        )

    # ── Write output ──────────────────────────────────────────────────────────
    out_path = paths.DERIVED / "wtd_channel.parquet"
    out.to_parquet(out_path, index=False)

    n_zs = out["wtd_channel_zs"].notna().sum()
    n_fan = out["wtd_channel_fan"].notna().sum()
    n_corr = out["wtd_corridor100_zs"].notna().sum()
    elapsed = time.time() - t0
    print(
        f"\nWrote {out_path}  ({len(out):,} rows, {out_path.stat().st_size/1e6:.1f} MB)"
        f"\n  wtd_channel_zs:      {n_zs:,} / {n_reaches:,} reaches"
        f"\n  wtd_channel_fan:     {n_fan:,} / {n_reaches:,} reaches"
        f"\n  wtd_corridor100_zs:  {n_corr:,} / {n_reaches:,} reaches"
        f"\n  spearman(ZS, Fan)={rho_zs_fan:.3f}  spearman(nearest, corridor)={rho_near_corr:.3f}"
        f"\n  Elapsed: {elapsed:.0f}s"
    )


if __name__ == "__main__":
    main()
