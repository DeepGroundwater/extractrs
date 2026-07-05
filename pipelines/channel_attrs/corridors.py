"""Buffered channel corridors from MERIT flowlines.

Two corridor products (spec: specs/2026-07-04-corridor-buffer-scaling.md):

- ``corridors_10m`` — fixed 10 m half-width floor: a minimum-width channel
  sample corridor kept as a sensitivity column.
- ``corridors_scaled`` — per-reach ``half_width = max(10 m, 1.5 * bankfull
  width)``. Orders 1–4 floor at 10 m (1.5 × 5.3 m < 10 m); orders 5–10 all
  clear the hinge, widening with channel size.

Width-estimate priority (spec §"Implementation rule"): observed channel width
(SWORD/GRWL) > modeled bankfull (Zarrabi et al. 2025) > order fallback
``w(omega) = 2 * 1.9^(omega-1)`` m. A 3 km cap is applied after priority
resolution to guard residual estuary/lake SWORD widths (paths.WIDTH_CAP_M).
"""
import argparse
import time

import numpy as np
import pandas as pd
import geopandas as gpd

from . import paths


# Module-level vectorized lookup into paths.WRF_HYDRO_BW_M (orders 1-10, clamped).
_bw_lookup = np.vectorize(
    lambda o: paths.WRF_HYDRO_BW_M[int(np.clip(round(o), 1, 10))],
    otypes=[float],
)


def order_bankfull_width_m(order):
    """Bottom channel width (m) from Strahler order via the WRF-Hydro CONUS lookup.

    Source: ``paths.WRF_HYDRO_BW_M`` (NCAR wrf_hydro_functions.py ``Mannings_Bw``,
    LR 7/01/2020). These are trapezoidal channel BOTTOM widths, narrower than
    bankfull surface width. Order is clamped to [1, 10].

    Hinge note: with the 1.5× rule and a 10 m floor, orders 1–4 are pinned
    (1.5 × 5.3 m = 7.95 m < 10 m); orders 5–10 all clear the hinge
    (order 5: 1.5 × 7.4 m = 11.1 m; order 10: 1.5 × 110 m = 165 m).
    """
    return _bw_lookup(np.asarray(order, dtype=float))


def scaled_half_width_m(width_est_m):
    """Corridor half-width (m) = ``max(positional-error floor, 1.5 * width)``."""
    return np.maximum(paths.E_POS_M, paths.ALPHA_HALF * np.asarray(width_est_m, dtype=float))


def _resolve_width_m(riv):
    """Per-reach bankfull-width estimate following the spec priority chain.

    Prefers observed then modeled width columns where present and finite; falls
    back to the order-derived width otherwise (currently every reach, until the
    Task-3 crosswalk supplies ``channel_width_obs`` / ``bankfull_width``).

    A cap of ``paths.WIDTH_CAP_M`` (3 km) is applied after priority resolution
    to guard against SWORD residual estuary/lake widths.
    """
    width = order_bankfull_width_m(riv["order"].to_numpy(dtype=float))
    # Lowest to highest priority, so later (better) sources overwrite earlier.
    for col in ("bankfull_width", "channel_width_obs"):
        if col in riv.columns:
            observed = riv[col].to_numpy(dtype=float)
            valid = np.isfinite(observed) & (observed > 0)
            width = np.where(valid, observed, width)
    # Cap: residual estuary/lake polygons in SWORD produce widths up to ~17 km
    # on CONUS bugfix1; a >3 km value is not a channel and must not drive corridor
    # sizing (would produce ~4.5 km half-width, blanketing entire flood plains).
    width = np.minimum(width, paths.WIDTH_CAP_M)
    return width


def _to_equal_area(riv):
    """Reproject to the equal-area CRS; a no-op when already there."""
    if riv.crs is None:
        raise ValueError(
            "GeoDataFrame has no CRS set; assign one before calling build_corridors "
            "or build_scaled_corridors"
        )
    if riv.crs != paths.CRS_EQUAL_AREA:
        riv = riv.to_crs(paths.CRS_EQUAL_AREA)
    return riv


def _merge_width_columns(riv: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Left-merge observed/modelled width columns; falls back gracefully when parquets absent.

    Only called when building ``corridors_scaled`` — not needed for the fixed-100 m set.
    Missing parquets are not an error: the width priority chain falls back to the
    order estimate for those reaches.
    """
    n_before = len(riv)
    for parquet, col in [
        (paths.DERIVED / "channel_width_obs.parquet", "channel_width_obs"),
        (paths.DERIVED / "bankfull.parquet", "bankfull_width"),
    ]:
        if not parquet.exists():
            print(f"  {parquet.name} not found — {col} falls back to order estimate")
            riv = riv.copy()
            riv[col] = np.nan
        else:
            df = pd.read_parquet(parquet, columns=["COMID", col])
            riv = riv.merge(df, on="COMID", how="left")
        assert len(riv) == n_before, f"Merge on {col} changed row count"

    obs_valid = riv["channel_width_obs"].notna() if "channel_width_obs" in riv.columns else pd.Series(False, index=riv.index)
    bf_valid = riv["bankfull_width"].notna() if "bankfull_width" in riv.columns else pd.Series(False, index=riv.index)
    print(
        f"  Width coverage: channel_width_obs={obs_valid.sum():,}  "
        f"bankfull_width={bf_valid.sum():,}  "
        f"order-fallback={(~(obs_valid | bf_valid)).sum():,}"
    )
    return riv


def build_corridors(riv: gpd.GeoDataFrame, half_width_m: float) -> gpd.GeoDataFrame:
    """Fixed-width corridors: buffer every reach by the same half-width."""
    riv = _to_equal_area(riv)
    out = riv[["COMID"]].copy()
    out["geometry"] = riv.geometry.buffer(half_width_m)
    return gpd.GeoDataFrame(out, geometry="geometry", crs=paths.CRS_EQUAL_AREA)


def build_scaled_corridors(riv: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Width-scaled corridors: ``half_width = max(10 m, 1.5 * bankfull width)``."""
    riv = _to_equal_area(riv)
    half_width = scaled_half_width_m(_resolve_width_m(riv))
    out = riv[["COMID"]].copy()
    out["half_width_m"] = half_width
    out["geometry"] = riv.geometry.buffer(half_width)
    return gpd.GeoDataFrame(out, geometry="geometry", crs=paths.CRS_EQUAL_AREA)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--only",
        choices=["scaled", "10m"],
        default=None,
        help="Build only one corridor set; omit to build both.",
    )
    args = parser.parse_args(argv)

    print(f"Reading MERIT flowlines from {paths.MERIT_RIV} ...")
    t0 = time.time()
    riv = gpd.read_file(paths.MERIT_RIV, columns=["COMID", "order"])
    print(f"  Read {len(riv):,} features in {time.time() - t0:.1f}s  (CRS: {riv.crs})")
    riv = riv.to_crs(paths.CRS_EQUAL_AREA)  # reproject once; build_* are no-ops below

    builders = [
        ("corridors_10m", lambda r: build_corridors(r, paths.E_POS_M)),
        ("corridors_scaled", build_scaled_corridors),
    ]
    riv_with_width = None  # loaded lazily; only needed for corridors_scaled
    for name, fn in builders:
        if args.only is not None and args.only != name.replace("corridors_", ""):
            print(f"Skipping {name} (--only {args.only})")
            continue
        build_riv = riv
        if name == "corridors_scaled":
            if riv_with_width is None:
                riv_with_width = _merge_width_columns(riv)
            build_riv = riv_with_width
        t1 = time.time()
        print(f"Building {name} ...")
        corr = fn(build_riv)
        out_path = paths.DERIVED / f"{name}.parquet"
        corr.to_parquet(out_path)
        elapsed = time.time() - t1
        size_mb = out_path.stat().st_size / 1e6
        print(f"  {name}: {len(corr):,} corridors, CRS={corr.crs.to_epsg()}, "
              f"{size_mb:.1f} MB  ({elapsed:.1f}s)")
        if name == "corridors_scaled":
            hw = corr["half_width_m"]
            print(f"  half_width_m: p50={hw.median():.1f}  p99={hw.quantile(0.99):.1f} "
                  f" max={hw.max():.1f}  count>{paths.E_POS_M:.0f}m="
                  f"{(hw > paths.E_POS_M).sum():,}")


if __name__ == "__main__":
    main()
