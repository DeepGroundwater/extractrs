"""Buffered channel corridors from MERIT flowlines.

Two corridor products (spec: specs/2026-07-04-corridor-buffer-scaling.md):

- ``corridors_100m`` — fixed 100 m half-width: the positional-error floor
  (Amatulli et al. 2022; StreamCat precedent, Hill et al. 2016). Valid for the
  majority of reaches (those below the width hinge) and kept as a sensitivity
  column.
- ``corridors_scaled`` — per-reach ``half_width = max(100 m, 1.5 * bankfull
  width)``. Below ~order 6 the positional error dominates and this equals the
  100 m set; large rivers widen with channel size so the corridor samples the
  channel itself, not a fixed strip across it.

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


def order_bankfull_width_m(order):
    """Fallback bankfull width (m) from Strahler order.

    ``w(omega) = 2 m * 1.9^(omega-1)``: order-1 width ~2 m (Downing et al.
    2012), growing ~1.9x per order (Leopold & Maddock 1953 hydraulic geometry
    composed with Horton area ratios). Used only where observed/modeled width is
    unavailable — order carries +/-1-order scatter, so it is the last resort.
    """
    return paths.ORDER_WIDTH_BASE_M * paths.ORDER_WIDTH_RATIO ** (np.asarray(order, dtype=float) - 1)


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
    """Reproject to the equal-area CRS; a no-op when already there.

    The CRS guard lets ``main()`` reproject the full frame once while keeping
    each ``build_*`` function self-contained (and testable on lon/lat input).
    """
    if riv.crs is None or riv.crs != paths.CRS_EQUAL_AREA:
        riv = riv.to_crs(paths.CRS_EQUAL_AREA)
    return riv


def build_corridors(riv: gpd.GeoDataFrame, half_width_m: float) -> gpd.GeoDataFrame:
    """Fixed-width corridors: buffer every reach by the same half-width."""
    riv = _to_equal_area(riv)
    out = riv[["COMID"]].copy()
    out["geometry"] = riv.geometry.buffer(half_width_m)
    return gpd.GeoDataFrame(out, geometry="geometry", crs=paths.CRS_EQUAL_AREA)


def build_scaled_corridors(riv: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Width-scaled corridors: ``half_width = max(100 m, 1.5 * bankfull width)``."""
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
        choices=["scaled", "100m"],
        default=None,
        help="Build only one corridor set; omit to build both.",
    )
    args = parser.parse_args(argv)

    print(f"Reading MERIT flowlines from {paths.MERIT_RIV} ...")
    t0 = time.time()
    riv = gpd.read_file(paths.MERIT_RIV, columns=["COMID", "order"])
    print(f"  Read {len(riv):,} features in {time.time() - t0:.1f}s  (CRS: {riv.crs})")

    # Merge real-width columns for the scaled-corridor priority chain.
    n_before = len(riv)
    cw = pd.read_parquet(paths.DERIVED / "channel_width_obs.parquet",
                         columns=["COMID", "channel_width_obs"])
    riv = riv.merge(cw, on="COMID", how="left")

    bf = pd.read_parquet(paths.DERIVED / "bankfull.parquet",
                         columns=["COMID", "bankfull_width"])
    riv = riv.merge(bf, on="COMID", how="left")
    assert len(riv) == n_before, (
        f"Merge changed row count: {n_before} → {len(riv)}"
    )
    n_obs = riv["channel_width_obs"].notna().sum()
    n_bf = riv["bankfull_width"].notna().sum()
    print(f"  Width coverage: channel_width_obs={n_obs:,}  bankfull_width={n_bf:,} "
          f"  order-fallback={n_before - max(n_obs, n_bf):,}")

    riv = riv.to_crs(paths.CRS_EQUAL_AREA)  # reproject once; build_* are no-ops below

    builders = [
        ("corridors_100m", lambda r: build_corridors(r, paths.E_POS_M)),
        ("corridors_scaled", build_scaled_corridors),
    ]
    for name, fn in builders:
        if args.only is not None and args.only != name.replace("corridors_", ""):
            print(f"Skipping {name} (--only {args.only})")
            continue
        t1 = time.time()
        print(f"Building {name} ...")
        corr = fn(riv)
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
