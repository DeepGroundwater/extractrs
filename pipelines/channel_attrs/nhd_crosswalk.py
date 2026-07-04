"""NHDPlusV2 -> MERIT crosswalk by buffered-geometry length-weighted matching
(mirrors Wade et al. 2025's MERIT-SWORD construction; spec §A0).

part_len = length of the NHD flowline clipped to the MERIT reach's matching
buffer (CROSSWALK_BUFFER_M = 300 m — the MERIT positional-error envelope).
match_frac = sum(part_len) / merit reach length, a per-COMID quality flag
(low values = headwater divergence between the networks; consumers should
NaN-out transfers below a threshold, default 0.3).

ID-space note (spec §A0): MERIT COMID and NHDPlus COMID are UNRELATED integer
namespaces that share the same column name. The foreign key (NHD COMID) is
stored as `foreign_id` throughout this module to prevent collisions.
"""
from __future__ import annotations

import time
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pyogrio

from . import paths


def build_crosswalk(
    merit: gpd.GeoDataFrame,
    nhd: gpd.GeoDataFrame,
    buffer_m: float,
    top_k: int,
) -> pd.DataFrame:
    """Return a (COMID, foreign_id, part_len, match_frac) crosswalk table.

    Parameters
    ----------
    merit : GeoDataFrame with columns ['COMID'] and linestring geometry.
    nhd   : GeoDataFrame with columns ['nhd_comid'] and linestring geometry.
    buffer_m : half-width of the matching envelope around each MERIT reach.
    top_k    : keep at most this many NHD partners per MERIT reach.
    """
    merit = merit.to_crs(paths.CRS_EQUAL_AREA)
    nhd = nhd.to_crs(paths.CRS_EQUAL_AREA)

    buf = merit.copy()
    buf["merit_len"] = merit.geometry.length
    buf["geometry"] = merit.geometry.buffer(buffer_m)

    joined = gpd.sjoin(
        nhd, buf[["COMID", "merit_len", "geometry"]], how="inner", predicate="intersects"
    )
    if joined.empty:
        return pd.DataFrame(columns=["COMID", "foreign_id", "part_len", "match_frac"])

    # Exact clipped length per (nhd, merit) candidate pair.
    buf_geo = buf.set_index("COMID")["geometry"]
    joined["part_len"] = [
        geom.intersection(buf_geo.loc[c]).length
        for geom, c in zip(joined.geometry, joined["COMID"])
    ]
    joined = joined[joined["part_len"] > 0].copy()

    out = (
        joined.rename(columns={"nhd_comid": "foreign_id"})[
            ["COMID", "foreign_id", "part_len", "merit_len"]
        ]
        .sort_values(["COMID", "part_len"], ascending=[True, False])
        .groupby("COMID", sort=False)
        .head(top_k)
    )

    frac = (
        out.groupby("COMID")
        .apply(
            lambda g: min(1.0, g["part_len"].sum() / g["merit_len"].iloc[0]),
            include_groups=False,
        )
        .rename("match_frac")
    )
    return out.drop(columns="merit_len").merge(frac, on="COMID")


def _merit_bounds(merit: gpd.GeoDataFrame) -> tuple[float, float, float, float]:
    """Return (minx, miny, maxx, maxy) in EPSG:4326 for the full MERIT dataset."""
    b = merit.to_crs("EPSG:4326").total_bounds
    return tuple(b)  # (minx, miny, maxx, maxy)


def main(
    gdb_path: str | Path | None = None,
    out_path: str | Path | None = None,
    tile_n: int = 8,
) -> None:
    """Build NHD→MERIT crosswalk at CONUS scale, tiled for memory efficiency.

    Parameters
    ----------
    gdb_path : path to the extracted NHD FileGDB directory. Defaults to the
               standard staging location.
    out_path : parquet output path. Defaults to paths.DERIVED / "nhd_merit_crosswalk.parquet".
    tile_n   : number of tiles per axis (tile_n x tile_n grid). Default 8 → 64 tiles.
               Use tile_n=1 for a single-tile timing test.
    """
    gdb_path = Path(gdb_path) if gdb_path else (
        Path("/mnt/ssd1/data/channel_attrs/raw/nhdplusv2/extracted")
        / "NHDPlusNationalData"
        / "NHDPlusV21_National_Seamless_Flattened_Lower48.gdb"
    )
    out_path = Path(out_path) if out_path else paths.DERIVED / "nhd_merit_crosswalk.parquet"

    print(f"Loading MERIT flowlines from {paths.MERIT_RIV}...")
    t0 = time.time()
    merit_full = gpd.read_file(
        paths.MERIT_RIV,
        columns=["COMID", "lengthkm"],
        engine="pyogrio",
    )
    print(f"  {len(merit_full):,} MERIT reaches loaded in {time.time()-t0:.1f}s")

    # Build a CONUS tile grid from MERIT bounds (EPSG:4326).
    merit_4326 = merit_full.to_crs("EPSG:4326")
    minx, miny, maxx, maxy = merit_4326.total_bounds
    # Small margin so edge reaches aren't missed.
    margin = 0.1  # degrees
    minx -= margin; miny -= margin; maxx += margin; maxy += margin

    xs = [minx + i * (maxx - minx) / tile_n for i in range(tile_n + 1)]
    ys = [miny + j * (maxy - miny) / tile_n for j in range(tile_n + 1)]

    OVERLAP_DEG = 0.01  # ~1 km — prevents seam losses

    all_pairs: list[pd.DataFrame] = []
    total_tiles = tile_n * tile_n
    tile_idx = 0

    print(f"Building crosswalk over {total_tiles} tiles ({tile_n}×{tile_n} grid)...")
    t_start = time.time()

    for i in range(tile_n):
        for j in range(tile_n):
            tile_idx += 1
            bbox = (
                xs[i] - OVERLAP_DEG,
                ys[j] - OVERLAP_DEG,
                xs[i + 1] + OVERLAP_DEG,
                ys[j + 1] + OVERLAP_DEG,
            )

            # Subset MERIT to tile.
            merit_tile = merit_4326.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]
            if merit_tile.empty:
                continue

            # Load NHD flowlines for this bbox from the GDB.
            try:
                nhd_tile = pyogrio.read_dataframe(
                    str(gdb_path),
                    layer="NHDFlowline_Network",
                    bbox=bbox,
                    columns=["COMID", "GNIS_NAME"],
                )
            except Exception as e:
                print(f"  tile {tile_idx}/{total_tiles} bbox={bbox}: NHD read error: {e}")
                continue

            if nhd_tile.empty:
                continue

            nhd_tile = nhd_tile.rename(columns={"COMID": "nhd_comid"})

            # Crosswalk on this tile.
            pairs = build_crosswalk(
                merit_tile.to_crs(paths.CRS_EQUAL_AREA),
                nhd_tile,
                buffer_m=paths.CROSSWALK_BUFFER_M,
                top_k=paths.CROSSWALK_TOP_K,
            )
            if not pairs.empty:
                # Carry GNIS_NAME for spot-check (merged back via nhd_comid lookup).
                all_pairs.append(pairs)

            elapsed = time.time() - t_start
            rate = tile_idx / elapsed if elapsed > 0 else 0
            remaining = (total_tiles - tile_idx) / rate if rate > 0 else 0
            print(
                f"  tile {tile_idx:3d}/{total_tiles} | "
                f"merit={len(merit_tile):6,} nhd={len(nhd_tile):6,} "
                f"pairs={len(pairs) if not pairs.empty else 0:6,} | "
                f"elapsed={elapsed/60:.1f}m ETA={remaining/60:.1f}m"
            )

    if not all_pairs:
        print("ERROR: no crosswalk pairs produced. Check GDB path and layer names.")
        return

    print("Concatenating and deduplicating...")
    combined = pd.concat(all_pairs, ignore_index=True)
    print(f"  {len(combined):,} pairs before dedup")

    # Dedupe (COMID, foreign_id) keeping max part_len (from the tile where it
    # falls most cleanly inside the buffer), then re-rank per COMID.
    combined = (
        combined.sort_values("part_len", ascending=False)
        .drop_duplicates(subset=["COMID", "foreign_id"])
        .sort_values(["COMID", "part_len"], ascending=[True, False])
    )

    # Recompute match_frac from the deduplicated pairs using MERIT lengths.
    merit_len = merit_full.copy()
    merit_len["merit_len"] = merit_full.to_crs(paths.CRS_EQUAL_AREA).geometry.length
    merit_len = merit_len[["COMID", "merit_len"]]

    combined = combined.merge(merit_len, on="COMID", how="left")
    frac = (
        combined.groupby("COMID")
        .apply(
            lambda g: min(1.0, g["part_len"].sum() / g["merit_len"].iloc[0]),
            include_groups=False,
        )
        .rename("match_frac")
    )
    # Drop per-pair match_frac from tiled runs (it was tile-local), attach global one.
    combined = combined.drop(columns=["match_frac", "merit_len"], errors="ignore")
    combined = combined.merge(frac, on="COMID")

    # Top-k re-ranking after dedup.
    combined = (
        combined.sort_values(["COMID", "part_len"], ascending=[True, False])
        .groupby("COMID", sort=False)
        .head(paths.CROSSWALK_TOP_K)
        .reset_index(drop=True)
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(out_path, index=False)

    total_pairs = len(combined)
    matched_comids = combined["COMID"].nunique()
    total_merit = len(merit_full)
    print(f"\nOutput: {out_path}")
    print(f"  Total pairs: {total_pairs:,}")
    print(f"  Matched MERIT COMIDs: {matched_comids:,} / {total_merit:,} "
          f"({matched_comids/total_merit:.1%})")

    # Coverage breakdown by uparea (need to reload with uparea).
    merit_up = gpd.read_file(paths.MERIT_RIV, columns=["COMID", "uparea"], engine="pyogrio")
    matched_set = set(combined["COMID"].unique())
    large = merit_up[merit_up["uparea"] > 100]
    large_matched = large[large["COMID"].isin(matched_set)]
    print(f"  uparea>100 km² matched: {len(large_matched):,} / {len(large):,} "
          f"({len(large_matched)/len(large):.1%})")

    mf = combined.groupby("COMID")["match_frac"].first()
    print(f"\nmatch_frac quartiles:")
    print(f"  p25={mf.quantile(0.25):.3f}  median={mf.quantile(0.50):.3f}  "
          f"p75={mf.quantile(0.75):.3f}  mean={mf.mean():.3f}")

    total_elapsed = time.time() - t_start
    print(f"\nTotal wall time: {total_elapsed/60:.1f} min")


if __name__ == "__main__":
    main()
