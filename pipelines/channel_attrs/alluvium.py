"""Alluvium fraction: Principal Aquifer polygon overlay on scaled corridors.

Deviation from plan Task 8: overlay is on ``corridors_scaled.parquet``
(not ``corridors_100m`` as the plan text says).  Rationale: large rivers
lie over alluvial aquifers and need the scaled corridor to actually cover
the channel + its alluvial floor; a fixed 100 m strip may miss the alluvial
body entirely for order-8/9 reaches.  Result is documented and the
100 m version is still available as a sensitivity run.

Alluvial classes chosen: ROCK_TYPE == 100 ("Unconsolidated sand and gravel
aquifers") from the USGS Principal Aquifers shapefile (ScienceBase item
10.5066/P9Y2HOUJ).  This class contains every explicitly alluvial/basin-fill
AQ_NAME (Mississippi Valley, Pecos, Basin-and-Range fill, Columbia Plateau
fill, etc.).  ROCK_TYPE 200 (Semiconsolidated sand, e.g. High Plains,
coastal-plain Cretaceous) is excluded: those are sedimentary bedrock, not
modern river alluvium.

Output: ``derived/alluvium_fraction.parquet``  (columns: COMID int64,
alluvium_fraction float64 ∈ [0,1]; 0 = no mapped alluvium, explicit fill).
"""
import time

import geopandas as gpd
import numpy as np
import pandas as pd

from . import paths

# ROCK_TYPE = 100: "Unconsolidated sand and gravel aquifers" — the alluvial class.
# Basis: all explicitly alluvial/basin-fill AQ_NAMEs fall exclusively in this
# class (verified by inspecting the shapefile: ROCK_TYPE × AQ_NAME cross-tab).
ALLUVIAL_ROCK_TYPES = (100,)


def load_alluvial_polys() -> gpd.GeoDataFrame:
    """Load and filter principal aquifer polygons to alluvial rock types."""
    aq = gpd.read_file(paths.RAW / "aquifers" / "us_aquifers.shp")
    alluvial = aq[aq["ROCK_TYPE"].isin(ALLUVIAL_ROCK_TYPES)].copy()
    alluvial = alluvial[["ROCK_TYPE", "ROCK_NAME", "AQ_NAME", "geometry"]]
    alluvial = alluvial.to_crs(paths.CRS_EQUAL_AREA)
    # Dissolve so the overlay sees one geometry per alluvial body, not 806 sliver polygons.
    alluvial = alluvial.dissolve().reset_index(drop=True)[["geometry"]]
    alluvial["is_alluvial"] = 1
    return alluvial


def overlay_chunk(
    corr_chunk: gpd.GeoDataFrame,
    alluvial: gpd.GeoDataFrame,
) -> pd.DataFrame:
    """Compute alluvium_fraction = intersect_area / corridor_area per COMID.

    Returns a DataFrame with columns COMID, alluvium_fraction.  Reaches with
    no intersection are NOT returned (caller fills with 0 after merging).
    """
    # Pre-compute corridor areas (EPSG:5070 → metres²).
    corr_chunk = corr_chunk.copy()
    corr_chunk["corridor_area"] = corr_chunk.geometry.area

    # Intersection overlay: returns only overlapping pieces.
    try:
        inter = gpd.overlay(corr_chunk[["COMID", "corridor_area", "geometry"]],
                            alluvial,
                            how="intersection",
                            keep_geom_type=True)
    except Exception:
        # If overlay produces no results (no alluvial coverage in this chunk).
        return pd.DataFrame(columns=["COMID", "alluvium_fraction"])

    if inter.empty:
        return pd.DataFrame(columns=["COMID", "alluvium_fraction"])

    inter["intersect_area"] = inter.geometry.area

    # Per-COMID: sum intersect area, take corridor_area (constant per COMID).
    agg = inter.groupby("COMID").agg(
        intersect_area=("intersect_area", "sum"),
        corridor_area=("corridor_area", "first"),
    ).reset_index()
    agg["alluvium_fraction"] = (agg["intersect_area"] / agg["corridor_area"]).clip(0.0, 1.0)
    return agg[["COMID", "alluvium_fraction"]]


def main() -> None:
    t_start = time.time()
    print(f"Loading corridors_scaled from {paths.DERIVED / 'corridors_scaled.parquet'} ...")
    corr = gpd.read_parquet(paths.DERIVED / "corridors_scaled.parquet")
    print(f"  {len(corr):,} corridors loaded")

    print("Loading alluvial aquifer polygons ...")
    alluvial = load_alluvial_polys()
    print(f"  {len(alluvial):,} dissolved alluvial polygon(s)  "
          f"(ROCK_TYPE {ALLUVIAL_ROCK_TYPES})")

    # Chunk by pfaf-4 prefix (first 4 digits of COMID string) for memory management.
    comid_str = corr["COMID"].astype(str)
    pfaf4 = comid_str.str[:4]
    chunks = sorted(pfaf4.unique())
    print(f"  Processing {len(chunks)} pfaf-4 chunks ...")

    results: list[pd.DataFrame] = []
    for i, p4 in enumerate(chunks):
        mask = pfaf4 == p4
        chunk = corr[mask].reset_index(drop=True)
        t0 = time.time()
        frac = overlay_chunk(chunk, alluvial)
        n_hit = len(frac)
        elapsed = time.time() - t0
        print(f"  [{i+1:02d}/{len(chunks)}] pfaf-4={p4}  {len(chunk):,} corridors  "
              f"→ {n_hit:,} with alluvium  ({elapsed:.1f}s)")
        results.append(frac)

    hits = pd.concat(results, ignore_index=True) if results else pd.DataFrame(
        columns=["COMID", "alluvium_fraction"]
    )

    # Left-join onto all COMIDs; fill 0 where no intersection.
    all_comids = pd.DataFrame({"COMID": corr["COMID"].values})
    out = all_comids.merge(hits, on="COMID", how="left")
    out["alluvium_fraction"] = out["alluvium_fraction"].fillna(0.0)

    # Verify coverage.
    assert len(out) == len(corr), f"Row count mismatch: {len(out)} != {len(corr)}"
    n_nonzero = (out["alluvium_fraction"] > 0).sum()
    total = len(out)
    print(f"\nAlluvium fraction summary:")
    print(f"  {n_nonzero:,}/{total:,} reaches have alluvium fraction > 0 ({n_nonzero/total:.1%})")
    print(out["alluvium_fraction"].describe())

    out_path = paths.DERIVED / "alluvium_fraction.parquet"
    out.to_parquet(out_path, index=False)
    print(f"\nWritten: {out_path}")
    print(f"Total elapsed: {time.time() - t_start:.1f}s")

    # Sanity checks.
    _sanity_checks(out, corr)


def _sanity_checks(out: pd.DataFrame, corr: gpd.GeoDataFrame) -> None:
    """Mississippi alluvial valley ≈ 1.0; Rocky Mountain headwaters ≈ 0."""
    import pyogrio
    # Mississippi alluvial valley: look at reaches near Vicksburg, MS (~32.3N, -90.9W)
    riv = pyogrio.read_dataframe(
        paths.MERIT_RIV, columns=["COMID"],
        bbox=(-91.5, 32.0, -90.5, 33.0),
    )
    miss_comids = set(riv["COMID"].tolist())
    miss_vals = out.loc[out["COMID"].isin(miss_comids), "alluvium_fraction"]
    if len(miss_vals):
        print(f"\nSanity Mississippi alluvial valley ({len(miss_vals)} reaches): "
              f"median={miss_vals.median():.3f}  (expect ~1.0)")
    else:
        print("\nSanity Mississippi: no reaches in bbox")

    # Rocky Mountain headwaters: Colorado Front Range (~40N, -105.5W)
    riv_co = pyogrio.read_dataframe(
        paths.MERIT_RIV, columns=["COMID"],
        bbox=(-106.0, 39.5, -105.0, 40.5),
    )
    co_comids = set(riv_co["COMID"].tolist())
    co_vals = out.loc[out["COMID"].isin(co_comids), "alluvium_fraction"]
    if len(co_vals):
        print(f"Sanity Rocky Mtn headwaters ({len(co_vals)} reaches): "
              f"median={co_vals.median():.3f}  (expect ~0)")
    else:
        print("Sanity Rocky Mtn: no reaches in bbox")


if __name__ == "__main__":
    main()
