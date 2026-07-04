"""SWORD/GRWL observed widths -> MERIT COMIDs via spatial join into MERIT
catchment polygons (bugfix1 version).

Background (2026-07-04): The Wade et al. 2025 (Zenodo 13152826) mb_to_sword
translation tables use a pre-bugfix MERIT-Basins COMID ordering that does NOT
match the COMIDs in the bugfix1 shapefile used by this project.  The same
integer COMID refers to completely different reaches in the two versions
(distances of 100–1100 km are common within a single pfaf-2 region), producing
zero rank correlation against bankfull widths when the table is used directly.

Fix: assign each SWORD NA reach to a MERIT bugfix1 COMID by spatial join of
the SWORD reach centerpoint (x, y) into the MERIT bugfix1 catchment polygons.
Only river reaches (lakeflag == 0) are used; lake/reservoir SWORD widths are
not representative of channel width.

Outputs derived/channel_width_obs.parquet (columns: COMID, channel_width_obs).
"""
import netCDF4
import numpy as np
import pandas as pd
import geopandas as gpd

from . import paths
from .transfer import weighted_transfer

_SWORD_NA = paths.RAW / "sword" / "netcdf" / "na_sword_v16.nc"

# CONUS pfaf-2 regions
CONUS_PFAF = {71, 72, 73, 74, 75, 76, 77, 78}


def load_sword_reaches() -> pd.DataFrame:
    """Load SWORD v16 NA reaches; return DataFrame with geometry columns.

    Only CONUS reaches (pfaf-2 regions 71-78) with positive width and
    lakeflag == 0 (river reaches, not lakes / reservoirs) are returned.
    Negative width is SWORD's missing-data sentinel and is excluded.
    """
    ds = netCDF4.Dataset(_SWORD_NA)
    rg = ds.groups["reaches"]
    reach_ids = rg["reach_id"][:].data
    xs = rg["x"][:].data
    ys = rg["y"][:].data
    widths = rg["width"][:].data.astype(float)
    reach_lengths = rg["reach_length"][:].data.astype(float)
    lakeflags = rg["lakeflag"][:].data
    ds.close()

    df = pd.DataFrame({
        "foreign_id": reach_ids,
        "x": xs,
        "y": ys,
        "width": widths,
        "reach_length": reach_lengths,
        "lakeflag": lakeflags,
    })

    pfaf = df["foreign_id"] // 1_000_000_000
    df = df[pfaf.isin(CONUS_PFAF) & (df["width"] > 0) & (df["lakeflag"] == 0)]
    return df.reset_index(drop=True)


def build_crosswalk(sword_df: pd.DataFrame) -> pd.DataFrame:
    """Spatial join SWORD reach centerpoints into MERIT bugfix1 catchments.

    Returns long DataFrame(COMID, foreign_id, part_len) where part_len is the
    SWORD reach_length (used as the weight in the subsequent transfer).
    """
    sword_gdf = gpd.GeoDataFrame(
        sword_df[["foreign_id", "reach_length"]],
        geometry=gpd.points_from_xy(sword_df["x"], sword_df["y"]),
        crs="EPSG:4326",
    )

    cat = gpd.read_file(paths.MERIT_CAT, columns=["COMID"])
    cat = cat.set_crs("EPSG:4326", allow_override=True)

    joined = gpd.sjoin(sword_gdf, cat[["COMID", "geometry"]], how="left", predicate="within")

    matched = joined[joined["COMID"].notna()][["foreign_id", "reach_length", "COMID"]].copy()
    unmatched = joined[joined["COMID"].isna()].copy()

    # Recover boundary reaches with nearest-catchment join (projected CRS).
    if len(unmatched) > 0:
        sword_proj = unmatched[["foreign_id", "reach_length", "geometry"]].to_crs(paths.CRS_EQUAL_AREA)
        cat_proj = cat[["COMID", "geometry"]].to_crs(paths.CRS_EQUAL_AREA)
        nearest = gpd.sjoin_nearest(sword_proj, cat_proj, how="left", max_distance=10_000)
        nearest = nearest[nearest["COMID"].notna()][["foreign_id", "reach_length", "COMID"]].copy()
        matched = pd.concat([matched, nearest], ignore_index=True)

    xwalk = matched[["COMID", "foreign_id", "reach_length"]].copy()
    xwalk["COMID"] = xwalk["COMID"].astype(int)
    xwalk = xwalk.rename(columns={"reach_length": "part_len"})
    return xwalk.reset_index(drop=True)


def main() -> None:
    print("Loading SWORD v16 NA river reaches (CONUS, lakeflag==0, width>0) ...")
    sword = load_sword_reaches()
    print(f"  SWORD CONUS river reaches: {len(sword):,}")

    print("Building MERIT-SWORD crosswalk via spatial join into bugfix1 catchments ...")
    xwalk = build_crosswalk(sword)
    print(f"  Crosswalk pairs: {len(xwalk):,}  unique COMIDs: {xwalk['COMID'].nunique():,}")

    print("Running weighted transfer ...")
    out = weighted_transfer(xwalk, sword[["foreign_id", "width"]], value_col="width")
    out = out.rename(columns={"width": "channel_width_obs"})

    covered = out["channel_width_obs"].notna()
    total = len(out)
    n_covered = covered.sum()
    print(f"\nCoverage: {n_covered:,} / {total:,} MERIT COMIDs "
          f"({100 * n_covered / total:.1f}%)")
    print(f"Median width (covered reaches): "
          f"{out.loc[covered, 'channel_width_obs'].median():.1f} m")
    print(f"Width distribution (covered):")
    print(out.loc[covered, "channel_width_obs"].describe().to_string())

    out_path = paths.DERIVED / "channel_width_obs.parquet"
    out.to_parquet(out_path, index=False)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
