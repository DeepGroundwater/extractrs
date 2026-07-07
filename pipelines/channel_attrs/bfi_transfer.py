"""BFI (baseflow index) + drainage density -> MERIT COMIDs.

BFI:
  Source: USGS Wolock BFI_CONUS.zip (NHD-indexed table, ScienceBase).
  Column: CAT_BFI — local catchment BFI in [0, 100] percent; divided by 100.
  Transfer: NHDPlus->MERIT crosswalk (Task 3), match_frac >= 0.3.

  NOTE: The Phase-A plan (Task 9, Step 1) described BFI transfer as a raster
  zonal-stats operation using bfi48grd. BFI_CONUS.zip is instead an NHD-indexed
  table (CAT_BFI already aggregated per NHD COMID by EPA), so we use the
  crosswalk approach directly without rioxarray. This deviates from the plan's
  raster path but is equivalent (both derive from the same Wolock BFI grid) and
  is simpler. Deviation recorded here.

Drainage density:
  lengthkm (km) from MERIT pfaf-7 shapefile divided by catchsize (km²) from
  merit_global_attributes_v2.nc.  Guard: catchsize <= 0 -> NaN.

Output: derived/basin_extras.parquet  (columns: COMID, bfi, drainage_density)
"""
import zipfile
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import pyogrio

from . import paths
from .transfer import weighted_transfer

MATCH_FRAC_MIN = 0.3


def load_bfi(zip_path: Path) -> pd.DataFrame:
    """Read BFI_CONUS.zip -> DataFrame(foreign_id int32, bfi float [0–100 %]).

    Negative COMIDs are coastal pseudo-COMIDs and are filtered out.
    """
    with zipfile.ZipFile(zip_path) as zf:
        txt_name = next(n for n in zf.namelist() if n.endswith(".txt"))
        with zf.open(txt_name) as f:
            df = pd.read_csv(f, usecols=["COMID", "CAT_BFI"])
    df = df[df["COMID"] > 0].copy()
    df = df.rename(columns={"COMID": "foreign_id", "CAT_BFI": "bfi"})
    df["foreign_id"] = df["foreign_id"].astype("int32")
    return df.reset_index(drop=True)


def compute_drainage_density(riv_path: Path, global_nc_path: Path) -> pd.DataFrame:
    """lengthkm / catchsize per MERIT COMID -> drainage_density (km⁻¹).

    catchsize <= 0 is guarded to NaN (no zero-area reaches should exist in
    MERIT, but the global NC includes non-CONUS COMIDs with placeholder values).
    """
    riv = pyogrio.read_dataframe(riv_path, columns=["COMID", "lengthkm"])

    ds = netCDF4.Dataset(global_nc_path)
    comids = np.array(ds["COMID"][:], dtype="int64")
    catchsize = np.array(ds["catchsize"][:], dtype="float64")
    ds.close()

    cat_df = pd.DataFrame({"COMID": comids, "catchsize": catchsize})
    merged = riv.merge(cat_df, on="COMID", how="left")
    merged["drainage_density"] = np.where(
        merged["catchsize"] > 0,
        merged["lengthkm"] / merged["catchsize"],
        np.nan,
    )
    return merged[["COMID", "drainage_density"]]


def main() -> None:
    bfi_raw = load_bfi(paths.RAW / "bfi" / "BFI_CONUS.zip")
    print(f"BFI records (positive COMID): {len(bfi_raw):,}")
    print(f"CAT_BFI range: {bfi_raw['bfi'].min():.2f} – {bfi_raw['bfi'].max():.2f} % (expect [0, 100])")

    xw = pd.read_parquet(paths.DERIVED / "nhd_merit_crosswalk.parquet")
    xw = xw[xw["match_frac"] >= MATCH_FRAC_MIN]

    bfi_out = weighted_transfer(xw, bfi_raw, value_col="bfi")
    bfi_out["bfi"] /= 100.0

    cov = bfi_out["bfi"].notna().sum()
    total = len(bfi_out)
    print(f"BFI coverage: {cov:,}/{total:,} ({cov/total:.1%})")
    print(bfi_out["bfi"].describe())

    dd = compute_drainage_density(paths.MERIT_RIV, paths.GLOBAL_NC)
    print(f"\nDrainage density: median={dd['drainage_density'].median():.3f} km⁻¹ (expect 0.1–2)")

    out = bfi_out.merge(dd, on="COMID", how="outer")
    out.to_parquet(paths.DERIVED / "basin_extras.parquet", index=False)
    print(f"Written: {paths.DERIVED / 'basin_extras.parquet'}")

    _sanity_checks(bfi_out)


def _sanity_checks(bfi_out: pd.DataFrame) -> None:
    """BFI in [0,1]; humid East > arid Southwest spot check."""
    import pyogrio

    assert bfi_out["bfi"].dropna().between(0, 1).all(), "BFI out of [0, 1] range"

    # Eastern reaches: bbox for Mid-Atlantic (Virginia / North Carolina)
    riv_east = pyogrio.read_dataframe(
        paths.MERIT_RIV, columns=["COMID", "COMID"],
        bbox=(-82.0, 35.0, -76.0, 38.0),
    )
    east_comids = set(riv_east["COMID"].tolist())
    east_bfi = bfi_out.loc[bfi_out["COMID"].isin(east_comids), "bfi"].dropna()

    # Southwest: Arizona / New Mexico
    riv_sw = pyogrio.read_dataframe(
        paths.MERIT_RIV, columns=["COMID", "COMID"],
        bbox=(-114.0, 33.0, -107.0, 37.0),
    )
    sw_comids = set(riv_sw["COMID"].tolist())
    sw_bfi = bfi_out.loc[bfi_out["COMID"].isin(sw_comids), "bfi"].dropna()

    print(f"\nSanity BFI East (Mid-Atlantic, {len(east_bfi)} reaches): median={east_bfi.median():.3f}")
    print(f"Sanity BFI SW (AZ/NM, {len(sw_bfi)} reaches): median={sw_bfi.median():.3f}")
    if len(east_bfi) and len(sw_bfi):
        print(f"  East > SW: {east_bfi.median() > sw_bfi.median()} (expect True)")


if __name__ == "__main__":
    main()
