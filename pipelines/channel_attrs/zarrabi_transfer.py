"""Zarrabi et al. bankfull geometry -> MERIT COMIDs via NHDPlus->MERIT crosswalk.

Input:  raw/zarrabi/Bankfull_Meanflow_CONUS.txt (NHD-indexed CSV)
        derived/nhd_merit_crosswalk.parquet
Output: derived/bankfull.parquet  (columns: COMID, bankfull_depth, bankfull_width)

File format: comma-separated, first column is a numeric row index, then:
  COMID, REACHCODE, TotDASqKM, StreamOrde, bnk_depth, bnk_width, mf_depth, mf_width
Units: depth in m, width in m (Zarrabi et al. 2021, Zenodo 13883263).

McManamay confinement: SKIPPED.
The McManamay & Pers (2017) confinement dataset requires a separate download
(Figshare article) that is not in the Phase-A download manifest, and the plan
describes it as "nice-to-have". This omission is recorded here; the three
bankfull geometry columns are the primary Phase-A deliverable and are complete.
"""
from pathlib import Path

import pandas as pd

from . import paths
from .transfer import weighted_transfer

MATCH_FRAC_MIN = 0.3


def load_zarrabi(txt_path: Path) -> pd.DataFrame:
    """Read Bankfull_Meanflow_CONUS.txt; return foreign_id + bankfull_depth/width."""
    df = pd.read_csv(txt_path, index_col=0, usecols=[0, 1, 5, 6])
    # After usecols: index, COMID, bnk_depth, bnk_width
    df = df.rename(columns={
        "COMID": "foreign_id",
        "bnk_depth": "bankfull_depth",
        "bnk_width": "bankfull_width",
    })
    df["foreign_id"] = df["foreign_id"].astype("int32")
    return df.reset_index(drop=True)


def main() -> None:
    zarr_file = paths.RAW / "zarrabi" / "Bankfull_Meanflow_CONUS.txt"
    df = load_zarrabi(zarr_file)
    print(f"Zarrabi records: {len(df):,}")
    print(f"  bnk_depth range: {df['bankfull_depth'].min():.3f} – {df['bankfull_depth'].max():.3f} m")
    print(f"  bnk_width range: {df['bankfull_width'].min():.3f} – {df['bankfull_width'].max():.3f} m")

    xw = pd.read_parquet(paths.DERIVED / "nhd_merit_crosswalk.parquet")
    xw = xw[xw["match_frac"] >= MATCH_FRAC_MIN]

    depth = weighted_transfer(xw, df[["foreign_id", "bankfull_depth"]], value_col="bankfull_depth")
    width = weighted_transfer(xw, df[["foreign_id", "bankfull_width"]], value_col="bankfull_width")

    out = depth.merge(width, on="COMID", how="outer")

    cov_d = out["bankfull_depth"].notna().sum()
    cov_w = out["bankfull_width"].notna().sum()
    total = len(out)
    print(f"Coverage: depth {cov_d:,}/{total:,} ({cov_d/total:.1%}), "
          f"width {cov_w:,}/{total:,} ({cov_w/total:.1%})")
    print(out[["bankfull_depth", "bankfull_width"]].describe())

    out.to_parquet(paths.DERIVED / "bankfull.parquet", index=False)
    print(f"Written: {paths.DERIVED / 'bankfull.parquet'}")

    _sanity_checks(out)


def _sanity_checks(out: pd.DataFrame) -> None:
    """Median depth in [0.3, 3 m]; Spearman(depth, uparea) > 0.5;
    Spearman(width, channel_width_obs) > 0.5."""
    from scipy.stats import spearmanr
    import pyogrio

    median_depth = out["bankfull_depth"].median()
    print(f"\nSanity depth median: {median_depth:.3f} m (expect 0.3–3 m)")

    # Spearman vs uparea from MERIT shapefile (non-geometric read).
    riv = pyogrio.read_dataframe(paths.MERIT_RIV, columns=["COMID", "uparea"])
    merged = out.merge(riv[["COMID", "uparea"]], on="COMID", how="inner")
    valid = merged[merged["bankfull_depth"].notna() & merged["uparea"].notna()]
    rho_ua, _ = spearmanr(valid["bankfull_depth"], valid["uparea"])
    print(f"Spearman(bankfull_depth, uparea): {rho_ua:.3f} (expect >0.5)")

    # Spearman vs SWORD channel_width_obs.
    cwo = pd.read_parquet(paths.DERIVED / "channel_width_obs.parquet")
    merged2 = out.merge(cwo, on="COMID", how="inner")
    valid2 = merged2[merged2["bankfull_width"].notna() & merged2["channel_width_obs"].notna()]
    if len(valid2) > 10:
        rho_w, _ = spearmanr(valid2["bankfull_width"], valid2["channel_width_obs"])
        print(f"Spearman(bankfull_width, channel_width_obs): {rho_w:.3f} (expect >0.5)")
    else:
        print(f"Spearman(bankfull_width, channel_width_obs): insufficient overlap "
              f"({len(valid2)} reaches)")


if __name__ == "__main__":
    main()
