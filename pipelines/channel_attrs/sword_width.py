"""SWORD/GRWL observed widths -> MERIT COMIDs via the published MERIT-SWORD
translation tables (Wade et al. 2025, Zenodo 13152826).

Adaptation notes (verified 2026-07-04):
- Translation table MERIT reach ID variable: `mb` (NOT `COMID` as the plan
  draft assumed). Dimension and coordinate variable both named `mb`.
  CONUS values like 71000001 (pfaf-2 prefix + 6-digit reach index).
- CONUS = pfaf-2 regions 71–78 (8 files: mb_to_sword_pfaf_7[1-8]_translate.nc).
- SWORD v16 NA file (`na_sword_v16.nc`) has a `/reaches` group with
  `reach_id` (int64) and `width` (float64, 38 048 NA reaches).
  xarray decode fails on that file's string variables; read with netCDF4.
- SWORD `width` range: -1.0 to 69 823 m; negative sentinel = missing data,
  filtered out before transfer.
"""
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr

from . import paths
from .transfer import weighted_transfer

# Extracted translation directory
_TRANS_DIR = paths.RAW / "merit_sword" / "ms_translate" / "mb_to_sword"
_SWORD_NA = paths.RAW / "sword" / "netcdf" / "na_sword_v16.nc"

# CONUS pfaf-2 regions (confirmed from zip listing)
CONUS_REGIONS = [71, 72, 73, 74, 75, 76, 77, 78]


def melt_translation(ds: xr.Dataset, comid_var: str = "mb") -> pd.DataFrame:
    """Melt a wide translation table into long (COMID, foreign_id, part_len) form.

    Translation files have sword_1..sword_40 (SWORD reach_id) and
    part_len_1..part_len_40 (overlap length in m) on a MERIT-reach dimension.
    Zero values are absent/no-match and are dropped.
    """
    frames = []
    for k in range(1, 41):
        s, p = f"sword_{k}", f"part_len_{k}"
        if s not in ds or p not in ds:
            break
        df = pd.DataFrame({
            "COMID": ds[comid_var].values,
            "foreign_id": ds[s].values,
            "part_len": ds[p].values,
        })
        frames.append(df[(df.foreign_id > 0) & (df.part_len > 0)])
    if not frames:
        return pd.DataFrame(columns=["COMID", "foreign_id", "part_len"])
    return pd.concat(frames, ignore_index=True)


def load_sword_widths() -> pd.DataFrame:
    """Load SWORD v16 NA reach widths; return DataFrame(foreign_id, width).

    Uses netCDF4 directly — xarray raises AttributeError on the string
    variable `river_name` in this file.
    Negative widths are SWORD's missing-data sentinel and are set to NaN.
    """
    ds = netCDF4.Dataset(_SWORD_NA)
    rg = ds.groups["reaches"]
    reach_ids = rg["reach_id"][:].data
    widths = rg["width"][:].data.astype(float)
    ds.close()
    widths[widths < 0] = np.nan
    return pd.DataFrame({"foreign_id": reach_ids, "width": widths})


def main() -> None:
    print("Loading SWORD v16 NA widths ...")
    sword = load_sword_widths()
    print(f"  SWORD NA reaches: {len(sword):,}  (non-null widths: {sword['width'].notna().sum():,})")

    print("Melting MERIT-SWORD translation tables for CONUS regions 71-78 ...")
    xwalk_frames = []
    for region in CONUS_REGIONS:
        nc_path = _TRANS_DIR / f"mb_to_sword_pfaf_{region}_translate.nc"
        ds = xr.open_dataset(nc_path)
        frame = melt_translation(ds, comid_var="mb")
        ds.close()
        print(f"  pfaf-{region}: {len(frame):,} (COMID, SWORD) pairs")
        xwalk_frames.append(frame)

    xwalk = pd.concat(xwalk_frames, ignore_index=True)
    print(f"Total crosswalk pairs: {len(xwalk):,}")

    print("Running weighted transfer ...")
    out = weighted_transfer(xwalk, sword, value_col="width")
    out = out.rename(columns={"width": "channel_width_obs"})

    # Coverage stats
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
