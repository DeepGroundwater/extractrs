"""Assemble Phase-A per-COMID attributes into merit_channel_attributes_v1.nc.

Global COMID vector from merit_global_attributes_v2.nc (2,939,404 int64).
Each parquet is reindexed onto that vector; NaN outside coverage.
Variables: one f64 per attribute on dim COMID — byte-schema-identical to
merit_global_attributes_v2.nc (no extra dims, no groups).
Stats JSON over finite values only, six-key format matching ddrs convention.
"""
import json
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

from . import paths


# (variable_name, parquet_path, column_in_parquet)
VARIABLES = [
    ("channel_width_obs",   paths.DERIVED / "channel_width_obs.parquet",   "channel_width_obs"),
    ("corridor_impervious", paths.DERIVED / "corridor_impervious.parquet",  "corridor_impervious"),
    ("bankfull_depth",      paths.DERIVED / "bankfull.parquet",             "bankfull_depth"),
    ("bankfull_width",      paths.DERIVED / "bankfull.parquet",             "bankfull_width"),
    ("wtd_channel_zs",      paths.DERIVED / "wtd_channel.parquet",          "wtd_channel_zs"),
    ("wtd_channel_fan",     paths.DERIVED / "wtd_channel.parquet",          "wtd_channel_fan"),
    ("channel_wtd_bed_rel", paths.DERIVED / "wtd_bedrel.parquet",           "channel_wtd_bed_rel"),
    ("losing_fraction",     paths.DERIVED / "wtd_bedrel.parquet",           "losing_fraction"),
    ("alluvium_fraction",   paths.DERIVED / "alluvium_fraction.parquet",    "alluvium_fraction"),
    ("bfi",                 paths.DERIVED / "basin_extras.parquet",         "bfi"),
    ("drainage_density",    paths.DERIVED / "basin_extras.parquet",         "drainage_density"),
]


def assemble(global_comids: np.ndarray, frames: dict, out: Path) -> None:
    """Write a netCDF4 file with one f64 var per attribute aligned to global_comids.

    Parameters
    ----------
    global_comids:
        1-D int64 array of all global COMID values (the coordinate axis).
    frames:
        Dict of variable_name -> DataFrame with columns [COMID, <variable_name>].
        COMIDs absent from a frame become NaN in the output.
    out:
        Output netCDF4 path (overwritten if exists).
    """
    n = len(global_comids)
    index = pd.Index(global_comids, name="COMID")

    with netCDF4.Dataset(out, "w") as nc:
        nc.createDimension("COMID", n)
        cv = nc.createVariable("COMID", "i8", ("COMID",))
        cv[:] = global_comids

        for name, df in frames.items():
            arr = df.set_index("COMID")[name].reindex(index).values.astype("float64")
            var = nc.createVariable(name, "f8", ("COMID",))
            var[:] = arr


def _load_frames() -> dict:
    """Read each parquet once (cached by path) and return the frames dict."""
    cache: dict[str, pd.DataFrame] = {}
    frames: dict[str, pd.DataFrame] = {}
    for var_name, parquet_path, col in VARIABLES:
        key = str(parquet_path)
        if key not in cache:
            cache[key] = pd.read_parquet(parquet_path)
        df = cache[key][["COMID", col]].rename(columns={col: var_name})
        frames[var_name] = df
    return frames


def _compute_stats(nc_path: Path) -> dict:
    """Compute per-variable stats over finite values; matches ddrs six-key format."""
    stats = {}
    with netCDF4.Dataset(nc_path) as nc:
        nc.set_auto_mask(False)
        for name in nc.variables:
            if name == "COMID":
                continue
            arr = np.array(nc[name][:], dtype="float64")
            finite = arr[np.isfinite(arr)]
            if len(finite) == 0:
                continue
            stats[name] = {
                "min":  float(np.min(finite)),
                "max":  float(np.max(finite)),
                "mean": float(np.mean(finite)),
                "std":  float(np.std(finite)),
                "p10":  float(np.percentile(finite, 10)),
                "p90":  float(np.percentile(finite, 90)),
            }
    return stats


def main() -> None:
    # 1. Global COMID vector (2,939,404 int64)
    with netCDF4.Dataset(paths.GLOBAL_NC) as nc:
        global_comids = np.array(nc["COMID"][:], dtype="int64")
    print(f"Global COMID vector: {len(global_comids):,} entries")

    # 2. Load parquets and assemble
    frames = _load_frames()
    assemble(global_comids, frames, paths.OUT_NC)
    print(f"Written: {paths.OUT_NC}")

    # 3. Stats JSON — same naming convention as existing ddrs stats files:
    #    merit_attribute_statistics_<nc_basename>.json
    stats = _compute_stats(paths.OUT_NC)
    stats_path = paths.STATS_DIR / f"merit_attribute_statistics_{paths.OUT_NC.name}.json"
    stats_path.write_text(json.dumps(stats, indent=2))
    print(f"Stats: {stats_path}")

    # 4. Validation per plan step 4
    print("\n--- Validation ---")
    with netCDF4.Dataset(paths.GLOBAL_NC) as a, netCDF4.Dataset(paths.OUT_NC) as b:
        a.set_auto_mask(False)
        b.set_auto_mask(False)
        n_global = len(a.dimensions["COMID"])
        n_channel = len(b.dimensions["COMID"])
        assert n_global == n_channel, f"dim mismatch: {n_global} vs {n_channel}"
        assert (np.array(a["COMID"][:100]) == np.array(b["COMID"][:100])).all(), \
            "COMID ordering diverges in first 100"
        for v in b.variables:
            if v == "COMID":
                continue
            assert b[v].dimensions == ("COMID",), f"{v}: unexpected dims"
            arr = np.array(b[v][:], dtype="float64")
            assert arr.dtype == np.float64, v
            finite = int(np.isfinite(arr).sum())
            print(f"  {v}: finite={finite:,} ({finite / n_channel:.1%})")
    print("PASS: dim length + COMID order + schema OK")


if __name__ == "__main__":
    main()
