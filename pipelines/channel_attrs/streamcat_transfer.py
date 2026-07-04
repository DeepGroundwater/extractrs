"""StreamCat pctimp2019catrp100 (100 m riparian buffer imperviousness, NHDPlusV2)
-> MERIT COMIDs via the Task-3 NHDPlus->MERIT crosswalk.

Input:  raw/streamcat/pctimp2019_Region*.json
        derived/nhd_merit_crosswalk.parquet
Output: derived/corridor_impervious.parquet  (columns: COMID, corridor_impervious)

API JSON structure: {"items": [{..., "comid": <int>, "pctimp2019catrp100": <float|null>, ...}]}
Values are in [0, 100] percent; divided by 100 -> [0, 1] fraction in the output.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import paths
from .transfer import weighted_transfer

MATCH_FRAC_MIN = 0.3


def _parse_region_json(path: Path) -> pd.DataFrame:
    """Stream-parse one region JSON file via ijson; tolerates truncated files.

    The API returns {"items": [...]}; each item has at minimum "comid" and
    "pctimp2019catrp100" (null is allowed). ijson emits complete item objects
    as they close, so a truncated final item is silently dropped rather than
    raising.
    """
    import ijson

    rows: list[dict] = []
    try:
        with path.open("rb") as f:
            for item in ijson.items(f, "items.item"):
                val = item.get("pctimp2019catrp100")
                rows.append({
                    "foreign_id": int(item["comid"]),
                    # ijson returns Decimal for JSON numbers; cast to float.
                    "corridor_impervious": float(val) if val is not None else None,
                })
    except (ijson.common.IncompleteJSONError, Exception):
        # Truncated file: keep whatever was parsed before the cut.
        pass
    if not rows:
        return pd.DataFrame(columns=["foreign_id", "corridor_impervious"])
    df = pd.DataFrame(rows)
    return df


def load_streamcat(raw_dir: Path) -> pd.DataFrame:
    """Parse all Region JSON files; return DataFrame with foreign_id + corridor_impervious.

    Handles truncated files gracefully via streaming parse (ijson): any
    partially-written final record is silently dropped.
    """
    frames = []
    for p in sorted(raw_dir.glob("pctimp2019_Region*.json")):
        df = _parse_region_json(p)
        n = len(df)
        if n == 0:
            print(f"  WARNING: {p.name} yielded 0 records (empty or unreadable)")
        else:
            print(f"  {p.name}: {n:,} records")
        frames.append(df)
    sc = pd.concat(frames, ignore_index=True)
    # Align dtype with crosswalk foreign_id (int32).
    sc["foreign_id"] = sc["foreign_id"].astype("int32")
    return sc


def main() -> None:
    sc = load_streamcat(paths.RAW / "streamcat")
    n_records = len(sc)
    n_nonnull = sc["corridor_impervious"].notna().sum()
    print(f"StreamCat records: {n_records:,}  non-null: {n_nonnull:,}")

    xw = pd.read_parquet(paths.DERIVED / "nhd_merit_crosswalk.parquet")
    xw = xw[xw["match_frac"] >= MATCH_FRAC_MIN]

    out = weighted_transfer(xw, sc, value_col="corridor_impervious")
    out["corridor_impervious"] /= 100.0

    coverage = out["corridor_impervious"].notna().sum()
    total = len(out)
    print(f"MERIT COMIDs with value: {coverage:,}/{total:,} ({coverage/total:.1%})")
    print(out["corridor_impervious"].describe())

    out.to_parquet(paths.DERIVED / "corridor_impervious.parquet", index=False)
    print(f"Written: {paths.DERIVED / 'corridor_impervious.parquet'}")

    # --- Sanity checks ---
    _sanity_checks(out)


def _sanity_checks(out: pd.DataFrame) -> None:
    """Spot-check LA River (urban, expect > 0.6) and rural Montana (expect < 0.05)."""
    import pyogrio
    riv = pyogrio.read_dataframe(
        paths.MERIT_RIV, columns=["COMID", "COMID"],
        bbox=(-118.5, 33.8, -117.9, 34.2),  # LA River bbox (lon_min, lat_min, lon_max, lat_max)
    )
    la_comids = set(riv["COMID"].tolist())
    la_vals = out.loc[out["COMID"].isin(la_comids), "corridor_impervious"].dropna()
    if len(la_vals):
        print(f"\nSanity LA River ({len(la_vals)} reaches): median={la_vals.median():.3f} "
              f"(expect >0.6)")
    else:
        print("\nSanity LA River: no matching MERIT reaches found in bbox")

    riv_mt = pyogrio.read_dataframe(
        paths.MERIT_RIV, columns=["COMID", "COMID"],
        bbox=(-115.0, 46.5, -114.0, 47.5),  # Rural Montana (Bitterroot Valley)
    )
    mt_comids = set(riv_mt["COMID"].tolist())
    mt_vals = out.loc[out["COMID"].isin(mt_comids), "corridor_impervious"].dropna()
    if len(mt_vals):
        print(f"Sanity rural Montana ({len(mt_vals)} reaches): median={mt_vals.median():.3f} "
              f"(expect <0.05)")
    else:
        print("Sanity rural Montana: no matching MERIT reaches found in bbox")


if __name__ == "__main__":
    main()
