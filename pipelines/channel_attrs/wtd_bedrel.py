"""Bed-relative channel head and losing-stream fraction.

``channel_wtd_bed_rel = wtd_channel_zs − bankfull_depth``

Sign convention (preserved through Phase C):
  - **Positive** = water table is BELOW the channel bed = losing-possible.
  - **Negative** = water table is above the channel bed = gaining (perennial).

``losing_fraction`` is the per-reach fraction of densified channel vertices
where ``wtd_zs > bankfull_depth`` (i.e., water table below bed at that
point).  It comes from the ``frac_below_zs`` column saved by Task 6 in
``derived/wtd_channel.parquet``; no re-sampling of the raster is needed.

Output: ``derived/wtd_bedrel.parquet``
  Columns: COMID, channel_wtd_bed_rel, losing_fraction

Spec-mandated regime checks (three COMIDs, recorded in the report):
  - Ogallala/High-Plains reach: expect strongly positive bed-relative WTD.
  - Appalachian perennial reach: expect negative bed-relative WTD.
  - LA River: positive WTD + high imperviousness > 0.6 (Phase-C falsification
    pair — impervious-sealed channel where the high WTD depth is driven by
    canal lining, not natural losing-stream dynamics).
"""
import numpy as np
import pandas as pd

from . import paths

# Spec-mandated regime-check COMIDs (identified from MERIT shapefile geometry +
# empirical ZS WTD sampling — verified against derived wtd_channel.parquet).
#
# Ogallala: headwater ephemeral reach in southern High Plains (TX/KS border,
# lat ~37°N lon ~-102°W), order 1, wtd_channel_zs=63.8 m, bankfull_depth=1.15 m
# → bed_rel = +62.6 m (strongly losing-possible, deep Ogallala drawdown).
COMID_OGALLALA = 74049945
# Appalachian: New River headwaters (WV), order 6, uparea ~104k km²
# → bed_rel = -7.87 m (strongly gaining — perennial, as expected).
COMID_APPALACHIAN = 74037194
# LA River: urban tributary, lat ~34.0°N lon ~-118.2°W, order 1,
# wtd_channel_zs=3.57 m, bankfull_depth=2.21 m → bed_rel = +1.37 m (losing),
# corridor_impervious=0.71 > 0.6 — Phase-C falsification pair.
COMID_LA_RIVER = 77030730


def bed_relative(
    wtd: pd.DataFrame,
    bankfull: pd.DataFrame,
) -> pd.DataFrame:
    """Compute bed-relative channel head from WTD and bankfull depth.

    Parameters
    ----------
    wtd
        DataFrame with columns ``COMID`` and ``wtd_channel_zs`` (positive
        below-surface metres after the Task-6 sign conversion).
    bankfull
        DataFrame with columns ``COMID`` and ``bankfull_depth`` (metres).

    Returns
    -------
    pd.DataFrame
        Columns: ``COMID``, ``channel_wtd_bed_rel``.
        Positive = water table below bed = losing-possible.
    """
    merged = wtd[["COMID", "wtd_channel_zs"]].merge(
        bankfull[["COMID", "bankfull_depth"]], on="COMID", how="left"
    )
    merged["channel_wtd_bed_rel"] = merged["wtd_channel_zs"] - merged["bankfull_depth"]
    return merged[["COMID", "channel_wtd_bed_rel"]]


def main() -> None:
    # ── Load inputs ───────────────────────────────────────────────────────────
    print("Reading wtd_channel.parquet ...")
    wtd = pd.read_parquet(paths.DERIVED / "wtd_channel.parquet")
    print(f"  {len(wtd):,} rows, columns: {wtd.columns.tolist()}")

    print("Reading bankfull.parquet ...")
    bnk = pd.read_parquet(paths.DERIVED / "bankfull.parquet", columns=["COMID", "bankfull_depth"])
    print(f"  {len(bnk):,} reaches with bankfull depth  (median={bnk['bankfull_depth'].median():.3f} m)")

    # ── bed-relative head ─────────────────────────────────────────────────────
    bedrel = bed_relative(wtd, bnk)
    bedrel["COMID"] = bedrel["COMID"].astype(int)

    # ── losing_fraction from Task-6 per-vertex frac_below ────────────────────
    if "frac_below_zs" in wtd.columns:
        frac = wtd[["COMID", "frac_below_zs"]].rename(
            columns={"frac_below_zs": "losing_fraction"}
        )
        frac["COMID"] = frac["COMID"].astype(int)
        out = bedrel.merge(frac, on="COMID", how="left")
        n_frac = out["losing_fraction"].notna().sum()
        print(f"  losing_fraction from Task-6 frac_below_zs: {n_frac:,} reaches")
    else:
        print(
            "  frac_below_zs not found in wtd_channel.parquet — "
            "losing_fraction will be NaN (re-run wtd_sample.main() to populate)"
        )
        out = bedrel.copy()
        out["losing_fraction"] = np.nan

    # ── Assertions on sign convention ─────────────────────────────────────────
    bedrel_finite = out["channel_wtd_bed_rel"].dropna()
    losing_frac_finite = out["losing_fraction"].dropna()

    assert len(bedrel_finite) > 0, "channel_wtd_bed_rel is entirely NaN"
    pct_positive = float(np.mean(bedrel_finite > 0) * 100)
    pct_negative = float(np.mean(bedrel_finite < 0) * 100)
    print(
        f"\nBed-relative distribution:  "
        f"p50={bedrel_finite.median():.2f} m  "
        f"positive={pct_positive:.1f}%  negative={pct_negative:.1f}%"
    )

    if len(losing_frac_finite) > 0:
        losing_global = float(np.mean(losing_frac_finite > 0.5) * 100)
        print(f"Losing fraction:  median={losing_frac_finite.median():.3f}  "
              f"majority-losing (>0.5): {losing_global:.1f}%")

    # ── Spec-mandated regime checks ───────────────────────────────────────────
    print("\nRegime checks:")
    checks = {
        "Ogallala/High-Plains": COMID_OGALLALA,
        "Appalachian perennial": COMID_APPALACHIAN,
        "LA River": COMID_LA_RIVER,
    }
    for label, comid in checks.items():
        row = out[out["COMID"] == comid]
        if row.empty:
            print(f"  {label} (COMID={comid}): NOT FOUND in output")
            continue
        brel = row["channel_wtd_bed_rel"].item()
        frac = row["losing_fraction"].item() if "losing_fraction" in row.columns else np.nan
        print(f"  {label} (COMID={comid}): channel_wtd_bed_rel={brel:.2f} m  "
              f"losing_fraction={frac:.3f}")

    # Ogallala should be strongly positive (deep water table, >10 m typical)
    ogallala_row = out[out["COMID"] == COMID_OGALLALA]
    if not ogallala_row.empty and pd.notna(ogallala_row["channel_wtd_bed_rel"].item()):
        val = ogallala_row["channel_wtd_bed_rel"].item()
        if val <= 0:
            print(
                f"  WARNING: Ogallala reach channel_wtd_bed_rel={val:.2f} m is non-positive; "
                "expected strongly positive for High-Plains losing stream."
            )

    # Appalachian should be negative (high rainfall, shallow WT)
    appal_row = out[out["COMID"] == COMID_APPALACHIAN]
    if not appal_row.empty and pd.notna(appal_row["channel_wtd_bed_rel"].item()):
        val = appal_row["channel_wtd_bed_rel"].item()
        if val >= 0:
            print(
                f"  WARNING: Appalachian reach channel_wtd_bed_rel={val:.2f} m is non-negative; "
                "expected negative for perennial gaining stream."
            )

    # LA River: pair with corridor_impervious > 0.6 (Phase-C falsification)
    la_row = out[out["COMID"] == COMID_LA_RIVER]
    if not la_row.empty:
        la_brel = la_row["channel_wtd_bed_rel"].item()
        try:
            imp = pd.read_parquet(
                paths.DERIVED / "corridor_impervious.parquet",
                columns=["COMID", "corridor_impervious"]
            )
            la_imp = imp[imp["COMID"] == COMID_LA_RIVER]
            la_imp_val = la_imp["corridor_impervious"].item() if not la_imp.empty else np.nan
        except Exception:
            la_imp_val = np.nan
        print(
            f"  LA River COMID={COMID_LA_RIVER}:  channel_wtd_bed_rel={la_brel:.2f} m  "
            f"corridor_impervious={la_imp_val:.2f}  "
            f"(Phase-C pair: expect positive WTD + impervious>0.6)"
        )
        if pd.notna(la_imp_val) and la_imp_val < 0.4:
            print(
                "  WARNING: LA River corridor_impervious below 0.4 — "
                "check COMID selection or StreamCat transfer."
            )

    # ── Write output ──────────────────────────────────────────────────────────
    out_path = paths.DERIVED / "wtd_bedrel.parquet"
    out.to_parquet(out_path, index=False)

    n_brel = out["channel_wtd_bed_rel"].notna().sum()
    n_frac_out = out["losing_fraction"].notna().sum()
    print(
        f"\nWrote {out_path}  ({len(out):,} rows, {out_path.stat().st_size/1e6:.1f} MB)"
        f"\n  channel_wtd_bed_rel:  {n_brel:,} non-NaN"
        f"\n  losing_fraction:      {n_frac_out:,} non-NaN"
    )


if __name__ == "__main__":
    main()
