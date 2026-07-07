"""Length-weighted attribute transfer over reach-ID crosswalk tables.

A crosswalk row says: `part_len` meters of MERIT reach COMID run along
foreign reach `foreign_id`. Transfer = sum(part_len * value) / sum(part_len)
over matched rows with non-null values.
"""
import numpy as np
import pandas as pd


def weighted_transfer(
    xwalk: pd.DataFrame, attrs: pd.DataFrame, value_col: str
) -> pd.DataFrame:
    m = xwalk.merge(attrs[["foreign_id", value_col]], on="foreign_id", how="left")
    m = m[m[value_col].notna() & (m["part_len"] > 0)]
    if m.empty:
        return pd.DataFrame({"COMID": xwalk["COMID"].unique(), value_col: np.nan})
    m["_wv"] = m["part_len"] * m[value_col]
    g = m.groupby("COMID").agg(_wv=("_wv", "sum"), _w=("part_len", "sum"))
    out = (g["_wv"] / g["_w"]).rename(value_col).reset_index()
    return out.merge(
        pd.DataFrame({"COMID": xwalk["COMID"].unique()}), on="COMID", how="right"
    )
