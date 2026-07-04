"""Buffered channel corridors from MERIT flowlines.

Buffer rationale (spec §3): MERIT lateral positional error is typically
100-300 m (Amatulli et al. 2022); 100 m half-width matches the StreamCat
precedent (Hill et al. 2016); reaches under broad GFPLAIN floodplains widen
to 200 m (flat-valley error tail).
"""
import time

import geopandas as gpd

from . import paths


def build_corridors(riv: gpd.GeoDataFrame, half_width_m: float) -> gpd.GeoDataFrame:
    riv = riv.to_crs(paths.CRS_EQUAL_AREA)
    out = riv[["COMID"]].copy()
    out["geometry"] = riv.geometry.buffer(half_width_m)
    return gpd.GeoDataFrame(out, geometry="geometry", crs=paths.CRS_EQUAL_AREA)


def main() -> None:
    print(f"Reading MERIT flowlines from {paths.MERIT_RIV} ...")
    t0 = time.time()
    riv = gpd.read_file(paths.MERIT_RIV, columns=["COMID"])
    print(f"  Read {len(riv):,} features in {time.time() - t0:.1f}s  (CRS: {riv.crs})")

    for name, hw in [
        ("corridors_100m", paths.CORRIDOR_HALF_WIDTH_M),
        ("corridors_wide", paths.CORRIDOR_WIDE_M),
    ]:
        t1 = time.time()
        print(f"Buffering {name} (half-width {hw:.0f} m) ...")
        corr = build_corridors(riv, hw)
        out_path = paths.DERIVED / f"{name}.parquet"
        corr.to_parquet(out_path)
        elapsed = time.time() - t1
        size_mb = out_path.stat().st_size / 1e6
        print(f"  {name}: {len(corr):,} corridors, CRS={corr.crs.to_epsg()}, "
              f"{size_mb:.1f} MB  ({elapsed:.1f}s)")


if __name__ == "__main__":
    main()
