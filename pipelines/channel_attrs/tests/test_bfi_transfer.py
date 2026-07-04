"""Unit tests for bfi_transfer focused on the drainage-density guard.

The weighted_transfer core is tested in test_transfer.py.
The load_bfi helper is tested via the zipfile path below.
"""
import io
import zipfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipelines.channel_attrs.bfi_transfer import compute_drainage_density, load_bfi


# ---------------------------------------------------------------------------
# compute_drainage_density guard: catchsize <= 0 -> NaN
# ---------------------------------------------------------------------------

class _FakeNC:
    """Minimal netCDF4.Dataset stand-in for testing."""

    def __init__(self, comids, catchsizes):
        self._data = {"COMID": np.array(comids, dtype="int64"),
                      "catchsize": np.array(catchsizes, dtype="float64")}

    def __getitem__(self, key):
        return self._data[key]

    def close(self):
        pass


def _make_riv_parquet(tmp_path: Path, comids: list[int], lengthkms: list[float]) -> Path:
    import geopandas as gpd
    from shapely.geometry import LineString

    # Build a minimal shapefile-like GeoDataFrame and write as parquet.
    # pyogrio reads shapefiles, so we need to write a temp shapefile.
    import tempfile, os
    shp = tmp_path / "riv.shp"
    gdf = gpd.GeoDataFrame(
        {"COMID": comids, "lengthkm": lengthkms},
        geometry=[LineString([(0, i), (1, i)]) for i in range(len(comids))],
        crs="EPSG:4326",
    )
    gdf.to_file(shp)
    return shp


def test_drainage_density_guard_zero_catchsize(tmp_path, monkeypatch):
    """catchsize == 0 must yield NaN drainage_density."""
    import geopandas as gpd
    from shapely.geometry import LineString
    import netCDF4

    # Minimal shapefile with 3 reaches.
    comids = [1, 2, 3]
    lengths = [2.0, 3.0, 1.0]
    catchsizes = [4.0, 0.0, -1.0]  # COMID 2 and 3 are bad (zero / negative)

    shp = tmp_path / "riv.shp"
    gpd.GeoDataFrame(
        {"COMID": comids, "lengthkm": lengths},
        geometry=[LineString([(0, i), (1, i)]) for i in range(3)],
        crs="EPSG:4326",
    ).to_file(shp)

    nc_path = tmp_path / "fake.nc"
    ds_w = netCDF4.Dataset(nc_path, "w")
    ds_w.createDimension("COMID", 3)
    v_comid = ds_w.createVariable("COMID", "i8", ("COMID",))
    v_comid[:] = comids
    v_cs = ds_w.createVariable("catchsize", "f8", ("COMID",))
    v_cs[:] = catchsizes
    ds_w.close()

    result = compute_drainage_density(shp, nc_path)

    row1 = result.loc[result["COMID"] == 1, "drainage_density"].item()
    row2 = result.loc[result["COMID"] == 2, "drainage_density"].item()
    row3 = result.loc[result["COMID"] == 3, "drainage_density"].item()

    assert abs(row1 - 2.0 / 4.0) < 1e-9, f"Expected 0.5, got {row1}"
    assert np.isnan(row2), f"Expected NaN for catchsize=0, got {row2}"
    assert np.isnan(row3), f"Expected NaN for catchsize<0, got {row3}"


def test_drainage_density_positive():
    """Normal case: lengthkm / catchsize."""
    import geopandas as gpd
    from shapely.geometry import LineString
    import netCDF4, tempfile

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        shp = tmp_path / "riv.shp"
        gpd.GeoDataFrame(
            {"COMID": [10, 20], "lengthkm": [5.0, 10.0]},
            geometry=[LineString([(0, 0), (1, 0)]), LineString([(0, 1), (1, 1)])],
            crs="EPSG:4326",
        ).to_file(shp)

        nc_path = tmp_path / "fake.nc"
        ds_w = netCDF4.Dataset(nc_path, "w")
        ds_w.createDimension("COMID", 2)
        v_c = ds_w.createVariable("COMID", "i8", ("COMID",))
        v_c[:] = [10, 20]
        v_cs = ds_w.createVariable("catchsize", "f8", ("COMID",))
        v_cs[:] = [2.5, 5.0]
        ds_w.close()

        result = compute_drainage_density(shp, nc_path)

    assert abs(result.loc[result["COMID"] == 10, "drainage_density"].item() - 2.0) < 1e-9
    assert abs(result.loc[result["COMID"] == 20, "drainage_density"].item() - 2.0) < 1e-9


# ---------------------------------------------------------------------------
# load_bfi: zip-reading and negative-COMID filtering
# ---------------------------------------------------------------------------

def _make_bfi_zip(tmp_path: Path, records: list[dict]) -> Path:
    csv_content = "COMID,CAT_BFI,NODATA\n"
    for r in records:
        csv_content += f"{r['COMID']},{r['CAT_BFI']},0.0\n"
    zip_path = tmp_path / "BFI_CONUS.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("BFI_CONUS.txt", csv_content)
    return zip_path


def test_load_bfi_filters_negative_comids(tmp_path):
    zip_path = _make_bfi_zip(tmp_path, [
        {"COMID": -1000, "CAT_BFI": 30.0},
        {"COMID": 500, "CAT_BFI": 45.0},
        {"COMID": 600, "CAT_BFI": 70.0},
    ])
    df = load_bfi(zip_path)
    assert len(df) == 2
    assert set(df["foreign_id"].tolist()) == {500, 600}


def test_load_bfi_dtype_and_column_names(tmp_path):
    zip_path = _make_bfi_zip(tmp_path, [{"COMID": 100, "CAT_BFI": 55.5}])
    df = load_bfi(zip_path)
    assert set(df.columns) == {"foreign_id", "bfi"}
    assert df["foreign_id"].dtype == "int32"
    assert df["bfi"].item() == 55.5
