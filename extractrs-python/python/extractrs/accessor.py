"""xarray accessors for extractrs zonal statistics."""

import warnings

import numpy as np
import xarray as xr

from extractrs._extractrs import CoverageCache, build_cache, apply_stat, apply_stat_batch


def _extract_grid_params(da):
    """Extract grid extent and cell size from a DataArray's coordinates.

    Looks for 1D coordinate arrays along the spatial dimensions (y, x)
    and computes the grid bounding box and cell spacing from cell centers.
    """
    # Find spatial dimension names
    y_name = None
    x_name = None
    for name in da.dims:
        low = name.lower()
        if low in ("y", "lat", "latitude"):
            y_name = name
        elif low in ("x", "lon", "longitude"):
            x_name = name

    if y_name is None or x_name is None:
        raise ValueError(
            f"Cannot identify spatial dimensions. Found dims: {list(da.dims)}. "
            "Expected 'y'/'lat' and 'x'/'lon'."
        )

    y_vals = da[y_name].values.astype(np.float64)
    x_vals = da[x_name].values.astype(np.float64)

    if len(x_vals) < 2 or len(y_vals) < 2:
        raise ValueError("Spatial dimensions must have at least 2 values.")

    dx = abs(float(x_vals[1] - x_vals[0]))
    dy = abs(float(y_vals[1] - y_vals[0]))

    # Cell centers → cell edges (half-cell padding)
    xmin = float(x_vals.min()) - dx / 2
    xmax = float(x_vals.max()) + dx / 2
    ymin = float(y_vals.min()) - dy / 2
    ymax = float(y_vals.max()) + dy / 2

    return xmin, ymin, xmax, ymax, dx, dy


def _get_crs(obj):
    """Try to get CRS from an xarray object via rioxarray."""
    try:
        return obj.rio.crs
    except Exception:
        return None


def _prepare_geometries(gdf, obj, id_col):
    """Prepare WKB list and ID list from a GeoDataFrame.

    Handles CRS reprojection if rioxarray is available and CRS differs.
    """
    raster_crs = _get_crs(obj)
    geom_crs = gdf.crs

    if raster_crs is not None and geom_crs is not None:
        if not geom_crs.equals(raster_crs):
            warnings.warn(
                f"extractrs: Reprojecting geometries from {geom_crs} to {raster_crs}",
                UserWarning,
                stacklevel=3,
            )
            gdf = gdf.to_crs(raster_crs)
    elif raster_crs is not None and geom_crs is None:
        warnings.warn(
            "extractrs: GeoDataFrame has no CRS set. Assuming it matches the raster CRS.",
            UserWarning,
            stacklevel=3,
        )
    elif raster_crs is None and geom_crs is not None:
        # No rioxarray or no CRS on raster — can't check
        pass

    wkb_list = [geom.wkb for geom in gdf.geometry]

    if id_col is not None:
        id_list = gdf[id_col].astype(np.int64).tolist()
    else:
        id_list = list(range(len(gdf)))

    return wkb_list, id_list, gdf


def _find_data_vars(ds, var=None, vars=None):
    """Determine which data variables to process."""
    if var is not None:
        if var not in ds.data_vars:
            raise ValueError(f"Variable '{var}' not found. Available: {list(ds.data_vars)}")
        return [var]
    if vars is not None:
        for v in vars:
            if v not in ds.data_vars:
                raise ValueError(f"Variable '{v}' not found. Available: {list(ds.data_vars)}")
        return list(vars)
    return list(ds.data_vars)


def _find_time_dim(da):
    """Find the time dimension name, or None if there isn't one."""
    for name in da.dims:
        if name.lower() in ("time", "t"):
            return name
    return None


def _get_nodata(da):
    """Get nodata value from a DataArray.

    xarray decodes _FillValue to NaN and moves it from attrs to encoding,
    so we check encoding first to recover the original numeric sentinel.
    """
    # Check encoding first — xarray moves _FillValue here after decoding
    for key in ("_FillValue", "missing_value"):
        if key in da.encoding:
            val = da.encoding[key]
            if val is not None and np.isfinite(float(val)):
                return float(val)

    # Try rioxarray (may return NaN if data is NaN-filled)
    try:
        nd = da.rio.nodata
        if nd is not None and np.isfinite(float(nd)):
            return float(nd)
    except Exception:
        pass

    # Try attrs (pre-decode or manually set)
    for key in ("_FillValue", "missing_value"):
        if key in da.attrs:
            val = da.attrs[key]
            if val is not None and np.isfinite(float(val)):
                return float(val)

    # Default: use a large sentinel
    return -9999.0


def _run_zonal(da, cache, id_list, id_col, stat):
    """Run zonal stats on a DataArray, returning an xarray DataArray."""
    nodata = _get_nodata(da)
    time_dim = _find_time_dim(da)
    id_coord_name = id_col if id_col is not None else "basin"

    if time_dim is not None:
        # 3D: (time, y, x) → apply_stat_batch → (time, basin)
        data = da.values.astype(np.float64)
        if data.ndim != 3:
            raise ValueError(
                f"Expected 3D array (time, y, x), got {data.ndim}D. "
                "Select a single variable first."
            )
        result = apply_stat_batch(cache, data, nodata, stat)
        return xr.DataArray(
            result,
            dims=[time_dim, id_coord_name],
            coords={
                time_dim: da[time_dim].values,
                id_coord_name: id_list,
            },
        )
    else:
        # 2D: (y, x) → apply_stat → (basin,)
        data = da.values.astype(np.float64)
        if data.ndim != 2:
            raise ValueError(f"Expected 2D array (y, x), got {data.ndim}D.")
        result = apply_stat(cache, data, nodata, stat)
        return xr.DataArray(
            np.asarray(result),
            dims=[id_coord_name],
            coords={id_coord_name: id_list},
        )


@xr.register_dataset_accessor("extrs")
class ExtractrsDatasetAccessor:
    """xarray Dataset accessor for extractrs zonal statistics."""

    def __init__(self, ds):
        self._ds = ds

    def zonal_stats(self, gdf, stat="mean", id_col=None, var=None, vars=None, method="exact"):
        """Compute zonal statistics for each geometry in a GeoDataFrame.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Polygons to compute statistics for.
        stat : str
            Statistic to compute: "mean", "sum", "count", "min", "max",
            "variance", "stdev".
        id_col : str, optional
            Column in gdf to use as basin IDs. If None, uses integer index.
        var : str, optional
            Single data variable to process.
        vars : list of str, optional
            Subset of data variables to process. If neither var nor vars
            is given, all data variables are processed.
        method : str, default "exact"
            Coverage method: "exact" for fractional pixel coverage,
            "center" for cell-center point-in-polygon (binary 0/1).

        Returns
        -------
        xarray.Dataset
            Dataset with spatial dims replaced by basin/id dimension.
        """
        data_vars = _find_data_vars(self._ds, var=var, vars=vars)

        # Use the first data variable to extract grid params and build cache
        ref_da = self._ds[data_vars[0]]
        xmin, ymin, xmax, ymax, dx, dy = _extract_grid_params(ref_da)
        wkb_list, id_list, gdf_proj = _prepare_geometries(gdf, ref_da, id_col)
        cache = build_cache(wkb_list, id_list, xmin, ymin, xmax, ymax, dx, dy, method)

        results = {}
        for vname in data_vars:
            da = self._ds[vname]
            results[vname] = _run_zonal(da, cache, id_list, id_col, stat)

        return xr.Dataset(results)


@xr.register_dataarray_accessor("extrs")
class ExtractrsDataArrayAccessor:
    """xarray DataArray accessor for extractrs zonal statistics."""

    def __init__(self, da):
        self._da = da

    def zonal_stats(self, gdf, stat="mean", id_col=None, method="exact"):
        """Compute zonal statistics for each geometry in a GeoDataFrame.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            Polygons to compute statistics for.
        stat : str
            Statistic to compute: "mean", "sum", "count", "min", "max",
            "variance", "stdev".
        id_col : str, optional
            Column in gdf to use as basin IDs. If None, uses integer index.
        method : str, default "exact"
            Coverage method: "exact" for fractional pixel coverage,
            "center" for cell-center point-in-polygon (binary 0/1).

        Returns
        -------
        xarray.DataArray
            DataArray with spatial dims replaced by basin/id dimension.
        """
        xmin, ymin, xmax, ymax, dx, dy = _extract_grid_params(self._da)
        wkb_list, id_list, gdf_proj = _prepare_geometries(gdf, self._da, id_col)
        cache = build_cache(wkb_list, id_list, xmin, ymin, xmax, ymax, dx, dy, method)
        return _run_zonal(self._da, cache, id_list, id_col, stat)
