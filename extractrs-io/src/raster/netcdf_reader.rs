use crate::raster::grid_builder::build_grid_from_coords;
use crate::IoError;
use extractrs_core::prelude::{BBox, Bounded, Grid};
use ndarray::Array2;
use std::path::Path;

/// A raster read from a Daymet NetCDF file.
pub struct DaymetRaster {
    /// 2D data array (rows, cols) — row 0 is northernmost.
    pub data: Array2<f64>,
    /// Fill/nodata value.
    pub nodata: f64,
    /// Grid in LCC coordinates (meters).
    pub grid: Grid<Bounded>,
}

/// Read a variable from a Daymet NetCDF file at a specific time index.
///
/// Reads the full spatial extent. For large files, use `read_daymet_var_bbox` instead.
pub fn read_daymet_var(
    nc_path: &Path,
    var_name: &str,
    time_index: usize,
) -> Result<DaymetRaster, IoError> {
    let file = netcdf::open(nc_path)?;

    let x_var = file
        .variable("x")
        .ok_or_else(|| IoError::Other("no 'x' variable in NetCDF".into()))?;
    let y_var = file
        .variable("y")
        .ok_or_else(|| IoError::Other("no 'y' variable in NetCDF".into()))?;

    let x_vals: Vec<f64> = read_coord_as_f64(&x_var)?;
    let y_vals: Vec<f64> = read_coord_as_f64(&y_var)?;

    let ny = y_vals.len();
    let nx = x_vals.len();
    let grid = build_grid_from_coords(&x_vals, &y_vals)?;

    let data_var = file
        .variable(var_name)
        .ok_or_else(|| IoError::Other(format!("no '{}' variable in NetCDF", var_name)))?;

    let nodata = get_nodata(&data_var);

    let data = read_slice(&data_var, var_name, time_index, 0, ny, 0, nx)?;

    Ok(DaymetRaster { data, nodata, grid })
}

/// Read a spatial subset of a Daymet variable, given an LCC bounding box.
///
/// This avoids loading the entire continental raster (~63M pixels) into memory.
pub fn read_daymet_var_bbox(
    nc_path: &Path,
    var_name: &str,
    time_index: usize,
    lcc_bbox: &BBox,
) -> Result<DaymetRaster, IoError> {
    let file = netcdf::open(nc_path)?;

    let x_var = file
        .variable("x")
        .ok_or_else(|| IoError::Other("no 'x' variable in NetCDF".into()))?;
    let y_var = file
        .variable("y")
        .ok_or_else(|| IoError::Other("no 'y' variable in NetCDF".into()))?;

    let x_all: Vec<f64> = read_coord_as_f64(&x_var)?;
    let y_all: Vec<f64> = read_coord_as_f64(&y_var)?;

    let dx = (x_all[1] - x_all[0]).abs();
    let dy = (y_all[0] - y_all[1]).abs();

    // Find index ranges for the bbox (with 1-cell padding)
    let col_start = x_all
        .iter()
        .position(|&x| x + dx / 2.0 > lcc_bbox.xmin)
        .unwrap_or(0)
        .saturating_sub(1);
    let col_end = x_all
        .iter()
        .rposition(|&x| x - dx / 2.0 < lcc_bbox.xmax)
        .map(|i| (i + 2).min(x_all.len()))
        .unwrap_or(x_all.len());

    // y is descending: y_all[0] is northernmost
    let row_start = y_all
        .iter()
        .position(|&y| y - dy / 2.0 < lcc_bbox.ymax)
        .unwrap_or(0)
        .saturating_sub(1);
    let row_end = y_all
        .iter()
        .rposition(|&y| y + dy / 2.0 > lcc_bbox.ymin)
        .map(|i| (i + 2).min(y_all.len()))
        .unwrap_or(y_all.len());

    let x_sub: Vec<f64> = x_all[col_start..col_end].to_vec();
    let y_sub: Vec<f64> = y_all[row_start..row_end].to_vec();

    if x_sub.len() < 2 || y_sub.len() < 2 {
        return Err(IoError::InvalidGrid(format!(
            "bbox doesn't intersect raster: x[{}..{}] y[{}..{}]",
            col_start, col_end, row_start, row_end
        )));
    }

    let ny = y_sub.len();
    let nx = x_sub.len();
    let grid = build_grid_from_coords(&x_sub, &y_sub)?;

    let data_var = file
        .variable(var_name)
        .ok_or_else(|| IoError::Other(format!("no '{}' variable in NetCDF", var_name)))?;

    let nodata = get_nodata(&data_var);

    let data = read_slice(&data_var, var_name, time_index, row_start, ny, col_start, nx)?;

    Ok(DaymetRaster { data, nodata, grid })
}

/// Read a coordinate variable, handling both f32 and f64 storage.
fn read_coord_as_f64(var: &netcdf::Variable) -> Result<Vec<f64>, IoError> {
    // Try f64 first, fall back to f32
    match var.get_values::<f64, _>(..) {
        Ok(vals) => Ok(vals),
        Err(_) => {
            let f32_vals: Vec<f32> = var
                .get_values(..)
                .map_err(|e| IoError::Other(format!("reading coordinate: {}", e)))?;
            Ok(f32_vals.iter().map(|&v| v as f64).collect())
        }
    }
}

fn get_nodata(data_var: &netcdf::Variable) -> f64 {
    data_var
        .attribute_value("_FillValue")
        .and_then(|r| r.ok())
        .and_then(|v| match v {
            netcdf::AttributeValue::Float(f) => Some(f as f64),
            netcdf::AttributeValue::Double(d) => Some(d),
            netcdf::AttributeValue::Short(s) => Some(s as f64),
            _ => None,
        })
        .unwrap_or(-9999.0)
}

fn read_slice(
    data_var: &netcdf::Variable,
    var_name: &str,
    time_index: usize,
    row_start: usize,
    ny: usize,
    col_start: usize,
    nx: usize,
) -> Result<Array2<f64>, IoError> {
    let dims = data_var.dimensions();
    if dims.len() == 3 {
        let flat: Vec<f32> = data_var
            .get_values([
                time_index..time_index + 1,
                row_start..row_start + ny,
                col_start..col_start + nx,
            ])
            .map_err(|e| IoError::Other(format!("reading {}: {}", var_name, e)))?;
        Array2::from_shape_vec((ny, nx), flat.iter().map(|&v| v as f64).collect())
            .map_err(|e| IoError::Other(format!("shape error: {}", e)))
    } else if dims.len() == 2 {
        let flat: Vec<f32> = data_var
            .get_values([row_start..row_start + ny, col_start..col_start + nx])
            .map_err(|e| IoError::Other(format!("reading {}: {}", var_name, e)))?;
        Array2::from_shape_vec((ny, nx), flat.iter().map(|&v| v as f64).collect())
            .map_err(|e| IoError::Other(format!("shape error: {}", e)))
    } else {
        Err(IoError::Other(format!(
            "expected 2D or 3D variable, got {}D",
            dims.len()
        )))
    }
}
