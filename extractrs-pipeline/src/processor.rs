use crate::coverage_cache::{CoverageCache, CoverageEntry};
use extractrs_io::raster::DaymetRaster;
use extractrs_stats::accumulator::{CellContribution, MeanAccum, StatAccumulator, StatValue};

/// Compute mean value for each basin using cached coverage fractions.
///
/// Returns (comid, mean_value) pairs. Basins with no valid data return NaN.
pub fn process_day(cache: &CoverageCache, raster: &DaymetRaster) -> Vec<(i64, f64)> {
    let (raster_rows, raster_cols) = raster.data.dim();
    let mut results = Vec::with_capacity(cache.entries.len());

    for entry in &cache.entries {
        let mean = compute_basin_mean(entry, raster, raster_rows, raster_cols);
        results.push((entry.comid, mean));
    }

    results
}

fn compute_basin_mean(
    entry: &CoverageEntry,
    raster: &DaymetRaster,
    raster_rows: usize,
    raster_cols: usize,
) -> f64 {
    let mut accum = MeanAccum::default();
    let (sub_rows, sub_cols) = entry.fractions.dim();

    for i in 0..sub_rows {
        let global_row = entry.row_offset + i;
        if global_row >= raster_rows {
            continue;
        }
        for j in 0..sub_cols {
            let global_col = entry.col_offset + j;
            if global_col >= raster_cols {
                continue;
            }

            let coverage = entry.fractions[[i, j]];
            if coverage <= 0.0 {
                continue;
            }

            let value = raster.data[[global_row, global_col]];
            let value_defined = (value - raster.nodata).abs() > 1e-6
                && value.is_finite();

            accum.process(&CellContribution {
                coverage,
                value,
                weight: 1.0,
                x: 0.0,
                y: 0.0,
                value_defined,
                weight_defined: true,
            });
        }
    }

    match accum.result() {
        StatValue::Float(v) => v,
        StatValue::OptFloat(Some(v)) => v,
        _ => f64::NAN,
    }
}
