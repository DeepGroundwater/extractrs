use extractrs_core::prelude::*;
use extractrs_stats::accumulator::*;
use extractrs_stats::variance::*;
use ndarray::Array2;
use numpy::{PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// A single polygon part: exterior ring + optional holes.
struct PolygonPart {
    exterior: Vec<Coord>,
    holes: Vec<Vec<Coord>>,
}

/// One basin's pre-computed coverage: sparse fractions + offset into the raster.
struct CacheEntry {
    id: i64,
    row_offset: usize,
    col_offset: usize,
    fractions: Array2<f32>,
}

/// Opaque coverage cache held on the Python side.
#[pyclass]
struct CoverageCache {
    entries: Vec<CacheEntry>,
    #[allow(dead_code)]
    grid: Grid<Bounded>,
}

/// Parse WKB bytes into polygon parts (exterior ring + holes).
/// Handles both Polygon and MultiPolygon WKB types.
fn parse_wkb(wkb: &[u8]) -> PyResult<Vec<PolygonPart>> {
    if wkb.len() < 5 {
        return Err(PyValueError::new_err("WKB too short"));
    }

    let big_endian = wkb[0] == 0;
    let geom_type = if big_endian {
        u32::from_be_bytes([wkb[1], wkb[2], wkb[3], wkb[4]])
    } else {
        u32::from_le_bytes([wkb[1], wkb[2], wkb[3], wkb[4]])
    };

    // Strip SRID flag and Z/M flags — only care about base type
    let base_type = geom_type & 0xFF;

    match base_type {
        3 => {
            // Polygon
            let part = parse_wkb_polygon(wkb, 5, big_endian)?;
            Ok(vec![part])
        }
        6 => {
            // MultiPolygon
            parse_wkb_multipolygon(wkb, big_endian)
        }
        _ => Err(PyValueError::new_err(format!(
            "Unsupported WKB geometry type: {}",
            geom_type
        ))),
    }
}

fn read_u32(data: &[u8], offset: usize, big_endian: bool) -> u32 {
    let bytes = [
        data[offset],
        data[offset + 1],
        data[offset + 2],
        data[offset + 3],
    ];
    if big_endian {
        u32::from_be_bytes(bytes)
    } else {
        u32::from_le_bytes(bytes)
    }
}

fn read_f64(data: &[u8], offset: usize, big_endian: bool) -> f64 {
    let bytes: [u8; 8] = data[offset..offset + 8].try_into().unwrap();
    if big_endian {
        f64::from_be_bytes(bytes)
    } else {
        f64::from_le_bytes(bytes)
    }
}

fn parse_ring(data: &[u8], offset: usize, big_endian: bool) -> PyResult<(Vec<Coord>, usize)> {
    let n_points = read_u32(data, offset, big_endian) as usize;
    let mut pos = offset + 4;
    let mut coords = Vec::with_capacity(n_points);
    for _ in 0..n_points {
        let x = read_f64(data, pos, big_endian);
        let y = read_f64(data, pos + 8, big_endian);
        coords.push(Coord::new(x, y));
        pos += 16;
    }
    Ok((coords, pos))
}

fn parse_wkb_polygon(data: &[u8], start: usize, big_endian: bool) -> PyResult<PolygonPart> {
    let n_rings = read_u32(data, start, big_endian) as usize;
    let mut pos = start + 4;

    if n_rings == 0 {
        return Ok(PolygonPart {
            exterior: vec![],
            holes: vec![],
        });
    }

    let (exterior, new_pos) = parse_ring(data, pos, big_endian)?;
    pos = new_pos;

    let mut holes = Vec::with_capacity(n_rings.saturating_sub(1));
    for _ in 1..n_rings {
        let (hole, new_pos) = parse_ring(data, pos, big_endian)?;
        pos = new_pos;
        holes.push(hole);
    }

    Ok(PolygonPart { exterior, holes })
}

fn parse_wkb_multipolygon(data: &[u8], big_endian: bool) -> PyResult<Vec<PolygonPart>> {
    let n_polys = read_u32(data, 5, big_endian) as usize;
    let mut pos = 9;
    let mut parts = Vec::with_capacity(n_polys);

    for _ in 0..n_polys {
        // Each sub-polygon has its own byte order + type header (5 bytes)
        let sub_be = data[pos] == 0;
        // Skip the 5-byte header
        pos += 5;
        let part = parse_wkb_polygon(data, pos, sub_be)?;

        // Advance past the polygon data
        let n_rings = read_u32(data, pos, sub_be) as usize;
        pos += 4;
        for _ in 0..n_rings {
            let n_pts = read_u32(data, pos, sub_be) as usize;
            pos += 4 + n_pts * 16;
        }

        if !part.exterior.is_empty() {
            parts.push(part);
        }
    }

    Ok(parts)
}

/// Create an accumulator from a stat name string.
fn make_accumulator(stat: &str) -> PyResult<Box<dyn StatAccumulator>> {
    match stat {
        "mean" => Ok(Box::new(MeanAccum::default())),
        "sum" => Ok(Box::new(SumAccum::default())),
        "count" => Ok(Box::new(CountAccum::default())),
        "min" => Ok(Box::new(MinAccum::default())),
        "max" => Ok(Box::new(MaxAccum::default())),
        "variance" => Ok(Box::new(VarianceAccum::default())),
        "stdev" => Ok(Box::new(StdevAccum::default())),
        _ => Err(PyValueError::new_err(format!(
            "Unknown stat '{}'. Supported: mean, sum, count, min, max, variance, stdev",
            stat
        ))),
    }
}

/// Extract a f64 result from a StatValue, returning NaN for None.
fn stat_to_f64(val: StatValue) -> f64 {
    match val {
        StatValue::Float(v) => v,
        StatValue::OptFloat(Some(v)) => v,
        StatValue::OptFloat(None) => f64::NAN,
        StatValue::Int(v) => v as f64,
        StatValue::FloatArray(_) => f64::NAN,
    }
}

#[pymethods]
impl CoverageCache {
    /// Number of basins in the cache.
    #[getter]
    fn n_basins(&self) -> usize {
        self.entries.len()
    }

    /// Basin IDs in cache order.
    fn ids<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        let ids: Vec<i64> = self.entries.iter().map(|e| e.id).collect();
        PyArray1::from_vec(py, ids)
    }
}

/// Build a coverage cache from WKB geometries and grid parameters.
///
/// Args:
///     wkb_list: list of WKB bytes (one per geometry)
///     id_list: list of integer IDs (same length as wkb_list)
///     xmin, ymin, xmax, ymax: grid extent
///     dx, dy: cell size
///
/// Returns:
///     CoverageCache object
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn build_cache(
    wkb_list: Vec<Vec<u8>>,
    id_list: Vec<i64>,
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
    dx: f64,
    dy: f64,
) -> PyResult<CoverageCache> {
    if wkb_list.len() != id_list.len() {
        return Err(PyValueError::new_err(
            "wkb_list and id_list must have the same length",
        ));
    }

    let grid = Grid::<Bounded>::new(BBox::new(xmin, ymin, xmax, ymax), dx, dy);
    let mut entries = Vec::with_capacity(wkb_list.len());

    for (wkb, id) in wkb_list.iter().zip(id_list.iter()) {
        let parts = parse_wkb(wkb)?;

        let mut part_results: Vec<(CoverageResult, usize, usize)> = Vec::new();

        for part in &parts {
            if part.exterior.is_empty() {
                continue;
            }

            let hole_slices: Vec<Vec<Coord>> = part.holes.clone();
            let result = raster_cell_intersection(&grid, &part.exterior, &hole_slices);

            if result.fractions.dim() == (0, 0) {
                continue;
            }

            let ro = grid.get_row(result.grid.extent().ymax - result.grid.dy() * 0.5);
            let co = grid.get_column(result.grid.extent().xmin + result.grid.dx() * 0.5);
            part_results.push((result, ro, co));
        }

        if part_results.is_empty() {
            continue;
        }

        // Single part fast path
        if part_results.len() == 1 {
            let (result, row_offset, col_offset) = part_results.into_iter().next().unwrap();
            entries.push(CacheEntry {
                id: *id,
                row_offset,
                col_offset,
                fractions: result.fractions,
            });
            continue;
        }

        // Multiple parts: merge into one array
        let mut min_row = usize::MAX;
        let mut min_col = usize::MAX;
        let mut max_row = 0usize;
        let mut max_col = 0usize;

        for (result, ro, co) in &part_results {
            let (rows, cols) = result.fractions.dim();
            min_row = min_row.min(*ro);
            min_col = min_col.min(*co);
            max_row = max_row.max(*ro + rows);
            max_col = max_col.max(*co + cols);
        }

        let total_rows = max_row - min_row;
        let total_cols = max_col - min_col;
        let mut combined = Array2::<f32>::zeros((total_rows, total_cols));

        for (result, ro, co) in &part_results {
            let (rows, cols) = result.fractions.dim();
            let dr = ro - min_row;
            let dc = co - min_col;
            for r in 0..rows {
                for c in 0..cols {
                    combined[[dr + r, dc + c]] += result.fractions[[r, c]];
                }
            }
        }

        entries.push(CacheEntry {
            id: *id,
            row_offset: min_row,
            col_offset: min_col,
            fractions: combined,
        });
    }

    Ok(CoverageCache { entries, grid })
}

/// Apply a statistic to a single 2D raster using the coverage cache.
///
/// Args:
///     cache: CoverageCache from build_cache()
///     raster: 2D numpy array (rows, cols) of raster values
///     nodata: nodata value to skip
///     stat: statistic name ("mean", "sum", "count", "min", "max", "variance", "stdev")
///
/// Returns:
///     1D numpy array of length n_basins with the stat result per basin
#[pyfunction]
fn apply_stat<'py>(
    py: Python<'py>,
    cache: &CoverageCache,
    raster: PyReadonlyArray2<'py, f64>,
    nodata: f64,
    stat: &str,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let raster = raster.as_array();
    let (raster_rows, raster_cols) = raster.dim();
    let nodata_is_nan = nodata.is_nan();
    let mut results = Vec::with_capacity(cache.entries.len());

    for entry in &cache.entries {
        let mut accum = make_accumulator(stat)?;
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

                let value = raster[[global_row, global_col]];
                let value_defined = if nodata_is_nan {
                    value.is_finite()
                } else {
                    (value - nodata).abs() > 1e-6 && value.is_finite()
                };

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

        results.push(stat_to_f64(accum.result()));
    }

    Ok(PyArray1::from_vec(py, results))
}

/// Apply a statistic to a 3D raster (time, rows, cols) using the coverage cache.
///
/// Args:
///     cache: CoverageCache from build_cache()
///     raster: 3D numpy array (time, rows, cols)
///     nodata: nodata value to skip
///     stat: statistic name
///
/// Returns:
///     2D numpy array (time, n_basins)
#[pyfunction]
fn apply_stat_batch<'py>(
    py: Python<'py>,
    cache: &CoverageCache,
    raster: numpy::PyReadonlyArray3<'py, f64>,
    nodata: f64,
    stat: &str,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let raster = raster.as_array();
    let (n_times, raster_rows, raster_cols) = raster.dim();
    let nodata_is_nan = nodata.is_nan();
    let n_basins = cache.entries.len();
    let mut results = Array2::<f64>::from_elem((n_times, n_basins), f64::NAN);

    for t in 0..n_times {
        let slice = raster.slice(ndarray::s![t, .., ..]);

        for (b, entry) in cache.entries.iter().enumerate() {
            let mut accum = make_accumulator(stat)?;
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

                    let value = slice[[global_row, global_col]];
                    let value_defined = if nodata_is_nan {
                        value.is_finite()
                    } else {
                        (value - nodata).abs() > 1e-6 && value.is_finite()
                    };

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

            results[[t, b]] = stat_to_f64(accum.result());
        }
    }

    Ok(PyArray2::from_owned_array(py, results))
}

/// extractrs: Fast exact zonal statistics — Rust engine.
#[pymodule]
fn _extractrs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<CoverageCache>()?;
    m.add_function(wrap_pyfunction!(build_cache, m)?)?;
    m.add_function(wrap_pyfunction!(apply_stat, m)?)?;
    m.add_function(wrap_pyfunction!(apply_stat_batch, m)?)?;
    Ok(())
}
