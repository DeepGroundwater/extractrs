use crate::bbox::BBox;
use crate::cell::Cell;
use crate::coord::Coord;
use crate::coverage;
use crate::floodfill::{FloodFill, EXTERIOR, FILLABLE, INTERIOR};
use crate::grid::{Bounded, Grid, Infinite};
use crate::measures;
use crate::side::Side;
use geo_types::{Coord as GeoCoord, LineString, Polygon};
use ndarray::Array2;

/// Result of intersecting a polygon with a raster grid.
#[derive(Debug)]
pub struct CoverageResult {
    /// 2D array of coverage fractions in [0.0, 1.0].
    pub fractions: Array2<f32>,
    /// The sub-grid these fractions correspond to.
    pub grid: Grid<Bounded>,
}

/// Compute exact fractional coverage of each grid cell by a polygon.
///
/// The polygon is specified as an exterior ring and zero or more interior
/// rings (holes). All rings must be closed (first == last coordinate).
///
/// Returns a `CoverageResult` with coverage fractions and the sub-grid
/// that was actually computed (cropped to the polygon's bounding box).
///
/// This is the top-level entry point for the core algorithm.
///
/// Port of exactextract's `RasterCellIntersection`.
pub fn raster_cell_intersection(
    raster_grid: &Grid<Bounded>,
    exterior_ring: &[Coord],
    interior_rings: &[Vec<Coord>],
) -> CoverageResult {
    // Compute bounding box of the geometry
    let geom_bbox = ring_bbox(exterior_ring);

    // Shrink grid to geometry extent
    let sub_grid = raster_grid.shrink_to_fit(&geom_bbox);
    if sub_grid.is_empty() {
        return CoverageResult {
            fractions: Array2::zeros((0, 0)),
            grid: sub_grid,
        };
    }

    let mut results = Array2::zeros((sub_grid.rows(), sub_grid.cols()));

    // Check for rectangular ring fast path
    if interior_rings.is_empty()
        && exterior_ring.len() == 5
        && (measures::area(exterior_ring) - geom_bbox.area()).abs() < 1e-10
    {
        process_rectangular_ring(&sub_grid, &geom_bbox, true, &mut results);
        return CoverageResult {
            fractions: results,
            grid: sub_grid,
        };
    }

    // Build geo::Polygon for PIP queries in flood fill
    let geo_polygon = build_geo_polygon(exterior_ring, interior_rings);

    // Process exterior ring
    process_ring(
        &sub_grid,
        exterior_ring,
        true,  // exterior
        true,  // areal
        &geo_polygon,
        &mut results,
    );

    // Process interior rings (holes)
    for hole in interior_rings {
        process_ring(
            &sub_grid,
            hole,
            false, // interior
            true,
            &geo_polygon,
            &mut results,
        );
    }

    CoverageResult {
        fractions: results,
        grid: sub_grid,
    }
}

/// Process a single ring (exterior or interior) and accumulate results.
fn process_ring(
    sub_grid: &Grid<Bounded>,
    ring: &[Coord],
    exterior: bool,
    areal: bool,
    geo_polygon: &Polygon<f64>,
    results: &mut Array2<f32>,
) {
    if ring.len() < 4 && areal {
        return; // Degenerate
    }

    let ring_bbox = ring_bbox(ring);
    let ring_sub = sub_grid.shrink_to_fit(&ring_bbox);
    if ring_sub.is_empty() {
        return;
    }

    let inf_grid = ring_sub.to_infinite();

    // Check single-cell fast path
    if inf_grid.rows() == 3 && inf_grid.cols() == 3 {
        let cell_bb = inf_grid.cell_bbox(1, 1);
        if cell_bb.contains_box(&ring_bbox) {
            let row = sub_grid.get_row(ring_sub.extent().ymax - ring_sub.dy() * 0.5);
            let col = sub_grid.get_column(ring_sub.extent().xmin + ring_sub.dx() * 0.5);
            let cell_area = sub_grid.cell_bbox(row, col).area();
            if cell_area > 0.0 {
                let ring_area = measures::area(ring);
                let frac = (ring_area / cell_area) as f32;
                if exterior {
                    results[[row, col]] += frac;
                } else {
                    results[[row, col]] -= frac;
                }
            }
            return;
        }
    }

    // Ensure CCW orientation for exterior rings
    let mut coords: Vec<Coord> = ring.to_vec();
    if areal && measures::area_signed(&coords) < 0.0 {
        coords.reverse();
    }

    // Phase 1: Traverse cells
    let mut cells: Vec<Vec<Option<Cell>>> =
        (0..inf_grid.rows())
            .map(|_| (0..inf_grid.cols()).map(|_| None).collect())
            .collect();

    traverse_cells(&mut cells, &mut coords, &inf_grid, areal);

    // Phase 2-3: Collect coverage fractions with flood fill
    let finite_rows = inf_grid.rows() - 2;
    let finite_cols = inf_grid.cols() - 2;
    let mut areas = Array2::from_elem((finite_rows, finite_cols), FILLABLE);
    let finite_grid = inf_grid.to_bounded();
    let ff = FloodFill::new(geo_polygon.clone(), finite_grid);

    for i in 0..finite_rows {
        for j in 0..finite_cols {
            let inf_i = i + 1;
            let inf_j = j + 1;
            if let Some(ref cell) = cells[inf_i][inf_j] {
                if !cell.determined() {
                    areas[[i, j]] = if ff.cell_is_inside(i, j) {
                        INTERIOR
                    } else {
                        EXTERIOR
                    };
                } else {
                    // Compute covered fraction from traversals
                    let traversal_coords: Vec<&[Coord]> = cell
                        .traversals()
                        .iter()
                        .filter(|t| t.traversed() || t.is_closed_ring())
                        .filter(|t| t.multiple_unique_coordinates())
                        .map(|t| t.coords())
                        .collect();

                    if traversal_coords.is_empty() {
                        areas[[i, j]] = if ff.cell_is_inside(i, j) {
                            INTERIOR
                        } else {
                            EXTERIOR
                        };
                    } else {
                        let cell_bb = ring_sub.cell_bbox(i, j);
                        areas[[i, j]] = coverage::covered_fraction(&cell_bb, &traversal_coords);
                    }
                }
            }
        }
    }

    // Flood fill remaining FILLABLE cells
    ff.flood(&mut areas);

    // Phase 4: Accumulate into results
    let factor: f32 = if exterior { 1.0 } else { -1.0 };
    let row_offset = sub_grid.get_row(ring_sub.extent().ymax - ring_sub.dy() * 0.5);
    let col_offset = sub_grid.get_column(ring_sub.extent().xmin + ring_sub.dx() * 0.5);

    for i in 0..finite_rows {
        for j in 0..finite_cols {
            let ri = row_offset + i;
            let ci = col_offset + j;
            if ri < results.dim().0 && ci < results.dim().1 {
                results[[ri, ci]] += factor * areas[[i, j]];
            }
        }
    }
}

/// Core traversal loop: walk coordinates through the grid, recording
/// traversals in each cell they pass through.
///
/// Port of exactextract's `traverse_cells()`.
fn traverse_cells(
    cells: &mut [Vec<Option<Cell>>],
    coords: &mut Vec<Coord>,
    grid: &Grid<Infinite>,
    areal: bool,
) {
    if coords.is_empty() {
        return;
    }

    let mut pos: usize = 0;
    let mut row = grid.get_row(coords[0].y);
    let mut col = grid.get_column(coords[0].x);
    let mut last_exit: Option<Coord> = None;

    while pos < coords.len() {
        // Ensure cell exists
        if cells[row][col].is_none() {
            cells[row][col] = Some(Cell::new(grid.cell_bbox(row, col)));
        }
        let cell = cells[row][col].as_mut().unwrap();

        // Feed coordinates to this cell until it exits
        loop {
            if pos >= coords.len() && last_exit.is_none() {
                break;
            }

            let next_coord = if let Some(exit) = last_exit {
                exit
            } else {
                coords[pos]
            };

            let prev_original = if pos > 0 {
                Some(&coords[pos - 1])
            } else {
                None
            };

            let consumed = cell.take(&next_coord, prev_original);

            if cell.last_traversal().exited() {
                let exit_coord = *cell.last_traversal().exit_coordinate();
                if exit_coord != next_coord {
                    last_exit = Some(exit_coord);
                } else {
                    last_exit = None;
                    if consumed {
                        // The coordinate that triggered exit was consumed
                    } else {
                        // Don't advance pos — coordinate still needs processing
                    }
                }
                break;
            } else {
                if last_exit.is_some() {
                    last_exit = None;
                } else {
                    pos += 1;
                }
            }
        }

        cell.force_exit();

        if cell.last_traversal().exited() {
            // For areal geometries starting mid-cell, re-queue coordinates
            if areal && !cell.last_traversal().traversed() {
                let requeue: Vec<Coord> = cell.last_traversal().coords().to_vec();
                coords.extend_from_slice(&requeue);
            }

            // Move to adjacent cell based on exit side
            match cell.last_traversal().exit_side() {
                Side::Top => {
                    if row == 0 {
                        break;
                    }
                    row -= 1;
                }
                Side::Bottom => {
                    row += 1;
                    if row >= grid.rows() {
                        break;
                    }
                }
                Side::Left => {
                    if col == 0 {
                        break;
                    }
                    col -= 1;
                }
                Side::Right => {
                    col += 1;
                    if col >= grid.cols() {
                        break;
                    }
                }
                Side::None => {
                    // Closed ring within cell — advance to next coordinate
                    if last_exit.is_none() {
                        pos += 1;
                    }
                }
            }
        } else {
            // Cell didn't exit — we're done
            break;
        }
    }
}

/// Fast path for axis-aligned rectangular polygons.
///
/// Computes coverage by simple box-box intersection for each cell.
fn process_rectangular_ring(
    grid: &Grid<Bounded>,
    ring_bbox: &BBox,
    _exterior: bool,
    results: &mut Array2<f32>,
) {
    for row in 0..grid.rows() {
        for col in 0..grid.cols() {
            let cell_bb = grid.cell_bbox(row, col);
            if !cell_bb.intersects(ring_bbox) {
                continue;
            }
            let isect = cell_bb.intersection(ring_bbox);
            let frac = isect.area() / cell_bb.area();
            results[[row, col]] = frac as f32;
        }
    }
}

/// Compute bounding box of a coordinate ring.
fn ring_bbox(ring: &[Coord]) -> BBox {
    let mut bbox = BBox::make_empty();
    for c in ring {
        if c.x < bbox.xmin {
            bbox.xmin = c.x;
        }
        if c.x > bbox.xmax {
            bbox.xmax = c.x;
        }
        if c.y < bbox.ymin {
            bbox.ymin = c.y;
        }
        if c.y > bbox.ymax {
            bbox.ymax = c.y;
        }
    }
    bbox
}

/// Convert our Coord arrays to a geo_types::Polygon for PIP queries.
fn build_geo_polygon(exterior: &[Coord], interiors: &[Vec<Coord>]) -> Polygon<f64> {
    let ext_ls: LineString<f64> = exterior
        .iter()
        .map(|c| GeoCoord { x: c.x, y: c.y })
        .collect();
    let int_ls: Vec<LineString<f64>> = interiors
        .iter()
        .map(|ring| ring.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
        .collect();
    Polygon::new(ext_ls, int_ls)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_grid(xmin: f64, ymin: f64, xmax: f64, ymax: f64, dx: f64, dy: f64) -> Grid<Bounded> {
        Grid::<Bounded>::new(BBox::new(xmin, ymin, xmax, ymax), dx, dy)
    }

    #[test]
    fn test_rectangular_polygon_full_coverage() {
        let grid = make_grid(0.0, 0.0, 4.0, 4.0, 1.0, 1.0);
        let exterior = vec![
            Coord::new(0.0, 0.0),
            Coord::new(4.0, 0.0),
            Coord::new(4.0, 4.0),
            Coord::new(0.0, 4.0),
            Coord::new(0.0, 0.0),
        ];

        let result = raster_cell_intersection(&grid, &exterior, &[]);
        assert_eq!(result.fractions.dim(), (4, 4));

        for row in 0..4 {
            for col in 0..4 {
                assert!(
                    (result.fractions[[row, col]] - 1.0).abs() < 1e-5,
                    "cell [{row},{col}] = {}, expected 1.0",
                    result.fractions[[row, col]]
                );
            }
        }
    }

    #[test]
    fn test_rectangular_polygon_half_coverage() {
        let grid = make_grid(0.0, 0.0, 4.0, 4.0, 1.0, 1.0);
        // Rectangle covering left half
        let exterior = vec![
            Coord::new(0.0, 0.0),
            Coord::new(2.0, 0.0),
            Coord::new(2.0, 4.0),
            Coord::new(0.0, 4.0),
            Coord::new(0.0, 0.0),
        ];

        let result = raster_cell_intersection(&grid, &exterior, &[]);

        // Left 2 columns should be 1.0, right 2 should be 0.0
        for row in 0..4 {
            assert!(
                (result.fractions[[row, 0]] - 1.0).abs() < 1e-5,
                "row {row} col 0"
            );
            assert!(
                (result.fractions[[row, 1]] - 1.0).abs() < 1e-5,
                "row {row} col 1"
            );
        }
    }

    #[test]
    fn test_rectangular_partial_cell() {
        let grid = make_grid(0.0, 0.0, 2.0, 2.0, 1.0, 1.0);
        // Rectangle covering 50% of cell (0,0) (which is top-left, y=1..2)
        let exterior = vec![
            Coord::new(0.0, 1.0),
            Coord::new(0.5, 1.0),
            Coord::new(0.5, 2.0),
            Coord::new(0.0, 2.0),
            Coord::new(0.0, 1.0),
        ];

        let result = raster_cell_intersection(&grid, &exterior, &[]);

        // The top-left cell should be 50%
        assert!(
            (result.fractions[[0, 0]] - 0.5).abs() < 1e-5,
            "got {}",
            result.fractions[[0, 0]]
        );
    }

    #[test]
    fn test_empty_intersection() {
        let grid = make_grid(0.0, 0.0, 4.0, 4.0, 1.0, 1.0);
        // Polygon completely outside grid
        let exterior = vec![
            Coord::new(10.0, 10.0),
            Coord::new(12.0, 10.0),
            Coord::new(12.0, 12.0),
            Coord::new(10.0, 12.0),
            Coord::new(10.0, 10.0),
        ];

        let result = raster_cell_intersection(&grid, &exterior, &[]);
        assert_eq!(result.fractions.dim(), (0, 0));
    }

    #[test]
    fn test_triangle_polygon() {
        let grid = make_grid(0.0, 0.0, 2.0, 2.0, 1.0, 1.0);
        // Right triangle: (0,0)-(2,0)-(0,2)-(0,0) = area 2.0 over 4 cells
        let exterior = vec![
            Coord::new(0.0, 0.0),
            Coord::new(2.0, 0.0),
            Coord::new(0.0, 2.0),
            Coord::new(0.0, 0.0),
        ];

        let result = raster_cell_intersection(&grid, &exterior, &[]);

        let (rows, cols) = result.fractions.dim();
        let mut total = 0.0f32;
        for r in 0..rows {
            for c in 0..cols {
                let f = result.fractions[[r, c]];
                eprintln!("  [{r},{c}] = {f:.6}");
                total += f;
            }
        }
        eprintln!("  TOTAL = {total:.6} (expected 2.000)");
        eprintln!("  ERROR = {:.6}", (total - 2.0).abs());

        assert!(
            (total - 2.0).abs() < 0.01,
            "total coverage = {total}, expected 2.0"
        );
    }
}
