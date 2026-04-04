use crate::bbox::BBox;
use crate::cell::Cell;
use crate::coord::Coord;
use crate::coverage;
use crate::floodfill::{FloodFill, EXTERIOR, FILLABLE, INTERIOR};
use crate::grid::{Bounded, Grid, Infinite};
use crate::measures;
use crate::side::Side;
use geo::Contains;
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

/// Compute cell-center coverage of each grid cell by a polygon.
///
/// For each cell in the polygon's bounding box, tests whether the cell
/// center falls inside the polygon. Cells whose center is inside get
/// coverage = 1.0; all others get 0.0. No fractional coverage is computed.
///
/// This matches the rasterization approach used by xarray-spatial and
/// other scanline-based zonal statistics tools.
pub fn raster_cell_center(
    raster_grid: &Grid<Bounded>,
    exterior_ring: &[Coord],
    interior_rings: &[Vec<Coord>],
) -> CoverageResult {
    let geom_bbox = ring_bbox(exterior_ring);
    let sub_grid = raster_grid.shrink_to_fit(&geom_bbox);
    if sub_grid.is_empty() {
        return CoverageResult {
            fractions: Array2::zeros((0, 0)),
            grid: sub_grid,
        };
    }

    let geo_polygon = build_geo_polygon(exterior_ring, interior_rings);
    let mut results = Array2::zeros((sub_grid.rows(), sub_grid.cols()));

    for row in 0..sub_grid.rows() {
        let y = sub_grid.y_for_row(row);
        for col in 0..sub_grid.cols() {
            let x = sub_grid.x_for_col(col);
            if geo_polygon.contains(&geo_types::Point(geo_types::Coord { x, y })) {
                results[[row, col]] = 1.0;
            }
        }
    }

    CoverageResult {
        fractions: results,
        grid: sub_grid,
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

    // --- raster_cell_center tests ---

    #[test]
    fn test_center_full_coverage() {
        // Polygon covers entire 4x4 grid — all cell centers inside.
        let grid = make_grid(0.0, 0.0, 4.0, 4.0, 1.0, 1.0);
        let exterior = vec![
            Coord::new(0.0, 0.0),
            Coord::new(4.0, 0.0),
            Coord::new(4.0, 4.0),
            Coord::new(0.0, 4.0),
            Coord::new(0.0, 0.0),
        ];

        let result = raster_cell_center(&grid, &exterior, &[]);
        assert_eq!(result.fractions.dim(), (4, 4));
        for row in 0..4 {
            for col in 0..4 {
                assert_eq!(
                    result.fractions[[row, col]], 1.0,
                    "cell [{row},{col}] should be 1.0"
                );
            }
        }
    }

    #[test]
    fn test_center_half_coverage() {
        // Polygon covers left half of a 4x4 grid.
        // shrink_to_fit clips the result to the polygon bbox [0,2]x[0,4] → 2 cols, 4 rows.
        // geo::Contains is strict (boundary excluded), so cell centers at
        // x=0.5 are inside but x=1.5 is also inside (< 2.0 boundary).
        let grid = make_grid(0.0, 0.0, 4.0, 4.0, 1.0, 1.0);
        let exterior = vec![
            Coord::new(0.0, 0.0),
            Coord::new(2.0, 0.0),
            Coord::new(2.0, 4.0),
            Coord::new(0.0, 4.0),
            Coord::new(0.0, 0.0),
        ];

        let result = raster_cell_center(&grid, &exterior, &[]);

        // Sub-grid is clipped to polygon bbox: 2 cols × 4 rows
        assert_eq!(result.fractions.dim(), (4, 2));

        for row in 0..4 {
            // col 0 center at x=0.5 → strictly inside
            assert_eq!(result.fractions[[row, 0]], 1.0, "row {row} col 0");
            // col 1 center at x=1.5 → strictly inside (boundary at x=2.0)
            assert_eq!(result.fractions[[row, 1]], 1.0, "row {row} col 1");
        }
    }

    #[test]
    fn test_center_binary_only() {
        // A polygon covering 50% of a cell should still produce 1.0
        // (center is inside) or 0.0 (center is outside) — never 0.5.
        let grid = make_grid(0.0, 0.0, 2.0, 2.0, 1.0, 1.0);
        // Small rectangle: covers top-left quadrant of cell (0,0).
        // Cell center at (0.5, 1.5) is inside this box → 1.0.
        let exterior = vec![
            Coord::new(0.0, 1.0),
            Coord::new(1.0, 1.0),
            Coord::new(1.0, 2.0),
            Coord::new(0.0, 2.0),
            Coord::new(0.0, 1.0),
        ];

        let result = raster_cell_center(&grid, &exterior, &[]);

        for row in 0..result.fractions.dim().0 {
            for col in 0..result.fractions.dim().1 {
                let f = result.fractions[[row, col]];
                assert!(
                    f == 0.0 || f == 1.0,
                    "cell [{row},{col}] = {f}, expected 0.0 or 1.0"
                );
            }
        }
    }

    #[test]
    fn test_center_empty_intersection() {
        // Polygon outside grid → empty result.
        let grid = make_grid(0.0, 0.0, 4.0, 4.0, 1.0, 1.0);
        let exterior = vec![
            Coord::new(10.0, 10.0),
            Coord::new(12.0, 10.0),
            Coord::new(12.0, 12.0),
            Coord::new(10.0, 12.0),
            Coord::new(10.0, 10.0),
        ];

        let result = raster_cell_center(&grid, &exterior, &[]);
        assert_eq!(result.fractions.dim(), (0, 0));
    }

    #[test]
    fn test_center_triangle() {
        // Right triangle: (0,0)-(4,0)-(0,4). On a 4x4 grid, cell centers
        // below the hypotenuse (y < 4-x) are inside.
        let grid = make_grid(0.0, 0.0, 4.0, 4.0, 1.0, 1.0);
        let exterior = vec![
            Coord::new(0.0, 0.0),
            Coord::new(4.0, 0.0),
            Coord::new(0.0, 4.0),
            Coord::new(0.0, 0.0),
        ];

        let result = raster_cell_center(&grid, &exterior, &[]);

        // Cell centers: (col+0.5, ymax - (row+0.5))
        // Hypotenuse: y = 4 - x, so inside when y < 4-x (strictly).
        // Row 0 (y=3.5): x<0.5 → col 0 only? No: 3.5 < 4-0.5=3.5 is false (boundary).
        // geo::Contains is strict: boundary points are NOT inside.
        // Row 3 (y=0.5): x<3.5 → cols 0,1,2 inside.
        let total: f32 = result.fractions.iter().sum();
        // Triangle area = 8, cells inside ≈ 6 (interior cells only, boundary excluded)
        assert!(
            total >= 3.0 && total <= 10.0,
            "total = {total}, expected between 3 and 10"
        );
        // Every cell is binary
        for &f in result.fractions.iter() {
            assert!(f == 0.0 || f == 1.0, "got {f}");
        }
    }

    #[test]
    fn test_center_with_hole() {
        // Square with a hole: exterior 0..6, hole 2..4.
        // On a 6x6 grid, cell centers inside the hole should be 0.0.
        let grid = make_grid(0.0, 0.0, 6.0, 6.0, 1.0, 1.0);
        let exterior = vec![
            Coord::new(0.0, 0.0),
            Coord::new(6.0, 0.0),
            Coord::new(6.0, 6.0),
            Coord::new(0.0, 6.0),
            Coord::new(0.0, 0.0),
        ];
        let hole = vec![
            Coord::new(2.0, 2.0),
            Coord::new(4.0, 2.0),
            Coord::new(4.0, 4.0),
            Coord::new(2.0, 4.0),
            Coord::new(2.0, 2.0),
        ];

        let result = raster_cell_center(&grid, &exterior, &[hole]);

        // All cells should be 1.0 except the 4 cells whose centers
        // fall inside the hole (x=2.5,3.5 and y=2.5,3.5).
        // Grid rows: row 0 → y=5.5, row 5 → y=0.5
        // Hole y=2..4 means rows with center y in (2,4): y=2.5 (row 3), y=3.5 (row 2)
        // Hole x=2..4 means cols with center x in (2,4): x=2.5 (col 2), x=3.5 (col 3)
        let total: f32 = result.fractions.iter().sum();
        assert_eq!(total, 32.0, "36 cells - 4 hole cells = 32");

        // Verify hole cells are 0
        assert_eq!(result.fractions[[2, 2]], 0.0, "hole cell [2,2]");
        assert_eq!(result.fractions[[2, 3]], 0.0, "hole cell [2,3]");
        assert_eq!(result.fractions[[3, 2]], 0.0, "hole cell [3,2]");
        assert_eq!(result.fractions[[3, 3]], 0.0, "hole cell [3,3]");

        // Verify non-hole cells are 1
        assert_eq!(result.fractions[[0, 0]], 1.0);
        assert_eq!(result.fractions[[5, 5]], 1.0);
    }

    #[test]
    fn test_center_agrees_with_exact_on_interior() {
        // For a polygon that aligns with the grid, both methods should
        // agree: interior cells get 1.0 in both.
        let grid = make_grid(0.0, 0.0, 4.0, 4.0, 1.0, 1.0);
        let exterior = vec![
            Coord::new(0.0, 0.0),
            Coord::new(4.0, 0.0),
            Coord::new(4.0, 4.0),
            Coord::new(0.0, 4.0),
            Coord::new(0.0, 0.0),
        ];

        let exact = raster_cell_intersection(&grid, &exterior, &[]);
        let center = raster_cell_center(&grid, &exterior, &[]);

        assert_eq!(exact.fractions.dim(), center.fractions.dim());
        for row in 0..4 {
            for col in 0..4 {
                assert_eq!(
                    exact.fractions[[row, col]], center.fractions[[row, col]],
                    "mismatch at [{row},{col}]"
                );
            }
        }
    }

    // --- raster_cell_intersection tests (original) ---

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
