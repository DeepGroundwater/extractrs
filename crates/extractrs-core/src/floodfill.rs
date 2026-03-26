use crate::grid::{Bounded, Grid};
use geo::Contains;
use geo_types::{Coord as GeoCoord, Point, Polygon};
use ndarray::Array2;
use std::collections::VecDeque;

/// Sentinel values for flood fill classification.
pub const EXTERIOR: f32 = 0.0;
pub const INTERIOR: f32 = 1.0;
pub const FILLABLE: f32 = -1.0;

/// Flood fill engine that classifies undetermined cells as interior
/// or exterior to a polygon.
///
/// Uses `geo::Contains` for point-in-polygon tests (replaces GEOS
/// PreparedContains in the C++ implementation).
///
/// Port of exactextract's `FloodFill` class.
pub struct FloodFill {
    polygon: Polygon<f64>,
    grid: Grid<Bounded>,
}

impl FloodFill {
    pub fn new(polygon: Polygon<f64>, grid: Grid<Bounded>) -> Self {
        Self { polygon, grid }
    }

    /// Test if a cell center is inside the polygon.
    pub fn cell_is_inside(&self, row: usize, col: usize) -> bool {
        let x = self.grid.x_for_col(col);
        let y = self.grid.y_for_row(row);
        self.polygon.contains(&Point(GeoCoord { x, y }))
    }

    /// Classify all FILLABLE cells as INTERIOR or EXTERIOR.
    ///
    /// For each FILLABLE cell encountered, performs a PIP test and then
    /// propagates the result to all connected FILLABLE cells using
    /// scanline flood fill.
    pub fn flood(&self, areas: &mut Array2<f32>) {
        let (nrows, ncols) = areas.dim();
        for i in 0..nrows {
            for j in 0..ncols {
                if areas[[i, j]] == FILLABLE {
                    let fill_value = if self.cell_is_inside(i, j) {
                        INTERIOR
                    } else {
                        EXTERIOR
                    };
                    flood_from_pixel(areas, i, j, fill_value);
                }
            }
        }
    }
}

/// BFS scanline flood fill from a seed pixel.
///
/// Fills all connected FILLABLE cells with `fill_value`.
/// Uses a queue-based approach (BFS) to avoid stack overflow
/// on large connected regions.
///
/// Port of exactextract's `flood_from_pixel()`.
fn flood_from_pixel(arr: &mut Array2<f32>, start_i: usize, start_j: usize, fill_value: f32) {
    let (nrows, ncols) = arr.dim();
    let mut queue = VecDeque::new();
    queue.push_back((start_i, start_j));

    while let Some((i, j)) = queue.pop_front() {
        if arr[[i, j]] == fill_value {
            continue; // Already filled
        }
        if arr[[i, j]] != FILLABLE {
            continue; // Not fillable (already determined)
        }

        // Check left neighbor
        if j > 0 && arr[[i, j - 1]] == FILLABLE {
            queue.push_back((i, j - 1));
        }

        let j0 = j;

        // Fill along this row to the right
        let mut jj = j;
        while jj < ncols && arr[[i, jj]] == FILLABLE {
            arr[[i, jj]] = fill_value;
            jj += 1;
        }
        let j1 = jj;

        // Initiate scanlines above
        if i > 0 {
            for jj in j0..j1 {
                if arr[[i - 1, jj]] == FILLABLE {
                    queue.push_back((i - 1, jj));
                }
            }
        }

        // Initiate scanlines below
        if i < nrows - 1 {
            for jj in j0..j1 {
                if arr[[i + 1, jj]] == FILLABLE {
                    queue.push_back((i + 1, jj));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bbox::BBox;
    use geo_types::{Coord as GeoCoord, LineString};

    fn make_polygon(coords: &[(f64, f64)]) -> Polygon<f64> {
        let exterior: LineString<f64> = coords
            .iter()
            .map(|(x, y)| GeoCoord { x: *x, y: *y })
            .collect();
        Polygon::new(exterior, vec![])
    }

    #[test]
    fn test_flood_simple() {
        let grid = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 3.0, 3.0), 1.0, 1.0);
        let poly = make_polygon(&[
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 3.0),
            (0.0, 3.0),
            (0.0, 0.0),
        ]);
        let ff = FloodFill::new(poly, grid);

        let mut areas = Array2::from_elem((3, 3), FILLABLE);
        ff.flood(&mut areas);

        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(areas[[i, j]], INTERIOR);
            }
        }
    }

    #[test]
    fn test_flood_partial() {
        // Polygon covers left half of a 4x2 grid.
        // In practice, the boundary column (col 2) would have a partial
        // coverage fraction from the traversal algorithm, acting as a
        // barrier between interior and exterior fillable regions.
        let grid = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 4.0, 2.0), 1.0, 1.0);
        let poly = make_polygon(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (0.0, 0.0),
        ]);
        let ff = FloodFill::new(poly, grid);

        // Pre-fill boundary column (col 2) with partial coverage as would
        // happen after the traversal phase — this separates the regions.
        let mut areas = Array2::from_elem((2, 4), FILLABLE);
        areas[[0, 2]] = 0.5; // boundary cell, already determined
        areas[[1, 2]] = 0.5; // boundary cell, already determined
        ff.flood(&mut areas);

        // Left 2 cols should be interior
        assert_eq!(areas[[0, 0]], INTERIOR);
        assert_eq!(areas[[0, 1]], INTERIOR);
        // Col 2 remains at its pre-filled value
        assert_eq!(areas[[0, 2]], 0.5);
        // Right col should be exterior
        assert_eq!(areas[[0, 3]], EXTERIOR);
    }
}
