use crate::bbox::BBox;

/// Marker trait for grid extent modes.
pub trait ExtentMode: std::fmt::Debug + Clone {
    const PADDING: usize;
}

/// Bounded grid: no padding cells. Coordinates must be within extent.
#[derive(Debug, Clone)]
pub struct Bounded;
impl ExtentMode for Bounded {
    const PADDING: usize = 0;
}

/// Infinite grid: 1 padding cell on each side. Coordinates outside the
/// extent land in sentinel cells extending to +/- infinity.
#[derive(Debug, Clone)]
pub struct Infinite;
impl ExtentMode for Infinite {
    const PADDING: usize = 1;
}

/// A regular rectangular grid with uniform cell spacing.
///
/// Port of exactextract's `Grid<extent_tag>` template.
///
/// Row 0 is at the top (y = ymax). Y increases downward in row indices.
#[derive(Debug, Clone)]
pub struct Grid<E: ExtentMode> {
    extent: BBox,
    dx: f64,
    dy: f64,
    num_rows: usize,
    num_cols: usize,
    _marker: std::marker::PhantomData<E>,
}

impl<E: ExtentMode> Grid<E> {
    pub fn new(extent: BBox, dx: f64, dy: f64) -> Self {
        let data_rows = if extent.ymax > extent.ymin {
            (((extent.ymax - extent.ymin) / dy).round()) as usize
        } else {
            0
        };
        let data_cols = if extent.xmax > extent.xmin {
            (((extent.xmax - extent.xmin) / dx).round()) as usize
        } else {
            0
        };
        Self {
            extent,
            dx,
            dy,
            num_rows: 2 * E::PADDING + data_rows,
            num_cols: 2 * E::PADDING + data_cols,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn empty(dx: f64, dy: f64) -> Self {
        Self {
            extent: BBox::make_empty(),
            dx,
            dy,
            num_rows: 2 * E::PADDING,
            num_cols: 2 * E::PADDING,
            _marker: std::marker::PhantomData,
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.num_rows <= 2 * E::PADDING && self.num_cols <= 2 * E::PADDING
    }

    #[inline]
    pub fn rows(&self) -> usize {
        self.num_rows
    }

    #[inline]
    pub fn cols(&self) -> usize {
        self.num_cols
    }

    #[inline]
    pub fn size(&self) -> usize {
        self.num_rows * self.num_cols
    }

    #[inline]
    pub fn dx(&self) -> f64 {
        self.dx
    }

    #[inline]
    pub fn dy(&self) -> f64 {
        self.dy
    }

    #[inline]
    pub fn extent(&self) -> &BBox {
        &self.extent
    }

    /// X coordinate of the center of column `col`.
    #[inline]
    pub fn x_for_col(&self, col: usize) -> f64 {
        self.extent.xmin + ((col as f64) - (E::PADDING as f64) + 0.5) * self.dx
    }

    /// Y coordinate of the center of row `row`.
    #[inline]
    pub fn y_for_row(&self, row: usize) -> f64 {
        self.extent.ymax - ((row as f64) - (E::PADDING as f64) + 0.5) * self.dy
    }
}

// --- Bounded-specific methods ---

impl Grid<Bounded> {
    /// Convert to an infinite grid (adds 1 padding cell on each side).
    pub fn to_infinite(&self) -> Grid<Infinite> {
        Grid::<Infinite>::new(self.extent, self.dx, self.dy)
    }

    /// Get the column index for an x coordinate (bounded: must be within extent).
    pub fn get_column(&self, x: f64) -> usize {
        if x == self.extent.xmax {
            return self.num_cols - 1;
        }
        ((x - self.extent.xmin) / self.dx).floor() as usize
    }

    /// Get the row index for a y coordinate (bounded: must be within extent).
    pub fn get_row(&self, y: f64) -> usize {
        if y == self.extent.ymin {
            return self.num_rows - 1;
        }
        ((self.extent.ymax - y) / self.dy).floor() as usize
    }

    /// Get the bounding box of a cell at (row, col).
    ///
    /// Rightmost and bottommost cells extend to the exact grid extent
    /// to prevent floating-point gaps.
    pub fn cell_bbox(&self, row: usize, col: usize) -> BBox {
        let xmin = self.extent.xmin + col as f64 * self.dx;
        let ymax = self.extent.ymax - row as f64 * self.dy;

        let xmax = if col == self.num_cols - 1 {
            self.extent.xmax
        } else {
            self.extent.xmin + (col + 1) as f64 * self.dx
        };

        let ymin = if row == self.num_rows - 1 {
            self.extent.ymin
        } else {
            self.extent.ymax - (row + 1) as f64 * self.dy
        };

        BBox::new(xmin, ymin, xmax, ymax)
    }

    /// Shrink the grid to fit a bounding box, snapping to grid lines.
    pub fn shrink_to_fit(&self, b: &BBox) -> Grid<Bounded> {
        if self.is_empty() || !self.extent.intersects(b) {
            return Grid::<Bounded>::empty(self.dx, self.dy);
        }

        let col_min = ((b.xmin - self.extent.xmin) / self.dx).floor().max(0.0) as usize;
        let col_max =
            (((b.xmax - self.extent.xmin) / self.dx).ceil() as usize).min(self.num_cols);
        let row_min = ((self.extent.ymax - b.ymax) / self.dy).floor().max(0.0) as usize;
        let row_max =
            (((self.extent.ymax - b.ymin) / self.dy).ceil() as usize).min(self.num_rows);

        if col_min >= col_max || row_min >= row_max {
            return Grid::<Bounded>::empty(self.dx, self.dy);
        }

        let new_extent = BBox::new(
            self.extent.xmin + col_min as f64 * self.dx,
            self.extent.ymax - row_max as f64 * self.dy,
            self.extent.xmin + col_max as f64 * self.dx,
            self.extent.ymax - row_min as f64 * self.dy,
        );
        Grid::<Bounded>::new(new_extent, self.dx, self.dy)
    }

    /// Subdivide this grid into chunks of at most `max_cells` cells.
    pub fn subdivide(&self, max_cells: usize) -> Vec<Grid<Bounded>> {
        if self.size() <= max_cells || self.is_empty() {
            return vec![self.clone()];
        }

        let rows_per_chunk = (max_cells / self.num_cols).max(1);
        let mut chunks = Vec::new();
        let mut row = 0;

        while row < self.num_rows {
            let chunk_rows = rows_per_chunk.min(self.num_rows - row);
            let ymax = self.extent.ymax - row as f64 * self.dy;
            let ymin = self.extent.ymax - (row + chunk_rows) as f64 * self.dy;
            let chunk_extent = BBox::new(self.extent.xmin, ymin, self.extent.xmax, ymax);
            chunks.push(Grid::<Bounded>::new(chunk_extent, self.dx, self.dy));
            row += chunk_rows;
        }

        chunks
    }
}

// --- Infinite-specific methods ---

impl Grid<Infinite> {
    /// Convert to a bounded grid (removes padding cells).
    pub fn to_bounded(&self) -> Grid<Bounded> {
        Grid::<Bounded>::new(self.extent, self.dx, self.dy)
    }

    /// Get the column index for an x coordinate.
    /// Coordinates outside the extent land in padding cells (0 or cols-1).
    pub fn get_column(&self, x: f64) -> usize {
        if x < self.extent.xmin {
            return 0;
        }
        if x > self.extent.xmax {
            return self.num_cols - 1;
        }
        if x == self.extent.xmax {
            return self.num_cols - 2;
        }
        let col =
            Infinite::PADDING + ((x - self.extent.xmin) / self.dx).floor() as usize;
        col.min(self.get_column_max())
    }

    /// Get the row index for a y coordinate.
    /// Coordinates outside the extent land in padding cells (0 or rows-1).
    pub fn get_row(&self, y: f64) -> usize {
        if y > self.extent.ymax {
            return 0;
        }
        if y < self.extent.ymin {
            return self.num_rows - 1;
        }
        if y == self.extent.ymin {
            return self.num_rows - 2;
        }
        let row =
            Infinite::PADDING + ((self.extent.ymax - y) / self.dy).floor() as usize;
        row.min(self.get_row_max())
    }

    fn get_column_max(&self) -> usize {
        self.num_cols.saturating_sub(2)
    }

    fn get_row_max(&self) -> usize {
        self.num_rows.saturating_sub(2)
    }

    /// Get the bounding box of a cell at (row, col).
    /// Padding cells extend to +/- infinity.
    pub fn cell_bbox(&self, row: usize, col: usize) -> BBox {
        let xmin = if col == 0 {
            f64::NEG_INFINITY
        } else if col == self.num_cols - 1 {
            self.extent.xmax
        } else {
            self.extent.xmin + (col as f64 - 1.0) * self.dx
        };

        let xmax = if col == self.num_cols - 1 {
            f64::INFINITY
        } else if col == 0 {
            self.extent.xmin
        } else {
            self.extent.xmin + col as f64 * self.dx
        };

        let ymax = if row == 0 {
            f64::INFINITY
        } else if row == self.num_rows - 1 {
            self.extent.ymin
        } else {
            self.extent.ymax - (row as f64 - 1.0) * self.dy
        };

        let ymin = if row == self.num_rows - 1 {
            f64::NEG_INFINITY
        } else if row == 0 {
            self.extent.ymax
        } else {
            self.extent.ymax - row as f64 * self.dy
        };

        BBox::new(xmin, ymin, xmax, ymax)
    }

    /// Row offset between this grid and another infinite grid.
    pub fn row_offset(&self, other: &Grid<Infinite>) -> usize {
        ((other.extent.ymax - self.extent.ymax) / self.dy).round() as usize
    }

    /// Column offset between this grid and another infinite grid.
    pub fn col_offset(&self, other: &Grid<Infinite>) -> usize {
        ((self.extent.xmin - other.extent.xmin) / self.dx).round() as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bounded_basic() {
        let g = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 10.0, 10.0), 1.0, 1.0);
        assert_eq!(g.rows(), 10);
        assert_eq!(g.cols(), 10);
        assert!(!g.is_empty());
    }

    #[test]
    fn test_bounded_get_column() {
        let g = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 3.0, 3.0), 1.0, 1.0);
        assert_eq!(g.get_column(0.0), 0);
        assert_eq!(g.get_column(0.5), 0);
        assert_eq!(g.get_column(1.0), 1);
        assert_eq!(g.get_column(3.0), 2); // xmax maps to last col
    }

    #[test]
    fn test_bounded_get_row() {
        let g = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 3.0, 3.0), 1.0, 1.0);
        assert_eq!(g.get_row(3.0), 0); // ymax = top row
        assert_eq!(g.get_row(2.5), 0);
        assert_eq!(g.get_row(2.0), 1);
        assert_eq!(g.get_row(0.0), 2); // ymin = bottom row
    }

    #[test]
    fn test_bounded_cell_bbox() {
        let g = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 3.0, 3.0), 1.0, 1.0);
        let bb = g.cell_bbox(0, 0);
        assert_eq!(bb, BBox::new(0.0, 2.0, 1.0, 3.0));
    }

    #[test]
    fn test_bounded_cell_bbox_last() {
        // Last row/col should snap to extent boundary
        let g = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 3.0, 3.0), 1.0, 1.0);
        let bb = g.cell_bbox(2, 2);
        assert_eq!(bb.xmax, 3.0);
        assert_eq!(bb.ymin, 0.0);
    }

    #[test]
    fn test_infinite_padding() {
        let g = Grid::<Infinite>::new(BBox::new(0.0, 0.0, 3.0, 3.0), 1.0, 1.0);
        assert_eq!(g.rows(), 5); // 3 data + 2 padding
        assert_eq!(g.cols(), 5);
    }

    #[test]
    fn test_infinite_get_column() {
        let g = Grid::<Infinite>::new(BBox::new(0.0, 0.0, 3.0, 3.0), 1.0, 1.0);
        assert_eq!(g.get_column(-1.0), 0); // outside → padding
        assert_eq!(g.get_column(0.0), 1); // xmin → first data col
        assert_eq!(g.get_column(0.5), 1);
        assert_eq!(g.get_column(3.0), 3); // xmax → last data col
        assert_eq!(g.get_column(4.0), 4); // outside → padding
    }

    #[test]
    fn test_infinite_cell_bbox_padding() {
        let g = Grid::<Infinite>::new(BBox::new(0.0, 0.0, 2.0, 2.0), 1.0, 1.0);
        let left_pad = g.cell_bbox(1, 0);
        assert_eq!(left_pad.xmin, f64::NEG_INFINITY);
        assert_eq!(left_pad.xmax, 0.0);
    }

    #[test]
    fn test_shrink_to_fit() {
        let g = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 10.0, 10.0), 1.0, 1.0);
        let shrunk = g.shrink_to_fit(&BBox::new(2.5, 3.5, 5.5, 7.5));
        assert_eq!(shrunk.cols(), 4); // cols 2..6
        assert_eq!(shrunk.rows(), 5); // rows covering y 3..8
    }

    #[test]
    fn test_subdivide() {
        let g = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 10.0, 10.0), 1.0, 1.0);
        let chunks = g.subdivide(30);
        assert!(chunks.len() > 1);
        let total_cells: usize = chunks.iter().map(|c| c.size()).sum();
        assert_eq!(total_cells, 100);
    }

    #[test]
    fn test_center_coordinates() {
        let g = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 4.0, 4.0), 2.0, 2.0);
        assert!((g.x_for_col(0) - 1.0).abs() < 1e-10);
        assert!((g.x_for_col(1) - 3.0).abs() < 1e-10);
        assert!((g.y_for_row(0) - 3.0).abs() < 1e-10);
        assert!((g.y_for_row(1) - 1.0).abs() < 1e-10);
    }
}
