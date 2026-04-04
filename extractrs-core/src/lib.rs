//! `extractrs-core` — Pure-Rust raster-polygon coverage engine.
//!
//! This crate computes how much of each raster grid cell is covered by
//! arbitrary polygons. It provides two coverage methods:
//!
//! - [`raster_cell_intersection`] — **Exact fractional coverage**. Computes
//!   the precise area of intersection between each cell and the polygon.
//!   Port of the algorithm from [exactextract](https://github.com/isciences/exactextract)
//!   with GEOS replaced by the `geo` crate.
//!
//! - [`raster_cell_center`] — **Cell-center point-in-polygon**. Tests whether
//!   each cell's center point falls inside the polygon (binary 1.0 / 0.0).
//!   Produces identical results to xarray-spatial and other scanline-based
//!   rasterization tools.
//!
//! # Exact Algorithm
//!
//! The exact algorithm proceeds in 4 phases:
//!
//! 1. **Ring Traversal**: Walk each polygon ring coordinate-by-coordinate,
//!    computing exact exit points via line-cell intersection.
//! 2. **Coverage Fraction**: For each touched cell, assemble closed rings
//!    from traversal chains + cell boundary, compute area via shoelace.
//! 3. **Flood Fill**: Classify untouched cells as interior/exterior via
//!    point-in-polygon test and scanline propagation.
//! 4. **Hole Subtraction**: Exterior rings add area, interior rings subtract.
//!
//! # Example
//!
//! ```
//! use extractrs_core::prelude::*;
//!
//! let grid = Grid::<Bounded>::new(BBox::new(0.0, 0.0, 10.0, 10.0), 1.0, 1.0);
//! let exterior = vec![
//!     Coord::new(2.0, 2.0),
//!     Coord::new(8.0, 2.0),
//!     Coord::new(8.0, 8.0),
//!     Coord::new(2.0, 8.0),
//!     Coord::new(2.0, 2.0),
//! ];
//!
//! // Exact fractional coverage
//! let result = raster_cell_intersection(&grid, &exterior, &[]);
//! assert_eq!(result.fractions.dim(), (6, 6));
//!
//! // Cell-center binary coverage
//! let result = raster_cell_center(&grid, &exterior, &[]);
//! assert_eq!(result.fractions.dim(), (6, 6));
//! ```

pub mod bbox;
pub mod cell;
pub mod coord;
pub mod coverage;
pub mod crossing;
pub mod floodfill;
pub mod grid;
pub mod intersection;
pub mod measures;
pub mod perimeter;
pub mod side;
pub mod traversal;

/// Convenient re-exports for common usage.
pub mod prelude {
    pub use crate::bbox::BBox;
    pub use crate::coord::Coord;
    pub use crate::grid::{Bounded, Grid, Infinite};
    pub use crate::intersection::{raster_cell_center, raster_cell_intersection, CoverageResult};
}
