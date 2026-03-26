//! `extractrs-core` — Pure-Rust exact fractional-pixel coverage fraction engine.
//!
//! This crate computes the exact fraction of each raster grid cell that is
//! covered by arbitrary polygons. It is a port of the core algorithm from
//! [exactextract](https://github.com/isciences/exactextract) with GEOS
//! replaced by the `geo` crate for point-in-polygon queries.
//!
//! # Algorithm
//!
//! The algorithm proceeds in 4 phases:
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
//! let result = raster_cell_intersection(&grid, &exterior, &[]);
//! // Result grid is shrunk to the polygon's bounding box (6x6 cells)
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
    pub use crate::intersection::{raster_cell_intersection, CoverageResult};
}
