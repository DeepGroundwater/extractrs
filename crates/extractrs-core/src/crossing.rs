use crate::coord::Coord;
use crate::side::Side;

/// A point where a line segment crosses a cell boundary.
#[derive(Debug, Clone, Copy)]
pub struct Crossing {
    pub side: Side,
    pub coord: Coord,
}

impl Crossing {
    #[inline]
    pub fn new(side: Side, x: f64, y: f64) -> Self {
        Self {
            side,
            coord: Coord::new(x, y),
        }
    }
}
