use crate::bbox::BBox;
use crate::coord::Coord;

/// Compute the perimeter distance of a point on a box boundary.
///
/// Distance is measured counter-clockwise from the bottom-left corner:
///
/// ```text
///     (xmin,ymax)──────────(xmax,ymax)
///         │  ←── TOP ────▶  │
///    LEFT │                  │ RIGHT
///     ▲   │                  │  │
///     │   │                  │  ▼
///         │  ◀── BOTTOM ──  │
///     (xmin,ymin)──────────(xmax,ymin)
///         ▲ start here
/// ```
///
/// - LEFT side (x=xmin):   d = y - ymin
/// - TOP side (y=ymax):    d = height + (x - xmin)
/// - RIGHT side (x=xmax):  d = height + width + (ymax - y)
/// - BOTTOM side (y=ymin):  d = 2*height + width + (xmax - x)
///
/// Port of exactextract's `perimeter_distance()`.
pub fn perimeter_distance(bbox: &BBox, c: &Coord) -> f64 {
    let h = bbox.height();
    let w = bbox.width();

    if c.x == bbox.xmin {
        // Left side
        c.y - bbox.ymin
    } else if c.y == bbox.ymax {
        // Top side
        h + (c.x - bbox.xmin)
    } else if c.x == bbox.xmax {
        // Right side
        h + w + (bbox.ymax - c.y)
    } else if c.y == bbox.ymin {
        // Bottom side
        2.0 * h + w + (bbox.xmax - c.x)
    } else {
        panic!(
            "Cannot compute perimeter distance for point ({}, {}) not on boundary of {:?}",
            c.x, c.y, bbox
        );
    }
}

/// CCW distance from measure1 to measure2 along the perimeter.
///
/// If measure2 <= measure1, the distance is simply measure1 - measure2.
/// Otherwise, we wrap around: perimeter + measure1 - measure2.
#[inline]
pub fn perimeter_distance_ccw(measure1: f64, measure2: f64, total_perimeter: f64) -> f64 {
    if measure2 <= measure1 {
        measure1 - measure2
    } else {
        total_perimeter + measure1 - measure2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_box() -> BBox {
        BBox::new(0.0, 0.0, 1.0, 1.0)
    }

    #[test]
    fn test_corners() {
        let b = unit_box();
        assert!((perimeter_distance(&b, &Coord::new(0.0, 0.0)) - 0.0).abs() < 1e-10);
        assert!((perimeter_distance(&b, &Coord::new(0.0, 1.0)) - 1.0).abs() < 1e-10);
        assert!((perimeter_distance(&b, &Coord::new(1.0, 1.0)) - 2.0).abs() < 1e-10);
        assert!((perimeter_distance(&b, &Coord::new(1.0, 0.0)) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_midpoints() {
        let b = unit_box();
        // Left side midpoint
        assert!((perimeter_distance(&b, &Coord::new(0.0, 0.5)) - 0.5).abs() < 1e-10);
        // Top side midpoint
        assert!((perimeter_distance(&b, &Coord::new(0.5, 1.0)) - 1.5).abs() < 1e-10);
        // Right side midpoint
        assert!((perimeter_distance(&b, &Coord::new(1.0, 0.5)) - 2.5).abs() < 1e-10);
        // Bottom side midpoint
        assert!((perimeter_distance(&b, &Coord::new(0.5, 0.0)) - 3.5).abs() < 1e-10);
    }

    #[test]
    fn test_ccw_distance() {
        assert!((perimeter_distance_ccw(3.0, 1.0, 4.0) - 2.0).abs() < 1e-10);
        assert!((perimeter_distance_ccw(1.0, 3.0, 4.0) - 2.0).abs() < 1e-10);
        assert!((perimeter_distance_ccw(1.0, 1.0, 4.0) - 0.0).abs() < 1e-10);
    }
}
