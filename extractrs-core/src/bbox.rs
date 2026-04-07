use crate::coord::Coord;
use crate::crossing::Crossing;
use crate::side::Side;

/// Axis-aligned bounding box.
///
/// Port of exactextract's `Box` class. The `crossing()` method is the core
/// geometric operation for computing where a line exits a cell boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BBox {
    pub xmin: f64,
    pub ymin: f64,
    pub xmax: f64,
    pub ymax: f64,
}

impl BBox {
    #[inline]
    pub fn new(xmin: f64, ymin: f64, xmax: f64, ymax: f64) -> Self {
        Self {
            xmin,
            ymin,
            xmax,
            ymax,
        }
    }

    pub fn make_empty() -> Self {
        Self {
            xmin: f64::INFINITY,
            ymin: f64::INFINITY,
            xmax: f64::NEG_INFINITY,
            ymax: f64::NEG_INFINITY,
        }
    }

    #[inline]
    pub fn width(&self) -> f64 {
        self.xmax - self.xmin
    }

    #[inline]
    pub fn height(&self) -> f64 {
        self.ymax - self.ymin
    }

    #[inline]
    pub fn area(&self) -> f64 {
        self.width() * self.height()
    }

    #[inline]
    pub fn perimeter(&self) -> f64 {
        2.0 * self.width() + 2.0 * self.height()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.xmin >= self.xmax || self.ymin >= self.ymax
    }

    #[inline]
    pub fn upper_left(&self) -> Coord {
        Coord::new(self.xmin, self.ymax)
    }

    #[inline]
    pub fn upper_right(&self) -> Coord {
        Coord::new(self.xmax, self.ymax)
    }

    #[inline]
    pub fn lower_left(&self) -> Coord {
        Coord::new(self.xmin, self.ymin)
    }

    #[inline]
    pub fn lower_right(&self) -> Coord {
        Coord::new(self.xmax, self.ymin)
    }

    /// Determine which side of the box a coordinate lies on.
    /// Uses exact floating-point comparison — this is intentional.
    pub fn side(&self, c: &Coord) -> Side {
        if c.x == self.xmin {
            return Side::Left;
        }
        if c.x == self.xmax {
            return Side::Right;
        }
        if c.y == self.ymin {
            return Side::Bottom;
        }
        if c.y == self.ymax {
            return Side::Top;
        }
        Side::None
    }

    /// Inclusive containment test (boundary counts as inside).
    #[inline]
    pub fn contains(&self, c: &Coord) -> bool {
        c.x >= self.xmin && c.x <= self.xmax && c.y >= self.ymin && c.y <= self.ymax
    }

    /// Strict containment test (boundary counts as outside).
    #[inline]
    pub fn strictly_contains(&self, c: &Coord) -> bool {
        c.x > self.xmin && c.x < self.xmax && c.y > self.ymin && c.y < self.ymax
    }

    /// Test if this box contains another box entirely.
    pub fn contains_box(&self, b: &BBox) -> bool {
        b.xmin >= self.xmin && b.xmax <= self.xmax && b.ymin >= self.ymin && b.ymax <= self.ymax
    }

    pub fn intersects(&self, other: &BBox) -> bool {
        self.xmin <= other.xmax
            && self.xmax >= other.xmin
            && self.ymin <= other.ymax
            && self.ymax >= other.ymin
    }

    pub fn intersection(&self, other: &BBox) -> BBox {
        BBox::new(
            self.xmin.max(other.xmin),
            self.ymin.max(other.ymin),
            self.xmax.min(other.xmax),
            self.ymax.min(other.ymax),
        )
    }

    pub fn expand_to_include(&self, other: &BBox) -> BBox {
        BBox::new(
            self.xmin.min(other.xmin),
            self.ymin.min(other.ymin),
            self.xmax.max(other.xmax),
            self.ymax.max(other.ymax),
        )
    }

    pub fn translate(&self, dx: f64, dy: f64) -> BBox {
        BBox::new(
            self.xmin + dx,
            self.ymin + dy,
            self.xmax + dx,
            self.ymax + dy,
        )
    }

    /// Compute where a line segment from c1 to c2 exits this box.
    ///
    /// This is the core geometric algorithm. c1 must be inside the box.
    /// Uses quadrant-based projection with clamping for numerical robustness.
    ///
    /// Port of exactextract's `Box::crossing()` (box.cpp:43-118).
    ///
    /// **Degenerate fallback:** the original exactextract implementation
    /// throws "Never get here" when a vertical/horizontal segment lies
    /// entirely outside the box (i.e., the line never intersects this box).
    /// This can be triggered by low-quality polygon data with collinear
    /// vertices that happen to fall outside the current cell — see
    /// isciences/exactextract#126. Rather than panic, we return the box
    /// side closest to the segment so the traversal can make progress.
    /// Callers should also avoid invoking `crossing()` in this state when
    /// possible; see `Cell::take()` for the upstream guard.
    pub fn crossing(&self, c1: &Coord, c2: &Coord) -> Crossing {
        // Vertical line (x constant)
        if c1.x == c2.x {
            if c2.y >= self.ymax {
                return Crossing::new(Side::Top, c1.x, self.ymax);
            } else if c2.y <= self.ymin {
                return Crossing::new(Side::Bottom, c1.x, self.ymin);
            } else {
                // Degenerate: vertical line entirely outside the box on x.
                // Return the closest side so the traversal can advance.
                let y = clamp(c2.y, self.ymin, self.ymax);
                if c1.x <= self.xmin {
                    return Crossing::new(Side::Left, self.xmin, y);
                } else {
                    return Crossing::new(Side::Right, self.xmax, y);
                }
            }
        }

        // Horizontal line (y constant)
        if c1.y == c2.y {
            if c2.x >= self.xmax {
                return Crossing::new(Side::Right, self.xmax, c1.y);
            } else if c2.x <= self.xmin {
                return Crossing::new(Side::Left, self.xmin, c1.y);
            } else {
                // Degenerate: horizontal line entirely outside the box on y.
                let x = clamp(c2.x, self.xmin, self.xmax);
                if c1.y <= self.ymin {
                    return Crossing::new(Side::Bottom, x, self.ymin);
                } else {
                    return Crossing::new(Side::Top, x, self.ymax);
                }
            }
        }

        // Diagonal line: compute absolute slope and branch on quadrant
        let m = ((c2.y - c1.y) / (c2.x - c1.x)).abs();
        let up = c2.y > c1.y;
        let right = c2.x > c1.x;

        if up && right {
            // 1st quadrant (up-right)
            let y2 = c1.y + m * (self.xmax - c1.x);
            if y2 < self.ymax {
                Crossing::new(Side::Right, self.xmax, clamp(y2, self.ymin, self.ymax))
            } else {
                let x2 = c1.x + (self.ymax - c1.y) / m;
                Crossing::new(Side::Top, clamp(x2, self.xmin, self.xmax), self.ymax)
            }
        } else if up && !right {
            // 2nd quadrant (up-left)
            let y2 = c1.y + m * (c1.x - self.xmin);
            if y2 < self.ymax {
                Crossing::new(Side::Left, self.xmin, clamp(y2, self.ymin, self.ymax))
            } else {
                let x2 = c1.x - (self.ymax - c1.y) / m;
                Crossing::new(Side::Top, clamp(x2, self.xmin, self.xmax), self.ymax)
            }
        } else if !up && right {
            // 4th quadrant (down-right)
            let y2 = c1.y - m * (self.xmax - c1.x);
            if y2 > self.ymin {
                Crossing::new(Side::Right, self.xmax, clamp(y2, self.ymin, self.ymax))
            } else {
                let x2 = c1.x + (c1.y - self.ymin) / m;
                Crossing::new(Side::Bottom, clamp(x2, self.xmin, self.xmax), self.ymin)
            }
        } else {
            // 3rd quadrant (down-left)
            let y2 = c1.y - m * (c1.x - self.xmin);
            if y2 > self.ymin {
                Crossing::new(Side::Left, self.xmin, clamp(y2, self.ymin, self.ymax))
            } else {
                let x2 = c1.x - (c1.y - self.ymin) / m;
                Crossing::new(Side::Bottom, clamp(x2, self.xmin, self.xmax), self.ymin)
            }
        }
    }
}

/// Clamp a value to [lo, hi]. Prevents floating-point overflow from
/// diagonal line projections.
#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_box() -> BBox {
        BBox::new(0.0, 0.0, 1.0, 1.0)
    }

    #[test]
    fn test_dimensions() {
        let b = BBox::new(1.0, 2.0, 4.0, 6.0);
        assert_eq!(b.width(), 3.0);
        assert_eq!(b.height(), 4.0);
        assert_eq!(b.area(), 12.0);
    }

    #[test]
    fn test_side_detection() {
        let b = unit_box();
        assert_eq!(b.side(&Coord::new(0.0, 0.5)), Side::Left);
        assert_eq!(b.side(&Coord::new(1.0, 0.5)), Side::Right);
        assert_eq!(b.side(&Coord::new(0.5, 0.0)), Side::Bottom);
        assert_eq!(b.side(&Coord::new(0.5, 1.0)), Side::Top);
        assert_eq!(b.side(&Coord::new(0.5, 0.5)), Side::None);
    }

    #[test]
    fn test_crossing_vertical_up() {
        let b = unit_box();
        let c = b.crossing(&Coord::new(0.5, 0.5), &Coord::new(0.5, 2.0));
        assert_eq!(c.side, Side::Top);
        assert_eq!(c.coord, Coord::new(0.5, 1.0));
    }

    #[test]
    fn test_crossing_vertical_down() {
        let b = unit_box();
        let c = b.crossing(&Coord::new(0.5, 0.5), &Coord::new(0.5, -1.0));
        assert_eq!(c.side, Side::Bottom);
        assert_eq!(c.coord, Coord::new(0.5, 0.0));
    }

    #[test]
    fn test_crossing_horizontal_right() {
        let b = unit_box();
        let c = b.crossing(&Coord::new(0.5, 0.5), &Coord::new(2.0, 0.5));
        assert_eq!(c.side, Side::Right);
        assert_eq!(c.coord, Coord::new(1.0, 0.5));
    }

    #[test]
    fn test_crossing_diagonal_up_right() {
        let b = unit_box();
        // 45-degree line from center to (2,2): hits corner (1,1) exactly.
        // Algorithm projects to right edge first: y2 = 0.5 + 1*(1-0.5) = 1.0.
        // Since y2 == ymax (not strictly <), falls through to Top exit.
        let c = b.crossing(&Coord::new(0.5, 0.5), &Coord::new(2.0, 2.0));
        assert_eq!(c.side, Side::Top);
        assert!((c.coord.y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_crossing_vertical_degenerate_left_of_box() {
        // Vertical segment whose x is strictly left of the box.
        // Pre-fix this would panic with "vertical line doesn't exit box".
        let b = unit_box();
        let c = b.crossing(&Coord::new(-2.0, 0.5), &Coord::new(-2.0, 0.5));
        assert_eq!(c.side, Side::Left);
        assert_eq!(c.coord.x, 0.0);
    }

    #[test]
    fn test_crossing_vertical_degenerate_right_of_box() {
        let b = unit_box();
        let c = b.crossing(&Coord::new(5.0, 0.5), &Coord::new(5.0, 0.5));
        assert_eq!(c.side, Side::Right);
        assert_eq!(c.coord.x, 1.0);
    }

    #[test]
    fn test_crossing_horizontal_degenerate_above_box() {
        let b = unit_box();
        let c = b.crossing(&Coord::new(0.5, 5.0), &Coord::new(0.5, 5.0));
        assert_eq!(c.side, Side::Top);
        assert_eq!(c.coord.y, 1.0);
    }

    #[test]
    fn test_crossing_horizontal_degenerate_below_box() {
        let b = unit_box();
        let c = b.crossing(&Coord::new(0.5, -2.0), &Coord::new(0.5, -2.0));
        assert_eq!(c.side, Side::Bottom);
        assert_eq!(c.coord.y, 0.0);
    }

    #[test]
    fn test_crossing_diagonal_steep_up_right() {
        let b = unit_box();
        // Steep slope: exits top, not right
        let c = b.crossing(&Coord::new(0.5, 0.5), &Coord::new(0.6, 5.0));
        assert_eq!(c.side, Side::Top);
        assert_eq!(c.coord.y, 1.0);
    }

    #[test]
    fn test_intersects() {
        let a = BBox::new(0.0, 0.0, 2.0, 2.0);
        let b = BBox::new(1.0, 1.0, 3.0, 3.0);
        assert!(a.intersects(&b));

        let c = BBox::new(3.0, 3.0, 4.0, 4.0);
        assert!(!a.intersects(&c));
    }

    #[test]
    fn test_intersection() {
        let a = BBox::new(0.0, 0.0, 2.0, 2.0);
        let b = BBox::new(1.0, 1.0, 3.0, 3.0);
        let i = a.intersection(&b);
        assert_eq!(i, BBox::new(1.0, 1.0, 2.0, 2.0));
    }

    #[test]
    fn test_containment() {
        let b = unit_box();
        assert!(b.contains(&Coord::new(0.0, 0.0))); // boundary inclusive
        assert!(!b.strictly_contains(&Coord::new(0.0, 0.0))); // boundary exclusive
        assert!(b.strictly_contains(&Coord::new(0.5, 0.5)));
    }
}
