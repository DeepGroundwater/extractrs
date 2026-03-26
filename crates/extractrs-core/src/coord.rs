/// A 2D coordinate with exact floating-point comparison semantics.
///
/// Exact `==` comparison is intentional and critical for boundary detection
/// in the traversal algorithm.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coord {
    pub x: f64,
    pub y: f64,
}

impl Coord {
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    #[inline]
    pub fn distance(&self, other: &Coord) -> f64 {
        let dx = other.x - self.x;
        let dy = other.y - self.y;
        (dx * dx + dy * dy).sqrt()
    }

    #[inline]
    pub fn equals(&self, other: &Coord, tol: f64) -> bool {
        (other.x - self.x).abs() < tol && (other.y - self.y).abs() < tol
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance() {
        let a = Coord::new(0.0, 0.0);
        let b = Coord::new(3.0, 4.0);
        assert!((a.distance(&b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_exact_equality() {
        let a = Coord::new(1.0, 2.0);
        let b = Coord::new(1.0, 2.0);
        assert_eq!(a, b);
    }

    #[test]
    fn test_tolerance_equality() {
        let a = Coord::new(1.0, 2.0);
        let b = Coord::new(1.0 + 1e-10, 2.0 - 1e-10);
        assert!(a.equals(&b, 1e-8));
        assert!(!a.equals(&b, 1e-12));
    }
}
