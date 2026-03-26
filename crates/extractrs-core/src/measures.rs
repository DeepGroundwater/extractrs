use crate::coord::Coord;

/// Signed area of a ring via the shoelace formula.
///
/// Positive = counter-clockwise, negative = clockwise.
/// Uses first coordinate as reference point to reduce rounding error.
///
/// Port of exactextract's `area_signed()` (measures.cpp).
pub fn area_signed(ring: &[Coord]) -> f64 {
    if ring.len() < 3 {
        return 0.0;
    }
    let x0 = ring[0].x;
    let mut sum = 0.0;
    for i in 1..ring.len() - 1 {
        let x = ring[i].x - x0;
        let y_next = ring[i + 1].y;
        let y_prev = ring[i - 1].y;
        // Standard shoelace: positive = CCW in math convention (y-up)
        sum += x * (y_next - y_prev);
    }
    sum / 2.0
}

/// Absolute area of a ring.
#[inline]
pub fn area(ring: &[Coord]) -> f64 {
    area_signed(ring).abs()
}

/// Total length of a coordinate sequence.
pub fn length(coords: &[Coord]) -> f64 {
    let mut sum = 0.0;
    for i in 1..coords.len() {
        sum += coords[i - 1].distance(&coords[i]);
    }
    sum
}

/// Check if a coordinate sequence has at least two distinct points.
pub fn has_multiple_unique_coordinates(coords: &[Coord]) -> bool {
    for i in 1..coords.len() {
        if coords[0] != coords[i] {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_area_ccw_unit_square() {
        // CCW unit square
        let ring = vec![
            Coord::new(0.0, 0.0),
            Coord::new(1.0, 0.0),
            Coord::new(1.0, 1.0),
            Coord::new(0.0, 1.0),
            Coord::new(0.0, 0.0),
        ];
        let a = area_signed(&ring);
        assert!((a - 1.0).abs() < 1e-10, "expected +1.0, got {a}");
    }

    #[test]
    fn test_area_cw_unit_square() {
        // CW unit square (reversed)
        let ring = vec![
            Coord::new(0.0, 0.0),
            Coord::new(0.0, 1.0),
            Coord::new(1.0, 1.0),
            Coord::new(1.0, 0.0),
            Coord::new(0.0, 0.0),
        ];
        let a = area_signed(&ring);
        assert!((a + 1.0).abs() < 1e-10, "expected -1.0, got {a}");
    }

    #[test]
    fn test_area_triangle() {
        let ring = vec![
            Coord::new(0.0, 0.0),
            Coord::new(2.0, 0.0),
            Coord::new(1.0, 2.0),
            Coord::new(0.0, 0.0),
        ];
        assert!((area(&ring) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_length() {
        let coords = vec![
            Coord::new(0.0, 0.0),
            Coord::new(3.0, 0.0),
            Coord::new(3.0, 4.0),
        ];
        assert!((length(&coords) - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_degenerate_area() {
        assert_eq!(area_signed(&[]), 0.0);
        assert_eq!(area_signed(&[Coord::new(0.0, 0.0)]), 0.0);
        assert_eq!(
            area_signed(&[Coord::new(0.0, 0.0), Coord::new(1.0, 0.0)]),
            0.0
        );
    }

    #[test]
    fn test_has_multiple_unique() {
        let same = vec![Coord::new(1.0, 1.0), Coord::new(1.0, 1.0)];
        assert!(!has_multiple_unique_coordinates(&same));

        let diff = vec![Coord::new(1.0, 1.0), Coord::new(2.0, 1.0)];
        assert!(has_multiple_unique_coordinates(&diff));
    }
}
