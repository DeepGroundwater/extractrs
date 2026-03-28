use crate::bbox::BBox;
use crate::coord::Coord;
use crate::measures;
use crate::side::Side;
use crate::traversal::Traversal;

/// A single grid cell that collects traversal information as polygon
/// rings pass through it.
///
/// Port of exactextract's `Cell` class.
#[derive(Debug)]
pub struct Cell {
    bbox: BBox,
    traversals: Vec<Traversal>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Location {
    Inside,
    Outside,
    Boundary,
}

/// Classify a coordinate relative to a bbox.
fn classify(bbox: &BBox, c: &Coord) -> Location {
    if bbox.strictly_contains(c) {
        Location::Inside
    } else if bbox.contains(c) {
        Location::Boundary
    } else {
        Location::Outside
    }
}

impl Cell {
    pub fn new(bbox: BBox) -> Self {
        Self {
            bbox,
            traversals: vec![Traversal::new()],
        }
    }

    #[inline]
    pub fn bbox(&self) -> &BBox {
        &self.bbox
    }

    #[inline]
    pub fn area(&self) -> f64 {
        self.bbox.area()
    }

    /// Ensure the last traversal is available for writing.
    /// If the last one is already exited or is a closed ring, push a new one.
    /// Returns the index of the active traversal.
    fn ensure_active_traversal(&mut self) -> usize {
        let last = self.traversals.last().unwrap();
        if last.exited() || last.is_closed_ring() {
            self.traversals.push(Traversal::new());
        }
        self.traversals.len() - 1
    }

    /// Feed a coordinate to this cell.
    ///
    /// Returns `true` if the coordinate was consumed (inside/boundary).
    /// Returns `false` if the coordinate was outside, triggering an exit.
    ///
    /// `prev_original`: the original (uninterpolated) previous coordinate,
    /// used for robustness when computing cell exit crossings.
    /// See exactextract regression test #7.
    pub fn take(&mut self, c: &Coord, prev_original: Option<&Coord>) -> bool {
        let idx = self.ensure_active_traversal();

        if self.traversals[idx].is_empty() {
            // Start a new traversal
            let side = self.bbox.side(c);
            self.traversals[idx].enter(*c, side);
            return true;
        }

        if classify(&self.bbox, c) != Location::Outside {
            // Coordinate is inside or on boundary
            self.traversals[idx].add(*c);

            if self.traversals[idx].is_closed_ring() {
                self.traversals[idx].force_exit(Side::None);
            }
            return true;
        }

        // Coordinate exits cell.
        // Use prev_original for crossing computation to avoid robustness issues
        // from interpolated coordinates (exactextract cell.cpp:126-129).
        let crossing = if let Some(prev) = prev_original {
            self.bbox.crossing(prev, c)
        } else {
            let last = *self.traversals[idx].last_coordinate();
            self.bbox.crossing(&last, c)
        };
        self.traversals[idx].exit(crossing.coord, crossing.side);

        false
    }

    /// Force the last traversal to exit if it's on a boundary.
    pub fn force_exit(&mut self) {
        let last_coord = {
            let t = self.traversals.last().unwrap();
            if t.exited() || t.is_empty() {
                return;
            }
            *t.last_coordinate()
        };
        if classify(&self.bbox, &last_coord) == Location::Boundary {
            let side = self.bbox.side(&last_coord);
            self.traversals.last_mut().unwrap().force_exit(side);
        }
    }

    /// A cell is "determined" if it has at least one traversal with
    /// multiple unique coordinates that has either traversed or forms
    /// a closed ring.
    pub fn determined(&self) -> bool {
        for t in &self.traversals {
            if !t.traversed() && !t.is_closed_ring() {
                continue;
            }
            if t.multiple_unique_coordinates() {
                return true;
            }
        }
        false
    }

    /// Total length of all traversals through this cell.
    pub fn traversal_length(&self) -> f64 {
        self.traversals
            .iter()
            .map(|t| measures::length(t.coords()))
            .sum()
    }

    /// Get the last traversal (read-only).
    pub fn last_traversal(&self) -> &Traversal {
        self.traversals.last().unwrap()
    }

    /// Get all traversals (for coverage fraction computation).
    pub fn traversals(&self) -> &[Traversal] {
        &self.traversals
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_cell() -> Cell {
        Cell::new(BBox::new(0.0, 0.0, 1.0, 1.0))
    }

    #[test]
    fn test_take_inside() {
        let mut cell = unit_cell();
        assert!(cell.take(&Coord::new(0.0, 0.5), None)); // boundary = inside
        assert!(cell.take(&Coord::new(0.5, 0.5), None)); // interior
        assert!(!cell.last_traversal().exited());
    }

    #[test]
    fn test_take_exits() {
        let mut cell = unit_cell();
        assert!(cell.take(&Coord::new(0.0, 0.5), None)); // enter from left
        assert!(cell.take(&Coord::new(0.5, 0.5), None)); // interior
        assert!(!cell.take(&Coord::new(2.0, 0.5), None)); // exits right
        assert!(cell.last_traversal().exited());
        assert_eq!(cell.last_traversal().exit_side(), Side::Right);
    }

    #[test]
    fn test_closed_ring_within_cell() {
        let mut cell = unit_cell();
        assert!(cell.take(&Coord::new(0.2, 0.2), None));
        assert!(cell.take(&Coord::new(0.8, 0.2), None));
        assert!(cell.take(&Coord::new(0.5, 0.8), None));
        assert!(cell.take(&Coord::new(0.2, 0.2), None)); // closes ring
        assert!(cell.last_traversal().is_closed_ring());
        // Closed rings are NOT "exited" — they have exit = Side::None.
        // They're detected via is_closed_ring() instead.
    }

    #[test]
    fn test_determined() {
        let mut cell = unit_cell();
        assert!(!cell.determined());

        cell.take(&Coord::new(0.0, 0.5), None);
        cell.take(&Coord::new(0.5, 0.5), None);
        cell.take(&Coord::new(2.0, 0.5), None);
        assert!(cell.determined());
    }

    #[test]
    fn test_multiple_traversals() {
        let mut cell = unit_cell();
        // First traversal: enter left, exit right
        cell.take(&Coord::new(0.0, 0.3), None);
        cell.take(&Coord::new(2.0, 0.3), None);
        assert!(cell.last_traversal().exited());

        // Second traversal: enter bottom, exit top
        cell.take(&Coord::new(0.5, 0.0), None);
        cell.take(&Coord::new(0.5, 2.0), None);
        assert_eq!(cell.traversals().len(), 2);
    }
}
