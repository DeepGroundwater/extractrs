use crate::coord::Coord;
use crate::side::Side;

/// Records the path of a polygon ring through a single grid cell.
///
/// A traversal captures the entry side, coordinates within the cell,
/// and exit side. A cell may have multiple traversals if the ring
/// enters and exits the cell multiple times.
///
/// Port of exactextract's `Traversal` class.
#[derive(Debug)]
pub struct Traversal {
    coords: Vec<Coord>,
    entry: Side,
    exit: Side,
}

impl Default for Traversal {
    fn default() -> Self {
        Self::new()
    }
}

impl Traversal {
    pub fn new() -> Self {
        Self {
            coords: Vec::new(),
            entry: Side::None,
            exit: Side::None,
        }
    }

    /// Mark entry into the cell from a boundary side.
    pub fn enter(&mut self, c: Coord, s: Side) {
        debug_assert!(self.coords.is_empty(), "Traversal already started");
        self.coords.push(c);
        self.entry = s;
    }

    /// Add an intermediate coordinate within the cell.
    #[inline]
    pub fn add(&mut self, c: Coord) {
        self.coords.push(c);
    }

    /// Mark exit from the cell through a boundary side.
    pub fn exit(&mut self, c: Coord, s: Side) {
        self.coords.push(c);
        self.exit = s;
    }

    /// Force-set the exit side without adding a coordinate.
    /// Used for closed rings that don't actually exit the cell.
    #[inline]
    pub fn force_exit(&mut self, s: Side) {
        self.exit = s;
    }

    /// True if this traversal forms a closed ring (first == last, size >= 3).
    pub fn is_closed_ring(&self) -> bool {
        self.coords.len() >= 3 && self.coords[0] == *self.coords.last().unwrap()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.coords.is_empty()
    }

    #[inline]
    pub fn entered(&self) -> bool {
        self.entry != Side::None
    }

    #[inline]
    pub fn exited(&self) -> bool {
        self.exit != Side::None
    }

    /// True if this traversal has both entered and exited the cell.
    #[inline]
    pub fn traversed(&self) -> bool {
        self.entered() && self.exited()
    }

    /// True if the coordinates contain at least two distinct points.
    pub fn multiple_unique_coordinates(&self) -> bool {
        for i in 1..self.coords.len() {
            if self.coords[0] != self.coords[i] {
                return true;
            }
        }
        false
    }

    #[inline]
    pub fn entry_side(&self) -> Side {
        self.entry
    }

    #[inline]
    pub fn exit_side(&self) -> Side {
        self.exit
    }

    #[inline]
    pub fn coords(&self) -> &[Coord] {
        &self.coords
    }

    #[inline]
    pub fn last_coordinate(&self) -> &Coord {
        self.coords.last().expect("empty traversal has no last coordinate")
    }

    pub fn exit_coordinate(&self) -> &Coord {
        debug_assert!(self.exited(), "no exit coordinate for incomplete traversal");
        self.last_coordinate()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_enter_add_exit() {
        let mut t = Traversal::new();
        assert!(t.is_empty());
        assert!(!t.entered());
        assert!(!t.exited());

        t.enter(Coord::new(0.0, 0.5), Side::Left);
        assert!(t.entered());
        assert!(!t.exited());

        t.add(Coord::new(0.5, 0.5));

        t.exit(Coord::new(1.0, 0.5), Side::Right);
        assert!(t.traversed());
        assert_eq!(t.exit_side(), Side::Right);
        assert_eq!(*t.exit_coordinate(), Coord::new(1.0, 0.5));
    }

    #[test]
    fn test_closed_ring() {
        let mut t = Traversal::new();
        t.enter(Coord::new(0.2, 0.2), Side::None);
        t.add(Coord::new(0.8, 0.2));
        t.add(Coord::new(0.5, 0.8));
        t.add(Coord::new(0.2, 0.2)); // close the ring
        assert!(t.is_closed_ring());
    }

    #[test]
    fn test_multiple_unique() {
        let mut t = Traversal::new();
        t.enter(Coord::new(0.5, 0.5), Side::None);
        t.add(Coord::new(0.5, 0.5));
        assert!(!t.multiple_unique_coordinates());

        t.add(Coord::new(0.6, 0.5));
        assert!(t.multiple_unique_coordinates());
    }
}
