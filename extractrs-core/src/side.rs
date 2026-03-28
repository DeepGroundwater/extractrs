/// Which side of a cell boundary a traversal enters/exits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Side {
    None,
    Left,
    Right,
    Top,
    Bottom,
}
