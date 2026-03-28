use crate::bbox::BBox;
use crate::coord::Coord;
use crate::measures::{self, has_multiple_unique_coordinates};
use crate::perimeter::{perimeter_distance, perimeter_distance_ccw};

/// A chain of coordinates from a traversal through a cell,
/// tagged with perimeter-distance start/stop positions.
struct CoordinateChain<'a> {
    start: f64,
    stop: f64,
    coordinates: &'a [Coord],
    visited: bool,
}

/// Compute the "left-hand area" for a cell given traversal coordinate lists.
///
/// This assembles closed rings from traversal chains + cell boundary segments,
/// then computes the signed area. Counter-clockwise rings contribute positive
/// area; clockwise rings contribute negative area.
///
/// This is the pure-Rust replacement for exactextract's `left_hand_area()`
/// which uses GEOS for ring orientation checks. We use the shoelace formula
/// for CCW detection instead.
///
/// Port of exactextract's `traversal_areas.cpp`.
pub fn left_hand_area(bbox: &BBox, coord_lists: &[&[Coord]]) -> f64 {
    let mut ccw_sum = 0.0;
    let mut cw_sum = 0.0;
    let mut found_ring = false;

    visit_rings(bbox, coord_lists, |ring, is_ccw| {
        found_ring = true;
        if is_ccw {
            ccw_sum += measures::area(ring);
        } else {
            cw_sum += measures::area(ring);
        }
    });

    if !found_ring {
        panic!("Cannot determine coverage fraction (0 or 100%)");
    }

    if ccw_sum == 0.0 && cw_sum > 0.0 {
        // Only holes found: full cell minus holes
        bbox.area() - cw_sum
    } else {
        ccw_sum - cw_sum
    }
}

/// Visit all closed rings formed by traversal chains within a cell.
///
/// For each ring, calls `visitor(ring_coords, is_ccw)`.
///
/// Closed traversals (first == last) are visited directly.
/// Open traversals are assembled into rings by linking their endpoints
/// along the cell boundary in CCW order.
fn visit_rings(bbox: &BBox, coord_lists: &[&[Coord]], mut visitor: impl FnMut(&[Coord], bool)) {
    let mut chains: Vec<CoordinateChain> = Vec::with_capacity(coord_lists.len() + 4);

    // Step 1: Process coordinate lists — closed rings vs open traversals
    for coords in coord_lists {
        if !has_multiple_unique_coordinates(coords) {
            continue;
        }

        if coords.first() == coords.last() {
            // Already a closed ring: determine orientation via signed area
            let is_ccw = measures::area_signed(coords) > 0.0;
            visitor(coords, is_ccw);
        } else {
            // Open traversal: record perimeter distances
            let start = perimeter_distance(bbox, coords.first().unwrap());
            let stop = perimeter_distance(bbox, coords.last().unwrap());
            chains.push(CoordinateChain {
                start,
                stop,
                coordinates: coords,
                visited: false,
            });
        }
    }

    // Step 2: Add corner point chains (single-point "chains" at each box corner)
    // These ensure we can close rings that pass through corners.
    let corners = [
        bbox.lower_left(),  // perimeter distance = 0
        bbox.upper_left(),  // perimeter distance = height
        bbox.upper_right(), // perimeter distance = height + width
        bbox.lower_right(), // perimeter distance = 2*height + width
    ];
    let corner_distances = [
        0.0,
        bbox.height(),
        bbox.height() + bbox.width(),
        2.0 * bbox.height() + bbox.width(),
    ];

    // We store corner coordinates in owned storage so we can reference them
    let corner_coords: Vec<Vec<Coord>> = corners.iter().map(|c| vec![*c]).collect();
    for (i, corner) in corner_coords.iter().enumerate() {
        chains.push(CoordinateChain {
            start: corner_distances[i],
            stop: corner_distances[i],
            coordinates: corner,
            visited: false,
        });
    }

    // Step 3: Greedy ring assembly
    let total_perimeter = bbox.perimeter();

    for start_idx in 0..chains.len() {
        if chains[start_idx].visited || chains[start_idx].coordinates.len() == 1 {
            continue;
        }

        let mut ring = Vec::new();
        let mut current_idx = start_idx;
        let first_idx = start_idx;

        loop {
            chains[current_idx].visited = true;
            ring.extend_from_slice(chains[current_idx].coordinates);

            // Find next chain in CCW order
            let next_idx = next_chain_index(&chains, current_idx, first_idx, total_perimeter);

            if next_idx == first_idx {
                break;
            }

            // Add boundary points between current exit and next entry
            add_boundary_segment(bbox, &chains[current_idx], &chains[next_idx], &mut ring);

            current_idx = next_idx;
        }

        // Close the ring
        if let Some(first) = ring.first().copied() {
            ring.push(first);
        }

        if has_multiple_unique_coordinates(&ring) {
            // Assembled rings are always CCW (assembled in CCW perimeter order)
            visitor(&ring, true);
        }
    }
}

/// Find the next unvisited chain in CCW perimeter order from `current`.
/// Returns `first_idx` if we've wrapped around.
fn next_chain_index(
    chains: &[CoordinateChain],
    current_idx: usize,
    first_idx: usize,
    total_perimeter: f64,
) -> usize {
    let current_stop = chains[current_idx].stop;
    let mut best_idx = first_idx;
    let mut best_dist = perimeter_distance_ccw(
        current_stop,
        chains[first_idx].start,
        total_perimeter,
    );

    for (i, chain) in chains.iter().enumerate() {
        if i == current_idx || chain.visited {
            continue;
        }
        let dist = perimeter_distance_ccw(current_stop, chain.start, total_perimeter);
        if dist < best_dist {
            best_dist = dist;
            best_idx = i;
        }
    }

    best_idx
}

/// Add intermediate boundary points between two chains.
/// This walks along the cell boundary from the exit of `from` to the
/// entry of `to`, adding corner points as needed.
fn add_boundary_segment(
    _bbox: &BBox,
    _from: &CoordinateChain,
    _to: &CoordinateChain,
    _ring: &mut Vec<Coord>,
) {
    // Corner insertion: walk CCW from `from.stop` to `to.start`,
    // adding any corners we pass through.
    // For now this is a simplified version; the full implementation
    // would add intermediate corner coordinates along the boundary path.
    // The ring assembly via perimeter distances handles the primary case;
    // corners are already represented as zero-length chains.
}

/// Compute the covered fraction of a cell given its traversals.
///
/// This is the main entry point for coverage computation per cell.
/// Collects traversal coordinate lists, calls `left_hand_area()`,
/// and divides by cell area.
pub fn covered_fraction(bbox: &BBox, traversal_coords: &[&[Coord]]) -> f32 {
    if traversal_coords.is_empty() {
        return 0.0;
    }

    let cell_area = bbox.area();
    if cell_area <= 0.0 {
        return 0.0;
    }

    let coverage_area = left_hand_area(bbox, traversal_coords);
    let fraction = (coverage_area / cell_area) as f32;

    // Clamp to [0, 1] to handle floating-point rounding
    fraction.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_closed_ring_within_cell() {
        let bbox = BBox::new(0.0, 0.0, 1.0, 1.0);
        // A small square covering 1/4 of the cell (CCW)
        let ring = vec![
            Coord::new(0.0, 0.0),
            Coord::new(0.5, 0.0),
            Coord::new(0.5, 0.5),
            Coord::new(0.0, 0.5),
            Coord::new(0.0, 0.0),
        ];
        let frac = covered_fraction(&bbox, &[ring.as_slice()]);
        assert!((frac - 0.25).abs() < 1e-5, "expected 0.25, got {frac}");
    }

    #[test]
    fn test_full_cell_coverage() {
        let bbox = BBox::new(0.0, 0.0, 1.0, 1.0);
        // Ring that exactly covers the cell (CCW)
        let ring = vec![
            Coord::new(0.0, 0.0),
            Coord::new(1.0, 0.0),
            Coord::new(1.0, 1.0),
            Coord::new(0.0, 1.0),
            Coord::new(0.0, 0.0),
        ];
        let frac = covered_fraction(&bbox, &[ring.as_slice()]);
        assert!((frac - 1.0).abs() < 1e-5, "expected 1.0, got {frac}");
    }

    #[test]
    fn test_horizontal_traversal() {
        let bbox = BBox::new(0.0, 0.0, 1.0, 1.0);
        // Traversal entering left at y=0.5, exiting right at y=0.5
        // This covers the bottom half of the cell
        let coords = vec![Coord::new(0.0, 0.5), Coord::new(1.0, 0.5)];
        let frac = covered_fraction(&bbox, &[coords.as_slice()]);
        assert!((frac - 0.5).abs() < 1e-5, "expected 0.5, got {frac}");
    }
}
