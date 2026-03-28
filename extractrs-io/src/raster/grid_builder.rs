use crate::IoError;
use extractrs_core::prelude::{BBox, Grid, Bounded};

/// Build a Grid<Bounded> from 1D coordinate arrays (cell centers).
///
/// Assumes regular spacing. Daymet x is ascending, y is descending (north to south).
pub fn build_grid_from_coords(x: &[f64], y: &[f64]) -> Result<Grid<Bounded>, IoError> {
    if x.len() < 2 || y.len() < 2 {
        return Err(IoError::InvalidGrid(format!(
            "need at least 2 coordinates, got x={} y={}",
            x.len(),
            y.len()
        )));
    }

    let dx = (x[1] - x[0]).abs();
    let dy = (y[0] - y[1]).abs(); // y[0] > y[1] for north-to-south

    // Verify uniform spacing
    for i in 1..x.len() {
        let step = (x[i] - x[i - 1]).abs();
        if (step - dx).abs() > 1.0 {
            return Err(IoError::InvalidGrid(format!(
                "non-uniform x spacing at index {}: {} vs expected {}",
                i, step, dx
            )));
        }
    }
    for i in 1..y.len() {
        let step = (y[i - 1] - y[i]).abs();
        if (step - dy).abs() > 1.0 {
            return Err(IoError::InvalidGrid(format!(
                "non-uniform y spacing at index {}: {} vs expected {}",
                i, step, dy
            )));
        }
    }

    // Grid extent from cell centers → cell edges
    let xmin = x[0] - dx / 2.0;
    let xmax = x[x.len() - 1] + dx / 2.0;

    // Determine y orientation
    let y_descending = y[0] > y[y.len() - 1];
    let (ymin, ymax) = if y_descending {
        (y[y.len() - 1] - dy / 2.0, y[0] + dy / 2.0)
    } else {
        (y[0] - dy / 2.0, y[y.len() - 1] + dy / 2.0)
    };

    let grid = Grid::<Bounded>::new(BBox::new(xmin, ymin, xmax, ymax), dx, dy);

    // Sanity checks
    if grid.cols() != x.len() {
        return Err(IoError::InvalidGrid(format!(
            "grid cols {} != x len {}",
            grid.cols(),
            x.len()
        )));
    }
    if grid.rows() != y.len() {
        return Err(IoError::InvalidGrid(format!(
            "grid rows {} != y len {}",
            grid.rows(),
            y.len()
        )));
    }

    Ok(grid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_grid() {
        let x = vec![500.0, 1500.0, 2500.0]; // 1000m spacing
        let y = vec![3500.0, 2500.0, 1500.0, 500.0]; // descending, 1000m spacing
        let grid = build_grid_from_coords(&x, &y).unwrap();
        assert_eq!(grid.cols(), 3);
        assert_eq!(grid.rows(), 4);
        assert!((grid.dx() - 1000.0).abs() < 1e-6);
        assert!((grid.dy() - 1000.0).abs() < 1e-6);
        assert!((grid.extent().xmin - 0.0).abs() < 1e-6);
        assert!((grid.extent().ymax - 4000.0).abs() < 1e-6);
    }
}
