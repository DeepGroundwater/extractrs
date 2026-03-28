use extractrs_core::prelude::{Bounded, CoverageResult, Grid, raster_cell_intersection};
use extractrs_io::proj::CrsTransformer;
use extractrs_io::vector::Basin;
use ndarray::Array2;

/// Pre-computed coverage fractions for one basin against a fixed raster grid.
pub struct CoverageEntry {
    pub comid: i64,
    pub row_offset: usize,
    pub col_offset: usize,
    pub fractions: Array2<f32>,
}

/// Cache of coverage fractions for all basins against a single raster grid.
pub struct CoverageCache {
    pub entries: Vec<CoverageEntry>,
    pub grid: Grid<Bounded>,
}

impl CoverageCache {
    /// Build cache by computing coverage fractions for each basin.
    ///
    /// Basins are in WGS84; the raster grid is in LCC. The transformer
    /// converts basin coords to LCC before intersection.
    pub fn build(
        basins: &[Basin],
        raster_grid: &Grid<Bounded>,
        transformer: &CrsTransformer,
    ) -> Result<Self, anyhow::Error> {
        let mut entries = Vec::with_capacity(basins.len());

        for (i, basin) in basins.iter().enumerate() {
            if i % 1000 == 0 && i > 0 {
                eprintln!("  coverage cache: {}/{} basins", i, basins.len());
            }

            // Compute coverage for each polygon part
            let mut part_results: Vec<(CoverageResult, usize, usize)> = Vec::new();

            for part in &basin.parts {
                let lcc_exterior = transformer.transform_coords(&part.exterior)?;
                let lcc_holes: Vec<Vec<_>> = part
                    .holes
                    .iter()
                    .map(|h| transformer.transform_coords(h))
                    .collect::<Result<_, _>>()?;

                let result: CoverageResult =
                    raster_cell_intersection(raster_grid, &lcc_exterior, &lcc_holes);

                if result.fractions.dim() == (0, 0) {
                    continue;
                }

                let ro = raster_grid
                    .get_row(result.grid.extent().ymax - result.grid.dy() * 0.5);
                let co = raster_grid
                    .get_column(result.grid.extent().xmin + result.grid.dx() * 0.5);

                part_results.push((result, ro, co));
            }

            if part_results.is_empty() {
                continue;
            }

            // Fast path: single part (common case)
            if part_results.len() == 1 {
                let (result, row_offset, col_offset) = part_results.into_iter().next().unwrap();
                entries.push(CoverageEntry {
                    comid: basin.comid,
                    row_offset,
                    col_offset,
                    fractions: result.fractions,
                });
                continue;
            }

            // Multiple parts: accumulate fractions into a single array
            let mut min_row = usize::MAX;
            let mut min_col = usize::MAX;
            let mut max_row = 0usize;
            let mut max_col = 0usize;

            for (result, ro, co) in &part_results {
                let (rows, cols) = result.fractions.dim();
                min_row = min_row.min(*ro);
                min_col = min_col.min(*co);
                max_row = max_row.max(*ro + rows);
                max_col = max_col.max(*co + cols);
            }

            let total_rows = max_row - min_row;
            let total_cols = max_col - min_col;
            let mut combined = Array2::<f32>::zeros((total_rows, total_cols));

            for (result, ro, co) in &part_results {
                let (rows, cols) = result.fractions.dim();
                let dr = ro - min_row;
                let dc = co - min_col;
                for r in 0..rows {
                    for c in 0..cols {
                        combined[[dr + r, dc + c]] += result.fractions[[r, c]];
                    }
                }
            }

            entries.push(CoverageEntry {
                comid: basin.comid,
                row_offset: min_row,
                col_offset: min_col,
                fractions: combined,
            });
        }

        eprintln!(
            "  coverage cache: {}/{} basins intersect raster",
            entries.len(),
            basins.len()
        );

        Ok(Self {
            entries,
            grid: raster_grid.clone(),
        })
    }
}
