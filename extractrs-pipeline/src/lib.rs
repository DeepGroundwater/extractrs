pub mod coverage_cache;
pub mod processor;

use chrono::NaiveDate;
use extractrs_core::prelude::BBox;
use std::path::PathBuf;

pub struct PipelineConfig {
    pub shapefile_path: PathBuf,
    pub download_dir: PathBuf,
    pub output_path: PathBuf,
    pub variable: String,
    pub start_date: NaiveDate,
    pub end_date: NaiveDate,
    /// Optional WGS84 bbox to filter basins spatially.
    pub bbox_filter: Option<BBox>,
    /// Max basins to process (for testing).
    pub max_basins: Option<usize>,
}
