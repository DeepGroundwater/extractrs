pub mod download;
pub mod output;
pub mod proj;
pub mod raster;
pub mod vector;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum IoError {
    #[error("GDAL error: {0}")]
    Gdal(#[from] gdal::errors::GdalError),

    #[error("NetCDF error: {0}")]
    Netcdf(#[from] netcdf::Error),

    #[error("Projection error: {0}")]
    Proj(#[from] ::proj::ProjError),

    #[error("Projection creation error: {0}")]
    ProjCreate(#[from] ::proj::ProjCreateError),

    #[error("HTTP error: {0}")]
    Http(#[from] reqwest::Error),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("CSV error: {0}")]
    Csv(#[from] csv::Error),

    #[error("Invalid grid: {0}")]
    InvalidGrid(String),

    #[error("{0}")]
    Other(String),
}
