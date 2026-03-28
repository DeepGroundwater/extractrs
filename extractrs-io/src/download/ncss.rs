use crate::IoError;
use extractrs_core::prelude::BBox;
use chrono::NaiveDate;
use std::path::{Path, PathBuf};

const DAYMET_NCSS_BASE: &str =
    "https://thredds.daac.ornl.gov/thredds/ncss/grid/ornldaac/1840";

pub struct NcssClient {
    client: reqwest::blocking::Client,
    download_dir: PathBuf,
}

impl NcssClient {
    pub fn new(download_dir: &Path) -> Self {
        std::fs::create_dir_all(download_dir).ok();
        Self {
            client: reqwest::blocking::Client::builder()
                .timeout(std::time::Duration::from_secs(300))
                .build()
                .expect("failed to build HTTP client"),
            download_dir: download_dir.to_path_buf(),
        }
    }

    /// Download a single day of a Daymet variable for a geographic bounding box.
    ///
    /// Returns the path to the downloaded NetCDF file.
    pub fn download_daymet(
        &self,
        var: &str,
        date: NaiveDate,
        bbox_wgs84: &BBox,
    ) -> Result<PathBuf, IoError> {
        let year = date.format("%Y");
        let date_str = date.format("%Y-%m-%d");

        let dataset = format!("daymet_v4_daily_na_{}_{}.nc", var, year);
        let url = format!(
            "{base}/{dataset}?var={var}\
             &north={north}&south={south}&east={east}&west={west}\
             &disableProjSubset=on&horizStride=1\
             &time_start={date}T12:00:00Z&time_end={date}T12:00:00Z\
             &timeStride=1&accept=netcdf",
            base = DAYMET_NCSS_BASE,
            dataset = dataset,
            var = var,
            north = bbox_wgs84.ymax,
            south = bbox_wgs84.ymin,
            east = bbox_wgs84.xmax,
            west = bbox_wgs84.xmin,
            date = date_str,
        );

        let out_path = self
            .download_dir
            .join(format!("daymet_{}_{}_{}.nc", var, date_str, "subset"));

        eprintln!("downloading {} {} ...", var, date_str);
        let resp = self.client.get(&url).send()?;

        if !resp.status().is_success() {
            return Err(IoError::Other(format!(
                "NCSS request failed: {} — {}",
                resp.status(),
                url
            )));
        }

        let bytes = resp.bytes()?;
        std::fs::write(&out_path, &bytes)?;
        eprintln!("  saved {} ({} bytes)", out_path.display(), bytes.len());

        Ok(out_path)
    }
}
