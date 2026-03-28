use crate::IoError;
use std::path::Path;

pub struct CsvOutput {
    writer: csv::Writer<std::fs::File>,
}

impl CsvOutput {
    pub fn new(path: &Path) -> Result<Self, IoError> {
        let writer = csv::Writer::from_path(path)?;
        Ok(Self { writer })
    }

    pub fn write_header(&mut self) -> Result<(), IoError> {
        self.writer
            .write_record(["comid", "date", "mean_value"])?;
        Ok(())
    }

    pub fn write_row(
        &mut self,
        comid: i64,
        date: &str,
        mean_value: f64,
    ) -> Result<(), IoError> {
        self.writer.write_record(&[
            comid.to_string(),
            date.to_string(),
            format!("{:.6}", mean_value),
        ])?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<(), IoError> {
        self.writer.flush()?;
        Ok(())
    }
}
