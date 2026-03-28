use anyhow::Result;
use chrono::NaiveDate;
use clap::Parser;
use extractrs_core::prelude::BBox;
use extractrs_io::download::NcssClient;
use extractrs_io::output::CsvOutput;
use extractrs_io::proj::CrsTransformer;
use extractrs_io::raster::read_daymet_var;
use extractrs_io::vector::read_basins_bbox;
use extractrs_pipeline::coverage_cache::CoverageCache;
use extractrs_pipeline::processor::process_day;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(name = "extractrs-pipeline")]
#[command(about = "Extract zonal statistics from Daymet for MERIT basins")]
struct Args {
    /// Path to MERIT basins shapefile
    #[arg(long)]
    shapefile: PathBuf,

    /// Output CSV path
    #[arg(long)]
    output: PathBuf,

    /// Daymet variable (prcp, tmin, tmax, srad, vp, swe, dayl)
    #[arg(long, default_value = "prcp")]
    var: String,

    /// Start date (YYYY-MM-DD)
    #[arg(long)]
    start_date: String,

    /// End date (YYYY-MM-DD)
    #[arg(long)]
    end_date: String,

    /// Directory for temporary downloads
    #[arg(long, default_value = "/tmp/extractrs_downloads")]
    download_dir: PathBuf,

    /// Bounding box filter: west,south,east,north
    #[arg(long)]
    bbox: Option<String>,

    /// Max basins to process (for testing)
    #[arg(long)]
    max_basins: Option<usize>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let t0 = Instant::now();

    let start = NaiveDate::parse_from_str(&args.start_date, "%Y-%m-%d")?;
    let end = NaiveDate::parse_from_str(&args.end_date, "%Y-%m-%d")?;

    let bbox = match &args.bbox {
        Some(s) => {
            let parts: Vec<f64> = s.split(',').map(|p| p.trim().parse().unwrap()).collect();
            BBox::new(parts[0], parts[1], parts[2], parts[3])
        }
        None => {
            eprintln!("no --bbox specified, using default test region (N Georgia)");
            BBox::new(-85.0, 34.0, -84.0, 35.0)
        }
    };

    // Step 1: Read basins
    eprintln!("reading basins from {} ...", args.shapefile.display());
    let t_read = Instant::now();
    let mut basins = read_basins_bbox(&args.shapefile, &bbox)?;
    if let Some(max) = args.max_basins {
        basins.truncate(max);
    }
    eprintln!(
        "  loaded {} basins in {:.1}s",
        basins.len(),
        t_read.elapsed().as_secs_f64()
    );

    if basins.is_empty() {
        eprintln!("no basins found in bbox — exiting");
        return Ok(());
    }

    // Step 2: Download first day to establish grid
    let ncss = NcssClient::new(&args.download_dir);
    let first_nc = ncss.download_daymet(&args.var, start, &bbox)?;
    let first_raster = read_daymet_var(&first_nc, &args.var, 0)?;
    eprintln!(
        "  raster grid: {}x{} cells, dx={:.0}m dy={:.0}m",
        first_raster.grid.cols(),
        first_raster.grid.rows(),
        first_raster.grid.dx(),
        first_raster.grid.dy(),
    );

    // Step 3: Build coverage cache
    eprintln!("building coverage cache ...");
    let t_cache = Instant::now();
    let transformer = CrsTransformer::wgs84_to_daymet_lcc()?;
    let cache = CoverageCache::build(&basins, &first_raster.grid, &transformer)?;
    let cache_time = t_cache.elapsed().as_secs_f64();
    eprintln!("  cache built in {:.2}s", cache_time);

    // Step 4: Process first day (already downloaded)
    let mut csv = CsvOutput::new(&args.output)?;
    csv.write_header()?;

    let t_proc = Instant::now();
    let results = process_day(&cache, &first_raster);
    let first_proc_time = t_proc.elapsed().as_secs_f64();
    let date_str = start.format("%Y-%m-%d").to_string();
    for &(comid, mean_val) in &results {
        if mean_val.is_finite() {
            csv.write_row(comid, &date_str, mean_val)?;
        }
    }
    eprintln!(
        "  day {} processed in {:.4}s ({} basins)",
        date_str,
        first_proc_time,
        results.len()
    );

    // Clean up first file
    std::fs::remove_file(&first_nc).ok();

    // Step 5: Process remaining days
    let mut day = start.succ_opt().unwrap();
    let mut day_count = 1u32;
    let mut total_proc_time = first_proc_time;

    while day <= end {
        day_count += 1;
        let nc_path = ncss.download_daymet(&args.var, day, &bbox)?;
        let raster = read_daymet_var(&nc_path, &args.var, 0)?;

        let t_day = Instant::now();
        let results = process_day(&cache, &raster);
        let day_proc_time = t_day.elapsed().as_secs_f64();
        total_proc_time += day_proc_time;

        let date_str = day.format("%Y-%m-%d").to_string();
        for &(comid, mean_val) in &results {
            if mean_val.is_finite() {
                csv.write_row(comid, &date_str, mean_val)?;
            }
        }

        // Clean up
        std::fs::remove_file(&nc_path).ok();

        if day_count % 10 == 0 {
            eprintln!(
                "  processed {}/{} days, avg proc {:.4}s/day",
                day_count,
                (end - start).num_days() + 1,
                total_proc_time / day_count as f64
            );
        }

        day = day.succ_opt().unwrap();
    }

    csv.flush()?;

    let total_time = t0.elapsed().as_secs_f64();
    eprintln!("\n=== Summary ===");
    eprintln!("  basins: {}", cache.entries.len());
    eprintln!("  days: {}", day_count);
    eprintln!("  cache build: {:.2}s", cache_time);
    eprintln!("  total compute: {:.2}s ({:.4}s/day avg)", total_proc_time, total_proc_time / day_count as f64);
    eprintln!("  total elapsed: {:.1}s", total_time);

    // 40-year projection
    let days_40yr = 365 * 40;
    let compute_40yr = total_proc_time / day_count as f64 * days_40yr as f64;
    eprintln!("\n=== 40-Year Projection ===");
    eprintln!("  compute only: {:.1}s ({:.1} min)", compute_40yr, compute_40yr / 60.0);
    eprintln!("  cache (one-time): {:.1}s", cache_time);
    eprintln!(
        "  with download (~30s/day): {:.1} hours",
        (30.0 * days_40yr as f64 + compute_40yr + cache_time) / 3600.0
    );

    Ok(())
}
