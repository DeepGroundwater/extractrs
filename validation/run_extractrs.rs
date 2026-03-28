// cargo run --release --package extractrs-pipeline --example validate

use anyhow::Result;
use extractrs_io::proj::CrsTransformer;
use extractrs_core::prelude::BBox;
use extractrs_io::raster::read_daymet_var;
use extractrs_io::vector::read_basins_bbox;
use extractrs_pipeline::coverage_cache::CoverageCache;
use extractrs_pipeline::processor::process_day;
use std::path::Path;
use std::time::Instant;

fn main() -> Result<()> {
    let nc_path_env = std::env::var("NC_PATH")
        .unwrap_or_else(|_| "/mnt/ssd1/data/daymet/daymet_v4_daily_na_tmin_2024.nc".into());
    let nc_path = Path::new(&nc_path_env);
    let var_name = std::env::var("VAR_NAME").unwrap_or_else(|_| "tmin".into());
    let time_idx: usize = std::env::var("TIME_IDX")
        .unwrap_or_else(|_| "182".into())
        .parse()?;
    let n_iters: usize = std::env::var("N_ITERS")
        .unwrap_or_else(|_| "1".into())
        .parse()?;

    let shp_path = Path::new(
        "/mnt/ssd1/data/merit/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp",
    );
    let out_path = Path::new("/tmp/extractrs_validation/extractrs_results.csv");
    std::fs::create_dir_all("/tmp/extractrs_validation")?;

    // Read CONUS basins
    let bbox_conus = BBox::new(-125.0, 24.0, -66.0, 50.0);
    let t0 = Instant::now();
    let basins = read_basins_bbox(shp_path, &bbox_conus)?;
    let read_time = t0.elapsed().as_secs_f64();
    eprintln!("read {} basins in {:.1}s", basins.len(), read_time);

    let transformer = CrsTransformer::wgs84_to_daymet_lcc()?;

    // Read full continental raster
    let t_read = Instant::now();
    eprintln!("reading full raster from {} ...", nc_path.display());
    let raster = read_daymet_var(nc_path, &var_name, time_idx)?;
    let raster_read_time = t_read.elapsed().as_secs_f64();
    eprintln!(
        "raster: {}x{} cells, dx={:.0}m, read in {:.1}s",
        raster.grid.cols(),
        raster.grid.rows(),
        raster.grid.dx(),
        raster_read_time
    );

    // Build coverage cache for ALL basins
    let t_cache = Instant::now();
    let cache = CoverageCache::build(&basins, &raster.grid, &transformer)?;
    let cache_time = t_cache.elapsed().as_secs_f64();
    eprintln!(
        "cache built in {:.1}s ({} entries)",
        cache_time,
        cache.entries.len()
    );

    // Process (with optional benchmark loop)
    let t_proc = Instant::now();
    let mut results = Vec::new();
    for _ in 0..n_iters {
        results = process_day(&cache, &raster);
    }
    let total_proc = t_proc.elapsed().as_secs_f64();
    let proc_per_iter = total_proc / n_iters as f64;
    eprintln!(
        "processed {} iters in {:.3}s ({:.3}s/iter)",
        n_iters, total_proc, proc_per_iter
    );

    // Write CSV
    let mut wtr = csv::Writer::from_path(out_path)?;
    wtr.write_record(["COMID", "mean"])?;
    let mut valid = 0;
    for &(comid, mean_val) in &results {
        if mean_val.is_finite() {
            wtr.write_record(&[comid.to_string(), format!("{:.6}", mean_val)])?;
            valid += 1;
        }
    }
    wtr.flush()?;
    eprintln!(
        "wrote {}/{} results to {}",
        valid,
        results.len(),
        out_path.display()
    );

    // Timing summary
    eprintln!("\n=== Timing ===");
    eprintln!("  basin read:  {:.1}s", read_time);
    eprintln!("  raster read: {:.1}s", raster_read_time);
    eprintln!("  cache build: {:.1}s", cache_time);
    eprintln!("  process/day: {:.3}s", proc_per_iter);
    eprintln!("  total:       {:.1}s", t0.elapsed().as_secs_f64());

    // 40-year projection
    let days_40yr: f64 = 365.0 * 40.0;
    eprintln!(
        "\n=== 40-Year Projection ({} basins) ===",
        cache.entries.len()
    );
    eprintln!("  cache (one-time):   {:.1}s", cache_time);
    eprintln!("  process per day:    {:.3}s", proc_per_iter);
    eprintln!(
        "  process 14600 days: {:.1}s ({:.1} min)",
        proc_per_iter * days_40yr,
        proc_per_iter * days_40yr / 60.0
    );
    eprintln!("  raster read/day:    {:.1}s", raster_read_time);
    eprintln!(
        "  total w/ reads:     {:.0}s ({:.1} min)",
        cache_time + (proc_per_iter + raster_read_time) * days_40yr,
        (cache_time + (proc_per_iter + raster_read_time) * days_40yr) / 60.0
    );

    Ok(())
}
