use crate::IoError;
use extractrs_core::prelude::{BBox, Coord};
use gdal::vector::LayerAccess;
use gdal::Dataset;
use std::path::Path;

/// A single polygon part (exterior ring + holes).
#[derive(Debug, Clone)]
pub struct PolygonPart {
    pub exterior: Vec<Coord>,
    pub holes: Vec<Vec<Coord>>,
}

/// A basin polygon with its COMID identifier.
/// MultiPolygon geometries store all sub-polygons in `parts`.
#[derive(Debug, Clone)]
pub struct Basin {
    pub comid: i64,
    pub parts: Vec<PolygonPart>,
    pub bbox: BBox,
}

/// Read all basins from a shapefile.
pub fn read_basins(shp_path: &Path) -> Result<Vec<Basin>, IoError> {
    let ds = Dataset::open(shp_path)?;
    let mut layer = ds.layer(0)?;

    // Find COMID field index
    let comid_idx = layer
        .defn()
        .fields()
        .position(|f| f.name() == "COMID")
        .ok_or_else(|| IoError::Other("no COMID field in shapefile".into()))?;

    let mut basins = Vec::with_capacity(layer.feature_count() as usize);
    for feature in layer.features() {
        if let Some(basin) = feature_to_basin(&feature, comid_idx)? {
            basins.push(basin);
        }
    }

    Ok(basins)
}

/// Read basins within a WGS84 bounding box.
pub fn read_basins_bbox(shp_path: &Path, bbox: &BBox) -> Result<Vec<Basin>, IoError> {
    let ds = Dataset::open(shp_path)?;
    let mut layer = ds.layer(0)?;

    let comid_idx = layer
        .defn()
        .fields()
        .position(|f| f.name() == "COMID")
        .ok_or_else(|| IoError::Other("no COMID field in shapefile".into()))?;

    layer.set_spatial_filter_rect(bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax);

    let mut basins = Vec::new();
    for feature in layer.features() {
        if let Some(basin) = feature_to_basin(&feature, comid_idx)? {
            basins.push(basin);
        }
    }

    Ok(basins)
}

fn feature_to_basin(
    feature: &gdal::vector::Feature,
    comid_idx: usize,
) -> Result<Option<Basin>, IoError> {
    let comid = match feature.field_as_integer64(comid_idx)? {
        Some(v) => v,
        None => return Ok(None),
    };

    let geom = match feature.geometry() {
        Some(g) => g,
        None => return Ok(None),
    };

    // Extract the polygon geometry (handle MultiPolygon by taking first)
    let geom_type = geom.geometry_type();
    let is_multi = matches!(
        geom_type,
        gdal::vector::OGRwkbGeometryType::wkbMultiPolygon
            | gdal::vector::OGRwkbGeometryType::wkbMultiPolygon25D
    );

    if is_multi {
        let n_parts = geom.geometry_count();
        if n_parts == 0 {
            return Ok(None);
        }

        let mut parts = Vec::with_capacity(n_parts);
        let mut bbox = BBox::make_empty();

        for i in 0..n_parts {
            let sub = unsafe { geom.get_unowned_geometry(i) };
            let ring_count = sub.geometry_count();
            if ring_count == 0 {
                continue;
            }

            let ext_ring = unsafe { sub.get_unowned_geometry(0) };
            let exterior = ring_to_coords(&ext_ring);
            if exterior.len() < 4 {
                continue;
            }

            let mut holes = Vec::new();
            for j in 1..ring_count {
                let int_ring = unsafe { sub.get_unowned_geometry(j) };
                let coords = ring_to_coords(&int_ring);
                if coords.len() >= 4 {
                    holes.push(coords);
                }
            }

            bbox = bbox.expand_to_include(&compute_bbox(&exterior));
            parts.push(PolygonPart { exterior, holes });
        }

        if parts.is_empty() {
            return Ok(None);
        }

        return Ok(Some(Basin { comid, parts, bbox }));
    }

    // Simple Polygon
    let ring_count = geom.geometry_count();
    if ring_count == 0 {
        return Ok(None);
    }

    let ext_ring = unsafe { geom.get_unowned_geometry(0) };
    let exterior = ring_to_coords(&ext_ring);
    if exterior.len() < 4 {
        return Ok(None);
    }

    let mut holes = Vec::new();
    for i in 1..ring_count {
        let int_ring = unsafe { geom.get_unowned_geometry(i) };
        let coords = ring_to_coords(&int_ring);
        if coords.len() >= 4 {
            holes.push(coords);
        }
    }

    let bbox = compute_bbox(&exterior);
    Ok(Some(Basin {
        comid,
        parts: vec![PolygonPart { exterior, holes }],
        bbox,
    }))
}

fn ring_to_coords(geom: &gdal::vector::Geometry) -> Vec<Coord> {
    let n = geom.point_count();
    let mut coords = Vec::with_capacity(n);
    for i in 0..n {
        let (x, y, _) = geom.get_point(i as i32);
        coords.push(Coord::new(x, y));
    }
    coords
}

fn compute_bbox(coords: &[Coord]) -> BBox {
    let mut bbox = BBox::make_empty();
    for c in coords {
        if c.x < bbox.xmin { bbox.xmin = c.x; }
        if c.x > bbox.xmax { bbox.xmax = c.x; }
        if c.y < bbox.ymin { bbox.ymin = c.y; }
        if c.y > bbox.ymax { bbox.ymax = c.y; }
    }
    bbox
}
