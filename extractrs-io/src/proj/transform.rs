use crate::IoError;
use extractrs_core::prelude::{BBox, Coord};
use proj::Proj;

const DAYMET_LCC: &str = "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs";
const WGS84: &str = "+proj=longlat +datum=WGS84 +no_defs";

pub struct CrsTransformer {
    proj: Proj,
}

impl CrsTransformer {
    /// Create a transformer from WGS84 (lon/lat) to Daymet LCC (meters).
    pub fn wgs84_to_daymet_lcc() -> Result<Self, IoError> {
        let proj = Proj::new_known_crs(WGS84, DAYMET_LCC, None)?;
        Ok(Self { proj })
    }

    /// Create a transformer from Daymet LCC (meters) to WGS84 (lon/lat).
    pub fn daymet_lcc_to_wgs84() -> Result<Self, IoError> {
        let proj = Proj::new_known_crs(DAYMET_LCC, WGS84, None)?;
        Ok(Self { proj })
    }

    /// Transform a slice of extractrs Coords.
    pub fn transform_coords(&self, coords: &[Coord]) -> Result<Vec<Coord>, IoError> {
        coords
            .iter()
            .map(|c| {
                let (x, y) = self.proj.convert((c.x, c.y))?;
                Ok(Coord::new(x, y))
            })
            .collect()
    }

    /// Transform a bounding box. Samples corners and midpoints for accuracy
    /// under non-affine projections (LCC).
    pub fn transform_bbox(&self, bbox: &BBox) -> Result<BBox, IoError> {
        let sample_points = vec![
            Coord::new(bbox.xmin, bbox.ymin),
            Coord::new(bbox.xmax, bbox.ymin),
            Coord::new(bbox.xmax, bbox.ymax),
            Coord::new(bbox.xmin, bbox.ymax),
            Coord::new((bbox.xmin + bbox.xmax) / 2.0, bbox.ymin),
            Coord::new((bbox.xmin + bbox.xmax) / 2.0, bbox.ymax),
            Coord::new(bbox.xmin, (bbox.ymin + bbox.ymax) / 2.0),
            Coord::new(bbox.xmax, (bbox.ymin + bbox.ymax) / 2.0),
        ];
        let transformed = self.transform_coords(&sample_points)?;
        let mut out = BBox::make_empty();
        for c in &transformed {
            if c.x < out.xmin {
                out.xmin = c.x;
            }
            if c.x > out.xmax {
                out.xmax = c.x;
            }
            if c.y < out.ymin {
                out.ymin = c.y;
            }
            if c.y > out.ymax {
                out.ymax = c.y;
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip() {
        let fwd = CrsTransformer::wgs84_to_daymet_lcc().unwrap();
        let inv = CrsTransformer::daymet_lcc_to_wgs84().unwrap();
        // Denver, CO: roughly -104.99, 39.74
        let orig = Coord::new(-104.99, 39.74);
        let lcc = fwd.transform_coords(&[orig]).unwrap();
        let back = inv.transform_coords(&lcc).unwrap();
        assert!((back[0].x - orig.x).abs() < 1e-6, "lon: {} vs {}", back[0].x, orig.x);
        assert!((back[0].y - orig.y).abs() < 1e-6, "lat: {} vs {}", back[0].y, orig.y);
    }

    #[test]
    fn test_bbox_transform() {
        let fwd = CrsTransformer::wgs84_to_daymet_lcc().unwrap();
        let bbox = BBox::new(-85.0, 34.0, -84.0, 35.0);
        let lcc_bbox = fwd.transform_bbox(&bbox).unwrap();
        // LCC coordinates should be in meters, roughly 1.3M east, -0.9M north
        assert!(lcc_bbox.xmin > 1_000_000.0, "xmin={}", lcc_bbox.xmin);
        assert!(lcc_bbox.width() > 50_000.0, "width={}", lcc_bbox.width());
    }
}
