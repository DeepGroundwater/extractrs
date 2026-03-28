use crate::accumulator::{CellContribution, StatAccumulator, StatValue};

/// Numerically stable online variance computation using West's (1979)
/// algorithm. Processes weighted observations one at a time.
///
/// Reference: D.H.D. West, "Updating Mean and Variance Estimates:
/// An Improved Method", Communications of the ACM, 22(9), 1979.
#[derive(Debug, Default)]
pub struct WestVariance {
    sum_w: f64,
    mean: f64,
    m2: f64,
}

impl WestVariance {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an observation with the given weight.
    pub fn process(&mut self, value: f64, weight: f64) {
        if weight <= 0.0 {
            return;
        }
        let new_sum_w = self.sum_w + weight;
        let delta = value - self.mean;
        let r = delta * weight / new_sum_w;
        self.m2 += self.sum_w * delta * r;
        self.mean += r;
        self.sum_w = new_sum_w;
    }

    /// Population variance.
    pub fn variance(&self) -> Option<f64> {
        if self.sum_w > 0.0 {
            Some(self.m2 / self.sum_w)
        } else {
            None
        }
    }

    /// Population standard deviation.
    pub fn stdev(&self) -> Option<f64> {
        self.variance().map(|v| v.sqrt())
    }

    /// Weighted mean.
    pub fn mean(&self) -> Option<f64> {
        if self.sum_w > 0.0 {
            Some(self.mean)
        } else {
            None
        }
    }
}

/// Coverage-weighted population variance.
#[derive(Debug, Default)]
pub struct VarianceAccum {
    west: WestVariance,
}

impl StatAccumulator for VarianceAccum {
    fn process(&mut self, cell: &CellContribution) {
        if cell.value_defined && cell.coverage > 0.0 {
            self.west.process(cell.value, cell.coverage as f64);
        }
    }

    fn result(&self) -> StatValue {
        StatValue::OptFloat(self.west.variance())
    }

    fn name(&self) -> &str {
        "variance"
    }
}

/// Coverage-weighted population standard deviation.
#[derive(Debug, Default)]
pub struct StdevAccum {
    west: WestVariance,
}

impl StatAccumulator for StdevAccum {
    fn process(&mut self, cell: &CellContribution) {
        if cell.value_defined && cell.coverage > 0.0 {
            self.west.process(cell.value, cell.coverage as f64);
        }
    }

    fn result(&self) -> StatValue {
        StatValue::OptFloat(self.west.stdev())
    }

    fn name(&self) -> &str {
        "stdev"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_west_variance_uniform() {
        let mut wv = WestVariance::new();
        // 5 equal-weight observations of value 10
        for _ in 0..5 {
            wv.process(10.0, 1.0);
        }
        assert!((wv.mean().unwrap() - 10.0).abs() < 1e-10);
        assert!((wv.variance().unwrap()).abs() < 1e-10);
    }

    #[test]
    fn test_west_variance_known() {
        let mut wv = WestVariance::new();
        // Values: 2, 4, 4, 4, 5, 5, 7, 9
        // Mean = 5, Variance = 4.0
        for v in [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0] {
            wv.process(v, 1.0);
        }
        assert!((wv.mean().unwrap() - 5.0).abs() < 1e-10);
        assert!((wv.variance().unwrap() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_west_variance_weighted() {
        let mut wv = WestVariance::new();
        // Value 0 with weight 3, value 1 with weight 1
        // Weighted mean = (0*3 + 1*1) / 4 = 0.25
        // Weighted var = (3*(0-0.25)^2 + 1*(1-0.25)^2) / 4 = (0.1875 + 0.5625)/4 = 0.1875
        wv.process(0.0, 3.0);
        wv.process(1.0, 1.0);
        assert!((wv.mean().unwrap() - 0.25).abs() < 1e-10);
        assert!((wv.variance().unwrap() - 0.1875).abs() < 1e-10);
    }
}
