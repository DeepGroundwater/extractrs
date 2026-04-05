/// A single cell's contribution to statistics.
#[derive(Debug, Clone, Copy)]
pub struct CellContribution {
    pub coverage: f32,
    pub value: f64,
    pub weight: f64,
    pub x: f64,
    pub y: f64,
    pub value_defined: bool,
    pub weight_defined: bool,
}

/// Result value from a stat computation.
#[derive(Debug, Clone)]
pub enum StatValue {
    Float(f64),
    Int(i64),
    OptFloat(Option<f64>),
    FloatArray(Vec<f64>),
}

/// Trait for all statistical accumulators.
///
/// Designed for incremental processing — cells are fed one at a time.
/// Accumulators must be Send + Sync for parallel processing.
pub trait StatAccumulator: Send + Sync {
    /// Feed one cell's contribution.
    fn process(&mut self, cell: &CellContribution);

    /// Extract the final result.
    fn result(&self) -> StatValue;

    /// Name of this statistic.
    fn name(&self) -> &str;
}

/// Coverage-weighted count: sum of coverage fractions.
#[derive(Debug, Default)]
pub struct CountAccum {
    sum_ci: f64,
}

impl StatAccumulator for CountAccum {
    fn process(&mut self, cell: &CellContribution) {
        if cell.value_defined {
            self.sum_ci += cell.coverage as f64;
        }
    }

    fn result(&self) -> StatValue {
        StatValue::Float(self.sum_ci)
    }

    fn name(&self) -> &str {
        "count"
    }
}

/// Coverage-weighted sum of values.
#[derive(Debug, Default)]
pub struct SumAccum {
    sum_xici: f64,
}

impl StatAccumulator for SumAccum {
    fn process(&mut self, cell: &CellContribution) {
        if cell.value_defined {
            self.sum_xici += cell.value * cell.coverage as f64;
        }
    }

    fn result(&self) -> StatValue {
        StatValue::Float(self.sum_xici)
    }

    fn name(&self) -> &str {
        "sum"
    }
}

/// Coverage-weighted mean of values.
#[derive(Debug, Default)]
pub struct MeanAccum {
    sum_ci: f64,
    sum_xici: f64,
}

impl StatAccumulator for MeanAccum {
    fn process(&mut self, cell: &CellContribution) {
        if cell.value_defined {
            let c = cell.coverage as f64;
            self.sum_ci += c;
            self.sum_xici += cell.value * c;
        }
    }

    fn result(&self) -> StatValue {
        if self.sum_ci > 0.0 {
            StatValue::Float(self.sum_xici / self.sum_ci)
        } else {
            StatValue::OptFloat(None)
        }
    }

    fn name(&self) -> &str {
        "mean"
    }
}

/// Minimum value (ignoring coverage weights).
#[derive(Debug)]
#[derive(Default)]
pub struct MinAccum {
    min: Option<f64>,
}


impl StatAccumulator for MinAccum {
    fn process(&mut self, cell: &CellContribution) {
        if cell.value_defined && cell.coverage > 0.0 {
            self.min = Some(match self.min {
                Some(v) => v.min(cell.value),
                None => cell.value,
            });
        }
    }

    fn result(&self) -> StatValue {
        StatValue::OptFloat(self.min)
    }

    fn name(&self) -> &str {
        "min"
    }
}

/// Maximum value (ignoring coverage weights).
#[derive(Debug)]
#[derive(Default)]
pub struct MaxAccum {
    max: Option<f64>,
}


impl StatAccumulator for MaxAccum {
    fn process(&mut self, cell: &CellContribution) {
        if cell.value_defined && cell.coverage > 0.0 {
            self.max = Some(match self.max {
                Some(v) => v.max(cell.value),
                None => cell.value,
            });
        }
    }

    fn result(&self) -> StatValue {
        StatValue::OptFloat(self.max)
    }

    fn name(&self) -> &str {
        "max"
    }
}

/// Weighted sum: Σ(x_i * c_i * w_i).
///
/// Port of exactextract's weighted_sum (raster_stats.h: m_sum_xiciwi).
/// Both value and weight must be defined for a cell to contribute.
#[derive(Debug, Default)]
pub struct WeightedSumAccum {
    sum_xiciwi: f64,
}

impl StatAccumulator for WeightedSumAccum {
    fn process(&mut self, cell: &CellContribution) {
        if cell.value_defined && cell.weight_defined {
            self.sum_xiciwi += cell.value * cell.coverage as f64 * cell.weight;
        }
    }

    fn result(&self) -> StatValue {
        StatValue::Float(self.sum_xiciwi)
    }

    fn name(&self) -> &str {
        "weighted_sum"
    }
}

/// Weighted mean: Σ(x_i * c_i * w_i) / Σ(c_i * w_i).
///
/// 1:1 port of exactextract's weighted_mean (raster_stats.h).
/// Accumulates the same two terms as the C++:
///   - sum_ciwi   = Σ(c_i * w_i)       (denominator / weighted_count)
///   - sum_xiciwi = Σ(x_i * c_i * w_i) (numerator / weighted_sum)
#[derive(Debug, Default)]
pub struct WeightedMeanAccum {
    sum_ciwi: f64,
    sum_xiciwi: f64,
}

impl StatAccumulator for WeightedMeanAccum {
    fn process(&mut self, cell: &CellContribution) {
        if cell.value_defined && cell.weight_defined {
            let ciwi = cell.coverage as f64 * cell.weight;
            self.sum_ciwi += ciwi;
            self.sum_xiciwi += cell.value * ciwi;
        }
    }

    fn result(&self) -> StatValue {
        if self.sum_ciwi > 0.0 {
            StatValue::Float(self.sum_xiciwi / self.sum_ciwi)
        } else {
            StatValue::OptFloat(None)
        }
    }

    fn name(&self) -> &str {
        "weighted_mean"
    }
}
