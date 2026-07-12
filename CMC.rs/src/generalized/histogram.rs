//! Histogram and logarithmic density-of-states storage.

use crate::generalized::GeneralizedError;

/// Fixed-size visit histogram with checked total-count accounting.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Histogram {
    counts: Vec<u64>,
    total: u64,
}

impl Histogram {
    pub fn new(bins: usize) -> Result<Self, GeneralizedError> {
        if bins == 0 {
            return Err(GeneralizedError::new(
                "histogram must contain at least one bin",
            ));
        }
        Ok(Self {
            counts: vec![0; bins],
            total: 0,
        })
    }

    #[inline]
    pub fn bins(&self) -> usize {
        self.counts.len()
    }

    #[inline]
    pub fn counts(&self) -> &[u64] {
        &self.counts
    }

    #[inline]
    pub const fn total(&self) -> u64 {
        self.total
    }

    #[inline]
    pub fn count(&self, bin: usize) -> u64 {
        self.counts[bin]
    }

    #[inline]
    pub fn visited_bins(&self) -> usize {
        self.counts.iter().filter(|&&count| count > 0).count()
    }

    pub fn record(&mut self, bin: usize) {
        let count = self
            .counts
            .get_mut(bin)
            .expect("histogram bin out of range");
        *count = count.checked_add(1).expect("histogram bin count overflow");
        self.total = self.total.checked_add(1).expect("histogram total overflow");
    }

    pub fn clear(&mut self) {
        self.counts.fill(0);
        self.total = 0;
    }

    /// Flatness over visited bins, with an explicit minimum represented fraction.
    pub fn is_flat(&self, minimum_fraction_of_mean: f64, minimum_visited_fraction: f64) -> bool {
        assert!(
            minimum_fraction_of_mean.is_finite() && (0.0..=1.0).contains(&minimum_fraction_of_mean)
        );
        assert!(
            minimum_visited_fraction.is_finite() && (0.0..=1.0).contains(&minimum_visited_fraction)
        );
        let visited = self.visited_bins();
        if visited == 0 {
            return false;
        }
        let required = (minimum_visited_fraction * self.bins() as f64).ceil() as usize;
        if visited < required.max(1) {
            return false;
        }
        let mean = self.total as f64 / visited as f64;
        self.counts
            .iter()
            .copied()
            .filter(|&count| count > 0)
            .all(|count| count as f64 >= minimum_fraction_of_mean * mean)
    }

    pub(crate) fn from_counts(counts: Vec<u64>) -> Result<Self, GeneralizedError> {
        if counts.is_empty() {
            return Err(GeneralizedError::new("histogram snapshot contains no bins"));
        }
        let total = counts.iter().try_fold(0_u64, |sum, &count| {
            sum.checked_add(count)
                .ok_or_else(|| GeneralizedError::new("histogram snapshot total overflow"))
        })?;
        Ok(Self { counts, total })
    }
}

/// Logarithmic density of states, defined up to an additive constant.
#[derive(Debug, Clone, PartialEq)]
pub struct LogDensityOfStates {
    values: Vec<f64>,
    visited: Vec<bool>,
}

impl LogDensityOfStates {
    pub fn new(bins: usize) -> Result<Self, GeneralizedError> {
        if bins == 0 {
            return Err(GeneralizedError::new(
                "density of states must contain at least one bin",
            ));
        }
        Ok(Self {
            values: vec![0.0; bins],
            visited: vec![false; bins],
        })
    }

    pub fn from_values(values: Vec<f64>, visited: Vec<bool>) -> Result<Self, GeneralizedError> {
        if values.is_empty() || values.len() != visited.len() {
            return Err(GeneralizedError::new(
                "density-of-states value and visited buffers must have equal non-zero length",
            ));
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(GeneralizedError::new(
                "density of states contains a non-finite value",
            ));
        }
        Ok(Self { values, visited })
    }

    #[inline]
    pub fn bins(&self) -> usize {
        self.values.len()
    }

    #[inline]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    #[inline]
    pub fn visited(&self) -> &[bool] {
        &self.visited
    }

    #[inline]
    pub fn value(&self, bin: usize) -> f64 {
        self.values[bin]
    }

    #[inline]
    pub fn is_visited(&self, bin: usize) -> bool {
        self.visited[bin]
    }

    #[inline]
    pub fn visited_bins(&self) -> usize {
        self.visited.iter().filter(|&&visited| visited).count()
    }

    pub fn add_visit_weight(&mut self, bin: usize, increment: f64) {
        assert!(
            increment.is_finite() && increment >= 0.0,
            "log-DOS increment must be finite and non-negative"
        );
        self.values[bin] += increment;
        assert!(self.values[bin].is_finite(), "log-DOS value overflow");
        self.visited[bin] = true;
    }

    /// Return values shifted so the largest visited value is zero.
    ///
    /// Unvisited bins remain `None`, preserving the distinction between an
    /// estimated zero and a macrostate never reached by the walk.
    pub fn shifted_to_max_zero(&self) -> Vec<Option<f64>> {
        let maximum = self
            .values
            .iter()
            .zip(&self.visited)
            .filter_map(|(&value, &visited)| visited.then_some(value))
            .max_by(f64::total_cmp);
        let Some(maximum) = maximum else {
            return vec![None; self.values.len()];
        };
        self.values
            .iter()
            .zip(&self.visited)
            .map(|(&value, &visited)| visited.then_some(value - maximum))
            .collect()
    }

    /// Recenter all bins so the largest visited estimate is zero.
    ///
    /// This is an additive gauge transformation and therefore leaves every
    /// transition probability and reweighted observable unchanged.
    pub fn normalize_max_zero(&mut self) {
        if let Some(maximum) = self
            .values
            .iter()
            .zip(&self.visited)
            .filter_map(|(&value, &visited)| visited.then_some(value))
            .max_by(f64::total_cmp)
        {
            self.shift(-maximum);
        }
    }

    /// Shift all values by a constant without changing physical predictions.
    pub fn shift(&mut self, offset: f64) {
        assert!(offset.is_finite(), "log-DOS shift must be finite");
        for value in &mut self.values {
            *value += offset;
            assert!(value.is_finite(), "log-DOS shift overflow");
        }
    }
}
