//! One-dimensional discrete and continuous macrostate axes.

use crate::generalized::GeneralizedError;

/// Maps a scalar macrostate value to a stable zero-based bin.
pub trait MacrostateAxis: Clone + Send + Sync {
    /// Number of bins.
    fn bins(&self) -> usize;

    /// Return the bin containing `value`, or `None` outside the represented range.
    fn bin(&self, value: f64) -> Option<usize>;

    /// Representative scalar value of one bin.
    fn center(&self, bin: usize) -> f64;

    /// Materialize all representative values in bin order.
    fn centers(&self) -> Vec<f64> {
        (0..self.bins()).map(|bin| self.center(bin)).collect()
    }
}

/// Uniform finite-width axis over the closed interval `[min, max]`.
///
/// Interior bins are half-open. The final bin includes `max`, avoiding a
/// special out-of-range result for an exactly represented upper boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BinnedAxis {
    min: f64,
    max: f64,
    width: f64,
    bins: usize,
}

impl BinnedAxis {
    /// Construct `bins` uniform bins spanning `[min, max]`.
    pub fn new(min: f64, max: f64, bins: usize) -> Result<Self, GeneralizedError> {
        if !min.is_finite() || !max.is_finite() || max <= min {
            return Err(GeneralizedError::new(
                "axis bounds must be finite with max greater than min",
            ));
        }
        if bins == 0 {
            return Err(GeneralizedError::new("axis must contain at least one bin"));
        }
        let width = (max - min) / bins as f64;
        if !width.is_finite() || width <= 0.0 {
            return Err(GeneralizedError::new(
                "axis width must be finite and positive",
            ));
        }
        Ok(Self {
            min,
            max,
            width,
            bins,
        })
    }

    /// Construct a uniform axis from a requested bin width.
    ///
    /// The span must contain an integral number of bins within floating-point
    /// tolerance; this prevents silently changing the requested upper bound.
    pub fn from_width(min: f64, max: f64, width: f64) -> Result<Self, GeneralizedError> {
        if !width.is_finite() || width <= 0.0 {
            return Err(GeneralizedError::new(
                "axis width must be finite and positive",
            ));
        }
        if !min.is_finite() || !max.is_finite() || max <= min {
            return Err(GeneralizedError::new(
                "axis bounds must be finite with max greater than min",
            ));
        }
        let raw_bins = (max - min) / width;
        let rounded = raw_bins.round();
        let tolerance = 64.0 * f64::EPSILON * raw_bins.abs().max(1.0);
        if (raw_bins - rounded).abs() > tolerance || rounded < 1.0 {
            return Err(GeneralizedError::new(
                "axis span must be an integral multiple of the requested width",
            ));
        }
        Self::new(min, max, rounded as usize)
    }

    #[inline]
    pub const fn min(self) -> f64 {
        self.min
    }

    #[inline]
    pub const fn max(self) -> f64 {
        self.max
    }

    #[inline]
    pub const fn width(self) -> f64 {
        self.width
    }
}

impl MacrostateAxis for BinnedAxis {
    #[inline]
    fn bins(&self) -> usize {
        self.bins
    }

    fn bin(&self, value: f64) -> Option<usize> {
        if !value.is_finite() || value < self.min || value > self.max {
            return None;
        }
        if value == self.max {
            return Some(self.bins - 1);
        }
        let index = ((value - self.min) / self.width).floor() as usize;
        // Rounding can place the greatest representable value below `max`
        // exactly on `bins`; the explicit bounds check above makes clamping safe.
        Some(index.min(self.bins - 1))
    }

    fn center(&self, bin: usize) -> f64 {
        assert!(bin < self.bins, "axis bin out of range");
        self.min + (bin as f64 + 0.5) * self.width
    }
}

/// Sorted set of physically allowed scalar macrostates.
#[derive(Debug, Clone, PartialEq)]
pub struct DiscreteAxis {
    values: Vec<f64>,
    absolute_tolerance: f64,
    relative_tolerance: f64,
}

impl DiscreteAxis {
    /// Construct a discrete axis with conservative energy-cache tolerances.
    pub fn new(values: Vec<f64>) -> Result<Self, GeneralizedError> {
        Self::with_tolerance(values, 1e-10, 1e-12)
    }

    /// Construct a discrete axis with explicit absolute and relative matching tolerances.
    pub fn with_tolerance(
        mut values: Vec<f64>,
        absolute_tolerance: f64,
        relative_tolerance: f64,
    ) -> Result<Self, GeneralizedError> {
        if values.is_empty() {
            return Err(GeneralizedError::new(
                "discrete axis must contain at least one value",
            ));
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(GeneralizedError::new(
                "discrete axis contains a non-finite value",
            ));
        }
        if !absolute_tolerance.is_finite()
            || absolute_tolerance < 0.0
            || !relative_tolerance.is_finite()
            || relative_tolerance < 0.0
        {
            return Err(GeneralizedError::new(
                "axis tolerances must be finite and non-negative",
            ));
        }
        values.sort_by(f64::total_cmp);
        for pair in values.windows(2) {
            if close(pair[0], pair[1], absolute_tolerance, relative_tolerance) {
                return Err(GeneralizedError::new(
                    "discrete axis contains duplicate values within tolerance",
                ));
            }
        }
        Ok(Self {
            values,
            absolute_tolerance,
            relative_tolerance,
        })
    }

    #[inline]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    #[inline]
    pub const fn absolute_tolerance(&self) -> f64 {
        self.absolute_tolerance
    }

    #[inline]
    pub const fn relative_tolerance(&self) -> f64 {
        self.relative_tolerance
    }
}

impl MacrostateAxis for DiscreteAxis {
    #[inline]
    fn bins(&self) -> usize {
        self.values.len()
    }

    fn bin(&self, value: f64) -> Option<usize> {
        if !value.is_finite() {
            return None;
        }
        match self
            .values
            .binary_search_by(|candidate| candidate.total_cmp(&value))
        {
            Ok(index) => Some(index),
            Err(index) => {
                let left = index.checked_sub(1);
                let right = (index < self.values.len()).then_some(index);
                [left, right]
                    .into_iter()
                    .flatten()
                    .min_by(|&first, &second| {
                        (self.values[first] - value)
                            .abs()
                            .total_cmp(&(self.values[second] - value).abs())
                    })
                    .filter(|&candidate| {
                        close(
                            self.values[candidate],
                            value,
                            self.absolute_tolerance,
                            self.relative_tolerance,
                        )
                    })
            }
        }
    }

    fn center(&self, bin: usize) -> f64 {
        self.values[bin]
    }
}

#[inline]
fn close(left: f64, right: f64, absolute: f64, relative: f64) -> bool {
    (left - right).abs() <= absolute + relative * left.abs().max(right.abs())
}
