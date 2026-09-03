//! Measurement accumulation with binning.
//!
//! This module provides [`Measurements`] for accumulating observable samples
//! during Monte Carlo simulations, with automatic binning to reduce memory usage.
//!
//! # Scalar and Array Observables
//!
//! [`Accumulator`] supports both scalar (`f64`) and array (`ndarray::ArrayD<f64>`)
//! observables. Shape is determined by the first sample added.
//!
//! # Binning
//!
//! Samples are grouped into bins of configurable size. Each bin stores only
//! the mean of its samples, reducing memory for large simulations while
//! preserving statistical accuracy.
//!
//! # Usage
//!
//! ```text
//! // Scalar observable
//! ctx.measure("Energy", energy_value);
//!
//! // Array observable (e.g. correlation function)
//! let corr: Vec<f64> = compute_correlation();
//! ctx.measure_array("Correlation", &corr);
//! ```
//!
//! # Finalization
//!
//! After simulation completes, [`Measurements::finalize()`] returns estimates
//! for each observable with mean and error.

use crate::estimate::ComplexEstimate;
use crate::Estimate;
use ndarray::{Array1, ArrayD};
use std::collections::HashMap;

#[cfg(feature = "hdf5")]
use hdf5::Group;

/// Complex number representation for observables.
/// Stored as separate real/imaginary parts matching Carlo.jl HDF5 format.
#[derive(Debug, Clone)]
pub struct ComplexValue {
    pub re: f64,
    pub im: f64,
}

impl ComplexValue {
    /// Create a complex observable sample from real and imaginary parts.
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
}

/// Single observable accumulator with binning.
///
/// Supports both scalar (shape `[]`) and array (shape `[N]`, `[N,M]`, ...) observables.
/// Shape is determined by the first sample added.
pub struct Accumulator {
    /// Completed bins (each bin is mean of bin_capacity samples).
    bins: Vec<Array1<f64>>,

    /// Current bin being filled (accumulating sum).
    current_bin: Array1<f64>,

    /// Number of samples per bin.
    bin_capacity: usize,

    /// Samples accumulated in current bin.
    current_filling: usize,

    /// Shape of the observable (empty = scalar).
    shape: Vec<usize>,

    /// Total sample count.
    total_count: usize,
}

impl Accumulator {
    /// Create a new scalar accumulator.
    pub fn new(bin_capacity: usize) -> Self {
        Self {
            bins: Vec::new(),
            current_bin: Array1::zeros(0),
            bin_capacity,
            current_filling: 0,
            shape: vec![],
            total_count: 0,
        }
    }

    /// Create an array accumulator with given shape.
    pub fn with_shape(bin_capacity: usize, shape: &[usize]) -> Self {
        let flat_size: usize = shape.iter().product();
        Self {
            bins: Vec::new(),
            current_bin: Array1::zeros(flat_size),
            bin_capacity,
            current_filling: 0,
            shape: shape.to_vec(),
            total_count: 0,
        }
    }

    /// Add a scalar sample.
    pub fn add(&mut self, value: f64) {
        // If this is the first sample and we're scalar-shaped, init
        if self.shape.is_empty() && self.current_bin.is_empty() {
            self.current_bin = Array1::zeros(1);
        }

        self.current_bin[0] += value;
        self.total_count += 1;
        self.current_filling += 1;

        if self.current_filling >= self.bin_capacity {
            let bin_mean = self.current_bin[0] / self.bin_capacity as f64;
            self.bins.push(Array1::from_vec(vec![bin_mean]));
            self.current_bin.fill(0.0);
            self.current_filling = 0;
        }
    }

    /// Add an array sample (flat slice).
    /// Shape is determined by the first call.
    pub fn add_array(&mut self, values: &[f64]) {
        if values.is_empty() {
            return;
        }

        if self.shape.is_empty() {
            // First sample: determine shape
            self.shape = vec![values.len()];
            self.current_bin = Array1::zeros(values.len());
        }

        debug_assert_eq!(
            values.len(),
            self.current_bin.len(),
            "Shape mismatch: expected {}, got {}",
            self.current_bin.len(),
            values.len()
        );

        for (dst, &src) in self.current_bin.iter_mut().zip(values.iter()) {
            *dst += src;
        }
        self.total_count += 1;
        self.current_filling += 1;

        if self.current_filling >= self.bin_capacity {
            // Normalize completed bin
            let scale = 1.0 / self.bin_capacity as f64;
            self.bins.push(self.current_bin.clone() * scale);
            self.current_bin.fill(0.0);
            self.current_filling = 0;
        }
    }

    /// Finalize and return estimate.
    pub fn finalize(&self) -> Estimate {
        let mut all_bins: Vec<f64> = Vec::new();

        // Completed bins: compute mean of flattened bin values
        for bin in &self.bins {
            if bin.is_empty() {
                all_bins.push(0.0);
            } else {
                all_bins.push(bin.sum() / bin.len() as f64);
            }
        }

        // Partial bin
        if self.current_filling > 0 && !self.current_bin.is_empty() {
            let partial_mean = self.current_bin.sum()
                / (self.current_filling as f64 * self.current_bin.len() as f64);
            all_bins.push(partial_mean);
        }

        Estimate::from_bins_with_autocorr(&all_bins, self.autocorr_time())
    }

    /// Check if any complete bins exist.
    pub fn has_complete_bins(&self) -> bool {
        !self.bins.is_empty()
    }

    /// Get bin capacity.
    pub fn bin_capacity(&self) -> usize {
        self.bin_capacity
    }

    /// Get completed bins.
    pub fn bins(&self) -> &[Array1<f64>] {
        &self.bins
    }

    /// Get shape of observable (empty = scalar).
    pub fn shape(&self) -> &[usize] {
        &self.shape
    }

    /// Get total sample count.
    pub fn total_count(&self) -> usize {
        self.total_count
    }

    /// Get number of bins.
    pub fn num_bins(&self) -> usize {
        self.bins.len()
    }

    /// Get total samples including partial bin.
    pub fn total_samples(&self) -> usize {
        self.total_count
    }

    /// Get rebin means as flat array for autocorrelation analysis.
    pub fn rebin_means(&self) -> Vec<f64> {
        self.bins.iter().map(|b| b.sum() / b.len() as f64).collect()
    }

    /// Get all bin data as Array2 for covariance estimation.
    /// Shape: (n_bins, flat_size)
    pub fn bin_matrix(&self) -> ArrayD<f64> {
        if self.bins.is_empty() {
            return ArrayD::zeros(vec![0]);
        }

        let n_bins = self.bins.len();
        let flat_size = self.bins[0].len();

        let mut data = Vec::with_capacity(n_bins * flat_size);
        for bin in &self.bins {
            data.extend(bin.iter().copied());
        }

        ArrayD::from_shape_vec(vec![n_bins, flat_size], data)
            .unwrap_or_else(|_| ArrayD::zeros(vec![0]))
    }

    /// Compute autocorrelation time for scalar observables.
    pub fn autocorr_time(&self) -> f64 {
        if self.bins.len() < 2 {
            return 1.0;
        }

        let means: Vec<f64> = self.bins.iter().map(|b| b.sum() / b.len() as f64).collect();

        let n = means.len() as f64;
        let mu = means.iter().sum::<f64>() / n;

        // Variance without rebinning
        let var_no_rebin: f64 = means.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / (n - 1.0);

        // Rebin: pair consecutive samples
        let rebin_n = n as usize / 2;
        if rebin_n < 2 {
            return 1.0;
        }

        let rebin_means: Vec<f64> = (0..rebin_n)
            .map(|i| (means[2 * i] + means[2 * i + 1]) / 2.0)
            .collect();
        let rebin_mu = rebin_means.iter().sum::<f64>() / rebin_n as f64;
        let var_rebin: f64 = rebin_means
            .iter()
            .map(|x| (x - rebin_mu).powi(2))
            .sum::<f64>()
            / (rebin_n as f64 - 1.0);

        // Expected variance without correlations: var_no_rebin / 2
        let expected_var = var_no_rebin / 2.0;
        if expected_var <= 0.0 {
            return 1.0;
        }

        0.5 * (var_rebin / expected_var - 1.0).max(0.0)
    }

    /// Compute covariance matrix from bins.
    /// Returns None for scalar observables or if not enough bins.
    pub fn covariance(&self) -> Option<ArrayD<f64>> {
        if self.shape.is_empty() || self.bins.len() < 3 {
            return None;
        }

        let n_bins = self.bins.len();
        let flat_size = self.shape.iter().product::<usize>();

        if flat_size == 0 {
            return None;
        }

        // Compute mean
        let mut mean = Array1::zeros(flat_size);
        for bin in &self.bins {
            mean += bin;
        }
        mean /= n_bins as f64;

        // Compute covariance: Σ = 1/(N-1) * Σ_i (x_i - μ)(x_i - μ)^T
        let mut cov = ArrayD::zeros(vec![flat_size, flat_size]);

        for bin in &self.bins {
            let diff = bin - &mean;
            for i in 0..flat_size {
                for j in 0..flat_size {
                    cov[[i, j]] += diff[i] * diff[j];
                }
            }
        }

        let scale = 1.0 / (n_bins as f64 - 1.0);
        cov *= scale;

        Some(cov)
    }

    /// Compute autocorrelation time from stored bins using the variance ratio method.
    /// Uses 2-sample rebinning for a basic estimate.
    pub fn autocorr_time_from_bins(&self) -> f64 {
        let n = self.bins.len();
        if n < 4 {
            return 1.0;
        }

        let means: Vec<f64> = self.bins.iter().map(|b| b[0]).collect();
        let n_f64 = n as f64;
        let mu = means.iter().sum::<f64>() / n_f64;

        // Variance without rebinning
        let var_no_rebin = means.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / (n_f64 - 1.0);

        // Rebin: pair consecutive samples
        let rebin_n = n / 2;
        if rebin_n < 2 {
            return 1.0;
        }

        let rebin_means: Vec<f64> = (0..rebin_n)
            .map(|i| (means[2 * i] + means[2 * i + 1]) / 2.0)
            .collect();
        let rebin_mu = rebin_means.iter().sum::<f64>() / rebin_n as f64;
        let var_rebin = rebin_means
            .iter()
            .map(|x| (x - rebin_mu).powi(2))
            .sum::<f64>()
            / (rebin_n as f64 - 1.0);

        let expected_var = var_no_rebin / 2.0;
        if expected_var <= 0.0 {
            return 0.0;
        }

        0.5 * (var_rebin / expected_var - 1.0).max(0.0)
    }
}

/// Complex observable accumulator with binning.
/// Stores real and imaginary parts as separate scalar accumulators,
/// matching Carlo.jl HDF5 format (re/im split).
pub struct ComplexAccumulator {
    /// Real part bins.
    re_bins: Vec<f64>,
    /// Imaginary part bins.
    im_bins: Vec<f64>,
    /// Sum of real parts in current bin.
    re_current: f64,
    /// Sum of imaginary parts in current bin.
    im_current: f64,
    /// Number of samples per bin.
    bin_capacity: usize,
    /// Samples in current bin.
    current_filling: usize,
    /// Total sample count.
    total_count: usize,
}

impl ComplexAccumulator {
    /// Create a new complex accumulator.
    pub fn new(bin_capacity: usize) -> Self {
        Self {
            re_bins: Vec::new(),
            im_bins: Vec::new(),
            re_current: 0.0,
            im_current: 0.0,
            bin_capacity,
            current_filling: 0,
            total_count: 0,
        }
    }

    /// Add a complex sample.
    pub fn add(&mut self, re: f64, im: f64) {
        self.re_current += re;
        self.im_current += im;
        self.total_count += 1;
        self.current_filling += 1;

        if self.current_filling >= self.bin_capacity {
            let scale = 1.0 / self.bin_capacity as f64;
            self.re_bins.push(self.re_current * scale);
            self.im_bins.push(self.im_current * scale);
            self.re_current = 0.0;
            self.im_current = 0.0;
            self.current_filling = 0;
        }
    }

    /// Finalize and return complex estimate.
    pub fn finalize(&self) -> crate::estimate::ComplexEstimate {
        use crate::estimate::ComplexEstimate;
        ComplexEstimate::new(
            Estimate::from_bins(&self.re_bins),
            Estimate::from_bins(&self.im_bins),
        )
    }

    /// Get number of completed bins.
    pub fn num_bins(&self) -> usize {
        self.re_bins.len()
    }

    /// Get total samples.
    pub fn total_count(&self) -> usize {
        self.total_count
    }

    /// Get rebin means for real and imaginary parts.
    pub fn rebin_means(&self) -> (Vec<f64>, Vec<f64>) {
        (self.re_bins.clone(), self.im_bins.clone())
    }
}

/// Measurement collector managing multiple observables.
pub struct Measurements {
    observables: HashMap<String, Accumulator>,
    complex_observables: HashMap<String, ComplexAccumulator>,
    default_binsize: usize,
}

impl Measurements {
    /// Create a collector with `default_binsize` for auto-registered observables.
    pub fn new(default_binsize: usize) -> Self {
        Self {
            observables: HashMap::new(),
            complex_observables: HashMap::new(),
            default_binsize,
        }
    }

    /// Add a scalar sample to an observable (auto-create if needed).
    pub fn add_sample(&mut self, name: &str, value: f64) {
        if !self.observables.contains_key(name) {
            self.observables
                .insert(name.to_string(), Accumulator::new(self.default_binsize));
        }
        self.observables
            .get_mut(name)
            .expect("observable just inserted")
            .add(value);
    }

    /// Add an array sample to an observable.
    pub fn add_sample_array(&mut self, name: &str, values: &[f64]) {
        if !self.observables.contains_key(name) {
            self.observables.insert(
                name.to_string(),
                Accumulator::with_shape(self.default_binsize, &[values.len()]),
            );
        }
        self.observables
            .get_mut(name)
            .expect("observable just inserted")
            .add_array(values);
    }

    /// Register a scalar observable with custom binsize.
    pub fn register(&mut self, name: &str, binsize: usize) {
        self.observables
            .insert(name.to_string(), Accumulator::new(binsize));
    }

    /// Register an array observable with custom binsize and shape.
    pub fn register_array(&mut self, name: &str, binsize: usize, shape: &[usize]) {
        self.observables
            .insert(name.to_string(), Accumulator::with_shape(binsize, shape));
    }

    /// Add a complex sample to an observable.
    pub fn add_sample_complex(&mut self, name: &str, re: f64, im: f64) {
        if !self.complex_observables.contains_key(name) {
            self.complex_observables.insert(
                name.to_string(),
                ComplexAccumulator::new(self.default_binsize),
            );
        }
        self.complex_observables
            .get_mut(name)
            .expect("observable just inserted")
            .add(re, im);
    }

    /// Register a complex observable with custom binsize.
    pub fn register_complex(&mut self, name: &str, binsize: usize) {
        self.complex_observables
            .insert(name.to_string(), ComplexAccumulator::new(binsize));
    }

    /// Finalize all observables and return estimates.
    pub fn finalize(&self) -> HashMap<String, Estimate> {
        self.observables
            .iter()
            .map(|(name, acc)| (name.clone(), acc.finalize()))
            .collect()
    }

    /// Finalize complex observables and return estimates.
    pub fn finalize_complex(&self) -> HashMap<String, ComplexEstimate> {
        self.complex_observables
            .iter()
            .map(|(name, acc)| (name.clone(), acc.finalize()))
            .collect()
    }

    /// Get observables map (for iteration).
    pub fn observables(&self) -> &HashMap<String, Accumulator> {
        &self.observables
    }

    /// Get complex observables map.
    pub fn complex_observables(&self) -> &HashMap<String, ComplexAccumulator> {
        &self.complex_observables
    }
}

#[cfg(feature = "hdf5")]
impl Accumulator {
    /// Write accumulator bins to HDF5 group.
    pub fn write_hdf5(&self, group: &mut Group, name: &str) -> Result<(), crate::CarloError> {
        let obs_group = group
            .create_group(name)
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot create group {}: {}", name, e),
            })?;

        obs_group
            .new_dataset_builder()
            .with_data(&[self.bin_capacity as u64])
            .create("bin_length")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot write bin_length: {}", e),
            })?;

        // Write shape
        let shape_u64: Vec<u64> = self.shape.iter().map(|&x| x as u64).collect();
        obs_group
            .new_dataset_builder()
            .with_data(&shape_u64)
            .create("shape")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot write shape: {}", e),
            })?;

        // Write bins as 2D array: (flat_size, n_bins)
        if !self.bins.is_empty() {
            let flat_size = self.bins[0].len();
            let mut data: Vec<f64> = Vec::with_capacity(flat_size * self.bins.len());
            for bin in &self.bins {
                data.extend(bin.iter().copied());
            }

            obs_group
                .new_dataset_builder()
                .with_data(&data)
                .create("samples")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot write samples: {}", e),
                })?;
        }

        Ok(())
    }

    /// Read accumulator from HDF5 group.
    pub fn read_hdf5(group: &Group) -> Result<Self, crate::CarloError> {
        let bin_length: u64 = group
            .dataset("bin_length")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot read bin_length: {}", e),
            })?
            .read_1d::<u64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot parse bin_length: {}", e),
            })?[0];

        let shape: Vec<usize> = group
            .dataset("shape")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot read shape: {}", e),
            })?
            .read_1d::<u64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot parse shape: {}", e),
            })?
            .iter()
            .map(|&x| x as usize)
            .collect();

        let bins: Vec<Array1<f64>> = if shape.is_empty() {
            // Scalar case
            group
                .dataset("samples")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot read samples: {}", e),
                })?
                .read_1d::<f64>()
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot parse samples: {}", e),
                })?
                .to_vec()
                .into_iter()
                .map(|v| Array1::from_vec(vec![v]))
                .collect()
        } else {
            // Array case
            let flat_size: usize = shape.iter().product();
            let data: Vec<f64> = group
                .dataset("samples")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot read samples: {}", e),
                })?
                .read_1d::<f64>()
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot parse samples: {}", e),
                })?
                .to_vec();

            data.chunks(flat_size)
                .map(|chunk| Array1::from_vec(chunk.to_vec()))
                .collect()
        };

        let total_count = bins.len() * bin_length as usize;

        Ok(Self {
            bins,
            current_bin: Array1::zeros(if shape.is_empty() {
                0
            } else {
                shape.iter().product()
            }),
            bin_capacity: bin_length as usize,
            current_filling: 0,
            shape,
            total_count,
        })
    }

    /// Write accumulator checkpoint (includes partial bin).
    pub fn write_checkpoint_hdf5(
        &self,
        group: &mut Group,
        name: &str,
    ) -> Result<(), crate::CarloError> {
        let obs_group = group
            .create_group(name)
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot create group {}: {}", name, e),
            })?;

        // Write bin capacity
        obs_group
            .new_dataset_builder()
            .with_data(&[self.bin_capacity as u64])
            .create("bin_capacity")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write bin_capacity: {}", e),
            })?;

        // Write shape
        let shape_u64: Vec<u64> = self.shape.iter().map(|&x| x as u64).collect();
        obs_group
            .new_dataset_builder()
            .with_data(&shape_u64)
            .create("shape")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write shape: {}", e),
            })?;

        // Write completed bins
        if !self.bins.is_empty() {
            let flat_size = self.bins[0].len();
            let mut data: Vec<f64> = Vec::with_capacity(flat_size * self.bins.len());
            for bin in &self.bins {
                data.extend(bin.iter().copied());
            }

            obs_group
                .new_dataset_builder()
                .with_data(&data)
                .create("bins")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot write bins: {}", e),
                })?;
        }

        // Write partial bin (current_bin, only up to current_filling)
        if self.current_filling > 0 && !self.current_bin.is_empty() {
            // current_bin contains sum of current_filling samples
            let scale = 1.0 / self.current_filling as f64;
            let partial_mean: Vec<f64> = self.current_bin.iter().map(|x| x * scale).collect();

            obs_group
                .new_dataset_builder()
                .with_data(&partial_mean)
                .create("partial_bin_mean")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot write partial_bin_mean: {}", e),
                })?;

            obs_group
                .new_dataset_builder()
                .with_data(&[self.current_filling as u64])
                .create("partial_bin_filling")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot write partial_bin_filling: {}", e),
                })?;
        }

        // Write total count
        obs_group
            .new_dataset_builder()
            .with_data(&[self.total_count as u64])
            .create("total_count")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write total_count: {}", e),
            })?;

        Ok(())
    }

    /// Read accumulator checkpoint (includes partial bin).
    pub fn read_checkpoint_hdf5(group: &Group) -> Result<Self, crate::CarloError> {
        let bin_capacity: usize = group
            .dataset("bin_capacity")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read bin_capacity: {}", e),
            })?
            .read_1d::<u64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse bin_capacity: {}", e),
            })?[0] as usize;

        let shape: Vec<usize> = group
            .dataset("shape")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read shape: {}", e),
            })?
            .read_1d::<u64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse shape: {}", e),
            })?
            .iter()
            .map(|&x| x as usize)
            .collect();

        let flat_size: usize = if shape.is_empty() {
            0
        } else {
            shape.iter().product()
        };

        let bins: Vec<Array1<f64>> = if let Ok(ds) = group.dataset("bins") {
            let data: Vec<f64> = ds
                .read_1d::<f64>()
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot parse bins: {}", e),
                })?
                .to_vec();

            let effective_flat = if flat_size > 0 {
                flat_size
            } else {
                // Scalar: each element is one bin
                data.len()
            };

            if effective_flat == 0 {
                Vec::new()
            } else {
                data.chunks(effective_flat)
                    .map(|chunk| Array1::from_vec(chunk.to_vec()))
                    .collect()
            }
        } else {
            Vec::new()
        };

        let total_count: usize = group
            .dataset("total_count")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read total_count: {}", e),
            })?
            .read_1d::<u64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse total_count: {}", e),
            })?[0] as usize;

        // Restore partial bin if present
        let (current_bin, current_filling) = if let Ok(ds) = group.dataset("partial_bin_mean") {
            let partial_mean: Vec<f64> = ds
                .read_1d::<f64>()
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot parse partial_bin_mean: {}", e),
                })?
                .to_vec();

            let filling: u64 = group
                .dataset("partial_bin_filling")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot read partial_bin_filling: {}", e),
                })?
                .read_1d::<u64>()
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot parse partial_bin_filling: {}", e),
                })?[0];

            // Convert mean back to sum
            let sum = partial_mean
                .iter()
                .map(|x| x * filling as f64)
                .collect::<Vec<f64>>();

            (Array1::from_vec(sum), filling as usize)
        } else {
            (Array1::zeros(flat_size), 0)
        };

        Ok(Self {
            bins,
            current_bin,
            bin_capacity,
            current_filling,
            shape,
            total_count,
        })
    }
}

#[cfg(feature = "hdf5")]
impl Measurements {
    /// Write all measurements to HDF5 file.
    pub fn write_hdf5(&self, file: &mut Group) -> Result<(), crate::CarloError> {
        let mut obs_group =
            file.create_group("observables")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot create observables group: {}", e),
                })?;
        for (name, acc) in &self.observables {
            acc.write_hdf5(&mut obs_group, name)?;
        }
        for (name, acc) in &self.complex_observables {
            acc.write_hdf5(&mut obs_group, name)?;
        }
        Ok(())
    }

    /// Read measurements from HDF5 file.
    pub fn read_hdf5(file: &Group) -> Result<Self, crate::CarloError> {
        let mut measurements = Self::new(100); // default binsize
        let obs_group =
            file.group("observables")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot open observables group: {}", e),
                })?;

        for name_result in
            obs_group
                .member_names()
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot list members: {}", e),
                })?
        {
            let name = name_result;
            if let Ok(obs) = obs_group.group(&name) {
                // Check if this is a complex observable
                let is_complex = obs
                    .dataset("is_complex")
                    .ok()
                    .and_then(|ds| ds.read_1d::<u64>().ok().map(|arr| arr[0] == 1))
                    .unwrap_or(false);

                if is_complex {
                    if let Ok(acc) = ComplexAccumulator::read_hdf5(&obs) {
                        measurements.complex_observables.insert(name, acc);
                    }
                } else if let Ok(acc) = Accumulator::read_hdf5(&obs) {
                    measurements.observables.insert(name, acc);
                }
            }
        }
        Ok(measurements)
    }

    /// Write measurements checkpoint (includes partial bins).
    pub fn write_checkpoint_hdf5(&self, group: &mut Group) -> Result<(), crate::CarloError> {
        // Write default binsize
        group
            .new_dataset_builder()
            .with_data(&[self.default_binsize as u64])
            .create("default_binsize")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot write default_binsize: {}", e),
            })?;

        // Write real observables with full state
        let mut obs_group =
            group
                .create_group("observables")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot create observables group: {}", e),
                })?;

        for (name, acc) in &self.observables {
            acc.write_checkpoint_hdf5(&mut obs_group, name)?;
        }

        // Write complex observables with full state
        let mut complex_group = group.create_group("complex_observables").map_err(|e| {
            crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot create complex_observables group: {}", e),
            }
        })?;

        for (name, acc) in &self.complex_observables {
            acc.write_hdf5(&mut complex_group, name)?;
        }

        Ok(())
    }

    /// Read measurements checkpoint (includes partial bins).
    pub fn read_checkpoint_hdf5(group: &Group) -> Result<Self, crate::CarloError> {
        let default_binsize: usize = group
            .dataset("default_binsize")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot read default_binsize: {}", e),
            })?
            .read_1d::<u64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "checkpoint".into(),
                reason: format!("Cannot parse default_binsize: {}", e),
            })?[0] as usize;

        let mut measurements = Self::new(default_binsize);

        let obs_group =
            group
                .group("observables")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "checkpoint".into(),
                    reason: format!("Cannot open observables group: {}", e),
                })?;

        for name_result in obs_group.member_names().unwrap_or_default() {
            let name = name_result;
            if let Ok(obs) = obs_group.group(&name) {
                if let Ok(acc) = Accumulator::read_checkpoint_hdf5(&obs) {
                    measurements.observables.insert(name, acc);
                }
            }
        }

        // Read complex observables
        if let Ok(complex_group) = group.group("complex_observables") {
            for name_result in complex_group.member_names().unwrap_or_default() {
                let name = name_result;
                if let Ok(obs) = complex_group.group(&name) {
                    if let Ok(acc) = ComplexAccumulator::read_hdf5(&obs) {
                        measurements.complex_observables.insert(name, acc);
                    }
                }
            }
        }

        Ok(measurements)
    }
}

#[cfg(feature = "hdf5")]
impl ComplexAccumulator {
    /// Write complex accumulator to HDF5 group.
    /// Stores re/im parts separately matching Carlo.jl format.
    pub fn write_hdf5(&self, group: &mut Group, name: &str) -> Result<(), crate::CarloError> {
        let obs_group = group
            .create_group(name)
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot create group {}: {}", name, e),
            })?;

        obs_group
            .new_dataset_builder()
            .with_data(&[self.bin_capacity as u64])
            .create("bin_length")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot write bin_length: {}", e),
            })?;

        // Mark as complex observable
        obs_group
            .new_dataset_builder()
            .with_data(&[1u64])
            .create("is_complex")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot write is_complex: {}", e),
            })?;

        // Write real part samples
        if !self.re_bins.is_empty() {
            obs_group
                .new_dataset_builder()
                .with_data(&self.re_bins)
                .create("samples_re")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot write samples_re: {}", e),
                })?;
            obs_group
                .new_dataset_builder()
                .with_data(&self.im_bins)
                .create("samples_im")
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot write samples_im: {}", e),
                })?;
        }

        Ok(())
    }

    /// Read complex accumulator from HDF5 group.
    pub fn read_hdf5(group: &Group) -> Result<Self, crate::CarloError> {
        let bin_length: u64 = group
            .dataset("bin_length")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot read bin_length: {}", e),
            })?
            .read_1d::<u64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot parse bin_length: {}", e),
            })?[0];

        let re_bins: Vec<f64> = group
            .dataset("samples_re")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot read samples_re: {}", e),
            })?
            .read_1d::<f64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot parse samples_re: {}", e),
            })?
            .to_vec();

        let im_bins: Vec<f64> = group
            .dataset("samples_im")
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot read samples_im: {}", e),
            })?
            .read_1d::<f64>()
            .map_err(|e| crate::CarloError::InvalidConfig {
                field: "hdf5".into(),
                reason: format!("Cannot parse samples_im: {}", e),
            })?
            .to_vec();

        let total_count = re_bins.len() * bin_length as usize;

        Ok(Self {
            re_bins,
            im_bins,
            re_current: 0.0,
            im_current: 0.0,
            bin_capacity: bin_length as usize,
            current_filling: 0,
            total_count,
        })
    }
}
