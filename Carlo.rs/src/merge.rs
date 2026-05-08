//! Result merging and rebinning analysis.
//!
//! After Monte Carlo simulations complete, this module provides tools for
//! combining results from multiple runs and estimating statistical errors.
//!
//! # Key Functions
//!
//! - [`calc_rebin_count()`]: Determine optimal rebin count
//! - [`calc_rebin_length()`]: Calculate rebin length from samples
//! - [`compute_regular_autocorr_time()`]: Estimate autocorrelation time
//! - [`cov_of_mean()`]: Compute covariance tensor of the mean
//! - [`list_meas_files()`]: List measurement files in a directory
//!
//! ## HDF5 Functions (requires `hdf5` feature)
//!
//! When compiled with `--features hdf5`:
//! - `merge_results()`: Merge HDF5 measurement files from a task directory
//! - `merge_results_from_files()`: Merge from explicit file list
//! - `iterate_measfile_observables()`: Iterate over observables in files
//!
//! # Rebinning
//!
//! Rebinning reduces correlations between samples by grouping consecutive
//! measurements. The rebin length is chosen automatically based on sample
//! count or can be specified manually via [`MergeOptions`].
//!
//! # Example
//!
//! ```rust
//! use carlo_rs::merge::{MergeOptions, calc_rebin_count, calc_rebin_length};
//!
//! let sample_count = 10000u64;
//! let rebin_count = calc_rebin_count(sample_count, 10);
//! let rebin_length = calc_rebin_length(sample_count, None);
//! ```

use ndarray::ArrayD;
use serde::Serialize;
use std::path::PathBuf;

#[cfg(feature = "hdf5")]
use std::collections::HashMap;

#[cfg(feature = "hdf5")]
use hdf5::File as Hdf5File;

/// Determine the number of bins in the rebin procedure.
/// Rebinning will not be performed if sample_count <= min_bin_count.
pub fn calc_rebin_count(sample_count: u64, min_bin_count: u64) -> u64 {
    if sample_count <= min_bin_count {
        sample_count
    } else {
        min_bin_count + ((sample_count - min_bin_count) as f64).cbrt().round() as u64
    }
}

/// Calculate the rebin length from total samples.
pub fn calc_rebin_length(total_sample_count: u64, rebin_length: Option<u64>) -> u64 {
    if total_sample_count == 0 {
        1
    } else if let Some(len) = rebin_length {
        len
    } else {
        total_sample_count / calc_rebin_count(total_sample_count, 10)
    }
}

/// Observable metadata from HDF5 measurement files.
#[derive(Debug, Clone)]
pub struct ObservableType<T> {
    /// Internal bin length from measurement.
    pub internal_bin_length: u64,

    /// Shape of the observable (excluding sample dimension).
    pub shape: Vec<usize>,

    /// Total sample count across all runs.
    pub total_sample_count: u64,

    /// Phantom data for type tracking.
    _type: std::marker::PhantomData<T>,
}

impl<T> ObservableType<T> {
    pub fn new(internal_bin_length: u64, shape: Vec<usize>, total_sample_count: u64) -> Self {
        Self {
            internal_bin_length,
            shape,
            total_sample_count,
            _type: std::marker::PhantomData,
        }
    }
}

/// Merged observable with statistics.
///
/// JSON serialization mirrors Carlo.jl's `JSON.lower(::ResultObservable)`:
/// `mean`, `error`, `rebin_len`, `autocorr_time` (scalar max), `rebin_count`,
/// `internal_bin_len`, and optionally `covariance`.
#[derive(Debug, Clone)]
pub struct ResultObservable<T> {
    /// Internal bin length from measurement.
    pub internal_bin_length: u64,

    /// Rebin length used in analysis.
    pub rebin_length: u64,

    /// Mean value.
    pub mean: ArrayD<T>,

    /// Standard error of the mean.
    pub error: ArrayD<T>,

    /// Covariance matrix (optional, for array observables).
    pub covariance: Option<ArrayD<T>>,

    /// Autocorrelation time estimate.
    pub autocorrelation_time: ArrayD<f64>,

    /// Rebin means for jackknife analysis.
    pub rebin_means: ArrayD<T>,
}

/// Number of rebinned bins = last dimension of rebin_means.
pub fn rebin_count<T>(obs: &ResultObservable<T>) -> u64 {
    obs.rebin_means.shape().last().copied().unwrap_or(0) as u64
}

impl<T: Serialize> Serialize for ResultObservable<T> {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;

        let autocorr_scalar = self
            .autocorrelation_time
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max);

        let n_fields = if self.covariance.is_some() { 7 } else { 6 };
        let mut st = s.serialize_struct("ResultObservable", n_fields)?;
        st.serialize_field("mean", &self.mean)?;
        st.serialize_field("error", &self.error)?;
        st.serialize_field("rebin_len", &self.rebin_length)?;
        st.serialize_field("autocorr_time", &autocorr_scalar)?;
        st.serialize_field("rebin_count", &rebin_count(self))?;
        st.serialize_field("internal_bin_len", &self.internal_bin_length)?;
        if let Some(ref cov) = self.covariance {
            st.serialize_field("covariance", cov)?;
        }
        st.end()
    }
}

// Deserialize is not implemented — ResultObservable is written in Carlo.jl
// JSON format and read back via `measurement_from_obs()` in resulttools.rs.
// Use the manual JSON reader there instead of serde Deserialize.

/// Compute regular autocorrelation time from variance ratio.
/// τ = 0.5 * ((σ_binned / σ_unbinned)^2 - 1)
///
/// This definition gives τ = 0 for uncorrelated data.
pub fn compute_regular_autocorr_time(_mu: f64, sigma: f64, sigma_no_rebin: f64) -> f64 {
    if sigma_no_rebin <= 0.0 {
        return 0.0;
    }
    let ratio = sigma / sigma_no_rebin;
    0.5 * (ratio * ratio - 1.0).max(0.0)
}

/// Compute covariance tensor of the mean.
///
/// For an observable with shape S and N bins, returns a tensor of shape [S..., S...].
/// For a 1D observable with n components, returns an n×n covariance matrix.
///
/// Formula: cov[i,j] = (1/(N*(N-1))) * Σ_k (X_i[k] - μ_i)(X_j[k] - μ_j)
pub fn cov_of_mean(bins: &ArrayD<f64>) -> ArrayD<f64> {
    let ndim = bins.ndim();
    if ndim == 0 {
        return ArrayD::zeros(ndarray::IxDyn(&[1, 1]));
    }

    let obs_shape: Vec<usize> = bins.shape()[..ndim - 1].to_vec();
    let n_bins = bins.shape()[ndim - 1];

    if n_bins < 2 {
        let total_size: usize = obs_shape.iter().product();
        return ArrayD::zeros(
            std::iter::once(total_size)
                .chain(std::iter::once(total_size))
                .collect::<Vec<_>>(),
        );
    }

    // Flatten observation dimensions
    let obs_dim: usize = obs_shape.iter().product();
    let mean = bins.mean_axis(ndarray::Axis(ndim - 1)).unwrap();
    let mean_flat: Vec<f64> = mean.iter().copied().collect();

    let mut cov = ndarray::Array2::<f64>::zeros((obs_dim, obs_dim));

    for k in 0..n_bins {
        let col = bins.index_axis(ndarray::Axis(ndim - 1), k);
        let col_flat: Vec<f64> = col.iter().copied().collect();

        for i in 0..obs_dim {
            let di = col_flat[i] - mean_flat[i];
            for j in 0..obs_dim {
                let dj = col_flat[j] - mean_flat[j];
                cov[[i, j]] += di * dj;
            }
        }
    }

    let n = n_bins as f64;
    cov.mapv_inplace(|v| v / (n * (n - 1.0)));

    let new_shape: Vec<usize> = obs_shape
        .iter()
        .cloned()
        .chain(obs_shape.iter().cloned())
        .collect();
    cov.into_shape_with_order(new_shape).unwrap()
}

/// Compute decorrelated autocorrelation time using eigenvalue decomposition.
///
/// Implements the whitening transform approach from Carlo.jl:
/// 1. Compute unbinned covariance from binned samples
/// 2. Eigenvalue decompose: Σ = Q Λ Q^T
/// 3. Apply whitening: T = Λ^(-1/2) Q^T
/// 4. Transform binned covariance to decorrelated basis
/// 5. Compute per-component correlation times
pub fn compute_decorrelated_autocorr_time(
    bins: &ArrayD<f64>,
    mu: &ArrayD<f64>,
    binned_cov: &ArrayD<f64>,
    _total_samples: usize,
) -> ArrayD<f64> {
    let ndim = bins.ndim();
    let obs_shape: Vec<usize> = bins.shape()[..ndim - 1].to_vec();
    let obs_dim: usize = obs_shape.iter().product();
    let n_bins = bins.shape()[ndim - 1];
    let m = n_bins as f64;

    if obs_dim == 0 || n_bins < 2 {
        return ArrayD::zeros(ndarray::IxDyn(&obs_shape));
    }

    // Flatten mu
    let mu_flat: Vec<f64> = mu.iter().copied().collect();

    // Compute unbinned covariance: Σ_unbinned = (1/(M-1)) * Σ_k (x_k - μ)(x_k - μ)^T
    let mut sigma_unbinned = ndarray::Array2::<f64>::zeros((obs_dim, obs_dim));

    for k in 0..n_bins {
        let x_k = bins.index_axis(ndarray::Axis(ndim - 1), k);
        let x_flat: Vec<f64> = x_k.iter().copied().collect();

        for i in 0..obs_dim {
            let di = x_flat[i] - mu_flat[i];
            for j in 0..obs_dim {
                sigma_unbinned[[i, j]] += di * (x_flat[j] - mu_flat[j]);
            }
        }
    }

    let m_f64 = n_bins as f64;
    if m_f64 > 1.0 {
        sigma_unbinned.mapv_inplace(|v| v / (m_f64 - 1.0));
    }

    // Eigenvalue decomposition using nalgebra (pure Rust, no BLAS dependency)
    use nalgebra::{DMatrix, SymmetricEigen};
    let eigen_matrix: DMatrix<f64> =
        DMatrix::from_row_iterator(obs_dim, obs_dim, sigma_unbinned.iter().cloned());
    let eigen_decomp = SymmetricEigen::new(eigen_matrix);
    let eigenvalues = eigen_decomp.eigenvalues;
    let eigenvectors = eigen_decomp.eigenvectors;

    // Whitening transform: T = Λ^(-1/2) Q^T
    let tolerance = 1e-10;
    let lambda_inv_sqrt: Vec<f64> = eigenvalues
        .iter()
        .map(|&lambda| {
            if lambda > tolerance {
                1.0 / lambda.sqrt()
            } else {
                0.0
            }
        })
        .collect();

    let transform: ndarray::Array2<f64> =
        ndarray::Array2::from_diag(&ndarray::Array1::from_vec(lambda_inv_sqrt)).dot(
            &ndarray::Array2::from_shape_fn((obs_dim, obs_dim), |(i, j)| eigenvectors[(j, i)]),
        );

    // Reshape binned_cov to 2D
    let binned_cov_2d: ndarray::Array2<f64> =
        if binned_cov.ndim() == 2 && binned_cov.shape()[0] == obs_dim {
            binned_cov
                .view()
                .into_shape_with_order((obs_dim, obs_dim))
                .unwrap()
                .to_owned()
        } else {
            ndarray::Array2::eye(obs_dim)
        };

    // Transform binned covariance: Σ_binned_decorr = T Σ_binned T^T
    let sigma_binned_decorr: ndarray::Array2<f64> =
        transform.dot(&binned_cov_2d).dot(&transform.t());

    // Extract diagonal (variances in decorrelated basis)
    let binned_variances_decorr: ndarray::Array1<f64> =
        ndarray::Array1::from_shape_fn(obs_dim, |i| sigma_binned_decorr[[i, i]].max(0.0));

    // τ = 0.5 * (M * σ²_decorr - 1)
    let correlation_times: ndarray::Array1<f64> =
        binned_variances_decorr.mapv(|v: f64| (0.5_f64 * (m * v - 1.0_f64)).max(0.0));

    correlation_times.into_shape_with_order(obs_shape).unwrap()
}

/// List measurement files in a task directory.
pub fn list_meas_files(taskdir: &PathBuf) -> Result<Vec<PathBuf>, crate::CarloError> {
    use std::fs;

    let mut files = Vec::new();
    if let Ok(entries) = fs::read_dir(taskdir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                if name.ends_with(".meas.h5") {
                    files.push(path);
                }
            }
        }
    }

    // Sort by name for consistent ordering
    files.sort();
    Ok(files)
}

/// Iterate over observables in HDF5 measurement files.
///
/// For each observable in each file, calls the provided function with:
/// - The observable name
/// - The HDF5 group for the observable
/// - Any accumulated state from previous calls
#[cfg(feature = "hdf5")]
pub fn iterate_measfile_observables<F, T>(
    filenames: &[PathBuf],
    mut f: F,
) -> Result<HashMap<String, T>, crate::CarloError>
where
    F: FnMut(&str, &hdf5::Group, Option<T>) -> Result<T, crate::CarloError>,
{
    let mut states = HashMap::new();

    for filename in filenames {
        let file = Hdf5File::open(filename).map_err(|e| crate::CarloError::InvalidConfig {
            field: "hdf5".into(),
            reason: format!("Cannot open {}: {}", filename.display(), e),
        })?;

        let obs_group =
            file.group("observables")
                .map_err(|_| crate::CarloError::InvalidConfig {
                    field: "observables".into(),
                    reason: format!("No observables group in {}", filename.display()),
                })?;

        for name_result in obs_group.member_names().unwrap_or_default() {
            if let Ok(obs_name) = name_result {
                if let Ok(obs) = obs_group.group(&obs_name) {
                    let state = states.remove(&obs_name);
                    let new_state = f(&obs_name, &obs, state)?;
                    states.insert(obs_name, new_state);
                }
            }
        }
    }

    Ok(states)
}

/// Options for merging results.
#[derive(Debug, Clone, Default)]
pub struct MergeOptions {
    /// Override rebin length (None = automatic).
    pub rebin_length: Option<u64>,

    /// Number of samples to skip at start.
    pub sample_skip: u64,

    /// Estimate covariance matrices.
    pub estimate_covariance: bool,
}

/// Merge results from a list of HDF5 measurement files.
#[cfg(feature = "hdf5")]
pub fn merge_results_from_files(
    filenames: &[PathBuf],
    options: &MergeOptions,
) -> Result<HashMap<String, ResultObservable<f64>>, crate::CarloError> {
    if filenames.is_empty() {
        return Ok(HashMap::new());
    }

    // First pass: collect observable types and counts
    let obs_types: HashMap<String, ObservableType<f64>> =
        iterate_measfile_observables(filenames, |name, group, state| {
            let bin_length: u64 = group.dataset("bin_length")?.read_1d::<u64>().map_err(|e| {
                crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot read bin_length for {}: {}", name, e),
                }
            })?[0];

            let samples =
                group
                    .dataset("samples")
                    .map_err(|e| crate::CarloError::InvalidConfig {
                        field: "hdf5".into(),
                        reason: format!("Cannot read samples for {}: {}", name, e),
                    })?;

            let shape = samples.shape()[..samples.shape().len() - 1].to_vec();
            let sample_count = samples.shape()[samples.shape().len() - 1] as u64;

            Ok(match state {
                Some(mut t) => {
                    t.total_sample_count += sample_count.saturating_sub(options.sample_skip);
                    t
                }
                None => ObservableType::new(
                    bin_length,
                    shape,
                    sample_count.saturating_sub(options.sample_skip),
                ),
            })
        })?;

    // Second pass: read actual samples and compute statistics
    let results: HashMap<String, ResultObservable<f64>> = obs_types
        .into_iter()
        .map(|(name, obs_type)| {
            let rebin_len = calc_rebin_length(obs_type.total_sample_count, options.rebin_length);
            match accumulate_and_compute(&name, &obs_type, filenames, options, rebin_len) {
                Ok(result) => (name, result),
                Err(e) => {
                    tracing::warn!("Failed to compute statistics for {}: {}", name, e);
                    ResultObservable {
                        internal_bin_length: obs_type.internal_bin_length,
                        rebin_length: rebin_len,
                        mean: ndarray::ArrayD::zeros(obs_type.shape.clone()),
                        error: ndarray::ArrayD::zeros(obs_type.shape.clone()),
                        covariance: None,
                        autocorrelation_time: ndarray::ArrayD::zeros(vec![1]),
                        rebin_means: ndarray::ArrayD::zeros(obs_type.shape),
                    }
                }
            }
        })
        .collect();

    Ok(results)
}

/// Accumulate samples from files and compute statistics for a single observable.
#[cfg(feature = "hdf5")]
fn accumulate_and_compute(
    name: &str,
    obs_type: &ObservableType<f64>,
    filenames: &[PathBuf],
    options: &MergeOptions,
    rebin_length: u64,
) -> Result<ResultObservable<f64>, crate::CarloError> {
    let mut state = AddSamplesState::new(
        obs_type.internal_bin_length,
        &obs_type.shape,
        rebin_length,
        options.estimate_covariance,
    );

    // Read samples from all files and accumulate rebin bins
    for filename in filenames {
        let file = Hdf5File::open(filename).map_err(|e| crate::CarloError::InvalidConfig {
            field: "hdf5".into(),
            reason: format!("Cannot open {}: {}", filename.display(), e),
        })?;

        let obs = file
            .group("observables")
            .map_err(|_| crate::CarloError::InvalidConfig {
                field: "observables".into(),
                reason: "No observables group".into(),
            })?
            .group(name)
            .map_err(|_| crate::CarloError::InvalidConfig {
                field: name.into(),
                reason: format!("Observable {} not found", name),
            })?;

        let samples: ArrayD<f64> =
            obs.dataset("samples")?
                .read()
                .map_err(|e| crate::CarloError::InvalidConfig {
                    field: "hdf5".into(),
                    reason: format!("Cannot read samples for {}: {}", name, e),
                })?;

        // Accumulate rebin bins from samples
        let n_samples = samples.shape().last().copied().unwrap_or(0);
        let offset = options.sample_skip as usize;
        let mut pos = offset;
        while pos + rebin_length as usize <= n_samples {
            state.add_rebin_bin(&samples, pos);
            pos += rebin_length as usize;
        }
    }

    if state.bin_count() == 0 {
        return Ok(ResultObservable {
            internal_bin_length: obs_type.internal_bin_length,
            rebin_length,
            mean: ndarray::ArrayD::zeros(obs_type.shape.clone()),
            error: ndarray::ArrayD::zeros(obs_type.shape.clone()),
            covariance: None,
            autocorrelation_time: ndarray::ArrayD::zeros(vec![1]),
            rebin_means: ndarray::ArrayD::zeros(obs_type.shape.clone()),
        });
    }

    let mu = state.mean();
    let sigma = state.std_of_mean();
    let sigma_no_rebin = state.std_of_mean_no_rebin();

    // Compute autocorrelation time
    let autocorrelation_time = compute_autocorrelation_time(&state, &mu, &sigma, &sigma_no_rebin);

    // Compute covariance if requested and observable is multi-component
    let covariance = if options.estimate_covariance
        && obs_type.shape.len() > 0
        && obs_type.shape.iter().product::<usize>() > 1
    {
        let rebin_bins = state.rebin_bins_array();
        if rebin_bins.is_some() {
            let bins_arr = rebin_bins.unwrap();
            let cov = cov_of_mean(&bins_arr);
            Some(cov)
        } else {
            None
        }
    } else {
        None
    };

    // Compute decorrelated autocorrelation time if covariance was computed
    let autocorrelation_time =
        if options.estimate_covariance && covariance.is_some() && obs_type.shape.len() > 0 {
            let rebin_bins = state.rebin_bins_array().unwrap();
            let cov = covariance.as_ref().unwrap();
            compute_decorrelated_autocorr_time(&rebin_bins, &mu, cov, state.bin_count())
        } else {
            autocorrelation_time
        };

    Ok(ResultObservable {
        internal_bin_length: obs_type.internal_bin_length,
        rebin_length,
        mean: mu,
        error: sigma,
        covariance,
        autocorrelation_time,
        rebin_means: state.rebin_means(),
    })
}

/// Compute autocorrelation time from accumulated statistics.
/// For scalar observables: uses variance ratio.
/// For array observables: computes per-element autocorrelation time.
#[allow(dead_code)]
fn compute_autocorrelation_time(
    _state: &AddSamplesState<f64>,
    mu: &ArrayD<f64>,
    sigma: &ArrayD<f64>,
    sigma_no_rebin: &ArrayD<f64>,
) -> ArrayD<f64> {
    // For scalar observables (shape == [] or [1]), return scalar wrapped in 1D array
    let shape = mu.shape();
    if shape.is_empty() || shape == [1] {
        let sigma_val = if shape.is_empty() {
            sigma[ndarray::IxDyn(&[])]
        } else {
            sigma[ndarray::IxDyn(&[0])]
        };
        let sigma_no_rebin_val = if sigma_no_rebin.shape().is_empty() {
            sigma_no_rebin[ndarray::IxDyn(&[])]
        } else {
            sigma_no_rebin[ndarray::IxDyn(&[0])]
        };
        let tau = compute_regular_autocorr_time(0.0, sigma_val, sigma_no_rebin_val);
        return ndarray::ArrayD::from_elem(ndarray::IxDyn(&[1]), tau);
    }

    // For array observables: compute per-element autocorrelation time
    sigma
        .iter()
        .zip(sigma_no_rebin.iter())
        .map(|(s, s0)| compute_regular_autocorr_time(0.0, *s, *s0))
        .collect::<ndarray::Array1<f64>>()
        .into_dyn()
}

/// Internal state for accumulating samples during merge.
///
/// Tracks both the raw sum and sum-of-squares to compute
/// both binned and unbinned statistics for autocorrelation analysis.
pub struct AddSamplesState<T> {
    #[allow(dead_code)]
    bin_length: u64,
    rebin_length: u64,
    shape: Vec<usize>,
    estimate_covariance: bool,
    /// Sum of bin means (for final mean computation)
    sum: ArrayD<T>,
    /// Sum of squared bin means (for unbinned variance)
    sum_sq: ArrayD<T>,
    /// Individual bin means stored for jackknife and covariance
    bin_means: Vec<ArrayD<T>>,
    /// Number of complete bins accumulated
    n_bins: usize,
}

impl AddSamplesState<f64> {
    pub fn new(
        bin_length: u64,
        shape: &[usize],
        rebin_length: u64,
        estimate_covariance: bool,
    ) -> Self {
        let shape_vec = shape.to_vec();
        Self {
            bin_length,
            rebin_length,
            shape: shape_vec.clone(),
            estimate_covariance,
            sum: ndarray::ArrayD::zeros(ndarray::IxDyn(&shape_vec)),
            sum_sq: ndarray::ArrayD::zeros(ndarray::IxDyn(&shape_vec)),
            bin_means: Vec::new(),
            n_bins: 0,
        }
    }

    /// Add one rebin bin's worth of samples.
    ///
    /// `samples` is the full sample array from the HDF5 file (last dim = samples).
    /// `offset` is the starting index.
    pub fn add_rebin_bin(&mut self, samples: &ArrayD<f64>, offset: usize) {
        let n = self.rebin_length as usize;
        let ndim = samples.ndim();
        // For scalar observables stored as 1D array [n_samples], samples.ndim() == 1
        // For array observables stored as [*shape, n_samples], samples.ndim() > 1
        // Handle 0-dim arrays (single scalar sample per bin)
        if ndim == 0 {
            let end = (offset + n).min(1);
            if offset < end {
                let bin_mean = samples.clone();
                self.sum += &bin_mean;
                self.sum_sq += &bin_mean.mapv(|v| v * v);
                if self.estimate_covariance {
                    self.bin_means.push(bin_mean);
                }
                self.n_bins += 1;
            }
            return;
        }

        let n_samples = samples.shape()[ndim - 1];
        let last_axis = ndarray::Axis(ndim - 1);
        let end = (offset + n).min(n_samples);

        // Compute mean of samples[offset..offset+n]
        let shape_dyn = ndarray::IxDyn(&self.shape);
        let mut bin_mean = ndarray::ArrayD::zeros(shape_dyn);

        for i in offset..end {
            let sample = samples.index_axis(last_axis, i);
            bin_mean += &sample;
        }

        let actual_n = (end - offset) as f64;
        if actual_n > 0.0 {
            bin_mean /= actual_n;

            // Accumulate sum and sum of squares
            self.sum += &bin_mean;
            self.sum_sq += &bin_mean.mapv(|v| v * v);

            if self.estimate_covariance {
                self.bin_means.push(bin_mean);
            }
            self.n_bins += 1;
        }
    }

    /// Compute mean from accumulated bins.
    pub fn mean(&self) -> ArrayD<f64> {
        if self.n_bins == 0 {
            return ndarray::ArrayD::zeros(ndarray::IxDyn(&self.shape));
        }
        self.sum.clone() / (self.n_bins as f64)
    }

    /// Compute std of mean from accumulated bins (binned estimate).
    pub fn std_of_mean(&self) -> ArrayD<f64> {
        if self.n_bins < 2 {
            return ndarray::ArrayD::zeros(ndarray::IxDyn(&self.shape));
        }
        let n = self.n_bins as f64;
        let mean = self.mean();
        let variance = self.sum_sq.clone() / n - mean.mapv(|v| v * v);
        variance.mapv(|v| (v.max(0.0)).sqrt() / (n - 1.0).sqrt())
    }

    /// Compute unbinned std of mean (no rebinning applied).
    ///
    /// Since we only store binned means, we approximate the unbinned variance
    /// by scaling the binned variance by rebin_length. This matches the
    /// Carlo.jl approach of comparing binned vs unbinned variance.
    pub fn std_of_mean_no_rebin(&self) -> ArrayD<f64> {
        if self.n_bins < 2 {
            return ndarray::ArrayD::zeros(ndarray::IxDyn(&self.shape));
        }
        let n = self.n_bins as f64;
        let mean = self.mean();
        let variance = self.sum_sq.clone() / n - mean.mapv(|v| v * v);
        let rebin_factor = self.rebin_length as f64;
        variance.mapv(|v| {
            let adjusted = v * rebin_factor;
            (adjusted.max(0.0)).sqrt() / (n * rebin_factor - 1.0).sqrt()
        })
    }

    /// Number of complete bins accumulated.
    pub fn bin_count(&self) -> usize {
        self.n_bins
    }

    /// Get rebin means array (shape: [*shape, bin_count]).
    pub fn rebin_means(&self) -> ArrayD<f64> {
        if self.bin_means.is_empty() {
            let mut result = ndarray::ArrayD::zeros(ndarray::IxDyn(&self.shape));
            // Store the mean if we have data
            if self.n_bins > 0 {
                result = self.mean();
            }
            return result;
        }

        // Stack all bin means into a single array with shape [*shape, n_bins]
        let n_bins = self.bin_means.len();
        let obs_dim: usize = self.shape.iter().product();
        let mut data = vec![0.0f64; obs_dim * n_bins];

        for (k, bin) in self.bin_means.iter().enumerate() {
            let slice = &mut data[k * obs_dim..(k + 1) * obs_dim];
            for (i, v) in bin.iter().enumerate() {
                slice[i] = *v;
            }
        }

        let mut new_shape = self.shape.clone();
        new_shape.push(n_bins);
        ndarray::ArrayD::from_shape_vec(ndarray::IxDyn(&new_shape), data).unwrap()
    }

    /// Get rebin bins as ArrayD for covariance computation.
    /// Returns None if covariance tracking is not enabled.
    pub fn rebin_bins_array(&self) -> Option<ArrayD<f64>> {
        if self.bin_means.is_empty() {
            return None;
        }

        let n_bins = self.bin_means.len();
        let obs_dim: usize = self.shape.iter().product();
        let mut data = vec![0.0f64; obs_dim * n_bins];

        for (k, bin) in self.bin_means.iter().enumerate() {
            let slice = &mut data[k * obs_dim..(k + 1) * obs_dim];
            for (i, v) in bin.iter().enumerate() {
                slice[i] = *v;
            }
        }

        let mut new_shape = self.shape.clone();
        new_shape.push(n_bins);
        Some(ndarray::ArrayD::from_shape_vec(ndarray::IxDyn(&new_shape), data).unwrap())
    }
}

/// Merge results from a task directory.
#[cfg(feature = "hdf5")]
pub fn merge_results(
    taskdir: &PathBuf,
    options: &MergeOptions,
) -> Result<HashMap<String, ResultObservable<f64>>, crate::CarloError> {
    let files = list_meas_files(taskdir)?;
    merge_results_from_files(&files, options)
}

/// Merge raw observables, then let the model register derived observables
/// (evaluables) via `register_fn`, mirroring Carlo.jl's
/// `merge_results(::Type{MC}, taskdir; parameters...)`.
///
/// The callback receives an [`Evaluator`](crate::evaluable::Evaluator)
/// pre-loaded with the merged observables. Use
/// [`Evaluator::evaluate()`](crate::evaluable::Evaluator::evaluate) to
/// register derived quantities (ratios, variances, etc.).
///
/// Evaluable means are appended to the result map. For full error bars
/// on evaluables, prefer using the `Evaluator` directly.
#[cfg(feature = "hdf5")]
pub fn merge_task_results<F>(
    taskdir: &PathBuf,
    options: &MergeOptions,
    register_fn: F,
) -> Result<HashMap<String, ResultObservable<f64>>, crate::CarloError>
where
    F: FnOnce(&mut crate::evaluable::Evaluator),
{
    let results = merge_results(taskdir, options)?;

    let estimate_covariance = options.estimate_covariance;
    let mut evaluator = crate::evaluable::Evaluator::new(results.clone(), estimate_covariance);

    register_fn(&mut evaluator);

    // Merge evaluable means into result map
    let mut results = results;
    for (name, mean_array) in evaluator.evaluables() {
        let obs = ResultObservable {
            internal_bin_length: 0, // evaluables don't have internal bins
            rebin_length: 0,
            mean: mean_array.clone(),
            error: ArrayD::zeros(mean_array.shape()),
            covariance: None,
            autocorrelation_time: ArrayD::zeros(mean_array.shape()),
            rebin_means: mean_array.clone(),
        };
        results.insert(name.clone(), obs);
    }

    Ok(results)
}
