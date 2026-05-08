use std::collections::HashMap;

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};

use crate::{CarloError, ComplexEstimate, Estimate};

/// Simulation results container.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Results {
    /// Observable estimates (scalar and array).
    #[serde(default)]
    estimates: HashMap<String, Estimate>,

    /// Complex observable estimates (stored as {re, im}).
    #[serde(default)]
    complex_estimates: HashMap<String, ComplexResult>,

    /// Run metadata.
    #[serde(default)]
    metadata: Metadata,
}

/// Complex result matching Carlo.jl JSON format: `{re, im}`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComplexResult {
    /// Real part with re/im breakdown (not just mean/error).
    pub re: Estimate,
    /// Imaginary part with re/im breakdown.
    pub im: Estimate,
}

impl ComplexResult {
    /// Create from complex estimate.
    pub fn from_estimate(est: &ComplexEstimate) -> Self {
        Self {
            re: est.re.clone(),
            im: est.im.clone(),
        }
    }

    /// Reconstruct as complex estimate.
    pub fn to_estimate(&self) -> ComplexEstimate {
        ComplexEstimate::new(self.re.clone(), self.im.clone())
    }
}

/// Run metadata for reproducibility.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Metadata {
    /// Carlo version.
    pub version: String,

    /// Timestamp.
    pub timestamp: DateTime<Utc>,

    /// Base seed used.
    pub base_seed: u64,

    /// Thermalization sweeps.
    pub thermalization_sweeps: u64,

    /// Measurement sweeps.
    pub measurement_sweeps: u64,

    /// Number of parallel tasks.
    pub n_tasks: usize,
}

impl Default for Metadata {
    fn default() -> Self {
        Self {
            version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp: Utc::now(),
            base_seed: 0,
            thermalization_sweeps: 0,
            measurement_sweeps: 0,
            n_tasks: 1,
        }
    }
}

impl Results {
    pub fn new() -> Self {
        Self {
            estimates: HashMap::new(),
            complex_estimates: HashMap::new(),
            metadata: Metadata::default(),
        }
    }

    /// Create from measurements finalize output.
    pub fn from_measurements(measurements: &HashMap<String, Estimate>) -> Self {
        Self {
            estimates: measurements.clone(),
            complex_estimates: HashMap::new(),
            metadata: Metadata::default(),
        }
    }

    /// Create from measurements including complex estimates.
    pub fn from_measurements_with_complex(
        measurements: &HashMap<String, Estimate>,
        complex: &HashMap<String, ComplexEstimate>,
    ) -> Self {
        Self {
            estimates: measurements.clone(),
            complex_estimates: complex
                .iter()
                .map(|(k, v)| (k.clone(), ComplexResult::from_estimate(v)))
                .collect(),
            metadata: Metadata::default(),
        }
    }

    /// Add an observable estimate.
    pub fn add(&mut self, name: &str, estimate: Estimate) {
        self.estimates.insert(name.to_string(), estimate);
    }

    /// Add a complex observable estimate.
    pub fn add_complex(&mut self, name: &str, estimate: ComplexEstimate) {
        self.complex_estimates
            .insert(name.to_string(), ComplexResult::from_estimate(&estimate));
    }

    /// Get an observable estimate.
    pub fn get(&self, name: &str) -> Option<&Estimate> {
        self.estimates.get(name)
    }

    /// Get a complex observable estimate.
    pub fn get_complex(&self, name: &str) -> Option<&ComplexResult> {
        self.complex_estimates.get(name)
    }

    /// Merge multiple results into one.
    /// Uses weighted averaging based on inverse variance weighting.
    pub fn merge(results: &[Results]) -> Self {
        if results.is_empty() {
            return Self::new();
        }

        if results.len() == 1 {
            return results[0].clone();
        }

        // Merge estimates using weighted average
        let mut merged_estimates: HashMap<String, Estimate> = HashMap::new();

        // Collect all observable names
        let all_names: std::collections::HashSet<String> = results
            .iter()
            .flat_map(|r| r.estimates.keys().cloned())
            .collect();

        for name in all_names {
            // Collect estimates for this observable
            let estimates: Vec<&Estimate> = results
                .iter()
                .filter_map(|r| r.estimates.get(&name))
                .collect();

            if estimates.is_empty() {
                continue;
            }

            // Weighted average: weight = n_bins (inverse variance proxy)
            let total_bins: usize = estimates.iter().map(|e| e.n_bins).sum();

            if total_bins == 0 {
                merged_estimates.insert(
                    name,
                    Estimate {
                        mean: 0.0,
                        stderr: 0.0,
                        autocorr_time: 1.0,
                        n_bins: 0,
                    },
                );
                continue;
            }

            // Weighted mean: sum(mean * n_bins) / total_bins
            let weighted_mean: f64 = estimates
                .iter()
                .map(|e| e.mean * e.n_bins as f64)
                .sum::<f64>()
                / total_bins as f64;

            // Combined variance estimate
            // For independent runs: variance_combined ≈ sum(var_i) / n_runs²
            // Using stderr² = variance / n_bins
            let combined_variance: f64 = if estimates.len() > 1 {
                estimates
                    .iter()
                    .map(|e| (e.stderr * (e.n_bins as f64).sqrt()).powi(2))
                    .sum::<f64>()
                    / (estimates.len() as f64).powi(2)
            } else {
                estimates[0].stderr.powi(2) * estimates[0].n_bins as f64
            };

            let combined_stderr = combined_variance.sqrt() / (total_bins as f64).sqrt();

            merged_estimates.insert(
                name,
                Estimate {
                    mean: weighted_mean,
                    stderr: combined_stderr,
                    autocorr_time: estimates.iter().map(|e| e.autocorr_time).sum::<f64>()
                        / estimates.len() as f64,
                    n_bins: total_bins,
                },
            );
        }

        // Merge metadata
        let total_measurement_sweeps = results.iter().map(|r| r.metadata.measurement_sweeps).sum();

        let metadata = Metadata {
            version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp: Utc::now(),
            base_seed: results[0].metadata.base_seed,
            thermalization_sweeps: results[0].metadata.thermalization_sweeps,
            measurement_sweeps: total_measurement_sweeps,
            n_tasks: results.len(),
        };

        // Merge complex estimates
        let merged_complex: HashMap<String, ComplexResult> = results
            .iter()
            .flat_map(|r| r.complex_estimates.keys().cloned())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .filter_map(|name| {
                let re_estimates: Vec<&Estimate> = results
                    .iter()
                    .filter_map(|r| r.complex_estimates.get(&name).map(|cr| &cr.re))
                    .collect();
                let im_estimates: Vec<&Estimate> = results
                    .iter()
                    .filter_map(|r| r.complex_estimates.get(&name).map(|cr| &cr.im))
                    .collect();

                if re_estimates.is_empty() || im_estimates.is_empty() {
                    return None;
                }

                let merge_real = |ests: &[&Estimate]| -> Estimate {
                    let total_bins: usize = ests.iter().map(|e| e.n_bins).sum();
                    if total_bins == 0 {
                        return Estimate::from_bins(&[]);
                    }
                    let weighted_mean: f64 =
                        ests.iter().map(|e| e.mean * e.n_bins as f64).sum::<f64>()
                            / total_bins as f64;
                    Estimate {
                        mean: weighted_mean,
                        stderr: ests[0].stderr,
                        autocorr_time: ests.iter().map(|e| e.autocorr_time).sum::<f64>()
                            / ests.len() as f64,
                        n_bins: total_bins,
                    }
                };

                Some((
                    name,
                    ComplexResult {
                        re: merge_real(&re_estimates),
                        im: merge_real(&im_estimates),
                    },
                ))
            })
            .collect();

        Self {
            estimates: merged_estimates,
            complex_estimates: merged_complex,
            metadata,
        }
    }

    /// Export to JSON string.
    /// Complex observables are serialized in {re, im} format matching Carlo.jl.
    pub fn to_json(&self) -> Result<String, CarloError> {
        #[derive(Serialize)]
        struct JsonResults<'a> {
            observables: &'a HashMap<String, Estimate>,
            #[serde(skip_serializing_if = "HashMap::is_empty")]
            complex_observables: &'a HashMap<String, ComplexResult>,
            metadata: &'a Metadata,
        }

        let json_data = JsonResults {
            observables: &self.estimates,
            complex_observables: &self.complex_estimates,
            metadata: &self.metadata,
        };

        serde_json::to_string_pretty(&json_data).map_err(|e| CarloError::InvalidConfig {
            field: "json_output".into(),
            reason: e.to_string(),
        })
    }

    /// Set metadata.
    pub fn set_metadata(&mut self, metadata: Metadata) {
        self.metadata = metadata;
    }

    /// Get metadata.
    pub fn metadata(&self) -> &Metadata {
        &self.metadata
    }

    /// Get all estimates.
    pub fn estimates(&self) -> &HashMap<String, Estimate> {
        &self.estimates
    }

    /// Get all complex estimates.
    pub fn complex_estimates(&self) -> &HashMap<String, ComplexResult> {
        &self.complex_estimates
    }
}

impl Default for Results {
    fn default() -> Self {
        Self::new()
    }
}
