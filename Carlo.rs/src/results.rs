//! Simulation result containers and metadata.

use std::collections::{BTreeMap, HashMap};

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use serde_json::Value as Json;

use crate::{ArrayEstimate, CarloError, ComplexEstimate, Estimate};

/// Simulation results container.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Results {
    /// Observable estimates (scalar).
    #[serde(default)]
    estimates: HashMap<String, Estimate>,

    /// Array observable estimates, element-preserving.
    ///
    /// Serialized as `array_observables` and omitted entirely when empty, so
    /// JSON written before this field existed round-trips unchanged.
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    array_estimates: HashMap<String, ArrayEstimate>,

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

    /// Extensible provenance attributes.
    ///
    /// A generic string-keyed map for domain adapters to record run
    /// provenance (model family, process policy, curve coordinates, ...).
    /// Carlo.rs itself never reads these values; it only stores, merges,
    /// and round-trips them. Omitted from JSON when empty so pre-existing
    /// results files stay byte-compatible.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub attributes: BTreeMap<String, Json>,
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
            attributes: BTreeMap::new(),
        }
    }
}

impl Metadata {
    /// Record one provenance attribute.
    ///
    /// Keys are domain-defined; Carlo.rs attaches no semantics to them.
    /// Later writes under the same key overwrite earlier ones.
    pub fn set_attribute(&mut self, key: impl Into<String>, value: impl Into<Json>) {
        self.attributes.insert(key.into(), value.into());
    }

    /// Look up one provenance attribute.
    pub fn attribute(&self, key: &str) -> Option<&Json> {
        self.attributes.get(key)
    }
}

impl Results {
    /// Create an empty result set with default metadata.
    pub fn new() -> Self {
        Self {
            estimates: HashMap::new(),
            array_estimates: HashMap::new(),
            complex_estimates: HashMap::new(),
            metadata: Metadata::default(),
        }
    }

    /// Create from measurements finalize output.
    pub fn from_measurements(measurements: &HashMap<String, Estimate>) -> Self {
        Self {
            estimates: measurements.clone(),
            array_estimates: HashMap::new(),
            complex_estimates: HashMap::new(),
            metadata: Metadata::default(),
        }
    }

    /// Create from split measurements output (scalar plus array observables).
    ///
    /// This is the finalize path that keeps array observables
    /// element-preserving; see
    /// [`Measurements::finalize_split`](crate::Measurements::finalize_split).
    pub fn from_measurements_with_arrays(
        measurements: &HashMap<String, Estimate>,
        arrays: &HashMap<String, ArrayEstimate>,
    ) -> Self {
        Self {
            estimates: measurements.clone(),
            array_estimates: arrays.clone(),
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
            array_estimates: HashMap::new(),
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

    /// Add an array observable estimate (element-preserving).
    pub fn add_array(&mut self, name: &str, estimate: ArrayEstimate) {
        self.array_estimates.insert(name.to_string(), estimate);
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

    /// Get an array observable estimate.
    pub fn get_array(&self, name: &str) -> Option<&ArrayEstimate> {
        self.array_estimates.get(name)
    }

    /// Get a complex observable estimate.
    pub fn get_complex(&self, name: &str) -> Option<&ComplexResult> {
        self.complex_estimates.get(name)
    }

    /// Merge multiple results into one.
    /// Uses weighted averaging based on inverse variance weighting.
    ///
    /// Array observables merge element-by-element with the same weighting;
    /// results disagreeing on an array observable's shape skip that
    /// observable with a warning instead of being silently flattened.
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

            merged_estimates.insert(name, merge_weighted(&estimates));
        }

        // Merge array estimates element-wise, preserving shape
        let merged_arrays = merge_array_estimates(results);

        // Merge metadata
        let total_measurement_sweeps = results.iter().map(|r| r.metadata.measurement_sweeps).sum();

        let metadata = Metadata {
            version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp: Utc::now(),
            base_seed: results[0].metadata.base_seed,
            thermalization_sweeps: results[0].metadata.thermalization_sweeps,
            measurement_sweeps: total_measurement_sweeps,
            n_tasks: results.len(),
            attributes: merge_metadata_attributes(results),
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

                Some((
                    name,
                    ComplexResult {
                        re: merge_weighted(&re_estimates),
                        im: merge_weighted(&im_estimates),
                    },
                ))
            })
            .collect();

        Self {
            estimates: merged_estimates,
            array_estimates: merged_arrays,
            complex_estimates: merged_complex,
            metadata,
        }
    }

    /// Export to JSON string.
    /// Complex observables are serialized in {re, im} format matching Carlo.jl.
    ///
    /// Array observables appear under `array_observables` with their shape
    /// and one element estimate per element; the key is omitted when there
    /// are none, so scalar-only output is byte-identical to the legacy
    /// format.
    pub fn to_json(&self) -> Result<String, CarloError> {
        #[derive(Serialize)]
        struct JsonResults<'a> {
            observables: &'a HashMap<String, Estimate>,
            #[serde(skip_serializing_if = "HashMap::is_empty")]
            array_observables: &'a HashMap<String, ArrayEstimate>,
            #[serde(skip_serializing_if = "HashMap::is_empty")]
            complex_observables: &'a HashMap<String, ComplexResult>,
            metadata: &'a Metadata,
        }

        let json_data = JsonResults {
            observables: &self.estimates,
            array_observables: &self.array_estimates,
            complex_observables: &self.complex_estimates,
            metadata: &self.metadata,
        };

        serde_json::to_string_pretty(&json_data).map_err(|e| CarloError::InvalidConfig {
            field: "json_output".into(),
            reason: e.to_string(),
        })
    }

    /// Parse results from a JSON string produced by [`Results::to_json`].
    ///
    /// Accepts legacy files without `array_observables` and metadata without
    /// `attributes`: both default to empty.
    pub fn from_json(json: &str) -> Result<Self, CarloError> {
        #[derive(Deserialize)]
        struct JsonResults {
            #[serde(default)]
            observables: HashMap<String, Estimate>,
            #[serde(default)]
            array_observables: HashMap<String, ArrayEstimate>,
            #[serde(default)]
            complex_observables: HashMap<String, ComplexResult>,
            #[serde(default)]
            metadata: Metadata,
        }

        let parsed: JsonResults =
            serde_json::from_str(json).map_err(CarloError::SerializationError)?;
        Ok(Self {
            estimates: parsed.observables,
            array_estimates: parsed.array_observables,
            complex_estimates: parsed.complex_observables,
            metadata: parsed.metadata,
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

    /// Get all array estimates.
    pub fn array_estimates(&self) -> &HashMap<String, ArrayEstimate> {
        &self.array_estimates
    }

    /// Get all complex estimates.
    pub fn complex_estimates(&self) -> &HashMap<String, ComplexResult> {
        &self.complex_estimates
    }
}

/// Merge one observable's estimates across runs.
///
/// Weighted mean by bin count; combined stderr from per-run variances;
/// autocorrelation times averaged; bin counts summed. This is exactly the
/// weighting [`Results::merge`] has always applied to scalars, factored out
/// so array elements merge identically.
fn merge_weighted(estimates: &[&Estimate]) -> Estimate {
    let total_bins: usize = estimates.iter().map(|e| e.n_bins).sum();

    if total_bins == 0 {
        return Estimate {
            mean: 0.0,
            stderr: 0.0,
            autocorr_time: 1.0,
            n_bins: 0,
        };
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

    Estimate {
        mean: weighted_mean,
        stderr: combined_stderr,
        autocorr_time: estimates.iter().map(|e| e.autocorr_time).sum::<f64>()
            / estimates.len() as f64,
        n_bins: total_bins,
    }
}

/// Merge array observables element-wise across runs.
///
/// An observable name is merged only when every result carrying it reports
/// the same shape; otherwise the observable is skipped with a warning (never
/// flattened or truncated). Observables missing from some results merge the
/// ones present, matching the scalar behavior.
fn merge_array_estimates(results: &[Results]) -> HashMap<String, ArrayEstimate> {
    let all_names: std::collections::HashSet<String> = results
        .iter()
        .flat_map(|r| r.array_estimates.keys().cloned())
        .collect();

    let mut merged = HashMap::new();
    for name in all_names {
        let entries: Vec<&ArrayEstimate> = results
            .iter()
            .filter_map(|r| r.array_estimates.get(&name))
            .collect();
        if entries.is_empty() {
            continue;
        }

        let shape = &entries[0].shape;
        if entries.iter().any(|e| e.shape != *shape) {
            tracing::warn!(
                observable = %name,
                "array observable shape differs across results; skipping in merge"
            );
            continue;
        }

        let n_elements = shape.iter().product::<usize>();
        let elements: Vec<Estimate> = (0..n_elements)
            .map(|i| {
                let column: Vec<&Estimate> = entries.iter().map(|e| &e.elements[i]).collect();
                merge_weighted(&column)
            })
            .collect();

        // Element count equals the shape product by construction, so this
        // cannot fail; map the impossible case to a skip defensively.
        if let Ok(estimate) = ArrayEstimate::new(shape.clone(), elements) {
            merged.insert(name, estimate);
        } else {
            tracing::warn!(observable = %name, "merged array estimate failed validation");
        }
    }
    merged
}

/// Merge provenance attributes across runs.
///
/// Attributes are kept when every result that has them agrees; conflicting
/// values are dropped with a warning rather than picking a winner silently.
fn merge_metadata_attributes(results: &[Results]) -> BTreeMap<String, Json> {
    let mut merged: BTreeMap<String, Json> = BTreeMap::new();
    let mut conflicts: Vec<String> = Vec::new();

    for result in results {
        for (key, value) in &result.metadata.attributes {
            match merged.get(key) {
                Some(existing) if existing != value => {
                    if !conflicts.contains(key) {
                        conflicts.push(key.clone());
                    }
                }
                _ => {
                    merged.insert(key.clone(), value.clone());
                }
            }
        }
    }

    for key in conflicts {
        tracing::warn!(attribute = %key, "metadata attribute differs across results; dropping it from the merge");
        merged.remove(&key);
    }
    merged
}

impl Default for Results {
    fn default() -> Self {
        Self::new()
    }
}
