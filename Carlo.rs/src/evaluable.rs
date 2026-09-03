//! Jackknife resampling for error propagation.
//!
//! Jackknife resampling allows computing errors on derived observables
//! (functions of primary observables) while properly propagating correlations.
//!
//! # Usage
//!
//! Use [`jackknife()`] to compute a derived observable with error propagation.
//!
//! # Theory
//!
//! Jackknife creates n resampled datasets by omitting one bin at a time,
//! computes the derived observable on each, and estimates error from
//! the variance of these jackknife estimates.

use crate::merge::ResultObservable;
use crate::CarloError;
use ndarray::ArrayD;
use std::collections::HashMap;

/// Result type for jackknife resampling: (mean, error, covariance)
pub type JackknifeResult = (ArrayD<f64>, ArrayD<f64>, Option<ArrayD<f64>>);

/// Perform jackknife resampling.
#[allow(clippy::type_complexity)]
pub fn jackknife<F, T>(
    func: F,
    sample_sets: &[ArrayD<T>],
    _estimate_covariance: bool,
) -> Result<JackknifeResult, CarloError>
where
    F: Fn(&[ArrayD<T>]) -> ArrayD<T>,
    T: ndarray::LinalgScalar
        + std::fmt::Debug
        + num_traits::ToPrimitive
        + num_traits::FromPrimitive
        + ndarray::ScalarOperand,
{
    let sample_count = sample_sets.iter().map(|s| s.len()).min().unwrap_or(0);
    if sample_count == 0 {
        return Err(CarloError::InvalidConfig {
            field: "samples".into(),
            reason: "Empty sample set".into(),
        });
    }

    // Compute means (sum / n for each sample set)
    let n = T::from_usize(sample_count).unwrap_or(T::one());
    let means: Vec<ArrayD<T>> = sample_sets
        .iter()
        .map(|s| {
            let sum = s.sum();
            ArrayD::from_elem(s.shape(), sum / n)
        })
        .collect();
    let complete_eval = func(&means);

    // Convert complete_eval to f64 for output
    let complete_f64 = complete_eval.mapv(|v| v.to_f64().unwrap_or(0.0));

    // Compute jackknife evaluations
    let jacked_evals: Vec<ArrayD<T>> = (0..sample_count)
        .map(|k| {
            let jacked_means: Vec<ArrayD<T>> = sample_sets
                .iter()
                .map(|s| {
                    let sum = s.sum();
                    let val_k = s.iter().nth(k).copied().unwrap_or(T::zero());
                    let n_minus_1 = T::from_usize(sample_count - 1).unwrap_or(T::one());
                    ArrayD::from_elem(s.shape(), (sum - val_k) / n_minus_1)
                })
                .collect();
            func(&jacked_means)
        })
        .collect();

    // Jackknife mean
    let zero_arr = ArrayD::zeros(complete_eval.raw_dim());
    let jacked_mean: ArrayD<T> = jacked_evals.iter().fold(zero_arr, |acc, e| acc + e.clone()) / n;
    let jacked_mean_f64 = jacked_mean.mapv(|v| v.to_f64().unwrap_or(0.0));

    // Bias-corrected mean: n * complete - (n-1) * jacked_mean
    let n_f64 = sample_count as f64;
    let n_minus_1_f64 = (sample_count - 1) as f64;
    let bias_corrected_mean = complete_f64 * n_f64 - jacked_mean_f64.clone() * n_minus_1_f64;

    // Error: sqrt((n-1)/n * Σ(jacked_eval - jacked_mean)^2)
    let zero_f64 = ArrayD::zeros(complete_eval.raw_dim());
    let sum_sq_diff = jacked_evals.iter().fold(zero_f64, |acc, e| {
        let e_f64 = e.mapv(|v| v.to_f64().unwrap_or(0.0));
        let diff = &e_f64 - &jacked_mean_f64;
        acc + diff.mapv(|v| v * v)
    });

    let scale = n_minus_1_f64 / n_f64;
    let error = sum_sq_diff.mapv(|v: f64| (scale * v).sqrt());

    Ok((bias_corrected_mean, error, None))
}

/// Evaluator for defining derived observables.
///
/// Uses `f64` as the scalar type for computed values.
pub struct Evaluator {
    observables: HashMap<String, ResultObservable<f64>>,
    evaluables: HashMap<String, ArrayD<f64>>,
    #[allow(dead_code)]
    estimate_covariance: bool,
}

impl Evaluator {
    pub fn new(
        observables: HashMap<String, ResultObservable<f64>>,
        estimate_covariance: bool,
    ) -> Self {
        Self {
            observables,
            evaluables: HashMap::new(),
            estimate_covariance,
        }
    }

    pub fn observables(&self) -> &HashMap<String, ResultObservable<f64>> {
        &self.observables
    }

    pub fn evaluables(&self) -> &HashMap<String, ArrayD<f64>> {
        &self.evaluables
    }

    /// Define an evaluable (derived observable) using jackknife resampling.
    ///
    /// # Arguments
    /// - `name`: Name of the evaluable
    /// - `ingredients`: Names of the observables to use
    /// - `func`: Function that computes the derived observable from ingredient means
    ///
    /// # Example
    /// ```text
    /// evaluator.evaluate("susceptibility", &["magnetization", "magnetization_sq"], |args| {
    ///     &args[1] - &(&args[0] * &args[0])  // χ = <M²> - <M>²
    /// });
    /// ```
    pub fn evaluate<F>(
        &mut self,
        name: &str,
        ingredients: &[&str],
        func: F,
    ) -> Result<(), CarloError>
    where
        F: Fn(&[ArrayD<f64>]) -> ArrayD<f64>,
    {
        // Check that all ingredients exist
        let missing: Vec<&&str> = ingredients
            .iter()
            .filter(|i| !self.observables.contains_key(**i))
            .collect();

        if !missing.is_empty() {
            tracing::warn!(
                "Evaluable '{}': ingredients {:?} not found. Skipping...",
                name,
                missing
            );
            return Ok(());
        }

        // Collect rebin_means from ingredients
        let sample_sets: Vec<ArrayD<f64>> = ingredients
            .iter()
            .map(|i| {
                self.observables
                    .get(*i)
                    .map(|obs| obs.rebin_means.clone())
                    .unwrap_or_else(|| ArrayD::zeros(vec![0]))
            })
            .collect();

        // Check that we have samples
        let bin_count = sample_sets.first().map(|s| s.len()).unwrap_or(0);
        if bin_count == 0 {
            return Ok(());
        }

        // Run jackknife
        let (mean, _error, _cov) = jackknife(func, &sample_sets, self.estimate_covariance)?;

        // Store result
        self.evaluables.insert(name.to_string(), mean);

        Ok(())
    }
}

/// MultiplexEvaluator for parallel tempering chains.
///
/// Allows registering multiple evaluation functions (one per PT chain)
/// and stacking the results together.
///
/// # Usage
///
/// ```text
/// let mut multi_eval = MultiplexEvaluator::new(4); // 4 PT chains
///
/// // Register evaluator for each chain
/// for chain_idx in 0..4 {
///     multi_eval.evaluate("order_parameter", &["magnetization"], move |args| {
///         // Chain-specific computation
///         args[0] * temperature_factors[chain_idx]
///     });
/// }
///
/// // Run evaluations and stack results
/// multi_eval.run_evaluations(&mut evaluator);
/// ```
pub struct MultiplexEvaluator {
    /// Number of PT chains.
    entry_count: usize,
    /// Registered evaluables: name -> (ingredients, functions).
    evals: HashMap<String, EvaluableEntry>,
}

/// An evaluable entry: (ingredient names, evaluation functions per PT chain).
type EvaluableEntry = (
    Vec<String>,
    Vec<Box<dyn Fn(&[ArrayD<f64>]) -> ArrayD<f64> + Send + Sync>>,
);

impl MultiplexEvaluator {
    /// Create a new MultiplexEvaluator for the given number of chains.
    pub fn new(entry_count: usize) -> Self {
        Self {
            entry_count,
            evals: HashMap::new(),
        }
    }

    /// Register an evaluation function for a chain.
    ///
    /// The function will be called with the ingredient arrays and should return
    /// the computed value for that chain.
    pub fn evaluate<F>(&mut self, name: &str, ingredients: &[&str], func: F)
    where
        F: Fn(&[ArrayD<f64>]) -> ArrayD<f64> + Send + Sync + 'static,
    {
        let entry = self.evals.entry(name.to_string()).or_insert_with(|| {
            (
                ingredients.iter().map(|s| s.to_string()).collect(),
                Vec::new(),
            )
        });

        // Check ingredient consistency
        let existing_ingredients: Vec<&str> = entry.0.iter().map(|s| s.as_str()).collect();
        if existing_ingredients != ingredients {
            tracing::error!(
                "Evaluable '{}' has inconsistent ingredients: {:?} != {:?}",
                name,
                existing_ingredients,
                ingredients
            );
            return;
        }

        entry.1.push(Box::new(func));
    }

    /// Run all registered evaluations and populate the evaluator.
    ///
    /// Results from all chains are stacked along a new first dimension.
    pub fn run_evaluations(&self, evaluator: &mut Evaluator) -> Result<(), CarloError> {
        for (name, (ingredients, funcs)) in &self.evals {
            if funcs.len() != self.entry_count {
                return Err(CarloError::InvalidConfig {
                    field: "multiplex_evaluator".into(),
                    reason: format!(
                        "Evaluable '{}': expected {} functions, got {}",
                        name,
                        self.entry_count,
                        funcs.len()
                    ),
                });
            }

            // Collect sample sets from ingredients
            let sample_sets: Vec<ArrayD<f64>> = ingredients
                .iter()
                .filter_map(|i| evaluator.observables.get(i))
                .map(|obs| obs.rebin_means.clone())
                .collect();

            if sample_sets.len() != ingredients.len() {
                tracing::warn!(
                    "Evaluable '{}': some ingredients not found. Skipping...",
                    name
                );
                continue;
            }

            let bin_count = sample_sets.first().map(|s| s.len()).unwrap_or(0);
            if bin_count == 0 {
                continue;
            }

            // Create stacked evaluation function
            let stacked_func = move |args: &[ArrayD<f64>]| -> ArrayD<f64> {
                use ndarray::{Array, IxDyn};

                // Compute result for each chain
                let results: Vec<ArrayD<f64>> = funcs.iter().map(|f| f(args)).collect();

                // Stack results along new first dimension
                if results.is_empty() {
                    return ArrayD::zeros(vec![0]);
                }

                // Get the shape of the result
                let result_shape = results[0].shape();
                let entry_count = results.len();

                // Build stacked array with shape (entry_count, ...result_shape)
                let mut stacked_shape = vec![entry_count];
                stacked_shape.extend(result_shape.iter().copied());

                let total_size: usize = stacked_shape.iter().product();
                let mut stacked_data = Vec::with_capacity(total_size);

                for result in &results {
                    stacked_data.extend(result.iter().copied());
                }

                Array::from_shape_vec(IxDyn(&stacked_shape), stacked_data)
                    .unwrap_or_else(|_| ArrayD::zeros(vec![0]))
            };

            // Run jackknife with stacked function
            let (mean, _error, _cov) =
                jackknife(stacked_func, &sample_sets, evaluator.estimate_covariance)?;

            evaluator.evaluables.insert(name.clone(), mean);
        }

        Ok(())
    }

    /// Get number of chains.
    pub fn entry_count(&self) -> usize {
        self.entry_count
    }

    /// Get number of registered evaluables.
    pub fn len(&self) -> usize {
        self.evals.len()
    }

    /// Check if no evaluables are registered.
    pub fn is_empty(&self) -> bool {
        self.evals.is_empty()
    }
}
