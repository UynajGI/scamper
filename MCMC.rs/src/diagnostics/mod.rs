mod ess;
mod rank;

use serde::{Deserialize, Serialize};

use crate::{McmcError, MemoryTrace, TraceStore};

/// Per-parameter convergence and Monte Carlo accuracy diagnostics.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ParameterDiagnostics {
    pub name: String,
    pub mean: f64,
    pub std_dev: f64,
    pub mcse: f64,
    pub rhat: f64,
    pub ess_bulk: f64,
    pub ess_tail: f64,
    pub median: f64,
    pub q05: f64,
    pub q95: f64,
}

/// Diagnostics aggregated across independent chains.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct MultiChainDiagnostics {
    pub parameters: Vec<ParameterDiagnostics>,
    pub total_divergences: u64,
    pub mean_acceptance: f64,
    pub chains: usize,
    pub draws_per_chain: usize,
    /// Per-chain energy Bayesian fraction of missing information.
    ///
    /// `None` means that a chain does not contain a complete, finite energy
    /// series or that the series has zero variance.
    #[serde(default)]
    pub chain_ebfmi: Vec<Option<f64>>,
    /// Number of retained transitions that exhausted the configured NUTS
    /// tree-depth limit.
    #[serde(default)]
    pub max_tree_depth_hits: u64,
}

/// Compute the energy Bayesian fraction of missing information (E-BFMI).
///
/// The estimate is the mean squared change in Hamiltonian energy divided by
/// the sample variance numerator of the energy series. Missing or non-finite
/// energies, fewer than three draws, and a constant energy series return
/// `None`.
pub fn energy_bfmi(energies: &[Option<f64>]) -> Option<f64> {
    if energies.len() < 3 {
        return None;
    }
    let values = energies.iter().copied().collect::<Option<Vec<_>>>()?;
    if values.iter().any(|value| !value.is_finite()) {
        return None;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let denominator = values
        .iter()
        .map(|value| {
            let centered = value - mean;
            centered * centered
        })
        .sum::<f64>();
    if !denominator.is_finite() || denominator <= 0.0 {
        return None;
    }
    let numerator = values
        .windows(2)
        .map(|pair| {
            let difference = pair[1] - pair[0];
            difference * difference
        })
        .sum::<f64>();
    let estimate = numerator / denominator;
    estimate.is_finite().then_some(estimate)
}

/// Compute rank-normalized split R-hat, bulk/tail ESS and MCSE.
pub fn diagnose(
    traces: &[MemoryTrace],
    parameter_names: &[String],
) -> Result<MultiChainDiagnostics, McmcError> {
    if traces.len() < 2 || traces.iter().any(|trace| trace.len() < 4) {
        return Err(McmcError::InsufficientDraws);
    }
    let dimension = traces[0].dimension();
    if traces.iter().any(|trace| trace.dimension() != dimension) {
        return Err(McmcError::InconsistentTraceDimension);
    }
    if !parameter_names.is_empty() && parameter_names.len() != dimension {
        return Err(McmcError::DimensionMismatch {
            expected: dimension,
            actual: parameter_names.len(),
        });
    }
    let draws_per_chain = traces.iter().map(TraceStore::len).min().unwrap_or(0);
    let mut parameters = Vec::with_capacity(dimension);
    for parameter in 0..dimension {
        let chains = traces
            .iter()
            .map(|trace| {
                trace
                    .parameter(parameter)
                    .expect("validated parameter index")
                    .take(draws_per_chain)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let pooled = chains.iter().flatten().copied().collect::<Vec<_>>();
        let mean = pooled.iter().sum::<f64>() / pooled.len() as f64;
        let variance = pooled
            .iter()
            .map(|value| {
                let delta = value - mean;
                delta * delta
            })
            .sum::<f64>()
            / (pooled.len() - 1) as f64;
        let std_dev = variance.sqrt();
        let mut sorted = pooled.clone();
        sorted.sort_by(f64::total_cmp);
        let q05 = rank::quantile_sorted(&sorted, 0.05);
        let median = rank::quantile_sorted(&sorted, 0.5);
        let q95 = rank::quantile_sorted(&sorted, 0.95);

        let split = rank::split_chains(&chains);
        let rank_normalized = rank::rank_normalize(&split);
        let folded = split
            .iter()
            .map(|chain| {
                chain
                    .iter()
                    .map(|value| (value - median).abs())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let folded_rank_normalized = rank::rank_normalize(&folded);
        let rhat = rank::rhat(&rank_normalized).max(rank::rhat(&folded_rank_normalized));
        let ess_bulk = ess::effective_sample_size(&rank_normalized);
        let lower_indicators = split
            .iter()
            .map(|chain| {
                chain
                    .iter()
                    .map(|value| if *value <= q05 { 1.0 } else { 0.0 })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let upper_indicators = split
            .iter()
            .map(|chain| {
                chain
                    .iter()
                    .map(|value| if *value >= q95 { 1.0 } else { 0.0 })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let ess_tail = ess::effective_sample_size(&lower_indicators)
            .min(ess::effective_sample_size(&upper_indicators));
        let mcse = if ess_bulk > 0.0 {
            std_dev / ess_bulk.sqrt()
        } else {
            f64::NAN
        };
        let name = parameter_names
            .get(parameter)
            .cloned()
            .unwrap_or_else(|| format!("x[{parameter}]"));
        parameters.push(ParameterDiagnostics {
            name,
            mean,
            std_dev,
            mcse,
            rhat,
            ess_bulk,
            ess_tail,
            median,
            q05,
            q95,
        });
    }

    let total_divergences = traces
        .iter()
        .flat_map(|trace| trace.divergences().iter().copied())
        .map(u64::from)
        .sum();
    let finite_acceptance = traces
        .iter()
        .flat_map(|trace| trace.acceptance_rates().iter().copied().flatten())
        .filter(|value| value.is_finite())
        .collect::<Vec<_>>();
    let mean_acceptance = if finite_acceptance.is_empty() {
        f64::NAN
    } else {
        finite_acceptance.iter().sum::<f64>() / finite_acceptance.len() as f64
    };

    let chain_ebfmi = traces
        .iter()
        .map(|trace| energy_bfmi(trace.energies()))
        .collect();
    let max_tree_depth_hits = traces
        .iter()
        .flat_map(|trace| trace.max_tree_depth_reached().iter().copied())
        .map(u64::from)
        .sum();

    Ok(MultiChainDiagnostics {
        parameters,
        total_divergences,
        mean_acceptance,
        chains: traces.len(),
        draws_per_chain,
        chain_ebfmi,
        max_tree_depth_hits,
    })
}
