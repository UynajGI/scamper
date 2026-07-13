use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::diagnostics::diagnose;
use crate::target::LogDensity;
use crate::{
    EuclideanState, McmcError, MemoryTrace, MultiChainDiagnostics, SamplingPhase, TraceStore,
    TransitionKernel,
};

/// Configuration for independent-chain execution.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct McmcConfig {
    pub chains: usize,
    pub warmup: u64,
    pub samples: u64,
    pub thinning: usize,
    pub base_seed: u64,
    pub parameter_names: Vec<String>,
}

impl Default for McmcConfig {
    fn default() -> Self {
        Self {
            chains: 4,
            warmup: 1_000,
            samples: 2_000,
            thinning: 1,
            base_seed: 42,
            parameter_names: Vec::new(),
        }
    }
}

impl McmcConfig {
    pub fn validate(&self, initial_positions: &[Vec<f64>]) -> Result<usize, McmcError> {
        if self.chains < 2 {
            return Err(McmcError::InvalidConfig(
                "multi-chain diagnostics require at least two chains".to_string(),
            ));
        }
        if self.samples < 4 {
            return Err(McmcError::InvalidConfig(
                "at least four production transitions are required".to_string(),
            ));
        }
        if self.thinning == 0 {
            return Err(McmcError::InvalidConfig(
                "thinning must be at least one".to_string(),
            ));
        }
        let retained = usize::try_from(self.samples)
            .unwrap_or(usize::MAX)
            .div_ceil(self.thinning);
        if retained < 4 {
            return Err(McmcError::InvalidConfig(
                "thinning must retain at least four draws per chain".to_string(),
            ));
        }
        if initial_positions.len() != self.chains {
            return Err(McmcError::DimensionMismatch {
                expected: self.chains,
                actual: initial_positions.len(),
            });
        }
        let dimension = initial_positions.first().map_or(0, Vec::len);
        if dimension == 0
            || initial_positions
                .iter()
                .any(|position| position.len() != dimension)
        {
            return Err(McmcError::InvalidConfig(
                "all initial positions must have the same positive dimension".to_string(),
            ));
        }
        if !self.parameter_names.is_empty() && self.parameter_names.len() != dimension {
            return Err(McmcError::DimensionMismatch {
                expected: dimension,
                actual: self.parameter_names.len(),
            });
        }
        Ok(dimension)
    }
}

/// Output of one independent chain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChainOutput {
    pub chain_id: usize,
    pub trace: MemoryTrace,
    pub final_position: Vec<f64>,
    pub final_log_density: f64,
}

/// Complete multi-chain posterior output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct McmcOutput {
    pub chains: Vec<ChainOutput>,
    pub diagnostics: MultiChainDiagnostics,
}

/// Run independent chains in Rayon and compute cross-chain diagnostics.
pub fn run_multichain<T, K, FT, FK>(
    target_factory: FT,
    kernel_factory: FK,
    initial_positions: Vec<Vec<f64>>,
    config: McmcConfig,
) -> Result<McmcOutput, McmcError>
where
    T: LogDensity<[f64]>,
    K: TransitionKernel<T>,
    FT: Fn(usize) -> T + Sync,
    FK: Fn(usize) -> K + Sync,
{
    let dimension = config.validate(&initial_positions)?;
    let retained = usize::try_from(config.samples)
        .unwrap_or(usize::MAX)
        .div_ceil(config.thinning);
    let outputs = initial_positions
        .into_par_iter()
        .enumerate()
        .map(|(chain_id, initial)| {
            let mut target = target_factory(chain_id);
            let mut state = EuclideanState::initialize(&mut target, initial)?;
            let mut kernel = kernel_factory(chain_id);
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(chain_seed(config.base_seed, chain_id));

            kernel.on_phase_start(&mut target, SamplingPhase::Warmup, &state)?;
            for _ in 0..config.warmup {
                let report =
                    kernel.transition(&mut target, &mut state, &mut rng, SamplingPhase::Warmup)?;
                report.validate()?;
            }
            kernel.on_phase_end(&mut target, SamplingPhase::Warmup, &state)?;
            kernel.on_phase_start(&mut target, SamplingPhase::Sampling, &state)?;

            let mut trace = MemoryTrace::new(dimension, config.thinning)?;
            trace.reserve_draws(retained);
            for _ in 0..config.samples {
                let report = kernel.transition(
                    &mut target,
                    &mut state,
                    &mut rng,
                    SamplingPhase::Sampling,
                )?;
                let _retained = trace.record(chain_id, &state, &report)?;
            }
            kernel.on_phase_end(&mut target, SamplingPhase::Sampling, &state)?;
            Ok(ChainOutput {
                chain_id,
                trace,
                final_position: state.position().clone(),
                final_log_density: state.log_density(),
            })
        })
        .collect::<Result<Vec<_>, McmcError>>()?;

    let traces = outputs
        .iter()
        .map(|output| output.trace.clone())
        .collect::<Vec<_>>();
    let diagnostics = diagnose(&traces, &config.parameter_names)?;
    Ok(McmcOutput {
        chains: outputs,
        diagnostics,
    })
}

pub(crate) fn chain_seed(base_seed: u64, chain_id: usize) -> u64 {
    let mut value = base_seed ^ (chain_id as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}
