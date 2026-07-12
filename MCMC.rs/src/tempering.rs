use rand::{Rng, RngExt};
use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::multichain::chain_seed;
use crate::target::{validate_log_density, LogDensity};
use crate::{EuclideanState, McmcError, MemoryTrace, SamplingPhase, TraceStore, TransitionKernel};

/// Configuration for local shared-memory replica exchange.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TemperingConfig {
    /// User-defined ladder values passed to `target_factory` and `kernel_factory`.
    pub ladder: Vec<f64>,
    pub warmup: u64,
    pub samples: u64,
    pub thinning: usize,
    /// Number of local transitions between neighboring exchange attempts.
    pub exchange_interval: u64,
    pub base_seed: u64,
    pub parameter_names: Vec<String>,
}

impl Default for TemperingConfig {
    fn default() -> Self {
        Self {
            ladder: vec![1.0, 0.5, 0.25, 0.125],
            warmup: 1_000,
            samples: 2_000,
            thinning: 1,
            exchange_interval: 10,
            base_seed: 42,
            parameter_names: Vec::new(),
        }
    }
}

impl TemperingConfig {
    pub fn validate(&self, initial_positions: &[Vec<f64>]) -> Result<usize, McmcError> {
        if self.ladder.len() < 2 {
            return Err(McmcError::InvalidConfig(
                "parallel tempering requires at least two ladder values".to_string(),
            ));
        }
        if self.ladder.iter().any(|value| !value.is_finite()) {
            return Err(McmcError::InvalidConfig(
                "tempering ladder values must be finite".to_string(),
            ));
        }
        if self.samples == 0 || self.thinning == 0 || self.exchange_interval == 0 {
            return Err(McmcError::InvalidConfig(
                "samples, thinning and exchange interval must be positive".to_string(),
            ));
        }
        if initial_positions.len() != self.ladder.len() {
            return Err(McmcError::DimensionMismatch {
                expected: self.ladder.len(),
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
                "all replica initial positions must have the same positive dimension".to_string(),
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

/// Exchange diagnostics for one neighboring ladder edge.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ExchangeEdgeDiagnostics {
    pub left_index: usize,
    pub right_index: usize,
    pub attempts: u64,
    pub acceptances: u64,
}

impl ExchangeEdgeDiagnostics {
    pub fn acceptance_rate(&self) -> Option<f64> {
        (self.attempts > 0).then(|| self.acceptances as f64 / self.attempts as f64)
    }
}

/// Posterior trace and final state for one fixed ladder slot.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TemperedChainOutput {
    pub slot: usize,
    pub ladder_value: f64,
    pub trace: MemoryTrace,
    pub final_position: Vec<f64>,
    pub final_log_density: f64,
}

/// Output of a local parallel-tempering run.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TemperingOutput {
    pub parameter_names: Vec<String>,
    pub replicas: Vec<TemperedChainOutput>,
    pub exchanges: Vec<ExchangeEdgeDiagnostics>,
}

struct Replica<T, K> {
    slot: usize,
    ladder_value: f64,
    target: T,
    kernel: K,
    state: EuclideanState,
    rng: Xoshiro256PlusPlus,
    trace: MemoryTrace,
}

/// Run fixed-slot replica exchange using Rayon for local transitions.
///
/// Each factory receives `(slot, ladder_value)`. The exchange ratio is computed
/// generically by evaluating each neighboring state under the other slot's
/// target, so callers may temper the likelihood only, the full posterior, or a
/// different model parameter without changing this runtime.
pub fn run_parallel_tempering<T, K, FT, FK>(
    target_factory: FT,
    kernel_factory: FK,
    initial_positions: Vec<Vec<f64>>,
    config: TemperingConfig,
) -> Result<TemperingOutput, McmcError>
where
    T: LogDensity<[f64]>,
    K: TransitionKernel<T>,
    FT: Fn(usize, f64) -> T + Sync,
    FK: Fn(usize, f64) -> K + Sync,
{
    let dimension = config.validate(&initial_positions)?;
    let retained = usize::try_from(config.samples)
        .unwrap_or(usize::MAX)
        .div_ceil(config.thinning);

    let mut replicas = initial_positions
        .into_par_iter()
        .enumerate()
        .map(|(slot, initial)| {
            let ladder_value = config.ladder[slot];
            let mut target = target_factory(slot, ladder_value);
            let state = EuclideanState::initialize(&mut target, initial)?;
            let mut trace = MemoryTrace::new(dimension, config.thinning)?;
            trace.reserve_draws(retained);
            Ok(Replica {
                slot,
                ladder_value,
                target,
                kernel: kernel_factory(slot, ladder_value),
                state,
                rng: Xoshiro256PlusPlus::seed_from_u64(chain_seed(config.base_seed, slot)),
                trace,
            })
        })
        .collect::<Result<Vec<_>, McmcError>>()?;

    let mut exchanges = (0..config.ladder.len() - 1)
        .map(|left_index| ExchangeEdgeDiagnostics {
            left_index,
            right_index: left_index + 1,
            attempts: 0,
            acceptances: 0,
        })
        .collect::<Vec<_>>();
    let mut exchange_rng =
        Xoshiro256PlusPlus::seed_from_u64(config.base_seed ^ 0xA076_1D64_78BD_642F);
    let mut exchange_round = 0_u64;
    let mut transitions_since_exchange = 0_u64;

    start_phase(&mut replicas, SamplingPhase::Warmup)?;
    run_phase(
        &mut replicas,
        SamplingPhase::Warmup,
        config.warmup,
        config.exchange_interval,
        &mut transitions_since_exchange,
        &mut exchange_round,
        &mut exchanges,
        &mut exchange_rng,
    )?;
    end_phase(&mut replicas, SamplingPhase::Warmup)?;

    start_phase(&mut replicas, SamplingPhase::Sampling)?;
    run_phase(
        &mut replicas,
        SamplingPhase::Sampling,
        config.samples,
        config.exchange_interval,
        &mut transitions_since_exchange,
        &mut exchange_round,
        &mut exchanges,
        &mut exchange_rng,
    )?;
    end_phase(&mut replicas, SamplingPhase::Sampling)?;

    let outputs = replicas
        .into_iter()
        .map(|replica| TemperedChainOutput {
            slot: replica.slot,
            ladder_value: replica.ladder_value,
            trace: replica.trace,
            final_position: replica.state.position().clone(),
            final_log_density: replica.state.log_density(),
        })
        .collect();
    Ok(TemperingOutput {
        parameter_names: config.parameter_names,
        replicas: outputs,
        exchanges,
    })
}

fn start_phase<T, K>(replicas: &mut [Replica<T, K>], phase: SamplingPhase) -> Result<(), McmcError>
where
    T: LogDensity<[f64]>,
    K: TransitionKernel<T>,
{
    for replica in replicas {
        replica
            .kernel
            .on_phase_start(&mut replica.target, phase, &replica.state)?;
    }
    Ok(())
}

fn end_phase<T, K>(replicas: &mut [Replica<T, K>], phase: SamplingPhase) -> Result<(), McmcError>
where
    T: LogDensity<[f64]>,
    K: TransitionKernel<T>,
{
    for replica in replicas {
        replica
            .kernel
            .on_phase_end(&mut replica.target, phase, &replica.state)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_phase<T, K>(
    replicas: &mut [Replica<T, K>],
    phase: SamplingPhase,
    transitions: u64,
    exchange_interval: u64,
    transitions_since_exchange: &mut u64,
    exchange_round: &mut u64,
    exchanges: &mut [ExchangeEdgeDiagnostics],
    exchange_rng: &mut Xoshiro256PlusPlus,
) -> Result<(), McmcError>
where
    T: LogDensity<[f64]>,
    K: TransitionKernel<T>,
{
    let mut remaining = transitions;
    while remaining > 0 {
        let until_exchange = exchange_interval - *transitions_since_exchange;
        let batch = remaining.min(until_exchange);
        replicas.par_iter_mut().try_for_each(|replica| {
            for _ in 0..batch {
                let report = replica.kernel.transition(
                    &mut replica.target,
                    &mut replica.state,
                    &mut replica.rng,
                    phase,
                )?;
                report.validate()?;
                if phase == SamplingPhase::Sampling {
                    let _retained = replica
                        .trace
                        .record(replica.slot, &replica.state, &report)?;
                }
            }
            Ok::<(), McmcError>(())
        })?;
        remaining -= batch;
        *transitions_since_exchange += batch;

        if *transitions_since_exchange == exchange_interval {
            attempt_neighbor_exchanges(replicas, *exchange_round, exchanges, exchange_rng)?;
            *exchange_round = (*exchange_round).saturating_add(1);
            *transitions_since_exchange = 0;
        }
    }
    Ok(())
}

fn attempt_neighbor_exchanges<T, K, R>(
    replicas: &mut [Replica<T, K>],
    exchange_round: u64,
    exchanges: &mut [ExchangeEdgeDiagnostics],
    rng: &mut R,
) -> Result<(), McmcError>
where
    T: LogDensity<[f64]>,
    K: TransitionKernel<T>,
    R: Rng + ?Sized,
{
    let offset = if replicas.len() == 2 {
        0
    } else {
        (exchange_round & 1) as usize
    };
    let mut left_index = offset;
    while left_index + 1 < replicas.len() {
        let right_index = left_index + 1;
        let (left_slice, right_slice) = replicas.split_at_mut(right_index);
        let left = &mut left_slice[left_index];
        let right = &mut right_slice[0];

        exchanges[left_index].attempts = exchanges[left_index].attempts.saturating_add(1);
        let left_cross = validate_log_density(left.target.log_density(right.state.position()))?;
        let right_cross = validate_log_density(right.target.log_density(left.state.position()))?;
        let log_acceptance =
            left_cross + right_cross - left.state.log_density() - right.state.log_density();
        if log_acceptance.is_nan() {
            return Err(McmcError::InvalidLogDensity {
                value: log_acceptance,
            });
        }
        let accepted = log_acceptance >= 0.0
            || rng.random::<f64>().max(f64::MIN_POSITIVE).ln() < log_acceptance;
        if accepted {
            left.state
                .exchange_position_with(&mut right.state, left_cross, right_cross);
            left.state.cache_mut().invalidate_gradient();
            right.state.cache_mut().invalidate_gradient();
            exchanges[left_index].acceptances = exchanges[left_index].acceptances.saturating_add(1);
        } else {
            left.state.mark_rejected_transition();
            right.state.mark_rejected_transition();
        }
        left_index += 2;
    }
    Ok(())
}
