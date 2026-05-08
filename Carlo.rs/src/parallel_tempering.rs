//! Parallel tempering Monte Carlo.
//!
//! Parallel tempering (also known as replica exchange Monte Carlo) runs
//! multiple simulations at different temperatures/parameters, exchanging
//! configurations to improve sampling at difficult parameter values.
//!
//! # How it works
//!
//! 1. Each MPI rank runs one chain at a different parameter value (e.g., temperature)
//! 2. Chains communicate to collect measurements across all temperatures
//! 3. At fixed intervals, neighboring chains attempt Metropolis exchanges
//! 4. The even-odd pairing alternates each exchange step to ensure all neighbors interact
//!
//! # Usage
//!
//! 1. Implement [`ParallelTemperingCompatible`] for your model
//! 2. Create [`ParallelTemperingConfig`] with parameter values
//! 3. Use [`run_parallel_tempering()`] for MPI-based PT simulation
//!
//! # MPI Requirements
//!
//! PT requires the `mpi` feature. The number of MPI ranks must equal
//! the number of parameter values in the config.

use crate::{Context, MonteCarlo};
use rand::Rng;
use rand::SeedableRng;

#[cfg(feature = "mpi")]
use rand::RngExt;

#[cfg(feature = "mpi")]
use crate::{CarloError, FromParams, Metadata, Params, Results};

#[cfg(feature = "mpi")]
use mpi::topology::{Communicator, SimpleCommunicator};

// MPI tags for PT communication
#[cfg(feature = "mpi")]
mod tags {
    pub const PT_WEIGHT_MSG: i32 = 4573792;
    pub const PT_SWITCH_MSG: i32 = 4573793;
    pub const PT_MEASUREMENTS_TAG: i32 = 4573794;
}

// ============================================================================
// Configuration
// ============================================================================

/// Configuration for parallel tempering runs.
#[derive(Debug, Clone)]
pub struct ParallelTemperingConfig {
    /// Name of the parameter to vary (e.g., "beta", "temperature").
    pub parameter: String,
    /// Values of the parameter for each chain.
    pub values: Vec<f64>,
    /// Interval (in sweeps) between exchange attempts.
    pub interval: u64,
}

// ============================================================================
// Compatible trait
// ============================================================================

/// Trait for Monte Carlo implementations compatible with parallel tempering.
///
/// Required methods:
/// - `log_weight_ratio`: Compute log(w(new_param)/w(current_param))
/// - `change_parameter`: Update the model to use a new parameter value
pub trait ParallelTemperingCompatible: MonteCarlo {
    /// Compute log ratio of weights for parameter change.
    ///
    /// Let W(x, p) be the weight of configuration x at parameter p,
    /// and W(x, p') be the weight after changing to new_value.
    /// Returns log(W(x, p') / W(x, p)).
    fn log_weight_ratio(&self, param: &str, new_value: f64) -> f64;

    /// Update the model to use a new value for the given parameter.
    fn change_parameter(&mut self, param: &str, new_value: f64);
}

// ============================================================================
// Parallel Tempering MC wrapper
// ============================================================================

/// Parallel tempering wrapper around a child MC implementation.
///
/// Each instance manages one chain of the parallel tempering simulation.
/// The `chain_idx` identifies which parameter value this chain uses.
pub struct ParallelTemperingMC<MC: ParallelTemperingCompatible> {
    /// Name of the tempered parameter.
    pub parameter_name: String,

    /// All parameter values across chains.
    pub parameter_values: Vec<f64>,

    /// Exchange attempt interval (sweeps).
    pub tempering_interval: u64,

    /// Which chain this instance represents.
    pub chain_idx: usize,

    /// Child MC implementation.
    pub child_mc: MC,

    /// Collected measurements (buffered between sync points).
    pt_measurements: PtMeasurements,
}

/// Buffered measurements for parallel tempering synchronization.
#[derive(Debug, Clone)]
struct PtMeasurements {
    queue: Vec<(String, f64)>,
}

impl PtMeasurements {
    fn new() -> Self {
        Self { queue: Vec::new() }
    }

    fn add(&mut self, name: &str, value: f64) {
        self.queue.push((name.to_string(), value));
    }

    #[cfg(feature = "mpi")]
    fn clear(&mut self) {
        self.queue.clear();
    }
}

impl<MC: ParallelTemperingCompatible> ParallelTemperingMC<MC> {
    /// Create a new PT wrapper.
    pub fn new(config: &ParallelTemperingConfig, chain_idx: usize, child_mc: MC) -> Self {
        Self {
            parameter_name: config.parameter.clone(),
            parameter_values: config.values.clone(),
            tempering_interval: config.interval,
            chain_idx,
            child_mc,
            pt_measurements: PtMeasurements::new(),
        }
    }

    /// Get the current parameter value for this chain.
    pub fn current_value(&self) -> f64 {
        self.parameter_values[self.chain_idx]
    }

    /// Get chain index.
    pub fn chain_idx(&self) -> usize {
        self.chain_idx
    }

    /// Change the current chain's parameter value.
    pub fn set_chain_idx(&mut self, new_idx: usize) {
        self.chain_idx = new_idx;
        self.child_mc
            .change_parameter(&self.parameter_name, self.parameter_values[self.chain_idx]);
    }

    /// Collect a measurement into the PT buffer.
    pub fn pt_measure(&mut self, name: &str, value: f64) {
        self.pt_measurements.add(name, value);
    }

    /// Finalize PT measurements into regular estimates.
    pub fn finalize_pt_measurements<R: Rng + SeedableRng>(&mut self, ctx: &mut Context<R>) {
        for (name, value) in self.pt_measurements.queue.drain(..) {
            ctx.measure(&name, value);
        }
    }
}

// ============================================================================
// MPI PT Exchange
// ============================================================================

/// Manages MPI-based parallel tempering exchanges between chains.
///
/// # Usage
///
/// ```rust,ignore
/// use carlo_rs::parallel_tempering::PtExchange;
///
/// let exchange = PtExchange::new(&comm, &config, &params, seed, binsize)?;
///
/// loop {
///     exchange.step(&mut ctx)?;
///     if exchange.should_exchange() {
///         exchange.try_exchange()?;
///     }
///     if exchange.is_complete() {
///         break;
///     }
/// }
///
/// let results = exchange.finalize();
/// ```
#[cfg(feature = "mpi")]
pub struct PtExchange<MC: ParallelTemperingCompatible, R: Rng + SeedableRng> {
    /// Communicator for all PT chains.
    comm: SimpleCommunicator,

    /// PT MC instance for this chain.
    mc: ParallelTemperingMC<MC>,

    /// Context for this chain.
    ctx: Context<R>,

    /// Target sweeps per chain.
    target_sweeps: u64,

    /// Completed sweeps.
    sweeps_done: u64,

    /// Base seed for reproducibility.
    base_seed: u64,
}

#[cfg(feature = "mpi")]
impl<MC: ParallelTemperingCompatible + FromParams<Rng = R>, R: Rng + SeedableRng + Send>
    PtExchange<MC, R>
{
    /// Create a new PT exchange instance.
    ///
    /// # Arguments
    /// * `comm` - MPI communicator containing all chains (one rank per chain)
    /// * `config` - PT configuration with parameter values
    /// * `params` - Base parameters (will be modified with chain-specific value)
    /// * `seed` - Base random seed
    /// * `binsize` - Binning size for measurements
    /// * `target_sweeps` - Number of sweeps to run
    pub fn new(
        comm: SimpleCommunicator,
        config: &ParallelTemperingConfig,
        params: &Params,
        seed: u64,
        binsize: usize,
        target_sweeps: u64,
    ) -> Result<Self, CarloError> {
        let rank = comm.rank();
        let n_chains = comm.size();

        if n_chains as usize != config.values.len() {
            return Err(CarloError::InvalidConfig {
                field: "pt_chains".into(),
                reason: format!(
                    "Number of MPI ranks ({}) != number of PT values ({})",
                    n_chains,
                    config.values.len()
                ),
            });
        }

        let chain_idx = rank as usize;

        // Create chain-specific parameters
        let mut chain_params = params.clone();
        chain_params.set(&config.parameter, config.values[chain_idx].to_string());

        // Create RNG with chain-specific seed
        let chain_seed = seed.wrapping_add(chain_idx as u64 * 100000);
        let rng = R::seed_from_u64(chain_seed);

        let mut ctx = Context::new_with_binsize(rng, 0, binsize);
        let child_mc = MC::from_params(&chain_params, &mut ctx.rng)?;

        let mc = ParallelTemperingMC::new(config, chain_idx, child_mc);

        Ok(Self {
            comm,
            mc,
            ctx,
            target_sweeps,
            sweeps_done: 0,
            base_seed: seed,
        })
    }

    /// Execute one PT step (sweep + optional exchange).
    pub fn step(&mut self) {
        self.mc.child_mc.sweep(&mut self.ctx);
        self.ctx.advance_sweep();

        // Collect measurements if thermalized
        if self.ctx.is_thermalized() {
            self.mc.child_mc.measure(&mut self.ctx);
        }

        // Check if we should attempt exchange
        if self.ctx.sweep_count() > 0 && self.ctx.sweep_count() % self.mc.tempering_interval == 0 {
            let _ = self.try_exchange();
        }

        self.sweeps_done += 1;
    }

    /// Try to exchange with a neighbor chain.
    ///
    /// Uses even-odd pairing: at even exchange steps, chain 0 pairs with 1,
    /// chain 2 with 3, etc. At odd steps, chain 1 pairs with 2, etc.
    pub fn try_exchange(&mut self) -> Result<bool, CarloError> {
        let n_chains = self.comm.size() as usize;
        if n_chains < 2 {
            return Ok(false);
        }

        let exchange_step = self.ctx.sweep_count() / self.mc.tempering_interval;
        let pairing_offset = exchange_step & 1;

        // Determine partner
        let my_chain_idx = self.mc.chain_idx;
        let partner_chain_idx = if my_chain_idx % 2 == pairing_offset as usize {
            // Try to pair with higher index
            if my_chain_idx + 1 < n_chains {
                my_chain_idx + 1
            } else {
                return Ok(false); // No partner available
            }
        } else {
            // Try to pair with lower index
            if my_chain_idx > 0 {
                my_chain_idx - 1
            } else {
                return Ok(false); // No partner available
            }
        };

        // Find partner's rank in communicator
        let partner_rank = self.find_rank_for_chain(partner_chain_idx)?;

        // Compute log weight ratio
        let w = self.mc.child_mc.log_weight_ratio(
            &self.mc.parameter_name,
            self.mc.parameter_values[partner_chain_idx],
        );

        // Exchange weights and decide
        let accept = if my_chain_idx % 2 == pairing_offset as usize {
            // Even/offset chains receive first
            let partner_w = self.recv_weight(partner_rank)?;
            let accept = self.ctx.rng.random::<f64>() < (w + partner_w).exp();
            self.send_switch(partner_rank, accept)?;
            accept
        } else {
            // Odd chains send weight first
            self.send_weight(partner_rank, w)?;
            self.recv_switch(partner_rank)?
        };

        if accept {
            self.mc.set_chain_idx(partner_chain_idx);
        }

        Ok(accept)
    }

    /// Synchronize measurements across all chains.
    ///
    /// Gathers PT-buffered measurements from all chains and records
    /// them with proper chain permutation tracking.
    pub fn synchronize_measurements(&mut self) -> Result<(), CarloError> {
        use mpi::traits::*;

        let rank = self.comm.rank();
        let n_chains = self.comm.size();

        // Serialize measurement queue to bytes for MPI transfer
        let local_data = serde_json::to_string(&self.mc.pt_measurements.queue).map_err(|e| {
            CarloError::InvalidConfig {
                field: "pt_measurements".into(),
                reason: format!("Failed to serialize measurements: {}", e),
            }
        })?;
        let local_bytes = local_data.as_bytes().to_vec();
        self.mc.pt_measurements.clear();

        if rank == 0 {
            // Root collects from all other ranks
            let mut all_queues: Vec<Vec<(String, f64)>> = Vec::with_capacity(n_chains as usize);
            // Root's own data
            let root_data: Vec<(String, f64)> =
                serde_json::from_slice(&local_bytes).map_err(|e| CarloError::InvalidConfig {
                    field: "pt_measurements".into(),
                    reason: format!("Failed to deserialize root measurements: {}", e),
                })?;
            all_queues.push(root_data);

            for src in 1..n_chains {
                let (bytes, _) = self
                    .comm
                    .process_at_rank(src)
                    .receive_vec_with_tag::<u8>(tags::PT_MEASUREMENTS_TAG);
                let queue: Vec<(String, f64)> =
                    serde_json::from_slice(&bytes).map_err(|e| CarloError::InvalidConfig {
                        field: "pt_measurements".into(),
                        reason: format!(
                            "Failed to deserialize measurements from rank {}: {}",
                            src, e
                        ),
                    })?;
                all_queues.push(queue);
            }

            // Verify all chains recorded measurements in the same order
            let ref_names: Vec<String> = all_queues[0].iter().map(|(n, _)| n.clone()).collect();
            for q in &all_queues {
                let names: Vec<String> = q.iter().map(|(n, _)| n.clone()).collect();
                if names != ref_names {
                    return Err(CarloError::InvalidConfig {
                        field: "pt_measurements".into(),
                        reason: "Measurement order differs between chains".into(),
                    });
                }
            }

            // For each observable, collect values across chains
            for (i, (name, _)) in all_queues[0].iter().enumerate() {
                let values: Vec<f64> = all_queues.iter().map(|q| q[i].1).collect();
                for &v in &values {
                    self.ctx.measure(name, v);
                }
            }
        } else {
            // Non-root ranks send their queues to root
            self.comm
                .process_at_rank(0)
                .send_with_tag(local_bytes.as_slice(), tags::PT_MEASUREMENTS_TAG);
            self.mc.pt_measurements.clear();
        }

        Ok(())
    }

    /// Check if simulation is complete.
    pub fn is_complete(&self) -> bool {
        self.sweeps_done >= self.target_sweeps
    }

    /// Get sweep count.
    pub fn sweep_count(&self) -> u64 {
        self.sweeps_done
    }

    /// Get current chain index.
    pub fn chain_idx(&self) -> usize {
        self.mc.chain_idx
    }

    /// Get current parameter value.
    pub fn current_value(&self) -> f64 {
        self.mc.current_value()
    }

    /// Get context reference.
    pub fn context(&self) -> &Context<R> {
        &self.ctx
    }

    /// Finalize and return results.
    pub fn finalize(mut self) -> Results {
        // Flush any remaining PT measurements
        self.mc.finalize_pt_measurements(&mut self.ctx);

        let estimates = self.ctx.finalize_measurements();
        let mut results = Results::from_measurements(&estimates);
        results.set_metadata(Metadata {
            version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp: chrono::Utc::now(),
            base_seed: self.base_seed,
            thermalization_sweeps: self.ctx.thermalization_sweeps(),
            measurement_sweeps: self.sweeps_done,
            n_tasks: self.comm.size() as usize,
        });
        results
    }

    // ---- MPI helper methods ----

    fn find_rank_for_chain(&self, target_chain_idx: usize) -> Result<i32, CarloError> {
        // In PT mode, each rank corresponds to one chain, so rank == chain_idx
        Ok(target_chain_idx as i32)
    }

    fn send_weight(&self, dest: i32, weight: f64) -> Result<(), CarloError> {
        use mpi::traits::*;
        self.comm
            .process_at_rank(dest)
            .send_with_tag(&weight, tags::PT_WEIGHT_MSG);
        Ok(())
    }

    fn recv_weight(&self, source: i32) -> Result<f64, CarloError> {
        use mpi::traits::*;
        let (weight, _) = self
            .comm
            .process_at_rank(source)
            .receive_with_tag(tags::PT_WEIGHT_MSG);
        Ok(weight)
    }

    fn send_switch(&self, dest: i32, accept: bool) -> Result<(), CarloError> {
        use mpi::traits::*;
        self.comm
            .process_at_rank(dest)
            .send_with_tag(&accept, tags::PT_SWITCH_MSG);
        Ok(())
    }

    fn recv_switch(&self, source: i32) -> Result<bool, CarloError> {
        use mpi::traits::*;
        let (accept, _) = self
            .comm
            .process_at_rank(source)
            .receive_with_tag(tags::PT_SWITCH_MSG);
        Ok(accept)
    }
}

// ============================================================================
// MPI Entry Point
// ============================================================================

/// Run parallel tempering simulation using MPI.
///
/// Each MPI rank runs one chain at a different parameter value.
/// Chains exchange configurations at fixed intervals to improve sampling.
///
/// # Requirements
/// - Number of MPI ranks == number of parameter values
/// - Requires `mpi` feature
#[cfg(feature = "mpi")]
pub fn run_parallel_tempering<MC, R>(
    config: &ParallelTemperingConfig,
    params: &Params,
    seed: u64,
    binsize: usize,
    target_sweeps: u64,
    thermalization_sweeps: u64,
) -> Result<Option<Results>, CarloError>
where
    MC: ParallelTemperingCompatible + FromParams<Rng = R>,
    R: Rng + SeedableRng + Send,
{
    use mpi::traits::*;

    let universe = mpi::initialize().ok_or(MpiError::InitFailed)?;
    let world = universe.world();

    // Use the world communicator for PT
    let pt_comm = world;

    // Override thermalization in context
    let exchange = PtExchange::<MC, R>::new(pt_comm, config, params, seed, binsize, target_sweeps)?;

    // Set thermalization in context
    // (PtExchange creates context with 0 thermalization; we adjust here)

    let n_ranks = exchange.comm.size();
    let rank = exchange.comm.rank();

    if rank == 0 {
        eprintln!(
            "Parallel tempering: {} chains, {} sweeps, interval = {}",
            n_ranks, target_sweeps, config.interval
        );
        for (i, &v) in config.values.iter().enumerate() {
            eprintln!("  Chain {}: {} = {}", i, config.parameter, v);
        }
    }

    let mut exchange = exchange;
    exchange.ctx = Context::new_with_binsize(exchange.ctx.rng, thermalization_sweeps, binsize);

    while !exchange.is_complete() {
        exchange.step();
    }

    let results = exchange.finalize();

    // Only rank 0 returns results
    if rank == 0 {
        Ok(Some(results))
    } else {
        Ok(None)
    }
}

// ============================================================================
// Backward compatibility: non-MPI stubs
// ============================================================================

/// MPI Error type (re-export from mpi module)
#[cfg(feature = "mpi")]
pub use crate::backend::MpiError;

#[cfg(not(feature = "mpi"))]
/// Stub error when MPI is not available
#[derive(Debug, thiserror::Error)]
pub enum MpiError {
    #[error("MPI feature not enabled")]
    NotAvailable,
}
