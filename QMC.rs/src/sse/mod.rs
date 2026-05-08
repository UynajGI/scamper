//! SSE (Stochastic Series Expansion) algorithm module.

mod diagonal;
mod engine;
mod estimators;
mod loop_;
mod measurements;
pub mod vertex_data;
pub mod vertex_list;

pub use engine::{OperatorSequence, SSEEngine, Vertex};
pub use estimators::CorrelationResult;
pub use vertex_data::VertexData;
pub use vertex_list::VertexList;
pub use crate::hilbert::OpType; // Re-export for convenience

use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use crate::hilbert::HilbertSpace;
use crate::lattice::{BondType, Lattice};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::HashMap;

/// Domain trait for lattice QMC methods.
/// Note: This trait does NOT require MonteCarlo - SSECore implements MonteCarlo.
pub trait LatticeQMC {
    /// Access lattice topology.
    fn lattice(&self) -> &Lattice;

    /// Total number of sites.
    fn n_sites(&self) -> usize {
        self.lattice().n_sites
    }
}

/// Method trait for SSE Monte Carlo.
/// Models implement this trait to define their physics.
pub trait SSEMonteCarlo: LatticeQMC {
    /// Associated HilbertSpace type.
    type HilbertSpace: HilbertSpace;

    /// Define operators on each bond type.
    /// Returns (OpType, coupling_constant) pairs.
    /// The coupling_constant should be the bare coupling (e.g., J for Heisenberg).
    /// The diagonal matrix element computation in the HilbertSpace handles the
    /// shifted form (s1*s2 + 1)/4 to ensure positive weights.
    fn bond_operators(&self, bond_type: BondType) -> Vec<(OpType, f64)>;

    /// Access HilbertSpace implementation.
    fn hilbert_space(&self) -> &Self::HilbertSpace;

    /// Simulation inverse temperature.
    fn beta(&self) -> f64;

    /// Total diagonal shift for the shifted Hamiltonian.
    /// H' = H + C, so E = -<n>/beta - C.
    fn diagonal_shift(&self) -> f64 {
        0.0
    }

    /// Loop parameter (e.g., XXZ anisotropy Δ) used by directed loop algorithm.
    /// Default is 1.0 (Heisenberg point, no bounce needed).
    fn loop_parameter(&self) -> f64 {
        1.0
    }
}

/// Trait for models that can be created from Params for use with Scheduler.
/// SSECore requires the inner model to implement this.
pub trait SSEMonteCarloFromParams: SSEMonteCarlo + FromParams {}

/// SSE core wrapper providing default MonteCarlo implementation.
pub struct SSECore<MC: SSEMonteCarlo> {
    /// SSE engine
    pub engine: SSEEngine<MC::HilbertSpace>,
    /// User model
    pub mc: MC,
    /// Specific heat accumulators: (<n>, <n^2>) accumulated per measurement sweep.
    n_sum: f64,
    n_sq_sum: f64,
    n_measurements: usize,
}

impl<MC: SSEMonteCarlo> SSECore<MC> {
    /// Create SSE core from user model.
    pub fn new(mc: MC) -> Self {
        let lattice = mc.lattice().clone();
        let beta = mc.beta();
        let n_sites = lattice.n_sites;

        // Estimate max_length: M ~ N x beta x (average operators per site)
        let max_length = (n_sites as f64 * beta * 2.0) as usize + 100;

        let hs = mc.hilbert_space().clone();

        // Build weights from bond_operators
        let weights = Self::build_weights(&mc, &lattice);

        let engine = SSEEngine::new(lattice, hs, max_length, weights, beta, mc.diagonal_shift(), mc.loop_parameter());

        SSECore {
            engine,
            mc,
            n_sum: 0.0,
            n_sq_sum: 0.0,
            n_measurements: 0,
        }
    }

    /// Collect bond weights from model.
    fn build_weights(mc: &MC, lattice: &Lattice) -> HashMap<BondType, f64> {
        let mut weights = HashMap::new();

        for site in &lattice.sites {
            for neighbor in site {
                if !weights.contains_key(&neighbor.bond_type) {
                    let ops = mc.bond_operators(neighbor.bond_type);
                    let diag_weight: f64 = ops
                        .iter()
                        .filter(|(op, _)| *op == OpType::Diagonal)
                        .map(|(_, w)| w)
                        .sum();
                    let offdiag_weight: f64 = ops
                        .iter()
                        .filter(|(op, _)| *op == OpType::OffDiagonal)
                        .map(|(_, w)| w)
                        .sum();
                    // Store diagonal weight, but if it's zero use off-diagonal weight.
                    // This ensures pure off-diagonal models (XY) have non-zero weights
                    // for the diagonal update to insert off-diagonal operators.
                    let weight = if diag_weight > 1e-10 {
                        diag_weight
                    } else {
                        offdiag_weight
                    };
                    weights.insert(neighbor.bond_type, weight);
                }
            }
        }

        weights
    }

    /// Compute specific heat from accumulated operator count statistics.
    ///
    /// C = (<n^2> - <n>^2 - <n>) / N
    ///
    /// Returns None if fewer than 2 measurements.
    pub fn compute_specific_heat(&self) -> Option<f64> {
        if self.n_measurements < 2 {
            return None;
        }
        let n_mean = self.n_sum / self.n_measurements as f64;
        let n_sq_mean = self.n_sq_sum / self.n_measurements as f64;
        let ns = self.engine.lattice.n_sites as f64;
        Some((n_sq_mean - n_mean * n_mean - n_mean) / ns)
    }

    /// Reset specific heat accumulators.
    pub fn reset_specific_heat(&mut self) {
        self.n_sum = 0.0;
        self.n_sq_sum = 0.0;
        self.n_measurements = 0;
    }
}

impl<MC: SSEMonteCarlo> MonteCarlo for SSECore<MC> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.engine.diagonal_update(&mut ctx.rng);
        self.engine.loopupdate(&mut ctx.rng);
        ctx.advance_sweep();
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        // Energy (diagnostic estimator)
        let energy = self.engine.compute_energy();
        ctx.measure("Energy", energy);

        // Magnetization
        let magnetization = self.engine.compute_magnetization();
        ctx.measure("Magnetization", magnetization);

        // Correlation functions and susceptibilities
        if let Some(result) = self.engine.measure_correlations() {
            // Staggered susceptibility
            ctx.measure("StaggeredSusceptibility", result.staggered_susceptibility);

            // Structure factor (array observable)
            ctx.measure_array("StructureFactor", &result.structure_factor);

            // Correlation function (array observable)
            ctx.measure_array("Correlation", &result.correlation);
        }

        // Specific heat from operator count fluctuations
        let n = self.engine.op_seq.n_operators as f64;
        ctx.measure("OperatorCount", n);
        self.n_sum += n;
        self.n_sq_sum += n * n;
        self.n_measurements += 1;

        // Report specific heat if we have enough data
        if let Some(cv) = self.compute_specific_heat() {
            ctx.measure("SpecificHeat", cv);
        }
    }

    fn name(&self) -> &'static str {
        "SSECore"
    }
}

impl<MC: SSEMonteCarlo + FromParams<Rng = Xoshiro256PlusPlus>> FromParams for SSECore<MC> {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        // Delegate to inner model's from_params
        let mc = MC::from_params(params, rng)?;
        Ok(Self::new(mc))
    }
}
