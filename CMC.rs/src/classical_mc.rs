//! Pre-built wrapper that composes Hamiltonian + Algorithm into a Carlo.rs [`MonteCarlo`] impl.

use crate::algorithm::Algorithm;
use crate::hamiltonian::{ClusterModel, Hamiltonian, Measurable, Proposable};
use crate::lattice::build_hypercubic;
use crate::observables::DefaultObservableSet;
use crate::system::System;
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, Params};

/// Pre-composed classical Monte Carlo simulation.
///
/// `H` = Hamiltonian (e.g. `IsingModel`), `A` = algorithm (e.g. `MetropolisCore`).
///
/// # Usage
///
/// ```ignore
/// type IsingMetro = ClassicalMC<IsingModel, MetropolisCore>;
/// let params = ...; // L=16, J=1, beta=0.5
/// let scheduler = Scheduler::new(backend, config);
/// let results = scheduler.run_one::<IsingMetro>(&params);
/// ```
pub struct ClassicalMC<H, A>
where
    H: Hamiltonian + Measurable,
    A: Algorithm<H>,
{
    pub system: System,
    pub model: H,
    pub algorithm: A,
    pub observables: DefaultObservableSet<H>,
}

impl<H, A> ClassicalMC<H, A>
where
    H: Hamiltonian + Measurable,
    A: Algorithm<H>,
{
    pub fn new(system: System, model: H, algorithm: A) -> Self {
        Self {
            system,
            model,
            algorithm,
            observables: DefaultObservableSet::new(),
        }
    }

    pub fn with_observables(system: System, model: H, algorithm: A, observables: DefaultObservableSet<H>) -> Self {
        Self {
            system,
            model,
            algorithm,
            observables,
        }
    }
}

// ── MonteCarlo impl ────────────────────────────────────────

impl<H, A> MonteCarlo for ClassicalMC<H, A>
where
    H: Hamiltonian + Measurable,
    A: Algorithm<H>,
{
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.algorithm
            .sweep(&mut self.system, &self.model, &mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        for obs in self.observables.iter() {
            ctx.measure(obs.name(), obs.measure(&self.system, &self.model));
        }
    }

    fn name(&self) -> &'static str {
        self.algorithm.name()
    }
}

// ── FromParams impl ─────────────────────────────────────────

/// Parse lattice dimensions from params.
///
/// - `Lx, Ly` → 2D square
/// - `Lx, Ly, Lz` → 3D cubic
/// - `L` → 1D chain (fallback)
fn parse_lattice(params: &Params) -> (Vec<usize>, Vec<crate::lattice::BondType>) {
    use crate::lattice::BondType;

    if let Some(lx) = params.get::<usize>("Lx") {
        let ly = params.get::<usize>("Ly").unwrap_or(lx);
        if let Some(lz) = params.get::<usize>("Lz") {
            return (
                vec![lx, ly, lz],
                vec![BondType::SquareX, BondType::SquareY, BondType::SquareZ],
            );
        }
        return (vec![lx, ly], vec![BondType::SquareX, BondType::SquareY]);
    }

    let l = params.get::<usize>("L").unwrap_or(10);
    (vec![l], vec![BondType::ChainX])
}

impl<H, A> FromParams for ClassicalMC<H, A>
where
    H: Hamiltonian + Measurable + Proposable + ClusterModel + FromHamiltonianParams,
    A: Algorithm<H> + Default,
{
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let pbc: bool = params
            .get::<String>("pbc")
            .map(|s| s == "true" || s == "1")
            .unwrap_or(true);
        let (dims, bond_types) = parse_lattice(params);
        let lattice = build_hypercubic(&dims, &bond_types, pbc);

        let beta = params.get::<f64>("beta").unwrap_or(1.0);
        let model = H::from_hamiltonian_params(params)?;
        let spin_dim = model.spin_dim();
        let algorithm = A::default();

        // Random initial spins
        let mut system = System::new(lattice, spin_dim, 0.0, beta);
        for site in 0..system.n_sites() {
            let spin = model.propose(rng);
            system.spin_at_mut(site, spin_dim).copy_from_slice(&spin);
        }
        // Compute initial energy
        system.energy = model.compute_total_energy(&system.spins, &system.lattice, beta);

        Ok(Self {
            system,
            model,
            algorithm,
            observables: DefaultObservableSet::new(),
        })
    }
}

/// Minimal per-model param parsing trait. Separates model params (J, q, ...) from
/// lattice params (L, dims) and temperature (β).
pub trait FromHamiltonianParams: Hamiltonian + Sized {
    fn from_hamiltonian_params(params: &Params) -> Result<Self, CarloError>;
}

// ── IsingModel FromHamiltonianParams ────────────────────────

impl FromHamiltonianParams for crate::models::IsingModel {
    fn from_hamiltonian_params(params: &Params) -> Result<Self, CarloError> {
        let j = params.get::<f64>("J").unwrap_or(1.0);
        Ok(Self::new(j))
    }
}

impl FromHamiltonianParams for crate::models::PottsModel {
    fn from_hamiltonian_params(params: &Params) -> Result<Self, CarloError> {
        let j = params.get::<f64>("J").unwrap_or(1.0);
        let q = params.get::<usize>("q").unwrap_or(3);
        Ok(Self::new(j, q))
    }
}

impl FromHamiltonianParams for crate::models::XYModel {
    fn from_hamiltonian_params(params: &Params) -> Result<Self, CarloError> {
        let j = params.get::<f64>("J").unwrap_or(1.0);
        Ok(Self::new(j))
    }
}

impl FromHamiltonianParams for crate::models::HeisenbergModel {
    fn from_hamiltonian_params(params: &Params) -> Result<Self, CarloError> {
        let j = params.get::<f64>("J").unwrap_or(1.0);
        Ok(Self::new(j))
    }
}

// ── Tests ───────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::MetropolisCore;
    use crate::models::IsingModel;
    use carlo_rs::{RayonBackend, RunConfig, Scheduler};

    #[test]
    fn test_classical_mc_end_to_end() {
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let config = RunConfig {
            thermalization_sweeps: 200,
            measurement_sweeps: 500,
            binsize: 50,
            base_seed: 42,
            ..Default::default()
        };

        let backend = RayonBackend::new(1);
        let scheduler = Scheduler::new(backend, config);
        let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

        let energy = results.get("Energy").expect("Energy observable missing");
        let mag = results
            .get("Magnetization")
            .expect("Magnetization observable missing");

        // At beta=1 (T=1), ferromagnetic Ising chain should have negative energy
        assert!(energy.mean < 0.0);
        assert!(energy.stderr > 0.0);
        assert!((0.0..=1.0).contains(&mag.mean));
    }
}
