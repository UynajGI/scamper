//! Pre-built wrapper that composes Model + Algorithm into a Carlo.rs [`MonteCarlo`] impl.

use crate::algorithm::Algorithm;
use crate::lattice::build_hypercubic;
use crate::model::Model;
use crate::system::System;
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, Params};

/// Pre-composed classical Monte Carlo simulation.
///
/// `M` = physics model (e.g. `IsingModel`), `A` = algorithm (e.g. `MetropolisCore`).
///
/// # Usage
///
/// ```ignore
/// type IsingMetro = ClassicalMC<IsingModel, MetropolisCore>;
/// let params = ...; // L=16, J=1, beta=0.5
/// let scheduler = Scheduler::new(backend, config);
/// let results = scheduler.run_one::<IsingMetro>(&params);
/// ```
pub struct ClassicalMC<M: Model, A: Algorithm<M>> {
    pub system: System,
    pub model: M,
    pub algorithm: A,
}

impl<M: Model, A: Algorithm<M>> ClassicalMC<M, A> {
    pub fn new(system: System, model: M, algorithm: A) -> Self {
        Self {
            system,
            model,
            algorithm,
        }
    }
}

// ── MonteCarlo impl ────────────────────────────────────────

impl<M: Model, A: Algorithm<M>> MonteCarlo for ClassicalMC<M, A> {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.algorithm
            .sweep(&mut self.system, &self.model, &mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure("Energy", self.system.energy);
        let mag = self.model.magnetization(&self.system.spins);
        ctx.measure("Magnetization", mag);
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

impl<M, A> FromParams for ClassicalMC<M, A>
where
    M: Model + FromModelParams,
    A: Algorithm<M> + Default,
{
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let pbc: bool = params
            .get::<String>("pbc")
            .map(|s| s == "true" || s == "1")
            .unwrap_or(true);
        let (dims, bond_types) = parse_lattice(params);
        let lattice = build_hypercubic(&dims, &bond_types, pbc);

        let model = M::from_model_params(params)?;
        let spin_dim = model.spin_dim();
        let algorithm = A::default();

        // Random initial spins
        let mut system = System::new(lattice, spin_dim, 0.0);
        for site in 0..system.n_sites() {
            let spin = model.random_spin(rng);
            system.spin_at_mut(site, spin_dim).copy_from_slice(&spin);
        }
        // Compute initial energy
        system.energy = model.compute_total_energy(&system.spins, &system.lattice);

        Ok(Self {
            system,
            model,
            algorithm,
        })
    }
}

/// Minimal per-model param parsing trait. Separates model params from lattice params.
pub trait FromModelParams: Model + Sized {
    fn from_model_params(params: &Params) -> Result<Self, CarloError>;
}

// ── IsingModel FromModelParams ──────────────────────────────

impl FromModelParams for crate::model::IsingModel {
    fn from_model_params(params: &Params) -> Result<Self, CarloError> {
        let j = params.get::<f64>("J").unwrap_or(1.0);
        let beta = params.get::<f64>("beta").unwrap_or(1.0);
        Ok(Self::new(j, beta))
    }
}

impl FromModelParams for crate::model::PottsModel {
    fn from_model_params(params: &Params) -> Result<Self, CarloError> {
        let j = params.get::<f64>("J").unwrap_or(1.0);
        let beta = params.get::<f64>("beta").unwrap_or(1.0);
        let q = params.get::<usize>("q").unwrap_or(3);
        Ok(Self::new(j, beta, q))
    }
}

impl FromModelParams for crate::model::XYModel {
    fn from_model_params(params: &Params) -> Result<Self, CarloError> {
        let j = params.get::<f64>("J").unwrap_or(1.0);
        let beta = params.get::<f64>("beta").unwrap_or(1.0);
        Ok(Self::new(j, beta))
    }
}

impl FromModelParams for crate::model::HeisenbergModel {
    fn from_model_params(params: &Params) -> Result<Self, CarloError> {
        let j = params.get::<f64>("J").unwrap_or(1.0);
        let beta = params.get::<f64>("beta").unwrap_or(1.0);
        Ok(Self::new(j, beta))
    }
}

// ── Tests ───────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::MetropolisCore;
    use crate::model::IsingModel;
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
