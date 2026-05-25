//! Pre-built wrapper that composes Hamiltonian + Algorithm into a Carlo.rs [`MonteCarlo`] impl.

use crate::algorithm::Algorithm;
use crate::hamiltonian::{ClusterModel, Hamiltonian, Measurable, Proposable};
use crate::lattice::build_hypercubic;
use crate::observables::DefaultObservableSet;
use crate::system::System;
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, ParallelTemperingCompatible, Params};
use serde_json::Value as Json;

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
// ── JSON snapshot save/load ──────────────────────────────────
impl<H, A> ClassicalMC<H, A>
where
    H: Hamiltonian + Measurable,
    A: Algorithm<H>,
{
    /// Save full simulation state as JSON value.
    pub fn save_snapshot(&self) -> Json {
        serde_json::json!({
            "spins": &self.system.spins,
            "energy": self.system.energy,
            "beta": self.system.beta,
            "spin_dim": self.model.spin_dim(),
            "n_sites": self.system.n_sites(),
            "offsets": &self.system.lattice.offsets,
            "neighbors": &self.system.lattice.neighbors,
        })
    }

    /// Load simulation state from JSON value.
    pub fn load_snapshot(&mut self, snapshot: &Json) -> Result<(), CarloError> {
        let spins: Vec<f64> = snapshot["spins"]
            .as_array()
            .ok_or_else(|| CarloError::InvalidConfig {
                field: "snapshot.spins".into(),
                reason: "missing or invalid".into(),
            })?
            .iter()
            .map(|v| v.as_f64().unwrap_or(0.0))
            .collect();
        let energy = snapshot["energy"].as_f64().unwrap_or(0.0);
        let beta = snapshot["beta"].as_f64().unwrap_or(1.0);

        if spins.len() != self.system.spins.len() {
            return Err(CarloError::InvalidConfig {
                field: "snapshot.spins".into(),
                reason: format!(
                    "spin count mismatch: expected {}, got {}",
                    self.system.spins.len(),
                    spins.len()
                ),
            });
        }

        self.system.spins.copy_from_slice(&spins);
        self.system.energy = energy;
        self.system.beta = beta;

        Ok(())
    }
}

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

    pub fn with_observables(
        system: System,
        model: H,
        algorithm: A,
        observables: DefaultObservableSet<H>,
    ) -> Self {
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
        let mut e = 0.0;
        let mut m = 0.0;
        for obs in self.observables.iter() {
            let v = obs.measure(&self.system, &self.model);
            match obs.name() {
                "Energy" => e = v,
                "Magnetization" => m = v,
                _ => {}
            }
            ctx.measure(obs.name(), v);
        }
        ctx.measure("E2", e * e);
        ctx.measure("M2", m * m);
        ctx.measure("M4", m * m * m * m);
    }

    fn name(&self) -> &'static str {
        self.algorithm.name()
    }
}

// ── ParallelTemperingCompatible impl ─────────────────────────

impl<H, A> ParallelTemperingCompatible for ClassicalMC<H, A>
where
    H: Hamiltonian + Measurable,
    A: Algorithm<H>,
{
    fn log_weight_ratio(&self, param: &str, new_value: f64) -> f64 {
        match param {
            "beta" => {
                // W(x, β) = exp(-βE), log(W'/W) = -(β'-β)E
                (self.system.beta - new_value) * self.system.energy
            }
            _ => panic!("unsupported PT param: {param}"),
        }
    }

    fn change_parameter(&mut self, param: &str, new_value: f64) {
        match param {
            "beta" => {
                self.system.beta = new_value;
                // Recompute energy (for most classical models E is β-independent,
                // but we recompute for correctness with future models).
                self.system.energy = self.model.compute_total_energy(
                    &self.system.spins,
                    &self.system.lattice,
                    new_value,
                );
            }
            _ => panic!("unsupported PT param: {param}"),
        }
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
    use carlo_rs::{ParallelTemperingCompatible, RayonBackend, RunConfig, Scheduler};
    use rand::SeedableRng;

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

    #[test]
    fn test_pt_log_weight_ratio() {
        let mut params = Params::new();
        params.set("L", 4usize);
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42);
        let mc = <ClassicalMC<IsingModel, MetropolisCore> as carlo_rs::FromParams>::from_params(
            &params, &mut rng,
        )
        .unwrap();

        let e = mc.system.energy;
        // log(W(β')/W(β)) = -(β'-β)E = (β - β')E
        let lr = mc.log_weight_ratio("beta", 2.0);
        let expected = (1.0 - 2.0) * e;
        assert!(
            (lr - expected).abs() < 1e-10,
            "lr = {}, expected = {}",
            lr,
            expected
        );
    }

    #[test]
    fn test_moment_observables_present() {
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let config = RunConfig {
            thermalization_sweeps: 100,
            measurement_sweeps: 200,
            binsize: 50,
            base_seed: 42,
            ..Default::default()
        };

        let backend = RayonBackend::new(1);
        let scheduler = Scheduler::new(backend, config);
        let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

        let e2 = results.get("E2").expect("E2 observable missing");
        let m2 = results.get("M2").expect("M2 observable missing");
        let m4 = results.get("M4").expect("M4 observable missing");

        assert!(e2.mean > 0.0, "E² should be positive");
        assert!((0.0..=1.0).contains(&m2.mean), "M² should be in [0,1]");
        assert!((0.0..=1.0).contains(&m4.mean), "M⁴ should be in [0,1]");
    }

    #[test]
    fn test_snapshot_round_trip() {
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42);
        let mut mc =
            <ClassicalMC<IsingModel, MetropolisCore> as carlo_rs::FromParams>::from_params(
                &params, &mut rng,
            )
            .unwrap();

        // Run a few sweeps
        let mut ctx = carlo_rs::Context::new(rng, 0);
        for _ in 0..10 {
            mc.sweep(&mut ctx);
        }

        let energy_before = mc.system.energy;
        let spins_before = mc.system.spins.clone();

        // Save snapshot
        let snapshot = mc.save_snapshot();

        // Mutate state
        mc.system.energy = 999.0;
        mc.system.spins.fill(0.0);

        // Load snapshot
        mc.load_snapshot(&snapshot).unwrap();

        assert!((mc.system.energy - energy_before).abs() < 1e-10);
        assert_eq!(mc.system.spins, spins_before);
    }

    #[test]
    fn test_pt_change_parameter() {
        let mut params = Params::new();
        params.set("L", 4usize);
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42);
        let mut mc =
            <ClassicalMC<IsingModel, MetropolisCore> as carlo_rs::FromParams>::from_params(
                &params, &mut rng,
            )
            .unwrap();

        mc.change_parameter("beta", 2.5);
        assert!((mc.system.beta - 2.5).abs() < 1e-10);
        // Energy should be recomputed (same config, same energy for Ising)
        let e = mc
            .model
            .compute_total_energy(&mc.system.spins, &mc.system.lattice, 2.5);
        assert!((mc.system.energy - e).abs() < 1e-10);
    }
}
