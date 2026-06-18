//! Pre-built wrapper that composes Hamiltonian + Algorithm into a Carlo.rs [`MonteCarlo`] impl.

use crate::algorithm::Algorithm;
use crate::hamiltonian::{ClusterModel, Hamiltonian, Measurable, Proposable};
use crate::lattice::{build_honeycomb, build_hypercubic, build_kagome, build_triangular};
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

/// Build a lattice from params.
///
/// Uses `lattice_type` param to select the builder:
/// - `"chain"` (default): 1D chain
/// - `"square"`: 2D square (hypercubic)
/// - `"cubic"`: 3D cubic (hypercubic)
/// - `"triangular"`: 2D triangular (PBC only)
/// - `"honeycomb"`: 2D honeycomb (PBC only)
/// - `"kagome"`: 2D kagome (PBC only)
///
/// The non-bravais lattices (triangular/honeycomb/kagome) currently only
/// support periodic boundaries; passing `pbc=false` for those returns an
/// error rather than silently producing a PBC lattice. The hypercubic family
/// honors `pbc`.
///
/// Dimensions via `L` (1D) or `Lx`, `Ly` (2D) or `Lx`, `Ly`, `Lz` (3D).
fn build_lattice_from_params(
    params: &Params,
    pbc: bool,
) -> Result<crate::lattice::CsrLattice, CarloError> {
    let lt = params
        .get::<String>("lattice_type")
        .unwrap_or_else(|| "chain".to_string());

    match lt.as_str() {
        "triangular" | "honeycomb" | "kagome" => {
            // These builders are PBC-only today; reject an explicit open
            // request instead of silently producing a periodic lattice.
            if !pbc {
                return Err(CarloError::InvalidConfig {
                    field: "pbc".into(),
                    reason: format!(
                        "lattice_type `{lt}` only supports periodic boundaries; \
                         set pbc=true (the default) or omit it"
                    ),
                });
            }
            let lx = params.get::<usize>("Lx").unwrap_or(match lt.as_str() {
                "kagome" => 2,
                _ => 4,
            });
            let ly = params.get::<usize>("Ly").unwrap_or(lx);
            Ok(match lt.as_str() {
                "triangular" => build_triangular(lx, ly),
                "honeycomb" => build_honeycomb(lx, ly),
                "kagome" => build_kagome(lx, ly),
                _ => unreachable!(),
            })
        }
        _ => {
            // hypercubic family: chain, square, cubic
            use crate::lattice::BondType;
            let (dims, bond_types) = if let Some(lx) = params.get::<usize>("Lx") {
                let ly = params.get::<usize>("Ly").unwrap_or(lx);
                if let Some(lz) = params.get::<usize>("Lz") {
                    (
                        vec![lx, ly, lz],
                        vec![BondType::SquareX, BondType::SquareY, BondType::SquareZ],
                    )
                } else {
                    (vec![lx, ly], vec![BondType::SquareX, BondType::SquareY])
                }
            } else {
                let l = params.get::<usize>("L").unwrap_or(10);
                (vec![l], vec![BondType::ChainX])
            };
            Ok(build_hypercubic(&dims, &bond_types, pbc))
        }
    }
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
        let lattice = build_lattice_from_params(params, pbc)?;

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

    #[test]
    fn test_lattice_type_triangular() {
        let mut params = Params::new();
        params.set("Lx", 3usize);
        params.set("Ly", 3usize);
        params.set("lattice_type", "triangular");
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42);
        let mc = <ClassicalMC<IsingModel, MetropolisCore> as carlo_rs::FromParams>::from_params(
            &params, &mut rng,
        )
        .unwrap();

        assert_eq!(mc.system.n_sites(), 9);
        assert_eq!(mc.system.lattice.n_bonds, 54);
        assert_eq!(mc.system.lattice.degree(0), 6);
    }

    #[test]
    fn test_lattice_type_honeycomb() {
        let mut params = Params::new();
        params.set("Lx", 4usize);
        params.set("Ly", 4usize);
        params.set("lattice_type", "honeycomb");
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42);
        let mc = <ClassicalMC<IsingModel, MetropolisCore> as carlo_rs::FromParams>::from_params(
            &params, &mut rng,
        )
        .unwrap();

        assert_eq!(mc.system.n_sites(), 16);
        assert_eq!(mc.system.lattice.n_bonds, 48);
        assert_eq!(mc.system.lattice.degree(0), 3);
    }

    #[test]
    fn test_lattice_type_kagome() {
        let mut params = Params::new();
        params.set("Lx", 2usize);
        params.set("Ly", 2usize);
        params.set("lattice_type", "kagome");
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(42);
        let mc = <ClassicalMC<IsingModel, MetropolisCore> as carlo_rs::FromParams>::from_params(
            &params, &mut rng,
        )
        .unwrap();

        assert_eq!(mc.system.n_sites(), 12);
        assert_eq!(mc.system.lattice.n_bonds, 48);
        assert_eq!(mc.system.lattice.degree(0), 4);
    }
}
