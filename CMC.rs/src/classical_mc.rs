//! Carlo.rs adapter for composing a physical model, update kernel and observables.

use crate::algorithms::{Algorithm, SimulationPhase};
use crate::lattice::graph::{
    build_chain, build_honeycomb, build_hypercubic, build_kagome, build_square, build_triangular,
    BondType, CsrLattice,
};
use crate::lattice::interaction::{Hamiltonian, Initializable, Measurable};
use crate::lattice::state::System;
use crate::observables::{DefaultObservableSet, ObservableSet};
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, ParallelTemperingCompatible, Params};
use serde_json::Value as Json;
use std::{fmt::Display, str::FromStr};

/// Pre-composed classical Monte Carlo chain.
///
/// Carlo.rs remains responsible for RNG ownership, thermalization/measurement
/// scheduling, accumulation, parallel backends and parallel tempering.  CMC.rs
/// supplies the physical state transition and observable implementations.
pub struct ClassicalMC<H, A, O = DefaultObservableSet<H>>
where
    H: Hamiltonian,
    A: Algorithm<H>,
    O: ObservableSet<H>,
{
    pub system: System,
    pub model: H,
    pub algorithm: A,
    pub observables: O,
}

impl<H, A> ClassicalMC<H, A, DefaultObservableSet<H>>
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
}

impl<H, A, O> ClassicalMC<H, A, O>
where
    H: Hamiltonian,
    A: Algorithm<H>,
    O: ObservableSet<H>,
{
    pub fn with_observables(system: System, model: H, algorithm: A, observables: O) -> Self {
        Self {
            system,
            model,
            algorithm,
            observables,
        }
    }

    /// Versioned JSON state snapshot.  Topology metadata is included for
    /// validation; the lattice itself remains constructor-owned.
    pub fn save_snapshot(&self) -> Json {
        let edges = self
            .system
            .lattice
            .edges
            .iter()
            .map(|edge| {
                serde_json::json!({
                    "source": edge.source,
                    "target": edge.target,
                    "kind": edge.kind.as_label(),
                    "weight": edge.weight,
                })
            })
            .collect::<Vec<_>>();
        serde_json::json!({
            "format": "cmc-rs-snapshot-v2",
            "spins": &self.system.spins,
            "beta": self.system.beta,
            "spin_dim": self.model.spin_dim(),
            "n_sites": self.system.n_sites(),
            "n_edges": self.system.lattice.n_edges(),
            "offsets": &self.system.lattice.offsets,
            "neighbors": &self.system.lattice.neighbors,
            "edge_ids": &self.system.lattice.edge_ids,
            "edges": edges,
        })
    }

    /// Restore a snapshot and verify it against the current model/lattice.
    /// Cached energy is recomputed rather than blindly trusted.
    pub fn load_snapshot(&mut self, snapshot: &Json) -> Result<(), CarloError> {
        let format = snapshot["format"].as_str().unwrap_or("");
        if format != "cmc-rs-snapshot-v2" {
            return Err(CarloError::CheckpointCorrupted {
                detail: format!(
                    "unknown snapshot format: expected 'cmc-rs-snapshot-v2', got '{format}'"
                ),
            });
        }

        let values = snapshot["spins"]
            .as_array()
            .ok_or_else(|| invalid_checkpoint("snapshot.spins", "missing or invalid array"))?;
        let mut spins = Vec::with_capacity(values.len());
        for value in values {
            let spin = value
                .as_f64()
                .ok_or_else(|| invalid_checkpoint("snapshot.spins", "contains a non-number"))?;
            if !spin.is_finite() {
                return Err(invalid_checkpoint(
                    "snapshot.spins",
                    "contains a non-finite value",
                ));
            }
            spins.push(spin);
        }

        if spins.len() != self.system.spins.len() {
            return Err(invalid_checkpoint(
                "snapshot.spins",
                format!(
                    "length mismatch: expected {}, got {}",
                    self.system.spins.len(),
                    spins.len()
                ),
            ));
        }
        if let Some(n_sites) = snapshot["n_sites"].as_u64() {
            if n_sites as usize != self.system.n_sites() {
                return Err(invalid_checkpoint("snapshot.n_sites", "topology mismatch"));
            }
        }
        if let Some(spin_dim) = snapshot["spin_dim"].as_u64() {
            if spin_dim as usize != self.model.spin_dim() {
                return Err(invalid_checkpoint("snapshot.spin_dim", "model mismatch"));
            }
        }
        if let Some(n_edges) = snapshot["n_edges"].as_u64() {
            if n_edges as usize != self.system.lattice.n_edges() {
                return Err(invalid_checkpoint("snapshot.n_edges", "topology mismatch"));
            }
        }
        validate_snapshot_usize_array(snapshot, "offsets", &self.system.lattice.offsets)?;
        validate_snapshot_usize_array(snapshot, "neighbors", &self.system.lattice.neighbors)?;
        validate_snapshot_usize_array(snapshot, "edge_ids", &self.system.lattice.edge_ids)?;
        validate_snapshot_edges(snapshot, &self.system.lattice)?;

        let beta = snapshot["beta"]
            .as_f64()
            .ok_or_else(|| invalid_checkpoint("snapshot.beta", "missing or invalid"))?;
        self.system
            .set_beta(beta)
            .map_err(|reason| invalid_checkpoint("snapshot.beta", reason))?;
        self.system.spins.copy_from_slice(&spins);
        self.system.recompute_energy(&self.model);
        Ok(())
    }
}

impl<H, A, O> MonteCarlo for ClassicalMC<H, A, O>
where
    H: Hamiltonian,
    A: Algorithm<H>,
    O: ObservableSet<H>,
{
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, context: &mut Context<Self::Rng>) {
        let phase = SimulationPhase::from_run_phase(context.phase());
        self.algorithm
            .sweep_with_phase(&mut self.system, &self.model, &mut context.rng, phase);
    }

    fn measure(&mut self, context: &mut Context<Self::Rng>) {
        self.observables
            .measure_all(&self.system, &self.model, context);
    }

    fn name(&self) -> &'static str {
        self.algorithm.name()
    }
}

impl<H, A, O> ParallelTemperingCompatible for ClassicalMC<H, A, O>
where
    H: Hamiltonian,
    A: Algorithm<H>,
    O: ObservableSet<H>,
{
    fn log_weight_ratio(&self, parameter: &str, new_value: f64) -> f64 {
        match parameter {
            "beta" => (self.system.beta - new_value) * self.system.energy,
            _ => panic!("unsupported classical PT parameter: {parameter}"),
        }
    }

    fn change_parameter(&mut self, parameter: &str, new_value: f64) {
        match parameter {
            "beta" => {
                self.system
                    .set_beta(new_value)
                    .expect("parallel-tempering beta must be finite and non-negative");
                // Physical energy is beta-independent, but exact recomputation
                // also repairs any accumulated cache drift before an exchange.
                self.system.recompute_energy(&self.model);
            }
            _ => panic!("unsupported classical PT parameter: {parameter}"),
        }
    }
}

fn validate_snapshot_usize_array(
    snapshot: &Json,
    name: &str,
    expected: &[usize],
) -> Result<(), CarloError> {
    let Some(values) = snapshot[name].as_array() else {
        // Version-1 snapshots did not carry complete topology metadata.
        return Ok(());
    };
    if values.len() != expected.len() {
        return Err(invalid_checkpoint(
            format!("snapshot.{name}"),
            "topology length mismatch",
        ));
    }
    for (index, (value, expected_value)) in values.iter().zip(expected).enumerate() {
        if value.as_u64().map(|number| number as usize) != Some(*expected_value) {
            return Err(invalid_checkpoint(
                format!("snapshot.{name}[{index}]"),
                "topology mismatch",
            ));
        }
    }
    Ok(())
}

fn validate_snapshot_edges(snapshot: &Json, lattice: &CsrLattice) -> Result<(), CarloError> {
    let Some(values) = snapshot["edges"].as_array() else {
        return Ok(());
    };
    if values.len() != lattice.edges.len() {
        return Err(invalid_checkpoint(
            "snapshot.edges",
            "topology length mismatch",
        ));
    }
    for (index, (value, expected)) in values.iter().zip(&lattice.edges).enumerate() {
        let source = value["source"].as_u64().map(|number| number as usize);
        let target = value["target"].as_u64().map(|number| number as usize);
        let kind = value["kind"].as_str().and_then(BondType::from_label);
        let weight = value["weight"].as_f64();
        if source != Some(expected.source)
            || target != Some(expected.target)
            || kind != Some(expected.kind)
            || weight.map(f64::to_bits) != Some(expected.weight.to_bits())
        {
            return Err(invalid_checkpoint(
                format!("snapshot.edges[{index}]"),
                "physical edge mismatch",
            ));
        }
    }
    Ok(())
}

fn invalid_checkpoint(field: impl Into<String>, reason: impl Into<String>) -> CarloError {
    CarloError::CheckpointCorrupted {
        detail: format!("{}: {}", field.into(), reason.into()),
    }
}

fn invalid(field: impl Into<String>, reason: impl Into<String>) -> CarloError {
    CarloError::InvalidConfig {
        field: field.into(),
        reason: reason.into(),
    }
}

pub(crate) fn parse_param<T>(params: &Params, name: &str) -> Result<Option<T>, CarloError>
where
    T: FromStr,
    T::Err: Display,
{
    if !params.contains(name) {
        return Ok(None);
    }
    let raw = params
        .get::<String>(name)
        .expect("Params values are stored as valid strings");
    raw.parse::<T>()
        .map(Some)
        .map_err(|error| invalid(name, format!("cannot parse `{raw}`: {error}")))
}

pub(crate) fn parse_bool(params: &Params, name: &str, default: bool) -> Result<bool, CarloError> {
    let Some(value) = params.get::<String>(name) else {
        return Ok(default);
    };
    match value.to_ascii_lowercase().as_str() {
        "true" | "1" | "yes" | "on" => Ok(true),
        "false" | "0" | "no" | "off" => Ok(false),
        _ => Err(invalid(name, format!("expected boolean, got `{value}`"))),
    }
}

fn positive_dimension(params: &Params, name: &str, default: usize) -> Result<usize, CarloError> {
    let value = parse_param::<usize>(params, name)?.unwrap_or(default);
    if value == 0 {
        Err(invalid(name, "dimension must be positive"))
    } else {
        Ok(value)
    }
}

/// Build a CMC arbitrary graph from standard lattice parameters.
pub fn build_lattice_from_params(params: &Params, pbc: bool) -> Result<CsrLattice, CarloError> {
    let inferred = if params.contains("Lx") {
        "hypercubic"
    } else {
        "chain"
    };
    let lattice_type = params
        .get::<String>("lattice_type")
        .unwrap_or_else(|| inferred.to_string())
        .to_ascii_lowercase();

    match lattice_type.as_str() {
        "chain" => {
            let length = if params.contains("L") {
                positive_dimension(params, "L", 10)?
            } else {
                positive_dimension(params, "Lx", 10)?
            };
            Ok(build_chain(length, pbc))
        }
        "square" => {
            let lx = positive_dimension(params, "Lx", 4)?;
            let ly = positive_dimension(params, "Ly", lx)?;
            Ok(build_square(lx, ly, pbc))
        }
        "cubic" => {
            let lx = positive_dimension(params, "Lx", 4)?;
            let ly = positive_dimension(params, "Ly", lx)?;
            let lz = positive_dimension(params, "Lz", lx)?;
            Ok(build_hypercubic(
                &[lx, ly, lz],
                &[BondType::CubicX, BondType::CubicY, BondType::CubicZ],
                pbc,
            ))
        }
        "hypercubic" => {
            if params.contains("Lx") {
                let lx = positive_dimension(params, "Lx", 4)?;
                let ly = positive_dimension(params, "Ly", lx)?;
                if params.contains("Lz") {
                    let lz = positive_dimension(params, "Lz", lx)?;
                    Ok(build_hypercubic(
                        &[lx, ly, lz],
                        &[BondType::CubicX, BondType::CubicY, BondType::CubicZ],
                        pbc,
                    ))
                } else {
                    Ok(build_square(lx, ly, pbc))
                }
            } else {
                Ok(build_chain(positive_dimension(params, "L", 10)?, pbc))
            }
        }
        "triangular" | "honeycomb" | "kagome" => {
            if !pbc {
                return Err(invalid(
                    "pbc",
                    format!("lattice_type `{lattice_type}` currently requires periodic boundaries"),
                ));
            }
            let default = if lattice_type == "kagome" { 2 } else { 4 };
            let lx = positive_dimension(params, "Lx", default)?;
            let ly = positive_dimension(params, "Ly", lx)?;
            match lattice_type.as_str() {
                "triangular" if lx < 2 || ly < 2 => Err(invalid(
                    "Lx/Ly",
                    "triangular lattice requires both dimensions >= 2",
                )),
                "honeycomb" if lx < 2 || ly < 2 || !lx.is_multiple_of(2) => Err(invalid(
                    "Lx/Ly",
                    "honeycomb lattice requires even Lx >= 2 and Ly >= 2",
                )),
                "kagome" if lx < 2 || ly < 2 => Err(invalid(
                    "Lx/Ly",
                    "kagome lattice requires both dimensions >= 2",
                )),
                "triangular" => Ok(build_triangular(lx, ly)),
                "honeycomb" => Ok(build_honeycomb(lx, ly)),
                "kagome" => Ok(build_kagome(lx, ly)),
                _ => unreachable!(),
            }
        }
        _ => Err(invalid(
            "lattice_type",
            format!(
                "unknown lattice `{lattice_type}`; expected chain, square, cubic, hypercubic, triangular, honeycomb, or kagome"
            ),
        )),
    }
}

impl<H, A> FromParams for ClassicalMC<H, A, DefaultObservableSet<H>>
where
    H: Hamiltonian + Initializable + Measurable + FromHamiltonianParams,
    A: Algorithm<H> + Default,
{
    fn validate_params(params: &Params) -> Result<(), CarloError> {
        let beta = parse_param::<f64>(params, "beta")?.unwrap_or(1.0);
        if !beta.is_finite() || beta < 0.0 {
            return Err(invalid("beta", "must be finite and non-negative"));
        }
        let pbc = parse_bool(params, "pbc", true)?;
        build_lattice_from_params(params, pbc)?;
        H::from_hamiltonian_params(params)?;
        Ok(())
    }

    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Self::validate_params(params)?;
        let pbc = parse_bool(params, "pbc", true)?;
        let lattice = build_lattice_from_params(params, pbc)?;
        let beta = parse_param::<f64>(params, "beta")?.unwrap_or(1.0);
        let model = H::from_hamiltonian_params(params)?;
        let spin_dim = model.spin_dim();
        let mut system = System::new(lattice, spin_dim, 0.0, beta);

        let initial_state = params
            .get::<String>("initial_state")
            .unwrap_or_else(|| "hot".to_string())
            .to_ascii_lowercase();
        for site in 0..system.n_sites() {
            let spin = match initial_state.as_str() {
                "hot" | "random" => model.random_spin(rng),
                "cold" | "ordered" => model.ordered_spin(),
                _ => {
                    return Err(invalid(
                        "initial_state",
                        "expected hot/random or cold/ordered",
                    ))
                }
            };
            if spin.len() != spin_dim {
                return Err(invalid(
                    "initial_state",
                    format!(
                        "model initializer returned {} components, expected {spin_dim}",
                        spin.len()
                    ),
                ));
            }
            system.spin_at_mut(site, spin_dim).copy_from_slice(&spin);
        }
        system.recompute_energy(&model);
        system
            .validate(&model)
            .map_err(|reason| invalid("system", reason))?;

        Ok(Self {
            system,
            model,
            algorithm: A::default(),
            observables: DefaultObservableSet::new(),
        })
    }
}

/// Parse model parameters independently of topology and Carlo scheduling.
pub trait FromHamiltonianParams: Hamiltonian + Sized {
    fn from_hamiltonian_params(params: &Params) -> Result<Self, CarloError>;
}

fn finite_coupling(params: &Params) -> Result<f64, CarloError> {
    let coupling = parse_param::<f64>(params, "J")?.unwrap_or(1.0);
    if coupling.is_finite() {
        Ok(coupling)
    } else {
        Err(invalid("J", "must be finite"))
    }
}

impl FromHamiltonianParams for crate::lattice::models::IsingModel {
    fn from_hamiltonian_params(params: &Params) -> Result<Self, CarloError> {
        Ok(Self::new(finite_coupling(params)?))
    }
}

impl FromHamiltonianParams for crate::lattice::models::PottsModel {
    fn from_hamiltonian_params(params: &Params) -> Result<Self, CarloError> {
        let q = parse_param::<usize>(params, "q")?.unwrap_or(3);
        if q < 2 {
            return Err(invalid("q", "Potts q must be >= 2"));
        }
        Ok(Self::new(finite_coupling(params)?, q))
    }
}

impl<const D: usize> FromHamiltonianParams for crate::lattice::models::ONModel<D> {
    fn from_hamiltonian_params(params: &Params) -> Result<Self, CarloError> {
        if D < 2 {
            return Err(invalid("spin_dim", "O(N) dimension must be >= 2"));
        }
        Ok(Self::new(finite_coupling(params)?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{IsingModel, MetropolisCore};
    use carlo_rs::{ParallelTemperingCompatible, RayonBackend, RunConfig, Scheduler};
    use rand::SeedableRng;

    #[test]
    fn carlo_scheduler_end_to_end() {
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 1.0);
        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 100,
                measurement_sweeps: 200,
                binsize: 50,
                base_seed: 42,
                ..Default::default()
            },
        )
        .run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);
        assert!(results.get("Energy").is_some());
        assert!(results.get("E2").is_some());
        assert!(results.get("M4").is_some());
    }

    #[test]
    fn snapshot_recomputes_energy() {
        let mut params = Params::new();
        params.set("L", 4usize);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
        let mut mc =
            <ClassicalMC<IsingModel, MetropolisCore> as FromParams>::from_params(&params, &mut rng)
                .unwrap();
        let snapshot = mc.save_snapshot();
        mc.system.energy = 123.0;
        mc.load_snapshot(&snapshot).unwrap();
        assert!(mc.system.energy_error(&mc.model).abs() < 1e-12);
    }

    #[test]
    fn snapshot_rejects_topology_mismatch() {
        let mut params = Params::new();
        params.set("L", 4usize);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(11);
        let mut mc =
            <ClassicalMC<IsingModel, MetropolisCore> as FromParams>::from_params(&params, &mut rng)
                .unwrap();
        let mut snapshot = mc.save_snapshot();
        snapshot["neighbors"][0] = serde_json::json!(999usize);
        assert!(mc.load_snapshot(&snapshot).is_err());
    }

    #[test]
    fn malformed_numeric_parameter_is_rejected() {
        let mut params = Params::new();
        params.set("L", "not-a-number");
        let error = build_lattice_from_params(&params, true).unwrap_err();
        assert!(matches!(error, CarloError::InvalidConfig { .. }));
    }

    #[test]
    fn unknown_lattice_is_rejected() {
        let mut params = Params::new();
        params.set("lattice_type", "unknown-geometry");
        let error = build_lattice_from_params(&params, true).unwrap_err();
        assert!(matches!(error, CarloError::InvalidConfig { .. }));
    }

    #[test]
    fn pt_uses_physical_energy() {
        let mut params = Params::new();
        params.set("L", 4usize);
        params.set("beta", 1.0);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2);
        let mc =
            <ClassicalMC<IsingModel, MetropolisCore> as FromParams>::from_params(&params, &mut rng)
                .unwrap();
        assert_eq!(
            mc.log_weight_ratio("beta", 2.0),
            (mc.system.beta - 2.0) * mc.system.energy
        );
    }
}
