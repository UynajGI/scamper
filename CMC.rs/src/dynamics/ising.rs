//! Conserved and continuous-time Ising dynamics.

use super::{DynamicsError, KineticRateLaw, RejectionFreeModel};
use crate::algorithms::{Algorithm, SimulationPhase};
use crate::audit::{audit_lattice_cache, should_audit_cache};
use crate::core::acceptance::MetropolisHastingsAcceptance;
use crate::core::cache::{BatchEnergyPatch, EnergyPatch};
use crate::core::r#move::BatchSpinMove;
use crate::core::trial::{metropolis_hastings_step, ProposedMove};
use crate::lattice::interaction::Hamiltonian;
use crate::lattice::models::IsingModel;
use crate::lattice::state::System;
use rand::{Rng, RngExt};
use serde_json::{json, Value as Json};

/// Reference direct-Gillespie Ising model: one possible flip event per site.
#[derive(Debug, Clone)]
pub struct KineticIsingModel {
    hamiltonian: IsingModel,
    rate_law: KineticRateLaw,
}

impl KineticIsingModel {
    pub fn new(coupling: f64, rate_law: KineticRateLaw) -> Result<Self, DynamicsError> {
        if !coupling.is_finite() {
            return Err(DynamicsError::new("kinetic Ising coupling must be finite"));
        }
        rate_law.validate()?;
        Ok(Self {
            hamiltonian: IsingModel::new(coupling),
            rate_law,
        })
    }

    #[inline]
    pub const fn hamiltonian(&self) -> &IsingModel {
        &self.hamiltonian
    }

    #[inline]
    pub const fn rate_law(&self) -> KineticRateLaw {
        self.rate_law
    }

    pub fn flip_delta_energy(&self, state: &System, site: usize) -> Result<f64, DynamicsError> {
        if site >= state.n_sites() {
            return Err(DynamicsError::new("kinetic Ising site out of range"));
        }
        let proposed = [-state.spins[site]];
        Ok(self
            .hamiltonian
            .delta_energy(&state.spins, &state.lattice, site, &proposed))
    }

    pub fn flip_rate(&self, state: &System, site: usize) -> Result<f64, DynamicsError> {
        self.rate_law
            .rate(state.beta, self.flip_delta_energy(state, site)?)
    }
}

impl RejectionFreeModel for KineticIsingModel {
    type State = System;
    type Patch = EnergyPatch;

    fn event_count(&self, state: &Self::State) -> usize {
        state.n_sites()
    }

    fn event_rate(&self, state: &Self::State, event: usize) -> Result<f64, DynamicsError> {
        self.flip_rate(state, event)
    }

    fn prepare_event(
        &self,
        state: &Self::State,
        event: usize,
        patch: &mut Self::Patch,
    ) -> Result<(), DynamicsError> {
        patch.delta_energy = self.flip_delta_energy(state, event)?;
        Ok(())
    }

    fn commit_event(&self, state: &mut Self::State, event: usize, patch: &Self::Patch) {
        state.spins[event] = -state.spins[event];
        state.energy += patch.delta_energy;
    }

    fn validate_state(&self, state: &Self::State) -> Result<(), DynamicsError> {
        state
            .validate(&self.hamiltonian)
            .map_err(DynamicsError::new)?;
        if state
            .spins
            .iter()
            .any(|spin| (*spin - 1.0).abs() > 1e-12 && (*spin + 1.0).abs() > 1e-12)
        {
            return Err(DynamicsError::new(
                "kinetic Ising state contains a non-Ising spin",
            ));
        }
        audit_lattice_cache(state, &self.hamiltonian).map_err(DynamicsError::new)
    }
}

/// One rejection-free BKL/n-fold-way spin-flip event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BklEvent {
    pub site: usize,
    pub delta_energy: f64,
    pub delta_time: f64,
    pub total_rate: f64,
}

/// Fenwick-tree BKL kernel for arbitrary weighted Ising graphs.
///
/// Traditional class tables are optimal on uniform lattices.  A rate Fenwick
/// tree is the equivalent rejection-free selection structure for arbitrary
/// weighted/multigraph lattices, with `O(log N)` selection and local updates.
pub struct BklIsingKernel {
    model: KineticIsingModel,
    state: System,
    rates: Vec<f64>,
    tree: FenwickRates,
    event_time: f64,
    events: u64,
    audit_interval: u64,
}

impl BklIsingKernel {
    pub fn new(
        model: KineticIsingModel,
        mut state: System,
        audit_interval: u64,
    ) -> Result<Self, DynamicsError> {
        state.recompute_energy(model.hamiltonian());
        model.validate_state(&state)?;
        let rates = (0..state.n_sites())
            .map(|site| model.flip_rate(&state, site))
            .collect::<Result<Vec<_>, _>>()?;
        let tree = FenwickRates::new(&rates)?;
        Ok(Self {
            model,
            state,
            rates,
            tree,
            event_time: 0.0,
            events: 0,
            audit_interval,
        })
    }

    #[inline]
    pub const fn model(&self) -> &KineticIsingModel {
        &self.model
    }

    #[inline]
    pub const fn state(&self) -> &System {
        &self.state
    }

    #[inline]
    pub const fn event_time(&self) -> f64 {
        self.event_time
    }

    #[inline]
    pub const fn events(&self) -> u64 {
        self.events
    }

    #[inline]
    pub fn total_rate(&self) -> f64 {
        self.tree.total()
    }

    pub fn validate(&self) -> Result<(), DynamicsError> {
        self.model.validate_state(&self.state)?;
        if self.rates.len() != self.state.n_sites() {
            return Err(DynamicsError::new("BKL rate-vector length mismatch"));
        }
        let exact = (0..self.state.n_sites())
            .map(|site| self.model.flip_rate(&self.state, site))
            .collect::<Result<Vec<_>, _>>()?;
        for (site, (&cached, &expected)) in self.rates.iter().zip(&exact).enumerate() {
            let tolerance = 1e-12 * (1.0 + expected.abs());
            if (cached - expected).abs() > tolerance {
                return Err(DynamicsError::new(format!(
                    "BKL cached rate mismatch at site {site}"
                )));
            }
        }
        self.tree.validate(&self.rates)?;
        if !self.event_time.is_finite() || self.event_time < 0.0 {
            return Err(DynamicsError::new("BKL event-time clock is invalid"));
        }
        Ok(())
    }

    /// Execute one event.  `Ok(None)` denotes an absorbing state.
    pub fn step(&mut self, rng: &mut impl Rng) -> Result<Option<BklEvent>, DynamicsError> {
        let total_rate = self.total_rate();
        if total_rate == 0.0 {
            return Ok(None);
        }
        if !total_rate.is_finite() || total_rate < 0.0 {
            return Err(DynamicsError::new("BKL total rate is invalid"));
        }
        let threshold = rng.random::<f64>() * total_rate;
        let site = self.tree.select(threshold)?;
        let delta_energy = self.model.flip_delta_energy(&self.state, site)?;
        let delta_time = super::gillespie::exponential_wait(total_rate, rng);
        self.state.spins[site] = -self.state.spins[site];
        self.state.energy += delta_energy;
        self.event_time += delta_time;
        self.events = self.events.saturating_add(1);

        self.refresh_local_rates(site)?;

        if should_audit_cache(self.events, self.audit_interval) {
            self.validate()?;
        }
        Ok(Some(BklEvent {
            site,
            delta_energy,
            delta_time,
            total_rate,
        }))
    }

    fn refresh_local_rates(&mut self, site: usize) -> Result<(), DynamicsError> {
        let mut affected = Vec::with_capacity(
            self.state.lattice.offsets[site + 1] - self.state.lattice.offsets[site] + 1,
        );
        affected.push(site);
        affected.extend(
            self.state.lattice.neighbors
                [self.state.lattice.offsets[site]..self.state.lattice.offsets[site + 1]]
                .iter()
                .copied(),
        );
        affected.sort_unstable();
        affected.dedup();
        for affected_site in affected {
            let rate = self.model.flip_rate(&self.state, affected_site)?;
            self.rates[affected_site] = rate;
            self.tree.set(affected_site, rate)?;
        }
        Ok(())
    }

    pub fn advance_by(&mut self, duration: f64, rng: &mut impl Rng) -> Result<u64, DynamicsError> {
        if !duration.is_finite() || duration <= 0.0 {
            return Err(DynamicsError::new(
                "BKL advance duration must be finite and positive",
            ));
        }
        let target = self.event_time + duration;
        let before = self.events;
        while self.event_time < target {
            let total_rate = self.total_rate();
            if total_rate == 0.0 {
                self.event_time = target;
                break;
            }
            if !total_rate.is_finite() || total_rate < 0.0 {
                return Err(DynamicsError::new("BKL total rate is invalid"));
            }
            let site = self.tree.select(rng.random::<f64>() * total_rate)?;
            let delta_time = super::gillespie::exponential_wait(total_rate, rng);
            let remaining = target - self.event_time;
            if delta_time > remaining {
                self.event_time = target;
                break;
            }
            let delta_energy = self.model.flip_delta_energy(&self.state, site)?;
            self.state.spins[site] = -self.state.spins[site];
            self.state.energy += delta_energy;
            self.event_time += delta_time;
            self.events = self.events.saturating_add(1);
            self.refresh_local_rates(site)?;
            if should_audit_cache(self.events, self.audit_interval) {
                self.validate()?;
            }
        }
        Ok(self.events.saturating_sub(before))
    }

    /// Versioned, model-validated JSON snapshot.  RNG state remains owned by Carlo.rs Context.
    pub fn save_snapshot(&self) -> Json {
        json!({
            "format": "cmc-rs-bkl-ising-v1",
            "model": {
                "coupling": self.model.hamiltonian().j,
                "beta": self.state.beta,
                "n_sites": self.state.n_sites(),
                "edges": self.state.lattice.edges.iter().map(|edge| json!({
                    "source": edge.source,
                    "target": edge.target,
                    "kind": edge.kind.as_label(),
                    "weight": edge.weight,
                })).collect::<Vec<_>>(),
                "rate_law": rate_law_json(self.model.rate_law()),
            },
            "state": {
                "spins": self.state.spins,
                "energy": self.state.energy,
                "event_time": self.event_time,
                "events": self.events,
            },
            "rate_cache": {
                "rates": self.rates,
                "fenwick_tree": self.tree.tree,
                "fenwick_values": self.tree.values,
            },
            "audit_interval": self.audit_interval,
        })
    }

    pub fn load_snapshot(&mut self, snapshot: &Json) -> Result<(), DynamicsError> {
        if snapshot["format"].as_str() != Some("cmc-rs-bkl-ising-v1") {
            return Err(DynamicsError::new("unknown BKL Ising snapshot format"));
        }
        let model = &snapshot["model"];
        require_same_f64(model, "coupling", self.model.hamiltonian().j)?;
        require_same_f64(model, "beta", self.state.beta)?;
        if required_usize(model, "n_sites")? != self.state.n_sites() {
            return Err(DynamicsError::new("BKL snapshot site-count mismatch"));
        }
        validate_edges(&self.state, &model["edges"])?;
        if model["rate_law"] != rate_law_json(self.model.rate_law()) {
            return Err(DynamicsError::new("BKL snapshot rate-law mismatch"));
        }
        let state = &snapshot["state"];
        let spins = state["spins"]
            .as_array()
            .ok_or_else(|| DynamicsError::new("BKL snapshot spins must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_f64()
                    .ok_or_else(|| DynamicsError::new("BKL snapshot contains a non-f64 spin"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        if spins.len() != self.state.n_sites() {
            return Err(DynamicsError::new("BKL snapshot spin-count mismatch"));
        }
        self.state.spins = spins;
        self.state.energy = required_f64(state, "energy")?;
        self.event_time = required_f64(state, "event_time")?;
        self.events = required_u64(state, "events")?;
        self.audit_interval = required_u64(snapshot, "audit_interval")?;
        let cache = &snapshot["rate_cache"];
        self.rates = required_f64_vec(cache, "rates")?;
        let fenwick_tree = required_f64_vec(cache, "fenwick_tree")?;
        let fenwick_values = required_f64_vec(cache, "fenwick_values")?;
        self.tree = FenwickRates::from_parts(fenwick_tree, fenwick_values)?;
        self.validate()
    }
}

/// Canonical nearest-bond spin-exchange (Kawasaki) dynamics.
#[derive(Debug, Clone)]
pub struct KawasakiCore {
    attempts_per_sweep: usize,
    patch: BatchEnergyPatch,
    attempts: u64,
    accepts: u64,
    last_attempts: u64,
    last_accepts: u64,
    sweeps: u64,
    audit_interval: u64,
}

impl KawasakiCore {
    pub fn new(attempts_per_sweep: usize) -> Self {
        Self {
            attempts_per_sweep,
            patch: BatchEnergyPatch::default(),
            attempts: 0,
            accepts: 0,
            last_attempts: 0,
            last_accepts: 0,
            sweeps: 0,
            audit_interval: 0,
        }
    }

    pub fn with_cache_audit_interval(mut self, interval: u64) -> Self {
        self.audit_interval = interval;
        self
    }

    #[inline]
    pub const fn attempts(&self) -> u64 {
        self.attempts
    }

    #[inline]
    pub const fn accepts(&self) -> u64 {
        self.accepts
    }

    #[inline]
    pub const fn last_attempts(&self) -> u64 {
        self.last_attempts
    }

    #[inline]
    pub const fn last_accepts(&self) -> u64 {
        self.last_accepts
    }
}

impl Default for KawasakiCore {
    fn default() -> Self {
        Self::new(0)
    }
}

impl Algorithm<IsingModel> for KawasakiCore {
    fn sweep_with_phase(
        &mut self,
        system: &mut System,
        model: &IsingModel,
        rng: &mut impl Rng,
        _phase: SimulationPhase,
    ) {
        let before_attempts = self.attempts;
        let before_accepts = self.accepts;
        let attempts = if self.attempts_per_sweep == 0 {
            system.n_sites()
        } else {
            self.attempts_per_sweep
        };
        let ensemble = system.canonical_ensemble();
        let acceptance = MetropolisHastingsAcceptance;
        for _ in 0..attempts {
            self.attempts = self.attempts.saturating_add(1);
            if system.lattice.edges.is_empty() {
                continue;
            }
            let edge = system.lattice.edges[rng.random_range(0..system.lattice.edges.len())];
            if edge.source == edge.target || system.spins[edge.source] == system.spins[edge.target]
            {
                continue;
            }
            let mut movement = BatchSpinMove::with_capacity(1, 2);
            movement.push(edge.source, &[system.spins[edge.target]]);
            movement.push(edge.target, &[system.spins[edge.source]]);
            let proposal = ProposedMove::symmetric(movement);
            let outcome = metropolis_hastings_step(
                system,
                model,
                &proposal,
                &ensemble,
                &acceptance,
                &mut self.patch,
                rng,
            );
            if outcome.accepted {
                self.accepts = self.accepts.saturating_add(1);
            }
        }
        self.last_attempts = self.attempts.saturating_sub(before_attempts);
        self.last_accepts = self.accepts.saturating_sub(before_accepts);
        self.sweeps = self.sweeps.wrapping_add(1);
        if should_audit_cache(self.sweeps, self.audit_interval) {
            audit_lattice_cache(system, model).expect("Kawasaki lattice cache audit failed");
        }
    }

    fn name(&self) -> &'static str {
        "Kawasaki exchange"
    }
}

#[derive(Debug, Clone)]
struct FenwickRates {
    tree: Vec<f64>,
    values: Vec<f64>,
}

impl FenwickRates {
    fn new(values: &[f64]) -> Result<Self, DynamicsError> {
        let mut result = Self {
            tree: vec![0.0; values.len() + 1],
            values: vec![0.0; values.len()],
        };
        for (index, &value) in values.iter().enumerate() {
            result.set(index, value)?;
        }
        Ok(result)
    }

    fn from_parts(tree: Vec<f64>, values: Vec<f64>) -> Result<Self, DynamicsError> {
        if tree.len() != values.len() + 1 {
            return Err(DynamicsError::new("Fenwick checkpoint length mismatch"));
        }
        if tree
            .iter()
            .chain(&values)
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(DynamicsError::new(
                "Fenwick checkpoint contains an invalid rate",
            ));
        }
        Ok(Self { tree, values })
    }

    fn set(&mut self, index: usize, value: f64) -> Result<(), DynamicsError> {
        if index >= self.values.len() {
            return Err(DynamicsError::new("Fenwick rate index out of range"));
        }
        if !value.is_finite() || value < 0.0 {
            return Err(DynamicsError::new(
                "Fenwick rate must be finite and non-negative",
            ));
        }
        let delta = value - self.values[index];
        self.values[index] = value;
        let mut cursor = index + 1;
        while cursor < self.tree.len() {
            self.tree[cursor] += delta;
            cursor += cursor.isolate_lowest_one();
        }
        Ok(())
    }

    fn total(&self) -> f64 {
        self.prefix_sum(self.values.len())
    }

    fn prefix_sum(&self, end: usize) -> f64 {
        let mut cursor = end;
        let mut sum = 0.0;
        while cursor > 0 {
            sum += self.tree[cursor];
            cursor &= cursor - 1;
        }
        sum
    }

    fn select(&self, threshold: f64) -> Result<usize, DynamicsError> {
        let total = self.total();
        if self.values.is_empty() || !threshold.is_finite() || threshold < 0.0 || threshold >= total
        {
            return Err(DynamicsError::new("Fenwick selection threshold is invalid"));
        }
        let mut index = 0usize;
        let mut accumulated = 0.0;
        let mut bit = 1usize;
        while bit < self.tree.len() {
            bit <<= 1;
        }
        bit >>= 1;
        while bit != 0 {
            let next = index + bit;
            if next < self.tree.len() && accumulated + self.tree[next] <= threshold {
                accumulated += self.tree[next];
                index = next;
            }
            bit >>= 1;
        }
        if index < self.values.len() && self.values[index] > 0.0 {
            Ok(index)
        } else {
            self.values
                .iter()
                .enumerate()
                .skip(index)
                .find_map(|(candidate, rate)| (*rate > 0.0).then_some(candidate))
                .ok_or_else(|| DynamicsError::new("Fenwick selection found no positive rate"))
        }
    }

    fn validate(&self, expected: &[f64]) -> Result<(), DynamicsError> {
        if self.values.len() != expected.len() {
            return Err(DynamicsError::new("Fenwick rate length mismatch"));
        }
        for (index, (&cached, &value)) in self.values.iter().zip(expected).enumerate() {
            let tolerance = 1e-12 * (1.0 + value.abs());
            if (cached - value).abs() > tolerance {
                return Err(DynamicsError::new(format!(
                    "Fenwick value mismatch at event {index}"
                )));
            }
        }
        let expected_total: f64 = expected.iter().sum();
        let tolerance = 1e-11 * (1.0 + expected_total.abs());
        if (self.total() - expected_total).abs() > tolerance {
            return Err(DynamicsError::new("Fenwick total-rate cache mismatch"));
        }
        Ok(())
    }
}

fn rate_law_json(rate_law: KineticRateLaw) -> Json {
    match rate_law {
        KineticRateLaw::Glauber { attempt_frequency } => json!({
            "kind": "glauber",
            "attempt_frequency": attempt_frequency,
        }),
        KineticRateLaw::Metropolis { attempt_frequency } => json!({
            "kind": "metropolis",
            "attempt_frequency": attempt_frequency,
        }),
    }
}

fn validate_edges(state: &System, value: &Json) -> Result<(), DynamicsError> {
    let edges = value
        .as_array()
        .ok_or_else(|| DynamicsError::new("BKL snapshot edges must be an array"))?;
    if edges.len() != state.lattice.edges.len() {
        return Err(DynamicsError::new("BKL snapshot edge-count mismatch"));
    }
    for (saved, edge) in edges.iter().zip(&state.lattice.edges) {
        if required_usize(saved, "source")? != edge.source
            || required_usize(saved, "target")? != edge.target
            || saved["kind"].as_str() != Some(edge.kind.as_label())
            || saved["weight"].as_f64().map(f64::to_bits) != Some(edge.weight.to_bits())
        {
            return Err(DynamicsError::new("BKL snapshot topology mismatch"));
        }
    }
    Ok(())
}

fn required_f64_vec(value: &Json, field: &str) -> Result<Vec<f64>, DynamicsError> {
    value[field]
        .as_array()
        .ok_or_else(|| DynamicsError::new(format!("BKL snapshot field `{field}` is invalid")))?
        .iter()
        .map(|entry| {
            entry
                .as_f64()
                .filter(|number| number.is_finite() && *number >= 0.0)
                .ok_or_else(|| {
                    DynamicsError::new(format!("BKL snapshot field `{field}` is invalid"))
                })
        })
        .collect()
}

fn required_f64(value: &Json, field: &str) -> Result<f64, DynamicsError> {
    value[field]
        .as_f64()
        .filter(|number| number.is_finite())
        .ok_or_else(|| DynamicsError::new(format!("BKL snapshot field `{field}` is invalid")))
}

fn required_u64(value: &Json, field: &str) -> Result<u64, DynamicsError> {
    value[field]
        .as_u64()
        .ok_or_else(|| DynamicsError::new(format!("BKL snapshot field `{field}` is invalid")))
}

fn required_usize(value: &Json, field: &str) -> Result<usize, DynamicsError> {
    usize::try_from(required_u64(value, field)?)
        .map_err(|_| DynamicsError::new(format!("BKL snapshot field `{field}` is too large")))
}

fn require_same_f64(value: &Json, field: &str, expected: f64) -> Result<(), DynamicsError> {
    if value[field].as_f64().map(f64::to_bits) == Some(expected.to_bits()) {
        Ok(())
    } else {
        Err(DynamicsError::new(format!(
            "BKL snapshot field `{field}` mismatch"
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::build_chain;

    #[test]
    fn fenwick_selection_and_updates_match_linear_weights() {
        let mut tree = FenwickRates::new(&[1.0, 2.0, 0.0, 4.0]).unwrap();
        assert_eq!(tree.select(0.5).unwrap(), 0);
        assert_eq!(tree.select(1.5).unwrap(), 1);
        assert_eq!(tree.select(3.5).unwrap(), 3);
        tree.set(1, 0.0).unwrap();
        tree.set(2, 3.0).unwrap();
        tree.validate(&[1.0, 0.0, 3.0, 4.0]).unwrap();
        assert_eq!(tree.select(1.5).unwrap(), 2);
    }

    #[test]
    fn kinetic_rate_ratio_obeys_detailed_balance() {
        let beta = 0.7;
        let delta = 2.3;
        for law in [
            KineticRateLaw::glauber(1.2).unwrap(),
            KineticRateLaw::metropolis(1.2).unwrap(),
        ] {
            let forward = law.rate(beta, delta).unwrap();
            let reverse = law.rate(beta, -delta).unwrap();
            assert!((forward / reverse - (-beta * delta).exp()).abs() < 1e-13);
        }
    }

    #[test]
    fn bkl_constructor_recomputes_energy() {
        let lattice = build_chain(4, true);
        let mut state = System::new(lattice, 1, 1.0, 0.4);
        state.energy = 99.0;
        let model = KineticIsingModel::new(1.0, KineticRateLaw::default()).unwrap();
        let kernel = BklIsingKernel::new(model, state, 0).unwrap();
        assert!((kernel.state().energy + 4.0).abs() < 1e-12);
    }
}
