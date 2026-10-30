//! Contact process on a static graph, driven by the direct Gillespie kernel.

use core::fmt;

use rand::{Rng, RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

use crate::dynamics::{DynamicsError, GillespieKernel, RejectionFreeModel};
use crate::{CsrLattice, GraphView, UndirectedGraphView};

/// One contact-process state: one activity bit per dense vertex.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ContactState {
    active: Vec<bool>,
}

impl ContactState {
    /// All-inactive state (the absorbing state) over `vertices` vertices.
    pub fn new(vertices: usize) -> Self {
        Self {
            active: vec![false; vertices],
        }
    }

    /// State from an explicit activity mask.
    pub fn from_mask(mask: &[bool]) -> Self {
        Self {
            active: mask.to_vec(),
        }
    }

    pub fn active(&self) -> &[bool] {
        &self.active
    }

    pub fn active_count(&self) -> usize {
        self.active.iter().filter(|&&active| active).count()
    }

    pub fn is_active(&self, vertex: usize) -> bool {
        self.active.get(vertex).copied().unwrap_or(false)
    }

    /// Set one vertex active or inactive; panics outside the vertex domain.
    pub fn set_active(&mut self, vertex: usize, active: bool) {
        self.active[vertex] = active;
    }

    pub fn vertices(&self) -> usize {
        self.active.len()
    }
}

/// Prepared patch for one selected event.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ContactPatch {
    #[default]
    Inert,
    Recovery(usize),
    Infection(usize),
}

/// Contact-process rates over a fixed graph.
///
/// Event catalog: vertices `0..V` are recovery events, edges `V..V+E` are
/// transmission events. A recovery event has rate `mu` iff the vertex is
/// active; a transmission event has rate `lambda` iff exactly one endpoint is
/// active, and it activates the inactive endpoint. Events with an
/// inapplicable state have rate zero and are never selected.
#[derive(Clone, Debug)]
pub struct ContactProcessModel {
    graph: CsrLattice,
    infection_rate: f64,
    recovery_rate: f64,
}

impl ContactProcessModel {
    pub fn new(
        graph: CsrLattice,
        infection_rate: f64,
        recovery_rate: f64,
    ) -> Result<Self, ContactError> {
        graph
            .validate()
            .map_err(|error| ContactError::Topology(error.to_string()))?;
        if graph.vertex_count() == 0 {
            return Err(ContactError::EmptyGraph);
        }
        validate_rate("infection", infection_rate)?;
        validate_rate("recovery", recovery_rate)?;
        Ok(Self {
            graph,
            infection_rate,
            recovery_rate,
        })
    }

    pub const fn graph(&self) -> &CsrLattice {
        &self.graph
    }

    pub const fn infection_rate(&self) -> f64 {
        self.infection_rate
    }

    pub const fn recovery_rate(&self) -> f64 {
        self.recovery_rate
    }

    fn edge_infection(&self, state: &ContactState, edge: usize) -> Option<usize> {
        let edge = self.graph.edge_id(edge).expect("dense edge index");
        let [left, right] = self.graph.edge_endpoints(edge);
        match (
            state.is_active(left.index()),
            state.is_active(right.index()),
        ) {
            (true, false) => Some(right.index()),
            (false, true) => Some(left.index()),
            _ => None,
        }
    }
}

fn validate_rate(name: &'static str, value: f64) -> Result<(), ContactError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(ContactError::Rate { name, value })
    }
}

impl RejectionFreeModel for ContactProcessModel {
    type State = ContactState;
    type Patch = ContactPatch;

    fn event_count(&self, _state: &Self::State) -> usize {
        self.graph.vertex_count() + self.graph.edge_count()
    }

    fn event_rate(&self, state: &Self::State, event: usize) -> Result<f64, DynamicsError> {
        let vertices = self.graph.vertex_count();
        if event < vertices {
            Ok(if state.is_active(event) {
                self.recovery_rate
            } else {
                0.0
            })
        } else {
            let edge = event - vertices;
            Ok(if self.edge_infection(state, edge).is_some() {
                self.infection_rate
            } else {
                0.0
            })
        }
    }

    fn prepare_event(
        &self,
        state: &Self::State,
        event: usize,
        patch: &mut Self::Patch,
    ) -> Result<(), DynamicsError> {
        let vertices = self.graph.vertex_count();
        *patch = if event < vertices {
            ContactPatch::Recovery(event)
        } else {
            let edge = event - vertices;
            self.edge_infection(state, edge)
                .map(ContactPatch::Infection)
                .unwrap_or_default()
        };
        Ok(())
    }

    fn commit_event(&self, state: &mut Self::State, _event: usize, patch: &Self::Patch) {
        match *patch {
            ContactPatch::Inert => {}
            ContactPatch::Recovery(vertex) => state.set_active(vertex, false),
            ContactPatch::Infection(vertex) => state.set_active(vertex, true),
        }
    }

    fn validate_state(&self, state: &Self::State) -> Result<(), DynamicsError> {
        if state.vertices() == self.graph.vertex_count() {
            Ok(())
        } else {
            Err(DynamicsError::new(format!(
                "contact state covers {} vertices; graph has {}",
                state.vertices(),
                self.graph.vertex_count()
            )))
        }
    }
}

/// One finite-time contact-process run from an explicit initial state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactRun {
    /// Whether any vertex is active at the end.
    pub survived: bool,
    /// Active vertices at the end.
    pub final_active: usize,
    /// `final_active / V`.
    pub final_density: f64,
    /// Accepted Gillespie events.
    pub events: u64,
    /// Physical time reached (equal to the requested duration unless the
    /// process was absorbed earlier).
    pub time: f64,
}

/// Advance one realization from `initial` through `duration` of event time.
///
/// The direct method recomputes the full event catalog before every event;
/// this is the reference order, not a scaled-time approximation.
pub fn contact_process_run(
    model: &ContactProcessModel,
    initial: &ContactState,
    duration: f64,
    rng: &mut impl Rng,
) -> Result<ContactRun, ContactError> {
    if !duration.is_finite() || duration <= 0.0 {
        return Err(ContactError::Duration { value: duration });
    }
    let vertices = model.graph.vertex_count();
    let mut kernel =
        GillespieKernel::new(model.clone(), initial.clone()).map_err(ContactError::Dynamics)?;
    kernel
        .advance_by(duration, rng)
        .map_err(ContactError::Dynamics)?;
    let state = kernel.state();
    let final_active = state.active_count();
    Ok(ContactRun {
        survived: final_active > 0,
        final_active,
        final_density: final_active as f64 / vertices as f64,
        events: kernel.events(),
        time: kernel.event_time(),
    })
}

/// Finite-time ensemble over independent initial-condition draws.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactEnsemble {
    pub replicas: usize,
    /// Replicas still active at the end.
    pub survived: usize,
    /// `survived / replicas`.
    pub survival_fraction: f64,
    /// Mean final density over all replicas (absorbed runs contribute zero).
    pub mean_density: f64,
}

/// Run `replicas` independent realizations, each starting from a fresh
/// per-vertex density draw (one `f64` per vertex in ascending order; active
/// iff `draw < density`). Deterministic given the seed.
pub fn contact_ensemble(
    model: &ContactProcessModel,
    initial_density: f64,
    replicas: usize,
    duration: f64,
    seed: u64,
) -> Result<ContactEnsemble, ContactError> {
    if !initial_density.is_finite() || !(0.0..=1.0).contains(&initial_density) {
        return Err(ContactError::Density {
            value: initial_density,
        });
    }
    if replicas == 0 {
        return Err(ContactError::EmptyEnsemble);
    }
    if !duration.is_finite() || duration <= 0.0 {
        return Err(ContactError::Duration { value: duration });
    }
    let vertices = model.graph.vertex_count();
    let mut survived = 0usize;
    let mut density_sum = 0.0;
    for replica in 0..replicas {
        let replica_seed = seed.wrapping_add((replica as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(replica_seed);
        let mut state = ContactState::new(vertices);
        for vertex in 0..vertices {
            state.set_active(vertex, rng.random::<f64>() < initial_density);
        }
        let run = contact_process_run(model, &state, duration, &mut rng)?;
        survived += usize::from(run.survived);
        density_sum += run.final_density;
    }
    Ok(ContactEnsemble {
        replicas,
        survived,
        survival_fraction: survived as f64 / replicas as f64,
        mean_density: density_sum / replicas as f64,
    })
}

/// Typed rejections for the temporal/X4 surface.
#[derive(Debug)]
pub enum ContactError {
    /// The graph failed topology validation.
    Topology(String),
    /// No vertices to host the process.
    EmptyGraph,
    /// A rate is not finite and strictly positive.
    Rate { name: &'static str, value: f64 },
    /// The requested duration is not finite and positive.
    Duration { value: f64 },
    /// The initial density is outside `[0, 1]`.
    Density { value: f64 },
    /// An ensemble needs at least one replica.
    EmptyEnsemble,
    /// The Gillespie kernel rejected the state or an event.
    Dynamics(DynamicsError),
}

impl fmt::Display for ContactError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Topology(reason) => {
                write!(formatter, "invalid contact-process topology: {reason}")
            }
            Self::EmptyGraph => formatter.write_str("contact process needs at least one vertex"),
            Self::Rate { name, value } => {
                write!(
                    formatter,
                    "{name} rate must be finite and positive, got {value}"
                )
            }
            Self::Duration { value } => {
                write!(
                    formatter,
                    "duration must be finite and positive, got {value}"
                )
            }
            Self::Density { value } => {
                write!(
                    formatter,
                    "initial density must be finite within [0, 1], got {value}"
                )
            }
            Self::EmptyEnsemble => formatter.write_str("ensemble needs at least one replica"),
            Self::Dynamics(source) => write!(formatter, "Gillespie kernel failure: {source}"),
        }
    }
}

impl std::error::Error for ContactError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Dynamics(source) => Some(source),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Bond, BondType};

    #[test]
    fn absorbing_state_has_zero_rates() {
        let graph = CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);
        let model = ContactProcessModel::new(graph, 0.7, 1.0).expect("model");
        let state = ContactState::new(2);
        assert_eq!(
            RejectionFreeModel::event_count(&model, &state),
            3,
            "two recovery + one transmission"
        );
        for event in 0..3 {
            assert_eq!(
                RejectionFreeModel::event_rate(&model, &state, event).expect("rate"),
                0.0
            );
        }
    }

    #[test]
    fn hand_event_rates_follow_activity() {
        let graph = CsrLattice::from_edges(
            3,
            vec![
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(1, 2, BondType::Generic, 1.0),
            ],
        );
        let model = ContactProcessModel::new(graph, 0.5, 2.0).expect("model");
        let state = ContactState::from_mask(&[true, false, false]);
        // Vertices 0..2 recovery, edges 2..4 transmission.
        let rates: Vec<f64> = (0..5)
            .map(|event| RejectionFreeModel::event_rate(&model, &state, event).expect("rate"))
            .collect();
        assert_eq!(rates, vec![2.0, 0.0, 0.0, 0.5, 0.0]);
    }

    #[test]
    fn deterministic_run_replays_bitwise() {
        let graph = CsrLattice::from_edges(
            4,
            vec![
                Bond::new(0, 1, BondType::Generic, 1.0),
                Bond::new(1, 2, BondType::Generic, 1.0),
                Bond::new(2, 3, BondType::Generic, 1.0),
            ],
        );
        let model = ContactProcessModel::new(graph, 2.0, 1.0).expect("model");
        let initial = ContactState::from_mask(&[true, true, false, false]);
        let mut first = Xoshiro256PlusPlus::seed_from_u64(7);
        let mut second = Xoshiro256PlusPlus::seed_from_u64(7);
        let left = contact_process_run(&model, &initial, 5.0, &mut first).expect("run");
        let right = contact_process_run(&model, &initial, 5.0, &mut second).expect("run");
        assert_eq!(left, right);
    }

    #[test]
    fn invalid_domains_are_typed() {
        let graph = || CsrLattice::from_edges(2, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);
        assert!(matches!(
            ContactProcessModel::new(graph(), 0.0, 1.0),
            Err(ContactError::Rate {
                name: "infection",
                value: 0.0
            })
        ));
        let model = ContactProcessModel::new(graph(), 0.5, 1.0).expect("model");
        assert!(matches!(
            contact_process_run(&model, &ContactState::new(3), 1.0, &mut rand::rng()),
            Err(ContactError::Dynamics(_))
        ));
        assert!(matches!(
            contact_ensemble(&model, 0.5, 0, 1.0, 1),
            Err(ContactError::EmptyEnsemble)
        ));
        assert!(matches!(
            contact_ensemble(&model, 1.5, 1, 1.0, 1),
            Err(ContactError::Density { value: 1.5 })
        ));
    }
}
