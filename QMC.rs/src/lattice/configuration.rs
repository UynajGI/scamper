//! Continuous-time interaction-expansion configurations and worldline links.

use rand::Rng;
use rand::RngExt;

use crate::local_space::{BasisState, LocalHilbertSpace, ParticleStatistics};

use super::error::LatticeQmcError;
use super::model::PositiveOperatorModel;
use super::vertex::{Event, Vertex};

/// Sampled continuous-time operator sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct LatticeConfiguration {
    beta: f64,
    initial_states: Vec<BasisState>,
    vertices: Vec<Vertex>,
}

impl LatticeConfiguration {
    /// Construct from explicit states at imaginary time zero.
    pub fn new(
        beta: f64,
        initial_states: Vec<BasisState>,
        model: &impl PositiveOperatorModel,
    ) -> Result<Self, LatticeQmcError> {
        if !beta.is_finite() || beta <= 0.0 {
            return Err(LatticeQmcError::parameter(
                "beta",
                format!("must be finite and positive, got {beta}"),
            ));
        }
        if model.space().statistics() == ParticleStatistics::Fermion {
            return Err(LatticeQmcError::InvalidModel(
                "fermionic local spaces require the reserved signed/determinant backend;                  the positive worldline engine must not silently discard exchange signs"
                    .into(),
            ));
        }
        if initial_states.len() != model.graph().site_count() {
            return Err(LatticeQmcError::InvalidConfiguration(format!(
                "expected {} initial states, got {}",
                model.graph().site_count(),
                initial_states.len()
            )));
        }
        for (site, &state) in initial_states.iter().enumerate() {
            model.space().validate_state(site, state)?;
        }
        Ok(Self {
            beta,
            initial_states,
            vertices: Vec::new(),
        })
    }

    /// Random product state in the local `Sz` basis.
    pub fn random<R: Rng + ?Sized>(
        beta: f64,
        model: &impl PositiveOperatorModel,
        rng: &mut R,
    ) -> Result<Self, LatticeQmcError> {
        let states = (0..model.graph().site_count())
            .map(|site| rng.random_range(0..model.space().dimension(site)) as BasisState)
            .collect();
        Self::new(beta, states, model)
    }

    /// Inverse temperature.
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// States on the time-zero boundary.
    pub fn initial_states(&self) -> &[BasisState] {
        &self.initial_states
    }

    /// Mutable time-zero states. Intended for update kernels.
    pub(crate) fn initial_states_mut(&mut self) -> &mut [BasisState] {
        &mut self.initial_states
    }

    /// Sampled vertices.
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }

    /// Mutable sampled vertices. Intended for update kernels.
    pub(crate) fn vertices_mut(&mut self) -> &mut Vec<Vertex> {
        &mut self.vertices
    }

    /// Expansion order.
    pub fn expansion_order(&self) -> usize {
        self.vertices.len()
    }

    /// Number of diagonal vertices.
    pub fn diagonal_order(&self, model: &impl PositiveOperatorModel) -> usize {
        self.vertices
            .iter()
            .filter(|vertex| model.term(vertex.term).kind(vertex.kind).is_diagonal())
            .count()
    }

    /// Number of off-diagonal vertices.
    pub fn offdiagonal_order(&self, model: &impl PositiveOperatorModel) -> usize {
        self.expansion_order() - self.diagonal_order(model)
    }

    /// Validate every worldline and local matrix element.
    pub fn validate(&self, model: &impl PositiveOperatorModel) -> Result<(), LatticeQmcError> {
        for (site, &state) in self.initial_states.iter().enumerate() {
            model.space().validate_state(site, state)?;
        }
        let index = WorldlineIndex::build(self, model)?;
        for leg in 0..index.leg_count() {
            let linked = index.linked_leg(leg);
            let state = index.state_on_leg(self, model, leg);
            let linked_state = index.state_on_leg(self, model, linked);
            if state != linked_state {
                return Err(LatticeQmcError::InvalidConfiguration(format!(
                    "worldline link {leg}<->{linked} carries states {state} and {linked_state}"
                )));
            }
        }
        for site in 0..model.graph().site_count() {
            if let Some(first) = index.first_event(site) {
                let state = index.state_on_leg(self, model, first.incoming_leg);
                if self.initial_states[site] != state {
                    return Err(LatticeQmcError::InvalidConfiguration(format!(
                        "time-zero state {} at site {site} differs from first incoming leg {state}",
                        self.initial_states[site]
                    )));
                }
            }
        }
        Ok(())
    }
}

/// Time-ordered event/link index derived from a configuration.
///
/// Vertices may remain in an unsorted packed vector. Rebuilding this index is
/// linearithmic in the expansion order and leaves diagonal add/remove blocks
/// free to use `swap_remove`.
#[derive(Debug, Clone, PartialEq)]
pub struct WorldlineIndex {
    linked: Vec<usize>,
    leg_vertex: Vec<usize>,
    leg_local: Vec<usize>,
    vertex_offsets: Vec<usize>,
    events: Vec<Vec<Event>>,
}

impl WorldlineIndex {
    /// Build the time-ordered worldline links.
    pub fn build(
        configuration: &LatticeConfiguration,
        model: &impl PositiveOperatorModel,
    ) -> Result<Self, LatticeQmcError> {
        let mut vertex_offsets = Vec::with_capacity(configuration.vertices.len() + 1);
        let mut leg_vertex = Vec::new();
        let mut leg_local = Vec::new();
        let mut events = vec![Vec::new(); model.graph().site_count()];
        vertex_offsets.push(0);

        for (vertex_id, vertex) in configuration.vertices.iter().enumerate() {
            if !vertex.tau.is_finite() || vertex.tau < 0.0 || vertex.tau >= configuration.beta {
                return Err(LatticeQmcError::InvalidConfiguration(format!(
                    "vertex {vertex_id} has invalid time {}",
                    vertex.tau
                )));
            }
            let Some(term) = model.terms().get(vertex.term) else {
                return Err(LatticeQmcError::InvalidConfiguration(format!(
                    "vertex {vertex_id} references missing term {}",
                    vertex.term
                )));
            };
            if vertex.kind >= term.kinds().len() {
                return Err(LatticeQmcError::InvalidConfiguration(format!(
                    "vertex {vertex_id} references missing kind {}",
                    vertex.kind
                )));
            }
            let offset = *vertex_offsets.last().expect("initial offset");
            for (local_site, &site) in term.sites().iter().enumerate() {
                let incoming_leg = offset + 2 * local_site;
                let outgoing_leg = incoming_leg + 1;
                events[site].push(Event {
                    tau: vertex.tau,
                    vertex: vertex_id,
                    local_site,
                    incoming_leg,
                    outgoing_leg,
                });
                leg_vertex.extend([vertex_id, vertex_id]);
                leg_local.extend([2 * local_site, 2 * local_site + 1]);
            }
            vertex_offsets.push(offset + 2 * term.sites().len());
        }

        let mut linked = vec![usize::MAX; leg_vertex.len()];
        for site_events in &mut events {
            site_events.sort_by(|left, right| {
                left.tau
                    .total_cmp(&right.tau)
                    .then(left.vertex.cmp(&right.vertex))
                    .then(left.local_site.cmp(&right.local_site))
            });
            if site_events.is_empty() {
                continue;
            }
            for position in 0..site_events.len() {
                let current = site_events[position];
                let next = site_events[(position + 1) % site_events.len()];
                linked[current.outgoing_leg] = next.incoming_leg;
                linked[next.incoming_leg] = current.outgoing_leg;
            }
        }
        if linked.contains(&usize::MAX) {
            return Err(LatticeQmcError::InvalidConfiguration(
                "not every local operator leg belongs to a worldline link".into(),
            ));
        }
        Ok(Self {
            linked,
            leg_vertex,
            leg_local,
            vertex_offsets,
            events,
        })
    }

    /// Total number of local legs.
    pub fn leg_count(&self) -> usize {
        self.linked.len()
    }

    /// Global leg offset of a vertex.
    pub fn vertex_offset(&self, vertex: usize) -> usize {
        self.vertex_offsets[vertex]
    }

    /// Vertex containing a global leg.
    pub fn vertex_of_leg(&self, leg: usize) -> usize {
        self.leg_vertex[leg]
    }

    /// Local leg index inside its term.
    pub fn local_leg(&self, leg: usize) -> usize {
        self.leg_local[leg]
    }

    /// Worldline partner of a leg.
    pub fn linked_leg(&self, leg: usize) -> usize {
        self.linked[leg]
    }

    /// Time-ordered events on a site.
    pub fn events(&self, site: usize) -> &[Event] {
        &self.events[site]
    }

    /// First event after the time-zero boundary.
    pub fn first_event(&self, site: usize) -> Option<Event> {
        self.events[site].first().copied()
    }

    /// Current state carried by a global leg.
    pub fn state_on_leg(
        &self,
        configuration: &LatticeConfiguration,
        model: &impl PositiveOperatorModel,
        leg: usize,
    ) -> BasisState {
        let vertex = &configuration.vertices[self.leg_vertex[leg]];
        model
            .term(vertex.term)
            .kind(vertex.kind)
            .state(self.leg_local[leg])
    }

    /// State immediately before a proposed time on one site.
    pub fn state_before(
        &self,
        configuration: &LatticeConfiguration,
        model: &impl PositiveOperatorModel,
        site: usize,
        tau: f64,
    ) -> BasisState {
        let events = &self.events[site];
        if events.is_empty() {
            return configuration.initial_states[site];
        }
        let position = events.partition_point(|event| event.tau < tau);
        if position < events.len() {
            return self.state_on_leg(configuration, model, events[position].incoming_leg);
        }
        self.state_on_leg(configuration, model, events[events.len() - 1].outgoing_leg)
    }

    /// Synchronize explicit time-zero states after a loop changes boundary segments.
    pub fn canonicalize_initial_states(
        &self,
        configuration: &mut LatticeConfiguration,
        model: &impl PositiveOperatorModel,
    ) {
        for site in 0..model.graph().site_count() {
            if let Some(first) = self.first_event(site) {
                let state = self.state_on_leg(configuration, model, first.incoming_leg);
                configuration.initial_states[site] = state;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::CsrGraph;
    use crate::lattice::model::SpinLatticeModel;

    #[test]
    fn empty_configuration_is_valid_for_arbitrary_spin() {
        let graph = CsrGraph::chain(5, true).expect("graph");
        let model = SpinLatticeModel::heisenberg(graph, 3, -1.0).expect("model");
        let configuration =
            LatticeConfiguration::new(8.0, vec![0, 1, 2, 3, 0], &model).expect("configuration");
        configuration.validate(&model).expect("valid");
    }
}
