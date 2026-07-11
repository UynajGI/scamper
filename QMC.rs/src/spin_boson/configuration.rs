//! Continuous-time operator configuration for one spin-1/2 impurity.

use super::error::SpinBosonError;
use super::model::SpinBosonModel;
use super::vertex::{Event, Spin, Vertex, LEGS_PER_VERTEX};

/// Sampled continuous-time retarded-vertex expansion.
#[derive(Debug, Clone, PartialEq)]
pub struct WormholeConfiguration {
    beta: f64,
    vertices: Vec<Vertex>,
    empty_spin: Spin,
}

impl WormholeConfiguration {
    /// Construct an empty operator string.
    pub fn new(beta: f64, empty_spin: Spin) -> Result<Self, SpinBosonError> {
        if !beta.is_finite() || beta <= 0.0 {
            return Err(SpinBosonError::parameter(
                "beta",
                format!("must be finite and positive, got {beta}"),
            ));
        }
        if !matches!(empty_spin, -1 | 1) {
            return Err(SpinBosonError::parameter("empty_spin", "must be -1 or +1"));
        }
        Ok(Self {
            beta,
            vertices: Vec::new(),
            empty_spin,
        })
    }

    /// Inverse temperature.
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// Set inverse temperature without changing the dimensionless time
    /// coordinates of existing vertices.
    pub fn set_beta_rescale(&mut self, beta: f64) -> Result<(), SpinBosonError> {
        if !beta.is_finite() || beta <= 0.0 {
            return Err(SpinBosonError::parameter(
                "beta",
                format!("must be finite and positive, got {beta}"),
            ));
        }
        let scale = beta / self.beta;
        for vertex in &mut self.vertices {
            vertex.tau_a *= scale;
            vertex.tau_b *= scale;
        }
        self.beta = beta;
        Ok(())
    }

    /// Retarded vertices.
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }

    /// Mutable retarded vertices, used by update kernels.
    pub(crate) fn vertices_mut(&mut self) -> &mut Vec<Vertex> {
        &mut self.vertices
    }

    /// Spin assigned to an empty operator string.
    pub fn empty_spin(&self) -> Spin {
        self.empty_spin
    }

    /// Change the empty-sector spin.
    pub(crate) fn set_empty_spin(&mut self, spin: Spin) {
        self.empty_spin = spin;
    }

    /// Expansion order.
    pub fn expansion_order(&self) -> usize {
        self.vertices.len()
    }

    /// Number of diagonal vertices.
    pub fn diagonal_order(&self, model: &SpinBosonModel) -> usize {
        self.vertices
            .iter()
            .filter(|vertex| {
                model
                    .interaction(vertex.interaction)
                    .kind(vertex.kind)
                    .is_diagonal()
            })
            .count()
    }

    /// Number of off-diagonal vertices.
    pub fn offdiagonal_order(&self, model: &SpinBosonModel) -> usize {
        self.expansion_order() - self.diagonal_order(model)
    }

    /// Time-ordered endpoint list.
    pub fn events(&self) -> Vec<Event> {
        let mut events = Vec::with_capacity(2 * self.vertices.len());
        for (vertex, sampled) in self.vertices.iter().enumerate() {
            events.push(Event {
                time: sampled.tau_a,
                vertex,
                endpoint: 0,
            });
            events.push(Event {
                time: sampled.tau_b,
                vertex,
                endpoint: 1,
            });
        }
        events.sort_by(|left, right| {
            left.time
                .total_cmp(&right.time)
                .then(left.vertex.cmp(&right.vertex))
                .then(left.endpoint.cmp(&right.endpoint))
        });
        events
    }

    /// Validate vertex bounds and full periodic worldline continuity.
    pub fn validate(&self, model: &SpinBosonModel) -> Result<(), SpinBosonError> {
        for vertex in &self.vertices {
            if !(0.0..self.beta).contains(&vertex.tau_a)
                || !(0.0..self.beta).contains(&vertex.tau_b)
            {
                return Err(SpinBosonError::InvalidConfiguration(
                    "vertex time outside [0,beta)".into(),
                ));
            }
            let interaction = model
                .interactions()
                .get(vertex.interaction)
                .ok_or_else(|| {
                    SpinBosonError::InvalidConfiguration("invalid interaction index".into())
                })?;
            if vertex.kind >= interaction.kinds().len() {
                return Err(SpinBosonError::InvalidConfiguration(
                    "invalid local vertex kind".into(),
                ));
            }
            if !vertex.omega.is_finite() || vertex.omega <= 0.0 {
                return Err(SpinBosonError::InvalidConfiguration(
                    "sampled frequency must be finite and positive".into(),
                ));
            }
        }

        let index = WorldlineIndex::build(self, model)?;
        index.validate_links(self, model)
    }
}

/// Cached sorted events and worldline leg links for one update.
#[derive(Debug, Clone, PartialEq)]
pub struct WorldlineIndex {
    events: Vec<Event>,
    links: Vec<usize>,
}

impl WorldlineIndex {
    /// Build the time-ordered circular worldline index.
    pub fn build(
        configuration: &WormholeConfiguration,
        model: &SpinBosonModel,
    ) -> Result<Self, SpinBosonError> {
        let events = configuration.events();
        if events.is_empty() {
            return Ok(Self {
                events,
                links: Vec::new(),
            });
        }
        let mut links = vec![usize::MAX; LEGS_PER_VERTEX * configuration.vertices.len()];
        for (position, event) in events.iter().enumerate() {
            let next = events[(position + 1) % events.len()];
            links[event.outgoing_leg()] = next.incoming_leg();
            links[next.incoming_leg()] = event.outgoing_leg();
        }
        if links.contains(&usize::MAX) {
            return Err(SpinBosonError::InvalidConfiguration(
                "incomplete worldline link construction".into(),
            ));
        }
        let index = Self { events, links };
        index.validate_links(configuration, model)?;
        Ok(index)
    }

    /// Sorted events.
    pub fn events(&self) -> &[Event] {
        &self.events
    }

    /// Involutive worldline partner of a global leg.
    pub fn linked_leg(&self, leg: usize) -> usize {
        self.links[leg]
    }

    /// Number of globally indexed legs.
    pub fn leg_count(&self) -> usize {
        self.links.len()
    }

    /// Spin immediately before `tau`.
    pub fn spin_before(
        &self,
        configuration: &WormholeConfiguration,
        model: &SpinBosonModel,
        tau: f64,
    ) -> Spin {
        if self.events.is_empty() {
            return configuration.empty_spin;
        }
        let first = self.events[0];
        let mut spin = event_spins(configuration, model, first).0;
        for event in &self.events {
            if event.time >= tau {
                break;
            }
            spin = event_spins(configuration, model, *event).1;
        }
        spin
    }

    /// Spin at `tau`, with times interpreted periodically.
    pub fn spin_at(
        &self,
        configuration: &WormholeConfiguration,
        model: &SpinBosonModel,
        tau: f64,
    ) -> Spin {
        if self.events.is_empty() {
            return configuration.empty_spin;
        }
        let periodic_tau = tau.rem_euclid(configuration.beta);
        let first = self.events[0];
        let mut spin = event_spins(configuration, model, first).0;
        for event in &self.events {
            if event.time > periodic_tau {
                break;
            }
            spin = event_spins(configuration, model, *event).1;
        }
        spin
    }

    fn validate_links(
        &self,
        configuration: &WormholeConfiguration,
        model: &SpinBosonModel,
    ) -> Result<(), SpinBosonError> {
        if self.events.is_empty() {
            if !matches!(configuration.empty_spin, -1 | 1) {
                return Err(SpinBosonError::InvalidConfiguration(
                    "invalid empty-sector spin".into(),
                ));
            }
            return Ok(());
        }

        let first = self.events[0];
        let (first_incoming, _) = event_spins(configuration, model, first);
        let mut propagated = first_incoming;
        for event in &self.events {
            let (incoming, outgoing) = event_spins(configuration, model, *event);
            if incoming != propagated {
                return Err(SpinBosonError::InvalidConfiguration(format!(
                    "worldline discontinuity at tau={}: expected {propagated}, found {incoming}",
                    event.time
                )));
            }
            propagated = outgoing;
        }
        if propagated != first_incoming {
            return Err(SpinBosonError::InvalidConfiguration(
                "worldline is not periodic".into(),
            ));
        }

        for (leg, partner) in self.links.iter().copied().enumerate() {
            if self.links[partner] != leg {
                return Err(SpinBosonError::InvalidConfiguration(
                    "worldline links are not involutive".into(),
                ));
            }
            let spin_left = global_leg_spin(configuration, model, leg);
            let spin_right = global_leg_spin(configuration, model, partner);
            if spin_left != spin_right {
                return Err(SpinBosonError::InvalidConfiguration(
                    "linked worldline legs carry different spins".into(),
                ));
            }
        }
        Ok(())
    }
}

/// Incoming and outgoing spins at one endpoint event.
pub fn event_spins(
    configuration: &WormholeConfiguration,
    model: &SpinBosonModel,
    event: Event,
) -> (Spin, Spin) {
    let vertex = &configuration.vertices[event.vertex];
    let kind = model.interaction(vertex.interaction).kind(vertex.kind);
    let base = 2 * event.endpoint;
    (kind.spin(base), kind.spin(base + 1))
}

/// Spin carried by one global leg.
pub fn global_leg_spin(
    configuration: &WormholeConfiguration,
    model: &SpinBosonModel,
    global_leg: usize,
) -> Spin {
    let vertex_id = global_leg / LEGS_PER_VERTEX;
    let local_leg = global_leg % LEGS_PER_VERTEX;
    let vertex = &configuration.vertices[vertex_id];
    model
        .interaction(vertex.interaction)
        .kind(vertex.kind)
        .spin(local_leg)
}

#[cfg(test)]
mod tests {
    use crate::spin_boson::bath::{Bath, SingleModeBath};
    use crate::spin_boson::model::SpinBosonModel;

    use super::*;

    #[test]
    fn empty_configuration_is_valid() {
        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = SpinBosonModel::xxz(bath, 0.5, 0.2, 0.0, None).expect("model");
        let configuration = WormholeConfiguration::new(4.0, 1).expect("configuration");
        configuration.validate(&model).expect("valid worldline");
    }
}
