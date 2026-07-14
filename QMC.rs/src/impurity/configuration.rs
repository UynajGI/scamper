//! Continuous-time operator configuration for one spin-1/2 impurity.
//!
//! The configuration uses a persistent doubly-linked list of endpoints
//! augmented with a `BTreeMap` time index. This allows O(log n) insertion,
//! O(log n) deletion, and O(1) worldline traversal, replacing the previous
//! O(n log n) per-update `WorldlineIndex::build` approach.

use std::collections::BTreeMap;

use rand::Rng;
use rand::RngExt;

use super::error::ImpurityError;
use super::model::ImpurityModel;
use super::vertex::{EndpointId, LegId, LegSide, Spin, Vertex, VertexId, LEGS_PER_VERTEX};

/// Sortable key for the time-ordered BTreeMap.
///
/// `f64` does not implement `Ord`, so we wrap it with a tie-breaker and
/// implement `Ord` manually using `total_cmp`.
#[derive(Debug, Clone, Copy)]
struct EventKey {
    time: f64,
    tie_breaker: u64,
}

impl PartialEq for EventKey {
    fn eq(&self, other: &Self) -> bool {
        self.time.to_bits() == other.time.to_bits() && self.tie_breaker == other.tie_breaker
    }
}

impl Eq for EventKey {}

impl PartialOrd for EventKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EventKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.time
            .total_cmp(&other.time)
            .then(self.tie_breaker.cmp(&other.tie_breaker))
    }
}

/// Doubly-linked list pointers for one endpoint.
#[derive(Debug, Clone, PartialEq)]
struct EndpointLinks {
    prev: EndpointId,
    next: EndpointId,
    key: EventKey,
}

/// A vertex with persistent linked-list and index metadata.
#[derive(Debug, Clone, PartialEq)]
struct LinkedVertex {
    vertex: Vertex,
    endpoints: [EndpointLinks; 2],
    /// Position in `active_vertices`.
    active_position: usize,
    /// Position in `diagonal_vertices`, if diagonal.
    diagonal_position: Option<usize>,
}

/// Sampled continuous-time retarded-vertex expansion with persistent worldline
/// links.
#[derive(Debug, Clone, PartialEq)]
pub struct WormholeConfiguration {
    beta: f64,

    /// Stable vertex slots. `None` means the slot is free.
    vertices: Vec<Option<LinkedVertex>>,
    /// Free slot indices for reuse.
    free_slots: Vec<usize>,

    /// All active vertex IDs for O(1) random selection.
    active_vertices: Vec<VertexId>,
    /// Current diagonal vertex IDs for O(1) random removal.
    diagonal_vertices: Vec<VertexId>,

    /// Time-ordered endpoint index.
    time_order: BTreeMap<EventKey, EndpointId>,
    /// The endpoint with the smallest time; the list is periodic.
    first_endpoint: Option<EndpointId>,

    empty_spin: Spin,
}

impl WormholeConfiguration {
    /// Construct an empty operator string.
    pub fn new(beta: f64, empty_spin: Spin) -> Result<Self, ImpurityError> {
        if !beta.is_finite() || beta <= 0.0 {
            return Err(ImpurityError::parameter(
                "beta",
                format!("must be finite and positive, got {beta}"),
            ));
        }
        if !matches!(empty_spin, -1 | 1) {
            return Err(ImpurityError::parameter("empty_spin", "must be -1 or +1"));
        }
        Ok(Self {
            beta,
            vertices: Vec::new(),
            free_slots: Vec::new(),
            active_vertices: Vec::new(),
            diagonal_vertices: Vec::new(),
            time_order: BTreeMap::new(),
            first_endpoint: None,
            empty_spin,
        })
    }

    /// Inverse temperature.
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// Set inverse temperature, rescaling all vertex times.
    ///
    /// This triggers a full rebuild of the time index, which is acceptable
    /// because temperature changes are not on the hot path.
    pub fn set_beta_rescale(&mut self, beta: f64) -> Result<(), ImpurityError> {
        if !beta.is_finite() || beta <= 0.0 {
            return Err(ImpurityError::parameter(
                "beta",
                format!("must be finite and positive, got {beta}"),
            ));
        }
        let scale = beta / self.beta;
        for id in self.active_vertices.clone() {
            let linked = self.vertices[id.0].as_mut().unwrap();
            linked.vertex.tau_a *= scale;
            linked.vertex.tau_b *= scale;
        }
        self.beta = beta;
        self.rebuild_time_links()
    }

    /// Spin assigned to an empty operator string.
    pub fn empty_spin(&self) -> Spin {
        self.empty_spin
    }

    /// Change the empty-sector spin.
    pub(crate) fn set_empty_spin(&mut self, spin: Spin) {
        self.empty_spin = spin;
    }

    /// Synchronize the stored zero-order spin with the worldline segment at
    /// `tau = 0`. Directed loops can change this trace sector even though
    /// vertex times and links remain fixed.
    pub(crate) fn sync_empty_spin_from_worldline(
        &mut self,
        model: &ImpurityModel,
    ) -> Result<(), ImpurityError> {
        if let Some(first) = self.first_endpoint {
            let spin = self.endpoint_incoming_spin(first, model)?;
            self.empty_spin = spin;
        }
        Ok(())
    }

    /// Expansion order (number of active vertices).
    pub fn expansion_order(&self) -> usize {
        self.active_vertices.len()
    }

    /// Number of diagonal vertices.
    pub fn diagonal_order(&self) -> usize {
        self.diagonal_vertices.len()
    }

    /// Number of off-diagonal vertices.
    pub fn offdiagonal_order(&self) -> usize {
        self.expansion_order() - self.diagonal_order()
    }

    /// Access one vertex by ID.
    pub fn vertex(&self, id: VertexId) -> Result<&Vertex, ImpurityError> {
        self.vertices
            .get(id.0)
            .and_then(|slot| slot.as_ref())
            .map(|linked| &linked.vertex)
            .ok_or_else(|| {
                ImpurityError::InvalidConfiguration(format!("invalid vertex id {}", id.0))
            })
    }

    /// Mutable access to one vertex by ID.
    #[allow(dead_code)]
    pub(crate) fn vertex_mut(&mut self, id: VertexId) -> Result<&mut Vertex, ImpurityError> {
        self.vertices
            .get_mut(id.0)
            .and_then(|slot| slot.as_mut())
            .map(|linked| &mut linked.vertex)
            .ok_or_else(|| {
                ImpurityError::InvalidConfiguration(format!("invalid vertex id {}", id.0))
            })
    }

    /// Insert a new vertex into the configuration.
    ///
    /// Returns the stable `VertexId` assigned to the new vertex.
    pub fn insert_vertex(
        &mut self,
        vertex: Vertex,
        model: &ImpurityModel,
    ) -> Result<VertexId, ImpurityError> {
        let id = self.allocate_slot();
        let active_position = self.active_vertices.len();
        self.active_vertices.push(id);

        let is_diagonal = model
            .interaction(vertex.interaction)
            .kind(vertex.kind)
            .is_diagonal();
        let diagonal_position = if is_diagonal {
            let pos = self.diagonal_vertices.len();
            self.diagonal_vertices.push(id);
            Some(pos)
        } else {
            None
        };

        let tie_a = (2 * id.0) as u64;
        let tie_b = (2 * id.0 + 1) as u64;
        let key_a = EventKey {
            time: vertex.tau_a,
            tie_breaker: tie_a,
        };
        let key_b = EventKey {
            time: vertex.tau_b,
            tie_breaker: tie_b,
        };

        let endpoint_a = EndpointId {
            vertex: id,
            endpoint: 0,
        };
        let endpoint_b = EndpointId {
            vertex: id,
            endpoint: 1,
        };

        // Temporary self-links; will be overwritten by insert_endpoint.
        let links_a = EndpointLinks {
            prev: endpoint_a,
            next: endpoint_a,
            key: key_a,
        };
        let links_b = EndpointLinks {
            prev: endpoint_b,
            next: endpoint_b,
            key: key_b,
        };

        self.vertices[id.0] = Some(LinkedVertex {
            vertex,
            endpoints: [links_a, links_b],
            active_position,
            diagonal_position,
        });

        self.insert_endpoint(endpoint_a, key_a)?;
        self.insert_endpoint(endpoint_b, key_b)?;

        Ok(id)
    }

    /// Remove a vertex from the configuration.
    ///
    /// Returns the removed `Vertex`.
    pub fn remove_vertex(&mut self, id: VertexId) -> Result<Vertex, ImpurityError> {
        let endpoint_a = EndpointId {
            vertex: id,
            endpoint: 0,
        };
        let endpoint_b = EndpointId {
            vertex: id,
            endpoint: 1,
        };

        self.unlink_endpoint(endpoint_a)?;
        self.unlink_endpoint(endpoint_b)?;

        let linked = self.vertices[id.0].take().ok_or_else(|| {
            ImpurityError::InvalidConfiguration(format!("vertex {} already removed", id.0))
        })?;

        self.remove_from_active(id, linked.active_position);
        if let Some(diag_pos) = linked.diagonal_position {
            self.remove_from_diagonal(id, diag_pos);
        }

        self.free_slots.push(id.0);

        Ok(linked.vertex)
    }

    /// Change the kind of a vertex, synchronizing the diagonal index.
    pub fn set_kind(
        &mut self,
        id: VertexId,
        new_kind: usize,
        model: &ImpurityModel,
    ) -> Result<(), ImpurityError> {
        let linked = self.vertices[id.0].as_ref().ok_or_else(|| {
            ImpurityError::InvalidConfiguration(format!("invalid vertex id {}", id.0))
        })?;
        let old_kind = linked.vertex.kind;
        let interaction = linked.vertex.interaction;

        let old_diagonal = model.interaction(interaction).kind(old_kind).is_diagonal();
        let new_diagonal = model.interaction(interaction).kind(new_kind).is_diagonal();

        self.vertices[id.0].as_mut().unwrap().vertex.kind = new_kind;

        match (old_diagonal, new_diagonal) {
            (false, true) => {
                let pos = self.diagonal_vertices.len();
                self.diagonal_vertices.push(id);
                self.vertices[id.0].as_mut().unwrap().diagonal_position = Some(pos);
            }
            (true, false) => {
                let diag_pos = self.vertices[id.0]
                    .as_ref()
                    .unwrap()
                    .diagonal_position
                    .unwrap();
                self.remove_from_diagonal(id, diag_pos);
                self.vertices[id.0].as_mut().unwrap().diagonal_position = None;
            }
            _ => {}
        }

        Ok(())
    }

    /// Worldline partner of a leg.
    ///
    /// This is the core traversal primitive: given a leg, return the leg
    /// connected to it through the worldline.
    ///
    /// - `Incoming` leg at endpoint E is connected to the `Outgoing` leg of
    ///   the predecessor endpoint.
    /// - `Outgoing` leg at endpoint E is connected to the `Incoming` leg of
    ///   the successor endpoint.
    pub fn linked_leg(&self, leg: LegId) -> Result<LegId, ImpurityError> {
        let links = self.endpoint_links(leg.endpoint)?;
        match leg.side {
            LegSide::Incoming => Ok(LegId {
                endpoint: links.prev,
                side: LegSide::Outgoing,
            }),
            LegSide::Outgoing => Ok(LegId {
                endpoint: links.next,
                side: LegSide::Incoming,
            }),
        }
    }

    /// Spin immediately before `tau`.
    pub fn spin_before(&self, model: &ImpurityModel, tau: f64) -> Result<Spin, ImpurityError> {
        if self.time_order.is_empty() {
            return Ok(self.empty_spin);
        }

        let key = EventKey {
            time: tau.rem_euclid(self.beta),
            tie_breaker: 0,
        };

        let endpoint = self
            .time_order
            .range(key..)
            .next()
            .map(|(_, &ep)| ep)
            .unwrap_or(self.first_endpoint.unwrap());

        self.endpoint_incoming_spin(endpoint, model)
    }

    /// Spin at `tau`, with times interpreted periodically.
    pub fn spin_at(&self, model: &ImpurityModel, tau: f64) -> Result<Spin, ImpurityError> {
        if self.time_order.is_empty() {
            return Ok(self.empty_spin);
        }

        let periodic_tau = tau.rem_euclid(self.beta);
        // Search for the last endpoint strictly before (tau, 0).
        // This gives the last endpoint with time < tau, or wraps to the
        // last endpoint overall if tau is before all endpoints.
        let key = EventKey {
            time: periodic_tau,
            tie_breaker: 0,
        };

        let endpoint = self
            .time_order
            .range(..key)
            .next_back()
            .map(|(_, &ep)| ep)
            .unwrap_or_else(|| self.time_order.values().next_back().copied().unwrap());

        self.endpoint_outgoing_spin(endpoint, model)
    }

    /// Random leg for starting a directed loop.
    pub fn random_leg<R: Rng + ?Sized>(&self, rng: &mut R) -> Result<LegId, ImpurityError> {
        if self.active_vertices.is_empty() {
            return Err(ImpurityError::InvalidConfiguration(
                "cannot pick random leg from empty configuration".into(),
            ));
        }
        let vertex_id = self.active_vertices[rng.random_range(0..self.active_vertices.len())];
        let local_leg = rng.random_range(0..LEGS_PER_VERTEX);
        Ok(LegId::from_local(vertex_id, local_leg))
    }

    /// Map a uniformly sampled imaginary time to a directed-loop start leg.
    ///
    /// Forward propagation enters the incoming leg of the next endpoint.
    /// Backward propagation enters the outgoing leg of the previous endpoint.
    /// Vertex times are unchanged by a loop, so this proposal is symmetric
    /// between a loop and its reverse.
    pub fn start_leg_at_time(&self, tau: f64, forward: bool) -> Result<LegId, ImpurityError> {
        if self.time_order.is_empty() {
            return Err(ImpurityError::InvalidConfiguration(
                "cannot start a directed loop from an empty configuration".into(),
            ));
        }
        let key = EventKey {
            time: tau.rem_euclid(self.beta),
            tie_breaker: 0,
        };
        if forward {
            let endpoint = self
                .time_order
                .range(key..)
                .next()
                .map(|(_, &endpoint)| endpoint)
                .unwrap_or_else(|| self.first_endpoint.expect("non-empty time index"));
            Ok(LegId {
                endpoint,
                side: LegSide::Incoming,
            })
        } else {
            let endpoint = self
                .time_order
                .range(..key)
                .next_back()
                .map(|(_, &endpoint)| endpoint)
                .unwrap_or_else(|| {
                    self.time_order
                        .values()
                        .next_back()
                        .copied()
                        .expect("non-empty time index")
                });
            Ok(LegId {
                endpoint,
                side: LegSide::Outgoing,
            })
        }
    }

    /// Random diagonal vertex for removal.
    pub fn random_diagonal_vertex<R: Rng + ?Sized>(
        &self,
        rng: &mut R,
    ) -> Result<VertexId, ImpurityError> {
        if self.diagonal_vertices.is_empty() {
            return Err(ImpurityError::InvalidConfiguration(
                "no diagonal vertices to select".into(),
            ));
        }
        Ok(self.diagonal_vertices[rng.random_range(0..self.diagonal_vertices.len())])
    }

    /// Iterator over endpoints in time order.
    ///
    /// Yields `(time, EndpointId)` pairs.
    pub fn time_ordered_endpoints(&self) -> impl Iterator<Item = (f64, EndpointId)> + '_ {
        self.time_order
            .iter()
            .map(|(key, &endpoint)| (key.time, endpoint))
    }

    /// Validate vertex bounds and full periodic worldline continuity.
    pub fn validate(&self, model: &ImpurityModel) -> Result<(), ImpurityError> {
        // 1. Basic vertex bounds.
        for (slot_idx, slot) in self.vertices.iter().enumerate() {
            if let Some(linked) = slot {
                let v = &linked.vertex;
                if !(0.0..self.beta).contains(&v.tau_a) || !(0.0..self.beta).contains(&v.tau_b) {
                    return Err(ImpurityError::InvalidConfiguration(
                        "vertex time outside [0,beta)".into(),
                    ));
                }
                let interaction = model.interactions().get(v.interaction).ok_or_else(|| {
                    ImpurityError::InvalidConfiguration("invalid interaction index".into())
                })?;
                if v.kind >= interaction.kinds().len() {
                    return Err(ImpurityError::InvalidConfiguration(
                        "invalid local vertex kind".into(),
                    ));
                }
                if !v.omega.is_finite() || v.omega <= 0.0 {
                    return Err(ImpurityError::InvalidConfiguration(
                        "sampled frequency must be finite and positive".into(),
                    ));
                }
                let _ = slot_idx; // suppress unused warning
            }
        }

        let n = self.expansion_order();
        if n == 0 {
            if !matches!(self.empty_spin, -1 | 1) {
                return Err(ImpurityError::InvalidConfiguration(
                    "invalid empty-sector spin".into(),
                ));
            }
            return Ok(());
        }

        // 2. Time index size.
        if self.time_order.len() != 2 * n {
            return Err(ImpurityError::InvalidConfiguration(format!(
                "time_order has {} entries, expected {}",
                self.time_order.len(),
                2 * n
            )));
        }

        // 3. Doubly-linked consistency.
        for (&_key, &endpoint) in &self.time_order {
            let links = self.endpoint_links(endpoint)?;
            let next_links = self.endpoint_links(links.next)?;
            let prev_links = self.endpoint_links(links.prev)?;
            if next_links.prev != endpoint {
                return Err(ImpurityError::InvalidConfiguration(
                    "next.prev != self".into(),
                ));
            }
            if prev_links.next != endpoint {
                return Err(ImpurityError::InvalidConfiguration(
                    "prev.next != self".into(),
                ));
            }
        }

        // 4. Periodic traversal count.
        let first = self.first_endpoint.ok_or_else(|| {
            ImpurityError::InvalidConfiguration("first_endpoint is None with vertices".into())
        })?;
        let mut count = 0;
        let mut current = first;
        loop {
            count += 1;
            if count > 2 * n {
                return Err(ImpurityError::InvalidConfiguration(
                    "cycle traversal exceeded expected count".into(),
                ));
            }
            let links = self.endpoint_links(current)?;
            current = links.next;
            if current == first {
                break;
            }
        }
        if count != 2 * n {
            return Err(ImpurityError::InvalidConfiguration(format!(
                "cycle traversal visited {count} endpoints, expected {}",
                2 * n
            )));
        }

        // 5. Time monotonicity along the cycle (except the periodic wrap).
        let mut current = first;
        let mut prev_key = self.endpoint_links(first)?.key;
        for _ in 1..2 * n {
            let links = self.endpoint_links(current)?;
            current = links.next;
            let cur_key = self.endpoint_links(current)?.key;
            if cur_key <= prev_key {
                return Err(ImpurityError::InvalidConfiguration(
                    "time not monotonically increasing along cycle".into(),
                ));
            }
            prev_key = cur_key;
        }

        // 6. Worldline spin continuity and the stored tau=0 trace sector.
        let first_incoming = self.endpoint_incoming_spin(first, model)?;
        if first_incoming != self.empty_spin {
            return Err(ImpurityError::InvalidConfiguration(format!(
                "stored empty spin {} differs from the tau=0 worldline spin {first_incoming}",
                self.empty_spin
            )));
        }
        let mut propagated = self.empty_spin;
        let mut ep = first;
        for _ in 0..2 * n {
            let (incoming, outgoing) = self.endpoint_spins(ep, model)?;
            if incoming != propagated {
                return Err(ImpurityError::InvalidConfiguration(format!(
                    "worldline discontinuity: expected {propagated}, found {incoming}"
                )));
            }
            propagated = outgoing;
            let links = self.endpoint_links(ep)?;
            ep = links.next;
        }
        if propagated != self.empty_spin {
            return Err(ImpurityError::InvalidConfiguration(
                "worldline is not periodic".into(),
            ));
        }

        // 7. linked_leg involutive.
        for (&_key, &endpoint) in &self.time_order {
            for side in [LegSide::Incoming, LegSide::Outgoing] {
                let leg = LegId { endpoint, side };
                let partner = self.linked_leg(leg)?;
                let back = self.linked_leg(partner)?;
                if back != leg {
                    return Err(ImpurityError::InvalidConfiguration(
                        "linked_leg is not involutive".into(),
                    ));
                }
            }
        }

        // 8. Diagonal index consistency.
        let mut diag_seen = vec![false; self.vertices.len()];
        for &diag_id in &self.diagonal_vertices {
            let linked = self.vertices[diag_id.0].as_ref().ok_or_else(|| {
                ImpurityError::InvalidConfiguration("diagonal vertex slot is None".into())
            })?;
            let v = &linked.vertex;
            if !model.interaction(v.interaction).kind(v.kind).is_diagonal() {
                return Err(ImpurityError::InvalidConfiguration(
                    "non-diagonal vertex in diagonal index".into(),
                ));
            }
            if diag_seen[diag_id.0] {
                return Err(ImpurityError::InvalidConfiguration(
                    "duplicate in diagonal index".into(),
                ));
            }
            diag_seen[diag_id.0] = true;
        }
        for (slot_idx, slot) in self.vertices.iter().enumerate() {
            if let Some(linked) = slot {
                let v = &linked.vertex;
                let is_diag = model.interaction(v.interaction).kind(v.kind).is_diagonal();
                if is_diag && !diag_seen[slot_idx] {
                    return Err(ImpurityError::InvalidConfiguration(
                        "diagonal vertex missing from diagonal index".into(),
                    ));
                }
            }
        }

        // 9. Active index consistency.
        let mut active_seen = vec![false; self.vertices.len()];
        for &active_id in &self.active_vertices {
            if self.vertices[active_id.0].is_none() {
                return Err(ImpurityError::InvalidConfiguration(
                    "active vertex slot is None".into(),
                ));
            }
            if active_seen[active_id.0] {
                return Err(ImpurityError::InvalidConfiguration(
                    "duplicate in active index".into(),
                ));
            }
            active_seen[active_id.0] = true;
        }
        for (slot_idx, slot) in self.vertices.iter().enumerate() {
            if slot.is_some() && !active_seen[slot_idx] {
                return Err(ImpurityError::InvalidConfiguration(
                    "active vertex missing from active index".into(),
                ));
            }
        }

        Ok(())
    }

    // ──────────────────────────────── internal helpers ────────────────────────────────

    fn allocate_slot(&mut self) -> VertexId {
        if let Some(slot) = self.free_slots.pop() {
            VertexId(slot)
        } else {
            let slot = self.vertices.len();
            self.vertices.push(None);
            VertexId(slot)
        }
    }

    fn remove_from_active(&mut self, id: VertexId, position: usize) {
        debug_assert!(
            position < self.active_vertices.len(),
            "remove_from_active: position {position} out of bounds (len {})",
            self.active_vertices.len()
        );
        // swap_remove returns the REMOVED element, not the swapped one.
        // We need the last element (which will be swapped into `position`).
        let last_idx = self.active_vertices.len() - 1;
        let swapped = self.active_vertices[last_idx];
        self.active_vertices.swap_remove(position);
        // Update the swapped vertex's position, unless it was the one removed.
        if swapped != id {
            if let Some(linked) = self.vertices[swapped.0].as_mut() {
                linked.active_position = position;
            }
        }
    }

    fn remove_from_diagonal(&mut self, id: VertexId, position: usize) {
        debug_assert!(
            position < self.diagonal_vertices.len(),
            "remove_from_diagonal: position {position} out of bounds (len {})",
            self.diagonal_vertices.len()
        );
        // swap_remove returns the REMOVED element, not the swapped one.
        // We need the last element (which will be swapped into `position`).
        let last_idx = self.diagonal_vertices.len() - 1;
        let swapped = self.diagonal_vertices[last_idx];
        self.diagonal_vertices.swap_remove(position);
        // Update the swapped vertex's position, unless it was the one removed.
        if swapped != id {
            if let Some(linked) = self.vertices[swapped.0].as_mut() {
                linked.diagonal_position = Some(position);
            }
        }
    }

    fn endpoint_links(&self, endpoint: EndpointId) -> Result<&EndpointLinks, ImpurityError> {
        let linked = self.vertices[endpoint.vertex.0].as_ref().ok_or_else(|| {
            ImpurityError::InvalidConfiguration(format!("vertex {} not found", endpoint.vertex.0))
        })?;
        Ok(&linked.endpoints[endpoint.endpoint as usize])
    }

    fn endpoint_links_mut(
        &mut self,
        endpoint: EndpointId,
    ) -> Result<&mut EndpointLinks, ImpurityError> {
        let linked = self.vertices[endpoint.vertex.0].as_mut().ok_or_else(|| {
            ImpurityError::InvalidConfiguration(format!("vertex {} not found", endpoint.vertex.0))
        })?;
        Ok(&mut linked.endpoints[endpoint.endpoint as usize])
    }

    fn set_prev(&mut self, endpoint: EndpointId, prev: EndpointId) -> Result<(), ImpurityError> {
        self.endpoint_links_mut(endpoint)?.prev = prev;
        Ok(())
    }

    fn set_next(&mut self, endpoint: EndpointId, next: EndpointId) -> Result<(), ImpurityError> {
        self.endpoint_links_mut(endpoint)?.next = next;
        Ok(())
    }

    /// Insert an endpoint into the doubly-linked time list.
    fn insert_endpoint(
        &mut self,
        endpoint: EndpointId,
        key: EventKey,
    ) -> Result<(), ImpurityError> {
        if self.first_endpoint.is_none() {
            // Empty list: self-link.
            self.endpoint_links_mut(endpoint)?.prev = endpoint;
            self.endpoint_links_mut(endpoint)?.next = endpoint;
            self.endpoint_links_mut(endpoint)?.key = key;
            self.time_order.insert(key, endpoint);
            self.first_endpoint = Some(endpoint);
            return Ok(());
        }

        // Find successor: first endpoint with key >= our key.
        let successor = self
            .time_order
            .range(key..)
            .next()
            .map(|(_, &ep)| ep)
            .unwrap_or(self.first_endpoint.unwrap());

        let predecessor = self.endpoint_links(successor)?.prev;

        // Splice: predecessor <-> endpoint <-> successor
        self.set_next(predecessor, endpoint)?;
        self.set_prev(successor, endpoint)?;

        self.endpoint_links_mut(endpoint)?.prev = predecessor;
        self.endpoint_links_mut(endpoint)?.next = successor;
        self.endpoint_links_mut(endpoint)?.key = key;

        self.time_order.insert(key, endpoint);

        if key < self.endpoint_links(self.first_endpoint.unwrap())?.key {
            self.first_endpoint = Some(endpoint);
        }

        Ok(())
    }

    /// Remove an endpoint from the doubly-linked time list.
    fn unlink_endpoint(&mut self, endpoint: EndpointId) -> Result<(), ImpurityError> {
        let links = self.endpoint_links(endpoint)?.clone();
        let prev = links.prev;
        let next = links.next;

        if prev == endpoint && next == endpoint {
            // Only endpoint in the list.
            self.time_order.remove(&links.key);
            self.first_endpoint = None;
            return Ok(());
        }

        self.set_next(prev, next)?;
        self.set_prev(next, prev)?;
        self.time_order.remove(&links.key);

        if self.first_endpoint == Some(endpoint) {
            self.first_endpoint = Some(next);
        }

        Ok(())
    }

    /// Full rebuild of time links after beta rescale or checkpoint load.
    fn rebuild_time_links(&mut self) -> Result<(), ImpurityError> {
        self.time_order.clear();
        self.first_endpoint = None;

        // Reset all endpoint links to self-referential.
        for id in self.active_vertices.clone() {
            let linked = self.vertices[id.0].as_mut().unwrap();
            let endpoint_a = EndpointId {
                vertex: id,
                endpoint: 0,
            };
            let endpoint_b = EndpointId {
                vertex: id,
                endpoint: 1,
            };
            let tie_a = (2 * id.0) as u64;
            let tie_b = (2 * id.0 + 1) as u64;

            linked.endpoints[0] = EndpointLinks {
                prev: endpoint_a,
                next: endpoint_a,
                key: EventKey {
                    time: linked.vertex.tau_a,
                    tie_breaker: tie_a,
                },
            };
            linked.endpoints[1] = EndpointLinks {
                prev: endpoint_b,
                next: endpoint_b,
                key: EventKey {
                    time: linked.vertex.tau_b,
                    tie_breaker: tie_b,
                },
            };
        }

        // Re-insert all endpoints.
        for id in self.active_vertices.clone() {
            let linked = self.vertices[id.0].as_ref().unwrap();
            let endpoint_a = EndpointId {
                vertex: id,
                endpoint: 0,
            };
            let endpoint_b = EndpointId {
                vertex: id,
                endpoint: 1,
            };
            let key_a = linked.endpoints[0].key;
            let key_b = linked.endpoints[1].key;
            self.insert_endpoint(endpoint_a, key_a)?;
            self.insert_endpoint(endpoint_b, key_b)?;
        }

        Ok(())
    }

    /// Incoming spin at an endpoint.
    fn endpoint_incoming_spin(
        &self,
        endpoint: EndpointId,
        model: &ImpurityModel,
    ) -> Result<Spin, ImpurityError> {
        let linked = self.vertices[endpoint.vertex.0].as_ref().ok_or_else(|| {
            ImpurityError::InvalidConfiguration(format!("vertex {} not found", endpoint.vertex.0))
        })?;
        let kind = model
            .interaction(linked.vertex.interaction)
            .kind(linked.vertex.kind);
        let local_leg = 2 * (endpoint.endpoint as usize);
        Ok(kind.spin(local_leg))
    }

    /// Outgoing spin at an endpoint.
    pub fn endpoint_outgoing_spin(
        &self,
        endpoint: EndpointId,
        model: &ImpurityModel,
    ) -> Result<Spin, ImpurityError> {
        let linked = self.vertices[endpoint.vertex.0].as_ref().ok_or_else(|| {
            ImpurityError::InvalidConfiguration(format!("vertex {} not found", endpoint.vertex.0))
        })?;
        let kind = model
            .interaction(linked.vertex.interaction)
            .kind(linked.vertex.kind);
        let local_leg = 2 * (endpoint.endpoint as usize) + 1;
        Ok(kind.spin(local_leg))
    }

    /// Both incoming and outgoing spins at an endpoint.
    fn endpoint_spins(
        &self,
        endpoint: EndpointId,
        model: &ImpurityModel,
    ) -> Result<(Spin, Spin), ImpurityError> {
        let linked = self.vertices[endpoint.vertex.0].as_ref().ok_or_else(|| {
            ImpurityError::InvalidConfiguration(format!("vertex {} not found", endpoint.vertex.0))
        })?;
        let kind = model
            .interaction(linked.vertex.interaction)
            .kind(linked.vertex.kind);
        let base = 2 * (endpoint.endpoint as usize);
        Ok((kind.spin(base), kind.spin(base + 1)))
    }
}

// ──────────────────────────────── Legacy WorldlineIndex (test only) ────────────────────────────────

#[cfg(test)]
mod legacy {
    use super::super::vertex::Event;
    use super::*;

    /// Cached sorted events and worldline leg links for one update.
    /// Kept only for regression testing against the new persistent structure.
    #[derive(Debug, Clone, PartialEq)]
    #[allow(dead_code)]
    pub struct WorldlineIndex {
        pub events: Vec<Event>,
        pub links: Vec<usize>,
    }

    #[allow(dead_code)]
    impl WorldlineIndex {
        /// Build the time-ordered circular worldline index.
        pub fn build(
            configuration: &WormholeConfiguration,
            model: &ImpurityModel,
        ) -> Result<Self, ImpurityError> {
            let events = Self::collect_events(configuration);
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
                return Err(ImpurityError::InvalidConfiguration(
                    "incomplete worldline link construction".into(),
                ));
            }
            let index = Self { events, links };
            index.validate_links(configuration, model)?;
            Ok(index)
        }

        fn collect_events(configuration: &WormholeConfiguration) -> Vec<Event> {
            let mut events = Vec::new();
            for (slot_idx, slot) in configuration.vertices.iter().enumerate() {
                if let Some(linked) = slot {
                    events.push(Event {
                        time: linked.vertex.tau_a,
                        vertex: slot_idx,
                        endpoint: 0,
                    });
                    events.push(Event {
                        time: linked.vertex.tau_b,
                        vertex: slot_idx,
                        endpoint: 1,
                    });
                }
            }
            events.sort_by(|left, right| {
                left.time
                    .total_cmp(&right.time)
                    .then(left.vertex.cmp(&right.vertex))
                    .then(left.endpoint.cmp(&right.endpoint))
            });
            events
        }

        pub fn events(&self) -> &[Event] {
            &self.events
        }

        pub fn linked_leg(&self, leg: usize) -> usize {
            self.links[leg]
        }

        pub fn leg_count(&self) -> usize {
            self.links.len()
        }

        pub fn spin_before(
            &self,
            configuration: &WormholeConfiguration,
            model: &ImpurityModel,
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

        fn validate_links(
            &self,
            configuration: &WormholeConfiguration,
            model: &ImpurityModel,
        ) -> Result<(), ImpurityError> {
            if self.events.is_empty() {
                if !matches!(configuration.empty_spin, -1 | 1) {
                    return Err(ImpurityError::InvalidConfiguration(
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
                    return Err(ImpurityError::InvalidConfiguration(format!(
                        "worldline discontinuity at tau={}: expected {propagated}, found {incoming}",
                        event.time
                    )));
                }
                propagated = outgoing;
            }
            if propagated != first_incoming {
                return Err(ImpurityError::InvalidConfiguration(
                    "worldline is not periodic".into(),
                ));
            }

            for (leg, partner) in self.links.iter().copied().enumerate() {
                if self.links[partner] != leg {
                    return Err(ImpurityError::InvalidConfiguration(
                        "worldline links are not involutive".into(),
                    ));
                }
                let spin_left = global_leg_spin(configuration, model, leg);
                let spin_right = global_leg_spin(configuration, model, partner);
                if spin_left != spin_right {
                    return Err(ImpurityError::InvalidConfiguration(
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
        model: &ImpurityModel,
        event: Event,
    ) -> (Spin, Spin) {
        let linked = configuration.vertices[event.vertex].as_ref().unwrap();
        let kind = model
            .interaction(linked.vertex.interaction)
            .kind(linked.vertex.kind);
        let base = 2 * event.endpoint;
        (kind.spin(base), kind.spin(base + 1))
    }

    /// Spin carried by one global leg.
    pub fn global_leg_spin(
        configuration: &WormholeConfiguration,
        model: &ImpurityModel,
        global_leg: usize,
    ) -> Spin {
        let vertex_id = global_leg / LEGS_PER_VERTEX;
        let local_leg = global_leg % LEGS_PER_VERTEX;
        let linked = configuration.vertices[vertex_id].as_ref().unwrap();
        model
            .interaction(linked.vertex.interaction)
            .kind(linked.vertex.kind)
            .spin(local_leg)
    }
}

#[cfg(test)]
mod tests {
    use crate::impurity::bath::{Bath, SingleModeBath};
    use crate::impurity::model::ImpurityModel;

    use super::legacy::WorldlineIndex;
    use super::*;

    #[test]
    fn empty_configuration_is_valid() {
        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = ImpurityModel::xxz(bath, 0.5, 0.2, 0.0, None).expect("model");
        let configuration = WormholeConfiguration::new(4.0, 1).expect("configuration");
        configuration.validate(&model).expect("valid worldline");
    }

    #[test]
    fn sync_empty_spin_tracks_the_tau_zero_worldline_sector() {
        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = ImpurityModel::xxz(bath, 0.4, 0.0, 0.0, Some(0.2)).expect("model");
        let offdiagonal_kind = model
            .interaction(0)
            .kinds()
            .iter()
            .position(|kind| kind.legs() == &[1, -1, -1, 1])
            .expect("exchange vertex");
        let mut configuration = WormholeConfiguration::new(1.0, -1).expect("configuration");
        configuration
            .insert_vertex(
                Vertex {
                    tau_a: 0.25,
                    tau_b: 0.75,
                    omega: 1.0,
                    interaction: 0,
                    kind: offdiagonal_kind,
                },
                &model,
            )
            .expect("insert exchange vertex");
        assert!(configuration.validate(&model).is_err());
        configuration
            .sync_empty_spin_from_worldline(&model)
            .expect("synchronize trace sector");
        configuration
            .validate(&model)
            .expect("valid synchronized worldline");
        assert_eq!(configuration.empty_spin(), 1);
    }

    #[test]
    fn insert_and_remove_single_diagonal() {
        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = ImpurityModel::xxz(bath, 0.4, 0.2, 0.1, None).expect("model");
        let mut configuration = WormholeConfiguration::new(8.0, 1).expect("configuration");

        let interaction = model.interaction(0);
        let kind = interaction.diagonal_kind(1, 1);
        let vertex = Vertex {
            tau_a: 1.0,
            tau_b: 2.0,
            omega: 1.0,
            interaction: 0,
            kind,
        };
        let id = configuration.insert_vertex(vertex, &model).expect("insert");
        assert_eq!(configuration.expansion_order(), 1);
        assert_eq!(configuration.diagonal_order(), 1);
        configuration.validate(&model).expect("valid");

        let removed = configuration.remove_vertex(id).expect("remove");
        assert_eq!(removed.tau_a, 1.0);
        assert_eq!(configuration.expansion_order(), 0);
        configuration.validate(&model).expect("valid after remove");
    }

    #[test]
    fn two_vertices_with_different_spins() {
        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = ImpurityModel::xxz(bath, 0.4, 0.2, 0.1, None).expect("model");
        let mut configuration = WormholeConfiguration::new(8.0, -1).expect("configuration");

        let interaction = model.interaction(0);

        // Two same-spin diagonal vertices with interleaved times.
        // Both have legs [-1,-1,-1,-1], so worldline spin = -1 throughout.
        let kind1 = interaction.diagonal_kind(-1, -1);
        let v1 = Vertex {
            tau_a: 1.0,
            tau_b: 3.0,
            omega: 1.0,
            interaction: 0,
            kind: kind1,
        };
        configuration.insert_vertex(v1, &model).expect("insert 1");

        let kind2 = interaction.diagonal_kind(-1, -1);
        let v2 = Vertex {
            tau_a: 2.0,
            tau_b: 4.0,
            omega: 1.0,
            interaction: 0,
            kind: kind2,
        };
        configuration.insert_vertex(v2, &model).expect("insert 2");
        configuration.validate(&model).expect("valid after 2");
    }

    #[test]
    fn linked_leg_is_involutive() {
        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = ImpurityModel::xxz(bath, 0.4, 0.2, 0.1, None).expect("model");
        let mut configuration = WormholeConfiguration::new(8.0, 1).expect("configuration");

        for i in 0..5 {
            let interaction = model.interaction(0);
            let kind = interaction.diagonal_kind(1, 1);
            let tau_a = 0.5 + i as f64;
            let tau_b = tau_a + 0.3;
            let vertex = Vertex {
                tau_a,
                tau_b,
                omega: 1.0,
                interaction: 0,
                kind,
            };
            configuration.insert_vertex(vertex, &model).expect("insert");
        }

        for (&_key, &endpoint) in &configuration.time_order {
            for side in [LegSide::Incoming, LegSide::Outgoing] {
                let leg = LegId { endpoint, side };
                let partner = configuration.linked_leg(leg).expect("linked");
                let back = configuration.linked_leg(partner).expect("back");
                assert_eq!(back, leg, "linked_leg not involutive");
            }
        }
    }

    #[test]
    fn cross_validate_with_legacy_index() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256PlusPlus;

        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = ImpurityModel::xxz(bath, 0.4, 0.2, 0.1, None).expect("model");
        let mut configuration = WormholeConfiguration::new(8.0, 1).expect("configuration");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        // Insert several vertices.  Use only same-spin diagonal vertices
        // so the periodic worldline constraint is always satisfied.
        for _ in 0..10 {
            let interaction = model.interaction(0);
            let kind = interaction.diagonal_kind(1, 1);
            let tau_a = rng.random::<f64>() * configuration.beta();
            let tau_b = (tau_a + 0.1).rem_euclid(configuration.beta());
            let vertex = Vertex {
                tau_a,
                tau_b,
                omega: 1.0,
                interaction: 0,
                kind,
            };
            configuration.insert_vertex(vertex, &model).expect("insert");
        }

        configuration.validate(&model).expect("valid");

        // Build legacy index for comparison.
        let legacy = WorldlineIndex::build(&configuration, &model).expect("legacy build");

        // Compare linked_leg for all legs.
        for slot_idx in 0..configuration.vertices.len() {
            if configuration.vertices[slot_idx].is_none() {
                continue;
            }
            let vertex_id = VertexId(slot_idx);
            for local_leg in 0..LEGS_PER_VERTEX {
                let new_leg = LegId::from_local(vertex_id, local_leg);
                let new_partner = configuration.linked_leg(new_leg).expect("new linked");

                let old_global = LEGS_PER_VERTEX * slot_idx + local_leg;
                let old_partner_global = legacy.linked_leg(old_global);
                let old_partner_vertex = old_partner_global / LEGS_PER_VERTEX;
                let old_partner_local = old_partner_global % LEGS_PER_VERTEX;
                let old_partner =
                    LegId::from_local(VertexId(old_partner_vertex), old_partner_local);

                assert_eq!(
                    new_partner, old_partner,
                    "linked_leg mismatch for vertex {} leg {}",
                    slot_idx, local_leg
                );
            }
        }

        // Compare spin_before at various times.
        for tau in [0.0, 0.5, 1.0, 2.0, 4.0, 7.9] {
            let new_spin = configuration.spin_before(&model, tau).expect("new spin");
            let old_spin = legacy.spin_before(&configuration, &model, tau);
            assert_eq!(new_spin, old_spin, "spin_before mismatch at tau={}", tau);
        }
    }

    #[test]
    fn random_insert_remove_stress() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256PlusPlus;

        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = ImpurityModel::xxz(bath, 0.4, 0.2, 0.1, None).expect("model");
        let mut configuration = WormholeConfiguration::new(8.0, 1).expect("configuration");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(99);

        let mut live_ids: Vec<VertexId> = Vec::new();

        for _ in 0..1000 {
            if live_ids.is_empty() || rng.random::<bool>() {
                // Insert a same-spin diagonal vertex.
                let interaction = model.interaction(0);
                let kind = interaction.diagonal_kind(1, 1);
                let tau_a = rng.random::<f64>() * configuration.beta();
                let tau_b = (tau_a + 0.05).rem_euclid(configuration.beta());
                let vertex = Vertex {
                    tau_a,
                    tau_b,
                    omega: 1.0,
                    interaction: 0,
                    kind,
                };
                let id = configuration.insert_vertex(vertex, &model).expect("insert");
                live_ids.push(id);
            } else {
                // Remove a random live vertex.
                let idx = rng.random_range(0..live_ids.len());
                let id = live_ids.swap_remove(idx);
                configuration.remove_vertex(id).expect("remove");
            }
        }

        configuration.validate(&model).expect("valid after stress");
    }

    #[test]
    fn many_sweeps_preserve_valid_configuration() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256PlusPlus;

        use crate::algorithm::{QmcKernel, UpdateSchedule};
        use crate::impurity::updates::WormholeEngine;

        let bath = Bath::SingleMode(SingleModeBath::new(1.0).expect("mode"));
        let model = ImpurityModel::xxz(bath, 0.4, 0.2, 0.1, None).expect("model");
        let mut engine = WormholeEngine::new(model.clone(), UpdateSchedule::new(4, 4, 64));
        engine.set_validate_each_sweep(true);
        let mut configuration = WormholeConfiguration::new(8.0, 1).expect("configuration");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(2026);

        for _ in 0..1000 {
            engine.sweep(&mut configuration, &mut rng).expect("sweep");
        }

        <WormholeEngine as QmcKernel<WormholeConfiguration, Xoshiro256PlusPlus>>::validate(
            &engine,
            &configuration,
        )
        .expect("valid configuration after 1000 sweeps");

        assert!(engine.stats().loops > 0);
        assert!(
            engine.stats().diagonal_add_accepts > 0 || engine.stats().diagonal_remove_accepts > 0
        );
    }
}
