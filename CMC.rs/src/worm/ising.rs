//! High-temperature graph representation of the ferromagnetic Ising model.
//!
//! For `K_e = beta J_e` and `t_e = tanh(K_e)`, the physical-sector weight is
//! proportional to `prod_e t_e^n_e` over even subgraphs. A worm state carries
//! two odd-degree defects at the ordered `(tail, head)` endpoints.

use super::{WormError, WormModel, WormSector, WormState, WormStepDelta, WormStepProposal};
use crate::lattice::graph::CsrLattice;
use rand::{Rng, RngExt};

/// Edge-occupation state and its incremental parity/log-weight caches.
#[derive(Debug, Clone, PartialEq)]
pub struct IsingGraphConfiguration {
    occupied: Vec<bool>,
    odd_parity: Vec<bool>,
    occupied_edges: usize,
    log_graph_weight: f64,
}

impl IsingGraphConfiguration {
    pub fn empty(n_sites: usize, n_edges: usize) -> Self {
        Self {
            occupied: vec![false; n_edges],
            odd_parity: vec![false; n_sites],
            occupied_edges: 0,
            log_graph_weight: 0.0,
        }
    }

    #[inline]
    pub fn occupied(&self) -> &[bool] {
        &self.occupied
    }

    #[inline]
    pub fn odd_parity(&self) -> &[bool] {
        &self.odd_parity
    }

    #[inline]
    pub const fn occupied_edges(&self) -> usize {
        self.occupied_edges
    }

    #[inline]
    pub const fn log_graph_weight(&self) -> f64 {
        self.log_graph_weight
    }

    pub fn from_occupied(
        model: &IsingGraphWormModel,
        occupied: Vec<bool>,
    ) -> Result<Self, WormError> {
        if occupied.len() != model.lattice.n_edges() {
            return Err(WormError::new("Ising graph checkpoint edge-count mismatch"));
        }
        let mut configuration = Self::empty(model.lattice.n_sites, occupied.len());
        for (edge_id, is_occupied) in occupied.into_iter().enumerate() {
            if is_occupied {
                let log_weight = model.log_edge_weights[edge_id];
                if !log_weight.is_finite() {
                    return Err(WormError::new(
                        "zero-weight Ising graph edge cannot be occupied",
                    ));
                }
                let edge = model.lattice.edges[edge_id];
                configuration.occupied[edge_id] = true;
                configuration.occupied_edges += 1;
                configuration.log_graph_weight += log_weight;
                configuration.odd_parity[edge.source] ^= true;
                configuration.odd_parity[edge.target] ^= true;
            }
        }
        Ok(configuration)
    }
}

/// One local graph-edge toggle from the current worm head.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IsingWormStep {
    pub edge_id: usize,
}

/// Reusable transactional patch for one graph-edge toggle.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct IsingGraphPatch {
    edge_id: usize,
    new_occupied: bool,
    log_weight_delta: f64,
}

/// Ferromagnetic Ising high-temperature graph model on an arbitrary multigraph.
#[derive(Debug, Clone)]
pub struct IsingGraphWormModel {
    lattice: CsrLattice,
    beta: f64,
    coupling: f64,
    edge_couplings: Vec<f64>,
    edge_tanh: Vec<f64>,
    log_edge_weights: Vec<f64>,
    occupied_energy_terms: Vec<f64>,
    base_energy: f64,
}

impl IsingGraphWormModel {
    pub fn new(lattice: CsrLattice, beta: f64, coupling: f64) -> Result<Self, WormError> {
        lattice.validate().map_err(WormError::new)?;
        ensure_single_component(&lattice)?;
        if !beta.is_finite() || beta < 0.0 {
            return Err(WormError::new(
                "Ising graph beta must be finite and non-negative",
            ));
        }
        if !coupling.is_finite() {
            return Err(WormError::new("Ising graph coupling must be finite"));
        }

        let mut edge_couplings = Vec::with_capacity(lattice.n_edges());
        let mut edge_tanh = Vec::with_capacity(lattice.n_edges());
        let mut log_edge_weights = Vec::with_capacity(lattice.n_edges());
        let mut occupied_energy_terms = Vec::with_capacity(lattice.n_edges());
        let mut base_energy = 0.0;

        for (edge_id, edge) in lattice.edges.iter().enumerate() {
            if edge.source == edge.target {
                return Err(WormError::new(format!(
                    "classical Ising worm does not support self-loop edge {edge_id}"
                )));
            }
            let edge_coupling = coupling * edge.weight;
            if !edge_coupling.is_finite() || edge_coupling < 0.0 {
                return Err(WormError::new(format!(
                    "classical Ising worm requires J * weight >= 0 on edge {edge_id}"
                )));
            }
            let t = (beta * edge_coupling).tanh();
            edge_couplings.push(edge_coupling);
            edge_tanh.push(t);
            log_edge_weights.push(if t == 0.0 { f64::NEG_INFINITY } else { t.ln() });
            base_energy -= edge_coupling * t;
            occupied_energy_terms.push(if t == 0.0 {
                0.0
            } else {
                -edge_coupling * (t.recip() - t)
            });
        }

        Ok(Self {
            lattice,
            beta,
            coupling,
            edge_couplings,
            edge_tanh,
            log_edge_weights,
            occupied_energy_terms,
            base_energy,
        })
    }

    #[inline]
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    #[inline]
    pub const fn beta(&self) -> f64 {
        self.beta
    }

    #[inline]
    pub const fn coupling(&self) -> f64 {
        self.coupling
    }

    #[inline]
    pub fn edge_couplings(&self) -> &[f64] {
        &self.edge_couplings
    }

    #[inline]
    pub fn edge_tanh(&self) -> &[f64] {
        &self.edge_tanh
    }

    pub fn empty_configuration(&self) -> IsingGraphConfiguration {
        IsingGraphConfiguration::empty(self.lattice.n_sites, self.lattice.n_edges())
    }

    /// Canonical energy estimator evaluated on a physical even graph.
    pub fn energy_estimator(&self, configuration: &IsingGraphConfiguration) -> f64 {
        self.base_energy
            + configuration
                .occupied
                .iter()
                .zip(&self.occupied_energy_terms)
                .filter_map(|(occupied, term)| occupied.then_some(*term))
                .sum::<f64>()
    }

    pub fn validate_configuration(
        &self,
        configuration: &IsingGraphConfiguration,
    ) -> Result<(), WormError> {
        if configuration.occupied.len() != self.lattice.n_edges() {
            return Err(WormError::new(
                "Ising graph occupation vector length mismatch",
            ));
        }
        if configuration.odd_parity.len() != self.lattice.n_sites {
            return Err(WormError::new("Ising graph parity vector length mismatch"));
        }

        let mut parity = vec![false; self.lattice.n_sites];
        let mut occupied_edges = 0usize;
        let mut log_graph_weight = 0.0;
        for (edge_id, &occupied) in configuration.occupied.iter().enumerate() {
            if !occupied {
                continue;
            }
            let log_weight = self.log_edge_weights[edge_id];
            if !log_weight.is_finite() {
                return Err(WormError::new(
                    "Ising graph contains an occupied zero-weight edge",
                ));
            }
            let edge = self.lattice.edges[edge_id];
            parity[edge.source] ^= true;
            parity[edge.target] ^= true;
            occupied_edges += 1;
            log_graph_weight += log_weight;
        }
        if parity != configuration.odd_parity {
            return Err(WormError::new("Ising graph parity cache mismatch"));
        }
        if occupied_edges != configuration.occupied_edges {
            return Err(WormError::new("Ising graph occupied-edge cache mismatch"));
        }
        let tolerance = 1e-12 * (1.0 + log_graph_weight.abs());
        if (log_graph_weight - configuration.log_graph_weight).abs() > tolerance {
            return Err(WormError::new("Ising graph log-weight cache mismatch"));
        }
        Ok(())
    }
}

impl WormModel for IsingGraphWormModel {
    type Configuration = IsingGraphConfiguration;
    type Defect = usize;
    type Step = IsingWormStep;
    type Patch = IsingGraphPatch;

    fn open_defect_count(&self, _configuration: &Self::Configuration) -> usize {
        self.lattice.n_sites
    }

    fn open_defect(
        &self,
        _configuration: &Self::Configuration,
        index: usize,
    ) -> Result<Self::Defect, WormError> {
        if index < self.lattice.n_sites {
            Ok(index)
        } else {
            Err(WormError::new("Ising worm opening site out of range"))
        }
    }

    fn propose_step(
        &self,
        state: &WormState<Self::Configuration, Self::Defect>,
        rng: &mut impl Rng,
    ) -> Result<Option<WormStepProposal<Self::Step>>, WormError> {
        let head = *state
            .head()
            .ok_or_else(|| WormError::new("Ising worm step requires a head"))?;
        let start = self.lattice.offsets[head];
        let end = self.lattice.offsets[head + 1];
        let degree = end - start;
        if degree == 0 {
            return Ok(None);
        }
        let incidence = rng.random_range(start..end);
        let edge_id = self.lattice.edge_ids[incidence];
        let neighbor = self.lattice.neighbors[incidence];
        let reverse_degree = self.lattice.offsets[neighbor + 1] - self.lattice.offsets[neighbor];
        let log_reverse_over_forward = (degree as f64).ln() - (reverse_degree as f64).ln();
        Ok(Some(WormStepProposal::new(
            IsingWormStep { edge_id },
            log_reverse_over_forward,
        )))
    }

    fn evaluate_step(
        &self,
        state: &WormState<Self::Configuration, Self::Defect>,
        step: &Self::Step,
        patch: &mut Self::Patch,
    ) -> Result<WormStepDelta<Self::Defect>, WormError> {
        let head = *state
            .head()
            .ok_or_else(|| WormError::new("Ising worm step requires a head"))?;
        let edge = self
            .lattice
            .edges
            .get(step.edge_id)
            .copied()
            .ok_or_else(|| WormError::new("Ising worm edge id out of range"))?;
        let new_head = edge
            .other(head)
            .ok_or_else(|| WormError::new("Ising worm edge is not incident to the head"))?;
        let new_occupied = !state.configuration().occupied[step.edge_id];
        let log_edge_weight = self.log_edge_weights[step.edge_id];
        let log_weight_ratio = if new_occupied {
            log_edge_weight
        } else {
            -log_edge_weight
        };
        *patch = IsingGraphPatch {
            edge_id: step.edge_id,
            new_occupied,
            log_weight_delta: log_weight_ratio,
        };
        Ok(WormStepDelta {
            new_head,
            log_weight_ratio,
        })
    }

    fn commit_step(
        &self,
        state: &mut WormState<Self::Configuration, Self::Defect>,
        step: &Self::Step,
        patch: &Self::Patch,
    ) {
        assert_eq!(patch.edge_id, step.edge_id, "stale Ising worm patch");
        assert!(
            patch.log_weight_delta.is_finite(),
            "accepted Ising worm patch must have finite weight"
        );
        let edge = self.lattice.edges[step.edge_id];
        let configuration = state.configuration_mut();
        assert_eq!(
            configuration.occupied[step.edge_id], !patch.new_occupied,
            "Ising worm patch does not match accepted state"
        );
        configuration.occupied[step.edge_id] = patch.new_occupied;
        if patch.new_occupied {
            configuration.occupied_edges += 1;
        } else {
            configuration.occupied_edges -= 1;
        }
        configuration.log_graph_weight += patch.log_weight_delta;
        configuration.odd_parity[edge.source] ^= true;
        configuration.odd_parity[edge.target] ^= true;
    }

    fn validate_state(
        &self,
        state: &WormState<Self::Configuration, Self::Defect>,
    ) -> Result<(), WormError> {
        state.validate_structure()?;
        self.validate_configuration(state.configuration())?;
        let mut expected = vec![false; self.lattice.n_sites];
        match state.sector() {
            WormSector::Physical => {}
            WormSector::Worm => {
                let head = *state
                    .head()
                    .ok_or_else(|| WormError::new("worm sector is missing its head"))?;
                let tail = *state
                    .tail()
                    .ok_or_else(|| WormError::new("worm sector is missing its tail"))?;
                if head >= self.lattice.n_sites || tail >= self.lattice.n_sites {
                    return Err(WormError::new("Ising worm endpoint out of range"));
                }
                if head != tail {
                    expected[head] = true;
                    expected[tail] = true;
                }
            }
        }
        if state.configuration().odd_parity != expected {
            return Err(WormError::new(
                "Ising graph parity does not match worm endpoints",
            ));
        }
        Ok(())
    }

    fn endpoint_bin_count(&self) -> usize {
        self.lattice.n_sites
    }

    fn endpoint_bin(&self, defect: &Self::Defect) -> Option<usize> {
        (*defect < self.lattice.n_sites).then_some(*defect)
    }
}

/// Reject multi-component lattices (including isolated sites).
///
/// The single defect pair of the worm diffuses within one connected
/// component of the bond graph. On a disconnected graph the other
/// components would keep their initial (empty-graph) occupation forever:
/// the walk would still look healthy while sampling a wrong ensemble — the
/// silent-garbage failure this check rules out at input time.
fn ensure_single_component(lattice: &CsrLattice) -> Result<(), WormError> {
    let mut seen = vec![false; lattice.n_sites];
    let mut queue = std::collections::VecDeque::new();
    seen[0] = true;
    queue.push_back(0usize);
    let mut reached = 1usize;
    while let Some(site) = queue.pop_front() {
        for (neighbor, _) in lattice.incidences(site) {
            if !seen[neighbor] {
                seen[neighbor] = true;
                reached += 1;
                queue.push_back(neighbor);
            }
        }
    }
    if reached != lattice.n_sites {
        return Err(WormError::new(format!(
            "classical Ising worm requires a connected (single-component) lattice: only {reached} \
             of {} sites are reachable from site 0, so the worm's defect pair would be confined \
             to one component and the remaining components would stay frozen at their initial \
             occupation, silently sampling a wrong ensemble",
            lattice.n_sites
        )));
    }
    Ok(())
}

/// Exact physical-sector high-temperature graph enumeration for small systems.
#[derive(Debug, Clone, PartialEq)]
pub struct ExactIsingGraphExpansion {
    pub log_reduced_partition: f64,
    pub mean_occupied_edges: f64,
    pub mean_energy: f64,
    pub edge_occupation_probabilities: Vec<f64>,
    pub physical_configurations: u64,
}

pub fn enumerate_ising_graph_expansion(
    model: &IsingGraphWormModel,
) -> Result<ExactIsingGraphExpansion, WormError> {
    let n_edges = model.lattice.n_edges();
    if n_edges > 24 {
        return Err(WormError::new(
            "exact Ising graph enumeration is limited to 24 edges",
        ));
    }
    let configurations = 1u64 << n_edges;
    let physical_log_weight = |mask: u64| {
        let mut parity = vec![false; model.lattice.n_sites];
        let mut log_weight = 0.0;
        for edge_id in 0..n_edges {
            if mask & (1u64 << edge_id) == 0 {
                continue;
            }
            let edge_log_weight = model.log_edge_weights[edge_id];
            if !edge_log_weight.is_finite() {
                return None;
            }
            let edge = model.lattice.edges[edge_id];
            parity[edge.source] ^= true;
            parity[edge.target] ^= true;
            log_weight += edge_log_weight;
        }
        parity.iter().all(|odd| !odd).then_some(log_weight)
    };

    let mut max_log_weight = f64::NEG_INFINITY;
    let mut physical_configurations = 0u64;
    for mask in 0..configurations {
        if let Some(log_weight) = physical_log_weight(mask) {
            max_log_weight = max_log_weight.max(log_weight);
            physical_configurations += 1;
        }
    }
    if physical_configurations == 0 {
        return Err(WormError::new(
            "exact Ising graph enumeration found no physical configuration",
        ));
    }

    let mut normalization = 0.0;
    let mut occupied_sum = 0.0;
    let mut energy_sum = 0.0;
    let mut edge_sums = vec![0.0; n_edges];
    for mask in 0..configurations {
        let Some(log_weight) = physical_log_weight(mask) else {
            continue;
        };
        let weight = (log_weight - max_log_weight).exp();
        normalization += weight;
        let occupied_edges = mask.count_ones() as usize;
        occupied_sum += weight * occupied_edges as f64;
        let mut energy = model.base_energy;
        for (edge_id, edge_sum) in edge_sums.iter_mut().enumerate() {
            if mask & (1u64 << edge_id) != 0 {
                *edge_sum += weight;
                energy += model.occupied_energy_terms[edge_id];
            }
        }
        energy_sum += weight * energy;
    }

    Ok(ExactIsingGraphExpansion {
        log_reduced_partition: max_log_weight + normalization.ln(),
        mean_occupied_edges: occupied_sum / normalization,
        mean_energy: energy_sum / normalization,
        edge_occupation_probabilities: edge_sums
            .into_iter()
            .map(|sum| sum / normalization)
            .collect(),
        physical_configurations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::build_chain;

    #[test]
    fn validation_detects_parity_and_log_weight_cache_corruption() {
        let model = IsingGraphWormModel::new(build_chain(4, true), 0.4, 1.0).unwrap();
        let mut configuration = model.empty_configuration();
        configuration.odd_parity[0] = true;
        assert!(model.validate_configuration(&configuration).is_err());

        let mut configuration = model.empty_configuration();
        configuration.log_graph_weight = 1.0;
        assert!(model.validate_configuration(&configuration).is_err());
    }
}
