//! Multi-component Ising graph worm ensemble.
//!
//! The Ising high-temperature graph expansion factorizes over the connected
//! components of the bond graph: `Z = prod_c Z_c` with independent even-subgraph
//! ensembles per component. A single two-defect worm diffuses within one
//! component only, so a disconnected lattice is sampled correctly by running
//! one independent worm (one defect pair, one domain-separated RNG stream) per
//! component and combining observables additively.
//!
//! Per-sweep RNG derivation: one `u64` salt per component is drawn from the
//! shared (checkpoint-owned) context stream, and the component's stream is
//! derived through [`carlo_rs::RngStreamKey`] with the component index in the
//! replica field. This keeps the streams domain-separated while leaving no RNG
//! state outside the checkpointed context — a restored run replays the exact
//! same per-component streams.

use super::{IsingGraphWormModel, WormConfig, WormError, WormKernel};
use crate::lattice::graph::{Bond, CsrLattice};
use rand::{Rng, RngExt};

/// One connected component's independent two-defect worm.
pub struct IsingComponentWorm {
    kernel: WormKernel<IsingGraphWormModel>,
    /// `global_sites[local]` is the global site index of local site `local`.
    global_sites: Vec<usize>,
}

impl IsingComponentWorm {
    #[inline]
    pub const fn kernel(&self) -> &WormKernel<IsingGraphWormModel> {
        &self.kernel
    }

    #[inline]
    pub fn kernel_mut(&mut self) -> &mut WormKernel<IsingGraphWormModel> {
        &mut self.kernel
    }

    #[inline]
    pub const fn model(&self) -> &IsingGraphWormModel {
        self.kernel.model()
    }

    /// Global site indices of this component, ascending.
    #[inline]
    pub fn global_sites(&self) -> &[usize] {
        &self.global_sites
    }

    /// Local (component) index of a global site, if it belongs here.
    #[inline]
    pub fn local_site(&self, global: usize) -> Option<usize> {
        self.global_sites.binary_search(&global).ok()
    }
}

/// Independent per-component worm dynamics on a possibly disconnected lattice.
pub struct IsingGraphWormEnsemble {
    /// The full lattice the ensemble was decomposed from.
    lattice: CsrLattice,
    beta: f64,
    coupling: f64,
    config: WormConfig,
    components: Vec<IsingComponentWorm>,
}

impl IsingGraphWormEnsemble {
    /// Decompose `lattice` into connected components and build one validated
    /// worm model/kernel pair per component.
    ///
    /// Isolated sites form trivial components whose only even subgraph is the
    /// empty graph; their worm harmlessly opens and closes in place. Every
    /// genuinely invalid input (empty lattice, non-finite `beta`/`coupling`,
    /// negative `coupling * weight`, self-loops) is still rejected loudly by
    /// the per-component model constructor.
    pub fn new(
        lattice: CsrLattice,
        beta: f64,
        coupling: f64,
        config: WormConfig,
    ) -> Result<Self, WormError> {
        config.validate()?;
        lattice.validate().map_err(WormError::new)?;
        let mut components = Vec::new();
        for sites in lattice.connected_components() {
            let sub_lattice = component_sublattice(&lattice, &sites);
            let model = IsingGraphWormModel::new(sub_lattice, beta, coupling)?;
            let configuration = model.empty_configuration();
            let kernel = WormKernel::new(model, configuration, config.clone())?;
            components.push(IsingComponentWorm {
                kernel,
                global_sites: sites,
            });
        }
        debug_assert_eq!(
            components
                .iter()
                .map(|component: &IsingComponentWorm| component.global_sites.len())
                .sum::<usize>(),
            lattice.n_sites,
            "components must partition all sites"
        );
        Ok(Self {
            lattice,
            beta,
            coupling,
            config,
            components,
        })
    }

    /// Number of connected components being sampled.
    #[inline]
    pub fn n_components(&self) -> usize {
        self.components.len()
    }

    /// The full lattice this ensemble samples.
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
    pub const fn config(&self) -> &WormConfig {
        &self.config
    }

    #[inline]
    pub fn components(&self) -> &[IsingComponentWorm] {
        &self.components
    }

    #[inline]
    pub fn components_mut(&mut self) -> &mut [IsingComponentWorm] {
        &mut self.components
    }

    /// The kernel of the sole component; use [`Self::components`] for
    /// multi-component lattices.
    #[inline]
    pub fn single_kernel(&self) -> &WormKernel<IsingGraphWormModel> {
        assert_eq!(
            self.components.len(),
            1,
            "kernel() serves a single component; a {}-component lattice has one kernel each — \
             use components()/components_mut()",
            self.components.len()
        );
        &self.components[0].kernel
    }

    /// Mutable counterpart of [`Self::single_kernel`].
    #[inline]
    pub fn single_kernel_mut(&mut self) -> &mut WormKernel<IsingGraphWormModel> {
        assert_eq!(
            self.components.len(),
            1,
            "kernel_mut() serves a single component; a {}-component lattice has one kernel each \
             — use components_mut()",
            self.components.len()
        );
        &mut self.components[0].kernel
    }

    /// Whether every component is currently in the physical sector.
    ///
    /// Because the joint extended-sector weight factorizes over components,
    /// conditioning on "all components physical" leaves the product of the
    /// physical even-subgraph ensembles — the correct full-graph measure for
    /// extensive observables.
    #[inline]
    pub fn all_physical(&self) -> bool {
        self.components
            .iter()
            .all(|component| component.kernel.state().is_physical())
    }

    /// Total energy and occupied-edge count, valid only when all components
    /// are physical ([`Self::all_physical`]).
    pub fn total_energy_and_occupied_edges(&self) -> Option<(f64, usize)> {
        if !self.all_physical() {
            return None;
        }
        let mut energy = 0.0;
        let mut occupied = 0usize;
        for component in &self.components {
            let state = component.kernel.state();
            energy += component
                .kernel
                .model()
                .energy_estimator(state.configuration());
            occupied += state.configuration().occupied_edges();
        }
        Some((energy, occupied))
    }

    /// One sweep of every component on domain-separated derived streams.
    ///
    /// Each component draws one salt from `rng` and derives its own stream via
    /// [`carlo_rs::RngStreamKey`] (component index in the replica field), so
    /// per-component randomness stays domain-separated while all state remains
    /// a deterministic function of the caller's stream — nothing is hidden
    /// from a checkpoint.
    pub fn sweep(&mut self, rng: &mut impl Rng) -> Result<(), WormError> {
        for (index, component) in self.components.iter_mut().enumerate() {
            let salt = rng.random::<u64>();
            let mut stream = carlo_rs::RngStreamKey::new(salt)
                .with_replica(index as u64)
                .seeded::<rand_xoshiro::Xoshiro256PlusPlus>();
            component.kernel.sweep(&mut stream)?;
        }
        Ok(())
    }

    /// Audit every component kernel.
    pub fn validate(&self) -> Result<(), WormError> {
        for component in &self.components {
            component.kernel.validate()?;
        }
        Ok(())
    }

    /// Worm-estimated two-point correlation `⟨s_tail s_head⟩`.
    ///
    /// Defined only when both sites belong to the same component: a defect
    /// pair never spans components. Cross-component Ising correlations
    /// factorize to `⟨s_i⟩⟨s_j⟩ = 0` and have no worm estimator.
    pub fn endpoint_correlation(&self, tail: usize, head: usize) -> Option<f64> {
        for component in &self.components {
            let Some(local_tail) = component.local_site(tail) else {
                continue;
            };
            let local_head = component.local_site(head)?;
            let histogram = component.kernel.endpoint_pairs()?;
            return histogram.correlation_ratio(local_tail, local_head);
        }
        None
    }
}

/// Extract one component's sub-lattice with local site labels `0..len`.
///
/// Edges are emitted in global edge order, so a single-component lattice maps
/// to itself with identical edge ids. Bond types and weights are preserved.
fn component_sublattice(lattice: &CsrLattice, sites: &[usize]) -> CsrLattice {
    let mut local_of = vec![usize::MAX; lattice.n_sites];
    for (local, &global) in sites.iter().enumerate() {
        local_of[global] = local;
    }
    let edges: Vec<Bond> = lattice
        .edges
        .iter()
        .filter(|edge| local_of[edge.source] != usize::MAX)
        .map(|edge| {
            Bond::new(
                local_of[edge.source],
                local_of[edge.target],
                edge.kind,
                edge.weight,
            )
        })
        .collect();
    CsrLattice::from_edges(sites.len(), edges)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::build_chain;
    use rand::SeedableRng;

    fn config() -> WormConfig {
        WormConfig {
            local_updates_per_sweep: 4,
            close_probability: 0.25,
            log_worm_fugacity: -1.0,
            track_endpoint_pairs: false,
            cache_audit_interval: 0,
        }
    }

    #[test]
    fn ensemble_runs_multi_component_lattices_per_component() {
        // 4-ring + 3-chain + isolated site: three independent components.
        let lattice = CsrLattice::from_edges(
            8,
            vec![
                crate::Bond::new(0, 1, crate::BondType::Generic, 1.0),
                crate::Bond::new(1, 2, crate::BondType::Generic, 1.0),
                crate::Bond::new(2, 3, crate::BondType::Generic, 1.0),
                crate::Bond::new(3, 0, crate::BondType::Generic, 1.0),
                crate::Bond::new(4, 5, crate::BondType::Generic, 1.0),
                crate::Bond::new(5, 6, crate::BondType::Generic, 1.0),
            ],
        );
        let mut ensemble = IsingGraphWormEnsemble::new(lattice, 0.44, 1.0, config()).unwrap();
        assert_eq!(ensemble.n_components(), 3);
        assert_eq!(ensemble.components()[0].global_sites(), &[0, 1, 2, 3]);
        assert_eq!(ensemble.components()[1].global_sites(), &[4, 5, 6]);
        assert_eq!(ensemble.components()[2].global_sites(), &[7]);
        // Every edge is carried by exactly one component sub-lattice.
        let edges: usize = ensemble
            .components()
            .iter()
            .map(|component| component.model().lattice().n_edges())
            .sum();
        assert_eq!(edges, 6);

        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0xE0);
        for _ in 0..2_000 {
            ensemble.sweep(&mut rng).unwrap();
        }
        ensemble.validate().unwrap();
        for component in ensemble.components() {
            assert!(component.kernel.statistics().completed_worms > 0);
        }
        // Cross-component correlation has no estimator; same-component has.
        assert!(ensemble.endpoint_correlation(0, 2).is_none()); // no histogram tracking
        assert!(ensemble.endpoint_correlation(0, 4).is_none()); // different components
    }

    #[test]
    fn ensemble_rejects_genuinely_invalid_input_loudly() {
        let disconnected = CsrLattice::from_edges(
            4,
            vec![
                crate::Bond::new(0, 1, crate::BondType::Generic, 1.0),
                crate::Bond::new(2, 3, crate::BondType::Generic, 1.0),
            ],
        );
        // Non-finite / negative beta and coupling must error, not panic.
        assert!(
            IsingGraphWormEnsemble::new(disconnected.clone(), f64::NAN, 1.0, config()).is_err()
        );
        assert!(IsingGraphWormEnsemble::new(disconnected.clone(), -0.1, 1.0, config()).is_err());
        assert!(
            IsingGraphWormEnsemble::new(disconnected.clone(), 0.4, f64::INFINITY, config())
                .is_err()
        );
        // Antiferromagnetic coupling on any edge is rejected per component.
        assert!(IsingGraphWormEnsemble::new(disconnected, 0.4, -1.0, config()).is_err());
        // Self-loops remain unsupported.
        let looped = CsrLattice::from_edges(
            3,
            vec![
                crate::Bond::new(0, 1, crate::BondType::Generic, 1.0),
                crate::Bond::new(2, 2, crate::BondType::Generic, 1.0),
            ],
        );
        assert!(IsingGraphWormEnsemble::new(looped, 0.4, 1.0, config()).is_err());
        // Zero sites cannot form a valid lattice at all.
        let empty = CsrLattice::try_from_edges(0, vec![]);
        assert!(empty.is_err());
    }

    #[test]
    fn ensemble_domain_separated_streams_are_deterministic() {
        // Same shared stream -> identical ensemble trajectory; a connected
        // lattice is a one-component special case.
        let lattice = build_chain(4, true);
        let mut left = IsingGraphWormEnsemble::new(lattice.clone(), 0.4, 1.0, config()).unwrap();
        let mut right = IsingGraphWormEnsemble::new(lattice, 0.4, 1.0, config()).unwrap();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x5EED);
        let mut replay = rng.clone();
        for _ in 0..500 {
            left.sweep(&mut rng).unwrap();
            right.sweep(&mut replay).unwrap();
        }
        let left_occupied: Vec<_> = left
            .components()
            .iter()
            .map(|c| c.kernel().state().configuration().occupied().to_vec())
            .collect();
        let right_occupied: Vec<_> = right
            .components()
            .iter()
            .map(|c| c.kernel().state().configuration().occupied().to_vec())
            .collect();
        assert_eq!(left_occupied, right_occupied);
    }
}
