//! Model-side abstractions for classical lattice Monte Carlo.
//!
//! The core contract is deliberately temperature independent: a Hamiltonian
//! returns physical energies, while algorithms/ensembles apply `beta` exactly
//! once.  Physical bonds come from [`crate::lattice::CsrLattice::edges`], so
//! weighted and parallel interactions are counted without hidden `/ 2` rules.

use crate::lattice::{Bond, CsrLattice};
use rand::Rng;
use smallvec::SmallVec;

/// Small-vector representation used by the compatibility API.
///
/// Dimensions above three are supported and spill to the heap.  Built-in
/// `ONModel<D>` therefore works for arbitrary `D`, while XY/Heisenberg remain
/// allocation-free.
pub type Spin = SmallVec<[f64; 3]>;

/// General physical-energy model.
///
/// Custom multi-site/factor models may implement this trait directly.  The
/// common onsite + pair-interaction case should implement [`PairInteraction`]
/// instead and receives an optimized, correct blanket implementation.
pub trait Hamiltonian: Send + Sync {
    fn spin_dim(&self) -> usize;
    fn coupling(&self) -> f64;

    /// Energy terms affected by replacing `site` with `proposed`.
    ///
    /// `beta` is retained for source compatibility.  The return value must be
    /// a physical energy and must not contain a temperature factor.
    fn local_energy(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        beta: f64,
        proposed: &[f64],
    ) -> f64;

    fn compute_total_energy(&self, spins: &[f64], lattice: &CsrLattice, beta: f64) -> f64;

    fn delta_energy(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        proposed: &[f64],
    ) -> f64 {
        let spin_dim = self.spin_dim();
        let old = &spins[site * spin_dim..(site + 1) * spin_dim];
        self.local_energy(spins, lattice, site, 1.0, proposed)
            - self.local_energy(spins, lattice, site, 1.0, old)
    }
}

/// Convenience interface for onsite plus physical pair-bond Hamiltonians.
pub trait PairInteraction: Send + Sync {
    fn spin_dim(&self) -> usize;
    fn coupling(&self) -> f64;
    fn bond_energy(&self, left: &[f64], right: &[f64], bond: &Bond) -> f64;

    fn onsite_energy(&self, _site: usize, _spin: &[f64]) -> f64 {
        0.0
    }
}

impl<T: PairInteraction> Hamiltonian for T {
    fn spin_dim(&self) -> usize {
        PairInteraction::spin_dim(self)
    }

    fn coupling(&self) -> f64 {
        PairInteraction::coupling(self)
    }

    fn local_energy(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        _beta: f64,
        proposed: &[f64],
    ) -> f64 {
        let spin_dim = PairInteraction::spin_dim(self);
        debug_assert_eq!(proposed.len(), spin_dim);
        let mut energy = self.onsite_energy(site, proposed);
        let mut previous_self_loop = None;

        for (_neighbor, edge_id) in lattice.incidences(site) {
            let edge = lattice.edges[edge_id];
            if edge.source == site && edge.target == site {
                if previous_self_loop == Some(edge_id) {
                    continue;
                }
                previous_self_loop = Some(edge_id);
                energy += self.bond_energy(proposed, proposed, &edge);
                continue;
            }

            let (left, right) = if edge.source == site {
                let base = edge.target * spin_dim;
                (proposed, &spins[base..base + spin_dim])
            } else {
                let base = edge.source * spin_dim;
                (&spins[base..base + spin_dim], proposed)
            };
            energy += self.bond_energy(left, right, &edge);
        }
        energy
    }

    fn compute_total_energy(&self, spins: &[f64], lattice: &CsrLattice, _beta: f64) -> f64 {
        let spin_dim = PairInteraction::spin_dim(self);
        assert_eq!(
            spins.len(),
            lattice.n_sites * spin_dim,
            "spin buffer length does not match lattice and spin dimension"
        );

        let onsite: f64 = (0..lattice.n_sites)
            .map(|site| {
                let base = site * spin_dim;
                self.onsite_energy(site, &spins[base..base + spin_dim])
            })
            .sum();

        onsite
            + lattice
                .edges
                .iter()
                .map(|edge| {
                    let left_base = edge.source * spin_dim;
                    let right_base = edge.target * spin_dim;
                    self.bond_energy(
                        &spins[left_base..left_base + spin_dim],
                        &spins[right_base..right_base + spin_dim],
                        edge,
                    )
                })
                .sum::<f64>()
    }
}

/// Independent initialization capability.
///
/// This is intentionally separate from Metropolis proposals: a model may be
/// initialized by Carlo.rs without pretending to support a particular update.
pub trait Initializable: Hamiltonian {
    fn random_spin(&self, rng: &mut impl Rng) -> Spin;

    fn ordered_spin(&self) -> Spin;
}

/// Symmetric model-level proposal used by [`StandardStrategy`](crate::StandardStrategy).
pub trait Proposable: Hamiltonian {
    /// Compatibility proposal that does not inspect the current spin.
    fn propose(&self, rng: &mut impl Rng) -> Spin;

    /// State-aware symmetric proposal.  Discrete models override this to
    /// avoid no-op proposals; existing custom models inherit `propose`.
    fn propose_from(&self, _current: &[f64], rng: &mut impl Rng) -> Spin {
        self.propose(rng)
    }

    fn normalize_spin(&self, _spin: &mut [f64]) {}
}

/// Magnetization or another model-native scalar order parameter.
pub trait Measurable: Hamiltonian {
    fn magnetization(&self, spins: &[f64]) -> f64;
}

/// Auxiliary state shared by cluster construction and transformation.
#[derive(Debug, Clone)]
pub enum ClusterAuxiliary {
    /// Discrete models do not need a global embedding variable for bonds.
    None,
    /// Assign every spin in a cluster to one discrete state.
    DiscreteTarget(f64),
    /// Reflection normal used by O(N) embedding clusters.
    Reflection(Spin),
    /// Leave this SW cluster unchanged.
    Identity,
}

/// Cluster-update policy implemented by models that support Wolff/SW.
///
/// Unlike the previous trait, no method has a panic default and bond
/// probabilities receive both endpoint spins, the physical bond and the
/// embedding auxiliary.  This is required for correct O(N) clusters.
pub trait ClusterModel: Hamiltonian {
    /// Auxiliary/target for one Wolff cluster, sampled from its seed spin.
    fn wolff_auxiliary(&self, seed_spin: &[f64], rng: &mut impl Rng) -> ClusterAuxiliary;

    /// Global auxiliary used while forming SW bonds.
    fn sw_bond_auxiliary(&self, rng: &mut impl Rng) -> ClusterAuxiliary;

    /// Independent transformation for one completed SW cluster.
    fn sw_cluster_auxiliary(
        &self,
        representative_spin: &[f64],
        bond_auxiliary: &ClusterAuxiliary,
        rng: &mut impl Rng,
    ) -> ClusterAuxiliary;

    /// Activation probability for one physical edge.
    fn cluster_bond_probability(
        &self,
        left: &[f64],
        right: &[f64],
        bond: &Bond,
        auxiliary: &ClusterAuxiliary,
        beta: f64,
    ) -> f64;

    /// Apply a completed cluster transformation to one site spin.
    fn transform_cluster_spin(&self, spin: &[f64], auxiliary: &ClusterAuxiliary) -> Spin;
}

/// Models with a linear local field, used by exact over-relaxation.
pub trait LocalFieldModel: Hamiltonian {
    fn local_field(&self, spins: &[f64], lattice: &CsrLattice, site: usize, output: &mut [f64]);
}

/// Exact discrete conditional sampler.
pub trait HeatBathable: Hamiltonian {
    fn heat_bath_sample_site(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        beta: f64,
        rng: &mut impl Rng,
    ) -> Spin;
}

/// Exact continuous conditional sampler.
pub trait ContinuousHeatBathable: Hamiltonian {
    fn heat_bath_sample_site(
        &self,
        spins: &[f64],
        lattice: &CsrLattice,
        site: usize,
        beta: f64,
        rng: &mut impl Rng,
    ) -> Spin;
}
