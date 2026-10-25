//! q = 2 Edwards-Sokal coupling: the actual Fortuin-Kasteleyn bond draw,
//! exposed and analyzed (roadmap §15 R3).
//!
//! The joint Edwards-Sokal measure couples a zero-field ferromagnetic Ising
//! spin configuration `s` with a bond subset `A` (an FK configuration):
//!
//! ```text
//! P(s, A) ~ prod_{e in A}     p_e * [s_u(e) = s_v(e)]
//!         prod_{e not in A}   (1 - p_e),        p_e = 1 - exp(-2 beta J_e),
//! ```
//!
//! with `J_e > 0` the edge coupling (`CsrLattice` [`Bond`](crate::Bond)
//! weight) and `beta >= 0` the inverse temperature. Alternating conditional
//! draws are Swendsen-Wang: given `s`, each **same-spin** edge opens
//! independently with probability `p_e` (opposite-spin edges are closed
//! with probability one); given `A`, each FK cluster receives one
//! independent uniform spin. The marginal of `A` is the q = 2 random-cluster
//! measure `P(A) ~ prod p prod (1 - p) 2^{k(A)}`, and the Edwards-Sokal
//! identity `P_FK(i connected to j) = <s_i s_j>` holds at every step.
//!
//! This module is the *sampler* side of the R track; R1
//! ([`BondSample`], [`analyze_bonds`]) is the analysis contract. One
//! `sweep()` performs one Edwards-Sokal step and stores the ACTUAL bond draw
//! as a [`BondSample`] — never a reconstruction — and the following
//! `measure()` runs the shared F3 [`analyze`] on that same sample through the
//! [`ObservablePlan`], exactly like [`StaticPercolationMC`]. A frozen-spins
//! mode fixes `s` forever, reducing the sweep to the pure conditional FK
//! draw (the exact-enumeration gates build on this).
//!
//! # The worm's even subgraph is NOT an FK configuration
//!
//! The worm module (`crate::worm`, `worm/ising.rs`) samples the Ising
//! high-temperature **graph representation**: edge subsets whose every vertex
//! has even degree, weighted by `prod_e tanh(beta J_e)^{n_e}` (two odd
//! defects in the worm sector). That is a different measure on a different
//! support — even subgraphs with tanh weights, no independent per-edge
//! Bernoulli law, no `2^{k(A)}` cluster factor — and it must never be
//! relabeled, re-exported, or analyzed as an FK configuration. The FK draw
//! here is unrestricted over same-spin-compatible subsets with
//! `1 - exp(-2 beta J_e)` weights and per-cluster spin symmetry. The
//! physics tests verify the two representations agree only at the level of
//! Ising observables (partition function, correlations), never
//! configuration by configuration.
//!
//! # RNG contract (dense order, bitwise documented)
//!
//! One unfrozen `sweep()` consumes variates in this exact order:
//!
//! 1. **FK bond pass** — physical edges in ascending dense `EdgeId` order.
//!    An edge whose endpoint spins differ is closed and consumes NOTHING.
//!    A same-spin edge with `0 < p_e < 1` consumes exactly one `f64`
//!    (`open iff u < p_e`). Exact endpoints skip the draw (mirroring the F2
//!    law contract): `p_e = 0` (e.g. `beta = 0`) is closed, `p_e = 1`
//!    (reachable only through `exp` underflow) is open.
//! 2. **Cluster spin pass** (skipped when frozen) — sites in ascending
//!    index order; the first site of each FK cluster (its lowest-index
//!    member, a cluster-set invariant) consumes exactly one `f64`
//!    (`+1 if u < 1/2 else -1`) for the whole cluster.
//!
//! # References
//!
//! Edwards & Sokal, Phys. Rev. D **38**, 2009 (1988) and Swendsen & Wang,
//! Phys. Rev. Lett. **58**, 86 (1987) pre-date arXiv, so no PDFs are
//! archivable (the Tarjan-1972 precedent); the exact double-sum identities —
//! `Z_Ising = exp(beta sum_e J_e) * sum_s sum_{A compatible} prod p prod
//! (1-p)` and `P_FK(i ~ j) = <s_i s_j>` — are the anchors, verified by
//! exhaustive enumeration in `tests/physics/edwards_sokal.rs`.

use core::cmp::Ordering;
use core::fmt;

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{
    analyze, analyze_with_labels, BondSample, ComponentWorkspace, ObservablePlan,
    ObservablePlanError, SampleView, VertexActivity,
};
use crate::{CsrLattice, EdgeId, GraphView, UndirectedGraphView, VertexId};

/// All-active vertex activity: the vertex half of bond percolation semantics.
///
/// An FK configuration is an edge subset, so every vertex of the reference
/// graph is tracked (isolates included) — exactly R1's vertex view, as a
/// standalone source so an [`EdwardsSokalMC`] (or any caller holding a
/// [`BondSample`]) can feed the shared analyzers without materializing a
/// dense vertex mask.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct AllActiveVertices {
    vertex_count: usize,
}

impl AllActiveVertices {
    /// One active vertex per vertex of the reference graph.
    pub const fn new(vertex_count: usize) -> Self {
        Self { vertex_count }
    }
}

impl VertexActivity for AllActiveVertices {
    #[inline]
    fn vertex_count(&self) -> usize {
        self.vertex_count
    }

    #[inline]
    fn vertex_active(&self, _vertex: VertexId) -> bool {
        true
    }
}

/// Edwards-Sokal / Swendsen-Wang adapter exposing its actual FK bond draw.
///
/// Owns a validated `CsrLattice` with strictly ferromagnetic couplings
/// (`J_e = Bond.weight > 0`, finite — typed rejection otherwise) and the
/// inverse temperature `beta >= 0` (finite). `sweep()` is one ES step whose
/// bond half is retained as [`Self::last_bond_sample`]; `measure()` analyzes
/// that sample with the preallocated workspace and records the planned F4
/// scalars. Sampling is a Markov chain (unlike the i.i.d. static family):
/// thermalize before measuring. Winding observables are a typed rejection —
/// FK winding analysis awaits an embedding-aware follow-up.
///
/// The frozen variant ([`Self::frozen`]) pins a caller-supplied spin
/// configuration, so every sweep redraws only the conditional FK bond law
/// (independent Bernoulli per same-spin edge) — the deterministic ensemble
/// the exact-enumeration gates exhaust.
pub struct EdwardsSokalMC {
    lattice: CsrLattice,
    beta: f64,
    /// Precomputed `p_e = 1 - exp(-2 beta J_e)` per dense physical edge.
    bond_probabilities: Vec<f64>,
    /// Current Ising spins, `+1`/`-1` per site.
    spins: Vec<i8>,
    /// Frozen mode: never reassign spins (pure conditional FK draws).
    frozen: bool,
    /// Union-find over FK clusters, rebuilt every bond pass.
    parent: Vec<usize>,
    rank: Vec<u8>,
    /// Cluster spin by root; `0` marks "not drawn this sweep".
    cluster_spin: Vec<i8>,
    /// Open edges of the current draw, ascending dense order (no duplicates).
    open_edges: Vec<EdgeId>,
    /// The exposed configuration: the actual draw, rebuilt every sweep.
    sample: BondSample,
    /// All-active vertex source paired with `sample` for the analyzers.
    vertices: AllActiveVertices,
    workspace: ComponentWorkspace,
    plan: ObservablePlan,
}

impl EdwardsSokalMC {
    /// Full Edwards-Sokal dynamics (Swendsen-Wang), spins initialized all-up.
    pub fn new(
        lattice: CsrLattice,
        beta: f64,
        plan: ObservablePlan,
    ) -> Result<Self, EdwardsSokalError> {
        let vertices = lattice.vertex_count();
        Self::build(lattice, beta, vec![1; vertices], false, plan)
    }

    /// Frozen-spins mode: `spins` (`+1`/`-1` per site) never changes, so each
    /// sweep is one draw of the conditional FK bond law alone.
    pub fn frozen(
        lattice: CsrLattice,
        beta: f64,
        spins: impl Into<Vec<i8>>,
        plan: ObservablePlan,
    ) -> Result<Self, EdwardsSokalError> {
        Self::build(lattice, beta, spins.into(), true, plan)
    }

    fn build(
        lattice: CsrLattice,
        beta: f64,
        spins: Vec<i8>,
        frozen: bool,
        plan: ObservablePlan,
    ) -> Result<Self, EdwardsSokalError> {
        lattice
            .validate()
            .map_err(|error| EdwardsSokalError::Topology(error.to_string()))?;
        if !beta.is_finite() || beta < 0.0 {
            return Err(EdwardsSokalError::InverseTemperature { beta });
        }
        let mut bond_probabilities = Vec::with_capacity(lattice.edge_count());
        for (edge, bond) in lattice.edges.iter().enumerate() {
            // The FK probability 1 - exp(-2 beta J) is a probability only for
            // finite J > 0; antiferromagnetic or vanishing couplings are a
            // different (or undefined) ensemble and are rejected loudly.
            let coupling = bond.weight;
            if !coupling.is_finite() || coupling <= 0.0 {
                return Err(EdwardsSokalError::Coupling {
                    edge,
                    weight: coupling,
                });
            }
            // beta, coupling finite and non-negative: -2 beta J is finite or
            // -inf (overflow), exp(...) in [0, 1), p in [0, 1].
            bond_probabilities.push(1.0 - (-2.0 * beta * coupling).exp());
        }
        if spins.len() != lattice.vertex_count() {
            return Err(EdwardsSokalError::SpinCount {
                expected: lattice.vertex_count(),
                actual: spins.len(),
            });
        }
        for (site, &spin) in spins.iter().enumerate() {
            if spin != -1 && spin != 1 {
                return Err(EdwardsSokalError::SpinValue { site, value: spin });
            }
        }
        if plan.requests_winding() {
            return Err(EdwardsSokalError::WindingNotSupported);
        }
        plan.validate(lattice.vertex_count())
            .map_err(EdwardsSokalError::Observable)?;
        for (index, query) in plan.boundary_queries().iter().enumerate() {
            if query.vertex_count() != lattice.vertex_count() {
                return Err(EdwardsSokalError::QueryDomain {
                    query: index,
                    expected: lattice.vertex_count(),
                    actual: query.vertex_count(),
                });
            }
        }
        let mut workspace = ComponentWorkspace::new();
        workspace
            .prepare(lattice.vertex_count(), plan.boundary_queries().len(), false)
            .map_err(|_| EdwardsSokalError::Capacity)?;
        let vertices = lattice.vertex_count();
        // Pre-sweep state: the legal all-closed sample of this edge domain
        // (not a default with a foreign domain).
        let sample = BondSample::from_open_edges(&lattice, std::iter::empty())
            .expect("empty open-edge set is always in-domain");
        Ok(Self {
            vertices: AllActiveVertices::new(vertices),
            lattice,
            beta,
            bond_probabilities,
            spins,
            frozen,
            parent: Vec::new(),
            rank: Vec::new(),
            cluster_spin: Vec::new(),
            open_edges: Vec::new(),
            sample,
            workspace,
            plan,
        })
    }

    /// The owned topology with its ferromagnetic couplings.
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// The inverse temperature (validated finite, `>= 0`).
    pub const fn beta(&self) -> f64 {
        self.beta
    }

    /// The per-edge FK open probabilities `1 - exp(-2 beta J_e)`, dense order.
    pub fn bond_probabilities(&self) -> &[f64] {
        &self.bond_probabilities
    }

    /// Whether this adapter runs in frozen-spins mode.
    pub const fn is_frozen(&self) -> bool {
        self.frozen
    }

    /// The current spin configuration (`+1`/`-1` per site).
    pub fn spins(&self) -> &[i8] {
        &self.spins
    }

    /// The observable plan recorded by `measure()`.
    pub const fn observable_plan(&self) -> &ObservablePlan {
        &self.plan
    }

    /// The ACTUAL FK bond draw of the most recent `sweep()`.
    ///
    /// This is the configuration `measure()` analyzes — the R3 contract: the
    /// exposed sample and the analyzed sample are the same object, never a
    /// reconstruction. Before the first sweep it is the empty sample of the
    /// correct edge domain.
    pub const fn last_bond_sample(&self) -> &BondSample {
        &self.sample
    }

    /// The conditional FK bond pass given the current spins (contract step 1)
    /// plus union-find cluster construction and sample exposure.
    fn draw_fk_bonds<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        let Self {
            lattice,
            spins,
            bond_probabilities,
            parent,
            rank,
            open_edges,
            ..
        } = self;
        let vertex_count = lattice.vertex_count();
        parent.clear();
        parent.extend(0..vertex_count);
        rank.clear();
        rank.resize(vertex_count, 0);
        open_edges.clear();
        for (index, edge) in lattice.edges.iter().enumerate() {
            let same_spin = spins[edge.source] == spins[edge.target];
            // Short-circuit: opposite-spin edges consume no variate; exact
            // endpoints (0 through beta = 0, 1 only through underflow)
            // consume none either — the F2 endpoint convention.
            let open = same_spin
                && match bond_probabilities[index] {
                    0.0 => false,
                    1.0 => true,
                    probability => rng.random::<f64>() < probability,
                };
            if open {
                open_edges.push(EdgeId::from_valid_index(index));
                union(parent, rank, edge.source, edge.target);
            }
        }
        // The draw is the configuration: ascending dense order is exactly the
        // sorted, deduplicated order BondSample stores.
        self.sample = BondSample::from_open_edges(&self.lattice, self.open_edges.iter().copied())
            .expect("open edges carry this lattice's dense edge IDs");
    }

    /// The conditional spin pass given the FK clusters (contract step 2).
    fn reassign_cluster_spins<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        if self.frozen {
            return;
        }
        let Self {
            lattice,
            spins,
            parent,
            cluster_spin,
            ..
        } = self;
        let vertex_count = lattice.vertex_count();
        cluster_spin.clear();
        cluster_spin.resize(vertex_count, 0);
        for (site, spin) in spins.iter_mut().enumerate() {
            let root = find(parent, site);
            if cluster_spin[root] == 0 {
                cluster_spin[root] = if rng.random::<f64>() < 0.5 { 1 } else { -1 };
            }
            *spin = cluster_spin[root];
        }
    }
}

impl MonteCarlo for EdwardsSokalMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.draw_fk_bonds(&mut ctx.rng);
        self.reassign_cluster_spins(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        // Destructure for the disjoint immutable/mutable field borrows.
        let Self {
            lattice,
            vertices,
            sample,
            workspace,
            plan,
            ..
        } = self;
        // The measured configuration is the exposed sample: the same
        // BondSample returned by `last_bond_sample()` after this sweep.
        let result = if plan.requests_labels() {
            analyze_with_labels(
                lattice,
                vertices,
                sample,
                plan.boundary_queries(),
                workspace,
            )
        } else {
            analyze(
                lattice,
                vertices,
                sample,
                plan.boundary_queries(),
                workspace,
            )
        }
        .expect("constructor validated FK-draw and plan domains");
        let view = SampleView::new(
            lattice.vertex_count(),
            result.summary(),
            result.query_outcomes(),
            None,
            result.labels(),
        );
        plan.record(&view, ctx);
    }

    fn name(&self) -> &'static str {
        "EdwardsSokal"
    }
}

/// Union-find root with full path compression (the SWCore convention).
fn find(parent: &mut [usize], site: usize) -> usize {
    let mut root = site;
    while parent[root] != root {
        root = parent[root];
    }
    let mut current = site;
    while parent[current] != root {
        let next = parent[current];
        parent[current] = root;
        current = next;
    }
    root
}

/// Union by rank; a self-loop or intra-cluster edge is a no-op.
fn union(parent: &mut [usize], rank: &mut [u8], left: usize, right: usize) {
    let left_root = find(parent, left);
    let right_root = find(parent, right);
    if left_root == right_root {
        return;
    }
    match rank[left_root].cmp(&rank[right_root]) {
        Ordering::Less => parent[left_root] = right_root,
        Ordering::Greater => parent[right_root] = left_root,
        Ordering::Equal => {
            parent[right_root] = left_root;
            rank[left_root] += 1;
        }
    }
}

/// Construction or validation failure of an Edwards-Sokal adapter.
#[derive(Debug)]
pub enum EdwardsSokalError {
    /// The lattice failed `CsrLattice::validate`.
    Topology(String),
    /// An edge coupling (`Bond.weight`) is non-finite or not strictly
    /// positive: the FK law `1 - exp(-2 beta J)` is defined for ferromagnetic
    /// `J > 0` only.
    Coupling { edge: usize, weight: f64 },
    /// The inverse temperature is non-finite or negative.
    InverseTemperature { beta: f64 },
    /// A frozen spin vector has the wrong length.
    SpinCount { expected: usize, actual: usize },
    /// A frozen spin value is not `+1` or `-1`.
    SpinValue { site: usize, value: i8 },
    /// The plan requests winding observables, which have no FK analysis yet.
    WindingNotSupported,
    /// The plan itself is invalid for this topology.
    Observable(ObservablePlanError),
    /// A boundary query vertex domain differs from the topology.
    QueryDomain {
        query: usize,
        expected: usize,
        actual: usize,
    },
    /// The component workspace could not be reserved.
    Capacity,
}

impl fmt::Display for EdwardsSokalError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Topology(reason) => write!(formatter, "invalid topology: {reason}"),
            Self::Coupling { edge, weight } => write!(
                formatter,
                "edge {edge} coupling {weight} is not finite and strictly positive; \
                 the Edwards-Sokal FK law requires ferromagnetic J > 0"
            ),
            Self::InverseTemperature { beta } => write!(
                formatter,
                "inverse temperature {beta} must be finite and non-negative"
            ),
            Self::SpinCount { expected, actual } => write!(
                formatter,
                "frozen spins have {actual} entries; topology has {expected} sites"
            ),
            Self::SpinValue { site, value } => write!(
                formatter,
                "spin at site {site} is {value}; expected +1 or -1"
            ),
            Self::WindingNotSupported => write!(
                formatter,
                "winding observables on the FK draw are not implemented yet; \
                 construct the plan without Wrapping/WindingRank"
            ),
            Self::Observable(source) => write!(formatter, "invalid observable plan: {source}"),
            Self::QueryDomain {
                query,
                expected,
                actual,
            } => write!(
                formatter,
                "boundary query {query} has V={actual}; adapter topology has V={expected}"
            ),
            Self::Capacity => formatter.write_str("could not reserve component workspace"),
        }
    }
}

impl std::error::Error for EdwardsSokalError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_chain, BoundaryQuery, StaticObservable};
    use rand::SeedableRng;

    fn plain_plan() -> ObservablePlan {
        ObservablePlan::new([StaticObservable::ActiveEdgeCount], Vec::new())
    }

    fn wedge() -> CsrLattice {
        CsrLattice::from_edges(
            4,
            vec![
                crate::Bond::new(0, 1, crate::BondType::Generic, 1.0),
                crate::Bond::new(1, 2, crate::BondType::Generic, 0.7),
                crate::Bond::new(2, 3, crate::BondType::Generic, 1.3),
                crate::Bond::new(0, 2, crate::BondType::Generic, 0.5),
            ],
        )
    }

    fn context() -> Context<Xoshiro256PlusPlus> {
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(1), 0, 1)
    }

    #[test]
    fn non_ferromagnetic_and_invalid_inputs_are_typed_rejections() {
        let coupling_lattice = |weight: f64| {
            CsrLattice::from_edges(
                3,
                vec![crate::Bond::new(0, 1, crate::BondType::Generic, weight)],
            )
        };
        assert!(matches!(
            EdwardsSokalMC::new(coupling_lattice(0.0), 0.5, plain_plan())
                .err()
                .expect("zero coupling is rejected"),
            EdwardsSokalError::Coupling {
                edge: 0,
                weight: 0.0
            }
        ));
        assert!(matches!(
            EdwardsSokalMC::new(coupling_lattice(-1.0), 0.5, plain_plan())
                .err()
                .expect("negative coupling is rejected"),
            EdwardsSokalError::Coupling {
                edge: 0,
                weight: -1.0
            }
        ));
        // Non-finite weights cannot even be built into a CsrLattice (the
        // constructor validates them), so the ferromagnetic gate's remaining
        // live surface is J <= 0.

        for beta in [-0.1, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                EdwardsSokalMC::new(wedge(), beta, plain_plan())
                    .err()
                    .expect("invalid beta is rejected"),
                EdwardsSokalError::InverseTemperature { .. }
            ));
        }

        assert!(matches!(
            EdwardsSokalMC::frozen(wedge(), 0.5, vec![1, 1, 1], plain_plan())
                .err()
                .expect("short spin vector is rejected"),
            EdwardsSokalError::SpinCount {
                expected: 4,
                actual: 3
            }
        ));
        assert!(matches!(
            EdwardsSokalMC::frozen(wedge(), 0.5, vec![1, 0, -1, 1], plain_plan())
                .err()
                .expect("zero spin is rejected"),
            EdwardsSokalError::SpinValue { site: 1, value: 0 }
        ));

        let winding = ObservablePlan::new([StaticObservable::Wrapping], Vec::new());
        assert!(matches!(
            EdwardsSokalMC::new(wedge(), 0.5, winding)
                .err()
                .expect("winding plan is rejected"),
            EdwardsSokalError::WindingNotSupported
        ));

        // A boundary query built for a different vertex domain is rejected
        // with the query's position and both counts.
        let small = build_chain(2, false);
        let foreign = BoundaryQuery::new(&small, &[0], &[1]).unwrap();
        let plan = ObservablePlan::new([StaticObservable::LargestSize], vec![foreign]);
        assert!(matches!(
            EdwardsSokalMC::new(wedge(), 0.5, plan)
                .err()
                .expect("foreign query domain is rejected"),
            EdwardsSokalError::QueryDomain {
                query: 0,
                expected: 4,
                actual: 2
            }
        ));
    }

    #[test]
    fn frozen_draw_replays_the_documented_rng_contract_bitwise() {
        let lattice = wedge();
        let beta = 0.9;
        let spins = vec![1, 1, -1, 1];
        let mut adapter =
            EdwardsSokalMC::frozen(lattice.clone(), beta, spins.clone(), plain_plan()).unwrap();
        let mut ctx = context();
        adapter.sweep(&mut ctx);

        // Independent replay: one f64 per same-spin edge in dense order,
        // nothing for opposite-spin edges (all p_e strictly inside (0, 1)).
        let mut reference = Xoshiro256PlusPlus::seed_from_u64(1);
        let expected: Vec<bool> = lattice
            .edges
            .iter()
            .map(|edge| {
                let same_spin = spins[edge.source] == spins[edge.target];
                let probability = 1.0 - (-2.0 * beta * edge.weight).exp();
                same_spin && reference.random::<f64>() < probability
            })
            .collect();
        assert_eq!(adapter.last_bond_sample().to_mask(), expected);
        assert_eq!(adapter.spins(), spins.as_slice(), "frozen spins never move");
        adapter.sweep(&mut ctx);
        assert_eq!(adapter.spins(), spins.as_slice());
    }

    #[test]
    fn unfrozen_sweep_replays_bonds_and_cluster_spins_bitwise() {
        let lattice = build_chain(5, false);
        let beta = 0.8;
        let mut adapter = EdwardsSokalMC::new(lattice.clone(), beta, plain_plan()).unwrap();
        let mut ctx = context();
        adapter.sweep(&mut ctx);

        let mut reference = Xoshiro256PlusPlus::seed_from_u64(1);
        let mut spins: Vec<i8> = vec![1; lattice.vertex_count()];
        for _ in 0..3 {
            // Bond pass.
            let mut mask = vec![false; lattice.edge_count()];
            let mut parent: Vec<usize> = (0..lattice.vertex_count()).collect();
            let mut rank = vec![0u8; lattice.vertex_count()];
            for (index, edge) in lattice.edges.iter().enumerate() {
                if spins[edge.source] == spins[edge.target] {
                    let probability = 1.0 - (-2.0 * beta * edge.weight).exp();
                    if reference.random::<f64>() < probability {
                        mask[index] = true;
                        union(&mut parent, &mut rank, edge.source, edge.target);
                    }
                }
            }
            // Spin pass: one f64 per cluster at its lowest-index site.
            let mut cluster_spin = vec![0i8; lattice.vertex_count()];
            for (site, spin) in spins.iter_mut().enumerate() {
                let root = find(&mut parent, site);
                if cluster_spin[root] == 0 {
                    cluster_spin[root] = if reference.random::<f64>() < 0.5 {
                        1
                    } else {
                        -1
                    };
                }
                *spin = cluster_spin[root];
            }
            assert_eq!(adapter.last_bond_sample().to_mask(), mask);
            assert_eq!(adapter.spins(), spins.as_slice());
            if adapter.last_bond_sample().open_count() == 0 {
                // Nothing more to learn from further all-closed steps.
                break;
            }
            adapter.sweep(&mut ctx);
        }
    }

    #[test]
    fn beta_endpoints_are_deterministic_and_consume_nothing() {
        // beta = 0: p_e = 0 exactly, every bond closed, no variates.
        let lattice = wedge();
        let mut cold = EdwardsSokalMC::new(lattice.clone(), 0.0, plain_plan()).unwrap();
        let mut ctx = context();
        cold.sweep(&mut ctx);
        assert_eq!(cold.last_bond_sample().open_count(), 0);

        // beta large enough that exp underflows: p_e = 1 exactly, every
        // same-spin edge open. All-up spins (the constructor default) make
        // every edge same-spin.
        let mut hot = EdwardsSokalMC::new(lattice.clone(), 1000.0, plain_plan()).unwrap();
        hot.sweep(&mut ctx);
        assert_eq!(hot.last_bond_sample().open_count(), lattice.edge_count());
    }

    #[test]
    fn measure_records_planned_observables_on_the_exposed_sample() {
        let lattice = build_chain(4, false);
        let query = BoundaryQuery::new(&lattice, &[0], &[3]).unwrap();
        let plan = ObservablePlan::new(
            [
                StaticObservable::ActiveEdgeCount,
                StaticObservable::ComponentCount,
                StaticObservable::LargestSize,
            ],
            vec![query],
        );
        let mut adapter = EdwardsSokalMC::new(lattice.clone(), 0.0, plan).unwrap();
        let mut ctx = context();
        for _ in 0..4 {
            adapter.sweep(&mut ctx);
            adapter.measure(&mut ctx);
        }
        let values = ctx.finalize_measurements();
        assert_eq!(values["ActiveEdgeCount"].mean, 0.0);
        assert_eq!(values["ComponentCount"].mean, 4.0);
        assert_eq!(values["LargestSize"].mean, 1.0);
        assert_eq!(values["BoundaryCrossing0"].mean, 0.0);
        assert_eq!(values["BoundaryCrossing0"].n_bins, 4);
    }

    #[test]
    fn accessors_expose_the_validated_model_state() {
        let lattice = wedge();
        let beta = 0.9;
        let adapter = EdwardsSokalMC::new(lattice.clone(), beta, plain_plan()).unwrap();
        assert!(!adapter.is_frozen());
        assert_eq!(adapter.beta(), beta);
        assert_eq!(adapter.lattice().edge_count(), lattice.edge_count());
        assert_eq!(adapter.spins(), vec![1; 4].as_slice());
        for (edge, &probability) in adapter.bond_probabilities().iter().enumerate() {
            let coupling = lattice.edges[edge].weight;
            assert_close(probability, 1.0 - (-2.0 * beta * coupling).exp());
        }
        // Empty pre-sweep sample still pins the correct edge domain.
        assert_eq!(
            adapter.last_bond_sample().edge_count(),
            lattice.edge_count()
        );
        assert_eq!(adapter.last_bond_sample().open_count(), 0);

        let frozen =
            EdwardsSokalMC::frozen(lattice, beta, vec![-1, 1, -1, 1], plain_plan()).unwrap();
        assert!(frozen.is_frozen());
    }

    fn assert_close(left: f64, right: f64) {
        assert!((left - right).abs() <= 1e-15, "{left:.17e} != {right:.17e}");
    }
}
