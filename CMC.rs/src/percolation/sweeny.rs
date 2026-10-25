//! Single-bond Sweeny sampler for the random-cluster model at arbitrary
//! real `q > 0` (roadmap v2 §6 Q1).
//!
//! The random-cluster (Fortuin-Kasteleyn) measure on a graph `G = (V, E)`
//! with edge parameter `p` and cluster weight `q` is
//!
//! ```text
//! phi_{p,q}(A) = p^{|A|} (1-p)^{E-|A|} q^{k(A)} / Z(p, q),   A subset of E,
//! v            = p / (1 - p),
//! w_{p,q}(A)   = v^{|A|} q^{k(A)}          (the same measure up to Z),
//! ```
//!
//! where `k(A)` is the number of connected components of `(V, A)` with every
//! vertex active (the platform's bond convention). This module samples
//! `phi_{p,q}` directly with Sweeny's single-bond heat-bath/Metropolis
//! dynamics: one physical edge is proposed at a time, and the acceptance
//! probability depends on the change in the number of open bonds `Δb` and
//! the change in the number of components `Δn`.
//!
//! # Acceptance formulas (Geyer 2015, Eqs. 4 and 5)
//!
//! With `v = p/(1-p)`, an insertion has `Δb = +1` and `Δn ∈ {0, -1}`
//! (internal edge versus edge amalgamating two clusters); a deletion has
//! `Δb = -1` and `Δn ∈ {0, +1}` (internal edge versus bridge). The heat-bath
//! acceptance used by Sweeny is
//!
//! ```text
//! p_acc^HB(Δb, Δn) = q^{Δn} v^{Δb} / (1 + q^{Δn} v^{Δb})      (Geyer Eq. 4)
//! ```
//!
//! and the Metropolis-Hastings variant is
//!
//! ```text
//! p_acc^MH(Δb, Δn) = min(1, q^{Δn} v^{Δb}).                    (Geyer Eq. 5)
//! ```
//!
//! Geyer's archived text (lines 259-265 of
//! `docs/ref/md/Geyer_2015_sweeny_random_cluster.md`) states the ratio bound
//!
//! ```text
//! 1/2 <= p_acc^HB / p_acc^MH = max(1, q^{Δn} v^{Δb}) / (1 + q^{Δn} v^{Δb}) < 1,
//! ```
//!
//! with the minimum `1/2` reached at `q^{Δn} v^{Δb} = 1` (the percolation
//! limit `q -> 1` at criticality `v = sqrt(q)`). Both acceptances are exposed
//! through [`SweenyRates`]; the sampler chooses one through
//! [`SweenyAcceptance`].
//!
//! # Dynamic connectivity: spanning forest with replacement search
//!
//! Sweeny's bottleneck is the `Δn` query. This module maintains a **spanning
//! forest** of the open bonds with parent/child pointers and a per-edge
//! tree/non-tree flag:
//!
//! * insertion of an edge joining two trees re-roots the **smaller** tree at
//!   its endpoint and links it through the actual edge (so every stored tree
//!   edge is a real open edge), and an edge inside one tree is a non-tree
//!   edge (`Δn = 0`); the link is `O(min(|C_u|, |C_v|) + depth)`;
//! * deletion of a non-tree edge is `O(1)` (`Δn = 0`);
//! * deletion of a tree edge runs an **interleaved breadth-first search**
//!   from both endpoints of the cut (one vertex expanded from each side per
//!   round, so the smaller side is exhausted first), then scans the incident
//!   edges of that smaller side for an open non-tree **replacement edge**.
//!   If one exists the sides are reconnected through it (`Δn = 0`);
//!   otherwise the component splits (`Δn = +1`). The cost is proportional to
//!   the volume of the smaller side, never a full component recomputation and
//!   never `O(E)` per proposal. The worst case for a bridge is `O(V)`.
//!
//! This is the spanning-forest/replacement-search scheme Geyer discusses
//! around his Table I (`docs/ref/md/Geyer_2015_sweeny_random_cluster.md`,
//! lines 548-575 and 585-605); it is not the asymptotically optimal
//! Holm-de Lichtenberg-Thorup structure, and the docs say so. Storing real
//! edges (instead of attaching a component root directly through a "virtual"
//! link) is what makes deletion exact: a later cut on the path cannot
//! invalidate a tree link. The price is the smaller-component re-root on
//! merging insertions, `O(min(|C_u|, |C_v|))`; `find` intentionally performs
//! no path compression because the parent/child adjacency must stay
//! consistent for cut and reroot, so a pathological sequence of re-roots can
//! produce deep trees. This is documented rather than hidden. The
//! connectivity query count is reported by [`SweenyMC::delta_n_queries`].
//!
//! # Unconditional decisions (Geyer's optimization)
//!
//! For a fixed proposal direction the two possible `Δn` values give two
//! acceptance probabilities `p_low <= p_high`. One uniform `u` is drawn per
//! proposal; `u < p_low` accepts and `u >= p_high` rejects **without** a
//! connectivity query, and only the interval `[p_low, p_high)` needs the
//! query. This is exactly Geyer's unconditional-accept/reject shortcut
//! (`docs/ref/md/Geyer_2015_sweeny_random_cluster.md`, lines 277-302); the
//! decision function is unchanged, so detailed balance is untouched.
//!
//! # RNG contract (documented and replayed by the tests)
//!
//! One `sweep()` proposes every physical edge exactly once, in ascending
//! dense [`EdgeId`] order. Each proposal consumes exactly one
//! `rng.random::<f64>()` uniform in `[0, 1)`; the move is accepted iff
//! `u < p_acc` for the chosen acceptance rule (heat-bath or Metropolis) and
//! the actual `Δn`. No other RNG call is made. Consequently:
//!
//! * **`q = 1`, heat-bath, all-closed start, one sweep** is bitwise the F2
//!   `BondBernoulli(p)` draw of R2 at the same seed: every edge is closed
//!   when proposed, the insertion acceptance is exactly `p`, and the
//!   comparison is the same `u < p`.
//! * The unconditional shortcut consumes the same single uniform; it only
//!   skips the `Δn` query, not a draw.
//!
//! # Allocation and cost
//!
//! All buffers (open mask, forest parent/child adjacency, cut scratch,
//! analyzer workspace) are reserved at construction. A steady-state
//! `sweep()` performs no heap allocation; only [`SweenyMC::realization`]
//! allocates on demand. `ComponentWorkspace` is prepared once. The only
//! expensive operation is the tree-edge replacement search described above;
//! [`SweenyMC::delta_n_queries`] counts how often it (or an insertion
//! connectivity query) was actually needed.
//!
//! # References
//!
//! M. Geyer, *Efficient simulation of the random-cluster model*, J. Stat.
//! Mech. (2015) P08010, arXiv:1307.6647 — Eqs. 2-5, the ratio bound, and the
//! connectivity-cost discussion (archived in `docs/ref/`). Sweeny's original
//! algorithm (Phys. Rev. B **27**, 4445 (1983)) predates arXiv; Geyer's
//! archived comparison is the operational anchor here. The exact FK weight
//! enumeration used by the gates is the Fortuin-Kasteleyn identity
//! (Physica **57**, 536 (1972), pre-arXiv).

use core::fmt;

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{
    analyze_bonds, analyze_bonds_with_labels, BondSample, ComponentWorkspace, ObservablePlan,
    ObservablePlanError, SampleView,
};
use crate::{CsrLattice, GraphView, UndirectedGraphView};

/// Sentinel for "no edge / no parent / no colour".
const NONE: usize = usize::MAX;

/// Acceptance rule exposed by the Sweeny sampler.
///
/// `HeatBath` uses Geyer Eq. 4, `Metropolis` Geyer Eq. 5. Both target the
/// same random-cluster measure; Geyer's experiments used Metropolis, his
/// historical algorithm used heat-bath.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SweenyAcceptance {
    HeatBath,
    Metropolis,
}

impl SweenyAcceptance {
    pub const fn name(self) -> &'static str {
        match self {
            Self::HeatBath => "heat-bath",
            Self::Metropolis => "Metropolis",
        }
    }
}

/// Validated `(q, p)` pair together with the odds `v = p/(1-p)`.
///
/// The acceptance formulas are pure functions of `q`, `v`, `Δb`, and `Δn`;
/// validating once here keeps every hot-path call panic-free and keeps the
/// formulas independently testable.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SweenyRates {
    q: f64,
    probability: f64,
    odds: f64,
    ln_q: f64,
    ln_odds: f64,
}

impl SweenyRates {
    /// Validate `q > 0` finite and `p` strictly inside `(0, 1)` finite.
    ///
    /// Endpoints are rejected: Sweeny's single-bond heat-bath needs a finite
    /// odds `v = p/(1-p)`, and the random-cluster measure at `p = 0` or
    /// `p = 1` is a deterministic (and separately trivial) configuration.
    pub fn new(q: f64, probability: f64) -> Result<Self, SweenyError> {
        if !q.is_finite() || q <= 0.0 {
            return Err(SweenyError::InvalidQ { q });
        }
        if !probability.is_finite() || probability <= 0.0 || probability >= 1.0 {
            return Err(SweenyError::InvalidProbability { probability });
        }
        let odds = probability / (1.0 - probability);
        if !odds.is_finite() || odds <= 0.0 {
            return Err(SweenyError::NonFiniteOdds { probability, odds });
        }
        Ok(Self {
            q,
            probability,
            odds,
            ln_q: q.ln(),
            ln_odds: odds.ln(),
        })
    }

    /// The cluster weight `q`.
    #[inline]
    pub const fn q(self) -> f64 {
        self.q
    }

    /// The bond probability `p`.
    #[inline]
    pub const fn probability(self) -> f64 {
        self.probability
    }

    /// The odds `v = p/(1-p)`.
    #[inline]
    pub const fn odds(self) -> f64 {
        self.odds
    }

    /// `log(q^{Δn} v^{Δb})`, the numerically stable core of both formulas.
    #[inline]
    pub fn log_ratio(self, delta_n: i32, delta_b: i32) -> f64 {
        f64::from(delta_n) * self.ln_q + f64::from(delta_b) * self.ln_odds
    }

    /// The weight ratio `q^{Δn} v^{Δb}`, saturated to `0`/`+inf` outside the
    /// f64 exponent range.
    pub fn ratio(self, delta_n: i32, delta_b: i32) -> f64 {
        let log_ratio = self.log_ratio(delta_n, delta_b);
        if log_ratio > 700.0 {
            f64::INFINITY
        } else if log_ratio < -745.0 {
            0.0
        } else {
            log_ratio.exp()
        }
    }

    /// Heat-bath acceptance, Geyer Eq. 4.
    ///
    /// Computed from the log ratio as `1/(1+e^{-l})` for `l >= 0` and
    /// `e^l/(1+e^l)` otherwise, so it never overflows even when the weight
    /// ratio itself would.
    pub fn heat_bath(self, delta_n: i32, delta_b: i32) -> f64 {
        let log_ratio = self.log_ratio(delta_n, delta_b);
        if log_ratio >= 0.0 {
            1.0 / (1.0 + (-log_ratio).exp())
        } else {
            let ratio = log_ratio.exp();
            ratio / (1.0 + ratio)
        }
    }

    /// Metropolis-Hastings acceptance, Geyer Eq. 5: `min(1, ratio)`.
    pub fn metropolis(self, delta_n: i32, delta_b: i32) -> f64 {
        let log_ratio = self.log_ratio(delta_n, delta_b);
        if log_ratio >= 0.0 {
            1.0
        } else {
            log_ratio.exp()
        }
    }

    /// The exact ratio `p_acc^HB / p_acc^MH` as the closed form
    /// `1/(1 + exp(-|log ratio|))`, in `(1/2, 1]`.
    ///
    /// For `|log ratio| > ~37` the f64 result rounds to exactly `1.0`; the
    /// mathematical inequality stays strict. The acceptance-bound gate
    /// exercises the non-saturated range and separately asserts the closed
    /// form.
    pub fn heat_bath_metropolis_ratio(self, delta_n: i32, delta_b: i32) -> f64 {
        let magnitude = self.log_ratio(delta_n, delta_b).abs();
        1.0 / (1.0 + (-magnitude).exp())
    }
}

/// One proposal record returned by [`SweenyMC::propose`].
///
/// `delta_n` is the actual component-count change of an accepted move, or the
/// queried change of a rejected ambiguous move; it is `None` only when the
/// move was accepted or rejected unconditionally and no query was needed.
/// `connectivity_queried` is true when the decision or the accepted state
/// update consulted the forest (an insertion link query or a deletion
/// replacement search).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SweenyProposal {
    pub edge: usize,
    pub insertion: bool,
    pub accepted: bool,
    pub delta_n: Option<i32>,
    pub draw: f64,
    pub connectivity_queried: bool,
}

/// Sweeny single-bond random-cluster sampler.
///
/// The chain starts all-closed. One [`MonteCarlo::sweep`] is one ascending
/// pass over every physical edge; the current bond configuration is exposed
/// as a dense mask ([`Self::mask`]) and as an R1 [`BondSample`]
/// ([`Self::realization`]), and [`MonteCarlo::measure`] analyzes that same
/// mask through the shared F3 analyzer and records the [`ObservablePlan`].
///
/// This is a Markov chain: thermalize before measuring. Winding observables
/// are a typed rejection (no embedding), like the other FK adapters.
pub struct SweenyMC {
    lattice: CsrLattice,
    rates: SweenyRates,
    acceptance: SweenyAcceptance,
    open: Vec<bool>,
    forest: SpanningForest,
    open_count: usize,
    component_count: usize,
    proposal_count: u64,
    delta_n_queries: u64,
    workspace: ComponentWorkspace,
    plan: ObservablePlan,
}

impl SweenyMC {
    /// Construct an all-closed Sweeny chain.
    ///
    /// Rejections: `q <= 0`/non-finite, `p` outside the open unit interval,
    /// non-finite odds, invalid topology, winding plans, invalid plans for
    /// this vertex count, and boundary queries with a foreign vertex domain.
    pub fn new(
        lattice: CsrLattice,
        q: f64,
        bond_probability: f64,
        acceptance: SweenyAcceptance,
        plan: ObservablePlan,
    ) -> Result<Self, SweenyError> {
        lattice
            .validate()
            .map_err(|error| SweenyError::Topology(error.to_string()))?;
        let rates = SweenyRates::new(q, bond_probability)?;
        if plan.requests_winding() {
            return Err(SweenyError::WindingUnsupported);
        }
        plan.validate(lattice.vertex_count())
            .map_err(SweenyError::Observable)?;
        for (query, boundary) in plan.boundary_queries().iter().enumerate() {
            if boundary.vertex_count() != lattice.vertex_count() {
                return Err(SweenyError::QueryDomain {
                    query,
                    expected: lattice.vertex_count(),
                    actual: boundary.vertex_count(),
                });
            }
        }
        let mut workspace = ComponentWorkspace::new();
        workspace
            .prepare(lattice.vertex_count(), plan.boundary_queries().len(), false)
            .map_err(|_| SweenyError::Capacity)?;
        let forest = SpanningForest::new(&lattice);
        let vertex_count = lattice.vertex_count();
        let edge_count = lattice.edge_count();
        Ok(Self {
            lattice,
            rates,
            acceptance,
            open: vec![false; edge_count],
            forest,
            open_count: 0,
            component_count: vertex_count,
            proposal_count: 0,
            delta_n_queries: 0,
            workspace,
            plan,
        })
    }

    /// The owned topology.
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// The validated `(q, p, v)` rates.
    pub const fn rates(&self) -> SweenyRates {
        self.rates
    }

    /// The selected acceptance rule.
    pub const fn acceptance(&self) -> SweenyAcceptance {
        self.acceptance
    }

    /// The observable plan recorded by `measure()`.
    pub const fn observable_plan(&self) -> &ObservablePlan {
        &self.plan
    }

    /// The dense bond mask of the current configuration.
    #[inline]
    pub fn mask(&self) -> &[bool] {
        &self.open
    }

    /// Number of open bonds `|A|`.
    #[inline]
    pub const fn open_edge_count(&self) -> usize {
        self.open_count
    }

    /// Number of connected components `k(A)` of the current configuration.
    #[inline]
    pub const fn component_count(&self) -> usize {
        self.component_count
    }

    /// Proposals executed since construction (or the last [`Self::reset`]).
    #[inline]
    pub const fn proposal_count(&self) -> u64 {
        self.proposal_count
    }

    /// `Δn` queries executed since construction (or the last reset): an
    /// insertion connectivity check or a deletion replacement search. With
    /// the unconditional shortcut this can be below the proposal count.
    #[inline]
    pub const fn delta_n_queries(&self) -> u64 {
        self.delta_n_queries
    }

    /// The current configuration as R1's owned open-edge record.
    pub fn realization(&self) -> BondSample {
        BondSample::from_mask(&self.lattice, &self.open)
            .expect("mask length pinned to the topology at construction")
    }

    /// Restore the all-closed start state and clear the proposal counters.
    ///
    /// The forest and analyzer workspace keep their reservations.
    pub fn reset(&mut self) {
        self.open.fill(false);
        self.open_count = 0;
        self.component_count = self.lattice.vertex_count();
        self.forest.reset(&self.lattice);
        self.proposal_count = 0;
        self.delta_n_queries = 0;
    }

    /// One single-bond proposal on `edge_index`, consuming exactly one
    /// uniform.
    ///
    /// This is the low-level step used by [`MonteCarlo::sweep`] and by the
    /// update-correctness fuzz; it is not normally called by users. The
    /// returned record documents whether the move was accepted, its actual
    /// `Δn` when known, and whether connectivity was consulted.
    pub fn propose<R: Rng + ?Sized>(
        &mut self,
        edge_index: usize,
        rng: &mut R,
    ) -> Result<SweenyProposal, SweenyError> {
        if edge_index >= self.lattice.edge_count() {
            return Err(SweenyError::EdgeIndex {
                index: edge_index,
                edge_count: self.lattice.edge_count(),
            });
        }
        let bond = self.lattice.edges[edge_index];
        let (source, target) = (bond.source, bond.target);
        let insertion = !self.open[edge_index];
        let delta_b = if insertion { 1 } else { -1 };

        // The two possible `Δn` values: insertion { -1, 0 }, deletion
        // { 0, +1 }. Their acceptance order depends on q: for q < 1 the
        // cluster-weight factor reverses it, so sort the two probabilities
        // explicitly before using them as the unconditional bounds.
        let (delta_n_first, delta_n_second) = if insertion { (-1, 0) } else { (0, 1) };
        let p_first = self.acceptance_value(delta_n_first, delta_b);
        let p_second = self.acceptance_value(delta_n_second, delta_b);
        let p_low = p_first.min(p_second);
        let p_high = p_first.max(p_second);

        let draw = rng.random::<f64>();
        self.proposal_count += 1;

        let mut delta_n = None;
        let mut connectivity_queried = false;
        // A cut plan computed for the acceptance decision is reused by the
        // accepted state update, so a single replacement search suffices.
        let mut pending_cut: Option<CutPlan> = None;
        let accepted = if draw < p_low {
            true
        } else if draw >= p_high {
            false
        } else {
            // The ambiguous interval genuinely needs Δn.
            connectivity_queried = true;
            self.delta_n_queries += 1;
            if insertion {
                let connected = self.forest.connected(source, target);
                let actual = if connected { 0 } else { -1 };
                delta_n = Some(actual);
                draw < self.acceptance_value(actual, delta_b)
            } else if self.forest.is_tree(edge_index) {
                let plan = self.forest.plan_cut(edge_index, &self.lattice, &self.open);
                let split = plan.replacement.is_none();
                let actual = if split { 1 } else { 0 };
                delta_n = Some(actual);
                let accepted = draw < self.acceptance_value(actual, delta_b);
                if accepted {
                    pending_cut = Some(plan);
                }
                accepted
            } else {
                delta_n = Some(0);
                draw < self.acceptance_value(0, delta_b)
            }
        };

        if accepted {
            if insertion {
                self.open[edge_index] = true;
                self.open_count += 1;
                // If the decision did not query connectivity, the state
                // update still must know whether to link or to store a
                // non-tree edge.
                let actual = match delta_n {
                    Some(value) => value,
                    None => {
                        connectivity_queried = true;
                        self.delta_n_queries += 1;
                        if self.forest.connected(source, target) {
                            0
                        } else {
                            -1
                        }
                    }
                };
                if actual == 0 {
                    self.forest.add_non_tree(edge_index);
                } else {
                    debug_assert_eq!(actual, -1);
                    self.forest.link(source, target, edge_index);
                    self.component_count -= 1;
                }
                delta_n = Some(actual);
            } else {
                self.open[edge_index] = false;
                self.open_count -= 1;
                let split = if self.forest.is_tree(edge_index) {
                    let plan = match pending_cut.take() {
                        Some(plan) => plan,
                        None => {
                            connectivity_queried = true;
                            self.delta_n_queries += 1;
                            self.forest.plan_cut(edge_index, &self.lattice, &self.open)
                        }
                    };
                    let split = plan.replacement.is_none();
                    self.forest.apply_cut(edge_index, plan);
                    split
                } else {
                    self.forest.remove_non_tree(edge_index);
                    false
                };
                if split {
                    self.component_count += 1;
                    delta_n = Some(1);
                } else {
                    delta_n = Some(0);
                }
            }
        }

        Ok(SweenyProposal {
            edge: edge_index,
            insertion,
            accepted,
            delta_n,
            draw,
            connectivity_queried,
        })
    }

    #[inline]
    fn acceptance_value(&self, delta_n: i32, delta_b: i32) -> f64 {
        match self.acceptance {
            SweenyAcceptance::HeatBath => self.rates.heat_bath(delta_n, delta_b),
            SweenyAcceptance::Metropolis => self.rates.metropolis(delta_n, delta_b),
        }
    }
}

impl MonteCarlo for SweenyMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        for edge in 0..self.lattice.edge_count() {
            self.propose(edge, &mut ctx.rng)
                .expect("dense edge index inside the topology");
        }
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let result = if self.plan.requests_labels() {
            analyze_bonds_with_labels(
                &self.lattice,
                &self.open,
                self.plan.boundary_queries(),
                &mut self.workspace,
            )
        } else {
            analyze_bonds(
                &self.lattice,
                &self.open,
                self.plan.boundary_queries(),
                &mut self.workspace,
            )
        }
        .expect("constructor prepared matching analysis domains");
        debug_assert_eq!(
            result.summary().component_count,
            self.component_count,
            "incremental component count must match the analyzer"
        );
        debug_assert_eq!(
            result.summary().active_edge_count,
            self.open_count,
            "incremental open count must match the analyzer"
        );
        let sample = SampleView::new(
            self.lattice.vertex_count(),
            result.summary(),
            result.query_outcomes(),
            None,
            result.labels(),
        );
        self.plan.record(&sample, ctx);
    }

    fn name(&self) -> &'static str {
        "Sweeny"
    }
}

/// A cut plan computed without mutating the forest.
///
/// `plan_cut` runs the interleaved search and the replacement scan; the
/// acceptance decision only needs `replacement.is_none()`, and an accepted
/// move applies the plan without repeating the expensive search.
struct CutPlan {
    child: usize,
    parent: usize,
    total: usize,
    side_is_a: bool,
    side_len: usize,
    replacement: Option<(usize, usize, usize)>,
}

/// Spanning forest over the open bonds with tree/non-tree classification.
///
/// Invariants: roots have `parent[root] == root` and `parent_edge[root] ==
/// NONE`; every non-root `v` appears once in `children[parent[v]]`; a tree
/// edge `e` has `tree_edge[e] == true` and `child_of_edge[e] == v` for the
/// unique vertex whose parent edge it is; `size[root]` is the component size.
/// Non-tree edges have `child_of_edge[e] == NONE`.
struct SpanningForest {
    parent: Vec<usize>,
    parent_edge: Vec<usize>,
    children: Vec<Vec<usize>>,
    child_of_edge: Vec<usize>,
    tree_edge: Vec<bool>,
    size: Vec<usize>,
    mark: Vec<u32>,
    stamp: u32,
    queue_a: Vec<usize>,
    queue_b: Vec<usize>,
    side_a: Vec<usize>,
    side_b: Vec<usize>,
    path_vertices: Vec<usize>,
    path_edges: Vec<usize>,
}

impl SpanningForest {
    fn new(lattice: &CsrLattice) -> Self {
        let vertex_count = lattice.vertex_count();
        let edge_count = lattice.edge_count();
        let mut children = Vec::with_capacity(vertex_count);
        for vertex in 0..vertex_count {
            // A spanning tree can use at most one incident edge per
            // neighbour, so the graph degree bounds the tree children.
            children.push(Vec::with_capacity(lattice.degree(vertex)));
        }
        Self {
            parent: (0..vertex_count).collect(),
            parent_edge: vec![NONE; vertex_count],
            children,
            child_of_edge: vec![NONE; edge_count],
            tree_edge: vec![false; edge_count],
            size: vec![1; vertex_count],
            mark: vec![0; vertex_count],
            stamp: 0,
            queue_a: Vec::with_capacity(vertex_count),
            queue_b: Vec::with_capacity(vertex_count),
            side_a: Vec::with_capacity(vertex_count),
            side_b: Vec::with_capacity(vertex_count),
            path_vertices: Vec::with_capacity(vertex_count),
            path_edges: Vec::with_capacity(vertex_count),
        }
    }

    fn reset(&mut self, lattice: &CsrLattice) {
        let vertex_count = lattice.vertex_count();
        self.parent.clear();
        self.parent.extend(0..vertex_count);
        self.parent_edge.fill(NONE);
        for children in &mut self.children {
            children.clear();
        }
        self.child_of_edge.fill(NONE);
        self.tree_edge.fill(false);
        self.size.fill(1);
        self.mark.fill(0);
        self.stamp = 0;
        self.queue_a.clear();
        self.queue_b.clear();
        self.side_a.clear();
        self.side_b.clear();
        self.path_vertices.clear();
        self.path_edges.clear();
    }

    #[inline]
    fn find(&self, mut vertex: usize) -> usize {
        while self.parent[vertex] != vertex {
            vertex = self.parent[vertex];
        }
        vertex
    }

    #[inline]
    fn connected(&self, left: usize, right: usize) -> bool {
        self.find(left) == self.find(right)
    }

    #[inline]
    fn is_tree(&self, edge: usize) -> bool {
        self.tree_edge[edge]
    }

    fn add_non_tree(&mut self, edge: usize) {
        debug_assert!(!self.tree_edge[edge]);
        self.child_of_edge[edge] = NONE;
    }

    fn remove_non_tree(&mut self, edge: usize) {
        debug_assert!(!self.tree_edge[edge]);
        self.child_of_edge[edge] = NONE;
    }

    /// Link two different components with the **actual** edge `source-target`.
    ///
    /// The smaller component is re-rooted at its endpoint of the edge and
    /// attached below the other endpoint, so every tree edge stored here is a
    /// real open edge. Storing virtual links (attaching the component root
    /// directly) would be cheaper, but a later cut elsewhere on the path could
    /// invalidate such a link; real-edge trees keep the deletion semantics
    /// exact. The price is a re-root of the smaller component on every
    /// component-merging insertion, documented in the module cost section.
    fn link(&mut self, source: usize, target: usize, edge: usize) {
        let source_root = self.find(source);
        let target_root = self.find(target);
        debug_assert_ne!(source_root, target_root);
        debug_assert!(!self.tree_edge[edge]);
        let source_size = self.size[source_root];
        let target_size = self.size[target_root];
        if source_size >= target_size {
            self.reroot(target);
            self.size[target] = target_size;
            self.parent[target] = source;
            self.parent_edge[target] = edge;
            self.children[source].push(target);
            self.child_of_edge[edge] = target;
            self.size[source_root] = source_size + target_size;
        } else {
            self.reroot(source);
            self.size[source] = source_size;
            self.parent[source] = target;
            self.parent_edge[source] = edge;
            self.children[target].push(source);
            self.child_of_edge[edge] = source;
            self.size[target_root] = source_size + target_size;
        }
        self.tree_edge[edge] = true;
    }

    /// Compute a cut plan (side split + replacement edge) without mutation.
    ///
    /// The deleted edge is still open and still a tree edge at this point;
    /// the replacement scan excludes it through the tree flag.
    fn plan_cut(&mut self, edge: usize, lattice: &CsrLattice, open: &[bool]) -> CutPlan {
        let child = self.child_of_edge[edge];
        debug_assert_ne!(child, NONE);
        let parent = self.parent[child];
        let total = self.size[self.find(parent)];
        let (side_is_a, side_len, stamp_side) = self.explore_sides(child, parent);
        let replacement = self.scan_replacement(side_is_a, stamp_side, lattice, open);
        CutPlan {
            child,
            parent,
            total,
            side_is_a,
            side_len,
            replacement,
        }
    }

    /// Apply a previously computed cut plan.
    ///
    /// Returns `true` when the component split (no replacement edge).
    fn apply_cut(&mut self, edge: usize, plan: CutPlan) -> bool {
        let position = self.children[plan.parent]
            .iter()
            .position(|&child| child == plan.child)
            .expect("spanning-forest parent/child pointers must agree");
        self.children[plan.parent].swap_remove(position);
        self.parent[plan.child] = plan.child;
        self.parent_edge[plan.child] = NONE;
        self.child_of_edge[edge] = NONE;
        self.tree_edge[edge] = false;

        if let Some((side_vertex, other_vertex, replacement)) = plan.replacement {
            self.reroot(side_vertex);
            self.size[side_vertex] = plan.side_len;
            let other_root = self.find(other_vertex);
            self.parent[side_vertex] = other_vertex;
            self.parent_edge[side_vertex] = replacement;
            self.children[other_vertex].push(side_vertex);
            self.child_of_edge[replacement] = side_vertex;
            self.tree_edge[replacement] = true;
            // The component is reconnected, so its size is unchanged. The
            // other side's stored size may still be the pre-cut total.
            self.size[other_root] = plan.total;
            false
        } else {
            let (side_root, other_root) = if plan.side_is_a {
                (plan.child, self.find(plan.parent))
            } else {
                (self.find(plan.parent), plan.child)
            };
            self.size[side_root] = plan.side_len;
            self.size[other_root] = plan.total - plan.side_len;
            true
        }
    }

    /// Interleaved BFS of the two sides of the cut `(child, parent)`.
    ///
    /// Side A is the child's subtree (children only); side B is everything
    /// else (children except `child`, plus the parent chain). The side with
    /// fewer vertices is exhausted first because each round expands one
    /// vertex per side, so the returned side is the smaller one (ties return
    /// A). The forest is not mutated; `parent[child]` is still the parent.
    ///
    /// Returns `(side_is_a, side_len, side_stamp)`. The *other* side is only
    /// partially explored when the search stops, so the replacement scan must
    /// identify it by "not in the smaller side" rather than by a visited
    /// stamp; an open edge leaving the smaller side can only enter the other
    /// side of the same component.
    fn explore_sides(&mut self, child: usize, parent: usize) -> (bool, usize, u32) {
        let (stamp_a, stamp_b) = self.next_stamps();
        self.queue_a.clear();
        self.queue_b.clear();
        self.side_a.clear();
        self.side_b.clear();
        self.mark[child] = stamp_a;
        self.queue_a.push(child);
        self.side_a.push(child);
        self.mark[parent] = stamp_b;
        self.queue_b.push(parent);
        self.side_b.push(parent);
        let mut index_a = 0usize;
        let mut index_b = 0usize;
        loop {
            if index_a == self.queue_a.len() {
                return (true, self.side_a.len(), stamp_a);
            }
            if index_b == self.queue_b.len() {
                return (false, self.side_b.len(), stamp_b);
            }
            let vertex_a = self.queue_a[index_a];
            index_a += 1;
            for index in 0..self.children[vertex_a].len() {
                let next = self.children[vertex_a][index];
                if self.mark[next] != stamp_a && self.mark[next] != stamp_b {
                    self.mark[next] = stamp_a;
                    self.queue_a.push(next);
                    self.side_a.push(next);
                }
            }
            let vertex_b = self.queue_b[index_b];
            index_b += 1;
            for index in 0..self.children[vertex_b].len() {
                let next = self.children[vertex_b][index];
                if next == child {
                    continue; // the cut edge
                }
                if self.mark[next] != stamp_a && self.mark[next] != stamp_b {
                    self.mark[next] = stamp_b;
                    self.queue_b.push(next);
                    self.side_b.push(next);
                }
            }
            let up = self.parent[vertex_b];
            if up != vertex_b && self.mark[up] != stamp_a && self.mark[up] != stamp_b {
                self.mark[up] = stamp_b;
                self.queue_b.push(up);
                self.side_b.push(up);
            }
        }
    }

    /// Find an open non-tree edge crossing from the smaller side.
    ///
    /// The smaller side is fully marked with `stamp_side`. Any open edge
    /// leaving it stays inside the same component, hence enters the other
    /// side of the cut; the other side need not be fully explored.
    fn scan_replacement(
        &self,
        side_is_a: bool,
        stamp_side: u32,
        lattice: &CsrLattice,
        open: &[bool],
    ) -> Option<(usize, usize, usize)> {
        let side = if side_is_a {
            &self.side_a
        } else {
            &self.side_b
        };
        for &vertex in side {
            for (neighbor, edge) in lattice.incidences(vertex) {
                if open[edge] && !self.tree_edge[edge] && self.mark[neighbor] != stamp_side {
                    return Some((vertex, neighbor, edge));
                }
            }
        }
        None
    }

    /// Re-root the tree containing `target` at `target`, reversing the
    /// parent/child links along the path to the old root.
    fn reroot(&mut self, target: usize) {
        self.path_vertices.clear();
        self.path_edges.clear();
        let mut current = target;
        self.path_vertices.push(current);
        while self.parent[current] != current {
            let next = self.parent[current];
            self.path_edges.push(self.parent_edge[current]);
            self.path_vertices.push(next);
            current = next;
        }
        let path_len = self.path_vertices.len();
        // `path_vertices[i + 1]` is the old parent of `path_vertices[i]`.
        for index in 0..path_len - 1 {
            let child = self.path_vertices[index];
            let parent = self.path_vertices[index + 1];
            let position = self.children[parent]
                .iter()
                .position(|&candidate| candidate == child)
                .expect("spanning-forest parent/child pointers must agree");
            self.children[parent].swap_remove(position);
        }
        self.parent[self.path_vertices[0]] = self.path_vertices[0];
        self.parent_edge[self.path_vertices[0]] = NONE;
        for index in 1..path_len {
            let parent = self.path_vertices[index - 1];
            let child = self.path_vertices[index];
            self.parent[child] = parent;
            self.parent_edge[child] = self.path_edges[index - 1];
            // The orientation of this tree edge reversed: the old parent is
            // now the child, so `child_of_edge` must follow.
            self.child_of_edge[self.path_edges[index - 1]] = child;
            self.children[parent].push(child);
        }
    }

    fn next_stamps(&mut self) -> (u32, u32) {
        if self.stamp > u32::MAX - 2 {
            self.mark.fill(0);
            self.stamp = 0;
        }
        let first = self.stamp + 1;
        let second = first + 1;
        self.stamp = second;
        (first, second)
    }
}

/// Validation or domain failure of a Sweeny sampler.
#[derive(Debug)]
pub enum SweenyError {
    /// `q` is not finite and strictly positive.
    InvalidQ { q: f64 },
    /// `p` is not finite and strictly inside `(0, 1)`.
    InvalidProbability { probability: f64 },
    /// `p/(1-p)` is not finite and positive.
    NonFiniteOdds { probability: f64, odds: f64 },
    /// The lattice failed `CsrLattice::validate`.
    Topology(String),
    /// The plan itself is invalid for this topology.
    Observable(ObservablePlanError),
    /// A boundary query vertex domain differs from the topology.
    QueryDomain {
        query: usize,
        expected: usize,
        actual: usize,
    },
    /// The plan requests winding observables, which have no FK analysis here.
    WindingUnsupported,
    /// A proposal index is outside the physical edge domain.
    EdgeIndex { index: usize, edge_count: usize },
    /// The component workspace could not be reserved.
    Capacity,
}

impl fmt::Display for SweenyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidQ { q } => write!(
                formatter,
                "random-cluster q must be finite and strictly positive, got {q}"
            ),
            Self::InvalidProbability { probability } => write!(
                formatter,
                "bond probability must be finite and strictly inside (0, 1), got {probability}"
            ),
            Self::NonFiniteOdds { probability, odds } => write!(
                formatter,
                "odds v = p/(1-p) = {odds} is not finite and positive for p = {probability}"
            ),
            Self::Topology(reason) => write!(formatter, "invalid topology: {reason}"),
            Self::Observable(source) => write!(formatter, "invalid observable plan: {source}"),
            Self::QueryDomain {
                query,
                expected,
                actual,
            } => write!(
                formatter,
                "boundary query {query} has V={actual}; adapter topology has V={expected}"
            ),
            Self::WindingUnsupported => write!(
                formatter,
                "winding observables are not part of the Sweeny adapter; call \
                 analyze_winding with an ExternalBondConfiguration instead"
            ),
            Self::EdgeIndex { index, edge_count } => write!(
                formatter,
                "proposal edge index {index} is outside the physical edge domain 0..{edge_count}"
            ),
            Self::Capacity => formatter.write_str("could not reserve component workspace"),
        }
    }
}

impl std::error::Error for SweenyError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_square, Bond, BondType, BoundaryQuery, StaticObservable};
    use rand::SeedableRng;

    fn plain_plan() -> ObservablePlan {
        ObservablePlan::new([StaticObservable::ActiveEdgeCount], Vec::new())
    }

    fn context(seed: u64) -> Context<Xoshiro256PlusPlus> {
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1)
    }

    fn independent_components(lattice: &CsrLattice, mask: &[bool]) -> (usize, Vec<usize>) {
        let vertex_count = lattice.vertex_count();
        let mut parent = (0..vertex_count).collect::<Vec<_>>();
        fn find(parent: &mut [usize], mut node: usize) -> usize {
            while parent[node] != node {
                node = parent[node];
            }
            node
        }
        for (edge, bond) in lattice.edges.iter().enumerate() {
            if mask[edge] {
                let left = find(&mut parent, bond.source);
                let right = find(&mut parent, bond.target);
                if left != right {
                    parent[left] = right;
                }
            }
        }
        let mut count = 0usize;
        let mut labels = Vec::with_capacity(vertex_count);
        for vertex in 0..vertex_count {
            let root = find(&mut parent, vertex);
            if root == vertex {
                count += 1;
            }
            labels.push(root);
        }
        // Normalize labels to the first member of each component for a
        // partition comparison independent of root choice.
        let mut canonical = vec![usize::MAX; vertex_count];
        for label in &mut labels {
            let root = *label;
            if canonical[root] == usize::MAX {
                canonical[root] = root;
            }
            *label = canonical[root];
        }
        (count, labels)
    }

    fn forest_partition(forest: &SpanningForest, vertex_count: usize) -> Vec<usize> {
        (0..vertex_count)
            .map(|vertex| forest.find(vertex))
            .collect()
    }

    fn assert_same_partition(left: &[usize], right: &[usize]) {
        assert_eq!(left.len(), right.len());
        for i in 0..left.len() {
            for j in 0..left.len() {
                assert_eq!(
                    left[i] == left[j],
                    right[i] == right[j],
                    "partitions disagree at ({i}, {j}): {left:?} vs {right:?}"
                );
            }
        }
    }

    fn assert_forest_invariants(mc: &SweenyMC) {
        let forest = &mc.forest;
        let vertex_count = mc.lattice.vertex_count();
        let mut tree_edges = 0usize;
        for (edge, &is_tree) in forest.tree_edge.iter().enumerate() {
            if is_tree {
                tree_edges += 1;
                assert!(mc.open[edge], "tree edge {edge} must be open");
                let child = forest.child_of_edge[edge];
                assert_ne!(child, NONE, "tree edge {edge} must own a child");
                assert_eq!(forest.parent_edge[child], edge);
            } else {
                assert_eq!(forest.child_of_edge[edge], NONE);
            }
        }
        assert_eq!(
            tree_edges,
            vertex_count - mc.component_count,
            "spanning forest must have V - k tree edges"
        );
        for vertex in 0..vertex_count {
            if forest.parent[vertex] == vertex {
                assert_eq!(forest.parent_edge[vertex], NONE);
            } else {
                let parent = forest.parent[vertex];
                assert!(
                    forest.children[parent].contains(&vertex),
                    "vertex {vertex} missing from children of {parent}"
                );
                assert!(forest.tree_edge[forest.parent_edge[vertex]]);
            }
            for &child in &forest.children[vertex] {
                assert_eq!(forest.parent[child], vertex);
            }
        }
    }

    fn fuzz_graph(seed: u64) -> CsrLattice {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        let vertex_count = 9usize;
        let edges = (0..16)
            .map(|_| {
                let source = rng.random_range(0..vertex_count);
                let target = rng.random_range(0..vertex_count);
                Bond::new(source, target, BondType::Generic, 1.0)
            })
            .collect::<Vec<_>>();
        CsrLattice::from_edges(vertex_count, edges)
    }

    #[test]
    fn rates_validate_the_open_domain_and_expose_both_formulas() {
        for q in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                SweenyRates::new(q, 0.5).expect_err("invalid q"),
                SweenyError::InvalidQ { .. }
            ));
        }
        for p in [0.0, 1.0, -0.1, 1.1, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                SweenyRates::new(2.0, p).expect_err("invalid p"),
                SweenyError::InvalidProbability { .. }
            ));
        }
        let rates = SweenyRates::new(2.0, 0.25).expect("valid rates");
        assert_eq!(rates.q(), 2.0);
        assert_eq!(rates.probability(), 0.25);
        assert_eq!(rates.odds(), 1.0 / 3.0);
        // q = 1 heat-bath reduces to the Bernoulli probabilities.
        let unit = SweenyRates::new(1.0, 0.3).expect("valid rates");
        assert!((unit.heat_bath(0, 1) - 0.3).abs() < 1e-15);
        assert!((unit.heat_bath(0, -1) - 0.7).abs() < 1e-15);
        assert!((unit.heat_bath(-1, 1) - 0.3).abs() < 1e-15);
        assert!((unit.heat_bath(1, -1) - 0.7).abs() < 1e-15);
    }

    #[test]
    fn acceptance_ratio_bound_holds_across_q_odds_and_moves() {
        for q in [0.5, 1.0, 2.0, 3.0, 6.0] {
            for probability in [0.05, 0.2, 0.5, 0.8, 0.95] {
                let rates = SweenyRates::new(q, probability).expect("valid rates");
                let odds = rates.odds();
                for delta_n in [-1, 0, 1] {
                    for delta_b in [-1, 1] {
                        let log_ratio = rates.log_ratio(delta_n, delta_b);
                        let heat_bath = rates.heat_bath(delta_n, delta_b);
                        let metropolis = rates.metropolis(delta_n, delta_b);
                        let ratio = rates.heat_bath_metropolis_ratio(delta_n, delta_b);
                        assert!(
                            (0.5 - 1e-12..=1.0).contains(&ratio),
                            "ratio {ratio} outside [1/2, 1] for q={q} p={probability} \
                             dn={delta_n} db={delta_b}"
                        );
                        let closed = if log_ratio >= 0.0 {
                            let x = log_ratio.exp();
                            x / (1.0 + x)
                        } else {
                            let x = log_ratio.exp();
                            1.0 / (1.0 + x)
                        };
                        assert!((ratio - closed).abs() < 1e-14);
                        assert!((heat_bath - metropolis * ratio).abs() < 1e-14);
                        if log_ratio == 0.0 {
                            assert_eq!(ratio, 0.5);
                        }
                        if log_ratio.abs() < 30.0 {
                            assert!(ratio < 1.0);
                        }
                        // Critical odds v = sqrt(q): the minimum moves sit at
                        // the percolation limit q -> 1.
                        if (odds - q.sqrt()).abs() < 1e-12 {
                            assert!(ratio >= 0.5 - 1e-12);
                        }
                    }
                }
            }
        }
        let unit = SweenyRates::new(1.0, 0.5).expect("valid rates");
        // At q = 1, v = 1 the internal insertion/deletion have x = 1 and
        // pHB/pMH = 1/2 exactly (Geyer's percolation-limit statement).
        assert_eq!(unit.heat_bath_metropolis_ratio(0, 1), 0.5);
        assert_eq!(unit.heat_bath_metropolis_ratio(0, -1), 0.5);
    }

    #[test]
    fn constructor_rejections_are_typed() {
        let lattice = build_square(3, 3, false);
        assert!(matches!(
            SweenyMC::new(
                lattice.clone(),
                -1.0,
                0.5,
                SweenyAcceptance::HeatBath,
                plain_plan()
            )
            .err()
            .expect("invalid q"),
            SweenyError::InvalidQ { .. }
        ));
        assert!(matches!(
            SweenyMC::new(
                lattice.clone(),
                2.0,
                1.0,
                SweenyAcceptance::HeatBath,
                plain_plan()
            )
            .err()
            .expect("invalid p"),
            SweenyError::InvalidProbability { .. }
        ));
        let winding = ObservablePlan::new([StaticObservable::Wrapping], Vec::new());
        assert!(matches!(
            SweenyMC::new(
                lattice.clone(),
                2.0,
                0.5,
                SweenyAcceptance::HeatBath,
                winding
            )
            .err()
            .expect("winding"),
            SweenyError::WindingUnsupported
        ));
        let foreign = crate::build_chain(5, false);
        let foreign_query = BoundaryQuery::new(&foreign, &[0], &[4]).expect("chain query");
        let foreign_plan =
            ObservablePlan::new([StaticObservable::LargestSize], vec![foreign_query]);
        assert!(matches!(
            SweenyMC::new(
                lattice.clone(),
                2.0,
                0.5,
                SweenyAcceptance::HeatBath,
                foreign_plan
            )
            .err()
            .expect("query domain"),
            SweenyError::QueryDomain {
                query: 0,
                expected: 9,
                actual: 5
            }
        ));
        let mut mc = SweenyMC::new(lattice, 2.0, 0.5, SweenyAcceptance::HeatBath, plain_plan())
            .expect("valid sampler");
        let mut ctx = context(1);
        assert!(matches!(
            mc.propose(mc.lattice.edge_count(), &mut ctx.rng)
                .expect_err("edge index"),
            SweenyError::EdgeIndex { .. }
        ));
    }

    #[test]
    fn update_fuzz_matches_full_recomputation() {
        for seed in 0..12u64 {
            let lattice = fuzz_graph(seed);
            let q = [0.5, 1.0, 2.0, 3.0][(seed as usize) % 4];
            let probability = [0.2, 0.4, 0.65, 0.85][(seed as usize + 1) % 4];
            let acceptance = if seed % 2 == 0 {
                SweenyAcceptance::HeatBath
            } else {
                SweenyAcceptance::Metropolis
            };
            let mut mc = SweenyMC::new(lattice.clone(), q, probability, acceptance, plain_plan())
                .expect("fuzz sampler");
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed ^ 0x5eed);
            for step in 0..2400usize {
                let edge = rng.random_range(0..lattice.edge_count());
                mc.propose(edge, &mut rng).expect("in-domain proposal");
                let (components, labels) = independent_components(&lattice, &mc.open);
                assert_eq!(mc.component_count, components, "seed {seed} step {step}");
                assert_eq!(
                    mc.open_count,
                    mc.open.iter().filter(|&&open| open).count(),
                    "seed {seed} step {step}"
                );
                let forest_labels = forest_partition(&mc.forest, lattice.vertex_count());
                assert_same_partition(&forest_labels, &labels);
                assert_forest_invariants(&mc);
            }
        }
    }

    #[test]
    fn reset_returns_to_the_all_closed_state() {
        let lattice = build_square(4, 4, false);
        let mut mc = SweenyMC::new(
            lattice.clone(),
            2.0,
            0.6,
            SweenyAcceptance::HeatBath,
            plain_plan(),
        )
        .expect("valid sampler");
        let mut ctx = context(7);
        for _ in 0..10 {
            mc.sweep(&mut ctx);
        }
        assert!(mc.open_edge_count() > 0);
        mc.reset();
        assert_eq!(mc.open_edge_count(), 0);
        assert_eq!(mc.component_count(), lattice.vertex_count());
        assert_eq!(mc.proposal_count(), 0);
        assert!(mc.mask().iter().all(|&open| !open));
    }
}
