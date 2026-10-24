//! Integer-`q` Chayes-Machta sampler: activate/sleep plus per-cluster colour
//! resampling (roadmap v2 §6 Q2).
//!
//! # The algorithm
//!
//! Chayes and Machta introduced a cluster dynamics for the random-cluster
//! model that generalizes Swendsen-Wang to non-integer `q` by splitting the
//! cluster weight into colour weights `q = sum_alpha q_alpha`. For **integer**
//! `q >= 1` each of the `q` colours has weight one, the joint measure over a
//! bond configuration `A` and a colouring `sigma` is
//!
//! ```text
//! mu_{p,q}(A, sigma) ∝ 1{sigma constant on each cluster of A}
//!                       * p^{|A|} (1-p)^{|E|-|A|}
//!                       * prod_alpha q_alpha^{k_alpha(sigma, A)},
//! ```
//!
//! and summing over the `q` colours of each of the `k(A)` clusters gives the
//! random-cluster marginal `q^{k(A)}`. Geyer's archived review names the
//! dynamics "the Chayes-Machta-Swendsen-Wang dynamics [12, 21]"
//! (`docs/ref/md/Geyer_2015_sweeny_random_cluster.md`, line 169), with
//! reference [21] = L. Chayes and J. Machta, Physica A **254**, 477 (1998).
//! Garoni, Ossola, Polin and Sokal (arXiv:1105.0373, Sec. 2) state the joint
//! measure Eq. 2.9 and the two half-steps:
//!
//! > "Given a bond configuration, choose a new colour configuration ... for
//! > each alpha in S update the random-cluster model on G[V_alpha] with
//! > parameter q_alpha" (Algorithm 1), and the colour resampling "randomly
//! > choose a value alpha in S with probability q_alpha/q and impose it on
//! > all vertices in that component."
//!
//! For integer `q` the single-colour update is the CM activate/sleep step:
//!
//! 1. **Colour resampling** — every current bond cluster draws one colour
//!    uniformly from `{0, ..., q-1}` (weight `q_alpha/q = 1/q`). Clusters of
//!    the selected **active** colour are awake; all others sleep.
//! 2. **Activate / sleep bond update** — remove every open edge whose two
//!    endpoints are awake, then re-add each edge between awake vertices
//!    independently with probability `p`. Sleeping edges are untouched.
//!
//! Repeating this with the active colour cycling through `0..q` is a Gibbs
//! sampler on the joint measure: step 1 is the exact colour conditional, step
//! 2 the exact bond conditional inside one colour class. Its bond marginal is
//! therefore the random-cluster measure (the same anchor as Q1), and at
//! `q = 2` it is the same marginal as R3's Edwards-Sokal/Swendsen-Wang chain.
//! At `q = 1` every vertex is awake, so a sweep redraws **every** edge as an
//! independent Bernoulli(p) draw — the R2 reduction, bitwise.
//!
//! # RNG contract (documented and replayed by the tests)
//!
//! One `sweep()` consumes variates in this order:
//!
//! 1. **Union-find rebuild** over the open edges (no variate).
//! 2. **Per-cluster colour pass** in ascending minimum-vertex order: the first
//!    vertex of each cluster draws one `usize` uniform in `0..q`; the whole
//!    cluster receives it. `q = 1` consumes no variate (the only colour).
//! 3. **Active colour**: sweep number `t` uses colour `t mod q` (no variate).
//! 4. **Awake-edge pass** in ascending dense `EdgeId` order: an edge with both
//!    endpoints in the active colour consumes exactly one `f64`
//!    (`rng.random::<f64>() < p`); every other edge consumes nothing.
//!
//! The constructor requires `p` strictly inside `(0, 1)`, so there are no
//! endpoint shortcuts to replay.
//!
//! # Cost and scope
//!
//! Unlike Q1's Sweeny chain this is a global cluster sweep: each `sweep()`
//! rebuilds a union-find over all open edges (`O(V + E)`) and re-draws the
//! awake subgraph. That is the standard Swendsen-Wang/CM cost and is not the
//! Q1 dynamic-connectivity problem. All buffers are reused; a steady-state
//! sweep allocates nothing, and [`ChayesMachtaMC::realization`] allocates
//! on demand.
//!
//! This is a Markov chain: thermalize before measuring. Winding observables
//! are a typed rejection (no embedding), like the other FK adapters.

use core::fmt;

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::{analyze_bonds, BondSample, ComponentWorkspace, ObservablePlan, ObservablePlanError};
use crate::{CsrLattice, GraphView, UndirectedGraphView};

const UNASSIGNED: usize = usize::MAX;

/// Integer-`q` Chayes-Machta adapter with activate/sleep cluster updates.
///
/// `q` is a positive integer (`q = 1` is the Bernoulli reduction), `p` is
/// strictly inside `(0, 1)`. The chain starts all-closed with every vertex in
/// colour `0`; the first sweep resamples colours and bonds. The current bond
/// mask is exposed through [`Self::mask`] and [`Self::realization`], and
/// [`MonteCarlo::measure`] analyzes that same mask through the shared F3
/// analyzer and the [`ObservablePlan`].
pub struct ChayesMachtaMC {
    lattice: CsrLattice,
    q: usize,
    probability: f64,
    open: Vec<bool>,
    colors: Vec<usize>,
    parent: Vec<usize>,
    rank: Vec<u8>,
    root_color: Vec<usize>,
    open_count: usize,
    component_count: usize,
    sweep_index: usize,
    workspace: ComponentWorkspace,
    plan: ObservablePlan,
}

impl ChayesMachtaMC {
    /// Construct an all-closed integer-`q` Chayes-Machta chain.
    ///
    /// Rejections: `q = 0`, `p` outside the open unit interval, invalid
    /// topology, winding plans, invalid plans for this vertex count, and
    /// boundary queries with a foreign vertex domain.
    pub fn new(
        lattice: CsrLattice,
        q: usize,
        probability: f64,
        plan: ObservablePlan,
    ) -> Result<Self, ChayesMachtaError> {
        lattice
            .validate()
            .map_err(|error| ChayesMachtaError::Topology(error.to_string()))?;
        if q == 0 {
            return Err(ChayesMachtaError::ZeroColors);
        }
        if !probability.is_finite() || probability <= 0.0 || probability >= 1.0 {
            return Err(ChayesMachtaError::InvalidProbability { probability });
        }
        if plan.requests_winding() {
            return Err(ChayesMachtaError::WindingUnsupported);
        }
        plan.validate(lattice.vertex_count())
            .map_err(ChayesMachtaError::Observable)?;
        for (query, boundary) in plan.boundary_queries().iter().enumerate() {
            if boundary.vertex_count() != lattice.vertex_count() {
                return Err(ChayesMachtaError::QueryDomain {
                    query,
                    expected: lattice.vertex_count(),
                    actual: boundary.vertex_count(),
                });
            }
        }
        let mut workspace = ComponentWorkspace::new();
        workspace
            .prepare(lattice.vertex_count(), plan.boundary_queries().len(), false)
            .map_err(|_| ChayesMachtaError::Capacity)?;
        let vertex_count = lattice.vertex_count();
        let edge_count = lattice.edge_count();
        Ok(Self {
            lattice,
            q,
            probability,
            open: vec![false; edge_count],
            colors: vec![0; vertex_count],
            parent: (0..vertex_count).collect(),
            rank: vec![0; vertex_count],
            root_color: vec![UNASSIGNED; vertex_count],
            open_count: 0,
            component_count: vertex_count,
            sweep_index: 0,
            workspace,
            plan,
        })
    }

    /// The owned topology.
    pub const fn lattice(&self) -> &CsrLattice {
        &self.lattice
    }

    /// The integer colour count `q`.
    pub const fn q(&self) -> usize {
        self.q
    }

    /// The bond probability `p`.
    pub const fn probability(&self) -> f64 {
        self.probability
    }

    /// The current dense bond mask.
    #[inline]
    pub fn mask(&self) -> &[bool] {
        &self.open
    }

    /// The current per-vertex colour (`0..q`).
    #[inline]
    pub fn colors(&self) -> &[usize] {
        &self.colors
    }

    /// Number of open bonds.
    #[inline]
    pub const fn open_edge_count(&self) -> usize {
        self.open_count
    }

    /// Number of connected components of the current bond configuration.
    #[inline]
    pub const fn component_count(&self) -> usize {
        self.component_count
    }

    /// The active colour of the next sweep.
    #[inline]
    pub const fn next_active_color(&self) -> usize {
        self.sweep_index
    }

    /// The observable plan recorded by `measure()`.
    pub const fn observable_plan(&self) -> &ObservablePlan {
        &self.plan
    }

    /// The current configuration as R1's owned open-edge record.
    pub fn realization(&self) -> BondSample {
        BondSample::from_mask(&self.lattice, &self.open)
            .expect("mask length pinned to the topology at construction")
    }

    /// One Chayes-Machta sweep (colour resampling + awake bond update).
    ///
    /// Exposed as an inherent method so unit tests can drive a sequence of
    /// sweeps without a scheduler context; [`MonteCarlo::sweep`] delegates
    /// here.
    pub fn step<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        self.rebuild_components();
        self.assign_colors(rng);
        let active = self.sweep_index;
        self.sweep_index = (self.sweep_index + 1) % self.q;
        for edge in 0..self.lattice.edge_count() {
            let bond = self.lattice.edges[edge];
            if self.colors[bond.source] == active && self.colors[bond.target] == active {
                // Awake edge: remove and re-add as an independent Bernoulli(p)
                // draw. One uniform per awake edge.
                self.open[edge] = rng.random::<f64>() < self.probability;
            }
        }
        self.rebuild_components();
    }

    /// Rebuild the union-find and the component count from the open mask.
    fn rebuild_components(&mut self) {
        let vertex_count = self.lattice.vertex_count();
        for vertex in 0..vertex_count {
            self.parent[vertex] = vertex;
            self.rank[vertex] = 0;
        }
        for edge in 0..self.lattice.edge_count() {
            if self.open[edge] {
                let bond = self.lattice.edges[edge];
                union(&mut self.parent, &mut self.rank, bond.source, bond.target);
            }
        }
        let mut components = 0usize;
        for vertex in 0..vertex_count {
            if find(&mut self.parent, vertex) == vertex {
                components += 1;
            }
        }
        self.component_count = components;
        self.open_count = self.open.iter().filter(|&&open| open).count();
    }

    /// Resample one colour per current component, in ascending minimum-vertex
    /// order, and mirror it into `colors`.
    fn assign_colors<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        self.root_color.fill(UNASSIGNED);
        let vertex_count = self.lattice.vertex_count();
        for vertex in 0..vertex_count {
            let root = find(&mut self.parent, vertex);
            if self.root_color[root] == UNASSIGNED {
                self.root_color[root] = if self.q == 1 {
                    0
                } else {
                    rng.random_range(0..self.q)
                };
            }
            self.colors[vertex] = self.root_color[root];
        }
    }
}

impl MonteCarlo for ChayesMachtaMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.step(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let result = analyze_bonds(
            &self.lattice,
            &self.open,
            self.plan.boundary_queries(),
            &mut self.workspace,
        )
        .expect("constructor prepared matching analysis domains");
        debug_assert_eq!(
            result.summary().component_count,
            self.component_count,
            "union-find count must match the analyzer"
        );
        self.plan.record(
            self.lattice.vertex_count(),
            result.summary(),
            result.query_outcomes(),
            None,
            ctx,
        );
    }

    fn name(&self) -> &'static str {
        "ChayesMachta"
    }
}

fn find(parent: &mut [usize], mut node: usize) -> usize {
    let mut root = node;
    while parent[root] != root {
        root = parent[root];
    }
    while parent[node] != root {
        let next = parent[node];
        parent[node] = root;
        node = next;
    }
    root
}

fn union(parent: &mut [usize], rank: &mut [u8], left: usize, right: usize) {
    let left_root = find(parent, left);
    let right_root = find(parent, right);
    if left_root == right_root {
        return;
    }
    match rank[left_root].cmp(&rank[right_root]) {
        core::cmp::Ordering::Less => parent[left_root] = right_root,
        core::cmp::Ordering::Greater => parent[right_root] = left_root,
        core::cmp::Ordering::Equal => {
            parent[right_root] = left_root;
            rank[left_root] += 1;
        }
    }
}

/// Construction or validation failure of a Chayes-Machta adapter.
#[derive(Debug)]
pub enum ChayesMachtaError {
    /// `q = 0`: at least one colour is required.
    ZeroColors,
    /// `p` is not finite and strictly inside `(0, 1)`.
    InvalidProbability { probability: f64 },
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
    /// The component workspace could not be reserved.
    Capacity,
}

impl fmt::Display for ChayesMachtaError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroColors => formatter.write_str("Chayes-Machta requires an integer q >= 1"),
            Self::InvalidProbability { probability } => write!(
                formatter,
                "bond probability must be finite and strictly inside (0, 1), got {probability}"
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
                "winding observables are not part of the Chayes-Machta adapter; call \
                 analyze_winding with an ExternalBondConfiguration instead"
            ),
            Self::Capacity => formatter.write_str("could not reserve component workspace"),
        }
    }
}

impl std::error::Error for ChayesMachtaError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_square, BoundaryQuery, StaticObservable};
    use rand::SeedableRng;

    fn plain_plan() -> ObservablePlan {
        ObservablePlan::new([StaticObservable::ActiveEdgeCount], Vec::new())
    }

    fn reference_components(lattice: &CsrLattice, mask: &[bool]) -> (usize, Vec<usize>) {
        let vertex_count = lattice.vertex_count();
        let mut parent = (0..vertex_count).collect::<Vec<_>>();
        fn find_local(parent: &mut [usize], mut node: usize) -> usize {
            while parent[node] != node {
                node = parent[node];
            }
            node
        }
        for (edge, bond) in lattice.edges.iter().enumerate() {
            if mask[edge] {
                let left = find_local(&mut parent, bond.source);
                let right = find_local(&mut parent, bond.target);
                if left != right {
                    parent[left] = right;
                }
            }
        }
        let mut labels = Vec::with_capacity(vertex_count);
        let mut count = 0usize;
        for vertex in 0..vertex_count {
            let root = find_local(&mut parent, vertex);
            if root == vertex {
                count += 1;
            }
            labels.push(root);
        }
        (count, labels)
    }

    #[test]
    fn constructor_rejections_are_typed() {
        let lattice = build_square(3, 3, false);
        assert!(matches!(
            ChayesMachtaMC::new(lattice.clone(), 0, 0.5, plain_plan())
                .err()
                .expect("zero colours"),
            ChayesMachtaError::ZeroColors
        ));
        for p in [0.0, 1.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                ChayesMachtaMC::new(lattice.clone(), 2, p, plain_plan())
                    .err()
                    .expect("invalid p"),
                ChayesMachtaError::InvalidProbability { .. }
            ));
        }
        let winding = ObservablePlan::new([StaticObservable::Wrapping], Vec::new());
        assert!(matches!(
            ChayesMachtaMC::new(lattice.clone(), 2, 0.5, winding)
                .err()
                .expect("winding"),
            ChayesMachtaError::WindingUnsupported
        ));
        let foreign = crate::build_chain(5, false);
        let foreign_query = BoundaryQuery::new(&foreign, &[0], &[4]).expect("chain query");
        let foreign_plan =
            ObservablePlan::new([StaticObservable::LargestSize], vec![foreign_query]);
        assert!(matches!(
            ChayesMachtaMC::new(lattice, 2, 0.5, foreign_plan)
                .err()
                .expect("query domain"),
            ChayesMachtaError::QueryDomain {
                query: 0,
                expected: 9,
                actual: 5
            }
        ));
    }

    #[test]
    fn sweeps_preserve_monochromatic_components_and_counts() {
        for seed in 0..8u64 {
            let lattice = build_square(4, 4, false);
            let q = 1 + (seed as usize % 4);
            let mut mc =
                ChayesMachtaMC::new(lattice.clone(), q, 0.4, plain_plan()).expect("valid adapter");
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed ^ 0xc4a5);
            for _ in 0..40 {
                mc.step(&mut rng);
                let (components, labels) = reference_components(&lattice, &mc.open);
                assert_eq!(mc.component_count, components);
                assert_eq!(mc.open_count, mc.open.iter().filter(|&&open| open).count());
                // Every open edge joins equal colours and every component is
                // monochromatic.
                for (edge, bond) in lattice.edges.iter().enumerate() {
                    if mc.open[edge] {
                        assert_eq!(
                            mc.colors[bond.source], mc.colors[bond.target],
                            "open edge {edge} joins different colours"
                        );
                    }
                }
                let mut colors_by_component = std::collections::HashMap::new();
                #[allow(clippy::needless_range_loop)]
                for vertex in 0..lattice.vertex_count() {
                    let root = labels[vertex];
                    match colors_by_component.entry(root) {
                        std::collections::hash_map::Entry::Vacant(entry) => {
                            entry.insert(mc.colors[vertex]);
                        }
                        std::collections::hash_map::Entry::Occupied(entry) => {
                            assert_eq!(*entry.get(), mc.colors[vertex]);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn q1_sweep_is_a_fresh_bernoulli_draw_in_dense_order() {
        let lattice = build_square(4, 5, false);
        let p = 0.37;
        let mut mc = ChayesMachtaMC::new(
            lattice.clone(),
            1,
            p,
            ObservablePlan::new([StaticObservable::ActiveEdgeCount], Vec::new()),
        )
        .expect("valid adapter");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x51ee);
        let expected = (0..lattice.edge_count())
            .map(|_| rng.random::<f64>() < p)
            .collect::<Vec<_>>();
        let mut actual = Xoshiro256PlusPlus::seed_from_u64(0x51ee);
        mc.step(&mut actual);
        assert_eq!(mc.open, expected);
        // A second q = 1 sweep is again a fresh draw, independent of the
        // previous configuration.
        let expected_second = (0..lattice.edge_count())
            .map(|_| rng.random::<f64>() < p)
            .collect::<Vec<_>>();
        mc.step(&mut actual);
        assert_eq!(mc.open, expected_second);
    }
}
