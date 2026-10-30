//! X2: two-layer mutual-component (cascading-failure) analysis.
//!
//! Model and definitional anchor: the interdependent-network order parameter
//! is the fraction of vertices in the **mutually connected largest
//! component**; cascading failures were introduced by Buldyrev et al.,
//! Nature 464, 1025 (2010) and are summarized in the archived Saberi review
//! (`docs/ref/md/Saberi_2015_percolation_review.md`, §V: "the order
//! parameter is taken as the fraction of nodes in the mutually connected
//! largest component", and the cascade mechanism). This module implements
//! the one-to-one dependency case at full initial occupancy:
//!
//! 1. Keep the largest component of layer A among the currently alive
//!    vertices.
//! 2. Keep the largest component of layer B among those survivors.
//! 3. Repeat until the alive set stops shrinking.
//!
//! Both layers share one dense vertex set; layer edges are independent. A
//! vertex survives exactly when it belongs to both layers' surviving giants.
//! Ties in the largest component use the analyzer's canonical rule (the
//! lowest vertex ID), so the result is deterministic. Random initial
//! failures and more than two layers are outside X2's first slice.

use carlo_rs::{Context, MonteCarlo};
use rand::{Rng, RngExt};
use rand_xoshiro::Xoshiro256PlusPlus;

use crate::percolation::{Probability, ProbabilityError};
use crate::{
    analyze_with_labels, Bond, BondType, ComponentWorkspace, CsrLattice, EdgeActivity, GraphView,
    UndirectedGraphView, VertexActivity,
};

/// Two network layers over one shared dense vertex set.
#[derive(Clone, Debug)]
pub struct MultiplexNetwork {
    layer_a: CsrLattice,
    layer_b: CsrLattice,
}

/// Fixed point of the two-layer cascade.
#[derive(Clone, Debug, PartialEq)]
pub struct MutualComponent {
    active: Vec<bool>,
    iterations: usize,
}

/// Typed multiplex rejection.
#[derive(Clone, Debug, PartialEq)]
pub enum MultiplexError {
    /// The two layers cover different vertex counts.
    VertexCount { first: usize, second: usize },
    /// A network with no vertices has no mutual component to analyze.
    EmptyNetwork,
    /// An input edge endpoint lies outside the dense vertex range.
    Endpoint {
        layer: &'static str,
        source: usize,
        target: usize,
        vertices: usize,
    },
    /// The monotone cascade failed to stabilize within `2V + 2` contractions.
    Cascade { steps: usize },
    /// An initial-activity mask does not match the vertex count.
    ActivityDomain { expected: usize, actual: usize },
    /// The initial survival probability is outside `[0, 1]`.
    Probability(ProbabilityError),
}

impl core::fmt::Display for MultiplexError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::VertexCount { first, second } => write!(
                formatter,
                "layers cover different vertex counts ({first} vs {second})"
            ),
            Self::EmptyNetwork => write!(formatter, "multiplex network has no vertices"),
            Self::Endpoint {
                layer,
                source,
                target,
                vertices,
            } => write!(
                formatter,
                "layer {layer}: edge {source}-{target} outside 0..{vertices}"
            ),
            Self::Cascade { steps } => write!(
                formatter,
                "two-layer cascade did not stabilize after {steps} contractions"
            ),
            Self::ActivityDomain { expected, actual } => write!(
                formatter,
                "initial activity covers {actual} vertices; expected {expected}"
            ),
            Self::Probability(source) => {
                write!(formatter, "invalid survival probability: {source}")
            }
        }
    }
}

impl std::error::Error for MultiplexError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Probability(source) => Some(source),
            _ => None,
        }
    }
}

impl MultiplexNetwork {
    /// Build the two-layer network; both layers must cover the same dense
    /// vertex set.
    pub fn new(layer_a: CsrLattice, layer_b: CsrLattice) -> Result<Self, MultiplexError> {
        if layer_a.vertex_count() != layer_b.vertex_count() {
            return Err(MultiplexError::VertexCount {
                first: layer_a.vertex_count(),
                second: layer_b.vertex_count(),
            });
        }
        if layer_a.vertex_count() == 0 {
            return Err(MultiplexError::EmptyNetwork);
        }
        Ok(Self { layer_a, layer_b })
    }

    pub fn vertices(&self) -> usize {
        self.layer_a.vertex_count()
    }

    pub const fn layer_a(&self) -> &CsrLattice {
        &self.layer_a
    }

    pub const fn layer_b(&self) -> &CsrLattice {
        &self.layer_b
    }

    /// Run the cascade to its fixed point and report the mutual component.
    ///
    /// ```
    /// use cmc_rs::multiplex_from_edges;
    ///
    /// let network = multiplex_from_edges(
    ///     4,
    ///     [(0, 1), (1, 2), (2, 3)],
    ///     [(0, 1), (2, 3)],
    /// )
    /// .expect("two layers");
    /// let mutual = network.mutual_component().expect("cascade");
    /// assert_eq!(mutual.active(), &[true, true, false, false]);
    /// assert_eq!(mutual.fraction(), 0.5);
    /// ```
    pub fn mutual_component(&self) -> Result<MutualComponent, MultiplexError> {
        let alive = vec![true; self.vertices()];
        self.mutual_component_from(&alive)
    }

    /// Run the cascade from an explicit initial-activity mask.
    ///
    /// `initial[i] == true` means vertex `i` survived the primary failure;
    /// vertices outside the mask are never revived. The mask must cover the
    /// shared vertex set.
    pub fn mutual_component_from(
        &self,
        initial: &[bool],
    ) -> Result<MutualComponent, MultiplexError> {
        let vertices = self.vertices();
        if initial.len() != vertices {
            return Err(MultiplexError::ActivityDomain {
                expected: vertices,
                actual: initial.len(),
            });
        }
        let mut workspace = ComponentWorkspace::new();
        let mut alive = initial.to_vec();
        let mut steps = 0usize;
        let max_steps = 2 * vertices + 2;
        loop {
            let before = alive.clone();
            for layer in [&self.layer_a, &self.layer_b] {
                let survivors = largest_component(layer, &alive, &mut workspace);
                alive = vec![false; vertices];
                for survivor in survivors {
                    alive[survivor] = true;
                }
                steps += 1;
            }
            if alive == before {
                return Ok(MutualComponent {
                    active: alive,
                    iterations: steps,
                }
                .validate(vertices));
            }
            if steps > max_steps {
                return Err(MultiplexError::Cascade { steps });
            }
        }
    }
}

impl MutualComponent {
    fn validate(self, vertices: usize) -> Self {
        debug_assert_eq!(self.active.len(), vertices);
        self
    }

    /// Vertices in the mutual component.
    pub fn active_vertices(&self) -> usize {
        self.active.iter().filter(|&&active| active).count()
    }

    /// `active_vertices / V`.
    pub fn fraction(&self) -> f64 {
        let vertices = self.active.len();
        if vertices == 0 {
            0.0
        } else {
            self.active_vertices() as f64 / vertices as f64
        }
    }

    /// Per-vertex membership (dense vertex order).
    pub fn active(&self) -> &[bool] {
        &self.active
    }

    /// Layer contractions performed before the fixed point.
    pub fn iterations(&self) -> usize {
        self.iterations
    }
}

#[derive(Clone, Copy)]
struct AliveVertices<'a>(&'a [bool]);

impl VertexActivity for AliveVertices<'_> {
    fn vertex_count(&self) -> usize {
        self.0.len()
    }

    fn vertex_active(&self, vertex: crate::VertexId) -> bool {
        self.0[vertex.index()]
    }
}

struct AliveLayerEdges<'a> {
    lattice: &'a CsrLattice,
    alive: &'a [bool],
}

impl EdgeActivity for AliveLayerEdges<'_> {
    fn edge_count(&self) -> usize {
        self.lattice.edge_count()
    }

    fn edge_active(&self, edge: crate::EdgeId) -> bool {
        let [source, target] = self.lattice.edge_endpoints(edge);
        self.alive[source.index()] && self.alive[target.index()]
    }
}

/// Dense indices of the layer's largest component among `alive` vertices.
fn largest_component(
    layer: &CsrLattice,
    alive: &[bool],
    workspace: &mut ComponentWorkspace,
) -> Vec<usize> {
    let result = analyze_with_labels(
        layer,
        &AliveVertices(alive),
        &AliveLayerEdges {
            lattice: layer,
            alive,
        },
        &[],
        workspace,
    )
    .expect("cascade activity domains match the layer topology");
    let Some(root) = result.summary().largest_component else {
        return Vec::new();
    };
    result
        .labels()
        .expect("labels were requested")
        .iter()
        .enumerate()
        .filter_map(|(vertex, label)| (*label == Some(root)).then_some(vertex))
        .collect()
}

/// Scheduler adapter: random primary failures followed by the two-layer
/// mutual-component cascade.
///
/// One `sweep()` draws one survival mask (one `f64` per vertex in ascending
/// dense order; vertex survives iff `draw < survival`), runs the cascade to
/// its fixed point, and stores the outcome. One `measure()` records
/// `MultiplexMutualFraction` (mutual vertices / total vertices),
/// `MultiplexMutualVertices`, `MultiplexInitialVertices`, and
/// `MultiplexCascadeIterations`.
///
/// Configure `thermalization_sweeps = 0`. The cascade allocates its working
/// vectors per sweep (the analysis layer predates a reusable cascade
/// workspace); this is documented rather than hidden.
pub struct MultiplexCascadeMC {
    network: MultiplexNetwork,
    survival: Probability,
    initial: Vec<bool>,
    mutual_vertices: usize,
    initial_vertices: usize,
    iterations: usize,
}

impl MultiplexCascadeMC {
    /// Build the adapter; the network must be non-empty and `survival` valid.
    pub fn new(network: MultiplexNetwork, survival: f64) -> Result<Self, MultiplexError> {
        let survival = Probability::new(survival).map_err(MultiplexError::Probability)?;
        let initial = vec![false; network.vertices()];
        Ok(Self {
            network,
            survival,
            initial,
            mutual_vertices: 0,
            initial_vertices: 0,
            iterations: 0,
        })
    }

    pub const fn network(&self) -> &MultiplexNetwork {
        &self.network
    }

    /// Survival probability of the primary failure step.
    pub const fn survival(&self) -> f64 {
        self.survival.get()
    }

    /// Mutual-component size of the last sweep.
    pub const fn mutual_vertices(&self) -> usize {
        self.mutual_vertices
    }

    /// Primary survivors of the last sweep.
    pub const fn initial_vertices(&self) -> usize {
        self.initial_vertices
    }

    /// Cascade contractions of the last sweep (`2` per full round).
    pub const fn iterations(&self) -> usize {
        self.iterations
    }

    /// Draw the primary failure mask and run one complete cascade.
    ///
    /// Draw contract: one `f64` per vertex in ascending dense order, one
    /// realization per call. Exposed for deterministic replay.
    pub fn realize<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        let survival = self.survival.get();
        for slot in self.initial.iter_mut() {
            *slot = rng.random::<f64>() < survival;
        }
        self.initial_vertices = self.initial.iter().filter(|&&alive| alive).count();
        let mutual = self
            .network
            .mutual_component_from(&self.initial)
            .expect("constructor validated the shared vertex domain");
        self.mutual_vertices = mutual.active_vertices();
        self.iterations = mutual.iterations();
    }
}

impl MonteCarlo for MultiplexCascadeMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.realize(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let vertices = self.network.vertices();
        ctx.measure(
            "MultiplexMutualFraction",
            self.mutual_vertices as f64 / vertices as f64,
        );
        ctx.measure("MultiplexMutualVertices", self.mutual_vertices as f64);
        ctx.measure("MultiplexInitialVertices", self.initial_vertices as f64);
        ctx.measure("MultiplexCascadeIterations", self.iterations as f64);
    }

    fn name(&self) -> &'static str {
        "MultiplexCascade"
    }
}

/// Convenience: two-layer network from raw edge lists.
///
/// Edges are undirected and emitted in input order; duplicate and self-loop
/// policy is the caller's (both are valid for the multiplex analysis).
pub fn multiplex_from_edges(
    vertices: usize,
    layer_a: impl IntoIterator<Item = (usize, usize)>,
    layer_b: impl IntoIterator<Item = (usize, usize)>,
) -> Result<MultiplexNetwork, MultiplexError> {
    let build = |edges: Vec<(usize, usize)>, layer: &'static str| {
        for &(source, target) in &edges {
            if source >= vertices || target >= vertices {
                return Err(MultiplexError::Endpoint {
                    layer,
                    source,
                    target,
                    vertices,
                });
            }
        }
        let bonds: Vec<Bond> = edges
            .into_iter()
            .map(|(source, target)| {
                Bond::new(
                    source.min(target),
                    source.max(target),
                    BondType::Generic,
                    1.0,
                )
            })
            .collect();
        Ok(CsrLattice::from_edges(vertices, bonds))
    };
    MultiplexNetwork::new(
        build(layer_a.into_iter().collect(), "A")?,
        build(layer_b.into_iter().collect(), "B")?,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn network(
        vertices: usize,
        layer_a: impl IntoIterator<Item = (usize, usize)>,
        layer_b: impl IntoIterator<Item = (usize, usize)>,
    ) -> MultiplexNetwork {
        multiplex_from_edges(vertices, layer_a, layer_b).expect("network")
    }

    fn active(network: &MultiplexNetwork) -> Vec<bool> {
        network
            .mutual_component()
            .expect("cascade")
            .active()
            .to_vec()
    }

    #[test]
    fn identical_chains_keep_every_vertex() {
        let network = network(
            5,
            [(0, 1), (1, 2), (2, 3), (3, 4)],
            [(0, 1), (1, 2), (2, 3), (3, 4)],
        );
        let mutual = network.mutual_component().expect("cascade");
        assert_eq!(mutual.active_vertices(), 5);
        assert_eq!(mutual.fraction(), 1.0);
        assert_eq!(mutual.iterations(), 2);
    }

    #[test]
    fn a_split_second_layer_cuts_the_chain() {
        // A is a path; B splits it into two pairs, keeping the canonical
        // lowest-ID largest component {0,1}.
        let network = network(4, [(0, 1), (1, 2), (2, 3)], [(0, 1), (2, 3)]);
        let mutual = network.mutual_component().expect("cascade");
        assert_eq!(mutual.active(), &[true, true, false, false]);
        assert_eq!(mutual.fraction(), 0.5);
    }

    #[test]
    fn a_star_and_a_disjoint_path_cascade_to_one_vertex() {
        // A: star centered at 0. B: 0 isolated, 1-2-3-4 a path. The cascade
        // passes through the full star, then the path, then the isolated
        // lowest-ID leaf.
        let network = network(
            5,
            [(0, 1), (0, 2), (0, 3), (0, 4)],
            [(1, 2), (2, 3), (3, 4)],
        );
        let mutual = network.mutual_component().expect("cascade");
        assert_eq!(mutual.active(), &[false, true, false, false, false]);
        assert_eq!(mutual.fraction(), 0.2);
    }

    /// Independent reference: union-find per layer on induced subgraphs, the
    /// same alternating map, no analyzer.
    fn reference(vertices: usize, a: &[(usize, usize)], b: &[(usize, usize)]) -> Vec<bool> {
        fn giant(vertices: usize, edges: &[(usize, usize)], alive: &[bool]) -> Vec<bool> {
            fn root(parent: &mut [usize], mut x: usize) -> usize {
                while parent[x] != x {
                    parent[x] = parent[parent[x]];
                    x = parent[x];
                }
                x
            }
            let mut parent: Vec<usize> = (0..vertices).collect();
            for &(left, right) in edges {
                if !alive[left] || !alive[right] {
                    continue;
                }
                let left = root(&mut parent, left);
                let right = root(&mut parent, right);
                if left != right {
                    parent[left.max(right)] = left.min(right);
                }
            }
            let mut sizes = vec![0usize; vertices];
            for (vertex, &is_alive) in alive.iter().enumerate() {
                if is_alive {
                    let root = root(&mut parent, vertex);
                    sizes[root] += 1;
                }
            }
            let best = (0..vertices)
                .filter(|&vertex| alive[vertex] && sizes[vertex] > 0)
                .max_by_key(|&vertex| (sizes[vertex], std::cmp::Reverse(vertex)));
            let Some(best) = best else {
                return vec![false; vertices];
            };
            (0..vertices)
                .map(|vertex| alive[vertex] && root(&mut parent, vertex) == best)
                .collect()
        }

        let mut alive = vec![true; vertices];
        loop {
            let before = alive.clone();
            alive = giant(vertices, a, &alive);
            alive = giant(vertices, b, &alive);
            if alive == before {
                return alive;
            }
        }
    }

    #[test]
    fn exact_two_layer_enumeration_matches_the_reference_map() {
        type Edges = [(usize, usize)];
        let cases: [(&Edges, &Edges); 4] = [
            (&[(0, 1), (1, 2), (2, 3), (3, 0)], &[(0, 2), (1, 3)]),
            (
                &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)],
                &[(0, 5), (1, 2), (3, 4)],
            ),
            (&[(0, 1), (2, 3), (4, 5)], &[(1, 2), (3, 4), (5, 0)]),
            (
                &[(0, 1), (0, 2), (0, 3), (0, 4)],
                &[(1, 2), (2, 3), (3, 4), (4, 1)],
            ),
        ];
        for (case, (a, b)) in cases.into_iter().enumerate() {
            let network = network(6, a.iter().copied(), b.iter().copied());
            assert_eq!(
                active(&network),
                reference(6, a, b),
                "case {case} disagrees with the independent reference"
            );
        }
    }

    #[test]
    fn rejections_are_typed() {
        let small = CsrLattice::from_edges(3, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);
        let large = CsrLattice::from_edges(4, vec![Bond::new(0, 1, BondType::Generic, 1.0)]);
        assert!(matches!(
            MultiplexNetwork::new(small, large),
            Err(MultiplexError::VertexCount {
                first: 3,
                second: 4
            })
        ));
        assert!(matches!(
            multiplex_from_edges(3, [(0, 3)], std::iter::empty()),
            Err(MultiplexError::Endpoint {
                layer: "A",
                source: 0,
                target: 3,
                vertices: 3
            })
        ));
    }
}
