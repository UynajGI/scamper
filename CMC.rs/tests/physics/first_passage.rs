//! F1 weighted shortest-path validation (roadmap v2 §3 F1).
//!
//! Gates, verbatim from the roadmap and stage brief:
//!
//! 1. **Bellman-Ford parity** — production Dijkstra vs a TEST-LOCAL
//!    independent Bellman-Ford (written below over an edge list, sharing no
//!    code with production) across chains, cycles, trees, and dense ER
//!    multigraphs with parallel edges and zero-weight edges, weights
//!    uniform [0,1) and with exact zeros mixed, plus single-source sweeps
//!    from EVERY vertex on the smallest graphs. Equality policy, stated
//!    precisely: on unique-path families (chains, trees) and on dyadic
//!    two-level weight fields ({0.0, 0.5, 1.0}, where every path sum is
//!    exact) distances assert BIT-equal; on uniform [0,1) fields equal-cost
//!    ties can pick different paths in the two algorithms, whose f64 sums
//!    may differ in final ulps, so those assert within 1e-12 relative
//!    (unreachable vertices must still agree exactly).
//! 2. **Hand cases** — a 9-vertex 11-edge multigraph hand-derived on paper
//!    from the documented policies BEFORE any assertion: parallel pair
//!    (lighter wins, pinned through the predecessor edge), a zero-weight
//!    triangle (cycle), equal-cost alternatives (tie observable through the
//!    predecessor choice), a self-loop, and a disconnected component (the
//!    `f64::INFINITY` unreachable convention), plus a pure zero-weight ring
//!    and a single-vertex graph.
//! 3. **Typed rejections** — negative weight (edge named), NaN, +infinity,
//!    wrong length, source domain, empty graph; every rejection fires at
//!    validation time so a rejected call leaves the previous run's results
//!    untouched and a valid graph then passes cleanly.
//! 4. **Determinism/reuse** — repeated calls from different sources on one
//!    prepared workspace are bitwise identical to fresh-workspace calls;
//!    reuse across graphs of growing (and shrinking) sizes is correct;
//!    same-size different-content graphs alternate on one workspace
//!    (dense IDs carry no graph provenance).

use cmc_rs::{
    build_chain, first_passage_run, Bond, BondType, CsrLattice, EdgeId, FirstPassageError,
    FirstPassageOutcome, FirstPassageWorkspace, GeodesicTracking, GraphView, Incidence,
    UndirectedGraphView, VertexId, NO_PREDECESSOR,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

// ─────────────────────────────────────────────────────────────────────────
// Test-local references and fixtures
// ─────────────────────────────────────────────────────────────────────────

/// Independent Bellman-Ford over an edge list: repeated full relaxation
/// rounds to the fixpoint (early exit when a round changes nothing).
/// Undirected: both directions of every edge relax. `INFINITY + w` stays
/// `INFINITY` and `INFINITY < INFINITY` is false, so unreachable vertices
/// never relax; non-negative weights make self-loops no-ops. Shares NO code
/// with the production Dijkstra.
fn bellman_ford(vertices: usize, edges: &[(usize, usize, f64)], source: usize) -> Vec<f64> {
    let mut dist = vec![f64::INFINITY; vertices];
    dist[source] = 0.0;
    for _ in 0..vertices {
        let mut changed = false;
        for &(u, v, w) in edges {
            if dist[u] + w < dist[v] {
                dist[v] = dist[u] + w;
                changed = true;
            }
            if dist[v] + w < dist[u] {
                dist[u] = dist[v] + w;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    dist
}

/// Build the production graph from an edge list, returning the lattice and
/// the weight slice read back out of `Bond.weight` in dense edge order —
/// the exact mapping the production rustdoc documents.
fn build(edges: &[(usize, usize, f64)], vertices: usize) -> (CsrLattice, Vec<f64>) {
    let bonds = edges
        .iter()
        .map(|&(source, target, weight)| Bond::new(source, target, BondType::Generic, weight))
        .collect();
    let lattice = CsrLattice::from_edges(vertices, bonds);
    let weights: Vec<f64> = lattice.edges.iter().map(|bond| bond.weight).collect();
    (lattice, weights)
}

fn chain_edges(n: usize, seed: u64, zero_fraction: f64) -> Vec<(usize, usize, f64)> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..n - 1)
        .map(|index| (index, index + 1, draw_weight(&mut rng, zero_fraction)))
        .collect()
}

fn ring_edges(n: usize, seed: u64, zero_fraction: f64) -> Vec<(usize, usize, f64)> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..n)
        .map(|index| (index, (index + 1) % n, draw_weight(&mut rng, zero_fraction)))
        .collect()
}

/// Random recursive tree: vertex i attaches to a uniform earlier vertex.
fn tree_edges(n: usize, seed: u64, zero_fraction: f64) -> Vec<(usize, usize, f64)> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut edges = Vec::with_capacity(n - 1);
    for index in 1..n {
        let parent = rng.random_range(0..index);
        edges.push((parent, index, draw_weight(&mut rng, zero_fraction)));
    }
    edges
}

/// Dense ER multigraph: `m` uniform random pairs (self-loops excluded,
/// duplicates kept — parallel edges are the point).
fn er_multigraph_edges(
    n: usize,
    m: usize,
    seed: u64,
    zero_fraction: f64,
) -> Vec<(usize, usize, f64)> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut edges = Vec::with_capacity(m + 3);
    while edges.len() < m {
        let u = rng.random_range(0..n);
        let v = rng.random_range(0..n);
        if u != v {
            edges.push((u, v, draw_weight(&mut rng, zero_fraction)));
        }
    }
    // Guarantee parallel edges independent of the birthday coincidence.
    for _ in 0..3 {
        let (u, v, w) = edges[0];
        edges.push((u, v, w));
    }
    edges
}

fn draw_weight(rng: &mut Xoshiro256PlusPlus, zero_fraction: f64) -> f64 {
    let draw = rng.random::<f64>();
    if draw < zero_fraction {
        0.0
    } else {
        draw
    }
}

/// Dyadic two-level field {0.0, 0.5, 1.0} (every path sum exact in f64).
fn dyadic_weight(rng: &mut Xoshiro256PlusPlus) -> f64 {
    [0.0, 0.5, 1.0][rng.random_range(0..3)]
}

/// A zero-vertex view for the `EmptyGraph` rejection.
struct EmptyView;

impl GraphView for EmptyView {
    fn vertex_count(&self) -> usize {
        0
    }
}

impl UndirectedGraphView for EmptyView {
    fn edge_count(&self) -> usize {
        0
    }

    fn edge_endpoints(&self, _edge: EdgeId) -> [VertexId; 2] {
        unreachable!("no edges exist")
    }

    fn incidences(&self, _vertex: VertexId) -> impl Iterator<Item = Incidence> + '_ {
        std::iter::empty()
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1: Bellman-Ford parity
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn dijkstra_matches_independent_bellman_ford() {
    let mut cases = 0_usize;
    let mut comparisons = 0_usize;
    let mut exact_cases = 0_usize;

    let check_exact = |label: &str,
                       edges: &[(usize, usize, f64)],
                       vertices: usize,
                       sources: &[usize],
                       cases: &mut usize,
                       exact_cases: &mut usize,
                       comparisons: &mut usize| {
        let (lattice, weights) = build(edges, vertices);
        for &source in sources {
            let mut workspace = FirstPassageWorkspace::new();
            let outcome = first_passage_run(
                &lattice,
                &weights,
                source,
                GeodesicTracking::Predecessors,
                &mut workspace,
            )
            .expect("exact family run");
            let reference = bellman_ford(vertices, edges, source);
            *cases += 1;
            *exact_cases += 1;
            assert_eq!(outcome.source, source);
            for (vertex, (&produced, &expected)) in workspace
                .distances()
                .iter()
                .zip(reference.iter())
                .enumerate()
            {
                *comparisons += 1;
                assert_eq!(
                    produced, expected,
                    "{label} V={vertices} source={source} vertex={vertex}: bit-exact parity"
                );
            }
        }
    };

    // Exact families: bit-equality. Chains and trees have unique paths
    // (both algorithms accumulate the same sum in path order); dyadic
    // weights make every path sum exact even where ties pick different
    // paths.
    for n in [1_usize, 2, 3, 5, 8, 13, 21] {
        check_exact(
            "chain",
            &chain_edges(n, 0x0E_0000 + n as u64, 0.25),
            n,
            &sources_for(n),
            &mut cases,
            &mut exact_cases,
            &mut comparisons,
        );
        check_exact(
            "ring",
            &ring_edges(n, 0x0E_1000 + n as u64, 0.25),
            n,
            &sources_for(n),
            &mut cases,
            &mut exact_cases,
            &mut comparisons,
        );
        check_exact(
            "tree",
            &tree_edges(n, 0x0E_2000 + n as u64, 0.25),
            n,
            &sources_for(n),
            &mut cases,
            &mut exact_cases,
            &mut comparisons,
        );
    }
    for n in [6_usize, 10, 14] {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0E_3000 + n as u64);
        let mut edges = Vec::new();
        while edges.len() < 4 * n {
            let u = rng.random_range(0..n);
            let v = rng.random_range(0..n);
            if u != v {
                edges.push((u, v, dyadic_weight(&mut rng)));
            }
        }
        for _ in 0..3 {
            edges.push(edges[0]); // guaranteed parallel edges
        }
        check_exact(
            "er-multigraph-dyadic",
            &edges,
            n,
            &sources_for(n),
            &mut cases,
            &mut exact_cases,
            &mut comparisons,
        );
    }

    // Unit-weight builder tori through the Bond.weight mapping (integer
    // path sums, massive ties).
    for lattice in [build_chain(12, true), build_chain(21, true)] {
        let vertices = lattice.vertex_count();
        let reference_edges: Vec<(usize, usize, f64)> = lattice
            .edges
            .iter()
            .map(|bond| (bond.source, bond.target, bond.weight))
            .collect();
        let weights: Vec<f64> = lattice.edges.iter().map(|bond| bond.weight).collect();
        for source in sources_for(vertices) {
            let mut workspace = FirstPassageWorkspace::new();
            first_passage_run(
                &lattice,
                &weights,
                source,
                GeodesicTracking::Off,
                &mut workspace,
            )
            .expect("builder torus run");
            let reference = bellman_ford(vertices, &reference_edges, source);
            cases += 1;
            exact_cases += 1;
            for (vertex, (&produced, &expected)) in workspace
                .distances()
                .iter()
                .zip(reference.iter())
                .enumerate()
            {
                comparisons += 1;
                assert_eq!(
                    produced, expected,
                    "unit torus V={vertices} source={source} vertex={vertex}"
                );
            }
        }
    }

    // Relative family: uniform [0,1) weights with exact zeros mixed — ties
    // can resolve to different paths in the two algorithms, so distances
    // assert within 1e-12 relative; reachability must still agree exactly.
    for (n, m, seed) in [
        (8_usize, 24_usize, 0x11_u64),
        (16, 64, 0x12),
        (24, 120, 0x13),
        (32, 256, 0x14),
        (48, 512, 0x15),
    ] {
        for zero_fraction in [0.0_f64, 0.25, 0.6] {
            let edges =
                er_multigraph_edges(n, m, seed + (zero_fraction * 100.0) as u64, zero_fraction);
            let (lattice, weights) = build(&edges, n);
            for source in sources_for(n) {
                let mut workspace = FirstPassageWorkspace::new();
                first_passage_run(
                    &lattice,
                    &weights,
                    source,
                    GeodesicTracking::Off,
                    &mut workspace,
                )
                .expect("relative family run");
                let reference = bellman_ford(n, &edges, source);
                cases += 1;
                for (vertex, (&produced, &expected)) in workspace
                    .distances()
                    .iter()
                    .zip(reference.iter())
                    .enumerate()
                {
                    comparisons += 1;
                    assert_eq!(
                        produced.is_finite(),
                        expected.is_finite(),
                        "reachability parity: multigraph n={n} m={m} zero_fraction={zero_fraction} \
                         source={source} vertex={vertex}"
                    );
                    if expected.is_finite() {
                        let tolerance = 1e-12 * expected.abs().max(1.0);
                        assert!(
                            (produced - expected).abs() <= tolerance,
                            "relative parity: multigraph n={n} m={m} source={source} \
                             vertex={vertex}: {produced} vs {expected}"
                        );
                    }
                }
            }
        }
    }

    eprintln!(
        "F1 BF parity: {cases} (graph, source) cases ({exact_cases} bit-exact), {comparisons} \
         vertex-distance comparisons, 0 mismatches"
    );
}

/// Every vertex as source on the smallest graphs, a spread otherwise.
fn sources_for(vertices: usize) -> Vec<usize> {
    if vertices <= 16 {
        (0..vertices).collect()
    } else {
        vec![0, vertices / 2, vertices - 1]
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 2: hand cases
// ─────────────────────────────────────────────────────────────────────────

/// The hand-worked multigraph: 9 vertices (0..8), 11 physical edges:
/// e0(0,1) w=5.0, e1(0,1) w=2.0 [parallel pair, lighter wins],
/// e2(1,2) w=0.0, e3(2,3) w=0.0, e4(3,1) w=0.0 [zero-weight triangle],
/// e5(0,4) w=1.0, e6(0,5) w=1.0, e7(4,6) w=1.0, e8(5,6) w=1.0
/// [equal-cost alternatives into 6], e9(6,6) w=0.0 [self-loop],
/// e10(7,8) w=1.0 [disconnected component]; vertices reached only from 0.
fn hand_graph() -> (CsrLattice, Vec<f64>) {
    build(
        &[
            (0, 1, 5.0),
            (0, 1, 2.0),
            (1, 2, 0.0),
            (2, 3, 0.0),
            (3, 1, 0.0),
            (0, 4, 1.0),
            (0, 5, 1.0),
            (4, 6, 1.0),
            (5, 6, 1.0),
            (6, 6, 0.0),
            (7, 8, 1.0),
        ],
        9,
    )
}

#[test]
fn hand_multigraph_distances_and_predecessors_match_derivation() {
    let (lattice, weights) = hand_graph();
    let mut workspace = FirstPassageWorkspace::new();
    let outcome = first_passage_run(
        &lattice,
        &weights,
        0,
        GeodesicTracking::Predecessors,
        &mut workspace,
    )
    .expect("hand graph run");

    // Settle trace, derived on paper from the documented tie policy
    // ((distance total_cmp, vertex id ascending), strict-improvement-only
    // predecessor updates) BEFORE running production:
    //   settle 0 (0.0): relax e0 -> dist1=5.0; e1 -> dist1=2.0 (lighter
    //     parallel bond wins); e5 -> dist4=1.0; e6 -> dist5=1.0.
    //   settle 4 (1.0, id 4 < 5 at equal distance): e7 -> dist6=2.0, pred
    //     (4, e7).
    //   settle 5 (1.0): e8 candidate 2.0 == 2.0, NOT strict -> tie keeps
    //     the predecessor set by the lower-id settler (4, e7).
    //   settle 1 (2.0, id 1 < 2, 3, 6 at equal distance): e2 -> dist2=2.0
    //     pred (1, e2); e4 -> dist3=2.0 pred (1, e4).
    //   settle 2 (2.0): e3 candidate 2.0 == 2.0, keeps (1, e4).
    //   settle 3 (2.0), then 6 (2.0, e9 self-loop 2.0 == 2.0, ignored);
    //   the 5.0 entry for 1 pops stale.
    let expected_distances = [
        0.0,
        2.0,
        2.0,
        2.0,
        1.0,
        1.0,
        2.0,
        f64::INFINITY,
        f64::INFINITY,
    ];
    let expected_pred_vertices = [
        NO_PREDECESSOR,
        0,
        1,
        1,
        0,
        0,
        4,
        NO_PREDECESSOR,
        NO_PREDECESSOR,
    ];
    let expected_pred_edges = [
        NO_PREDECESSOR,
        1, // the lighter parallel bond e1, not e0
        2,
        4,
        5,
        6,
        7, // equal-cost tie: the lower-id settler (4) kept the predecessor
        NO_PREDECESSOR,
        NO_PREDECESSOR,
    ];
    assert_eq!(workspace.distances(), expected_distances);
    assert_eq!(workspace.predecessor_vertices(), expected_pred_vertices);
    assert_eq!(workspace.predecessor_edges(), expected_pred_edges);
    assert_eq!(
        outcome,
        FirstPassageOutcome {
            source: 0,
            reached: 7
        }
    );
    for vertex in 0..7 {
        assert!(workspace.is_reachable(vertex), "vertex {vertex} reached");
    }
    for vertex in [7_usize, 8] {
        assert!(
            !workspace.is_reachable(vertex),
            "disconnected vertex {vertex} unreachable (INFINITY sentinel)"
        );
    }

    // Reverse direction: from inside the disconnected component everything
    // else is the sentinel.
    let outcome = first_passage_run(
        &lattice,
        &weights,
        7,
        GeodesicTracking::Predecessors,
        &mut workspace,
    )
    .expect("hand graph run from 7");
    let mut expected = vec![f64::INFINITY; 9];
    expected[7] = 0.0;
    expected[8] = 1.0;
    assert_eq!(workspace.distances(), expected);
    assert_eq!(outcome.reached, 2);
}

#[test]
fn zero_weight_ring_terminates_with_all_zero_distances() {
    let (lattice, weights) = build(&[(0, 1, 0.0), (1, 2, 0.0), (2, 3, 0.0), (3, 0, 0.0)], 4);
    let mut workspace = FirstPassageWorkspace::new();
    let outcome = first_passage_run(&lattice, &weights, 0, GeodesicTracking::Off, &mut workspace)
        .expect("zero ring run");
    assert_eq!(workspace.distances(), [0.0, 0.0, 0.0, 0.0]);
    assert_eq!(outcome.reached, 4);
}

#[test]
fn single_vertex_graph_has_trivial_distances() {
    let (lattice, weights) = build(&[], 1);
    let mut workspace = FirstPassageWorkspace::new();
    let outcome = first_passage_run(&lattice, &weights, 0, GeodesicTracking::Off, &mut workspace)
        .expect("single vertex run");
    assert_eq!(workspace.distances(), [0.0]);
    assert_eq!(outcome.reached, 1);
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3: typed rejections
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn typed_rejections_fire_at_validation_and_leave_results_untouched() {
    let (lattice, weights) = hand_graph();
    let mut workspace = FirstPassageWorkspace::new();
    let outcome = first_passage_run(
        &lattice,
        &weights,
        0,
        GeodesicTracking::Predecessors,
        &mut workspace,
    )
    .expect("valid baseline run");
    let baseline = workspace.distances().to_vec();
    let baseline_pred = workspace.predecessor_vertices().to_vec();
    assert_eq!(outcome.reached, 7);

    // Negative weight names the edge.
    let mut bad = weights.clone();
    bad[4] = -0.25;
    let error = first_passage_run(&lattice, &bad, 0, GeodesicTracking::Off, &mut workspace)
        .expect_err("negative weight rejected");
    assert!(matches!(
        error,
        FirstPassageError::NegativeWeight { edge: 4, value } if value == -0.25
    ));
    assert!(
        error.to_string().contains("edge 4"),
        "the rejection must name the edge: {error}"
    );

    // NaN and +infinity weights are non-finite rejections (separate edges).
    for (edge, value) in [
        (7_usize, f64::NAN),
        (0, f64::INFINITY),
        (2, f64::NEG_INFINITY),
    ] {
        let mut bad = weights.clone();
        bad[edge] = value;
        let error = first_passage_run(&lattice, &bad, 0, GeodesicTracking::Off, &mut workspace)
            .expect_err("non-finite weight rejected");
        assert!(
            matches!(error, FirstPassageError::WeightNotFinite { edge: e, .. } if e == edge),
            "edge {edge} with {value}: {error}"
        );
    }

    // Wrong weight-field length.
    let error = first_passage_run(
        &lattice,
        &weights[..10],
        0,
        GeodesicTracking::Off,
        &mut workspace,
    )
    .expect_err("wrong length rejected");
    assert!(matches!(
        error,
        FirstPassageError::WeightLengthMismatch {
            expected: 11,
            actual: 10
        }
    ));

    // Source outside the dense range.
    let error = first_passage_run(&lattice, &weights, 9, GeodesicTracking::Off, &mut workspace)
        .expect_err("source domain rejected");
    assert!(matches!(
        error,
        FirstPassageError::SourceDomain {
            source: 9,
            vertices: 9
        }
    ));

    // Zero-vertex graph.
    let error = first_passage_run(&EmptyView, &[], 0, GeodesicTracking::Off, &mut workspace)
        .expect_err("empty graph rejected");
    assert!(matches!(error, FirstPassageError::EmptyGraph));

    // Every rejection above left the valid run's results untouched.
    assert_eq!(workspace.distances(), baseline);
    assert_eq!(workspace.predecessor_vertices(), baseline_pred);

    // And a valid graph then passes cleanly on the same workspace.
    let outcome = first_passage_run(
        &lattice,
        &weights,
        0,
        GeodesicTracking::Predecessors,
        &mut workspace,
    )
    .expect("valid rerun after rejections");
    assert_eq!(outcome.reached, 7);
    assert_eq!(workspace.distances(), baseline);
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 4: determinism and workspace reuse
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn reused_workspace_matches_fresh_runs_bitwise_across_sources() {
    let edges = er_multigraph_edges(12, 36, 0x4001, 0.3);
    let (lattice, weights) = build(&edges, 12);
    let mut shared = FirstPassageWorkspace::new();
    let mut comparisons = 0_usize;

    // Fresh reference per source.
    let mut fresh = Vec::new();
    for source in 0..12 {
        let mut workspace = FirstPassageWorkspace::new();
        first_passage_run(
            &lattice,
            &weights,
            source,
            GeodesicTracking::Predecessors,
            &mut workspace,
        )
        .expect("fresh run");
        fresh.push((
            workspace.distances().to_vec(),
            workspace.predecessor_vertices().to_vec(),
            workspace.predecessor_edges().to_vec(),
        ));
    }

    // Shared workspace: cycle through all sources forward, then backward,
    // interleaving tracking-off runs; every result must reproduce the fresh
    // run bitwise.
    for source in (0..12).chain((0..12).rev()) {
        first_passage_run(
            &lattice,
            &weights,
            source,
            GeodesicTracking::Predecessors,
            &mut shared,
        )
        .expect("shared run");
        assert_eq!(
            shared.distances(),
            fresh[source].0,
            "distances source {source}"
        );
        assert_eq!(
            shared.predecessor_vertices(),
            fresh[source].1,
            "predecessor vertices source {source}"
        );
        assert_eq!(
            shared.predecessor_edges(),
            fresh[source].2,
            "predecessor edges source {source}"
        );
        comparisons += 1;
        // A tracking-off run in between must not disturb the buffers the
        // next tracking-on run rewrites; its own distances stay valid.
        first_passage_run(
            &lattice,
            &weights,
            11 - source,
            GeodesicTracking::Off,
            &mut shared,
        )
        .expect("shared tracking-off run");
        assert_eq!(shared.distances(), fresh[11 - source].0);
        comparisons += 1;
    }
    eprintln!("F1 reuse determinism: {comparisons} shared-vs-fresh bitwise comparisons, 0 diffs");
}

#[test]
fn workspace_reuse_across_growing_and_shrinking_graphs() {
    let mut workspace = FirstPassageWorkspace::new();

    let chain = chain_edges(8, 0x5001, 0.2);
    let (chain_lattice, chain_weights) = build(&chain, 8);
    for source in 0..8 {
        first_passage_run(
            &chain_lattice,
            &chain_weights,
            source,
            GeodesicTracking::Off,
            &mut workspace,
        )
        .expect("chain run");
        assert_eq!(
            workspace.distances(),
            bellman_ford(8, &chain, source).as_slice()
        );
    }
    assert_eq!(workspace.prepared_sizes(), (8, 7));

    let tree = tree_edges(24, 0x5002, 0.2);
    let (tree_lattice, tree_weights) = build(&tree, 24);
    for source in [0_usize, 11, 23] {
        first_passage_run(
            &tree_lattice,
            &tree_weights,
            source,
            GeodesicTracking::Off,
            &mut workspace,
        )
        .expect("tree run");
        assert_eq!(
            workspace.distances(),
            bellman_ford(24, &tree, source).as_slice()
        );
    }
    assert_eq!(workspace.prepared_sizes(), (24, 23));

    let dense = er_multigraph_edges(64, 320, 0x5003, 0.2);
    let (dense_lattice, dense_weights) = build(&dense, 64);
    for source in [0_usize, 31, 63] {
        first_passage_run(
            &dense_lattice,
            &dense_weights,
            source,
            GeodesicTracking::Predecessors,
            &mut workspace,
        )
        .expect("dense run");
        assert_eq!(
            workspace.distances(),
            bellman_ford(64, &dense, source).as_slice()
        );
    }
    assert_eq!(workspace.prepared_sizes(), (64, 323));

    // Shrinking back re-uses retained capacity and stays correct.
    for source in 0..8 {
        first_passage_run(
            &chain_lattice,
            &chain_weights,
            source,
            GeodesicTracking::Off,
            &mut workspace,
        )
        .expect("chain rerun");
        assert_eq!(
            workspace.distances(),
            bellman_ford(8, &chain, source).as_slice()
        );
    }
    assert_eq!(workspace.prepared_sizes(), (8, 7));
}

#[test]
fn same_size_different_graphs_alternate_on_one_workspace() {
    // Dense IDs carry no graph provenance: one workspace serves any graph
    // with the same (V, E), positionally reinterpreted.
    let edges_a = er_multigraph_edges(16, 48, 0x6001, 0.3);
    let edges_b = er_multigraph_edges(16, 48, 0x6002, 0.3);
    let (lattice_a, weights_a) = build(&edges_a, 16);
    let (lattice_b, weights_b) = build(&edges_b, 16);
    let mut workspace = FirstPassageWorkspace::new();

    first_passage_run(
        &lattice_a,
        &weights_a,
        0,
        GeodesicTracking::Predecessors,
        &mut workspace,
    )
    .expect("graph A first run");
    let first = (
        workspace.distances().to_vec(),
        workspace.predecessor_edges().to_vec(),
    );
    assert_eq!(
        workspace.distances(),
        bellman_ford(16, &edges_a, 0).as_slice()
    );

    first_passage_run(
        &lattice_b,
        &weights_b,
        0,
        GeodesicTracking::Predecessors,
        &mut workspace,
    )
    .expect("graph B run");
    assert_eq!(
        workspace.distances(),
        bellman_ford(16, &edges_b, 0).as_slice()
    );

    first_passage_run(
        &lattice_a,
        &weights_a,
        0,
        GeodesicTracking::Predecessors,
        &mut workspace,
    )
    .expect("graph A rerun");
    assert_eq!(workspace.distances(), first.0, "graph A reproduces bitwise");
    assert_eq!(workspace.predecessor_edges(), first.1);
}
