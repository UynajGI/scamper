//! G3 invasion percolation validation (roadmap §13).
//!
//! Gates, verbatim from the roadmap:
//!
//! 1. **fixed-weight hand-worked graphs** — a 7-vertex, 9-edge multigraph
//!    (parallel pair + self-loop) with hand-assigned thresholds, the complete
//!    invasion sequence, the frontier evolution, and the invasion tree
//!    verified BY HAND against production, in both site and bond mode;
//! 2. **binary-heap implementation versus scanning the complete frontier
//!    each step** — production against an independent test-side reference
//!    that rebuilds the whole frontier from scratch before every acceptance,
//!    across hundreds of random realizations (continuous AND tie-dense
//!    quantized thresholds, single and multiple seeds, all stop policies,
//!    trapping and non-trapping, multigraphs and two-component graphs);
//! 3. **tree uniqueness** — the invasion graph is a forest: exactly
//!    `invaded - seeds` accepted edges, no cycles (independent union-find
//!    replay), one parent per invaded non-seed vertex, one tree per seed;
//! 4. **non-trapping bond invasion agrees with the applicable Prim/MST
//!    prefix relation when thresholds are unique** — stated precisely:
//!    with a single seed, the accepted-edge sequence is exactly Prim's
//!    growth order from the seed and the terminal tree is the unique MST of
//!    the seed's connected component (verified against an independent
//!    test-side Prim AND an independent Kruskal on random unique-weight
//!    graphs); with `k` seeds it is Prim's growth order from a virtual
//!    super-root joined to every seed by `-inf` edges, and the terminal
//!    forest is the min-weight spanning forest with `k` trees, one seed
//!    each; a `Steps(n)` run is the exact `n`-edge Prim prefix.
//!
//! Plus the explicit-policy gates the roadmap asks to "make explicit": the
//! tie policy (equal thresholds resolve by the lowest dense ID — vertex ID
//! in site mode, edge ID in bond mode; `-0.0 < +0.0` under `total_cmp` are
//! distinct keys, not a tie), the stop policies (`GrowAll`,
//! `FirstSpanning{from,to}` at the BFS-verified spanning step, `Steps(n)`
//! exact), the Wilkinson-Willemsen trapping variant (hand-worked, including
//! the invaded-exit and retired-bond-entry cases), the adapter's bitwise
//! RNG draw contract and C0 history compression, the driven G1 substrate
//! (trajectory steps are the invasion times), typed rejections, and
//! workspace reuse across runs.

use std::cmp::Ordering;

use cmc_rs::{
    build_chain, build_square, invasion_run, Bond, BondType, BorrowedUndirectedCsr, CsrLattice,
    GraphView, InvasionError, InvasionMode, InvasionOutcome, InvasionPercolationMC, InvasionStop,
    InvasionStopReason, InvasionTrapping, InvasionWorkspace,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

// ─────────────────────────────────────────────────────────────────────────
// Fixtures
// ─────────────────────────────────────────────────────────────────────────

fn bond(source: usize, target: usize) -> Bond {
    Bond::new(source, target, BondType::Generic, 1.0)
}

/// The hand-worked graph: 7 vertices, 9 physical edges —
/// e0(0,1), e1(0,2), e2(1,3), e3(1,3) [parallel to e2], e4(2,3),
/// e5(3,3) [self-loop], e6(3,4), e7(4,5), e8(4,6).
fn hand_graph() -> CsrLattice {
    CsrLattice::try_from_edges(
        7,
        vec![
            bond(0, 1),
            bond(0, 2),
            bond(1, 3),
            bond(1, 3),
            bond(2, 3),
            bond(3, 3),
            bond(3, 4),
            bond(4, 5),
            bond(4, 6),
        ],
    )
    .expect("hand graph is valid")
}

/// Hand-assigned bond thresholds (all nine distinct; the self-loop carries
/// the largest and must never be accepted).
const HAND_BOND_THRESHOLDS: [f64; 9] = [0.50, 0.20, 0.70, 0.35, 0.60, 0.90, 0.10, 0.45, 0.80];

/// Hand-assigned site thresholds (the seed's entry is never read).
const HAND_SITE_THRESHOLDS: [f64; 7] = [0.00, 0.55, 0.30, 0.15, 0.65, 0.25, 0.85];

/// Star: e0(0,1), e1(0,2), e2(0,3) — the tie-policy fixture.
fn star() -> CsrLattice {
    CsrLattice::try_from_edges(4, vec![bond(0, 1), bond(0, 2), bond(0, 3)]).expect("star is valid")
}

/// Chain 0-1-2-3-4 — the trapping fixture.
fn chain5() -> CsrLattice {
    build_chain(5, false)
}

/// Two components (triangle + path) plus an isolated vertex.
fn two_components() -> CsrLattice {
    CsrLattice::try_from_edges(
        7,
        vec![bond(0, 1), bond(1, 2), bond(2, 0), bond(3, 4), bond(4, 5)],
    )
    .expect("two-component graph is valid")
}

/// Random multigraph with self-loops, parallel edges, and a guaranteed
/// connecting chain.
fn random_multigraph(vertices: usize, seed: u64) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut edges = Vec::new();
    for left in 0..vertices {
        // `right` from `left` admits self-loops; two windows admit parallels.
        for right in left..vertices {
            let roll = rng.random::<f64>();
            if roll < 0.22 {
                edges.push(bond(left, right));
            }
            if (0.70..0.78).contains(&roll) {
                edges.push(bond(left, right));
            }
        }
    }
    for vertex in 1..vertices {
        edges.push(bond(vertex - 1, vertex));
    }
    CsrLattice::try_from_edges(vertices, edges).expect("random multigraph is valid")
}

const NON_TRAPPING: InvasionTrapping = InvasionTrapping::NonTrapping;

fn trapping(exit: &[usize]) -> InvasionTrapping {
    InvasionTrapping::Trapping {
        exit: exit.to_vec(),
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Independent reference: scan the complete frontier every step
// ─────────────────────────────────────────────────────────────────────────

/// One frontier entry of the reference, exactly as a scan sees it.
#[derive(Clone, Copy, Debug, PartialEq)]
struct NaiveEntry {
    threshold: f64,
    /// Tie key: vertex index (site mode) or edge index (bond mode).
    id: usize,
    /// The un-invaded vertex acceptance would add.
    vertex: usize,
    /// The edge acceptance would record.
    edge: usize,
}

/// The reference's stop policy (mirrors the production API values).
#[derive(Clone, Copy, Debug)]
enum NaiveStop {
    GrowAll,
    FirstSpanning { from: usize, to: usize },
    Steps(usize),
}

impl NaiveStop {
    fn as_production(self) -> InvasionStop {
        match self {
            Self::GrowAll => InvasionStop::GrowAll,
            Self::FirstSpanning { from, to } => InvasionStop::FirstSpanning { from, to },
            Self::Steps(n) => InvasionStop::Steps(n),
        }
    }
}

/// Result of one reference realization: the full per-step trace.
struct NaiveTrace {
    order: Vec<usize>,
    accepted_edges: Vec<usize>,
    accepted_thresholds: Vec<f64>,
    /// `usize::MAX` un-invaded, the seed itself for seeds, the invader else.
    parent: Vec<usize>,
    /// Live frontier (sorted by the documented key) after planting (index 0)
    /// and after every accepted growth step.
    frontiers: Vec<Vec<NaiveEntry>>,
    /// Sites retired as trapped, in retirement order (trapping only).
    trapped: Vec<usize>,
    stop_reason: InvasionStopReason,
    stop_step: usize,
}

/// The independent reference: before every acceptance it REBUILDS the
/// complete frontier by scanning all edges (bond mode: every edge with
/// exactly one invaded endpoint, self-loops excluded) or all vertices
/// (site mode: every un-invaded vertex adjacent to the invaded set), sorts
/// it by the documented key `(total_cmp(threshold), dense id)`, and accepts
/// the first entry — the roadmap's "scanning the complete frontier each
/// step" alternative to the production binary heap.
///
/// Site-mode tree edges replicate the documented push rule: a site is
/// queued by the FIRST invader to reach it (seeds in caller order, then one
/// vertex per accepted step) and carries that invader's first incidence.
fn naive_invasion(
    graph: &CsrLattice,
    mode: InvasionMode,
    thresholds: &[f64],
    seeds: &[usize],
    stop: NaiveStop,
    exit: Option<&[usize]>,
) -> NaiveTrace {
    let vertices = graph.vertex_count();
    let endpoints = |edge: usize| [graph.edges[edge].source, graph.edges[edge].target];

    let mut parent = vec![usize::MAX; vertices];
    let mut carrying: Vec<Option<usize>> = vec![None; vertices]; // site mode
    let mut retired = vec![false; vertices];
    let mut order = Vec::new();
    let mut accepted_edges = Vec::new();
    let mut accepted_thresholds = Vec::new();
    let mut frontiers = Vec::new();
    let mut trapped = Vec::new();
    let invaded = |parent: &[usize], vertex: usize| parent[vertex] != usize::MAX;

    // Site-mode provenance: record the carrying edge exactly when the first
    // invader reaches the site (invasion events in production order).
    // (`CsrLattice::incidences` is the inherent (neighbor, edge_id) scan.)
    fn note_invasion_event(
        graph: &CsrLattice,
        parent: &[usize],
        carrying: &mut [Option<usize>],
        vertex: usize,
    ) {
        for (neighbor, edge) in graph.incidences(vertex) {
            if parent[neighbor] == usize::MAX && carrying[neighbor].is_none() {
                carrying[neighbor] = Some(edge);
            }
        }
    }

    for &seed in seeds {
        parent[seed] = seed;
        order.push(seed);
        if matches!(mode, InvasionMode::Site) {
            note_invasion_event(graph, &parent, &mut carrying, seed);
        }
    }
    if let NaiveStop::FirstSpanning { from, to } = stop {
        if invaded(&parent, from) && invaded(&parent, to) {
            return NaiveTrace {
                order,
                accepted_edges,
                accepted_thresholds,
                parent,
                frontiers,
                trapped,
                stop_reason: InvasionStopReason::SpanningStep,
                stop_step: 0,
            };
        }
    }
    if let NaiveStop::Steps(0) = stop {
        return NaiveTrace {
            order,
            accepted_edges,
            accepted_thresholds,
            parent,
            frontiers,
            trapped,
            stop_reason: InvasionStopReason::StepBudget,
            stop_step: 0,
        };
    }

    // Wilkinson-Willemsen escape test: BFS over the un-invaded (defending)
    // phase; an un-invaded exit vertex inside the region provides escape.
    fn escapes(graph: &CsrLattice, parent: &[usize], exit: Option<&[usize]>, start: usize) -> bool {
        let is_exit = |vertex: usize| exit.is_some_and(|exit| exit.contains(&vertex));
        if is_exit(start) {
            return true;
        }
        let mut seen = vec![false; graph.vertex_count()];
        seen[start] = true;
        let mut stack = vec![start];
        while let Some(vertex) = stack.pop() {
            for (neighbor, _) in graph.incidences(vertex) {
                if parent[neighbor] != usize::MAX || seen[neighbor] {
                    continue;
                }
                if is_exit(neighbor) {
                    return true;
                }
                seen[neighbor] = true;
                stack.push(neighbor);
            }
        }
        false
    }

    loop {
        // Rebuild the complete frontier from scratch.
        let mut frontier = Vec::new();
        match mode {
            InvasionMode::Bond => {
                for (edge, &threshold) in thresholds.iter().enumerate() {
                    let [left, right] = endpoints(edge);
                    if left == right {
                        continue; // self-loops never enter the frontier
                    }
                    match (invaded(&parent, left), invaded(&parent, right)) {
                        (true, false) => frontier.push(NaiveEntry {
                            threshold,
                            id: edge,
                            vertex: right,
                            edge,
                        }),
                        (false, true) => frontier.push(NaiveEntry {
                            threshold,
                            id: edge,
                            vertex: left,
                            edge,
                        }),
                        _ => {}
                    }
                }
            }
            InvasionMode::Site => {
                for vertex in 0..vertices {
                    // `carrying` keeps its entry after acceptance; a fresh
                    // scan must exclude already-invaded sites.
                    if !invaded(&parent, vertex) {
                        if let Some(edge) = carrying[vertex] {
                            frontier.push(NaiveEntry {
                                threshold: thresholds[vertex],
                                id: vertex,
                                vertex,
                                edge,
                            });
                        }
                    }
                }
            }
        }
        frontier.sort_by(|a, b| {
            a.threshold
                .total_cmp(&b.threshold)
                .then_with(|| a.id.cmp(&b.id))
        });
        frontiers.push(frontier.clone());

        // Walk the sorted frontier in key order (mirrors the production pop
        // loop: retired sites consume their entries, trapped sites retire
        // permanently, the first surviving entry is accepted).
        let mut accepted = None;
        for entry in &frontier {
            if exit.is_some() {
                if retired[entry.vertex] {
                    continue;
                }
                if !escapes(graph, &parent, exit, entry.vertex) {
                    retired[entry.vertex] = true;
                    trapped.push(entry.vertex);
                    continue;
                }
            }
            accepted = Some(*entry);
            break;
        }

        let Some(entry) = accepted else {
            let stop_step = accepted_edges.len();
            return NaiveTrace {
                order,
                accepted_edges,
                accepted_thresholds,
                parent,
                frontiers,
                trapped,
                stop_reason: InvasionStopReason::FrontierExhausted,
                stop_step,
            };
        };

        let [left, right] = endpoints(entry.edge);
        let invader = if left == entry.vertex { right } else { left };
        parent[entry.vertex] = invader;
        order.push(entry.vertex);
        accepted_edges.push(entry.edge);
        accepted_thresholds.push(entry.threshold);
        if matches!(mode, InvasionMode::Site) {
            note_invasion_event(graph, &parent, &mut carrying, entry.vertex);
        }

        let step = accepted_edges.len();
        match stop {
            NaiveStop::GrowAll => {}
            NaiveStop::Steps(budget) => {
                if budget == step {
                    return NaiveTrace {
                        order,
                        accepted_edges,
                        accepted_thresholds,
                        parent,
                        frontiers,
                        trapped,
                        stop_reason: InvasionStopReason::StepBudget,
                        stop_step: step,
                    };
                }
            }
            NaiveStop::FirstSpanning { from, to } => {
                if invaded(&parent, from) && invaded(&parent, to) {
                    return NaiveTrace {
                        order,
                        accepted_edges,
                        accepted_thresholds,
                        parent,
                        frontiers,
                        trapped,
                        stop_reason: InvasionStopReason::SpanningStep,
                        stop_step: step,
                    };
                }
            }
        }
    }
}

/// Run production and the scan reference side by side and assert complete
/// equality of sequence, tree, frontier evolution, stop behavior, and
/// outcome scalars; then assert the forest invariants. Returns the
/// reference trace for further gate-specific assertions.
fn assert_production_matches_naive(
    context: &str,
    graph: &CsrLattice,
    mode: InvasionMode,
    thresholds: &[f64],
    seeds: &[usize],
    stop: NaiveStop,
    exit: Option<&[usize]>,
) -> NaiveTrace {
    let reference = naive_invasion(graph, mode, thresholds, seeds, stop, exit);
    let trapping = match exit {
        None => NON_TRAPPING,
        Some(exit) => trapping(exit),
    };
    let mut workspace = InvasionWorkspace::new();
    let outcome = invasion_run(
        graph,
        mode,
        thresholds,
        seeds,
        stop.as_production(),
        &trapping,
        &mut workspace,
    )
    .unwrap_or_else(|error| panic!("{context}: production rejected a valid run: {error}"));

    assert_eq!(
        workspace.invaded_order(),
        reference.order.as_slice(),
        "{context}: invasion order"
    );
    assert_eq!(
        workspace.accepted_edges(),
        reference.accepted_edges.as_slice(),
        "{context}: accepted edges"
    );
    assert_eq!(
        workspace.accepted_thresholds(),
        reference.accepted_thresholds.as_slice(),
        "{context}: accepted thresholds"
    );
    assert_eq!(
        workspace.parent(),
        reference.parent.as_slice(),
        "{context}: parent tree"
    );
    assert_eq!(
        outcome.stop_reason, reference.stop_reason,
        "{context}: stop reason"
    );
    assert_eq!(
        outcome.stop_step, reference.stop_step,
        "{context}: stop step"
    );
    assert_eq!(
        outcome.accepted_total,
        reference.accepted_edges.len(),
        "{context}: accepted total"
    );
    assert_eq!(
        outcome.invaded_vertices,
        reference.order.len(),
        "{context}: invaded vertices"
    );
    assert_eq!(
        outcome.trapped_sites,
        reference.trapped.len(),
        "{context}: trapped sites"
    );
    assert_eq!(outcome.seed_count, seeds.len(), "{context}: seed count");
    assert_eq!(
        outcome.max_accepted_threshold,
        reference
            .accepted_thresholds
            .iter()
            .copied()
            .reduce(f64::max)
            .unwrap_or(0.0),
        "{context}: max accepted threshold"
    );
    assert_forest_invariants(context, graph, &outcome, workspace.parent(), seeds);
    reference
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3: tree uniqueness (forest invariants)
// ─────────────────────────────────────────────────────────────────────────

/// The invasion graph is a forest rooted at the seeds: exactly
/// `invaded - seeds` accepted edges, no cycle (independent union-find
/// replay of the accepted edges), one parent per invaded non-seed vertex,
/// and one tree per seed (trees never merge).
fn assert_forest_invariants(
    context: &str,
    graph: &CsrLattice,
    outcome: &InvasionOutcome,
    parent: &[usize],
    seeds: &[usize],
) {
    assert_eq!(
        outcome.accepted_total + seeds.len(),
        outcome.invaded_vertices,
        "{context}: edge count must equal invaded minus seeds"
    );
    for &seed in seeds {
        assert_eq!(parent[seed], seed, "{context}: a seed is its own parent");
    }
    let mut classes: Vec<usize> = (0..graph.vertex_count()).collect();
    let root = |classes: &mut Vec<usize>, mut node: usize| {
        while classes[node] != node {
            node = classes[node];
        }
        node
    };
    let mut tree_edges = 0_usize;
    for vertex in 0..graph.vertex_count() {
        if parent[vertex] == usize::MAX || parent[vertex] == vertex {
            continue;
        }
        let invader = parent[vertex];
        assert_ne!(invader, vertex, "{context}: no self parent");
        assert!(
            parent[invader] != usize::MAX,
            "{context}: the parent must already be invaded"
        );
        let (left, right) = (root(&mut classes, vertex), root(&mut classes, invader));
        assert_ne!(
            left, right,
            "{context}: accepted edges must not close a cycle"
        );
        classes[left] = right;
        tree_edges += 1;
    }
    assert_eq!(
        tree_edges, outcome.accepted_total,
        "{context}: parent edges"
    );

    // One tree per seed: BFS from each seed over parent edges partitions the
    // invaded vertices into exactly `seeds.len()` trees.
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); graph.vertex_count()];
    for vertex in 0..graph.vertex_count() {
        if parent[vertex] != usize::MAX && parent[vertex] != vertex {
            adjacency[vertex].push(parent[vertex]);
            adjacency[parent[vertex]].push(vertex);
        }
    }
    let mut visited = vec![false; graph.vertex_count()];
    let mut trees = 0_usize;
    let mut invaded_seen = 0_usize;
    for &seed in seeds {
        assert!(!visited[seed], "{context}: each seed starts its own tree");
        visited[seed] = true;
        invaded_seen += 1;
        let mut stack = vec![seed];
        while let Some(vertex) = stack.pop() {
            for &neighbor in &adjacency[vertex] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    invaded_seen += 1;
                    stack.push(neighbor);
                }
            }
        }
        trees += 1;
    }
    assert_eq!(trees, seeds.len(), "{context}: one tree per seed");
    assert_eq!(
        invaded_seen, outcome.invaded_vertices,
        "{context}: the forest spans exactly the invaded vertices"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1: fixed-weight hand-worked graphs
// ─────────────────────────────────────────────────────────────────────────

/// Hand-worked BOND trace on the multigraph, seed {0}, GrowAll:
///
/// - plant 0 → frontier {(e0,.50)→1, (e1,.20)→2};
/// - step 1 accepts e1 (.20): 2 joins; frontier {(e0,.50)→1, (e4,.60)→3};
/// - step 2 accepts e0 (.50): 1 joins; now THREE entries target 3 — the
///   parallel e3 (.35) beats e4 (.60) and e2 (.70);
/// - step 3 accepts e3 (.35): 3 joins via the PARALLEL edge; e2/e4 go stale;
/// - step 4 accepts e6 (.10): 4 joins; the self-loop e5 (.90) is never even
///   a candidate;
/// - step 5 accepts e7 (.45): 5 joins; step 6 accepts e8 (.80): 6 joins;
/// - stale e4, e2 are popped and discarded; frontier exhausts.
#[test]
fn hand_worked_bond_invasion_matches_manual_trace() {
    let graph = hand_graph();
    let reference = assert_production_matches_naive(
        "hand bond",
        &graph,
        InvasionMode::Bond,
        &HAND_BOND_THRESHOLDS,
        &[0],
        NaiveStop::GrowAll,
        None,
    );

    let mut workspace = InvasionWorkspace::new();
    let outcome = invasion_run(
        &graph,
        InvasionMode::Bond,
        &HAND_BOND_THRESHOLDS,
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("hand bond run is valid");

    assert_eq!(workspace.invaded_order(), &[0, 2, 1, 3, 4, 5, 6]);
    assert_eq!(workspace.accepted_edges(), &[1, 0, 3, 6, 7, 8]);
    assert_eq!(
        workspace.accepted_thresholds(),
        &[0.20, 0.50, 0.35, 0.10, 0.45, 0.80]
    );
    // The tree, by hand: 2 and 1 hang off 0, 3 off 1 (through the parallel
    // e3), 4 off 3, 5 and 6 off 4.
    let expected_parent = [0, 0, 0, 1, 3, 4, 4];
    assert_eq!(workspace.parent(), &expected_parent[..]);
    // The self-loop e5 and both stale entries (e2, e4) are never accepted.
    for forbidden in [2_usize, 4, 5] {
        assert!(
            !workspace.accepted_edges().contains(&forbidden),
            "edge {forbidden} must not be accepted"
        );
    }
    assert_eq!(
        outcome,
        InvasionOutcome {
            seed_count: 1,
            accepted_total: 6,
            invaded_vertices: 7,
            trapped_sites: 0,
            max_accepted_threshold: 0.80,
            stop_reason: InvasionStopReason::FrontierExhausted,
            stop_step: 6,
        }
    );

    // Frontier evolution, by hand: after planting (row 0) and after every
    // step (rows 1..=6), as (threshold, id, vertex, edge) sorted by the key.
    // The scan reference must produce exactly these live frontiers — the
    // stale heap entries (e2, e4 after step 3) never appear.
    let expected_frontiers: [&[(f64, usize, usize, usize)]; 7] = [
        &[(0.20, 1, 2, 1), (0.50, 0, 1, 0)],
        &[(0.50, 0, 1, 0), (0.60, 4, 3, 4)],
        &[(0.35, 3, 3, 3), (0.60, 4, 3, 4), (0.70, 2, 3, 2)],
        &[(0.10, 6, 4, 6)],
        &[(0.45, 7, 5, 7), (0.80, 8, 6, 8)],
        &[(0.80, 8, 6, 8)],
        &[],
    ];
    assert_eq!(reference.frontiers.len(), expected_frontiers.len());
    for (step, expected) in expected_frontiers.iter().enumerate() {
        let actual: Vec<(f64, usize, usize, usize)> = reference.frontiers[step]
            .iter()
            .map(|entry| (entry.threshold, entry.id, entry.vertex, entry.edge))
            .collect();
        assert_eq!(
            &actual[..],
            *expected,
            "hand bond frontier after step {step}"
        );
    }
    eprintln!("hand bond invasion: sequence, frontier, and tree verified");
}

/// Hand-worked SITE trace, seed {0}, GrowAll: the sequence differs from
/// bond mode (the site threshold field drives selection) —
///
/// - plant 0 → frontier {(v1,.55 via e0), (v2,.30 via e1)};
/// - step 1 accepts v2 (.30) through e1; frontier gains (v3,.15 via e4);
/// - step 2 accepts v3 (.15) through e4; frontier gains (v4,.65 via e6);
/// - step 3 accepts v1 (.55): its CARRYING edge is e0 — the first incidence
///   that queued it — not either parallel edge to 3 (both of whose endpoints
///   were invaded by then anyway);
/// - steps 4-6 accept v4 (.65), v5 (.25), v6 (.85).
#[test]
fn hand_worked_site_invasion_matches_manual_trace() {
    let graph = hand_graph();
    let reference = assert_production_matches_naive(
        "hand site",
        &graph,
        InvasionMode::Site,
        &HAND_SITE_THRESHOLDS,
        &[0],
        NaiveStop::GrowAll,
        None,
    );

    let mut workspace = InvasionWorkspace::new();
    let outcome = invasion_run(
        &graph,
        InvasionMode::Site,
        &HAND_SITE_THRESHOLDS,
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("hand site run is valid");

    assert_eq!(workspace.invaded_order(), &[0, 2, 3, 1, 4, 5, 6]);
    assert_eq!(workspace.accepted_edges(), &[1, 4, 0, 6, 7, 8]);
    assert_eq!(
        workspace.accepted_thresholds(),
        &[0.30, 0.15, 0.55, 0.65, 0.25, 0.85]
    );
    // By hand: 2 off 0 (e1), 3 off 2 (e4), 1 off 0 through its carrying
    // edge e0, 4 off 3, 5 and 6 off 4.
    let expected_parent = [0, 0, 0, 2, 3, 4, 4];
    assert_eq!(workspace.parent(), &expected_parent[..]);
    assert_eq!(outcome.accepted_total, 6);
    assert_eq!(outcome.max_accepted_threshold, 0.85);
    assert_eq!(outcome.stop_reason, InvasionStopReason::FrontierExhausted);
    assert_eq!(outcome.stop_step, 6);

    let expected_frontiers: [&[(f64, usize, usize, usize)]; 7] = [
        &[(0.30, 2, 2, 1), (0.55, 1, 1, 0)],
        &[(0.15, 3, 3, 4), (0.55, 1, 1, 0)],
        &[(0.55, 1, 1, 0), (0.65, 4, 4, 6)],
        &[(0.65, 4, 4, 6)],
        &[(0.25, 5, 5, 7), (0.85, 6, 6, 8)],
        &[(0.85, 6, 6, 8)],
        &[],
    ];
    assert_eq!(reference.frontiers.len(), expected_frontiers.len());
    for (step, expected) in expected_frontiers.iter().enumerate() {
        let actual: Vec<(f64, usize, usize, usize)> = reference.frontiers[step]
            .iter()
            .map(|entry| (entry.threshold, entry.id, entry.vertex, entry.edge))
            .collect();
        assert_eq!(
            &actual[..],
            *expected,
            "hand site frontier after step {step}"
        );
    }
    eprintln!("hand site invasion: sequence, frontier, and tree verified");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 2: binary heap versus complete-frontier scan
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn heap_invasion_matches_full_frontier_scan_reference() {
    let graphs = [
        ("chain", build_chain(12, false)),
        ("square", build_square(5, 5, false)),
        ("random-multigraph", random_multigraph(9, 4)),
        ("two-components", two_components()),
    ];
    let mut comparisons = 0_usize;
    for (name, graph) in &graphs {
        let vertices = graph.vertex_count();
        let edge_count = graph.n_edges();
        for mode in [InvasionMode::Site, InvasionMode::Bond] {
            for seed in 0..12_u64 {
                let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0300 + seed);
                let field_len = match mode {
                    InvasionMode::Site => vertices,
                    InvasionMode::Bond => edge_count,
                };
                // Even seeds: continuous thresholds. Odd seeds: quantized to
                // five levels so ties saturate the (threshold, id) key.
                let thresholds: Vec<f64> = (0..field_len)
                    .map(|_| {
                        let draw = rng.random::<f64>();
                        if seed % 2 == 0 {
                            draw
                        } else {
                            (draw * 5.0).floor() * 0.2
                        }
                    })
                    .collect();
                let seed_sets: [&[usize]; 2] = [&[0], &[0, vertices / 2]];
                for seeds in seed_sets {
                    for stop in [
                        NaiveStop::GrowAll,
                        NaiveStop::Steps(3),
                        NaiveStop::FirstSpanning {
                            from: 0,
                            to: vertices - 1,
                        },
                    ] {
                        let context =
                            format!("{name} {mode:?} seed {seed} seeds {seeds:?} stop {stop:?}");
                        assert_production_matches_naive(
                            &context,
                            graph,
                            mode,
                            &thresholds,
                            seeds,
                            stop,
                            None,
                        );
                        comparisons += 1;
                    }
                }
            }
        }
        eprintln!("heap-vs-scan parity: {name} ok");
    }
    assert!(comparisons >= 300, "expected a dense comparison grid");
}

#[test]
fn trapping_invasion_matches_full_frontier_scan_reference() {
    let graphs = [
        ("chain", build_chain(10, false)),
        ("random-multigraph", random_multigraph(8, 9)),
    ];
    for (name, graph) in &graphs {
        let vertices = graph.vertex_count();
        let edge_count = graph.n_edges();
        for mode in [InvasionMode::Site, InvasionMode::Bond] {
            for seed in 0..6_u64 {
                let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x7401 + seed);
                let field_len = match mode {
                    InvasionMode::Site => vertices,
                    InvasionMode::Bond => edge_count,
                };
                let thresholds: Vec<f64> = (0..field_len)
                    .map(|_| (rng.random::<f64>() * 4.0).floor() * 0.25)
                    .collect();
                let context = format!("{name} {mode:?} seed {seed} trapping exit {{last}}");
                assert_production_matches_naive(
                    &context,
                    graph,
                    mode,
                    &thresholds,
                    &[0],
                    NaiveStop::GrowAll,
                    Some(&[vertices - 1]),
                );
            }
        }
        eprintln!("trapping heap-vs-scan parity: {name} ok");
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3 (named): tree uniqueness on random realizations
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn invasion_graph_is_a_forest_with_one_tree_per_seed() {
    let graph = random_multigraph(10, 0x7ee);
    let vertices = graph.vertex_count();
    let edge_count = graph.n_edges();
    for mode in [InvasionMode::Site, InvasionMode::Bond] {
        for seed in 0..8_u64 {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
            let field_len = match mode {
                InvasionMode::Site => vertices,
                InvasionMode::Bond => edge_count,
            };
            let thresholds: Vec<f64> = (0..field_len)
                .map(|_| (rng.random::<f64>() * 3.0).floor() * 0.5)
                .collect();
            // Two seeds in the same component still grow disjoint trees.
            assert_production_matches_naive(
                "forest gate",
                &graph,
                mode,
                &thresholds,
                &[0, 5],
                NaiveStop::GrowAll,
                None,
            );
        }
    }
    eprintln!("forest invariants: verified on 16 random two-seed realizations");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 4: Prim/MST prefix relation (unique thresholds, bond, non-trapping)
// ─────────────────────────────────────────────────────────────────────────

/// Independent Prim: growth from a virtual super-root joined to every seed.
/// Each step scans ALL edges for cut-crossers (exactly one invaded
/// endpoint, self-loops excluded) and takes the minimum threshold — with
/// unique thresholds the cut minimum is unique, so this is Prim's order.
fn prim_growth_from_super_root(
    graph: &CsrLattice,
    thresholds: &[f64],
    seeds: &[usize],
) -> Vec<usize> {
    let vertices = graph.vertex_count();
    let mut active = vec![false; vertices];
    for &seed in seeds {
        active[seed] = true;
    }
    let mut accepted = Vec::new();
    loop {
        let mut best: Option<(f64, usize, usize)> = None; // (threshold, edge, vertex)
        for (edge, &threshold) in thresholds.iter().enumerate() {
            let [left, right] = [graph.edges[edge].source, graph.edges[edge].target];
            if left == right {
                continue;
            }
            match (active[left], active[right]) {
                (true, false) | (false, true) => {
                    let vertex = if active[left] { right } else { left };
                    let candidate = (threshold, edge, vertex);
                    if best.is_none() || candidate.0.total_cmp(&best.unwrap().0) == Ordering::Less {
                        best = Some(candidate);
                    }
                }
                _ => {}
            }
        }
        match best {
            Some((_, edge, vertex)) => {
                active[vertex] = true;
                accepted.push(edge);
            }
            None => return accepted,
        }
    }
}

/// Independent Kruskal: the MST edge set of the seed's connected component.
fn kruskal_mst_of_seed_component(
    graph: &CsrLattice,
    thresholds: &[f64],
    seed: usize,
) -> Vec<usize> {
    let vertices = graph.vertex_count();
    // Component membership by flood fill.
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertices];
    for edge in 0..graph.n_edges() {
        let [left, right] = [graph.edges[edge].source, graph.edges[edge].target];
        if left != right {
            adjacency[left].push(right);
            adjacency[right].push(left);
        }
    }
    let mut in_component = vec![false; vertices];
    in_component[seed] = true;
    let mut stack = vec![seed];
    while let Some(vertex) = stack.pop() {
        for &neighbor in &adjacency[vertex] {
            if !in_component[neighbor] {
                in_component[neighbor] = true;
                stack.push(neighbor);
            }
        }
    }
    // Kruskal restricted to edges inside the component.
    let mut order: Vec<usize> = (0..graph.n_edges())
        .filter(|&edge| {
            let [left, right] = [graph.edges[edge].source, graph.edges[edge].target];
            left != right && in_component[left] && in_component[right]
        })
        .collect();
    order.sort_by(|&a, &b| thresholds[a].total_cmp(&thresholds[b]).then(a.cmp(&b)));
    let mut classes: Vec<usize> = (0..vertices).collect();
    let root = |classes: &mut Vec<usize>, mut node: usize| {
        while classes[node] != node {
            node = classes[node];
        }
        node
    };
    let mut tree = Vec::new();
    for edge in order {
        let [left, right] = [graph.edges[edge].source, graph.edges[edge].target];
        let (left_root, right_root) = (root(&mut classes, left), root(&mut classes, right));
        if left_root != right_root {
            classes[left_root] = right_root;
            tree.push(edge);
        }
    }
    tree.sort_unstable();
    tree
}

/// Dyadic unique thresholds: a seeded shuffle of distinct dyadic values so
/// every comparison is exact and every threshold is unique.
fn unique_thresholds(edge_count: usize, seed: u64) -> Vec<f64> {
    let mut ranks: Vec<usize> = (0..edge_count).collect();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    for index in (1..ranks.len()).rev() {
        let swap_with = rng.random_range(0..=index);
        ranks.swap(index, swap_with);
    }
    ranks
        .iter()
        .map(|&rank| 0.03125 * (rank as f64 + 1.0))
        .collect()
}

/// The precise relation, verified: with unique thresholds a non-trapping
/// bond `GrowAll` realization accepts the Prim growth order (single seed:
/// Prim from the seed; multiple seeds: Prim from a virtual super-root
/// joined to the seeds), the terminal tree/forest is the unique MST of the
/// seed's component (single seed, Kruskal-checked), and a `Steps(n)` run is
/// the exact n-edge prefix.
#[test]
fn nontrapping_bond_invasion_is_prim_growth_order() {
    let graphs = [
        ("chain", build_chain(10, false)),
        ("square", build_square(5, 5, false)),
        ("random-multigraph", random_multigraph(9, 0x057)),
    ];
    for (name, graph) in &graphs {
        let edge_count = graph.n_edges();
        for seed in 0..4_u64 {
            let thresholds = unique_thresholds(edge_count, 0x9e11 + seed);
            let mut workspace = InvasionWorkspace::new();
            invasion_run(
                graph,
                InvasionMode::Bond,
                &thresholds,
                &[0],
                InvasionStop::GrowAll,
                &NON_TRAPPING,
                &mut workspace,
            )
            .expect("unique-threshold bond run is valid");

            // Accepted sequence is exactly the Prim growth order from 0.
            let prim = prim_growth_from_super_root(graph, &thresholds, &[0]);
            assert_eq!(
                workspace.accepted_edges(),
                &prim[..],
                "{name} seed {seed}: bond invasion must be Prim's growth order"
            );

            // Terminal tree is the unique MST of the seed's component.
            let mut produced = workspace.accepted_edges().to_vec();
            produced.sort_unstable();
            assert_eq!(
                produced,
                kruskal_mst_of_seed_component(graph, &thresholds, 0),
                "{name} seed {seed}: invasion tree must be the unique MST"
            );

            // Steps(n) is the exact n-edge prefix (when the frontier lasts).
            let prefix = prim.len().min(3);
            let mut stopped = InvasionWorkspace::new();
            let stopped_outcome = invasion_run(
                graph,
                InvasionMode::Bond,
                &thresholds,
                &[0],
                InvasionStop::Steps(prefix),
                &NON_TRAPPING,
                &mut stopped,
            )
            .expect("prefix run is valid");
            assert_eq!(
                stopped.accepted_edges(),
                &prim[..prefix],
                "{name} seed {seed}: Steps({prefix}) must be the Prim prefix"
            );
            assert_eq!(stopped_outcome.stop_reason, InvasionStopReason::StepBudget);

            // Multiple seeds: Prim from the virtual super-root.
            let mut multi = InvasionWorkspace::new();
            invasion_run(
                graph,
                InvasionMode::Bond,
                &thresholds,
                &[0, graph.vertex_count() / 2],
                InvasionStop::GrowAll,
                &NON_TRAPPING,
                &mut multi,
            )
            .expect("multi-seed run is valid");
            assert_eq!(
                multi.accepted_edges(),
                &prim_growth_from_super_root(graph, &thresholds, &[0, graph.vertex_count() / 2])[..],
                "{name} seed {seed}: multi-seed bond invasion must be super-root Prim"
            );
        }
        eprintln!("Prim/MST relation: {name} ok");
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Explicit tie policy
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn tie_policy_resolves_by_lowest_dense_id() {
    let graph = star();

    // Site mode: three equal site thresholds -> vertex-ID order 1, 2, 3.
    let site_tied = [0.0, 0.5, 0.5, 0.5];
    let mut workspace = InvasionWorkspace::new();
    invasion_run(
        &graph,
        InvasionMode::Site,
        &site_tied,
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("site tie run is valid");
    assert_eq!(workspace.invaded_order(), &[0, 1, 2, 3]);
    assert_eq!(workspace.accepted_edges(), &[0, 1, 2]);

    // Bond mode: three equal edge thresholds -> edge-ID order e0, e1, e2.
    let bond_tied = [0.5, 0.5, 0.5];
    invasion_run(
        &graph,
        InvasionMode::Bond,
        &bond_tied,
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("bond tie run is valid");
    assert_eq!(workspace.invaded_order(), &[0, 1, 2, 3]);
    assert_eq!(workspace.accepted_edges(), &[0, 1, 2]);

    // A two-way tie with a distinct third entry: only the tied pair swaps —
    // e1 (id 1) beats e2 (id 2) despite equal thresholds; site symmetric.
    let bond_partial = [0.9, 0.5, 0.5];
    invasion_run(
        &graph,
        InvasionMode::Bond,
        &bond_partial,
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("bond partial-tie run is valid");
    assert_eq!(workspace.invaded_order(), &[0, 2, 3, 1]);
    assert_eq!(workspace.accepted_edges(), &[1, 2, 0]);

    let site_partial = [0.0, 0.9, 0.5, 0.5];
    invasion_run(
        &graph,
        InvasionMode::Site,
        &site_partial,
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("site partial-tie run is valid");
    assert_eq!(workspace.invaded_order(), &[0, 2, 3, 1]);
    assert_eq!(workspace.accepted_edges(), &[1, 2, 0]);

    // Documented `total_cmp` fine print: -0.0 < +0.0 are DISTINCT keys, so
    // the -0.0 site wins even though its vertex ID (2) is larger than the
    // +0.0 site's (1) — an ID tie-break would give the opposite order.
    let signed_zeros = [1.0, 0.0, -0.0, 1.0];
    invasion_run(
        &graph,
        InvasionMode::Site,
        &signed_zeros,
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("signed-zero run is valid");
    assert_eq!(
        workspace.invaded_order(),
        &[0, 2, 1, 3],
        "-0.0 must precede +0.0 under total_cmp (distinct keys, not a tie)"
    );
    assert_eq!(workspace.accepted_thresholds(), &[-0.0, 0.0, 1.0]);
}

// ─────────────────────────────────────────────────────────────────────────
// Explicit stop policies
// ─────────────────────────────────────────────────────────────────────────

/// BFS over the tree prefix of the first `steps` accepted edges: is `goal`
/// reachable from `start`?
fn bfs_reaches_over_prefix(
    graph: &CsrLattice,
    accepted: &[usize],
    steps: usize,
    start: usize,
    goal: usize,
) -> bool {
    let vertices = graph.vertex_count();
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertices];
    for &edge in &accepted[..steps] {
        let [left, right] = [graph.edges[edge].source, graph.edges[edge].target];
        adjacency[left].push(right);
        adjacency[right].push(left);
    }
    let mut visited = vec![false; vertices];
    visited[start] = true;
    let mut stack = vec![start];
    while let Some(vertex) = stack.pop() {
        if vertex == goal {
            return true;
        }
        for &neighbor in &adjacency[vertex] {
            if !visited[neighbor] {
                visited[neighbor] = true;
                stack.push(neighbor);
            }
        }
    }
    false
}

#[test]
fn stop_policies_are_exact() {
    let graph = hand_graph();
    let mut workspace = InvasionWorkspace::new();

    // GrowAll reference (hand-verified elsewhere): full 6 steps.
    let full = invasion_run(
        &graph,
        InvasionMode::Bond,
        &HAND_BOND_THRESHOLDS,
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("grow-all run is valid");
    assert_eq!(full.stop_reason, InvasionStopReason::FrontierExhausted);
    assert_eq!(full.stop_step, 6);
    let full_edges = workspace.accepted_edges().to_vec();

    // Steps: exact budgets, including zero; an over-large budget falls
    // through to frontier exhaustion.
    for (budget, expected) in [(0, &[][..]), (3, &full_edges[..3]), (6, &full_edges[..6])] {
        let outcome = invasion_run(
            &graph,
            InvasionMode::Bond,
            &HAND_BOND_THRESHOLDS,
            &[0],
            InvasionStop::Steps(budget),
            &NON_TRAPPING,
            &mut workspace,
        )
        .expect("steps run is valid");
        assert_eq!(
            workspace.accepted_edges(),
            expected,
            "Steps({budget}) edges"
        );
        assert_eq!(outcome.accepted_total, expected.len());
        assert_eq!(outcome.stop_step, budget);
        assert_eq!(outcome.stop_reason, InvasionStopReason::StepBudget);
    }
    let over = invasion_run(
        &graph,
        InvasionMode::Bond,
        &HAND_BOND_THRESHOLDS,
        &[0],
        InvasionStop::Steps(100),
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("over-budget run is valid");
    assert_eq!(over.stop_reason, InvasionStopReason::FrontierExhausted);
    assert_eq!(over.accepted_total, 6);

    // FirstSpanning: for every target the stop must fire at exactly the
    // BFS-verified spanning step — the first prefix of accepted edges over
    // which `to` is reachable from the seed.
    for to in 1..7_usize {
        let outcome = invasion_run(
            &graph,
            InvasionMode::Bond,
            &HAND_BOND_THRESHOLDS,
            &[0],
            InvasionStop::FirstSpanning { from: 0, to },
            &NON_TRAPPING,
            &mut workspace,
        )
        .expect("spanning run is valid");
        assert_eq!(
            outcome.stop_reason,
            InvasionStopReason::SpanningStep,
            "to = {to}"
        );
        let step = outcome.stop_step;
        assert!(
            bfs_reaches_over_prefix(&graph, &full_edges, step, 0, to),
            "to = {to}: spanning must hold at the stop step"
        );
        assert!(
            step == 0 || !bfs_reaches_over_prefix(&graph, &full_edges, step - 1, 0, to),
            "to = {to}: spanning must fail one step earlier"
        );
        assert_eq!(workspace.accepted_edges(), &full_edges[..step]);
    }

    // Two seeds satisfying the policy at planting: step-zero spanning stop.
    let seeds = [0_usize, 3];
    let zero = invasion_run(
        &graph,
        InvasionMode::Bond,
        &HAND_BOND_THRESHOLDS,
        &seeds,
        InvasionStop::FirstSpanning { from: 0, to: 3 },
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("seed spanning run is valid");
    assert_eq!(zero.stop_reason, InvasionStopReason::SpanningStep);
    assert_eq!(zero.stop_step, 0);
    assert_eq!(zero.accepted_total, 0);
    assert_eq!(zero.invaded_vertices, 2);

    // An unreachable endpoint exhausts the frontier instead.
    let split = two_components();
    let exhausted = invasion_run(
        &split,
        InvasionMode::Bond,
        &[0.1, 0.2, 0.3, 0.4, 0.5],
        &[0],
        InvasionStop::FirstSpanning { from: 0, to: 4 },
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("unreachable spanning run is valid");
    assert_eq!(exhausted.stop_reason, InvasionStopReason::FrontierExhausted);
    assert_eq!(exhausted.accepted_total, 2, "only the triangle invades");
    eprintln!("stop policies: GrowAll/Steps/FirstSpanning exact");
}

// ─────────────────────────────────────────────────────────────────────────
// Wilkinson-Willemsen trapping (hand-worked)
// ─────────────────────────────────────────────────────────────────────────

/// Chain 0-1-2-3-4, seed {1}, site thresholds [.1, _, .2, .3, .4]:
///
/// - exit {0}: v0 (.1) is itself the outlet and invades; then v2's
///   defending region {2,3,4} cannot reach the (now invaded) outlet —
///   v2 is trapped and retired permanently, and the frontier empties.
/// - exit {4}: v0 (.1) is trapped FIRST (its region {0} reaches nothing);
///   the chain 2-3-4 invades normally, the outlet 4 last.
#[test]
fn trapping_hand_worked_chain_matches_manual_trace() {
    let graph = chain5();
    let thresholds = [0.1, 0.0, 0.2, 0.3, 0.4];
    let mut workspace = InvasionWorkspace::new();

    // Contrast: without trapping everything invades.
    let free = invasion_run(
        &graph,
        InvasionMode::Site,
        &thresholds,
        &[1],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("non-trapping chain run is valid");
    assert_eq!(workspace.invaded_order(), &[1, 0, 2, 3, 4]);
    assert_eq!(free.trapped_sites, 0);

    // Exit {0}: invade the outlet, then the branch is sealed.
    let sealed = invasion_run(
        &graph,
        InvasionMode::Site,
        &thresholds,
        &[1],
        InvasionStop::GrowAll,
        &trapping(&[0]),
        &mut workspace,
    )
    .expect("sealed trapping run is valid");
    assert_eq!(workspace.invaded_order(), &[1, 0]);
    assert_eq!(workspace.accepted_edges(), &[0]);
    assert_eq!(
        sealed,
        InvasionOutcome {
            seed_count: 1,
            accepted_total: 1,
            invaded_vertices: 2,
            trapped_sites: 1,
            max_accepted_threshold: 0.1,
            stop_reason: InvasionStopReason::FrontierExhausted,
            stop_step: 1,
        }
    );

    // Exit {4}: the near side is trapped, the far side drains to the outlet.
    let drained = invasion_run(
        &graph,
        InvasionMode::Site,
        &thresholds,
        &[1],
        InvasionStop::GrowAll,
        &trapping(&[4]),
        &mut workspace,
    )
    .expect("draining trapping run is valid");
    assert_eq!(workspace.invaded_order(), &[1, 2, 3, 4]);
    assert_eq!(workspace.accepted_edges(), &[1, 2, 3]);
    assert_eq!(drained.trapped_sites, 1);
    assert_eq!(drained.max_accepted_threshold, 0.4);

    assert_production_matches_naive(
        "trapping chain exit 0",
        &graph,
        InvasionMode::Site,
        &thresholds,
        &[1],
        NaiveStop::GrowAll,
        Some(&[0]),
    );
    assert_production_matches_naive(
        "trapping chain exit 4",
        &graph,
        InvasionMode::Site,
        &thresholds,
        &[1],
        NaiveStop::GrowAll,
        Some(&[4]),
    );
    eprintln!("trapping chain: both outlet placements hand-verified");
}

/// Bond-mode trapping with a parallel edge: retiring v2 must consume BOTH
/// frontier entries that target it (e1 .15 and e2 .20) without growth.
#[test]
fn trapping_bond_mode_retires_parallel_entries() {
    // Chain 1-2 doubled: e0(0,1), e1(1,2), e2(1,2), e3(2,3), e4(3,4).
    let graph = CsrLattice::try_from_edges(
        5,
        vec![bond(0, 1), bond(1, 2), bond(1, 2), bond(2, 3), bond(3, 4)],
    )
    .expect("doubled chain is valid");
    let thresholds = [0.1, 0.15, 0.20, 0.30, 0.40];
    let mut workspace = InvasionWorkspace::new();
    let outcome = invasion_run(
        &graph,
        InvasionMode::Bond,
        &thresholds,
        &[1],
        InvasionStop::GrowAll,
        &trapping(&[0]),
        &mut workspace,
    )
    .expect("bond trapping run is valid");
    // e0 invades the outlet 0; the -0.-side has no un-invaded outlet left:
    // e1's check retires v2, e2's pop must then skip the retired vertex.
    assert_eq!(workspace.invaded_order(), &[1, 0]);
    assert_eq!(workspace.accepted_edges(), &[0]);
    assert_eq!(outcome.trapped_sites, 1);
    assert_eq!(outcome.stop_reason, InvasionStopReason::FrontierExhausted);

    assert_production_matches_naive(
        "bond trapping parallel",
        &graph,
        InvasionMode::Bond,
        &thresholds,
        &[1],
        NaiveStop::GrowAll,
        Some(&[0]),
    );
    eprintln!("bond trapping: retired site consumes parallel entries");
}

// ─────────────────────────────────────────────────────────────────────────
// G1 substrate drive
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn substrate_trajectory_steps_are_invasion_times() {
    let graph = hand_graph();
    let mut workspace = InvasionWorkspace::new();
    invasion_run(
        &graph,
        InvasionMode::Bond,
        &HAND_BOND_THRESHOLDS,
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("substrate drive run is valid");

    let trajectory = workspace.substrate().trajectory();
    // Seeds record step 0; growth step i records the vertex and the edge,
    // both stamped with the 1-based invasion time.
    let expected_steps: Vec<usize> = std::iter::once(0)
        .chain((1..=6).flat_map(|step| [step, step]))
        .collect();
    assert_eq!(trajectory.steps(), &expected_steps[..]);
    let terminal = *trajectory.summaries().last().expect("nonempty trajectory");
    assert_eq!(terminal.active_vertex_count, 7);
    assert_eq!(terminal.active_edge_count, 6);
    assert_eq!(terminal.component_count, 1);
    assert_eq!(terminal.largest_component_size, 7);
    // Online: after the i-th COMPLETE growth (vertex op then edge op —
    // trajectory indices 2i and 2i+1) the active counts are 1 + i vertices
    // and i edges; the vertex-op summary inside step i still sees i - 1
    // edges because the edge activates after its vertex.
    for (index, summary) in trajectory.summaries().iter().enumerate() {
        if index % 2 == 0 {
            let growth = index / 2;
            assert_eq!(summary.active_vertex_count, 1 + growth);
            assert_eq!(summary.active_edge_count, growth);
        }
    }

    // Two seeds in one component grow two trees that never merge: the
    // substrate (which unions by accepted edges) sees two components.
    // Hand trace: frontier {e6 .10, e1 .20, e3 .35, e7+, e0 .50, e4 .60,
    // e2 .70, e8 .80} — 4 via e6, 2 via e1 (e4 then goes stale against the
    // seed 3), 1 via the parallel e3 (.35 beats e0 .50 from seed 0), then
    // e7, e8.
    let mut multi = InvasionWorkspace::new();
    let outcome = invasion_run(
        &graph,
        InvasionMode::Bond,
        &HAND_BOND_THRESHOLDS,
        &[0, 3],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut multi,
    )
    .expect("multi-seed substrate run is valid");
    assert_eq!(multi.invaded_order(), &[0, 3, 4, 2, 1, 5, 6]);
    assert_eq!(multi.accepted_edges(), &[6, 1, 3, 7, 8]);
    assert_eq!(outcome.accepted_total, 5);
    let terminal = *multi.substrate().trajectory().summaries().last().unwrap();
    assert_eq!(terminal.active_vertex_count, 7);
    assert_eq!(terminal.active_edge_count, 5);
    assert_eq!(
        terminal.component_count, 2,
        "the invasion forest has two trees"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Adapter: RNG draw contract, C0 history, and reuse
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn adapter_draw_contract_is_bitwise_reproducible() {
    let graph = hand_graph();
    for (mode, field_len) in [(InvasionMode::Bond, 9_usize), (InvasionMode::Site, 7_usize)] {
        let mut adapter = InvasionPercolationMC::new(
            graph.clone(),
            mode,
            vec![0],
            InvasionStop::GrowAll,
            NON_TRAPPING,
            6,
        )
        .expect("adapter constructs");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x63_33);
        adapter.realize(&mut rng);

        // Bitwise replay: the same RNG stream, drawn one f64 per element in
        // ascending dense order, fed to the RNG-free core.
        let mut replay = Xoshiro256PlusPlus::seed_from_u64(0x63_33);
        let thresholds: Vec<f64> = (0..field_len).map(|_| replay.random::<f64>()).collect();
        let mut workspace = InvasionWorkspace::new();
        let outcome = invasion_run(
            &graph,
            mode,
            &thresholds,
            &[0],
            InvasionStop::GrowAll,
            &NON_TRAPPING,
            &mut workspace,
        )
        .expect("replayed core run is valid");
        assert_eq!(adapter.last_outcome(), outcome, "{mode:?}: outcome");
        assert_eq!(
            adapter.last_order(),
            workspace.invaded_order(),
            "{mode:?}: order"
        );
        assert_eq!(
            adapter.last_accepted_edges(),
            workspace.accepted_edges(),
            "{mode:?}: edges"
        );
        assert_eq!(
            adapter.last_accepted_thresholds(),
            workspace.accepted_thresholds(),
            "{mode:?}: thresholds"
        );
    }
}

#[test]
fn adapter_history_is_one_fixed_shape_vector_sample() {
    let graph = hand_graph();

    // Shorter than the realization: truncate to the first `steps`.
    let mut short = InvasionPercolationMC::new(
        graph.clone(),
        InvasionMode::Bond,
        vec![0],
        InvasionStop::GrowAll,
        NON_TRAPPING,
        3,
    )
    .expect("short adapter constructs");
    short.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x31_57));
    let full_thresholds = short.last_accepted_thresholds().to_vec();
    assert_eq!(full_thresholds.len(), 6);
    assert_eq!(short.last_history(), &full_thresholds[..3]);
    assert_eq!(short.steps(), 3);

    // Longer: repeat the final accepted threshold to hold the shape.
    let mut long = InvasionPercolationMC::new(
        graph.clone(),
        InvasionMode::Bond,
        vec![0],
        InvasionStop::GrowAll,
        NON_TRAPPING,
        10,
    )
    .expect("long adapter constructs");
    long.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x31_57));
    let mut expected = full_thresholds.clone();
    expected.resize(10, *full_thresholds.last().expect("six accepted"));
    assert_eq!(long.last_history(), &expected[..]);
    assert_eq!(long.last_outcome().accepted_total, 6);

    // Nothing accepted: all-zero history (Steps(0) accepts nothing).
    let mut empty = InvasionPercolationMC::new(
        graph,
        InvasionMode::Bond,
        vec![0],
        InvasionStop::Steps(0),
        NON_TRAPPING,
        4,
    )
    .expect("empty adapter constructs");
    empty.realize(&mut Xoshiro256PlusPlus::seed_from_u64(0x31_57));
    assert_eq!(empty.last_history(), &[0.0, 0.0, 0.0, 0.0]);
    assert_eq!(empty.last_outcome().accepted_total, 0);
    assert_eq!(empty.last_outcome().max_accepted_threshold, 0.0);
}

#[test]
fn adapter_sweep_and_measure_record_through_context() {
    use carlo_rs::{Context, MonteCarlo};
    let graph = hand_graph();
    let mut adapter = InvasionPercolationMC::new(
        graph,
        InvasionMode::Site,
        vec![0],
        InvasionStop::GrowAll,
        NON_TRAPPING,
        5,
    )
    .expect("adapter constructs");
    let mut context =
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x5_33d), 0, 1 << 16);
    for _ in 0..4 {
        adapter.sweep(&mut context);
        adapter.measure(&mut context);
        assert_eq!(adapter.last_history().len(), 5);
    }
}

#[test]
fn workspace_reuse_across_runs_and_policies_leaves_no_stale_state() {
    let graph = hand_graph();
    let split = two_components();
    let mut workspace = InvasionWorkspace::new();

    // Identical inputs on a reused workspace give identical results.
    let mut first = Vec::new();
    for _ in 0..2 {
        let outcome = invasion_run(
            &graph,
            InvasionMode::Site,
            &HAND_SITE_THRESHOLDS,
            &[0],
            InvasionStop::GrowAll,
            &NON_TRAPPING,
            &mut workspace,
        )
        .expect("reuse site run is valid");
        first.push((
            outcome,
            workspace.invaded_order().to_vec(),
            workspace.accepted_edges().to_vec(),
        ));
    }
    assert_eq!(first[0], first[1], "reused workspace must be deterministic");

    // A trapping run followed by a non-trapping run on the SAME workspace:
    // the exit marks planted by the first must not leak into the second.
    invasion_run(
        &graph,
        InvasionMode::Site,
        &HAND_SITE_THRESHOLDS,
        &[1],
        InvasionStop::GrowAll,
        &trapping(&[0]),
        &mut workspace,
    )
    .expect("trapping run on reused workspace is valid");
    let after = invasion_run(
        &graph,
        InvasionMode::Site,
        &HAND_SITE_THRESHOLDS,
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("post-trapping run is valid");
    assert_eq!(after, first[0].0, "trapping marks must reset between runs");

    // Re-prepare to a different topology: fresh, correct behavior.
    // By hand: seed 0, frontier {e2 (.2)->2, e0 (.3)->1}; e2 first, then e1
    // (.1) beats e0 (.3) for vertex 1.
    let outcome = invasion_run(
        &split,
        InvasionMode::Bond,
        &[0.3, 0.1, 0.2, 0.4, 0.5],
        &[0],
        InvasionStop::GrowAll,
        &NON_TRAPPING,
        &mut workspace,
    )
    .expect("re-prepared run is valid");
    assert_eq!(outcome.invaded_vertices, 3);
    assert_eq!(workspace.invaded_order(), &[0, 2, 1]);
    assert_eq!(workspace.accepted_edges(), &[2, 1]);
}

// ─────────────────────────────────────────────────────────────────────────
// Typed rejections
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn invalid_inputs_are_typed_rejections() {
    let graph = hand_graph();
    let mut workspace = InvasionWorkspace::new();
    let mut run = |thresholds: &[f64],
                   seeds: &[usize],
                   stop: InvasionStop,
                   trapping: &InvasionTrapping|
     -> Result<InvasionOutcome, InvasionError> {
        invasion_run(
            &graph,
            InvasionMode::Site,
            thresholds,
            seeds,
            stop,
            trapping,
            &mut workspace,
        )
    };

    // Threshold field: wrong length, NaN, and infinity.
    assert!(matches!(
        run(&[0.1; 9], &[0], InvasionStop::GrowAll, &NON_TRAPPING),
        Err(InvasionError::ThresholdLengthMismatch {
            expected: 7,
            actual: 9
        })
    ));
    let mut nan_field = HAND_SITE_THRESHOLDS.to_vec();
    nan_field[2] = f64::NAN;
    assert!(matches!(
        run(&nan_field, &[0], InvasionStop::GrowAll, &NON_TRAPPING),
        Err(InvasionError::ThresholdNotFinite { index: 2, .. })
    ));
    let mut inf_field = HAND_SITE_THRESHOLDS.to_vec();
    inf_field[0] = f64::INFINITY;
    assert!(matches!(
        run(&inf_field, &[0], InvasionStop::GrowAll, &NON_TRAPPING),
        Err(InvasionError::ThresholdNotFinite { index: 0, .. })
    ));

    // Seeds: empty, out of domain, duplicated.
    assert!(matches!(
        run(
            &HAND_SITE_THRESHOLDS,
            &[],
            InvasionStop::GrowAll,
            &NON_TRAPPING
        ),
        Err(InvasionError::EmptySeeds)
    ));
    assert!(matches!(
        run(
            &HAND_SITE_THRESHOLDS,
            &[7],
            InvasionStop::GrowAll,
            &NON_TRAPPING
        ),
        Err(InvasionError::SeedDomain {
            seed: 7,
            vertices: 7
        })
    ));
    assert!(matches!(
        run(
            &HAND_SITE_THRESHOLDS,
            &[0, 0],
            InvasionStop::GrowAll,
            &NON_TRAPPING
        ),
        Err(InvasionError::DuplicateSeed { seed: 0 })
    ));

    // Spanning endpoints: out of domain.
    assert!(matches!(
        run(
            &HAND_SITE_THRESHOLDS,
            &[0],
            InvasionStop::FirstSpanning { from: 0, to: 7 },
            &NON_TRAPPING
        ),
        Err(InvasionError::SpanningDomain {
            vertex: 7,
            vertices: 7
        })
    ));

    // Trapping exits: empty, out of domain, duplicated.
    assert!(matches!(
        run(
            &HAND_SITE_THRESHOLDS,
            &[0],
            InvasionStop::GrowAll,
            &trapping(&[])
        ),
        Err(InvasionError::EmptyExitSet)
    ));
    assert!(matches!(
        run(
            &HAND_SITE_THRESHOLDS,
            &[0],
            InvasionStop::GrowAll,
            &trapping(&[7])
        ),
        Err(InvasionError::ExitDomain {
            vertex: 7,
            vertices: 7
        })
    ));
    assert!(matches!(
        run(
            &HAND_SITE_THRESHOLDS,
            &[0],
            InvasionStop::GrowAll,
            &trapping(&[1, 1])
        ),
        Err(InvasionError::DuplicateExit { vertex: 1 })
    ));

    // Empty graph: a zero-vertex borrowed view is a typed rejection.
    let empty = BorrowedUndirectedCsr::new(&[0], &[], &[], &[])
        .expect("zero-vertex borrowed view is valid");
    assert!(matches!(
        invasion_run(
            &empty,
            InvasionMode::Site,
            &[],
            &[0],
            InvasionStop::GrowAll,
            &NON_TRAPPING,
            &mut workspace
        ),
        Err(InvasionError::EmptyGraph)
    ));

    // Adapter: a zero-length recorded array is rejected at construction.
    assert!(matches!(
        InvasionPercolationMC::new(
            graph,
            InvasionMode::Site,
            vec![0],
            InvasionStop::GrowAll,
            NON_TRAPPING,
            0,
        ),
        Err(InvasionError::NoSteps)
    ));
}
