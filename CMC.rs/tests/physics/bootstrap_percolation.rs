//! E2 production validation: bootstrap percolation (activation-threshold
//! closure).
//!
//! Gates (roadmap §12 E2), each against independent test-side references
//! (no production analysis code):
//!
//! 1. **Exhaustive 2^V initial states**: path, cycle, star, and complete
//!    graphs; every seed mask x every homogeneous threshold r in
//!    `0..=max_degree + 1` x several seeded heterogeneous threshold
//!    vectors; both clocks (`SynchronousRounds`, `QueueClosure`) against
//!    a repeated-full-scan closure reference, including per-round
//!    trajectories against a simultaneous-update sweep reference.
//! 2. **`r = 0` activates all vertices** (both clocks; threshold-0
//!    vertices activate at the first update step).
//! 3. **`r > max_degree` preserves only the initial seeds**.
//! 4. **Order independence**: the FIFO queue closure, a
//!    one-random-qualifying-vertex-at-a-time closure, and a
//!    shuffled-sweep closure (both shuffled by a test RNG) reach the
//!    identical fixed point on random multigraphs with random seeds and
//!    thresholds, and equal the synchronous-rounds closure.
//! 5. **Hand-worked heterogeneous example**: path 0-1-2-3-4, thresholds
//!    `[2,1,1,2,3]`, seeds `{0,4}` — closure trajectory `[2,3,4,5]` in 3
//!    rounds / 3 events, hand-derived in the header of this file's
//!    reference section.
//! 6. **Adapter contracts**: fixed seeds reproduce the core bitwise; the
//!    documented Density draw contract (one f64 per vertex, ascending)
//!    replays bitwise; `measure` records ONE C0 array sample of fixed
//!    shape per realization plus the closure scalars.
//! 7. **Literature anchor** (archived-paper policy): Holroyd, "Sharp
//!    metastability threshold for two-dimensional bootstrap percolation",
//!    Probab. Theory Relat. Fields 125, 195-224 (2003)
//!    (arXiv:math/0206132, archived in `docs/ref/`): 2D r = 2 bootstrap on
//!    the periodic L x L torus. Theorem 1 gives `I(L,p) -> 1` iff
//!    `p ln L > pi^2/18`, so the finite-size seed threshold obeys
//!    `p_c(L) ~ pi^2/(18 ln L)`. Holroyd's convergence is notoriously slow
//!    (finite-L thresholds sit far below the asymptotic constant), so all
//!    gates are statistical brackets with real error bars, never exact
//!    finite-L threshold claims: (a) his own observable `I(L,p)` (full
//!    torus activation) is near-certain at `p+ = pi^2/(18 ln L)` and rare
//!    at `p- = p+/4` for L in {32, 64, 128}; (b) the measured
//!    half-majority density `p_hat(L)` decays monotonically with L (the
//!    1/ln L law), keeps `p_hat ln L` strictly inside `(lambda/4,
//!    lambda)`, and the metastable window is sharper than a factor 1.5 in
//!    p at every measured size.

use carlo_rs::{Context, Estimate, MonteCarlo};
use cmc_rs::{
    build_chain, build_square, Bond, BondType, BootstrapDynamics, BootstrapError,
    BootstrapPercolationMC, BootstrapSeeds, BootstrapThreshold, BootstrapWorkspace, CsrLattice,
    GraphView, VertexActivity, VertexId,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Realizations per (L, p) anchor measurement (env-tunable; the default
/// keeps the full gate under ~120 s in debug builds).
const ANCHOR_REALIZATIONS: usize = 400;

fn anchor_realizations() -> usize {
    std::env::var("SCUTTLE_BOOTSTRAP_SEEDS")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(ANCHOR_REALIZATIONS)
        .max(40)
}

// ─────────────────────────────────────────────────────────────────────────
// Independent references (hand-written; no production code).
// ─────────────────────────────────────────────────────────────────────────

/// Dense byte mask implementing [`VertexActivity`] for driving production.
struct MaskSeeds {
    active: Vec<u8>,
}

impl MaskSeeds {
    fn new(bits: &[u8]) -> Self {
        Self {
            active: bits.to_vec(),
        }
    }
}

impl VertexActivity for MaskSeeds {
    fn vertex_count(&self) -> usize {
        self.active.len()
    }

    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.active[vertex.index()] != 0
    }
}

/// Repeated full-scan closure: sweep the whole vertex list until nothing
/// changes (the roadmap's independent reference dynamics).
fn full_scan_closure(adjacency: &[Vec<usize>], seeds: &[bool], thresholds: &[u32]) -> Vec<bool> {
    let mut active = seeds.to_vec();
    let mut changed = true;
    while changed {
        changed = false;
        for (vertex, neighbors) in adjacency.iter().enumerate() {
            if !active[vertex] {
                let count = neighbors
                    .iter()
                    .filter(|&&neighbor| active[neighbor])
                    .count() as u32;
                if count >= thresholds[vertex] {
                    active[vertex] = true;
                    changed = true;
                }
            }
        }
    }
    active
}

/// Simultaneous-update sweep reference returning the per-round active
/// counts `[|A_0|, |A_1|, ...]` (round `t + 1` activates every inactive
/// vertex that meets its threshold against the round-`t` state).
fn full_scan_rounds(adjacency: &[Vec<usize>], seeds: &[bool], thresholds: &[u32]) -> Vec<usize> {
    let mut active = seeds.to_vec();
    let mut trajectory = vec![active.iter().filter(|&&flag| flag).count()];
    loop {
        let mut newly: Vec<usize> = Vec::new();
        for (vertex, neighbors) in adjacency.iter().enumerate() {
            if !active[vertex] {
                let count = neighbors
                    .iter()
                    .filter(|&&neighbor| active[neighbor])
                    .count() as u32;
                if count >= thresholds[vertex] {
                    newly.push(vertex);
                }
            }
        }
        if newly.is_empty() {
            return trajectory;
        }
        for vertex in newly {
            active[vertex] = true;
        }
        trajectory.push(active.iter().filter(|&&flag| flag).count());
    }
}

/// Order-randomized closure: activate ONE uniformly chosen qualifying
/// vertex per step (Fisher-Yates on the qualifying list with the test
/// RNG). The fixed point must not depend on this order.
fn one_at_a_time_closure(
    adjacency: &[Vec<usize>],
    seeds: &[bool],
    thresholds: &[u32],
    rng: &mut Xoshiro256PlusPlus,
) -> Vec<bool> {
    let mut active = seeds.to_vec();
    loop {
        let mut candidates: Vec<usize> = (0..adjacency.len())
            .filter(|&vertex| {
                !active[vertex]
                    && adjacency[vertex]
                        .iter()
                        .filter(|&&neighbor| active[neighbor])
                        .count() as u32
                        >= thresholds[vertex]
            })
            .collect();
        if candidates.is_empty() {
            return active;
        }
        // Fisher-Yates shuffle, then activate the first candidate.
        for position in (1..candidates.len()).rev() {
            let pick = rng.random_range(0..=position);
            candidates.swap(position, pick);
        }
        active[candidates[0]] = true;
    }
}

/// Shuffled-sweep closure: full scans whose vertex order is reshuffled by
/// the test RNG every sweep.
fn shuffled_sweep_closure(
    adjacency: &[Vec<usize>],
    seeds: &[bool],
    thresholds: &[u32],
    rng: &mut Xoshiro256PlusPlus,
) -> Vec<bool> {
    let mut active = seeds.to_vec();
    let mut order: Vec<usize> = (0..adjacency.len()).collect();
    loop {
        for position in (1..order.len()).rev() {
            let pick = rng.random_range(0..=position);
            order.swap(position, pick);
        }
        let mut changed = false;
        for &vertex in &order {
            if !active[vertex] {
                let count = adjacency[vertex]
                    .iter()
                    .filter(|&&neighbor| active[neighbor])
                    .count() as u32;
                if count >= thresholds[vertex] {
                    active[vertex] = true;
                    changed = true;
                }
            }
        }
        if !changed {
            return active;
        }
    }
}

fn active_set(workspace: &BootstrapWorkspace, vertices: usize) -> Vec<bool> {
    (0..vertices)
        .map(|index| workspace.is_active(index))
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────
// Hand-built topologies: reference adjacency lists (from scratch) paired
// with production lattices over the same edge sets.
// ─────────────────────────────────────────────────────────────────────────

fn path_parts(n: usize) -> (Vec<Vec<usize>>, CsrLattice) {
    let mut adjacency = vec![Vec::new(); n];
    let mut edges = Vec::new();
    for site in 0..n.saturating_sub(1) {
        adjacency[site].push(site + 1);
        adjacency[site + 1].push(site);
        edges.push(Bond::new(site, site + 1, BondType::Generic, 1.0));
    }
    (adjacency, CsrLattice::from_edges(n, edges))
}

fn cycle_parts(n: usize) -> (Vec<Vec<usize>>, CsrLattice) {
    let mut adjacency = vec![Vec::new(); n];
    for site in 0..n {
        let next = (site + 1) % n;
        adjacency[site].push(next);
        adjacency[next].push(site);
    }
    // Chain bonds plus the single wrap bond (a double edge for n = 2, a
    // self-loop for n = 1 — the multigraph C_n in each degenerate case).
    let mut edges: Vec<Bond> = (0..n.saturating_sub(1))
        .map(|site| Bond::new(site, site + 1, BondType::Generic, 1.0))
        .collect();
    if n > 0 {
        edges.push(Bond::new(n - 1, 0, BondType::Generic, 1.0));
    }
    (adjacency, CsrLattice::from_edges(n, edges))
}

fn star_parts(leaves: usize) -> (Vec<Vec<usize>>, CsrLattice) {
    let n = leaves + 1;
    let mut adjacency = vec![Vec::new(); n];
    let mut edges = Vec::new();
    for leaf in 1..=leaves {
        adjacency[0].push(leaf);
        adjacency[leaf].push(0);
        edges.push(Bond::new(0, leaf, BondType::Generic, 1.0));
    }
    (adjacency, CsrLattice::from_edges(n, edges))
}

fn complete_parts(n: usize) -> (Vec<Vec<usize>>, CsrLattice) {
    let mut adjacency = vec![Vec::new(); n];
    let mut edges = Vec::new();
    for left in 0..n {
        for right in (left + 1)..n {
            adjacency[left].push(right);
            adjacency[right].push(left);
            edges.push(Bond::new(left, right, BondType::Generic, 1.0));
        }
    }
    (adjacency, CsrLattice::from_edges(n, edges))
}

/// A seeded random simple graph plus a parallel-edge multigraph variant
/// for the order-independence gate.
fn random_graph_parts(
    vertices: usize,
    edge_probability: f64,
    extra_parallel: usize,
    seed: u64,
) -> (Vec<Vec<usize>>, CsrLattice) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut adjacency = vec![Vec::new(); vertices];
    let mut edges = Vec::new();
    for left in 0..vertices {
        for right in (left + 1)..vertices {
            if rng.random::<f64>() < edge_probability {
                for _ in 0..=extra_parallel {
                    adjacency[left].push(right);
                    adjacency[right].push(left);
                    edges.push(Bond::new(left, right, BondType::Generic, 1.0));
                }
            }
        }
    }
    (adjacency, CsrLattice::from_edges(vertices, edges))
}

/// Enumerate every seed mask as bools.
fn all_seed_masks(vertices: usize) -> Vec<Vec<bool>> {
    (0..1usize << vertices)
        .map(|mask| (0..vertices).map(|bit| mask & (1 << bit) != 0).collect())
        .collect()
}

fn seeds_from_mask(mask: &[bool]) -> MaskSeeds {
    MaskSeeds::new(&mask.iter().map(|&flag| u8::from(flag)).collect::<Vec<_>>())
}

/// The family graphs of the roadmap gate (paths, cycles, stars, complete
/// graphs), each as (reference adjacency, production lattice).
fn family_graphs() -> Vec<(&'static str, Vec<Vec<usize>>, CsrLattice)> {
    let (path_adjacency, path_lattice) = path_parts(5);
    let (cycle_adjacency, cycle_lattice) = cycle_parts(5);
    let (star_adjacency, star_lattice) = star_parts(4);
    let (complete4_adjacency, complete4_lattice) = complete_parts(4);
    let (complete5_adjacency, complete5_lattice) = complete_parts(5);
    vec![
        ("path5", path_adjacency, path_lattice),
        ("cycle5", cycle_adjacency, cycle_lattice),
        ("star4", star_adjacency, star_lattice),
        ("complete4", complete4_adjacency, complete4_lattice),
        ("complete5", complete5_adjacency, complete5_lattice),
    ]
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1: exhaustive 2^V initial states vs repeated full scans.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn exhaustive_initial_states_match_full_scan_on_family_graphs() {
    let cases = family_graphs();
    let mut workspace = BootstrapWorkspace::new();
    let mut trajectory = Vec::new();
    for (name, adjacency, lattice) in &cases {
        let vertices = adjacency.len();
        let max_degree = adjacency.iter().map(Vec::len).max().unwrap();
        // Homogeneous thresholds 0..=max_degree + 1 include both degenerate
        // gates (r = 0 and r > max_degree) in the same exhaustive sweep.
        let mut checked = 0usize;
        for r in 0..=(max_degree as u32 + 1) {
            let thresholds = vec![r; vertices];
            for mask in all_seed_masks(vertices) {
                let reference = full_scan_closure(adjacency, &mask, &thresholds);
                let seeds = seeds_from_mask(&mask);
                let rule = BootstrapThreshold::Homogeneous { value: r };
                let synchronous = workspace
                    .synchronous_rounds(lattice, &seeds, &rule, &mut trajectory)
                    .unwrap();
                assert_eq!(
                    active_set(&workspace, vertices),
                    reference,
                    "{name}: synchronous closure mismatch at r={r}, mask={mask:?}"
                );
                let queue = workspace
                    .queue_closure(lattice, &seeds, &rule, &mut trajectory)
                    .unwrap();
                assert_eq!(
                    active_set(&workspace, vertices),
                    reference,
                    "{name}: queue closure mismatch at r={r}, mask={mask:?}"
                );
                // Scalar bookkeeping invariants.
                let seed_count = mask.iter().filter(|&&flag| flag).count();
                let final_active = reference.iter().filter(|&&flag| flag).count();
                assert_eq!(synchronous.seed_count, seed_count);
                assert_eq!(queue.seed_count, seed_count);
                assert_eq!(synchronous.final_active, final_active);
                assert_eq!(queue.final_active, final_active);
                assert_eq!(
                    synchronous.activated_total,
                    final_active - seed_count,
                    "{name}: activated_total = final - seeds"
                );
                assert_eq!(
                    queue.activated_total,
                    final_active - seed_count,
                    "{name}: activated_total = final - seeds"
                );
                assert_eq!(queue.events, final_active - seed_count);
                checked += 1;
            }
        }
        assert!(checked >= 2usize.pow(vertices as u32), "{name}: exhausted");
    }
}

#[test]
fn exhaustive_synchronous_trajectories_match_simultaneous_sweeps() {
    let cases = family_graphs();
    let mut workspace = BootstrapWorkspace::new();
    let mut trajectory = Vec::new();
    for (name, adjacency, lattice) in &cases {
        let vertices = adjacency.len();
        let max_degree = adjacency.iter().map(Vec::len).max().unwrap();
        for r in 0..=(max_degree as u32 + 1) {
            let thresholds = vec![r; vertices];
            for mask in all_seed_masks(vertices) {
                let reference = full_scan_rounds(adjacency, &mask, &thresholds);
                let seeds = seeds_from_mask(&mask);
                let outcome = workspace
                    .synchronous_rounds(
                        lattice,
                        &seeds,
                        &BootstrapThreshold::Homogeneous { value: r },
                        &mut trajectory,
                    )
                    .unwrap();
                assert_eq!(
                    trajectory, reference,
                    "{name}: per-round trajectory mismatch at r={r}, mask={mask:?}"
                );
                assert_eq!(
                    outcome.rounds,
                    reference.len() - 1,
                    "{name}: rounds = sweeps that changed something"
                );
            }
        }
    }
}

#[test]
fn exhaustive_heterogeneous_thresholds_match_full_scan() {
    let cases = family_graphs();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00e2);
    let mut workspace = BootstrapWorkspace::new();
    let mut trajectory = Vec::new();
    for (name, adjacency, lattice) in &cases {
        let vertices = adjacency.len();
        // Three seeded threshold vectors per graph; values cover
        // 0..=degree+1 so zero thresholds and unreachable thresholds both
        // appear.
        for _ in 0..3 {
            let thresholds: Vec<u32> = (0..vertices)
                .map(|vertex| rng.random_range(0..=(adjacency[vertex].len() as u32 + 1)))
                .collect();
            for mask in all_seed_masks(vertices) {
                let reference = full_scan_closure(adjacency, &mask, &thresholds);
                let seeds = seeds_from_mask(&mask);
                let rule = BootstrapThreshold::Heterogeneous {
                    values: thresholds.clone(),
                };
                workspace
                    .synchronous_rounds(lattice, &seeds, &rule, &mut trajectory)
                    .unwrap();
                assert_eq!(
                    active_set(&workspace, vertices),
                    reference,
                    "{name}: heterogeneous synchronous mismatch, thresholds={thresholds:?}, mask={mask:?}"
                );
                workspace
                    .queue_closure(lattice, &seeds, &rule, &mut trajectory)
                    .unwrap();
                assert_eq!(
                    active_set(&workspace, vertices),
                    reference,
                    "{name}: heterogeneous queue mismatch, thresholds={thresholds:?}, mask={mask:?}"
                );
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gates 2 and 3: the degenerate-threshold semantics, stated verbatim.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn zero_threshold_activates_all_vertices() {
    let (path_adjacency, path_lattice) = path_parts(6);
    let (cycle_adjacency, cycle_lattice) = cycle_parts(6);
    let (star_adjacency, star_lattice) = star_parts(5);
    let (complete5_adjacency, complete5_lattice) = complete_parts(5);
    let cases: Vec<(&str, Vec<Vec<usize>>, CsrLattice)> = vec![
        ("path6", path_adjacency, path_lattice),
        ("cycle6", cycle_adjacency, cycle_lattice),
        ("star5", star_adjacency, star_lattice),
        ("complete5", complete5_adjacency, complete5_lattice),
    ];
    let mut workspace = BootstrapWorkspace::new();
    let mut trajectory = Vec::new();
    for (name, adjacency, lattice) in &cases {
        let vertices = adjacency.len();
        for mask in all_seed_masks(vertices) {
            let seeds = seeds_from_mask(&mask);
            let rule = BootstrapThreshold::Homogeneous { value: 0 };
            let synchronous = workspace
                .synchronous_rounds(lattice, &seeds, &rule, &mut trajectory)
                .unwrap();
            assert!(
                (0..vertices).all(|index| workspace.is_active(index)),
                "{name}: r = 0 must activate every vertex (mask={mask:?})"
            );
            assert_eq!(synchronous.final_active, vertices);
            // Zero thresholds qualify at the first update step.
            assert_eq!(synchronous.rounds, usize::from(mask.iter().any(|&f| !f)));
            let queue = workspace
                .queue_closure(lattice, &seeds, &rule, &mut trajectory)
                .unwrap();
            assert_eq!(queue.final_active, vertices);
            assert!(
                (0..vertices).all(|index| workspace.is_active(index)),
                "{name}: r = 0 queue closure must activate every vertex"
            );
        }
    }
}

#[test]
fn threshold_above_max_degree_preserves_only_seeds() {
    let (path_adjacency, path_lattice) = path_parts(6);
    let (cycle_adjacency, cycle_lattice) = cycle_parts(6);
    let (star_adjacency, star_lattice) = star_parts(5);
    let (complete5_adjacency, complete5_lattice) = complete_parts(5);
    let cases: Vec<(&str, Vec<Vec<usize>>, CsrLattice)> = vec![
        ("path6", path_adjacency, path_lattice),
        ("cycle6", cycle_adjacency, cycle_lattice),
        ("star5", star_adjacency, star_lattice),
        ("complete5", complete5_adjacency, complete5_lattice),
    ];
    let mut workspace = BootstrapWorkspace::new();
    let mut trajectory = Vec::new();
    for (name, adjacency, lattice) in &cases {
        let vertices = adjacency.len();
        let max_degree = adjacency.iter().map(Vec::len).max().unwrap() as u32;
        for mask in all_seed_masks(vertices) {
            let expected: Vec<bool> = mask.clone();
            let seeds = seeds_from_mask(&mask);
            let rule = BootstrapThreshold::Homogeneous {
                value: max_degree + 1,
            };
            let synchronous = workspace
                .synchronous_rounds(lattice, &seeds, &rule, &mut trajectory)
                .unwrap();
            assert_eq!(
                active_set(&workspace, vertices),
                expected,
                "{name}: r > max_degree must preserve only the seeds"
            );
            assert_eq!(synchronous.rounds, 0);
            assert_eq!(synchronous.activated_total, 0);
            workspace
                .queue_closure(lattice, &seeds, &rule, &mut trajectory)
                .unwrap();
            assert_eq!(
                active_set(&workspace, vertices),
                expected,
                "{name}: r > max_degree must preserve only the seeds (queue)"
            );
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 4: order independence of the final closure.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn queue_closure_fixed_point_is_order_independent() {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004f);
    let graphs = [
        random_graph_parts(8, 0.3, 0, 0x0051),
        random_graph_parts(8, 0.5, 0, 0x0052),
        random_graph_parts(9, 0.35, 1, 0x0053),
        path_parts(7),
        cycle_parts(7),
        star_parts(6),
        complete_parts(5),
    ];
    let mut workspace = BootstrapWorkspace::new();
    let mut trajectory = Vec::new();
    for (case, (adjacency, lattice)) in graphs.iter().enumerate() {
        let vertices = adjacency.len();
        for _ in 0..24 {
            // Random seeds and random heterogeneous thresholds.
            let mask: Vec<bool> = (0..vertices).map(|_| rng.random::<f64>() < 0.4).collect();
            let thresholds: Vec<u32> = (0..vertices)
                .map(|vertex| rng.random_range(0..=(adjacency[vertex].len() as u32 + 1)))
                .collect();
            let reference = full_scan_closure(adjacency, &mask, &thresholds);
            let one_at_a_time = one_at_a_time_closure(adjacency, &mask, &thresholds, &mut rng);
            let shuffled_sweep = shuffled_sweep_closure(adjacency, &mask, &thresholds, &mut rng);
            assert_eq!(one_at_a_time, reference, "case {case}: randomized order");
            assert_eq!(shuffled_sweep, reference, "case {case}: shuffled sweeps");
            let seeds = seeds_from_mask(&mask);
            let rule = BootstrapThreshold::Heterogeneous {
                values: thresholds.clone(),
            };
            let queue = workspace
                .queue_closure(lattice, &seeds, &rule, &mut trajectory)
                .unwrap();
            assert_eq!(
                active_set(&workspace, vertices),
                reference,
                "case {case}: FIFO queue closure = order-independent fixed point"
            );
            assert_eq!(queue.events, queue.activated_total);
            let synchronous = workspace
                .synchronous_rounds(lattice, &seeds, &rule, &mut trajectory)
                .unwrap();
            assert_eq!(
                active_set(&workspace, vertices),
                reference,
                "case {case}: synchronous and queue closures agree"
            );
            assert_eq!(synchronous.final_active, queue.final_active);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 5: hand-worked heterogeneous example.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn heterogeneous_hand_worked_example() {
    // Path 0-1-2-3-4, thresholds [2,1,1,2,3], seeds {0,4}:
    //   round 1: v1 (1 >= 1);
    //   round 2: v2 (1 >= 1 via v1);
    //   round 3: v3 (2 >= 2 via v2 and v4);
    // closure = all five vertices, trajectory [2,3,4,5].
    let (adjacency, lattice) = path_parts(5);
    let thresholds = vec![2, 1, 1, 2, 3];
    let mask = vec![true, false, false, false, true];
    assert_eq!(
        full_scan_rounds(&adjacency, &mask, &thresholds),
        vec![2, 3, 4, 5],
        "reference must reproduce the hand derivation"
    );
    let mut workspace = BootstrapWorkspace::new();
    let seeds = seeds_from_mask(&mask);
    let rule = BootstrapThreshold::Heterogeneous {
        values: thresholds.clone(),
    };
    let mut trajectory = Vec::new();
    let synchronous = workspace
        .synchronous_rounds(&lattice, &seeds, &rule, &mut trajectory)
        .unwrap();
    assert_eq!(trajectory, vec![2, 3, 4, 5]);
    assert_eq!(synchronous.rounds, 3);
    assert_eq!(synchronous.seed_count, 2);
    assert_eq!(synchronous.activated_total, 3);
    let queue = workspace
        .queue_closure(&lattice, &seeds, &rule, &mut trajectory)
        .unwrap();
    assert_eq!(queue.events, 3);
    assert_eq!(queue.final_active, 5);
    assert_eq!(active_set(&workspace, 5), vec![true; 5]);
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 6: adapter contracts (seeds, draw contract, C0 recording).
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn adapter_fixed_seeds_match_core_bitwise() {
    let lattice = build_square(6, 6, true);
    let vertices = lattice.vertex_count();
    let rule = BootstrapThreshold::Homogeneous { value: 2 };
    let seeds = BootstrapSeeds::Fixed {
        vertices: vec![0, 7, 20, 35],
    };
    for dynamics in [
        BootstrapDynamics::SynchronousRounds,
        BootstrapDynamics::QueueClosure,
    ] {
        let mut adapter =
            BootstrapPercolationMC::new(lattice.clone(), rule.clone(), seeds.clone(), dynamics, 8)
                .unwrap();
        let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(1), 0, 1 << 16);
        adapter.sweep(&mut ctx);
        adapter.measure(&mut ctx);
        // Core replay with the same explicit seeds.
        let mask: Vec<u8> = (0..vertices)
            .map(|index| u8::from([0, 7, 20, 35].contains(&index)))
            .collect();
        let core_seeds = MaskSeeds::new(&mask);
        let mut workspace = BootstrapWorkspace::new();
        let mut trajectory = Vec::new();
        let outcome = match dynamics {
            BootstrapDynamics::SynchronousRounds => {
                workspace.synchronous_rounds(&lattice, &core_seeds, &rule, &mut trajectory)
            }
            BootstrapDynamics::QueueClosure => {
                workspace.queue_closure(&lattice, &core_seeds, &rule, &mut trajectory)
            }
        }
        .unwrap();
        assert_eq!(adapter.last_outcome(), outcome);
        assert_eq!(adapter.last_trajectory(), trajectory.as_slice());
        // Fixed history: the recorded array pads the trajectory with the
        // final count up to the fixed shape.
        let final_active = *trajectory.last().unwrap() as f64;
        let expected: Vec<f64> = (0..8)
            .map(|step| {
                if step < trajectory.len() {
                    trajectory[step] as f64
                } else {
                    final_active
                }
            })
            .collect();
        assert_eq!(adapter.last_history(), expected.as_slice());
    }
}

#[test]
fn adapter_density_draw_contract_is_bitwise_reproducible() {
    let lattice = build_square(8, 8, true);
    let vertices = lattice.vertex_count();
    let density = 0.37;
    let seed = 0x0b00;
    let rule = BootstrapThreshold::Homogeneous { value: 2 };
    let mut adapter = BootstrapPercolationMC::new(
        lattice.clone(),
        rule.clone(),
        BootstrapSeeds::Density { p: density },
        BootstrapDynamics::QueueClosure,
        16,
    )
    .unwrap();
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1 << 16);
    adapter.sweep(&mut ctx);

    // Independent replay of the documented contract: one f64 draw per
    // vertex, ascending dense index, seed iff draw < p.
    let mut replay_rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mask: Vec<u8> = (0..vertices)
        .map(|_| u8::from(replay_rng.random::<f64>() < density))
        .collect();
    let replay_seeds = MaskSeeds::new(&mask);
    let mut workspace = BootstrapWorkspace::new();
    let mut trajectory = Vec::new();
    let outcome = workspace
        .queue_closure(&lattice, &replay_seeds, &rule, &mut trajectory)
        .unwrap();
    assert_eq!(adapter.last_outcome(), outcome);
    assert_eq!(adapter.last_trajectory(), trajectory.as_slice());
}

#[test]
fn adapter_records_one_vector_sample_plus_closure_scalars() {
    let lattice = build_square(6, 6, true);
    let steps = 10;
    let mut adapter = BootstrapPercolationMC::new(
        lattice,
        BootstrapThreshold::Homogeneous { value: 2 },
        BootstrapSeeds::Fixed {
            vertices: vec![0, 35],
        },
        BootstrapDynamics::SynchronousRounds,
        steps,
    )
    .unwrap();
    let sweeps = 16usize;
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(7), 0, 64);
    for _ in 0..sweeps {
        adapter.sweep(&mut ctx);
        adapter.measure(&mut ctx);
    }
    let (scalars, arrays) = ctx.finalize_measurements_split();
    let scalar: Option<&Estimate> = scalars.get("BootstrapFinalActive");
    let activated: Option<&Estimate> = scalars.get("BootstrapActivatedTotal");
    let rounds: Option<&Estimate> = scalars.get("BootstrapRounds");
    assert!(!scalars.contains_key("BootstrapEvents"));
    for estimate in [scalar, activated, rounds].into_iter().flatten() {
        assert_eq!(estimate.n_bins, 1, "binsize 64, {sweeps} samples: one bin");
    }
    let history = arrays
        .get("BootstrapActiveCount")
        .expect("BootstrapActiveCount array recorded");
    assert_eq!(history.shape, vec![steps]);
    assert_eq!(history.elements.len(), steps);
    // Deterministic seeds: every element mean equals the single realized
    // trajectory exactly, with zero error.
    let recorded = adapter.last_history();
    for (element, value) in history.elements.iter().zip(recorded) {
        assert_eq!(element.mean, *value);
        assert_eq!(element.stderr, 0.0);
    }
    assert_eq!(
        scalar.unwrap().mean as usize,
        recorded.last().copied().unwrap() as usize
    );
    // The closure invariant: activated = final - seeds (2 fixed seeds).
    assert_eq!(
        activated.unwrap().mean as usize + 2,
        scalar.unwrap().mean as usize
    );
    assert_eq!(rounds.unwrap().mean, adapter.last_outcome().rounds as f64);
}

#[test]
fn adapter_rejects_invalid_domains() {
    let lattice = build_chain(4, false);
    let invalid_density = BootstrapPercolationMC::new(
        lattice.clone(),
        BootstrapThreshold::Homogeneous { value: 2 },
        BootstrapSeeds::Density { p: 1.5 },
        BootstrapDynamics::SynchronousRounds,
        4,
    );
    assert!(matches!(
        invalid_density,
        Err(BootstrapError::Probability(_))
    ));
    let short_thresholds = BootstrapPercolationMC::new(
        lattice,
        BootstrapThreshold::Heterogeneous {
            values: vec![2, 2, 2],
        },
        BootstrapSeeds::Density { p: 0.5 },
        BootstrapDynamics::SynchronousRounds,
        4,
    );
    assert!(matches!(
        short_thresholds,
        Err(BootstrapError::ThresholdLengthMismatch {
            expected: 4,
            actual: 3
        })
    ));
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 7: Holroyd (2003) finite-size metastability anchor.
// ─────────────────────────────────────────────────────────────────────────

/// Fraction of realizations whose closure activates the majority of the
/// L x L torus, with a binomial standard error.
fn majority_probability(width: usize, density: f64, realizations: usize, seed: u64) -> (f64, f64) {
    let vertices = width * width;
    let mut adapter = BootstrapPercolationMC::new(
        build_square(width, width, true),
        BootstrapThreshold::Homogeneous { value: 2 },
        BootstrapSeeds::Density { p: density },
        BootstrapDynamics::SynchronousRounds,
        8,
    )
    .unwrap();
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1 << 20);
    let mut hits = 0usize;
    for _ in 0..realizations {
        adapter.sweep(&mut ctx);
        if adapter.last_outcome().final_active * 2 > vertices {
            hits += 1;
        }
    }
    let mean = hits as f64 / realizations as f64;
    let stderr = (mean * (1.0 - mean) / realizations as f64).sqrt();
    (mean, stderr)
}

#[test]
fn holroyd_sharp_metastability_bracket() {
    // Holroyd, PTRF 125, 195 (2003), Theorem 1: I(L,p) -> 1 iff
    // p ln L > lambda = pi^2/18 — the finite-size seed threshold decays
    // like 1/ln L. His convergence is notoriously slow (finite sizes sit
    // far below the asymptotic constant; the same effect that made the
    // pre-Holroyd numerics of Adler-Stauffer-Aharony misleading), so the
    // gates below are statistical brackets, never an exact threshold
    // claim: the measured half-majority density p_hat(L) must decay
    // monotonically in L (the 1/ln L decay), keep p_hat(L) ln L strictly
    // inside (lambda/4, lambda) (approaching the asymptotic constant from
    // below), and the metastable window must be sharper than a factor 1.5
    // in p at every measured size.
    let lambda = std::f64::consts::PI * std::f64::consts::PI / 18.0;
    let realizations = anchor_realizations();
    let mut p_hat = [0.0f64; 3];
    for (index, width) in [32usize, 64, 128].iter().enumerate() {
        let width = *width;
        // Bisect for the half-majority density; stop once the measured
        // probability is statistically compatible with 1/2.
        let mut low = 0.2 * lambda / (width as f64).ln();
        let mut high = 0.6 * lambda / (width as f64).ln();
        let mut converged = false;
        for _ in 0..9 {
            let mid = 0.5 * (low + high);
            let (mean, stderr) =
                majority_probability(width, mid, realizations, 0x0b00 + width as u64);
            println!(
                "HOLROYD L={width} bisect p={mid:.5} P(majority) = {mean:.4} +/- {stderr:.4}                  (p ln L / lambda = {:.3})",
                mid * (width as f64).ln() / lambda
            );
            if mean + 2.0 * stderr < 0.5 {
                low = mid;
            } else if mean - 2.0 * stderr > 0.5 {
                high = mid;
            } else {
                low = mid;
                high = mid;
                converged = true;
                break;
            }
        }
        assert!(converged, "L={width}: bisection must localize the crossing");
        let crossing = 0.5 * (low + high);
        p_hat[index] = crossing;
        // Sharpness of the metastable window at this size.
        let (low_win, low_err) =
            majority_probability(width, crossing / 1.5, realizations, 0x0b01 + width as u64);
        let (high_win, high_err) =
            majority_probability(width, 1.5 * crossing, realizations, 0x0b02 + width as u64);
        println!(
            "HOLROYD L={width} p_hat = {crossing:.5} +/- ~{:.5} (p_hat ln L = {:.4}, lambda = {lambda:.4});              window: P(p_hat/1.5) = {low_win:.4} +/- {low_err:.4}, P(1.5 p_hat) = {high_win:.4} +/- {high_err:.4}",
            0.5 * (high - low),
            crossing * (width as f64).ln()
        );
        let rescaled = crossing * (width as f64).ln();
        assert!(
            rescaled > lambda / 4.0 && rescaled < lambda,
            "L={width}: p_hat ln L = {rescaled:.4} must lie strictly inside (lambda/4, lambda) = ({:.4}, {lambda:.4})",
            lambda / 4.0
        );
        assert!(
            low_win + 3.0 * low_err < 0.1,
            "L={width}: majority activation must be rare at p_hat/1.5"
        );
        assert!(
            high_win - 3.0 * high_err > 0.9,
            "L={width}: majority activation must be near-certain at 1.5 p_hat"
        );
    }
    // The 1/ln L decay: monotonically decreasing threshold density, with
    // the extremes separated beyond the bisection uncertainty.
    assert!(
        p_hat[0] > p_hat[1] && p_hat[1] > p_hat[2],
        "p_hat must decay with L (measured {:?})",
        p_hat
    );
    let separation = 0.6 * lambda / 32f64.ln() - 0.2 * lambda / 128f64.ln();
    assert!(
        p_hat[0] - p_hat[2] > 0.1 * separation,
        "p_hat(32) - p_hat(128) = {:.5} must be a clear separation",
        p_hat[0] - p_hat[2]
    );
}

#[test]
fn holroyd_full_activation_bracket() {
    // The paper's own observable I(L,p): the probability that the whole
    // L x L torus activates. Statistical brackets at formula-derived
    // densities: near-certain at the asymptotic threshold
    // p+ = lambda / ln L (the p ln L > lambda side), rare at
    // p- = lambda / (4 ln L) (the p ln L < lambda side). Measured values
    // recorded with real error bars; no exact finite-L threshold claim.
    let lambda = std::f64::consts::PI * std::f64::consts::PI / 18.0;
    let realizations = anchor_realizations();
    for width in [32usize, 64, 128] {
        let vertices = width * width;
        let log_width = (width as f64).ln();
        for (label, density, seed) in [
            ("below", lambda / (4.0 * log_width), 0x0f01 + width as u64),
            ("above", lambda / log_width, 0x0f02 + width as u64),
        ] {
            let mut adapter = BootstrapPercolationMC::new(
                build_square(width, width, true),
                BootstrapThreshold::Homogeneous { value: 2 },
                BootstrapSeeds::Density { p: density },
                BootstrapDynamics::QueueClosure,
                8,
            )
            .unwrap();
            let mut ctx =
                Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1 << 20);
            let mut hits = 0usize;
            for _ in 0..realizations {
                adapter.sweep(&mut ctx);
                if adapter.last_outcome().final_active == vertices {
                    hits += 1;
                }
            }
            let mean = hits as f64 / realizations as f64;
            let stderr = (mean * (1.0 - mean) / realizations as f64)
                .sqrt()
                .max(0.5 / realizations as f64);
            println!(
                "HOLROYD_FULL L={width} {label} p={density:.5} (p ln L / lambda = {:.3}) \
                 P(full activation) = {mean:.4} +/- {stderr:.4}",
                density * log_width / lambda
            );
            match label {
                "below" => assert!(
                    mean + 3.0 * stderr < 0.05,
                    "L={width}: full activation must be rare at lambda/(4 ln L)"
                ),
                _ => assert!(
                    mean - 3.0 * stderr > 0.95,
                    "L={width}: full activation must be near-certain at lambda/ln L"
                ),
            }
        }
    }
}
