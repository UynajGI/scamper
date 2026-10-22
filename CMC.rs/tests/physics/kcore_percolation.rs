//! E3 production validation: k-core percolation.
//!
//! Gates (roadmap §12 E3), each anchored independently:
//!
//! 1. **Exhaustive initial subgraphs vs repeated full scans**: every `2^V`
//!    initial activity subset on path/cycle/star/complete graphs, uniform
//!    `k in {0,1,2,3,6}` plus heterogeneous threshold patterns, against an
//!    independently written repeated-full-scan reference (adjacency lists,
//!    simultaneous-deficiency removal per scan — also the independent
//!    check of the synchronous-round accounting).
//! 2. **`k = 0` removes nothing** (the core is exactly the initial
//!    activity) and **`k > max_degree` removes everything** in one
//!    synchronous round.
//! 3. **Queue-order independence**: random vertex relabelings (with
//!    shuffled edge insertion order) genuinely shuffle the deletion/queue
//!    order; the final core set and the complete summary stay identical —
//!    the mathematical uniqueness of the maximal min-residual-degree
//!    subgraph, verified empirically.
//! 4. **Hand-worked heterogeneous thresholds** with full wave accounting.
//! 5. **Classic deterministic facts**: the 1-core of any connected graph
//!    (no isolates) is the whole graph; the k-core of a tree is empty for
//!    `k >= 2`.
//! 6. **Adapter contract**: `KCoreMC` at `keep = 1` is bitwise the direct
//!    full-activity decomposition; its initial activity is bitwise the F2
//!    site law on the same stream; one realization records exactly one
//!    vector sample; the recorded-round budget truncates/pads without
//!    changing the exact scalars.
//! 7. **Literature anchor** — Janson & Luczak, Random Structures
//!    Algorithms 30, 50-62 (2007) (arXiv:math/0508453, archived in
//!    `docs/ref/`), Theorem 1.1 (after Pittel, Spencer & Wormald): on
//!    G(n, m) with mean degree `lambda = 2m/n`, the k-core (k >= 3) is
//!    empty whp below `lambda_k = min_{mu>0} mu/psi_{k-1}(mu)` and giant
//!    above it with `v(Core_k)/n -> psi_k(mu_k(lambda))`, where
//!    `psi_j(mu) = P(Poisson(mu) >= j)` and `mu_k(lambda)` is the largest
//!    solution of `mu = lambda psi_{k-1}(mu)`. An own solver computes the
//!    fixed-point system; G(n = 5e4, m) realizations bracket the 3-core
//!    emergence statistically (real error bars, never an exact threshold
//!    claim).

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    kcore_decompose, Bond, BondType, CoreThresholds, CsrLattice, GraphView, KCoreMC, KCoreSummary,
    KCoreWorkspace, ObservablePlan, StaticLaw, StaticObservable, StaticPercolationMC,
    VertexActivity,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

use super::common::{erdos_renyi_graph, star_lattice};

fn path_lattice(vertices: usize) -> CsrLattice {
    CsrLattice::from_edges(
        vertices,
        (0..vertices - 1)
            .map(|site| Bond::new(site, site + 1, BondType::Generic, 1.0))
            .collect(),
    )
}

fn cycle_lattice(vertices: usize) -> CsrLattice {
    CsrLattice::from_edges(
        vertices,
        (0..vertices)
            .map(|site| Bond::new(site, (site + 1) % vertices, BondType::Generic, 1.0))
            .collect(),
    )
}

fn complete_lattice(vertices: usize) -> CsrLattice {
    let mut edges = Vec::new();
    for left in 0..vertices {
        for right in (left + 1)..vertices {
            edges.push(Bond::new(left, right, BondType::Generic, 1.0));
        }
    }
    CsrLattice::from_edges(vertices, edges)
}

fn adjacency_of(lattice: &CsrLattice) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); lattice.n_sites];
    for bond in &lattice.edges {
        adjacency[bond.source].push(bond.target);
        adjacency[bond.target].push(bond.source);
    }
    adjacency
}

fn subsets(vertices: usize) -> impl Iterator<Item = Vec<bool>> {
    (0..1usize << vertices).map(move |mask| (0..vertices).map(|v| mask & (1 << v) != 0).collect())
}

// ─────────────────────────────────────────────────────────────────────────
// Independent reference: repeated full scans (no E1/E3 production code).
// ─────────────────────────────────────────────────────────────────────────

/// Repeated full-scan peeling: one scan finds every deficient vertex
/// against the current alive set; all of them are removed simultaneously
/// (one synchronous round); repeat to the fixed point. Returns the core
/// membership and the same summary fields the production path reports.
fn reference_kcore(
    adjacency: &[Vec<usize>],
    edges: &[(usize, usize)],
    thresholds: &[u32],
    active: &[bool],
) -> (Vec<bool>, KCoreSummary) {
    let vertices = adjacency.len();
    let mut alive = active.to_vec();
    let mut rounds = 0usize;
    loop {
        let mut deficient = Vec::new();
        for vertex in 0..vertices {
            if !alive[vertex] {
                continue;
            }
            let residual = adjacency[vertex]
                .iter()
                .filter(|&&neighbor| alive[neighbor])
                .count() as u32;
            if residual < thresholds[vertex] {
                deficient.push(vertex);
            }
        }
        if deficient.is_empty() {
            break;
        }
        rounds += 1;
        for vertex in deficient {
            alive[vertex] = false;
        }
    }
    let core_vertex_count = alive.iter().filter(|&&a| a).count();
    let core_edge_count = edges
        .iter()
        .filter(|&&(left, right)| alive[left] && alive[right])
        .count();
    // Largest core component by flood fill.
    let mut visited = vec![false; vertices];
    let mut largest = 0usize;
    for seed in 0..vertices {
        if !alive[seed] || visited[seed] {
            continue;
        }
        let mut stack = vec![seed];
        visited[seed] = true;
        let mut size = 0usize;
        while let Some(vertex) = stack.pop() {
            size += 1;
            for &neighbor in &adjacency[vertex] {
                if alive[neighbor] && !visited[neighbor] {
                    visited[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        largest = largest.max(size);
    }
    let initial: usize = active.iter().filter(|&&a| a).count();
    (
        alive,
        KCoreSummary {
            core_vertex_count,
            core_edge_count,
            largest_core_size: largest,
            peeling_rounds: rounds,
            deleted_total: initial - core_vertex_count,
        },
    )
}

fn edge_pairs(lattice: &CsrLattice) -> Vec<(usize, usize)> {
    lattice
        .edges
        .iter()
        .map(|bond| (bond.source, bond.target))
        .collect()
}

fn uniform_thresholds(k: usize) -> Vec<u32> {
    vec![k as u32; 16]
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 1: exhaustive initial subgraphs vs repeated full scans.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn exhaustive_initial_subgraphs_match_repeated_full_scan() {
    struct Case {
        name: &'static str,
        lattice: CsrLattice,
    }
    let cases = [
        Case {
            name: "path5",
            lattice: path_lattice(5),
        },
        Case {
            name: "cycle6",
            lattice: cycle_lattice(6),
        },
        Case {
            name: "star6",
            lattice: star_lattice(5),
        },
        Case {
            name: "complete4",
            lattice: complete_lattice(4),
        },
        Case {
            name: "complete5",
            lattice: complete_lattice(5),
        },
    ];

    let mut checked = 0usize;
    let mut workspace = KCoreWorkspace::default();
    for case in &cases {
        let vertices = case.lattice.n_sites;
        let adjacency = adjacency_of(&case.lattice);
        let edges = edge_pairs(&case.lattice);
        // Uniform k in {0, 1, 2, 3, 6} (6 > every max degree here) plus
        // two heterogeneous threshold patterns.
        let mut threshold_sets: Vec<(String, Vec<u32>)> = [0usize, 1, 2, 3, 6]
            .iter()
            .map(|&k| (format!("k={k}"), uniform_thresholds(k)))
            .collect();
        threshold_sets.push((
            "heterogeneous-period3".into(),
            (0..16).map(|v| (v % 3 + 1) as u32).collect(),
        ));
        threshold_sets.push((
            "heterogeneous-capped".into(),
            (0..16).map(|v| v.min(3) as u32).collect(),
        ));

        for (label, thresholds) in &threshold_sets {
            let core_thresholds = CoreThresholds::PerVertex(thresholds[..vertices].to_vec());
            for active in subsets(vertices) {
                let summary =
                    kcore_decompose(&case.lattice, &core_thresholds, &active, &mut workspace)
                        .unwrap_or_else(|error| panic!("{}/{}: {error}", case.name, label));
                let (reference_alive, reference) =
                    reference_kcore(&adjacency, &edges, &thresholds[..vertices], &active);
                for (vertex, &alive) in reference_alive.iter().enumerate() {
                    assert_eq!(
                        workspace.in_core(vertex),
                        alive,
                        "{}/{} subset {:?} vertex {vertex}: core membership diverged",
                        case.name,
                        label,
                        active
                    );
                }
                assert_eq!(
                    summary, reference,
                    "{}/{} subset {:?}: summary diverged",
                    case.name, label, active
                );
                assert_eq!(
                    workspace.round_survivors().len(),
                    summary.peeling_rounds + 1,
                    "survivor curve length must be rounds + 1"
                );
                checked += 1;
            }
        }
        eprintln!(
            "[E3 exhaustive] {}: all initial subgraphs match the reference",
            case.name
        );
    }
    eprintln!("[E3 exhaustive] {checked} (graph, threshold, subset) runs matched");
    assert!(checked > 1000, "the exhaustive corpus must stay sizable");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 2: k = 0 removes nothing; k > max_degree removes everything.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn k_zero_removes_nothing_and_large_k_removes_everything() {
    let lattices = [
        ("path7", path_lattice(7), 2),
        ("cycle6", cycle_lattice(6), 2),
        ("star6", star_lattice(5), 5),
        ("complete5", complete_lattice(5), 4),
    ];
    let mut workspace = KCoreWorkspace::default();
    let mut subsets_checked = 0;
    for (name, lattice, max_degree) in &lattices {
        let vertices = lattice.n_sites;
        // k = 0: the core is exactly the initial activity, no rounds.
        for active in subsets(vertices) {
            let summary = kcore_decompose(
                lattice,
                &CoreThresholds::uniform(0).unwrap(),
                &active,
                &mut workspace,
            )
            .unwrap();
            let initial = active.iter().filter(|&&a| a).count();
            assert_eq!(
                summary.core_vertex_count, initial,
                "{name}/k=0: nothing removed"
            );
            assert_eq!(summary.peeling_rounds, 0, "{name}/k=0: no peeling rounds");
            assert_eq!(summary.deleted_total, 0);
            assert!(workspace.deletion_order().is_empty());
            for (vertex, &is_active) in active.iter().enumerate() {
                assert_eq!(workspace.in_core(vertex), is_active);
            }
            subsets_checked += 1;
        }

        // k > max degree: everything is initially deficient; one
        // synchronous round removes the whole active graph.
        let full = vec![true; vertices];
        let summary = kcore_decompose(
            lattice,
            &CoreThresholds::uniform(max_degree + 1).unwrap(),
            &full,
            &mut workspace,
        )
        .unwrap();
        assert_eq!(summary.core_vertex_count, 0, "{name}/k>max: empty core");
        assert_eq!(summary.peeling_rounds, 1, "{name}/k>max: single round");
        assert_eq!(summary.deleted_total, vertices);
        assert_eq!(summary.core_edge_count, 0);
        assert_eq!(
            workspace.deletion_order(),
            &(0..vertices).collect::<Vec<_>>(),
            "{name}/k>max: ascending deletion order"
        );
        assert_eq!(workspace.round_survivors(), &[vertices as f64, 0.0]);
    }
    eprintln!("[E3 endpoints] k=0 verified on {subsets_checked} subsets; k>max_degree verified on 4 graphs");
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 3: queue-order independence (shuffled relabeling + edge order).
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn final_core_is_independent_of_queue_order() {
    // Relabeled copy of `base` under `permutation`, with a shuffled edge
    // insertion order so both the initial scan order and the CSR incidence
    // order (hence the queue arrival order) change.
    let relabeled = |base: &CsrLattice, permutation: &[usize], rng: &mut Xoshiro256PlusPlus| {
        let mut bonds: Vec<_> = base
            .edges
            .iter()
            .map(|bond| {
                Bond::new(
                    permutation[bond.source],
                    permutation[bond.target],
                    BondType::Generic,
                    1.0,
                )
            })
            .collect();
        for position in (1..bonds.len()).rev() {
            let pick = rng.random_range(0..=position);
            bonds.swap(position, pick);
        }
        CsrLattice::from_edges(base.n_sites, bonds)
    };

    let lattices = [
        ("path8", path_lattice(8)),
        ("cycle7", cycle_lattice(7)),
        ("star7", star_lattice(6)),
        ("complete5", complete_lattice(5)),
    ];
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xE30C_DE0F);
    let mut workspace = KCoreWorkspace::default();
    let mut distinct_orders: Vec<Vec<usize>> = Vec::new();

    for (name, lattice) in &lattices {
        let vertices = lattice.n_sites;
        // A handful of seeded initial subsets (including the full graph).
        let mut subsets_by_k: Vec<Vec<bool>> = vec![vec![true; vertices]];
        for _ in 0..2 {
            subsets_by_k.push((0..vertices).map(|_| rng.random::<bool>()).collect());
        }
        for k in [1usize, 2, 3] {
            for active in &subsets_by_k {
                let base_summary = kcore_decompose(
                    lattice,
                    &CoreThresholds::uniform(k).unwrap(),
                    active,
                    &mut workspace,
                )
                .unwrap();
                let base_core: Vec<bool> = (0..vertices)
                    .map(|vertex| workspace.in_core(vertex))
                    .collect();

                for _ in 0..12 {
                    let mut permutation: Vec<usize> = (0..vertices).collect();
                    for position in (1..vertices).rev() {
                        let pick = rng.random_range(0..=position);
                        permutation.swap(position, pick);
                    }
                    let shuffled = relabeled(lattice, &permutation, &mut rng);
                    // Original vertex v lives at relabeled index
                    // permutation[v]; carry the activity along the same
                    // mapping.
                    let mut relabeled_active = vec![false; vertices];
                    for (vertex, &target) in permutation.iter().enumerate() {
                        relabeled_active[target] = active[vertex];
                    }
                    let summary = kcore_decompose(
                        &shuffled,
                        &CoreThresholds::uniform(k).unwrap(),
                        &relabeled_active,
                        &mut workspace,
                    )
                    .unwrap();
                    assert_eq!(
                        summary, base_summary,
                        "{name}/k={k}: summary changed under relabeling"
                    );
                    for vertex in 0..vertices {
                        assert_eq!(
                            workspace.in_core(permutation[vertex]),
                            base_core[vertex],
                            "{name}/k={k}: core membership changed under relabeling"
                        );
                    }
                    // Map the deletion order back to original labels and
                    // register it as an order variant.
                    let mapped: Vec<usize> = workspace
                        .deletion_order()
                        .iter()
                        .map(|&v| permutation.iter().position(|&p| p == v).unwrap())
                        .collect();
                    if !distinct_orders.contains(&mapped) {
                        distinct_orders.push(mapped);
                    }
                }
            }
        }
    }
    eprintln!(
        "[E3 order] {} distinct deletion orders observed; all cores identical",
        distinct_orders.len()
    );
    assert!(
        distinct_orders.len() >= 2,
        "the relabelings must genuinely shuffle the queue order"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 4: hand-worked heterogeneous thresholds.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn heterogeneous_thresholds_hand_worked_examples() {
    let mut workspace = KCoreWorkspace::default();

    // Path 0-1-2-3-4, thresholds [2,3,1,3,2]: initially deficient are the
    // leaves (1 < 2) AND vertices 1, 3 (2 < 3); round 1 deletes all four,
    // vertex 2 drops to residual 0 < 1 and dies in round 2. Survivors
    // 5 -> 1 -> 0.
    let path5 = path_lattice(5);
    let summary = kcore_decompose(
        &path5,
        &CoreThresholds::PerVertex(vec![2, 3, 1, 3, 2]),
        &[true; 5],
        &mut workspace,
    )
    .unwrap();
    assert_eq!(summary.peeling_rounds, 2);
    assert_eq!(summary.deleted_total, 5);
    assert_eq!(workspace.round_survivors(), &[5.0, 1.0, 0.0]);
    assert_eq!(workspace.deletion_round(0), Some(1));
    assert_eq!(workspace.deletion_round(1), Some(1));
    assert_eq!(workspace.deletion_round(2), Some(2));
    assert_eq!(workspace.deletion_round(4), Some(1));

    // Star (center 0, five leaves), thresholds: center 6, leaves 1. The
    // center is deficient (5 < 6) and dies alone in round 1; the leaves
    // drop to residual 0 < 1 and die together in round 2. Survivors
    // 6 -> 5 -> 0.
    let star6 = star_lattice(5);
    let summary = kcore_decompose(
        &star6,
        &CoreThresholds::PerVertex(vec![6, 1, 1, 1, 1, 1]),
        &[true; 6],
        &mut workspace,
    )
    .unwrap();
    assert_eq!(summary.peeling_rounds, 2);
    assert_eq!(workspace.round_survivors(), &[6.0, 5.0, 0.0]);
    assert_eq!(workspace.deletion_round(0), Some(1));
    assert_eq!(workspace.deletion_round(1), Some(2));

    // Path 0-1-2-3, thresholds [2,2,1,2]: leaves die in round 1, vertex 1
    // in round 2 (1 < 2), vertex 2 in round 3 (0 < 1).
    let path4 = path_lattice(4);
    let summary = kcore_decompose(
        &path4,
        &CoreThresholds::PerVertex(vec![2, 2, 1, 2]),
        &[true; 4],
        &mut workspace,
    )
    .unwrap();
    assert_eq!(summary.peeling_rounds, 3);
    assert_eq!(workspace.round_survivors(), &[4.0, 2.0, 1.0, 0.0]);
    assert_eq!(workspace.deletion_round(2), Some(3));

    // A vertex whose threshold is 0 always survives the peeling: a
    // two-path with thresholds [0, 1] keeps both vertices even though the
    // leaf is initially deficient... no: the leaf has residual 1 >= 1, so
    // nothing is deleted at all.
    let path2 = path_lattice(2);
    let summary = kcore_decompose(
        &path2,
        &CoreThresholds::PerVertex(vec![0, 1]),
        &[true; 2],
        &mut workspace,
    )
    .unwrap();
    assert_eq!(summary.core_vertex_count, 2);
    assert_eq!(summary.peeling_rounds, 0);
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 5: classic deterministic facts.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn one_core_of_connected_graph_is_everything_and_trees_have_no_k2_core() {
    let mut workspace = KCoreWorkspace::default();
    let connected = [
        ("path8", path_lattice(8)),
        ("cycle7", cycle_lattice(7)),
        ("star7", star_lattice(6)),
        ("complete5", complete_lattice(5)),
    ];
    for (name, lattice) in &connected {
        let vertices = lattice.n_sites;
        let full = vec![true; vertices];
        // A connected graph on >= 2 vertices has no isolated vertex, so
        // its 1-core is the whole graph: nobody is ever deficient.
        let summary = kcore_decompose(
            lattice,
            &CoreThresholds::uniform(1).unwrap(),
            &full,
            &mut workspace,
        )
        .unwrap();
        assert_eq!(
            summary.core_vertex_count, vertices,
            "{name}: 1-core = whole"
        );
        assert_eq!(summary.core_edge_count, lattice.n_edges());
        assert_eq!(summary.largest_core_size, vertices);
        assert_eq!(summary.peeling_rounds, 0, "{name}: nothing to peel at k=1");

        // Trees (path and star are trees) always expose a leaf, so the
        // k-core is empty for every k >= 2.
        if name.starts_with("path") || name.starts_with("star") {
            for k in [2usize, 3] {
                let summary = kcore_decompose(
                    lattice,
                    &CoreThresholds::uniform(k).unwrap(),
                    &full,
                    &mut workspace,
                )
                .unwrap();
                assert_eq!(summary.core_vertex_count, 0, "{name}: k={k} core of a tree");
                assert_eq!(summary.deleted_total, vertices);
            }
        }
    }

    // Isolates are the only thing k = 1 removes: a path plus a detached
    // vertex keeps the path and drops the isolate in round 1.
    let lattice = CsrLattice::from_edges(
        4,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
        ],
    );
    let summary = kcore_decompose(
        &lattice,
        &CoreThresholds::uniform(1).unwrap(),
        &[true; 4],
        &mut workspace,
    )
    .unwrap();
    assert_eq!(summary.core_vertex_count, 3);
    assert_eq!(summary.peeling_rounds, 1);
    assert_eq!(workspace.deletion_order(), &[3]);
    assert_eq!(workspace.round_survivors(), &[4.0, 3.0]);
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 6: adapter contract.
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn kcore_mc_full_keep_reduces_to_direct_decomposition() {
    let lattice = path_lattice(16);
    let mut mc = KCoreMC::new(
        lattice.clone(),
        1.0,
        CoreThresholds::uniform(2).unwrap(),
        64,
    )
    .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x004b_434f_5245);
    mc.draw_and_peel(&mut rng);
    // keep = 1 exercises the law's exact shortcut: no RNG draws consumed.
    assert_eq!(mc.initial_activity(), vec![true; 16]);

    let mut workspace = KCoreWorkspace::default();
    let direct = kcore_decompose(
        &lattice,
        &CoreThresholds::uniform(2).unwrap(),
        &[true; 16],
        &mut workspace,
    )
    .unwrap();
    assert_eq!(mc.summary(), direct);
    assert_eq!(mc.deletion_order(), workspace.deletion_order());
    let curve = workspace.round_survivors();
    // The recorded curve is the true curve (budget 64 covers all 9
    // entries) plus the final-size padding out to the fixed shape.
    assert_eq!(&mc.round_survivors()[..curve.len()], curve);
    assert!(mc.round_survivors()[curve.len()..]
        .iter()
        .all(|&value| value == direct.core_vertex_count as f64));
}

#[test]
fn kcore_mc_initial_activity_is_bitwise_site_law() {
    let lattice = star_lattice(9);
    let mut mc = KCoreMC::new(
        lattice.clone(),
        0.61,
        CoreThresholds::uniform(2).unwrap(),
        32,
    )
    .unwrap();
    let mut adapter = StaticPercolationMC::new(
        lattice.clone(),
        StaticLaw::site(0.61).unwrap(),
        ObservablePlan::new([StaticObservable::LargestSize], Vec::new()),
    )
    .unwrap();

    let mut mc_ctx =
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x5173_4b43), 0, 8);
    let mut adapter_ctx =
        Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x5173_4b43), 0, 8);
    for sweep in 0..32 {
        mc.sweep(&mut mc_ctx);
        adapter.sweep(&mut adapter_ctx);
        for index in 0..lattice.n_sites {
            let vertex = lattice.vertex_id(index).unwrap();
            assert_eq!(
                mc.initial_activity()[index],
                adapter.configuration().vertex_active(vertex),
                "activity diverged at sweep {sweep}, vertex {index}"
            );
        }
    }
}

#[test]
fn kcore_mc_records_one_vector_sample_and_exact_scalars() {
    let lattice = cycle_lattice(32);
    let mut mc = KCoreMC::new(
        lattice.clone(),
        0.5,
        CoreThresholds::uniform(2).unwrap(),
        16,
    )
    .unwrap();
    let sweeps = 16;
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(11), 0, 1);
    for _ in 0..sweeps {
        mc.sweep(&mut ctx);
        mc.measure(&mut ctx);
    }
    let (scalars, arrays) = ctx.finalize_measurements_split();
    let array = &arrays["KCoreSurvivorCount"];
    assert_eq!(array.shape(), &[17], "fixed recorded-round budget shape");
    for element in array.elements() {
        assert_eq!(element.n_bins, sweeps, "one sample per realization");
    }
    for name in ["KCoreFinalSize", "KCoreRounds", "KCoreDeletedTotal"] {
        assert_eq!(
            scalars[name].n_bins, sweeps,
            "{name}: one sample per realization"
        );
    }
    // The survivor curve is non-increasing in every realization by
    // construction; the mean must be too.
    for window in array.elements().windows(2) {
        assert!(
            window[0].mean >= window[1].mean,
            "survivor counts are non-increasing"
        );
    }
}

#[test]
fn recorded_round_budget_truncates_and_pads_without_losing_scalars() {
    // path6 / k=2 needs 3 rounds; a budget of 1 records only the initial
    // count and round 1 while the scalars stay exact.
    let lattice = path_lattice(6);
    let mut mc =
        KCoreMC::new(lattice.clone(), 1.0, CoreThresholds::uniform(2).unwrap(), 1).unwrap();
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(3), 0, 8);
    mc.sweep(&mut ctx);
    mc.measure(&mut ctx);
    assert_eq!(mc.round_survivors(), &[6.0, 4.0]);
    assert_eq!(mc.summary().peeling_rounds, 3, "peeling ran to completion");
    assert_eq!(mc.summary().core_vertex_count, 0);
    assert_eq!(ctx.finalize_measurements_split().0["KCoreRounds"].mean, 3.0);

    // complete5 / k=1 never peels: the whole budget pads with the final
    // core size.
    let complete = complete_lattice(5);
    let mut quiet = KCoreMC::new(complete, 1.0, CoreThresholds::uniform(1).unwrap(), 8).unwrap();
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(4), 0, 8);
    quiet.sweep(&mut ctx);
    quiet.measure(&mut ctx);
    assert_eq!(quiet.round_survivors(), &[5.0; 9]);
}

#[test]
fn adapter_rejections_are_typed() {
    let lattice = path_lattice(4);
    assert!(matches!(
        KCoreMC::new(lattice.clone(), 0.5, CoreThresholds::uniform(1).unwrap(), 0),
        Err(cmc_rs::KCoreError::NoRecordedRounds)
    ));
    assert!(matches!(
        KCoreMC::new(lattice, 1.5, CoreThresholds::uniform(1).unwrap(), 4),
        Err(cmc_rs::KCoreError::Probability(_))
    ));
    assert!(matches!(
        KCoreMC::new(
            path_lattice(5),
            0.5,
            CoreThresholds::PerVertex(vec![1; 4]),
            4
        ),
        Err(cmc_rs::KCoreError::ThresholdLength {
            expected: 5,
            actual: 4
        })
    ));
}

// ─────────────────────────────────────────────────────────────────────────
// Gate 7: Janson-Luczak emergence threshold on G(n, m).
// ─────────────────────────────────────────────────────────────────────────

/// psi_j(mu) = P(Poisson(mu) >= j) = 1 - e^-mu sum_{i<j} mu^i/i!.
fn psi(j: u32, mu: f64) -> f64 {
    let mut term = 1.0; // mu^0 / 0!
    let mut cumulative = 1.0;
    for i in 1..j as usize {
        term *= mu / i as f64;
        cumulative += term;
    }
    1.0 - (-mu).exp() * cumulative
}

fn realization_count(default: usize) -> usize {
    std::env::var("SCUTTLE_KCORE_REALIZATIONS")
        .ok()
        .and_then(|raw| raw.trim().parse().ok())
        .unwrap_or(default)
}

#[test]
fn er_emergence_brackets_janson_luczak_threshold() {
    const K: u32 = 3;
    const N: usize = 50_000;

    // Own solver for the paper's fixed-point system: scan mu on a log grid
    // for the minimum of mu / psi_{k-1}(mu), then refine by golden section.
    let objective = |mu: f64| mu / psi(K - 1, mu);
    let mut best_mu = 0.01_f64;
    let mut best = f64::INFINITY;
    let mut mu = 0.01_f64;
    while mu < 20.0 {
        let value = objective(mu);
        if value < best {
            best = value;
            best_mu = mu;
        }
        mu *= 1.0005;
    }
    let mut low = best_mu / 1.05;
    let mut high = best_mu * 1.05;
    let phi = 0.618_033_988_749_895;
    let mut left = high - phi * (high - low);
    let mut right = low + phi * (high - low);
    for _ in 0..200 {
        if objective(left) < objective(right) {
            high = right;
            right = left;
            left = high - phi * (high - low);
        } else {
            low = left;
            left = right;
            right = low + phi * (high - low);
        }
    }
    let lambda_k = objective((low + high) / 2.0);
    eprintln!(
        "[E3 literature] solver: lambda_3 = {lambda_k:.6} at mu* = {:.4} (Janson-Luczak/PSW: 3.350919)",
        (low + high) / 2.0
    );
    // Deterministic identity of the paper's formula (not an MC claim).
    assert!(
        (lambda_k - 3.350_919).abs() < 2.0e-4,
        "lambda_3 solver value {lambda_k} disagrees with 3.350919"
    );

    // Above-threshold core fraction: mu_k(lambda) is the LARGEST solution
    // of mu = lambda psi_{k-1}(mu); iterate from mu = lambda.
    let theory_core_fraction = |mean_degree: f64| {
        let mut mu = mean_degree;
        for _ in 0..100_000 {
            let next = mean_degree * psi(K - 1, mu);
            if (next - mu).abs() < 1.0e-14 {
                mu = next;
                break;
            }
            mu = next;
        }
        psi(K, mu)
    };

    let realizations = realization_count(12).max(2);
    let measure = |mean_degree: f64, base_seed: u64| -> (f64, f64, f64) {
        let edges = (mean_degree * N as f64 / 2.0) as usize;
        let mut fractions = Vec::with_capacity(realizations);
        let mut largest = Vec::with_capacity(realizations);
        let mut workspace = KCoreWorkspace::default();
        for realization in 0..realizations {
            let lattice = erdos_renyi_graph(N, edges, base_seed + 7919 * realization as u64);
            let summary = kcore_decompose(
                &lattice,
                &CoreThresholds::uniform(K as usize).unwrap(),
                &[true; N],
                &mut workspace,
            )
            .unwrap();
            fractions.push(summary.core_vertex_count as f64 / N as f64);
            largest.push(summary.largest_core_size as f64 / N as f64);
        }
        let stats = |values: &[f64]| {
            let n = values.len() as f64;
            let mean = values.iter().sum::<f64>() / n;
            let variance = values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
            (mean, (variance / n).sqrt())
        };
        let (mean, stderr) = stats(&fractions);
        let (largest_mean, _) = stats(&largest);
        (mean, stderr, largest_mean)
    };

    // Below the threshold: the 3-core is extinguished (at most rare small
    // cores, never a giant one). A statistical bracket, not an exact claim.
    let below_c = lambda_k - 0.15;
    let (below_mean, below_stderr, below_largest) = measure(below_c, 0xE3_10_00);
    eprintln!(
        "[E3 literature] below: c = {below_c:.4}, core fraction = {below_mean:.5} ± {below_stderr:.5}, \
         largest-core fraction = {below_largest:.5} ({realizations} realizations)"
    );
    assert!(
        below_mean < 4.0 * below_stderr + 0.004,
        "below-threshold 3-core must vanish: {below_mean:.5} ± {below_stderr:.5}"
    );
    assert!(
        below_largest < 0.05,
        "below threshold no giant 3-core may appear: {below_largest}"
    );

    // Above the threshold: a giant 3-core exists and its fraction agrees
    // with psi_k(mu_k(c)). Slack 0.01 absorbs finite-N corrections.
    let above_c = lambda_k + 0.15;
    let theory = theory_core_fraction(above_c);
    let (above_mean, above_stderr, above_largest) = measure(above_c, 0xE3_20_00);
    eprintln!(
        "[E3 literature] above: c = {above_c:.4}, core fraction = {above_mean:.5} ± {above_stderr:.5} vs \
         psi_3(mu_3(c)) = {theory:.5}, largest-core fraction = {above_largest:.5}"
    );
    assert!(
        above_mean > 0.3,
        "above-threshold 3-core must be giant: {above_mean:.5}"
    );
    assert!(
        (above_mean - theory).abs() < 4.0 * above_stderr + 0.01,
        "above-threshold 3-core fraction {above_mean:.5} ± {above_stderr:.5} disagrees with theory {theory:.5}"
    );
    assert!(
        above_largest > 0.8 * theory,
        "the giant core component must carry most of the core: {above_largest:.5} vs {theory:.5}"
    );
}
