//! G2 Achlioptas process validation (roadmap §13).
//!
//! The roadmap gates, verbatim, and how each is realized here:
//!
//! 1. **Deterministic replay of a fixed candidate trace**: a test-side
//!    driver reimplements the documented RNG draw contract (initialization
//!    consumes nothing; per step `k` draws; without replacement: partial
//!    Fisher-Yates `random_range(0..L-j)` with the documented winner
//!    compaction) and drives the process through the public
//!    `step_with_candidates` API. With the same seed it reproduces the
//!    drawing sweep's accepted sequence BIT-TO-BIT; replaying the recorded
//!    trace again (no RNG at all) reproduces it a second time.
//! 2. **Every step against naive BFS**: after EVERY accepted edge of every
//!    realization, component count / largest size / second moment / active
//!    edge count are recomputed from scratch with an independent flood
//!    fill and compared with the G1 trajectory summaries (both sampling
//!    modes, both rules, chain/square/random/multigraph topologies); the
//!    `steps_to_giant` scalar is re-derived from the reference.
//! 3. **Candidate count one reduces to ordinary random edge addition**:
//!    m=1 without replacement draws one uniform edge from the remaining
//!    pool per step — bitwise-identical to a test-side uniform edge
//!    permutation under the same seed, and statistically indistinguishable
//!    from the Z1 `scan_bond` activation order (mean largest size at fixed
//!    occupation counts, combined error bars). m=1 WITH replacement reduces
//!    statistically to the same law (the first-appearance order of iid
//!    draws is a uniform permutation) with exact coupon-collector
//!    bookkeeping and mean-step statistics.
//! 4. **The terminal state contains every permitted edge**: exactly the
//!    universe (full and subset universes) is occupied at termination.
//!
//! Plus: the hand-worked product-vs-sum split (scores computed by hand),
//! replacement semantics (null steps recorded, never filtered), the C0
//! recording shape, and the Riordan & Warnke continuity anchor (Science
//! 335, 1189 (2012), arXiv:1102.5306): the product rule's largest-component
//! curve shows no O(1) jump and the finite-cluster susceptibility at the
//! transition GROWS with N — the signature of a continuous transition.

use std::collections::HashSet;

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    build_chain, build_square, scan_bond, AchlioptasGrid, AchlioptasProcessMC, AchlioptasRule,
    Bond, BondType, CandidateSampling, CsrLattice, EdgeUniverse, GraphView, NewmanZiffWorkspace,
    TiePolicy, UndirectedGraphView,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

// ── topology helpers ──────────────────────────────────────────────────────

fn random_graph(seed: u64) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let n = 10;
    let mut edges = Vec::new();
    for left in 0..n {
        for right in (left + 1)..n {
            if rng.random::<f64>() < 0.35 {
                edges.push(Bond::new(left, right, BondType::Generic, 1.0));
            }
        }
    }
    CsrLattice::try_from_edges(n, edges).expect("random graph is valid")
}

/// 3 vertices, 4 edges: self-loop (0,0), parallel pair (0,1) x2, (1,2).
fn multigraph() -> CsrLattice {
    CsrLattice::try_from_edges(
        3,
        vec![
            Bond::new(0, 0, BondType::Generic, 1.0),
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
        ],
    )
    .expect("multigraph is valid")
}

/// G(n, m) Erdos-Renyi graph via distinct uniform pairs (the bench-support
/// convention).
fn er_graph(vertices: usize, edge_count: usize, seed: u64) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut pairs: Vec<(usize, usize)> = Vec::with_capacity(edge_count);
    let mut seen: HashSet<(usize, usize)> = HashSet::with_capacity(edge_count);
    while pairs.len() < edge_count {
        let left = rng.random_range(0..vertices);
        let right = rng.random_range(0..vertices);
        if left != right {
            let pair = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            if seen.insert(pair) {
                pairs.push(pair);
            }
        }
    }
    let edges = pairs
        .into_iter()
        .map(|(left, right)| Bond::new(left, right, BondType::Generic, 1.0))
        .collect();
    CsrLattice::try_from_edges(vertices, edges).expect("er graph is valid")
}

// ── independent reference ─────────────────────────────────────────────────

/// From-scratch flood-fill facts of the subgraph occupied by the first
/// `count` accepted edges (all vertices active — the Achlioptas convention):
/// `(components, largest size, raw second moment, active edges)`.
fn naive_facts(
    graph: &CsrLattice,
    accepted: &[usize],
    count: usize,
) -> (usize, usize, u128, usize) {
    let vertices = graph.vertex_count();
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertices];
    let mut active_edges = 0_usize;
    for &edge in &accepted[..count] {
        let [left, right] = graph
            .edge_endpoints(graph.edge_id(edge).expect("dense edge"))
            .map(|vertex| vertex.index());
        // Both endpoints are always active; self-loops never merge.
        active_edges += 1;
        if left != right {
            adjacency[left].push(right);
            adjacency[right].push(left);
        }
    }
    let mut seen = vec![false; vertices];
    let mut components = 0_usize;
    let mut largest = 0_usize;
    let mut second_moment = 0_u128;
    let mut stack = Vec::new();
    for seed in 0..vertices {
        if seen[seed] {
            continue;
        }
        seen[seed] = true;
        stack.push(seed);
        let mut size = 0_usize;
        while let Some(vertex) = stack.pop() {
            size += 1;
            for &neighbor in &adjacency[vertex] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        components += 1;
        second_moment += (size as u128) * (size as u128);
        largest = largest.max(size);
    }
    (components, largest, second_moment, active_edges)
}

fn build_process(
    lattice: CsrLattice,
    rule: AchlioptasRule,
    sampling: CandidateSampling,
    m: usize,
    tie: TiePolicy,
    universe: EdgeUniverse,
) -> AchlioptasProcessMC {
    AchlioptasProcessMC::new(
        lattice,
        rule,
        sampling,
        m,
        tie,
        universe,
        AchlioptasGrid::Full,
    )
    .expect("test process is valid")
}

// ── gate 2: every step vs naive BFS ───────────────────────────────────────

#[test]
fn every_step_matches_naive_bfs_recomputation() {
    let cases = [
        ("chain", build_chain(8, false)),
        ("square", build_square(4, 4, false)),
        ("random", random_graph(4)),
        ("multigraph", multigraph()),
    ];
    let configurations = [
        (
            AchlioptasRule::Product,
            CandidateSampling::WithoutReplacement,
            2,
        ),
        (
            AchlioptasRule::Sum,
            CandidateSampling::WithoutReplacement,
            3,
        ),
        (
            AchlioptasRule::Product,
            CandidateSampling::WithReplacement,
            2,
        ),
    ];
    for (name, graph) in cases {
        let vertices = graph.vertex_count();
        for (rule, sampling, m) in configurations {
            let mut process = build_process(
                graph.clone(),
                rule,
                sampling,
                m,
                TiePolicy::LowestEdgeIndex,
                EdgeUniverse::AllEdges,
            );
            let edges = process.permitted_edge_count();
            for seed in 0..3_u64 {
                let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x00A2_0000 | seed);
                process.run_realization(&mut rng);
                let accepted = process.accepted_edges().to_vec();
                assert_eq!(accepted.len(), edges, "{name}: terminal occupation");
                let summaries = process.workspace().trajectory().summaries();
                let steps = process.workspace().trajectory().steps();
                // V vertex-initialization entries, then one per accepted edge.
                assert_eq!(summaries.len(), vertices + edges);
                assert!(steps[..vertices].iter().all(|&step| step == 0));
                let mut giant_reference = None;
                for count in 1..=edges {
                    let (components, largest, moment, active) =
                        naive_facts(&graph, &accepted, count);
                    let context =
                        format!("{name} {rule:?}/{sampling:?} m={m} seed {seed} t={count}");
                    let summary = summaries[vertices + count - 1];
                    assert_eq!(steps[vertices + count - 1], count, "{context}: step coord");
                    assert_eq!(summary.active_vertex_count, vertices, "{context}: vertices");
                    assert_eq!(summary.active_edge_count, active, "{context}: edges");
                    assert_eq!(summary.component_count, components, "{context}: components");
                    assert_eq!(
                        summary.largest_component_size, largest,
                        "{context}: largest"
                    );
                    assert_eq!(summary.raw_second_moment, moment, "{context}: moment");
                    if giant_reference.is_none() && 2 * largest > vertices {
                        giant_reference = Some(count);
                    }
                }
                assert_eq!(
                    process.steps_to_giant(),
                    giant_reference,
                    "{name} {rule:?}/{sampling:?} m={m} seed {seed}: giant step"
                );
                assert_eq!(
                    process.steps() - process.null_steps(),
                    edges,
                    "{name}: bookkeeping identity"
                );
                if matches!(sampling, CandidateSampling::WithoutReplacement) {
                    assert_eq!(process.null_steps(), 0, "{name}: no null steps");
                    assert_eq!(process.steps(), edges, "{name}: one step per edge");
                }
                // The largest curve mirrors the trajectory at every
                // checkpoint (the Full grid records every count).
                let curve = process.largest_curve();
                assert_eq!(curve.len(), edges + 1);
                assert_eq!(curve[0], 1.0, "{name}: V active isolates");
                for count in 1..=edges {
                    let context = format!("{name} {rule:?}/{sampling:?} m={m} seed {seed}");
                    assert_eq!(
                        curve[count],
                        summaries[vertices + count - 1].largest_component_size as f64,
                        "{context}: curve entry {count}"
                    );
                }
            }
        }
        eprintln!("achlioptas BFS parity: {name} ok");
    }
}

// ── gate 1: deterministic replay of a fixed candidate trace ───────────────

/// Test-side driver of the DOCUMENTED RNG contract (module docs of
/// `achlioptas.rs`), feeding `step_with_candidates` — it never calls
/// `run_realization`. Returns the recorded candidate trace and the accepted
/// sequence. Without replacement, the driver maintains its own pool with the
/// documented partial Fisher-Yates draws and the documented winner
/// compaction; that is what makes bit-identical reproduction of the sweep's
/// draws possible.
fn drive_by_contract(
    process: &mut AchlioptasProcessMC,
    rng: &mut Xoshiro256PlusPlus,
    record_trace: bool,
) -> (Vec<Vec<usize>>, Vec<usize>) {
    let permitted = process.permitted_edge_count();
    let m = process.candidates_per_step();
    let without_replacement = matches!(process.sampling(), CandidateSampling::WithoutReplacement);
    assert_eq!(
        permitted,
        process.lattice().edge_count(),
        "the contract driver covers the AllEdges universe"
    );
    let mut pool: Vec<usize> = (0..permitted).collect();
    let mut pool_len = permitted;
    let mut trace = Vec::new();
    let mut accepted = Vec::new();
    let mut added = 0_usize;
    while added < permitted {
        let count = if without_replacement {
            m.min(pool_len)
        } else {
            m
        };
        let mut candidates = Vec::with_capacity(count);
        for draw in 0..count {
            if without_replacement {
                let index = rng.random_range(0..pool_len - draw);
                candidates.push(pool[index]);
                pool.swap(index, pool_len - 1 - draw);
            } else {
                candidates.push(pool[rng.random_range(0..permitted)]);
            }
        }
        let record = process
            .step_with_candidates(&candidates)
            .expect("contract driver proposes legal candidates");
        if record.added_edge {
            added += 1;
            accepted.push(record.selected_edge);
            if without_replacement {
                // Documented compaction: the winner (draw position r, the
                // first occurrence — see the tie policy) leaves via swap
                // with the LAST pool slot L-1.
                let draw_position = candidates
                    .iter()
                    .position(|&edge| edge == record.selected_edge)
                    .expect("winner is a candidate");
                pool.swap(pool_len - 1 - draw_position, pool_len - 1);
                pool_len -= 1;
            }
        }
        if record_trace {
            trace.push(candidates);
        }
    }
    (trace, accepted)
}

#[test]
fn fixed_candidate_trace_replays_deterministically() {
    let cases = [
        ("chain", build_chain(10, false)),
        ("square", build_square(4, 4, false)),
        ("random", random_graph(9)),
    ];
    let configurations = [
        (
            AchlioptasRule::Product,
            CandidateSampling::WithoutReplacement,
            2,
        ),
        (
            AchlioptasRule::Sum,
            CandidateSampling::WithoutReplacement,
            4,
        ),
        (
            AchlioptasRule::Product,
            CandidateSampling::WithReplacement,
            3,
        ),
    ];
    for (name, graph) in cases {
        for (rule, sampling, m) in configurations {
            let configuration = format!("{name} {rule:?}/{sampling:?} m={m}");

            // (a) Two sweeps with the same seed are identical.
            let mut sweep_one = build_process(
                graph.clone(),
                rule,
                sampling,
                m,
                TiePolicy::LowestEdgeIndex,
                EdgeUniverse::AllEdges,
            );
            let mut sweep_two = build_process(
                graph.clone(),
                rule,
                sampling,
                m,
                TiePolicy::LowestEdgeIndex,
                EdgeUniverse::AllEdges,
            );
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0052_4550);
            sweep_one.run_realization(&mut rng);
            sweep_two.run_realization(&mut Xoshiro256PlusPlus::seed_from_u64(0x0052_4550));
            assert_eq!(
                sweep_one.accepted_edges(),
                sweep_two.accepted_edges(),
                "{configuration}: frozen-seed sweep determinism"
            );

            // (b) The contract driver (public step API only) reproduces the
            // sweep bit-for-bit under the same seed.
            let mut driven = build_process(
                graph.clone(),
                rule,
                sampling,
                m,
                TiePolicy::LowestEdgeIndex,
                EdgeUniverse::AllEdges,
            );
            let (trace, accepted) = drive_by_contract(
                &mut driven,
                &mut Xoshiro256PlusPlus::seed_from_u64(0x0052_4550),
                true,
            );
            assert_eq!(
                accepted,
                sweep_one.accepted_edges(),
                "{configuration}: contract driver vs sweep"
            );
            assert_eq!(
                driven.largest_curve(),
                sweep_one.largest_curve(),
                "{configuration}: replayed curve"
            );
            assert_eq!(
                driven.steps_to_giant(),
                sweep_one.steps_to_giant(),
                "{configuration}: replayed giant step"
            );
            assert_eq!(
                driven.workspace().summary(),
                sweep_one.workspace().summary(),
                "{configuration}: replayed terminal summary"
            );

            // (c) Replaying the recorded trace on a fresh adapter, with no
            // RNG anywhere, reproduces the sequence a second time.
            let mut replay = build_process(
                graph.clone(),
                rule,
                sampling,
                m,
                TiePolicy::LowestEdgeIndex,
                EdgeUniverse::AllEdges,
            );
            let mut replayed = Vec::new();
            for step in &trace {
                let record = replay.step_with_candidates(step).expect("recorded step");
                if record.added_edge {
                    replayed.push(record.selected_edge);
                }
            }
            assert_eq!(
                replayed,
                sweep_one.accepted_edges(),
                "{configuration}: pure-trace replay"
            );
            eprintln!("achlioptas trace replay: {configuration} ok");
        }
    }
}

// ── gate 3: candidate count one reduces to random edge addition ──────────

#[test]
fn candidate_count_one_without_replacement_bitwise_matches_uniform_permutation() {
    // The documented m=1 without-replacement contract IS a uniform edge
    // permutation: each step draws random_range(0..remaining) from a pool
    // and removes the drawn edge. The test-side permutation implements
    // exactly that; identical seeds must give identical edge orders.
    let cases = [
        ("chain", build_chain(16, false)),
        ("square", build_square(5, 5, false)),
        ("random", random_graph(12)),
        ("multigraph", multigraph()),
    ];
    for (name, graph) in cases {
        let mut process = build_process(
            graph.clone(),
            AchlioptasRule::Product,
            CandidateSampling::WithoutReplacement,
            1,
            TiePolicy::LowestEdgeIndex,
            EdgeUniverse::AllEdges,
        );
        let edges = process.permitted_edge_count();
        let seed = 0x004D_3120;
        process.run_realization(&mut Xoshiro256PlusPlus::seed_from_u64(seed));

        let mut pool: Vec<usize> = (0..edges).collect();
        let mut remaining = edges;
        let mut reference = Vec::with_capacity(edges);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        while remaining > 0 {
            let index = rng.random_range(0..remaining);
            reference.push(pool[index]);
            pool.swap(index, remaining - 1);
            remaining -= 1;
        }
        assert_eq!(
            process.accepted_edges(),
            reference.as_slice(),
            "{name}: m=1 accepted order == test-side uniform permutation"
        );
        // The candidate IS the added edge: scoring cannot reroute it.
        assert_eq!(process.steps(), edges);
        assert_eq!(process.null_steps(), 0);
        eprintln!("achlioptas m=1 bitwise reduction: {name} ok ({edges} edges)");
    }
}

#[test]
fn candidate_count_one_is_statistically_ordinary_edge_addition() {
    // Both m=1 sampling modes reduce to the uniform random edge process in
    // distribution; compare mean largest sizes at fixed occupation counts
    // against the Z1 bond scan (a fused, independent implementation).
    let graph = build_square(6, 6, false);
    let vertices = graph.vertex_count();
    let edges = graph.edge_count();
    let probes = [edges / 8, edges / 4, edges / 2, 3 * edges / 4];
    let seeds = 40_u64;
    for (rule, sampling) in [
        (
            AchlioptasRule::Product,
            CandidateSampling::WithoutReplacement,
        ),
        (AchlioptasRule::Sum, CandidateSampling::WithReplacement),
    ] {
        let mut process = build_process(
            graph.clone(),
            rule,
            sampling,
            1,
            TiePolicy::LowestEdgeIndex,
            EdgeUniverse::AllEdges,
        );
        let mut sums = vec![0_f64; probes.len()];
        let mut squares = vec![0_f64; probes.len()];
        for seed in 0..seeds {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x005A_3131 + seed);
            process.run_realization(&mut rng);
            let summaries = process.workspace().trajectory().summaries();
            for (position, &count) in probes.iter().enumerate() {
                let value = summaries[vertices + count - 1].largest_component_size as f64;
                sums[position] += value;
                squares[position] += value * value;
            }
        }
        // scan_bond reference over independent seeds.
        let mut scan = NewmanZiffWorkspace::new();
        let mut scan_sums = vec![0_f64; probes.len()];
        let mut scan_squares = vec![0_f64; probes.len()];
        for seed in 0..seeds {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x005A_3131 + 1_000 + seed);
            let curve = scan_bond(&graph, &mut scan, &mut rng).expect("bond scan");
            for (position, &count) in probes.iter().enumerate() {
                let value = curve
                    .point(count)
                    .expect("occupation point")
                    .largest_component_size as f64;
                scan_sums[position] += value;
                scan_squares[position] += value * value;
            }
        }
        for (position, &count) in probes.iter().enumerate() {
            let mean = sums[position] / seeds as f64;
            let scan_mean = scan_sums[position] / seeds as f64;
            let variance = squares[position] / seeds as f64 - mean * mean;
            let scan_variance = scan_squares[position] / seeds as f64 - scan_mean * scan_mean;
            let stderr = (variance / seeds as f64 + scan_variance / seeds as f64)
                .sqrt()
                .max(1e-9);
            let z = (mean - scan_mean) / stderr;
            assert!(
                z.abs() < 4.0,
                "{rule:?}/{sampling:?} occupation {count}: m=1 mean {mean} vs scan {scan_mean}, \
                 z = {z:.2}"
            );
            eprintln!(
                "achlioptas m=1 {rule:?}/{sampling:?} occupation {count}: \
                 means {mean:.2} vs {scan_mean:.2}, z = {z:.2}"
            );
        }
    }
}

#[test]
fn with_replacement_m_one_matches_coupon_collector_statistics() {
    // m=1 with replacement: total steps follow the coupon-collector law
    // (mean U * H_U); accepted count is exactly U; every step adds at most
    // one edge. Statistical bracket with real error bars, never exact.
    let edges = 24_usize;
    let graph = er_graph(14, edges, 0x0043_4F55);
    let mut process = build_process(
        graph,
        AchlioptasRule::Product,
        CandidateSampling::WithReplacement,
        1,
        TiePolicy::LowestEdgeIndex,
        EdgeUniverse::AllEdges,
    );
    let realizations = 60_u64;
    let mut steps_total = 0_f64;
    let mut squares = 0_f64;
    for seed in 0..realizations {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0043_4346 + seed);
        process.run_realization(&mut rng);
        assert_eq!(process.added_edge_count(), edges, "terminal occupation");
        assert_eq!(
            process.steps() - process.null_steps(),
            edges,
            "bookkeeping identity"
        );
        let steps = process.steps() as f64;
        steps_total += steps;
        squares += steps * steps;
    }
    let mean = steps_total / realizations as f64;
    let variance = squares / realizations as f64 - mean * mean;
    let stderr = (variance / realizations as f64).sqrt();
    // H_24 = sum_{k=1..24} 1/k; the coupon-collector mean.
    let harmonic: f64 = (1..=edges).map(|k| 1.0 / k as f64).sum();
    let expected = edges as f64 * harmonic;
    let z = (mean - expected) / stderr.max(1e-9);
    eprintln!(
        "achlioptas with-replacement m=1 steps: mean {mean:.1} vs coupon-collector \
         {expected:.1} (std err {stderr:.1}, z = {z:.2})"
    );
    assert!(
        z.abs() < 5.0,
        "coupon-collector mean not reproduced: {mean} vs {expected} (z = {z:.2})"
    );
}

// ── gate 4: terminal state = every permitted edge ─────────────────────────

#[test]
fn terminal_state_contains_every_permitted_edge() {
    let graph = build_square(4, 5, false);
    let total_edges = graph.edge_count();
    let subset: Vec<usize> = (0..total_edges).step_by(2).collect();
    let universes = [EdgeUniverse::AllEdges, EdgeUniverse::Subset(subset)];
    for universe in &universes {
        let expected: Vec<usize> = match universe {
            EdgeUniverse::AllEdges => (0..total_edges).collect(),
            EdgeUniverse::Subset(edges) => edges.clone(),
        };
        for sampling in [
            CandidateSampling::WithoutReplacement,
            CandidateSampling::WithReplacement,
        ] {
            let mut process = build_process(
                graph.clone(),
                AchlioptasRule::Product,
                sampling,
                2,
                TiePolicy::LowestEdgeIndex,
                match universe {
                    EdgeUniverse::AllEdges => EdgeUniverse::AllEdges,
                    EdgeUniverse::Subset(edges) => EdgeUniverse::Subset(edges.clone()),
                },
            );
            let permitted = process.permitted_edge_count();
            assert_eq!(permitted, expected.len());
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x0054_4552);
            process.run_realization(&mut rng);
            let summary = process.workspace().summary();
            assert_eq!(
                summary.active_edge_count, permitted,
                "every permitted edge is live at termination"
            );
            let mut accepted = process.accepted_edges().to_vec();
            accepted.sort_unstable();
            assert_eq!(accepted, expected, "accepted set == universe set");
            assert_eq!(process.remaining_edges().len(), 0);
            assert_eq!(
                process.workspace().trajectory().summaries().len(),
                graph.vertex_count() + permitted
            );
        }
    }
    eprintln!("achlioptas terminal-state gate ok (full + subset universes)");
}

// ── hand-worked product vs sum split ──────────────────────────────────────

#[test]
fn product_and_sum_rules_accept_different_edges_by_hand() {
    // Components of sizes 5, 2, 3, 1 (built by hand below); candidates
    // edge 7 = (4,10) joining (5,1) and edge 8 = (6,7) joining (2,3):
    //   product scores: 5*1 = 5 vs 2*3 = 6  -> product rule adds edge 7;
    //   sum scores:     5+1 = 6 vs 2+3 = 5  -> sum rule adds edge 8.
    let lattice = CsrLattice::try_from_edges(
        11,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0), // 0  -> {0..4} size 5
            Bond::new(1, 2, BondType::Generic, 1.0), // 1
            Bond::new(2, 3, BondType::Generic, 1.0), // 2
            Bond::new(3, 4, BondType::Generic, 1.0), // 3
            Bond::new(5, 6, BondType::Generic, 1.0), // 4  -> {5,6} size 2
            Bond::new(7, 8, BondType::Generic, 1.0), // 5  -> {7,8,9} size 3
            Bond::new(8, 9, BondType::Generic, 1.0), // 6
            Bond::new(4, 10, BondType::Generic, 1.0), // 7: joins (5,1)
            Bond::new(6, 7, BondType::Generic, 1.0), // 8: joins (2,3)
        ],
    )
    .expect("rule-split graph is valid");
    let mut product = build_process(
        lattice.clone(),
        AchlioptasRule::Product,
        CandidateSampling::WithoutReplacement,
        2,
        TiePolicy::LowestEdgeIndex,
        EdgeUniverse::AllEdges,
    );
    let mut sum = build_process(
        lattice,
        AchlioptasRule::Sum,
        CandidateSampling::WithoutReplacement,
        2,
        TiePolicy::LowestEdgeIndex,
        EdgeUniverse::AllEdges,
    );
    for process in [&mut product, &mut sum] {
        for edge in 0..=6 {
            process.step_with_candidates(&[edge]).unwrap();
        }
        let (components, largest, _, _) = naive_facts(process.lattice(), &[0, 1, 2, 3, 4, 5, 6], 7);
        assert_eq!(
            (components, largest),
            (4, 5),
            "component sizes 5, 2, 3, 1 by hand"
        );
    }
    assert_eq!(
        product.step_with_candidates(&[7, 8]).unwrap().selected_edge,
        7,
        "product: 5*1 < 2*3"
    );
    assert_eq!(
        sum.step_with_candidates(&[8, 7]).unwrap().selected_edge,
        8,
        "sum: 2+3 < 5+1"
    );
    eprintln!("achlioptas product-vs-sum hand-worked split ok");
}

// ── replacement semantics ─────────────────────────────────────────────────

#[test]
fn with_replacement_proposes_occupied_edges_as_null_steps() {
    // Documented decision: an already-added candidate competes normally; if
    // it wins, the step is a null step (draws consumed, nothing added, no
    // filtering/redraw). The occupied edge CAN win: propose it alone.
    let graph = build_chain(6, false);
    let mut process = build_process(
        graph.clone(),
        AchlioptasRule::Product,
        CandidateSampling::WithReplacement,
        2,
        TiePolicy::LowestEdgeIndex,
        EdgeUniverse::AllEdges,
    );
    let record = process.step_with_candidates(&[0]).unwrap();
    assert!(record.added_edge);
    let summary_after = process.workspace().summary();
    for _ in 0..3 {
        let null = process.step_with_candidates(&[0]).unwrap();
        assert_eq!(null.selected_edge, 0);
        assert!(!null.added_edge, "occupied winner adds nothing");
    }
    assert_eq!(process.null_steps(), 3);
    assert_eq!(process.steps(), 4);
    assert_eq!(process.added_edge_count(), 1);
    assert_eq!(
        process.workspace().summary(),
        summary_after,
        "null steps leave the state untouched"
    );
    // A fresh sweep from the same adapter still terminates complete.
    process.run_realization(&mut Xoshiro256PlusPlus::seed_from_u64(0x004E_554C));
    assert_eq!(process.added_edge_count(), graph.edge_count());
    let mut accepted = process.accepted_edges().to_vec();
    accepted.sort_unstable();
    accepted.dedup();
    assert_eq!(accepted.len(), graph.edge_count());
}

// ── C0 recording shape ────────────────────────────────────────────────────

#[test]
fn adapter_records_one_vector_sample_plus_scalars() {
    let lattice = build_square(4, 4, false);
    let edges = lattice.edge_count();
    let step = 4;
    let mut mc = AchlioptasProcessMC::new(
        lattice,
        AchlioptasRule::Product,
        CandidateSampling::WithoutReplacement,
        2,
        TiePolicy::LowestEdgeIndex,
        EdgeUniverse::AllEdges,
        AchlioptasGrid::Every { step },
    )
    .unwrap();
    let sweeps = 8;
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(0x0043_3041), 0, 1);
    for _ in 0..sweeps {
        mc.sweep(&mut ctx);
        mc.measure(&mut ctx);
    }
    let (scalars, arrays) = ctx.finalize_measurements_split();
    // (0..=24).step_by(4) ends exactly at 24 = edges: 7 checkpoints.
    let expected_checkpoints = edges / step + 1;
    let array = &arrays["AchlioptasLargestSize"];
    assert_eq!(
        array.shape(),
        &[expected_checkpoints],
        "Every-grid shape: E/step + 1 checkpoints when step divides E"
    );
    for element in array.elements() {
        assert_eq!(element.n_bins, sweeps, "one vector sample per realization");
    }
    let coordinates = &arrays["AchlioptasOccupiedEdges"];
    assert_eq!(coordinates.shape(), array.shape());
    let expected_occupied: Vec<f64> = (0..=edges)
        .step_by(step)
        .map(|count| count as f64)
        .collect();
    for (element, &value) in coordinates.elements().iter().zip(&expected_occupied) {
        assert!((element.mean - value).abs() < 1e-12, "constant coordinates");
        assert_eq!(element.n_bins, sweeps);
    }
    for name in [
        "AchlioptasGiantReached",
        "AchlioptasSteps",
        "AchlioptasNullSteps",
        "AchlioptasStepsToGiant",
    ] {
        assert_eq!(scalars[name].n_bins, sweeps, "{name}: one per realization");
    }
    assert_eq!(scalars["AchlioptasSteps"].mean, edges as f64);
    assert_eq!(scalars["AchlioptasNullSteps"].mean, 0.0);
    assert_eq!(scalars["AchlioptasGiantReached"].mean, 1.0);
    // Provenance strings are stable for Results metadata attributes.
    assert_eq!(AchlioptasRule::Product.as_provenance(), "product");
    assert_eq!(AchlioptasRule::Sum.as_provenance(), "sum");
    assert_eq!(
        CandidateSampling::WithoutReplacement.as_provenance(),
        "without-replacement"
    );
    assert_eq!(
        CandidateSampling::WithReplacement.as_provenance(),
        "with-replacement"
    );
    assert_eq!(
        TiePolicy::LowestEdgeIndex.as_provenance(),
        "lowest-edge-index"
    );
    assert_eq!(TiePolicy::EarliestDraw.as_provenance(), "earliest-draw");
}

// ── literature anchor: Riordan & Warnke (2012) continuity ─────────────────

/// One realization's largest-size and finite-cluster-susceptibility curves
/// over the whole occupation axis, read from the G1 trajectory (vertex
/// initialization entries skipped).
fn realization_curves(process: &mut AchlioptasProcessMC, seed: u64) -> (Vec<f64>, Vec<f64>) {
    let vertices = process.lattice().vertex_count();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    process.run_realization(&mut rng);
    let summaries = process.workspace().trajectory().summaries();
    let mut largest = Vec::with_capacity(summaries.len() - vertices);
    let mut susceptibility = Vec::with_capacity(summaries.len() - vertices);
    for summary in &summaries[vertices..] {
        let size = summary.largest_component_size as f64;
        largest.push(size);
        let finite = vertices as f64 - size;
        susceptibility.push(if finite > 0.0 {
            (summary.raw_second_moment as f64 - size * size) / finite
        } else {
            f64::NAN
        });
    }
    (largest, susceptibility)
}

fn mean_and_stderr(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let variance = values.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1.0);
    (mean, (variance / n).sqrt())
}

#[test]
fn riordan_warnke_continuity_anchor() {
    // Riordan & Warnke, Science 335, 1189 (2012): Achlioptas process phase
    // transitions (including the product rule) are CONTINUOUS. Statistical
    // signatures on fixed sparse ER graphs (mean degree 4, product rule,
    // m = 2, without replacement), N = 2^12 vs 2^15:
    //   (a) the finite-cluster susceptibility at the transition GROWS with
    //       N (a discontinuous transition would keep it O(1)) — the
    //       headline signature;
    //   (b) at each N the largest-component curve has no O(1) jump: the
    //       max single-window increase of S_max/N SHRINKS as the window is
    //       refined (a genuinely discontinuous transition would saturate at
    //       a fixed jump as the window shrinks), and stays below 1 in
    //       absolute terms.
    // Statistical gates with real error bars; never exact claims.
    let realizations = [48_usize, 24_usize];
    let mut peaks: Vec<(f64, f64)> = Vec::new(); // (peak susceptibility, stderr)
    let mut fine_jumps: Vec<f64> = Vec::new();
    let mut shrink_ratios: Vec<f64> = Vec::new();
    for (exponent, runs) in [12_usize, 15_usize].iter().zip(realizations) {
        let vertices = 1_usize << exponent;
        let mean_degree = 4.0;
        let edges = (vertices as f64 * mean_degree / 2.0) as usize;
        let mut process = build_process(
            er_graph(vertices, edges, 0x0052_5731 + *exponent as u64),
            AchlioptasRule::Product,
            CandidateSampling::WithoutReplacement,
            2,
            TiePolicy::LowestEdgeIndex,
            EdgeUniverse::AllEdges,
        );
        // Transition window 0.4 .. 1.6 edges per vertex, sampled on a fine
        // grid (windows of E/256 edges); coarser windows are subsamples of
        // the same mean curve.
        let stride = (edges / 256).max(1);
        let checkpoints: Vec<usize> = (0..=edges)
            .step_by(stride)
            .filter(|&count| {
                let density = count as f64 / vertices as f64;
                (0.4..=1.6).contains(&density)
            })
            .collect();
        assert!(
            checkpoints.len() >= 8,
            "enough windows inside the transition"
        );

        let mut per_realization_peaks = Vec::with_capacity(runs);
        let mut curve_sum = vec![0_f64; checkpoints.len()];
        for run in 0..runs {
            let seed = 0x0052_4143 + *exponent as u64 * 100 + run as u64;
            let (largest, susceptibility) = realization_curves(&mut process, seed);
            let low = (0.4 * vertices as f64) as usize;
            let high = ((1.6 * vertices as f64) as usize).min(susceptibility.len());
            per_realization_peaks.push(
                susceptibility[low..high]
                    .iter()
                    .copied()
                    .fold(f64::NAN, f64::max),
            );
            for (position, &count) in checkpoints.iter().enumerate() {
                curve_sum[position] += largest[count] / vertices as f64;
            }
        }
        let (peak, peak_stderr) = mean_and_stderr(&per_realization_peaks);
        peaks.push((peak, peak_stderr));
        let mean_curve: Vec<f64> = curve_sum.iter().map(|&v| v / runs as f64).collect();
        // Max increase of S_max/N over windows of (window + 1) checkpoints.
        let jump = |window: usize| {
            mean_curve
                .windows(window + 1)
                .map(|slice| slice[window] - slice[0])
                .fold(0.0_f64, f64::max)
        };
        let fine = jump(1);
        let coarse = jump(4);
        fine_jumps.push(fine);
        shrink_ratios.push(fine / coarse);
        eprintln!(
            "achlioptas anchor N=2^{exponent}: susceptibility peak {peak:.1} +/- \
             {peak_stderr:.1}; max S_max/N jump over {stride}-edge windows {fine:.3}, \
             over {}-edge windows {coarse:.3} (shrink ratio {:.3})",
            4 * stride,
            fine / coarse
        );
    }
    let (peak_small, stderr_small) = peaks[0];
    let (peak_large, stderr_large) = peaks[1];
    let ratio = peak_large / peak_small;
    let peak_difference = peak_large - peak_small;
    let difference_error = (stderr_small * stderr_small + stderr_large * stderr_large).sqrt();
    eprintln!(
        "achlioptas anchor: susceptibility peak grows {peak_small:.1} -> {peak_large:.1} \
         (x{ratio:.2}); difference {peak_difference:.1} +/- {difference_error:.1}"
    );
    assert!(
        peak_difference > 3.0 * difference_error,
        "finite-cluster susceptibility at the transition must grow with N \
         (continuous-transition signature): {peak_large} vs {peak_small} \
         (+/- {difference_error})"
    );
    assert!(
        ratio > 1.25,
        "susceptibility growth x{ratio:.2} is inconsistent with a bounded \
         (discontinuous) susceptibility"
    );
    for (exponent, (&jump, &shrink)) in [12_usize, 15_usize]
        .iter()
        .zip(fine_jumps.iter().zip(shrink_ratios.iter()))
    {
        assert!(
            jump < 1.0,
            "N=2^{exponent}: single-window jump {jump:.3} in S_max/N is O(1)-like"
        );
        assert!(
            shrink < 0.9,
            "N=2^{exponent}: refining the window does not shrink the max jump \
             ({shrink:.3} >= 0.9) — the signature of a discontinuous jump"
        );
    }
}
