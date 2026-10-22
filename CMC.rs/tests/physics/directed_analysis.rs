//! D2 production validation: static directed analysis.
//!
//! Gates (roadmap §11 D2):
//! 1. **Exhaustive definitional anchor**: every simple digraph on 3 and 4
//!    vertices (2^6 and 2^12 arc subsets) — SCCs and weak components match
//!    a brute-force transitive-closure ground truth computed from the
//!    mathematical definitions (u ~ v iff u reaches v and v reaches u).
//! 2. **Tarjan vs Kosaraju**: an independently written Kosaraju
//!    implementation (test-side, own code) produces identical SCC
//!    partitions on random digraphs across many seeds, plus the
//!    exhaustive small graphs.
//! 3. **Weak components vs undirected flood-fill**: test-side BFS over
//!    arcs-as-undirected matches on the same corpora.
//! 4. **SCC vs graph transpose**: SCC partitions are invariant under arc
//!    reversal (verified on the random corpus).
//! 5. **Closed forms**: directed chain, out-tree, in/out-star, ring, and
//!    complete DAG — component counts, acyclicity, largest sizes, and
//!    reachability relations checked against hand-derived values.
//! 6. **Reachability**: early-exit BFS equals the transitive closure on
//!    the exhaustive corpus; `reachable(v, v)` true; stamped workspace
//!    reuse across many queries stays correct.

use cmc_rs::{
    analyze_strong_components, analyze_weak_components, reachable, Arc, DirectedCsr, GraphView,
    ReachabilityWorkspace, StrongComponentWorkspace, WeakComponentWorkspace,
};

// ── Independent references (no production code beyond construction) ──────

/// Transitive closure by repeated squaring of the adjacency matrix
/// (boolean Floyd-Warshall-style ground truth).
fn closure(vertices: usize, arcs: &[(usize, usize)]) -> Vec<Vec<bool>> {
    let mut reach = vec![vec![false; vertices]; vertices];
    for (vertex, row) in reach.iter_mut().enumerate() {
        row[vertex] = true;
    }
    for &(source, target) in arcs {
        reach[source][target] = true;
    }
    #[allow(clippy::needless_range_loop)]
    for k in 0..vertices {
        let row_k = reach[k].clone();
        for row in reach.iter_mut() {
            if row[k] {
                for (j, flag) in row_k.iter().enumerate() {
                    if *flag {
                        row[j] = true;
                    }
                }
            }
        }
    }
    reach
}

/// Ground-truth SCC canonical labels from the closure definition.
fn truth_scc_labels(vertices: usize, arcs: &[(usize, usize)]) -> Vec<usize> {
    let reach = closure(vertices, arcs);
    let mut labels = vec![usize::MAX; vertices];
    for vertex in 0..vertices {
        if labels[vertex] != usize::MAX {
            continue;
        }
        let canonical = (0..vertices)
            .filter(|&other| reach[vertex][other] && reach[other][vertex])
            .min()
            .expect("vertex is reachable from itself");
        for other in 0..vertices {
            if reach[vertex][other] && reach[other][vertex] {
                labels[other] = canonical;
            }
        }
    }
    labels
}

/// Ground-truth weak-component canonical labels by flood fill.
fn truth_weak_labels(vertices: usize, arcs: &[(usize, usize)]) -> Vec<usize> {
    let mut adjacency = vec![Vec::new(); vertices];
    for &(source, target) in arcs {
        adjacency[source].push(target);
        adjacency[target].push(source);
    }
    let mut labels = vec![usize::MAX; vertices];
    for seed in 0..vertices {
        if labels[seed] != usize::MAX {
            continue;
        }
        let mut component = vec![seed];
        labels[seed] = seed;
        let mut head = 0;
        while head < component.len() {
            let vertex = component[head];
            head += 1;
            for &neighbor in &adjacency[vertex] {
                if labels[neighbor] == usize::MAX {
                    labels[neighbor] = seed;
                    component.push(neighbor);
                }
            }
        }
    }
    labels
}

/// Independent Kosaraju SCC labels (canonical lowest member).
fn kosaraju_labels(vertices: usize, arcs: &[(usize, usize)]) -> Vec<usize> {
    let mut out = vec![Vec::new(); vertices];
    let mut inc = vec![Vec::new(); vertices];
    for &(source, target) in arcs {
        out[source].push(target);
        inc[target].push(source);
    }
    // Finish order on G.
    let mut seen = vec![false; vertices];
    let mut order = Vec::new();
    for seed in 0..vertices {
        if seen[seed] {
            continue;
        }
        let mut stack = vec![(seed, 0)];
        seen[seed] = true;
        while let Some(&mut (vertex, ref mut cursor)) = stack.last_mut() {
            if *cursor < out[vertex].len() {
                let next = out[vertex][*cursor];
                *cursor += 1;
                if !seen[next] {
                    seen[next] = true;
                    stack.push((next, 0));
                }
            } else {
                stack.pop();
                order.push(vertex);
            }
        }
    }
    // Sweep the transpose in reverse finish order.
    let mut labels = vec![usize::MAX; vertices];
    for &seed in order.iter().rev() {
        if labels[seed] != usize::MAX {
            continue;
        }
        let mut component = vec![seed];
        labels[seed] = seed;
        let mut head = 0;
        while head < component.len() {
            let vertex = component[head];
            head += 1;
            for &neighbor in &inc[vertex] {
                if labels[neighbor] == usize::MAX {
                    labels[neighbor] = seed;
                    component.push(neighbor);
                }
            }
        }
        let canonical = *component.iter().min().expect("non-empty");
        for &member in &component {
            labels[member] = canonical;
        }
    }
    labels
}

fn production_scc_labels(graph: &DirectedCsr) -> Vec<usize> {
    let mut workspace = StrongComponentWorkspace::default();
    analyze_strong_components(graph, &mut workspace)
        .unwrap()
        .labels()
        .iter()
        .map(|label| label.index())
        .collect()
}

fn production_weak_labels(graph: &DirectedCsr) -> Vec<usize> {
    let mut workspace = WeakComponentWorkspace::default();
    analyze_weak_components(graph, &mut workspace)
        .unwrap()
        .labels()
        .iter()
        .map(|label| label.index())
        .collect()
}

fn partitions_equal(left: &[usize], right: &[usize]) -> bool {
    left == right
}

// ── Gate 1 + 2 + 3 + 6: exhaustive small digraphs ────────────────────────

#[test]
fn exhaustive_small_digraphs_match_definitions_and_kosaraju() {
    for vertices in [3_usize, 4] {
        let possible: Vec<(usize, usize)> = (0..vertices)
            .flat_map(|source| (0..vertices).map(move |target| (source, target)))
            .collect();
        let total = 1_usize << possible.len();
        for mask in 0..total {
            let arcs: Vec<(usize, usize)> = possible
                .iter()
                .enumerate()
                .filter(|(bit, _)| mask & (1 << bit) != 0)
                .map(|(_, arc)| *arc)
                .collect();
            let graph = DirectedCsr::from_arcs(
                vertices,
                arcs.iter().map(|&(s, t)| Arc::new(s, t)).collect(),
            )
            .unwrap();

            let truth = truth_scc_labels(vertices, &arcs);
            let produced = production_scc_labels(&graph);
            assert!(
                partitions_equal(&truth, &produced),
                "SCC mismatch on V={vertices} arcs={arcs:?}: truth {truth:?} got {produced:?}"
            );
            let kosaraju = kosaraju_labels(vertices, &arcs);
            assert!(
                partitions_equal(&truth, &kosaraju),
                "Kosaraju mismatch on arcs={arcs:?}"
            );

            let weak_truth = truth_weak_labels(vertices, &arcs);
            let weak_produced = production_weak_labels(&graph);
            assert!(
                partitions_equal(&weak_truth, &weak_produced),
                "weak mismatch on arcs={arcs:?}"
            );

            // Reachability equals the closure ground truth.
            let reach_truth = closure(vertices, &arcs);
            let mut reach_workspace = ReachabilityWorkspace::default();
            #[allow(clippy::needless_range_loop)]
            for source in 0..vertices {
                #[allow(clippy::needless_range_loop)]
                for target in 0..vertices {
                    let produced = reachable(
                        &graph,
                        graph.vertex_id(source).unwrap(),
                        graph.vertex_id(target).unwrap(),
                        &mut reach_workspace,
                    )
                    .unwrap();
                    assert_eq!(
                        produced, reach_truth[source][target],
                        "reachability {source}->{target} on arcs={arcs:?}"
                    );
                }
            }
        }
        eprintln!("[D2 exhaustive] all {total} digraphs on {vertices} vertices agree");
    }
}

// ── Gate 2 + 4: random digraphs, Tarjan vs Kosaraju, transpose invariance ─

#[test]
fn random_digraphs_match_kosaraju_and_transpose_invariance() {
    let mut rng = cmc_rs_test_rng(0xd2_4a);
    for round in 0..40 {
        let vertices = 5 + (rng() % 40) as usize;
        let arc_count = (rng() % 200) as usize;
        let mut arcs = Vec::with_capacity(arc_count);
        for _ in 0..arc_count {
            let source = (rng() % vertices as u64) as usize;
            let target = (rng() % vertices as u64) as usize;
            arcs.push((source, target));
        }
        let graph = DirectedCsr::from_arcs(
            vertices,
            arcs.iter().map(|&(s, t)| Arc::new(s, t)).collect(),
        )
        .unwrap();

        let truth = truth_scc_labels(vertices, &arcs);
        let produced = production_scc_labels(&graph);
        assert!(
            partitions_equal(&truth, &produced),
            "round {round}: Tarjan vs definition, arcs={arcs:?}"
        );
        let kosaraju = kosaraju_labels(vertices, &arcs);
        assert!(
            partitions_equal(&truth, &kosaraju),
            "round {round}: Kosaraju vs definition"
        );

        let weak_truth = truth_weak_labels(vertices, &arcs);
        let weak_produced = production_weak_labels(&graph);
        assert!(
            partitions_equal(&weak_truth, &weak_produced),
            "round {round}"
        );

        // Transpose invariance: reversing every arc leaves the SCC
        // partition unchanged.
        let transpose = DirectedCsr::from_arcs(
            vertices,
            arcs.iter().map(|&(s, t)| Arc::new(t, s)).collect(),
        )
        .unwrap();
        let transposed_labels = production_scc_labels(&transpose);
        assert!(
            partitions_equal(&produced, &transposed_labels),
            "round {round}: SCC partition must be transpose-invariant"
        );
    }
}

/// Deterministic xorshift test RNG (independent of production RNG).
fn cmc_rs_test_rng(mut state: u64) -> impl FnMut() -> u64 {
    move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    }
}

// ── Gate 5: closed forms ──────────────────────────────────────────────────

#[test]
fn directed_chains_trees_stars_rings_match_closed_forms() {
    let mut strong = StrongComponentWorkspace::default();
    let mut weak = WeakComponentWorkspace::default();

    // Directed chain 0 -> 1 -> ... -> 9: one weak component, ten singleton
    // SCCs, acyclic; reachability is exactly the chain order.
    let chain_arcs: Vec<Arc> = (0..9).map(|v| Arc::new(v, v + 1)).collect();
    let chain = DirectedCsr::from_arcs(10, chain_arcs.clone()).unwrap();
    let weak_summary = analyze_weak_components(&chain, &mut weak)
        .unwrap()
        .summary();
    assert_eq!(weak_summary.component_count, 1);
    assert_eq!(weak_summary.largest_component_size, 10);
    let strong_result = analyze_strong_components(&chain, &mut strong).unwrap();
    assert_eq!(strong_result.summary().component_count, 10);
    assert!(strong_result.summary().acyclic);
    assert_eq!(strong_result.summary().largest_scc_size, 1);
    let mut reach = ReachabilityWorkspace::default();
    for source in 0..10 {
        for target in 0..10 {
            assert_eq!(
                reachable(
                    &chain,
                    chain.vertex_id(source).unwrap(),
                    chain.vertex_id(target).unwrap(),
                    &mut reach
                )
                .unwrap(),
                source <= target,
                "chain reachability {source}->{target}"
            );
        }
    }

    // Out-star (hub -> leaves) vs in-star (leaves -> hub).
    let out_star =
        DirectedCsr::from_arcs(5, (1..5).map(|leaf| Arc::new(0, leaf)).collect()).unwrap();
    let out_summary = analyze_strong_components(&out_star, &mut strong)
        .unwrap()
        .summary();
    assert_eq!(out_summary.component_count, 5);
    assert!(out_summary.acyclic);
    let in_star =
        DirectedCsr::from_arcs(5, (1..5).map(|leaf| Arc::new(leaf, 0)).collect()).unwrap();
    let in_summary = analyze_strong_components(&in_star, &mut strong)
        .unwrap()
        .summary();
    assert_eq!(in_summary.component_count, 5);
    assert!(in_summary.acyclic);

    // Ring 0 -> 1 -> ... -> n-1 -> 0: one SCC of size n, cyclic.
    let ring =
        DirectedCsr::from_arcs(6, (0..6).map(|v| Arc::new(v, (v + 1) % 6)).collect()).unwrap();
    let ring_result = analyze_strong_components(&ring, &mut strong).unwrap();
    assert_eq!(ring_result.summary().component_count, 1);
    assert_eq!(ring_result.summary().largest_scc_size, 6);
    assert!(!ring_result.summary().acyclic);

    // Complete DAG on 5 vertices (i -> j iff i < j): acyclic, singleton
    // SCCs, everything reaches everything downstream.
    let dag_arcs: Vec<Arc> = (0..5)
        .flat_map(|i| ((i + 1)..5).map(move |j| Arc::new(i, j)))
        .collect();
    let dag = DirectedCsr::from_arcs(5, dag_arcs).unwrap();
    let dag_result = analyze_strong_components(&dag, &mut strong).unwrap();
    assert!(dag_result.summary().acyclic);
    assert_eq!(dag_result.summary().component_count, 5);
    for source in 0..5 {
        for target in 0..5 {
            assert_eq!(
                reachable(
                    &dag,
                    dag.vertex_id(source).unwrap(),
                    dag.vertex_id(target).unwrap(),
                    &mut reach
                )
                .unwrap(),
                source <= target
            );
        }
    }

    // A single self-loop is a cycle: acyclic false, SCC still singleton.
    let looped = DirectedCsr::from_arcs(1, vec![Arc::new(0, 0)]).unwrap();
    let looped_result = analyze_strong_components(&looped, &mut strong).unwrap();
    assert_eq!(looped_result.summary().component_count, 1);
    assert!(!looped_result.summary().acyclic);
}
