//! R1 external bond-configuration analysis validation.
//!
//! Decoupling contract under test: a bond configuration arrives as a raw mask
//! from ANY producer and is analyzed through the shared F3/W2 capabilities —
//! no sampler knowledge exists on the analysis side. Gates:
//!
//! 1. `ExternalBondConfiguration` (and the `analyze_bonds` convenience)
//!    matches an independent adjacency-BFS flood fill on the open-edge
//!    subgraph across chain/square/random topologies, random masks at
//!    p in {0.1, 0.5, 0.9}, multiple seeds: summary (component count,
//!    largest, second moment), canonical labels, and query outcomes.
//! 2. Typed rejections: wrong mask length from both entry points, empty mask
//!    on a non-empty graph, and `BondSample` buffer-length rejection. An
//!    all-`false` mask of correct length is legal (every vertex an isolate).
//!    Multigraphs with parallel edges and self-loops analyze at parity
//!    (an open self-loop is a counted union no-op per F3).
//! 3. Winding bridge: `ExternalBondConfiguration` feeds `analyze_winding`
//!    identically to a `StaticConfiguration` holding the same bonds
//!    (2x2 torus, random draws through both paths).
//! 4. Sampler decoupling: two synthetic producers (uniform Bernoulli and a
//!    checkerboard pattern) analyze identically through the same API.

use cmc_rs::{
    analyze_bonds, analyze_winding, analyze_with_labels, build_chain, build_square,
    build_square_with_embedding, Bond, BondBernoulli, BondConfigurationError, BondSample, BondType,
    BoundaryQuery, ComponentSummary, ComponentWorkspace, CsrLattice, EdgeActivity,
    ExternalBondConfiguration, GraphView, Probability, StaticConfiguration, UndirectedGraphView,
    VertexId, WindingWorkspace,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Independent reference on the open-edge subgraph: every vertex active,
/// adjacency BFS restricted to mask-open edges. Shares no code with the
/// production union-find.
fn reference<G: UndirectedGraphView>(
    graph: &G,
    mask: &[bool],
    queries: &[BoundaryQuery],
) -> (ComponentSummary, Vec<Option<VertexId>>, Vec<bool>) {
    // Every open edge counts as active (F3); a self-loop then unions
    // nothing, which the BFS below reproduces structurally.
    let active_edge_count = mask.iter().filter(|&&open| open).count();

    let mut seen = vec![false; graph.vertex_count()];
    let mut labels: Vec<Option<VertexId>> = vec![None; graph.vertex_count()];
    let mut component_count = 0usize;
    let mut largest_component = None;
    let mut largest_component_size = 0usize;
    let mut raw_second_moment = 0u128;
    for seed_index in 0..graph.vertex_count() {
        if seen[seed_index] {
            continue;
        }
        seen[seed_index] = true;
        let seed = graph.vertex_id(seed_index).expect("dense vertex ID");
        let mut stack = vec![seed];
        let mut members = Vec::new();
        let mut identity = seed;
        while let Some(vertex) = stack.pop() {
            members.push(vertex);
            identity = identity.min(vertex);
            for incidence in graph.incidences(vertex) {
                if !mask[incidence.edge.index()] {
                    continue;
                }
                let neighbor = incidence.neighbor.index();
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(incidence.neighbor);
                }
            }
        }
        for member in &members {
            labels[member.index()] = Some(identity);
        }
        let size = members.len();
        component_count += 1;
        raw_second_moment += (size as u128) * (size as u128);
        if size > largest_component_size
            || (size == largest_component_size
                && largest_component.is_none_or(|current: VertexId| identity < current))
        {
            largest_component = Some(identity);
            largest_component_size = size;
        }
    }

    let outcomes = queries
        .iter()
        .map(|query| {
            query.from().iter().any(|&from| {
                query.to().iter().any(|&to| {
                    labels[from.index()].is_some() && labels[from.index()] == labels[to.index()]
                })
            })
        })
        .collect();

    (
        ComponentSummary {
            active_vertex_count: graph.vertex_count(),
            active_edge_count,
            component_count,
            largest_component,
            largest_component_size,
            raw_second_moment,
        },
        labels,
        outcomes,
    )
}

fn queries<G: GraphView>(graph: &G) -> Vec<BoundaryQuery> {
    let end = graph.vertex_count().saturating_sub(1);
    let before_end = end.saturating_sub(1);
    vec![
        BoundaryQuery::new(graph, &[0], &[end]).expect("valid endpoint query"),
        BoundaryQuery::new(graph, &[0, 1], &[end, before_end]).expect("valid overlap query"),
        BoundaryQuery::new(graph, &[], &[end]).expect("valid empty query"),
    ]
}

/// Compare production analysis of one mask against the reference, through
/// both entry points (direct capability composition and the convenience).
fn assert_mask_at_parity<G: UndirectedGraphView>(
    graph: &G,
    mask: &[bool],
    queries: &[BoundaryQuery],
    workspace: &mut ComponentWorkspace,
    context: &str,
) {
    let (expected_summary, expected_labels, expected_outcomes) = reference(graph, mask, queries);
    let configuration = ExternalBondConfiguration::new(graph, mask)
        .unwrap_or_else(|error| panic!("{context}: construction failed: {error}"));
    let actual = analyze_with_labels(graph, &configuration, &configuration, queries, workspace)
        .unwrap_or_else(|error| panic!("{context}: analysis failed: {error}"));
    assert_eq!(actual.summary(), expected_summary, "{context}: summary");
    assert_eq!(
        actual.labels(),
        Some(expected_labels.as_slice()),
        "{context}: labels"
    );
    assert_eq!(
        actual.query_outcomes(),
        expected_outcomes,
        "{context}: query outcomes"
    );

    let convenience = analyze_bonds(graph, mask, queries, workspace)
        .unwrap_or_else(|error| panic!("{context}: analyze_bonds failed: {error}"));
    assert_eq!(
        convenience.summary(),
        expected_summary,
        "{context}: analyze_bonds summary"
    );
    assert_eq!(
        convenience.query_outcomes(),
        expected_outcomes,
        "{context}: analyze_bonds query outcomes"
    );
}

fn random_mask(edge_count: usize, p: f64, rng: &mut Xoshiro256PlusPlus) -> Vec<bool> {
    (0..edge_count).map(|_| rng.random::<f64>() < p).collect()
}

fn random_multigraph(seed: u64) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut edges = Vec::new();
    for left in 0..12 {
        for right in left..12 {
            if rng.random::<f64>() < 0.25 {
                edges.push(Bond::new(left, right, BondType::Generic, 1.0));
                // Parallel duplicate and (for left == right) a repeated
                // self-loop: both are distinct physical edges.
                if left == right || rng.random::<f64>() < 0.3 {
                    edges.push(Bond::new(left, right, BondType::Generic, 1.0));
                }
            }
        }
    }
    CsrLattice::from_edges(14, edges)
}

#[test]
fn random_masks_match_independent_bfs_reference_across_topologies_and_p() {
    let graphs: Vec<CsrLattice> = vec![
        build_chain(16, false),
        build_square(5, 5, false),
        random_multigraph(0x4633_5241_4e44),
    ];
    for (topology, graph) in graphs.iter().enumerate() {
        let queries = queries(graph);
        let mut workspace = ComponentWorkspace::new();
        workspace
            .prepare(graph.vertex_count(), queries.len(), true)
            .unwrap();
        for p in [0.1f64, 0.5, 0.9] {
            for seed in [1u64, 7, 0x5eed] {
                let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed ^ (p.to_bits() >> 40));
                for draw in 0..16 {
                    let mask = random_mask(graph.edge_count(), p, &mut rng);
                    assert_mask_at_parity(
                        graph,
                        &mask,
                        &queries,
                        &mut workspace,
                        &format!("topology={topology}/p={p}/seed={seed}/draw={draw}"),
                    );
                }
            }
        }
    }
}

#[test]
fn typed_rejections_and_legal_empty_semantics() {
    let graph = build_chain(4, false);
    let queries = queries(&graph);

    // Wrong length from the configuration constructor.
    assert_eq!(
        ExternalBondConfiguration::new(&graph, &[true; 2]).unwrap_err(),
        BondConfigurationError::MaskLength {
            expected: 3,
            actual: 2
        }
    );
    assert_eq!(
        BondSample::from_mask(&graph, &[true; 5]).unwrap_err(),
        BondConfigurationError::MaskLength {
            expected: 3,
            actual: 5
        }
    );
    // Empty mask on a non-empty graph: same typed rejection. It is NOT the
    // all-inactive configuration — that is a full-length all-false mask.
    assert_eq!(
        ExternalBondConfiguration::new(&graph, &[]).unwrap_err(),
        BondConfigurationError::MaskLength {
            expected: 3,
            actual: 0
        }
    );
    // And through the convenience, it is the analyzer's own domain error.
    let mut workspace = ComponentWorkspace::new();
    let error = analyze_bonds(&graph, &[], &queries, &mut workspace).unwrap_err();
    assert!(error.to_string().contains("edge activity length"));

    // Legal all-false full-length mask: nothing open, every vertex isolated.
    let result = analyze_bonds(&graph, &[false; 3], &queries, &mut workspace).unwrap();
    assert_eq!(
        result.summary(),
        ComponentSummary {
            active_vertex_count: 4,
            active_edge_count: 0,
            component_count: 4,
            largest_component: graph.vertex_id(0),
            largest_component_size: 1,
            raw_second_moment: 4,
        }
    );

    // BondSample mask-view buffer of the wrong length is typed too.
    let sample = BondSample::from_mask(&graph, &[true, false, false]).unwrap();
    assert_eq!(sample.open_count(), 1);
    assert_eq!(
        sample.write_mask(&mut [false; 2]).unwrap_err(),
        BondConfigurationError::MaskLength {
            expected: 3,
            actual: 2
        }
    );
}

#[test]
fn multigraph_parallel_edges_and_self_loops_analyze_at_parity() {
    // Self-loop at 0, parallel 1-2 pair, and a bridge 2-3.
    let graph = CsrLattice::from_edges(
        4,
        vec![
            Bond::new(0, 0, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
            Bond::new(2, 3, BondType::Generic, 1.0),
        ],
    );
    let queries = queries(&graph);
    let mut workspace = ComponentWorkspace::new();
    // Exhaustive over all 2^4 masks: parallel edges are independent mask
    // entries; an open self-loop is a counted union no-op.
    for mask_bits in 0..(1usize << graph.edge_count()) {
        let mask: Vec<bool> = (0..graph.edge_count())
            .map(|index| mask_bits & (1 << index) != 0)
            .collect();
        assert_mask_at_parity(
            &graph,
            &mask,
            &queries,
            &mut workspace,
            &format!("mask={mask_bits:?}"),
        );
    }
}

#[test]
fn winding_bridge_matches_static_configuration_with_same_bonds() {
    let (lattice, embedding) = build_square_with_embedding(2, 2, true).expect("2x2 torus");
    let bond_law = BondBernoulli::new(Probability::new(0.5).unwrap().into());
    let mut configuration = StaticConfiguration::new(lattice.vertex_count(), lattice.edge_count());
    let mut static_workspace = WindingWorkspace::new();
    let mut external_workspace = WindingWorkspace::new();
    let mut wrapping_seen = false;
    let mut contractible_seen = false;
    for seed in 0..64u64 {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5749_4e44 ^ seed);
        bond_law
            .sample(&lattice, &mut configuration, &mut rng)
            .expect("matching dimensions");
        // The external mask IS the StaticConfiguration's bond draw, read back
        // entry by entry — the two paths analyze literally the same bonds.
        let mask: Vec<bool> = (0..lattice.edge_count())
            .map(|index| configuration.edge_active(lattice.edge_id(index).unwrap()))
            .collect();
        let external =
            ExternalBondConfiguration::new(&lattice, &mask).expect("mask length matches");
        let static_result = analyze_winding(
            &lattice,
            &embedding,
            &configuration,
            &configuration,
            &mut static_workspace,
        )
        .unwrap();
        let external_result = analyze_winding(
            &lattice,
            &embedding,
            &external,
            &external,
            &mut external_workspace,
        )
        .unwrap();

        assert_eq!(
            static_result.summary(),
            external_result.summary(),
            "seed={seed}: summaries diverge"
        );
        assert_eq!(static_result.dimensions(), external_result.dimensions());
        for (component, (static_info, external_info)) in static_result
            .components()
            .iter()
            .zip(external_result.components())
            .enumerate()
        {
            assert_eq!(
                (
                    static_info.component(),
                    static_info.size(),
                    static_info.rank()
                ),
                (
                    external_info.component(),
                    external_info.size(),
                    external_info.rank()
                ),
                "seed={seed}: component {component} diverges"
            );
            assert_eq!(
                static_result.generators(component).collect::<Vec<_>>(),
                external_result.generators(component).collect::<Vec<_>>(),
                "seed={seed}: component {component} generators diverge"
            );
        }
        wrapping_seen |= external_result.summary().wrapping_component_count > 0;
        contractible_seen |= external_result.summary().wrapping_component_count == 0;
    }
    assert!(
        wrapping_seen && contractible_seen,
        "64 draws at p=0.5 on the 2x2 torus must exercise both wrapping and contractible cases"
    );
}

#[test]
fn different_synthetic_samplers_analyze_identically_through_one_api() {
    let cols = 6;
    let graph = build_square(cols, cols, false);
    let queries = queries(&graph);
    let mut workspace = ComponentWorkspace::new();
    workspace
        .prepare(graph.vertex_count(), queries.len(), true)
        .unwrap();

    // Producer A: independent bond Bernoulli at p = 0.5.
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4245_5250);
    let bernoulli_mask = random_mask(graph.edge_count(), 0.5, &mut rng);

    // Producer B: deterministic geometric checkerboard — edge open iff the
    // (row + column) parity of its lower-index endpoint is even. A completely
    // different draw mechanism that must flow through the same analysis.
    let checkerboard_mask: Vec<bool> = (0..graph.edge_count())
        .map(|index| {
            let edge = graph.edge_id(index).expect("dense edge ID");
            let [left, right] = graph.edge_endpoints(edge);
            let vertex = left.index().min(right.index());
            (vertex / cols + vertex % cols).is_multiple_of(2)
        })
        .collect();

    assert_ne!(
        bernoulli_mask, checkerboard_mask,
        "the two synthetic samplers must produce different configurations"
    );

    // Both producers hit exactly the same decoupled API.
    assert_mask_at_parity(
        &graph,
        &bernoulli_mask,
        &queries,
        &mut workspace,
        "bernoulli",
    );
    assert_mask_at_parity(
        &graph,
        &checkerboard_mask,
        &queries,
        &mut workspace,
        "checkerboard",
    );

    // And the BondSample record round-trips either producer's mask into the
    // same analysis (R3-readiness: an owned open-edge list is enough).
    for (label, mask) in [
        ("bernoulli", &bernoulli_mask),
        ("checkerboard", &checkerboard_mask),
    ] {
        let sample = BondSample::from_mask(&graph, mask).expect("mask length matches");
        let mut rebuilt = vec![false; graph.edge_count()];
        sample
            .write_mask(&mut rebuilt)
            .expect("buffer length matches");
        assert_eq!(rebuilt, *mask, "{label}: BondSample round trip");
        assert_mask_at_parity(&graph, &rebuilt, &queries, &mut workspace, label);
    }
}
