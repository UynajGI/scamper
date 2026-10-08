//! F3 topology-generic component analysis validation.

use cmc_rs::{
    analyze, analyze_with_labels, build_chain, build_honeycomb, build_hypercubic, build_kagome,
    build_square, build_triangular, Bond, BondType, BorrowedUndirectedCsr, BoundaryQuery,
    ComponentLabel, ComponentSummary, ComponentWorkspace, CsrLattice, EdgeActivity, EdgeId,
    GraphView, UndirectedGraphView, VertexActivity, VertexId,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

#[derive(Clone, Debug)]
struct VertexFlags(Vec<bool>);

impl VertexActivity for VertexFlags {
    fn vertex_count(&self) -> usize {
        self.0.len()
    }

    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.0[vertex.index()]
    }
}

#[derive(Clone, Debug)]
struct EdgeFlags(Vec<bool>);

impl EdgeActivity for EdgeFlags {
    fn edge_count(&self) -> usize {
        self.0.len()
    }

    fn edge_active(&self, edge: EdgeId) -> bool {
        self.0[edge.index()]
    }
}

#[derive(Debug, Eq, PartialEq)]
struct ReferenceResult {
    summary: ComponentSummary,
    labels: Vec<ComponentLabel>,
    outcomes: Vec<bool>,
}

fn reference<G: UndirectedGraphView>(
    graph: &G,
    vertices: &VertexFlags,
    edges: &EdgeFlags,
    queries: &[BoundaryQuery],
) -> ReferenceResult {
    let mut active_edge_count = 0usize;
    for edge_index in 0..graph.edge_count() {
        let edge = graph.edge_id(edge_index).expect("dense edge ID");
        let [left, right] = graph.edge_endpoints(edge);
        active_edge_count += usize::from(
            edges.edge_active(edge)
                && vertices.vertex_active(left)
                && vertices.vertex_active(right),
        );
    }

    let mut seen = vec![false; graph.vertex_count()];
    let mut labels = vec![None; graph.vertex_count()];
    let mut component_count = 0usize;
    let mut largest_component = None;
    let mut largest_component_size = 0usize;
    let mut raw_second_moment = 0u128;
    for seed_index in 0..graph.vertex_count() {
        let seed = graph.vertex_id(seed_index).expect("dense vertex ID");
        if seen[seed_index] || !vertices.vertex_active(seed) {
            continue;
        }
        seen[seed_index] = true;
        let mut stack = vec![seed];
        let mut members = Vec::new();
        let mut identity = seed;
        while let Some(vertex) = stack.pop() {
            members.push(vertex);
            identity = identity.min(vertex);
            for incidence in graph.incidences(vertex) {
                if !edges.edge_active(incidence.edge) || !vertices.vertex_active(incidence.neighbor)
                {
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
                labels[from.index()].is_some()
                    && query
                        .to()
                        .iter()
                        .any(|&to| labels[from.index()] == labels[to.index()])
            })
        })
        .collect();
    ReferenceResult {
        summary: ComponentSummary {
            active_vertex_count: vertices.0.iter().filter(|&&active| active).count(),
            active_edge_count,
            component_count,
            largest_component,
            largest_component_size,
            raw_second_moment,
        },
        labels,
        outcomes,
    }
}

fn queries<G: GraphView>(graph: &G) -> Vec<BoundaryQuery> {
    let end = graph.vertex_count().saturating_sub(1);
    vec![
        BoundaryQuery::new(graph, &[0, 0], &[end, end]).expect("valid endpoint query"),
        BoundaryQuery::new(graph, &[0, end], &[end]).expect("valid overlap query"),
        BoundaryQuery::new(graph, &[], &[end]).expect("valid empty query"),
    ]
}

fn compare<G: UndirectedGraphView>(
    graph: &G,
    vertices: &VertexFlags,
    edges: &EdgeFlags,
    queries: &[BoundaryQuery],
    workspace: &mut ComponentWorkspace,
) {
    let expected = reference(graph, vertices, edges, queries);
    let actual = analyze_with_labels(graph, vertices, edges, queries, workspace)
        .expect("matching dimensions");
    assert_eq!(actual.summary(), expected.summary);
    assert_eq!(actual.query_outcomes(), expected.outcomes);
    assert_eq!(actual.labels(), Some(expected.labels.as_slice()));
}

fn random_flags(
    graph: &impl UndirectedGraphView,
    law: usize,
    rng: &mut Xoshiro256PlusPlus,
) -> (VertexFlags, EdgeFlags) {
    let vertices = VertexFlags(
        (0..graph.vertex_count())
            .map(|_| law == 1 || rng.random::<f64>() < 0.57)
            .collect(),
    );
    let edges = EdgeFlags(
        (0..graph.edge_count())
            .map(|_| law == 0 || rng.random::<f64>() < 0.49)
            .collect(),
    );
    (vertices, edges)
}

fn random_graph() -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4633_5241_4e44);
    let mut edges = Vec::new();
    for left in 0..12 {
        for right in left..12 {
            if rng.random::<f64>() < 0.22 {
                edges.push(Bond::new(left, right, BondType::Generic, 1.0));
                if left != right && rng.random::<f64>() < 0.08 {
                    edges.push(Bond::new(left, right, BondType::Generic, 1.0));
                }
            }
        }
    }
    CsrLattice::from_edges(14, edges)
}

#[test]
fn exhaustive_small_site_bond_and_mixed_match_flood_fill() {
    let graph = build_square(2, 2, false);
    let queries = queries(&graph);
    let mut workspace = ComponentWorkspace::new();
    workspace.prepare(4, queries.len(), true).unwrap();
    for law in 0..3 {
        let bits = match law {
            0 => graph.vertex_count(),
            1 => graph.edge_count(),
            _ => graph.vertex_count() + graph.edge_count(),
        };
        for mask in 0..(1usize << bits) {
            let vertices = VertexFlags(
                (0..graph.vertex_count())
                    .map(|index| law == 1 || mask & (1 << index) != 0)
                    .collect(),
            );
            let edges = EdgeFlags(
                (0..graph.edge_count())
                    .map(|index| {
                        law == 0
                            || mask & (1 << (index + usize::from(law == 2) * graph.vertex_count()))
                                != 0
                    })
                    .collect(),
            );
            compare(&graph, &vertices, &edges, &queries, &mut workspace);
        }
    }
}

#[test]
fn all_topologies_and_laws_match_independent_reference() {
    let graphs = vec![
        build_chain(8, false),
        build_square(3, 3, false),
        build_hypercubic(
            &[2, 2, 2],
            &[BondType::CubicX, BondType::CubicY, BondType::CubicZ],
            false,
        ),
        build_triangular(2, 2),
        build_honeycomb(2, 2),
        build_kagome(2, 2),
        random_graph(),
    ];
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4633_5041_5249);
    let mut workspace = ComponentWorkspace::new();
    for graph in &graphs {
        let queries = queries(graph);
        workspace
            .prepare(graph.vertex_count(), queries.len(), true)
            .unwrap();
        for law in 0..3 {
            for _ in 0..128 {
                let (vertices, edges) = random_flags(graph, law, &mut rng);
                compare(graph, &vertices, &edges, &queries, &mut workspace);
            }
        }
    }
}

#[test]
fn owned_and_borrowed_csr_results_are_identical() {
    let graph = random_graph();
    let endpoints: Vec<_> = graph
        .edges
        .iter()
        .map(|edge| [edge.source, edge.target])
        .collect();
    let borrowed = BorrowedUndirectedCsr::new(
        &graph.offsets,
        &graph.neighbors,
        &graph.edge_ids,
        &endpoints,
    )
    .unwrap();
    let vertices = VertexFlags(
        (0..graph.vertex_count())
            .map(|index| index % 3 != 0)
            .collect(),
    );
    let edges = EdgeFlags(
        (0..graph.edge_count())
            .map(|index| index % 2 == 0)
            .collect(),
    );
    let owned_queries = queries(&graph);
    let borrowed_queries = queries(&borrowed);
    let mut owned_workspace = ComponentWorkspace::new();
    let mut borrowed_workspace = ComponentWorkspace::new();
    let owned = analyze_with_labels(
        &graph,
        &vertices,
        &edges,
        &owned_queries,
        &mut owned_workspace,
    )
    .unwrap();
    let borrowed_result = analyze_with_labels(
        &borrowed,
        &vertices,
        &edges,
        &borrowed_queries,
        &mut borrowed_workspace,
    )
    .unwrap();
    assert_eq!(owned.summary(), borrowed_result.summary());
    assert_eq!(owned.query_outcomes(), borrowed_result.query_outcomes());
    assert_eq!(owned.labels(), borrowed_result.labels());
}

#[test]
fn adversarial_multigraph_empty_and_activity_semantics() {
    let graph = CsrLattice::from_edges(
        6,
        vec![
            Bond::new(0, 0, BondType::Generic, 1.0),
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(2, 3, BondType::Generic, 1.0),
        ],
    );
    let vertices = VertexFlags(vec![true, true, false, true, true, false]);
    let edges = EdgeFlags(vec![true, true, true, true]);
    let queries = vec![
        BoundaryQuery::new(&graph, &[0, 0], &[0, 2]).unwrap(),
        BoundaryQuery::new(&graph, &[1], &[3]).unwrap(),
        BoundaryQuery::new(&graph, &[], &[]).unwrap(),
    ];
    let mut workspace = ComponentWorkspace::new();
    let result = analyze_with_labels(&graph, &vertices, &edges, &queries, &mut workspace).unwrap();
    assert_eq!(result.query_outcomes(), &[true, false, false]);
    assert_eq!(
        result.summary(),
        ComponentSummary {
            active_vertex_count: 4,
            active_edge_count: 3,
            component_count: 3,
            largest_component: graph.vertex_id(0),
            largest_component_size: 2,
            raw_second_moment: 6,
        }
    );
    assert_eq!(
        result.labels().unwrap(),
        &[
            graph.vertex_id(0),
            graph.vertex_id(0),
            None,
            graph.vertex_id(3),
            graph.vertex_id(4),
            None,
        ]
    );

    let empty_offsets = [0usize];
    let empty = BorrowedUndirectedCsr::new(&empty_offsets, &[], &[], &[]).unwrap();
    let query = BoundaryQuery::new(&empty, &[], &[]).unwrap();
    let empty_result = analyze(
        &empty,
        &VertexFlags(vec![]),
        &EdgeFlags(vec![]),
        &[query],
        &mut workspace,
    )
    .unwrap();
    assert_eq!(
        empty_result.summary(),
        ComponentSummary {
            active_vertex_count: 0,
            active_edge_count: 0,
            component_count: 0,
            largest_component: None,
            largest_component_size: 0,
            raw_second_moment: 0,
        }
    );
    assert_eq!(empty_result.query_outcomes(), &[false]);
    assert_eq!(empty_result.labels(), None);
}

#[test]
fn identity_is_canonical_minimum_and_ties_are_edge_order_independent() {
    let edges = vec![
        Bond::new(5, 4, BondType::Generic, 1.0),
        Bond::new(4, 3, BondType::Generic, 1.0),
        Bond::new(2, 1, BondType::Generic, 1.0),
        Bond::new(1, 0, BondType::Generic, 1.0),
    ];
    let reversed = edges.iter().cloned().rev().collect();
    for graph in [
        CsrLattice::from_edges(6, edges.clone()),
        CsrLattice::from_edges(6, reversed),
    ] {
        let mut workspace = ComponentWorkspace::new();
        let result = analyze(
            &graph,
            &VertexFlags(vec![true; 6]),
            &EdgeFlags(vec![true; 4]),
            &[],
            &mut workspace,
        )
        .unwrap();
        assert_eq!(result.summary().largest_component, graph.vertex_id(0));
        assert_eq!(result.summary().largest_component_size, 3);
    }
}

#[test]
fn query_validation_and_dimension_errors_are_typed_and_precede_mutation() {
    let graph = build_chain(4, false);
    let error = BoundaryQuery::new(&graph, &[0, 4], &[3]).unwrap_err();
    assert!(error.to_string().contains("position 1"));
    let query = BoundaryQuery::new(&graph, &[0], &[3]).unwrap();
    let mut workspace = ComponentWorkspace::new();
    workspace.prepare(9, 5, true).unwrap();
    let capacities = (
        workspace.vertex_capacity(),
        workspace.query_capacity(),
        workspace.label_capacity(),
    );
    let error = analyze(
        &graph,
        &VertexFlags(vec![true; 3]),
        &EdgeFlags(vec![true; 3]),
        &[query],
        &mut workspace,
    )
    .unwrap_err();
    assert!(error.to_string().contains("vertex activity length"));
    assert_eq!(
        capacities,
        (
            workspace.vertex_capacity(),
            workspace.query_capacity(),
            workspace.label_capacity(),
        )
    );

    let other = build_chain(5, false);
    let foreign_query = BoundaryQuery::new(&other, &[0], &[4]).unwrap();
    let error = analyze(
        &graph,
        &VertexFlags(vec![true; 4]),
        &EdgeFlags(vec![true; 3]),
        &[foreign_query],
        &mut workspace,
    )
    .unwrap_err();
    assert!(error.to_string().contains("boundary query 0"));
}

#[test]
fn edge_length_error_preserves_reusable_workspace_results() {
    let graph = build_chain(5, false);
    let queries = vec![
        BoundaryQuery::new(&graph, &[0], &[4]).unwrap(),
        BoundaryQuery::new(&graph, &[1], &[2]).unwrap(),
    ];
    let vertices = VertexFlags(vec![true, true, true, false, true]);
    let edges = EdgeFlags(vec![true, true, true, true]);
    let mut workspace = ComponentWorkspace::new();
    let before = analyze_with_labels(&graph, &vertices, &edges, &queries, &mut workspace).unwrap();
    let expected_summary = before.summary();
    let expected_outcomes = before.query_outcomes().to_vec();
    let expected_labels = before.labels().unwrap().to_vec();
    let capacities = (
        workspace.vertex_capacity(),
        workspace.query_capacity(),
        workspace.label_capacity(),
    );

    let error = analyze(
        &graph,
        &vertices,
        &EdgeFlags(vec![true; 3]),
        &queries,
        &mut workspace,
    )
    .unwrap_err();
    assert!(error.to_string().contains("edge activity length"));
    assert_eq!(
        capacities,
        (
            workspace.vertex_capacity(),
            workspace.query_capacity(),
            workspace.label_capacity(),
        )
    );

    let after = analyze_with_labels(&graph, &vertices, &edges, &queries, &mut workspace).unwrap();
    assert_eq!(after.summary(), expected_summary);
    assert_eq!(after.query_outcomes(), expected_outcomes);
    assert_eq!(after.labels(), Some(expected_labels.as_slice()));
}

#[test]
fn workspace_reuse_across_sizes_queries_and_labels_has_no_stale_output() {
    let large = build_chain(6, false);
    let large_vertices = VertexFlags(vec![true, true, false, true, true, true]);
    let large_edges = EdgeFlags(vec![true; 5]);
    let large_queries = vec![
        BoundaryQuery::new(&large, &[0], &[1]).unwrap(),
        BoundaryQuery::new(&large, &[0], &[5]).unwrap(),
    ];
    let mut workspace = ComponentWorkspace::new();

    let first = analyze_with_labels(
        &large,
        &large_vertices,
        &large_edges,
        &large_queries,
        &mut workspace,
    )
    .unwrap();
    assert_eq!(first.query_outcomes(), &[true, false]);
    assert_eq!(
        first.labels().unwrap(),
        &[
            large.vertex_id(0),
            large.vertex_id(0),
            None,
            large.vertex_id(3),
            large.vertex_id(3),
            large.vertex_id(3),
        ]
    );

    let small = build_chain(3, false);
    let small_result = analyze(
        &small,
        &VertexFlags(vec![true; 3]),
        &EdgeFlags(vec![true; 2]),
        &[],
        &mut workspace,
    )
    .unwrap();
    assert!(small_result.query_outcomes().is_empty());
    assert_eq!(small_result.labels(), None);
    assert_eq!(small_result.summary().largest_component_size, 3);

    let reordered_queries = vec![
        BoundaryQuery::new(&large, &[0], &[5]).unwrap(),
        BoundaryQuery::new(&large, &[4], &[5]).unwrap(),
        BoundaryQuery::new(&large, &[5], &[5]).unwrap(),
    ];
    let without_labels = analyze(
        &large,
        &large_vertices,
        &large_edges,
        &reordered_queries,
        &mut workspace,
    )
    .unwrap();
    assert_eq!(without_labels.query_outcomes(), &[false, true, true]);
    assert_eq!(without_labels.labels(), None);

    let empty_offsets = [0usize];
    let empty = BorrowedUndirectedCsr::new(&empty_offsets, &[], &[], &[]).unwrap();
    let empty_query = BoundaryQuery::new(&empty, &[], &[]).unwrap();
    let empty_result = analyze_with_labels(
        &empty,
        &VertexFlags(vec![]),
        &EdgeFlags(vec![]),
        &[empty_query],
        &mut workspace,
    )
    .unwrap();
    assert_eq!(empty_result.query_outcomes(), &[false]);
    assert_eq!(empty_result.labels(), Some([].as_slice()));

    let final_result = analyze_with_labels(
        &large,
        &large_vertices,
        &large_edges,
        &large_queries,
        &mut workspace,
    )
    .unwrap();
    assert_eq!(final_result.query_outcomes(), &[true, false]);
    assert_eq!(
        final_result.labels().unwrap(),
        &[
            large.vertex_id(0),
            large.vertex_id(0),
            None,
            large.vertex_id(3),
            large.vertex_id(3),
            large.vertex_id(3),
        ]
    );
}
