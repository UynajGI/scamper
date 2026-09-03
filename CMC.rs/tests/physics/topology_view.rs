use std::mem::size_of;

use cmc_rs::{
    Bond, BondType, BorrowedUndirectedCsr, EdgeId, GraphView, TopologyError, UndirectedGraphView,
    VertexId,
};

struct UndirectedFixture {
    owned: cmc_rs::CsrLattice,
    offsets: Vec<usize>,
    neighbors: Vec<usize>,
    edge_ids: Vec<usize>,
    endpoints: Vec<[usize; 2]>,
}

fn fixture() -> UndirectedFixture {
    let owned = cmc_rs::CsrLattice::from_edges(
        6,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(0, 1, BondType::Generic, 2.0),
            Bond::new(1, 1, BondType::Generic, 3.0),
            Bond::new(3, 4, BondType::Generic, 4.0),
        ],
    );
    let endpoints = owned
        .edges
        .iter()
        .map(|edge| [edge.source, edge.target])
        .collect();
    UndirectedFixture {
        owned: owned.clone(),
        offsets: owned.offsets,
        neighbors: owned.neighbors,
        edge_ids: owned.edge_ids,
        endpoints,
    }
}

#[test]
fn owned_and_borrowed_views_have_identical_topology() {
    let fixture = fixture();
    let borrowed = BorrowedUndirectedCsr::new(
        &fixture.offsets,
        &fixture.neighbors,
        &fixture.edge_ids,
        &fixture.endpoints,
    )
    .unwrap();
    let owned = fixture.owned;

    assert_eq!(owned.vertex_count(), borrowed.vertex_count());
    assert_eq!(owned.edge_count(), borrowed.edge_count());
    for edge_index in 0..owned.edge_count() {
        let owned_edge = owned.edge_id(edge_index).unwrap();
        let borrowed_edge = borrowed.edge_id(edge_index).unwrap();
        assert_eq!(
            UndirectedGraphView::edge_endpoints(&owned, owned_edge),
            borrowed.edge_endpoints(borrowed_edge)
        );
    }
    for vertex_index in 0..owned.vertex_count() {
        let owned_vertex = owned.vertex_id(vertex_index).unwrap();
        let borrowed_vertex = borrowed.vertex_id(vertex_index).unwrap();
        assert_eq!(
            UndirectedGraphView::incidences(&owned, owned_vertex).collect::<Vec<_>>(),
            borrowed.incidences(borrowed_vertex).collect::<Vec<_>>()
        );
    }

    assert_eq!(
        borrowed.incidences(borrowed.vertex_id(1).unwrap()).count(),
        4
    );
    assert_eq!(
        borrowed.incidences(borrowed.vertex_id(2).unwrap()).count(),
        0
    );
    assert_eq!(
        borrowed.incidences(borrowed.vertex_id(5).unwrap()).count(),
        0
    );
}

#[test]
fn incidence_iterators_exhaust_without_changing_results() {
    let fixture = fixture();
    let borrowed = BorrowedUndirectedCsr::new(
        &fixture.offsets,
        &fixture.neighbors,
        &fixture.edge_ids,
        &fixture.endpoints,
    )
    .unwrap();
    let owned = fixture.owned;
    let mut owned_iter = UndirectedGraphView::incidences(&owned, owned.vertex_id(1).unwrap());
    let mut borrowed_iter = borrowed.incidences(borrowed.vertex_id(1).unwrap());

    assert_eq!(owned_iter.by_ref().count(), 4);
    assert_eq!(borrowed_iter.by_ref().count(), 4);
    assert_eq!(owned_iter.next(), None);
    assert_eq!(borrowed_iter.next(), None);
}

#[test]
fn empty_borrowed_graph_is_valid() {
    let graph = BorrowedUndirectedCsr::new(&[0], &[], &[], &[]).unwrap();
    assert_eq!(graph.vertex_count(), 0);
    assert_eq!(graph.edge_count(), 0);
    assert_eq!(graph.vertex_id(0), None);
    assert_eq!(graph.edge_id(0), None);
}

#[test]
fn view_produced_ids_are_bounded_dense_transparent_indices() {
    let graph = cmc_rs::build_chain(3, false);
    assert_eq!(size_of::<VertexId>(), size_of::<usize>());
    assert_eq!(size_of::<EdgeId>(), size_of::<usize>());
    assert_eq!(graph.vertex_id(2).map(VertexId::index), Some(2));
    assert_eq!(graph.vertex_id(3), None);
    assert_eq!(graph.edge_id(1).map(EdgeId::index), Some(1));
    assert_eq!(graph.edge_id(2), None);
}

#[test]
fn borrowed_csr_rejects_invalid_offsets_and_lengths() {
    assert_eq!(
        BorrowedUndirectedCsr::new(&[], &[], &[], &[]).unwrap_err(),
        TopologyError::OffsetLength {
            expected: 1,
            actual: 0
        }
    );
    assert_eq!(
        BorrowedUndirectedCsr::new(&[1], &[], &[], &[]).unwrap_err(),
        TopologyError::OffsetStart { actual: 1 }
    );
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0, 2, 1], &[0], &[0], &[[0, 0]]).unwrap_err(),
        TopologyError::OffsetsNotMonotone { index: 1 }
    );
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0, 1], &[], &[], &[]).unwrap_err(),
        TopologyError::OffsetEnd {
            expected: 0,
            actual: 1
        }
    );
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0], &[], &[0], &[]).unwrap_err(),
        TopologyError::IncidenceLength {
            neighbors: 0,
            edge_ids: 1
        }
    );
}

#[test]
fn borrowed_csr_rejects_out_of_range_values_including_usize_max() {
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0, 1], &[1], &[0], &[[0, 0]]).unwrap_err(),
        TopologyError::NeighborOutOfRange {
            incidence: 0,
            neighbor: 1,
            vertex_count: 1
        }
    );
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0, 1], &[0], &[usize::MAX], &[[0, 0]]).unwrap_err(),
        TopologyError::EdgeIdOutOfRange {
            incidence: 0,
            edge: usize::MAX,
            edge_count: 1
        }
    );
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0], &[], &[], &[[0, 0]]).unwrap_err(),
        TopologyError::EndpointOutOfRange {
            edge: 0,
            endpoint: 0,
            vertex_count: 0
        }
    );
}

#[test]
fn borrowed_csr_rejects_mismatched_asymmetric_and_excess_incidences() {
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0, 1, 2], &[0, 1], &[0, 0], &[[0, 1]]).unwrap_err(),
        TopologyError::IncidenceMismatch {
            incidence: 0,
            vertex: 0,
            neighbor: 0,
            edge: 0
        }
    );
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0, 0], &[], &[], &[[0, 0]]).unwrap_err(),
        TopologyError::IncidenceMultiplicity { edge: 0, actual: 0 }
    );
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0, 2, 2], &[1, 1], &[0, 0], &[[0, 1]]).unwrap_err(),
        TopologyError::EndpointIncidenceMultiplicity {
            edge: 0,
            source_seen: true,
            target_seen: false
        }
    );
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0, 3], &[0, 0, 0], &[0, 0, 0], &[[0, 0]]).unwrap_err(),
        TopologyError::IncidenceMultiplicity { edge: 0, actual: 3 }
    );
}

#[test]
fn high_degree_parallel_and_self_loop_fixture_validates() {
    let mut edges = (1..2048)
        .map(|target| Bond::new(0, target, BondType::Generic, 1.0))
        .collect::<Vec<_>>();
    edges.push(Bond::new(0, 1, BondType::Generic, 2.0));
    edges.push(Bond::new(0, 0, BondType::Generic, 3.0));
    let owned = cmc_rs::CsrLattice::from_edges(2048, edges);
    let endpoints = owned
        .edges
        .iter()
        .map(|edge| [edge.source, edge.target])
        .collect::<Vec<_>>();
    let borrowed = BorrowedUndirectedCsr::new(
        &owned.offsets,
        &owned.neighbors,
        &owned.edge_ids,
        &endpoints,
    )
    .unwrap();

    assert_eq!(borrowed.edge_count(), 2049);
    assert_eq!(
        borrowed.incidences(borrowed.vertex_id(0).unwrap()).count(),
        2050
    );
}

#[test]
fn borrowed_csr_rejects_large_unreferenced_edge_table() {
    let endpoints = vec![[0, 0]; 4096];
    assert_eq!(
        BorrowedUndirectedCsr::new(&[0, 0], &[], &[], &endpoints).unwrap_err(),
        TopologyError::IncidenceMultiplicity { edge: 0, actual: 0 }
    );
}
