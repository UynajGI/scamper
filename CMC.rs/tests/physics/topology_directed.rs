//! D1 production validation: directed topology capabilities.
//!
//! Gates (roadmap §11 D1, mirroring the F1 gate structure):
//! 1. Owned/borrowed parity for both directions of adjacency.
//! 2. Empty graphs, isolated vertices, parallel arcs, and self-loops are
//!    first-class.
//! 3. Every malformed representation is rejected with a typed error:
//!    offset shape, incidence lengths, neighbor/arc ranges, per-arc
//!    source/target consistency, and arc multiplicity (missing or
//!    duplicated) on both the out and transpose sides.
//! 4. `DirectedCsr::validate` catches mutated public fields.

use cmc_rs::{
    Arc, BorrowedDirectedCsr, BorrowedDirectedOutCsr, DirectedCsr, DirectedGraphView,
    DirectedInGraphView, DirectedOutGraphView, GraphView, TopologyError, VertexId,
};

fn vertex<G: GraphView>(graph: &G, index: usize) -> VertexId {
    graph.vertex_id(index).expect("dense index")
}

fn sample_arcs() -> Vec<Arc> {
    // 0->1, 0->2, 1->2, 2->0 (ring), 3->3 (self-loop), 1->2 again
    // (parallel with arc 1), 4 isolated.
    vec![
        Arc::new(0, 1),
        Arc::new(0, 2),
        Arc::new(1, 2),
        Arc::new(2, 0),
        Arc::new(3, 3),
        Arc::new(1, 2),
    ]
}

#[test]
fn owned_directed_csr_exposes_both_directions() {
    let graph = DirectedCsr::from_arcs(5, sample_arcs()).unwrap();
    assert_eq!(graph.vertex_count(), 5);
    assert_eq!(graph.arc_count(), 6);

    let out: Vec<Vec<usize>> = (0..5)
        .map(|v| {
            graph
                .out_arcs(vertex(&graph, v))
                .map(|incidence| incidence.neighbor.index())
                .collect()
        })
        .collect();
    assert_eq!(out, vec![vec![1, 2], vec![2, 2], vec![0], vec![3], vec![]]);

    let incoming: Vec<Vec<usize>> = (0..5)
        .map(|v| {
            graph
                .in_arcs(vertex(&graph, v))
                .map(|incidence| incidence.neighbor.index())
                .collect()
        })
        .collect();
    assert_eq!(
        incoming,
        vec![vec![2], vec![0], vec![0, 1, 1], vec![3], vec![]]
    );

    // Arc identity: parallel arcs keep distinct IDs and the correct sense.
    let parallel: Vec<usize> = graph
        .out_arcs(vertex(&graph, 1))
        .map(|incidence| incidence.arc.index())
        .collect();
    assert_eq!(parallel, vec![2, 5]);
    let self_loop = graph.arc_id(4).expect("arc id");
    assert_eq!(
        graph.arc_endpoints(self_loop),
        [vertex(&graph, 3), vertex(&graph, 3)]
    );
    graph.validate().unwrap();
}

#[test]
fn empty_graph_and_isolated_vertices_are_valid() {
    let empty = DirectedCsr::from_arcs(0, Vec::new()).unwrap();
    assert_eq!(empty.vertex_count(), 0);
    assert_eq!(empty.arc_count(), 0);
    empty.validate().unwrap();

    let view = BorrowedDirectedCsr::new(&[0], &[], &[], &[], &[0], &[], &[]).unwrap();
    assert_eq!(view.vertex_count(), 0);
}

#[test]
fn borrowed_views_match_owned_parity() {
    let owned = DirectedCsr::from_arcs(5, sample_arcs()).unwrap();
    let endpoints: Vec<[usize; 2]> = owned
        .arcs
        .iter()
        .map(|arc| [arc.source, arc.target])
        .collect();

    let out_only = BorrowedDirectedOutCsr::new(
        &owned.out_offsets,
        &owned.out_neighbors,
        &owned.out_arc_ids,
        &endpoints,
    )
    .unwrap();
    let full = BorrowedDirectedCsr::new(
        &owned.out_offsets,
        &owned.out_neighbors,
        &owned.out_arc_ids,
        &endpoints,
        &owned.in_offsets,
        &owned.in_neighbors,
        &owned.in_arc_ids,
    )
    .unwrap();

    fn check_out_parity<Owned, View>(owned: &Owned, view: &View)
    where
        Owned: DirectedOutGraphView,
        View: DirectedOutGraphView,
    {
        for v in 0..5 {
            let owned_out: Vec<(usize, usize)> = owned
                .out_arcs(vertex(owned, v))
                .map(|i| (i.neighbor.index(), i.arc.index()))
                .collect();
            let borrowed_out: Vec<(usize, usize)> = view
                .out_arcs(vertex(view, v))
                .map(|i| (i.neighbor.index(), i.arc.index()))
                .collect();
            assert_eq!(owned_out, borrowed_out, "vertex {v} out parity");
        }
    }
    check_out_parity(&owned, &out_only);
    check_out_parity(&owned, &full);
    for v in 0..5 {
        let owned_in: Vec<(usize, usize)> = owned
            .in_arcs(vertex(&owned, v))
            .map(|i| (i.neighbor.index(), i.arc.index()))
            .collect();
        let borrowed_in: Vec<(usize, usize)> = full
            .in_arcs(vertex(&full, v))
            .map(|i| (i.neighbor.index(), i.arc.index()))
            .collect();
        assert_eq!(owned_in, borrowed_in, "vertex {v} in parity");
    }
}

#[test]
fn malformed_representations_are_rejected_with_typed_errors() {
    // A valid base: 0->1.
    let offsets = [0usize, 1, 1];
    let neighbors = [1usize];
    let arc_ids = [0usize];
    let endpoints = [[0usize, 1]];
    assert!(BorrowedDirectedOutCsr::new(&offsets, &neighbors, &arc_ids, &endpoints).is_ok());

    // Incidence length mismatch.
    assert_eq!(
        BorrowedDirectedOutCsr::new(&offsets, &neighbors[1..], &arc_ids, &endpoints).unwrap_err(),
        TopologyError::IncidenceLength {
            neighbors: 0,
            edge_ids: 1
        }
    );
    // Offsets: bad start / decreasing / bad end / wrong length.
    assert_eq!(
        BorrowedDirectedOutCsr::new(&[1, 1, 1], &neighbors, &arc_ids, &endpoints).unwrap_err(),
        TopologyError::OffsetStart { actual: 1 }
    );
    assert_eq!(
        BorrowedDirectedOutCsr::new(&[0, 2, 1], &neighbors, &arc_ids, &endpoints).unwrap_err(),
        TopologyError::OffsetsNotMonotone { index: 1 }
    );
    assert_eq!(
        BorrowedDirectedOutCsr::new(&[0, 1, 2], &neighbors, &arc_ids, &endpoints).unwrap_err(),
        TopologyError::OffsetEnd {
            expected: 1,
            actual: 2
        }
    );
    // Neighbor range.
    assert_eq!(
        BorrowedDirectedOutCsr::new(&offsets, &[7], &arc_ids, &endpoints).unwrap_err(),
        TopologyError::NeighborOutOfRange {
            incidence: 0,
            neighbor: 7,
            vertex_count: 2
        }
    );
    // Arc ID range.
    assert_eq!(
        BorrowedDirectedOutCsr::new(&offsets, &neighbors, &[9], &endpoints).unwrap_err(),
        TopologyError::EdgeIdOutOfRange {
            incidence: 0,
            edge: 9,
            edge_count: 1
        }
    );
    // Arc endpoint range in the table itself.
    assert_eq!(
        BorrowedDirectedOutCsr::new(&offsets, &neighbors, &arc_ids, &[[0, 5]]).unwrap_err(),
        TopologyError::EndpointOutOfRange {
            edge: 0,
            endpoint: 5,
            vertex_count: 2
        }
    );
    // Sense violation: incidence claims 1 -> 1 for arc 0->1.
    assert_eq!(
        BorrowedDirectedOutCsr::new(&[0, 0, 1], &neighbors, &arc_ids, &endpoints).unwrap_err(),
        TopologyError::OutArcMismatch {
            incidence: 0,
            vertex: 1,
            neighbor: 1,
            arc: 0
        }
    );
    // Neighbor inconsistent with the arc table (0 -> 0 claimed).
    assert_eq!(
        BorrowedDirectedOutCsr::new(&offsets, &[0], &arc_ids, &endpoints).unwrap_err(),
        TopologyError::OutArcMismatch {
            incidence: 0,
            vertex: 0,
            neighbor: 0,
            arc: 0
        }
    );
    // Duplicated arc in out incidences (two slots both carry arc 0).
    assert_eq!(
        BorrowedDirectedOutCsr::new(&[0, 2, 2], &[1, 1], &[0, 0], &endpoints).unwrap_err(),
        TopologyError::ArcOutMultiplicity { arc: 0, actual: 2 }
    );
    // Missing arc: two arcs in the table, only one out incidence.
    assert_eq!(
        BorrowedDirectedOutCsr::new(&offsets, &neighbors, &arc_ids, &endpoints_with_second())
            .unwrap_err(),
        TopologyError::ArcOutMultiplicity { arc: 1, actual: 0 }
    );

    // Transpose side: valid base for 0->1 is in-CSR with 1 -> {0}.
    let in_offsets = [0usize, 0, 1];
    let in_neighbors = [0usize];
    let in_arc_ids = [0usize];
    assert!(BorrowedDirectedCsr::new(
        &offsets,
        &neighbors,
        &arc_ids,
        &endpoints,
        &in_offsets,
        &in_neighbors,
        &in_arc_ids,
    )
    .is_ok());
    // Wrong-side incidence: arc 0 placed at vertex 0 (its source).
    assert_eq!(
        BorrowedDirectedCsr::new(
            &offsets,
            &neighbors,
            &arc_ids,
            &endpoints,
            &[0, 1, 1],
            &[1],
            &[0],
        )
        .unwrap_err(),
        TopologyError::InArcMismatch {
            incidence: 0,
            vertex: 0,
            neighbor: 1,
            arc: 0
        }
    );
    // Transpose vertex-count mismatch.
    assert_eq!(
        BorrowedDirectedCsr::new(
            &offsets,
            &neighbors,
            &arc_ids,
            &endpoints,
            &[0, 1],
            &[0],
            &[0],
        )
        .unwrap_err(),
        TopologyError::OffsetLength {
            expected: 3,
            actual: 2
        }
    );
    // Duplicated arc on the transpose side.
    assert_eq!(
        BorrowedDirectedCsr::new(
            &offsets,
            &neighbors,
            &arc_ids,
            &endpoints,
            &[0, 0, 2],
            &[0, 0],
            &[0, 0],
        )
        .unwrap_err(),
        TopologyError::ArcInMultiplicity { arc: 0, actual: 2 }
    );
}

fn endpoints_with_second() -> Vec<[usize; 2]> {
    vec![[0, 1], [1, 0]]
}

#[test]
fn validate_catches_mutated_public_fields() {
    let mut graph = DirectedCsr::from_arcs(3, vec![Arc::new(0, 1), Arc::new(1, 2)]).unwrap();
    graph.validate().unwrap();

    // Corrupt the out incidence's arc id: the slot at vertex 1 now points
    // at arc 0 (which runs 0 -> 1), so the sense check fires before the
    // multiplicity scan ever sees the duplicate.
    graph.out_arc_ids[1] = 0;
    assert_eq!(
        graph.validate().unwrap_err(),
        TopologyError::OutArcMismatch {
            incidence: 1,
            vertex: 1,
            neighbor: 2,
            arc: 0
        }
    );

    // A sense-consistent corruption instead trips multiplicity: two
    // parallel 0 -> 1 arcs, both out slots of vertex 0 claim arc 0.
    let mut graph = DirectedCsr::from_arcs(2, vec![Arc::new(0, 1), Arc::new(0, 1)]).unwrap();
    graph.out_arc_ids[1] = 0;
    assert_eq!(
        graph.validate().unwrap_err(),
        TopologyError::ArcOutMultiplicity { arc: 0, actual: 2 }
    );

    let mut graph = DirectedCsr::from_arcs(3, vec![Arc::new(0, 1), Arc::new(1, 2)]).unwrap();
    // Incidence 0 sits at vertex 1 (target of arc 0) and claims neighbor 2
    // instead of arc 0's source 0.
    graph.in_neighbors[0] = 2;
    assert_eq!(
        graph.validate().unwrap_err(),
        TopologyError::InArcMismatch {
            incidence: 0,
            vertex: 1,
            neighbor: 2,
            arc: 0
        }
    );
}
