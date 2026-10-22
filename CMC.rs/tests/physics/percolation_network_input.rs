//! Network-input production validation (N-track input PR).
//!
//! Gates (roadmap §9 "Network Input"):
//! 1. Edge lists build owned CSR graphs that are structurally identical to
//!    builder-made lattices (same V/E, same component structure).
//! 2. Every rejection is typed and carries the offending 1-based line:
//!    malformed lines, out-of-range endpoints, duplicates and self-loops
//!    under Reject policies.
//! 3. Policies behave as named: parallel edges survive AllowParallel as
//!    distinct physical edges; self-loops survive Allow.
//! 4. Comments and blank lines are skipped.
//! 5. Parsed graphs flow through the N2/N3 layers (profile + facade) and
//!    through the F1 zero-copy borrowed view with identical analysis.

use cmc_rs::{
    analyze, csr_from_edge_list, BorrowedUndirectedCsr, ComponentWorkspace, DuplicateEdgePolicy,
    EdgeListSpec, GraphView, NetworkInputError, NetworkProfile, ObservablePlan, RandomFailureMC,
    SelfLoopPolicy, StaticConfiguration, UndirectedGraphView,
};
use rand::SeedableRng as _;

#[test]
fn parsed_edge_list_matches_builder_lattice() {
    // Star with 4 leaves, written as an edge list with noise lines.
    let text = "# star graph\n0 1\n\n0 2\n  \n0 3\n# center last\n0 4\n";
    let parsed = csr_from_edge_list(EdgeListSpec::simple(5), text).expect("parse");
    let built = cmc_rs::CsrLattice::from_edges(
        5,
        vec![
            cmc_rs::Bond::new(0, 1, cmc_rs::BondType::Generic, 1.0),
            cmc_rs::Bond::new(0, 2, cmc_rs::BondType::Generic, 1.0),
            cmc_rs::Bond::new(0, 3, cmc_rs::BondType::Generic, 1.0),
            cmc_rs::Bond::new(0, 4, cmc_rs::BondType::Generic, 1.0),
        ],
    );

    assert_eq!(parsed.vertex_count(), built.vertex_count());
    assert_eq!(parsed.edge_count(), built.edge_count());
    // Original-degree profile (N2 substrate) is identical.
    let parsed_profile = NetworkProfile::from_graph(&parsed).unwrap();
    let built_profile = NetworkProfile::from_graph(&built).unwrap();
    assert_eq!(parsed_profile.bins(), built_profile.bins());

    // Identical component structure under a fully-active configuration.
    let mut parsed_workspace = ComponentWorkspace::new();
    let mut built_workspace = ComponentWorkspace::new();
    let parsed_all = StaticConfiguration::new(parsed.vertex_count(), parsed.edge_count());
    let built_all = StaticConfiguration::new(built.vertex_count(), built.edge_count());
    let parsed_summary = analyze(
        &parsed,
        &parsed_all,
        &parsed_all,
        &[],
        &mut parsed_workspace,
    )
    .unwrap()
    .summary();
    let built_summary = analyze(&built, &built_all, &built_all, &[], &mut built_workspace)
        .unwrap()
        .summary();
    assert_eq!(parsed_summary, built_summary);
}

#[test]
fn rejections_are_typed_with_line_numbers() {
    let spec = EdgeListSpec::simple(10);

    // Malformed: one token.
    assert_eq!(
        csr_from_edge_list(spec, "0 1\n2\n").unwrap_err(),
        NetworkInputError::Parse {
            line: 2,
            reason: "expected exactly two whitespace-separated integer endpoints"
        }
    );
    // Malformed: three tokens.
    assert_eq!(
        csr_from_edge_list(spec, "0 1 2\n").unwrap_err(),
        NetworkInputError::Parse {
            line: 1,
            reason: "expected exactly two whitespace-separated integer endpoints"
        }
    );
    // Malformed: non-integer (negative endpoints are not usize).
    assert_eq!(
        csr_from_edge_list(spec, "0 -1\n").unwrap_err(),
        NetworkInputError::Parse {
            line: 1,
            reason: "endpoints must be non-negative decimal integers"
        }
    );
    // Endpoint range, reported with the offending value and line.
    assert_eq!(
        csr_from_edge_list(spec, "0 1\n2 10\n").unwrap_err(),
        NetworkInputError::Endpoint {
            line: 2,
            value: 10,
            vertex_count: 10
        }
    );
    // Duplicate under Reject: the second occurrence reports.
    assert_eq!(
        csr_from_edge_list(spec, "1 2\n3 4\n2 1\n").unwrap_err(),
        NetworkInputError::DuplicateEdge {
            line: 3,
            left: 1,
            right: 2
        }
    );
    // Self-loop under Reject.
    assert_eq!(
        csr_from_edge_list(spec, "3 3\n").unwrap_err(),
        NetworkInputError::SelfLoop { line: 1, vertex: 3 }
    );

    // Display strings carry the line number for every variant.
    let rendered = NetworkInputError::DuplicateEdge {
        line: 7,
        left: 1,
        right: 2,
    }
    .to_string();
    assert!(rendered.contains("line 7"), "{rendered}");
}

#[test]
fn policies_keep_parallel_edges_and_self_loops_when_allowed() {
    // Same list as the rejection test, but the multigraph contract accepts
    // it; parallel edges stay distinct physical edges.
    let text = "1 2\n3 4\n2 1\n3 3\n";
    let strict = csr_from_edge_list(EdgeListSpec::simple(5), text);
    assert!(strict.is_err());

    let multi = csr_from_edge_list(EdgeListSpec::multigraph(5), text).expect("multigraph");
    assert_eq!(multi.edge_count(), 4);
    // The self-loop contributes two incidences at vertex 3 (CSR
    // convention), so its original degree is 3, not 2.
    let profile = NetworkProfile::from_graph(&multi).unwrap();
    assert_eq!(profile.degree(3), 3);

    // Mixed policies are expressible: parallel allowed, self-loops not.
    let mixed = EdgeListSpec {
        vertex_count: 5,
        duplicate_policy: DuplicateEdgePolicy::AllowParallel,
        self_loop_policy: SelfLoopPolicy::Reject,
    };
    assert_eq!(
        csr_from_edge_list(mixed, text).unwrap_err(),
        NetworkInputError::SelfLoop { line: 4, vertex: 3 }
    );
}

#[test]
fn parsed_graph_feeds_n2_facade_and_borrowed_view() {
    // A 10-vertex chain as text.
    let text: String = (0..9).map(|v| format!("{v} {}\n", v + 1)).collect();
    let lattice = csr_from_edge_list(EdgeListSpec::simple(10), &text).expect("chain");

    // N2 facade constructs and runs on the parsed graph.
    let mut mc =
        RandomFailureMC::new(lattice.clone(), 0.6, ObservablePlan::all(Vec::new())).unwrap();
    let mut ctx = carlo_rs::Context::new_with_binsize(
        rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1),
        0,
        4,
    );
    use carlo_rs::MonteCarlo as _;
    for _ in 0..4 {
        mc.sweep(&mut ctx);
        mc.measure(&mut ctx);
    }
    let (scalars, _) = ctx.finalize_measurements_split();
    assert!(scalars.contains_key("GiantFraction"));

    // F1 zero-copy view over the parsed lattice's complete CSR arrays
    // analyzes identically.
    let endpoints: Vec<[usize; 2]> = lattice
        .edges
        .iter()
        .map(|bond| [bond.source, bond.target])
        .collect();
    let view = BorrowedUndirectedCsr::new(
        &lattice.offsets,
        &lattice.neighbors,
        &lattice.edge_ids,
        &endpoints,
    )
    .expect("validated CSR arrays");
    assert_eq!(view.vertex_count(), lattice.vertex_count());
    assert_eq!(view.edge_count(), lattice.edge_count());

    let all_active = StaticConfiguration::new(lattice.vertex_count(), lattice.edge_count());
    let mut owned_workspace = ComponentWorkspace::new();
    let mut view_workspace = ComponentWorkspace::new();
    let owned = analyze(
        &lattice,
        &all_active,
        &all_active,
        &[],
        &mut owned_workspace,
    )
    .unwrap()
    .summary();
    let borrowed = analyze(&view, &all_active, &all_active, &[], &mut view_workspace)
        .unwrap()
        .summary();
    assert_eq!(owned, borrowed);
}
