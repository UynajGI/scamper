//! G1 incremental component process validation (roadmap §13).
//!
//! The Z1 prefix-parity convention applied to the generic substrate, plus the
//! substrate-specific gates:
//!
//! 1. every prefix of a random interleaved vertex/edge activation sequence
//!    matches a fresh F3 `analyze` pass over the same active sets (chain,
//!    square, random, multigraph topologies);
//! 2. exhaustive activation orderings on small multigraphs match an
//!    independent from-scratch flood-fill reference after EVERY prefix —
//!    online observables and order independence in one gate;
//! 3. the canonical largest-component identity matches F3 labels;
//! 4. driving the substrate with random vertex activations (every edge
//!    activated first) reproduces `scan_site` prefix summaries — a Z1
//!    compatibility demonstration, deliberately a test rather than a
//!    production adapter;
//! 5. reuse: prepare/prepare cycles with different sizes and reset cycles
//!    leave no stale state.

use cmc_rs::{
    analyze, analyze_with_labels, build_chain, build_square, scan_site, Bond, BondType,
    ComponentWorkspace, CsrLattice, EdgeActivity, EdgeId, GraphView, IncrementalSummary,
    IncrementalUnionWorkspace, NewmanZiffWorkspace, UndirectedGraphView, VertexActivity, VertexId,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// One incremental operation: activate a vertex or a physical edge.
#[derive(Clone, Copy)]
enum Op {
    Vertex(usize),
    Edge(usize),
}

impl Op {
    fn all(vertices: usize, edges: usize) -> Vec<Self> {
        (0..vertices)
            .map(Self::Vertex)
            .chain((0..edges).map(Self::Edge))
            .collect()
    }

    fn shuffled(vertices: usize, edges: usize, rng: &mut Xoshiro256PlusPlus) -> Vec<Self> {
        let mut ops = Self::all(vertices, edges);
        for index in (1..ops.len()).rev() {
            let swap_with = rng.random_range(0..=index);
            ops.swap(index, swap_with);
        }
        ops
    }
}

/// Mirrors the active sets of a G1 prefix (F3 activity contract).
struct Occupancy {
    sites: Vec<bool>,
    bonds: Vec<bool>,
}

impl Occupancy {
    fn new(graph: &CsrLattice) -> Self {
        Self {
            sites: vec![false; graph.vertex_count()],
            bonds: vec![false; graph.edge_count()],
        }
    }

    fn apply(&mut self, op: Op) {
        match op {
            Op::Vertex(vertex) => self.sites[vertex] = true,
            Op::Edge(edge) => self.bonds[edge] = true,
        }
    }
}

impl VertexActivity for Occupancy {
    fn vertex_count(&self) -> usize {
        self.sites.len()
    }

    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.sites[vertex.index()]
    }
}

impl EdgeActivity for Occupancy {
    fn edge_count(&self) -> usize {
        self.bonds.len()
    }

    fn edge_active(&self, edge: EdgeId) -> bool {
        self.bonds[edge.index()]
    }
}

/// Independent from-scratch reference: flood fill over the active subgraph.
/// Returns (components, largest size, canonical largest identity, raw second
/// moment, active edges) with the F3 conventions (self-loops count once,
/// parallel edges separately, lowest-index identity wins size ties).
fn naive_facts(
    graph: &CsrLattice,
    occupancy: &Occupancy,
) -> (usize, usize, Option<usize>, u128, usize) {
    let vertices = graph.vertex_count();
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); vertices];
    let mut active_edges = 0_usize;
    for index in 0..graph.edge_count() {
        if !occupancy.bonds[index] {
            continue;
        }
        let [left, right] = graph.edge_endpoints(graph.edge_id(index).expect("dense edge"));
        let (left, right) = (left.index(), right.index());
        if !occupancy.sites[left] || !occupancy.sites[right] {
            continue;
        }
        active_edges += 1;
        if left != right {
            adjacency[left].push(right);
            adjacency[right].push(left);
        }
    }

    let mut seen = vec![false; vertices];
    let mut components = 0_usize;
    let mut largest_size = 0_usize;
    let mut largest_identity = None;
    let mut second_moment = 0_u128;
    let mut stack = Vec::new();
    for seed in 0..vertices {
        if !occupancy.sites[seed] || seen[seed] {
            continue;
        }
        seen[seed] = true;
        stack.push(seed);
        let mut size = 0_usize;
        let mut identity = seed;
        while let Some(vertex) = stack.pop() {
            size += 1;
            identity = identity.min(vertex);
            for &neighbor in &adjacency[vertex] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        components += 1;
        second_moment += (size as u128) * (size as u128);
        if size > largest_size
            || (size == largest_size && largest_identity.is_none_or(|id| identity < id))
        {
            largest_size = size;
            largest_identity = Some(identity);
        }
    }
    (
        components,
        largest_size,
        largest_identity,
        second_moment,
        active_edges,
    )
}

fn f3_summary(graph: &CsrLattice, occupancy: &Occupancy) -> cmc_rs::ComponentSummary {
    let mut workspace = ComponentWorkspace::new();
    analyze(graph, occupancy, occupancy, &[], &mut workspace)
        .expect("occupancy domains match the graph")
        .summary()
}

fn assert_summary_matches_reference(
    context: &str,
    summary: IncrementalSummary,
    components: usize,
    largest_size: usize,
    largest_identity: Option<usize>,
    second_moment: u128,
    active_edges: usize,
) {
    assert_eq!(summary.component_count, components, "{context}: components");
    assert_eq!(
        summary.largest_component_size, largest_size,
        "{context}: largest size"
    );
    assert_eq!(
        summary.largest_component, largest_identity,
        "{context}: canonical largest identity"
    );
    assert_eq!(
        summary.raw_second_moment, second_moment,
        "{context}: second moment"
    );
    assert_eq!(summary.active_edge_count, active_edges, "{context}: edges");
}

/// Layer 1: random interleaved activations, every prefix against a fresh F3
/// `analyze` pass.
#[test]
fn random_activation_prefixes_match_fresh_analysis() {
    let cases = [
        ("chain", build_chain(8, false)),
        ("square", build_square(4, 4, false)),
        ("random", random_graph()),
        ("multigraph", multigraph()),
    ];
    for (name, graph) in cases {
        let vertices = graph.vertex_count();
        let edges = graph.edge_count();
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(vertices, edges).expect("prepare sizes");
        for seed in 0..4_u64 {
            workspace.reset();
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
            let ops = Op::shuffled(vertices, edges, &mut rng);
            let mut occupancy = Occupancy::new(&graph);
            let mut recorded = 0_usize;
            for (position, &op) in ops.iter().enumerate() {
                // The caller-defined step coordinate is the op position here;
                // any monotone or non-monotone law coordinate is legal.
                let step = position * 3 + 7;
                let summary = match op {
                    Op::Vertex(vertex) => workspace
                        .activate_vertex(&graph, graph.vertex_id(vertex).expect("vertex"), step)
                        .unwrap(),
                    Op::Edge(edge) => workspace
                        .activate_edge(&graph, graph.edge_id(edge).expect("edge"), step)
                        .unwrap(),
                };
                occupancy.apply(op);
                recorded += 1;
                let expected = f3_summary(&graph, &occupancy);
                let context = format!("{name} seed {seed} prefix {position}");
                assert_eq!(
                    summary.active_vertex_count, expected.active_vertex_count,
                    "{context}: vertices"
                );
                assert_eq!(
                    summary.component_count, expected.component_count,
                    "{context}: components"
                );
                assert_eq!(
                    summary.largest_component_size, expected.largest_component_size,
                    "{context}: largest"
                );
                assert_eq!(
                    summary.largest_component,
                    expected.largest_component.map(|vertex| vertex.index()),
                    "{context}: canonical identity"
                );
                assert_eq!(
                    summary.raw_second_moment, expected.raw_second_moment,
                    "{context}: second moment"
                );
                assert_eq!(
                    summary.active_edge_count, expected.active_edge_count,
                    "{context}: edges"
                );
                assert_eq!(workspace.summary(), summary, "{context}: live summary");
                // Trajectory bookkeeping: one entry per accepted op, in order.
                assert_eq!(
                    workspace.trajectory().len(),
                    recorded,
                    "{context}: trajectory len"
                );
                assert_eq!(
                    workspace.trajectory().steps()[recorded - 1],
                    step,
                    "{context}: trajectory step"
                );
                assert_eq!(
                    workspace.trajectory().summaries()[recorded - 1],
                    summary,
                    "{context}: trajectory summary"
                );
            }
        }
        eprintln!("incremental random-prefix parity: {name} ok");
    }
}

/// Layer 2: EVERY activation ordering of small multigraphs, every prefix
/// against the independent flood-fill reference.
#[test]
fn exhaustive_orderings_match_from_scratch_counts() {
    let cases = [
        ("multigraph", multigraph()),         // V=3, E=4 (self-loop + parallels)
        ("star", star()),                     // V=4, E=3
        ("path-with-loop", path_with_loop()), // V=3, E=3
    ];
    for (name, graph) in cases {
        let vertices = graph.vertex_count();
        let edges = graph.edge_count();
        let ops = Op::all(vertices, edges);
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(vertices, edges).expect("prepare sizes");
        let mut occupancy = Occupancy::new(&graph);
        let mut order: Vec<Op> = ops.clone();
        let mut permutations = 0_usize;
        enumerate_permutations(&mut order, 0, &mut |order| {
            workspace.reset();
            occupancy.sites.fill(false);
            occupancy.bonds.fill(false);
            for (position, &op) in order.iter().enumerate() {
                match op {
                    Op::Vertex(vertex) => {
                        workspace
                            .activate_vertex(
                                &graph,
                                graph.vertex_id(vertex).expect("vertex"),
                                position,
                            )
                            .unwrap();
                    }
                    Op::Edge(edge) => {
                        workspace
                            .activate_edge(&graph, graph.edge_id(edge).expect("edge"), position)
                            .unwrap();
                    }
                }
                occupancy.apply(op);
                let (components, largest, identity, moment, active_edges) =
                    naive_facts(&graph, &occupancy);
                let context = format!("{name} perm {permutations} prefix {position}");
                assert_summary_matches_reference(
                    &context,
                    workspace.summary(),
                    components,
                    largest,
                    identity,
                    moment,
                    active_edges,
                );
            }
            permutations += 1;
        });
        let expected_permutations = (1..=ops.len()).product::<usize>();
        assert_eq!(
            permutations, expected_permutations,
            "{name}: permutation count"
        );
        eprintln!("incremental exhaustive orderings: {name} ok ({permutations} orders)");
    }
}

/// Layer 3: the canonical identity matches the F3 `analyze_with_labels`
/// labels (a component label IS its lowest vertex ID) at every prefix.
#[test]
fn canonical_identity_matches_f3_labels() {
    let graph = random_graph();
    let vertices = graph.vertex_count();
    let edges = graph.edge_count();
    let mut workspace = IncrementalUnionWorkspace::new();
    workspace.prepare(vertices, edges).expect("prepare sizes");
    let mut f3_workspace = ComponentWorkspace::new();
    for seed in 0..4_u64 {
        workspace.reset();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(100 + seed);
        let ops = Op::shuffled(vertices, edges, &mut rng);
        let mut occupancy = Occupancy::new(&graph);
        for &op in &ops {
            match op {
                Op::Vertex(vertex) => {
                    workspace
                        .activate_vertex(&graph, graph.vertex_id(vertex).expect("vertex"), vertex)
                        .unwrap();
                }
                Op::Edge(edge) => {
                    workspace
                        .activate_edge(&graph, graph.edge_id(edge).expect("edge"), edge)
                        .unwrap();
                }
            }
            occupancy.apply(op);
            let result =
                analyze_with_labels(&graph, &occupancy, &occupancy, &[], &mut f3_workspace)
                    .expect("occupancy domains match");
            let summary = workspace.summary();
            let expected = result.summary();
            match expected.largest_component {
                None => {
                    // Only edges activated so far: nothing is active.
                    assert_eq!(
                        summary.largest_component, None,
                        "seed {seed}: no vertex active"
                    );
                    assert_eq!(summary.component_count, 0);
                    continue;
                }
                Some(largest) => {
                    let identity = largest.index();
                    assert_eq!(
                        summary.largest_component,
                        Some(identity),
                        "seed {seed}: G1 identity vs F3 selected largest"
                    );
                    // The F3 label of the identity vertex is the identity
                    // itself, and the whole component carries that label.
                    let labels = result.labels().expect("labels were requested");
                    assert_eq!(labels[identity], Some(largest));
                    let labeled: usize = (0..vertices)
                        .filter(|&index| labels[index] == Some(largest))
                        .count();
                    assert_eq!(labeled, expected.largest_component_size);
                    assert_eq!(
                        workspace
                            .component_size(graph.vertex_id(identity).expect("identity"))
                            .expect("identity vertex is active"),
                        summary.largest_component_size,
                        "component_size agrees with the selected component"
                    );
                }
            }
        }
    }
}

/// Layer 4: Z1 compatibility — random vertex activations over pre-activated
/// edges reproduce `scan_site` prefix summaries exactly (test-only
/// demonstration; Z1 keeps its own fused scan in production).
#[test]
fn random_vertex_replay_reproduces_scan_site() {
    let cases = [
        ("chain", build_chain(10, false)),
        ("square", build_square(4, 4, false)),
        ("random", random_graph()),
        ("multigraph", multigraph()),
    ];
    for (name, graph) in cases {
        let vertices = graph.vertex_count();
        let edges = graph.edge_count();
        let mut scan_workspace = NewmanZiffWorkspace::new();
        let mut workspace = IncrementalUnionWorkspace::new();
        workspace.prepare(vertices, edges).expect("prepare sizes");
        for seed in 0..4_u64 {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(200 + seed);
            let curve = scan_site(&graph, &mut scan_workspace, &mut rng).expect("scan");

            workspace.reset();
            for index in 0..edges {
                workspace
                    .activate_edge(&graph, graph.edge_id(index).expect("edge"), usize::MAX)
                    .expect("edges activate first");
            }
            assert_eq!(workspace.summary().active_edge_count, 0, "no endpoints yet");
            for (position, &vertex) in curve.activation_order().iter().enumerate() {
                let summary = workspace
                    .activate_vertex(
                        &graph,
                        graph.vertex_id(vertex).expect("vertex"),
                        position + 1,
                    )
                    .expect("each vertex activates once");
                let point = curve.point(position + 1).expect("point per prefix");
                let context = format!("{name} seed {seed} occupation {}", position + 1);
                assert_eq!(
                    summary.active_vertex_count,
                    position + 1,
                    "{context}: vertices"
                );
                assert_eq!(
                    summary.active_edge_count, point.active_edge_count,
                    "{context}: edges"
                );
                assert_eq!(
                    summary.component_count, point.component_count,
                    "{context}: components"
                );
                assert_eq!(
                    summary.largest_component_size, point.largest_component_size,
                    "{context}: largest"
                );
                assert_eq!(
                    summary.raw_second_moment, point.raw_second_moment,
                    "{context}: second moment"
                );
            }
            let terminal = workspace.summary();
            let last = curve.point(curve.trials()).expect("terminal point");
            assert_eq!(terminal.component_count, last.component_count);
            assert_eq!(terminal.raw_second_moment, last.raw_second_moment);
            // All V+E accepted operations fit the exact trajectory budget.
            assert_eq!(workspace.trajectory().len(), vertices + edges);
        }
        eprintln!("incremental Z1 site replay: {name} ok");
    }
}

/// Layer 5: reuse — reset cycles and prepare cycles with different sizes.
#[test]
fn prepare_cycles_and_resets_leave_no_stale_state() {
    let graph = random_graph();
    let vertices = graph.vertex_count();
    let edges = graph.edge_count();
    let mut workspace = IncrementalUnionWorkspace::new();

    // Grow into the random graph, run partially.
    workspace.prepare(vertices, edges).expect("prepare sizes");
    for index in 0..vertices / 2 {
        workspace
            .activate_vertex(&graph, graph.vertex_id(index).expect("vertex"), index)
            .unwrap();
    }
    assert!(!workspace.trajectory().is_empty());

    // Shrink to a small topology: everything resets.
    let small = multigraph();
    workspace.prepare(3, 4).expect("re-prepare smaller");
    assert_eq!(
        workspace.summary(),
        IncrementalSummary {
            active_vertex_count: 0,
            active_edge_count: 0,
            component_count: 0,
            largest_component: None,
            largest_component_size: 0,
            raw_second_moment: 0,
        }
    );
    assert!(workspace.trajectory().is_empty());

    // Two full realizations on the small graph in different random orders
    // reach the same terminal summary and a fresh trajectory each time.
    let mut terminals = Vec::new();
    for seed in 0..8_u64 {
        workspace.reset();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(300 + seed);
        let ops = Op::shuffled(3, 4, &mut rng);
        let mut occupancy = Occupancy::new(&small);
        for (position, &op) in ops.iter().enumerate() {
            match op {
                Op::Vertex(vertex) => {
                    workspace
                        .activate_vertex(&small, small.vertex_id(vertex).expect("vertex"), position)
                        .unwrap();
                }
                Op::Edge(edge) => {
                    workspace
                        .activate_edge(&small, small.edge_id(edge).expect("edge"), position)
                        .unwrap();
                }
            }
            occupancy.apply(op);
        }
        assert_eq!(workspace.trajectory().len(), 7);
        terminals.push(workspace.summary());
        let (components, largest, identity, moment, active_edges) = naive_facts(&small, &occupancy);
        assert_summary_matches_reference(
            "small terminal",
            workspace.summary(),
            components,
            largest,
            identity,
            moment,
            active_edges,
        );
    }
    terminals.dedup();
    assert_eq!(terminals.len(), 1, "order independence across seeds");
    let terminal = terminals[0];
    assert_eq!(terminal.active_vertex_count, 3);
    assert_eq!(terminal.active_edge_count, 4);
    assert_eq!(terminal.component_count, 1);

    // Grow back to the large graph: no stale state from the small one.
    workspace
        .prepare(vertices, edges)
        .expect("re-prepare larger");
    assert_eq!(workspace.summary().component_count, 0);
    assert!(workspace.trajectory().is_empty());
    workspace
        .activate_vertex(&graph, graph.vertex_id(0).expect("vertex"), 0)
        .unwrap();
    assert_eq!(workspace.summary().component_count, 1);
    assert_eq!(workspace.summary().largest_component, Some(0));
}

fn random_graph() -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(4);
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

fn star() -> CsrLattice {
    CsrLattice::try_from_edges(
        4,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(0, 2, BondType::Generic, 1.0),
            Bond::new(0, 3, BondType::Generic, 1.0),
        ],
    )
    .expect("star is valid")
}

fn path_with_loop() -> CsrLattice {
    CsrLattice::try_from_edges(
        3,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
            Bond::new(2, 2, BondType::Generic, 1.0),
        ],
    )
    .expect("path with loop is valid")
}

/// Invoke `visit` on every permutation of `items` (plain recursion; the
/// exhaustive gate is deliberately small).
fn enumerate_permutations(items: &mut [Op], start: usize, visit: &mut dyn FnMut(&[Op])) {
    if start == items.len() {
        visit(items);
        return;
    }
    for index in start..items.len() {
        items.swap(start, index);
        enumerate_permutations(items, start + 1, visit);
        items.swap(start, index);
    }
}
