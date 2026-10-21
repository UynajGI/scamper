use cmc_rs::{
    analyze, build_chain, build_hypercubic, build_square, Bond, BondType, BoundaryQuery,
    ComponentSummary, ComponentWorkspace, CsrLattice, EdgeActivity, EdgeId, GraphView,
    StaticConfiguration, UndirectedGraphView, VertexActivity, VertexId,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

pub const P_SITE: f64 = 0.5927;
pub const P_BOND: f64 = 0.5;

pub struct Case {
    pub name: &'static str,
    pub lattice: CsrLattice,
    pub from: Vec<usize>,
    pub to: Vec<usize>,
}

fn sparse_er(n: usize, edge_count: usize, seed: u64) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut pairs = Vec::with_capacity(edge_count);
    while pairs.len() < edge_count {
        let left = rng.random_range(0..n);
        let right = rng.random_range(0..n);
        if left != right {
            let pair = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            if !pairs.contains(&pair) {
                pairs.push(pair);
            }
        }
    }
    let edges = pairs
        .into_iter()
        .map(|(left, right)| Bond::new(left, right, BondType::Generic, 1.0))
        .collect();
    CsrLattice::from_edges(n, edges)
}

fn assert_degree_pool(lattice: &CsrLattice, endpoints: &[usize]) {
    assert_eq!(
        endpoints.len(),
        lattice.n_edges() * 2,
        "attachment pool must contain two entries per edge"
    );
    let mut multiplicity = vec![0usize; lattice.n_sites];
    for &vertex in endpoints {
        multiplicity[vertex] += 1;
    }
    for (vertex, count) in multiplicity.into_iter().enumerate() {
        assert_eq!(
            count,
            lattice.degree(vertex),
            "attachment pool multiplicity must equal degree for vertex {vertex}"
        );
    }
}

fn power_law(n: usize, edges_per_vertex: usize, seed: u64) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut edges = Vec::with_capacity((n - 1) * edges_per_vertex);
    let mut endpoints = vec![0usize, 1];
    edges.push(Bond::new(0, 1, BondType::Generic, 1.0));
    for vertex in 2..n {
        let mut targets = Vec::with_capacity(edges_per_vertex);
        while targets.len() < edges_per_vertex.min(vertex) {
            let target = endpoints[rng.random_range(0..endpoints.len())];
            if !targets.contains(&target) {
                targets.push(target);
            }
        }
        for target in targets {
            edges.push(Bond::new(vertex, target, BondType::Generic, 1.0));
            endpoints.push(vertex);
            endpoints.push(target);
        }
    }
    let lattice = CsrLattice::from_edges(n, edges);
    assert_degree_pool(&lattice, &endpoints);
    lattice
}

pub fn case(name: &str) -> Option<Case> {
    match name {
        "chain" => Some(Case {
            name: "chain",
            lattice: build_chain(4096, false),
            from: vec![0],
            to: vec![4095],
        }),
        "square" => Some(Case {
            name: "square",
            lattice: build_square(64, 64, false),
            from: (0..64).map(|row| row * 64).collect(),
            to: (0..64).map(|row| row * 64 + 63).collect(),
        }),
        "cubic" => Some(Case {
            name: "cubic",
            lattice: build_hypercubic(
                &[16, 16, 16],
                &[BondType::CubicX, BondType::CubicY, BondType::CubicZ],
                false,
            ),
            from: (0..4096).filter(|site| site % 16 == 0).collect(),
            to: (0..4096).filter(|site| site % 16 == 15).collect(),
        }),
        "sparse-er" => Some(Case {
            name: "sparse-er",
            lattice: sparse_er(4096, 16_384, 0x4552),
            from: vec![0],
            to: vec![4095],
        }),
        "power-law" => Some(Case {
            name: "power-law",
            lattice: power_law(4096, 4, 0x0050_4f57_4552),
            from: vec![0],
            to: vec![4095],
        }),
        _ => None,
    }
}

pub fn cases() -> Vec<Case> {
    ["chain", "square", "cubic", "sparse-er", "power-law"]
        .into_iter()
        .map(|name| case(name).expect("known benchmark case"))
        .collect()
}

#[derive(Clone, Copy)]
pub enum ReferenceMode {
    Site,
    Bond,
    Mixed,
}

impl ReferenceMode {
    // Each Criterion bench target is a separate crate, so items shared across
    // benches are dead in the ones that do not use them.
    #[allow(dead_code)]
    pub const ALL: [Self; 3] = [Self::Site, Self::Bond, Self::Mixed];

    #[allow(dead_code)]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Site => "site",
            Self::Bond => "bond",
            Self::Mixed => "mixed",
        }
    }
}

pub struct ReferenceConfiguration {
    mode: ReferenceMode,
    sites: Vec<bool>,
    bonds: Vec<bool>,
}

impl ReferenceConfiguration {
    pub fn new(case: &Case, mode: ReferenceMode) -> Self {
        Self {
            mode,
            sites: vec![matches!(mode, ReferenceMode::Bond); case.lattice.vertex_count()],
            bonds: vec![matches!(mode, ReferenceMode::Site); case.lattice.edge_count()],
        }
    }

    pub fn sample(&mut self, rng: &mut Xoshiro256PlusPlus) {
        if !matches!(self.mode, ReferenceMode::Bond) {
            self.sites
                .iter_mut()
                .for_each(|active| *active = rng.random::<f64>() < P_SITE);
        }
        if !matches!(self.mode, ReferenceMode::Site) {
            self.bonds
                .iter_mut()
                .for_each(|active| *active = rng.random::<f64>() < P_BOND);
        }
    }

    #[allow(dead_code)]
    pub fn storage_bytes(&self) -> usize {
        self.sites.capacity() + self.bonds.capacity()
    }
}

impl VertexActivity for ReferenceConfiguration {
    fn vertex_count(&self) -> usize {
        self.sites.len()
    }

    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.sites[vertex.index()]
    }
}

impl EdgeActivity for ReferenceConfiguration {
    fn edge_count(&self) -> usize {
        self.bonds.len()
    }

    fn edge_active(&self, edge: EdgeId) -> bool {
        self.bonds[edge.index()]
    }
}

pub fn sample_and_analyze(
    case: &Case,
    occupancy: &mut ReferenceConfiguration,
    rng: &mut Xoshiro256PlusPlus,
) -> ComponentSummary {
    occupancy.sample(rng);
    let queries = [BoundaryQuery::new(&case.lattice, &case.from, &case.to).unwrap()];
    let mut workspace = ComponentWorkspace::new();
    analyze(
        &case.lattice,
        occupancy,
        occupancy,
        &queries,
        &mut workspace,
    )
    .unwrap()
    .summary()
}

/// Periodic torus benchmark case for the W2 winding analyzer.
// Shared bench-support items are dead in bench targets that do not use them.
#[allow(dead_code)]
pub struct TorusCase {
    pub name: &'static str,
    pub lattice: CsrLattice,
    pub embedding: cmc_rs::LatticeEmbedding,
}

/// Tori with builder embeddings: 4096-site ring, 64x64 square, 16x16x8
/// cubic (D = 1, 2, 3).
#[allow(dead_code)]
pub fn torus_cases() -> Vec<TorusCase> {
    fn torus(
        name: &'static str,
        built: Result<(CsrLattice, cmc_rs::LatticeEmbedding), cmc_rs::EmbeddingError>,
    ) -> TorusCase {
        let (lattice, embedding) = built.expect("torus embedding");
        TorusCase {
            name,
            lattice,
            embedding,
        }
    }
    vec![
        torus("ring", cmc_rs::build_chain_with_embedding(4096, true)),
        torus("square", cmc_rs::build_square_with_embedding(64, 64, true)),
        torus(
            "cubic",
            cmc_rs::build_hypercubic_with_embedding(
                &[16, 16, 8],
                &[BondType::CubicX, BondType::CubicY, BondType::CubicZ],
                true,
            ),
        ),
    ]
}

/// Compact finite-cover lift-BFS winding reference for bench parity checks
/// (same algorithm as `tests/physics/percolation_winding.rs`, trimmed to
/// rank-only output): per component `(canonical vertex, size, rank)`.
///
/// The 5^D cover clips steps beyond radius 2; `clipped` reports whether any
/// step was dropped, in which case the detected ranks are a lower bound and
/// callers must skip parity assertions for that sample. On D = 1 rings every
/// simple cycle is the full ring winding exactly once, so clipping cannot
/// occur there.
#[allow(dead_code)]
pub fn lift_winding_reference(
    lattice: &CsrLattice,
    embedding: &cmc_rs::LatticeEmbedding,
    configuration: &StaticConfiguration,
) -> (Vec<(usize, usize, usize)>, bool) {
    use cmc_rs::{EdgeActivity as _, PeriodicEmbedding, VertexActivity as _};

    let dimensions = embedding.periodic_dimensions();
    let mut adjacency = vec![Vec::<(usize, Vec<i32>)>::new(); lattice.vertex_count()];
    for index in 0..lattice.edge_count() {
        let edge = lattice.edge_id(index).expect("dense edge index");
        if !configuration.edge_active(edge) {
            continue;
        }
        let [source, target] = lattice.edge_endpoints(edge);
        if !configuration.vertex_active(source) || !configuration.vertex_active(target) {
            continue;
        }
        let forward: Vec<i32> = embedding
            .edge_displacement(edge, source)
            .expect("endpoint")
            .iter()
            .collect();
        let backward: Vec<i32> = forward.iter().map(|&component| -component).collect();
        adjacency[source.index()].push((target.index(), forward));
        adjacency[target.index()].push((source.index(), backward));
    }

    let mut seen = vec![false; lattice.vertex_count()];
    let mut clipped = false;
    let mut components = Vec::new();
    for seed in 0..lattice.vertex_count() {
        let seed_id = lattice.vertex_id(seed).expect("dense vertex ID");
        if seen[seed] || !configuration.vertex_active(seed_id) {
            continue;
        }
        seen[seed] = true;
        let mut size = 1usize;
        let mut stack = vec![seed];
        while let Some(vertex) = stack.pop() {
            for &(neighbor, _) in &adjacency[vertex] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    size += 1;
                    stack.push(neighbor);
                }
            }
        }

        let side = 5usize;
        let cell_index = |cell: &[i32]| -> usize {
            cell.iter()
                .enumerate()
                .map(|(axis, &offset)| (offset + 2) as usize * side.pow(axis as u32))
                .sum()
        };
        let total_cells = side.pow(dimensions as u32);
        let mut lifted_seen = vec![false; lattice.vertex_count() * total_cells];
        let mut offsets = Vec::new();
        lifted_seen[seed * total_cells + cell_index(&vec![0; dimensions])] = true;
        let mut queue = vec![(seed, vec![0i32; dimensions])];
        let mut head = 0;
        while head < queue.len() {
            let (vertex, cell) = queue[head].clone();
            head += 1;
            if vertex == seed && cell.iter().any(|&offset| offset != 0) {
                offsets.push(cell.clone());
            }
            for (neighbor, displacement) in &adjacency[vertex] {
                let mut next_cell = cell.clone();
                let mut inside = true;
                for (axis, &step) in displacement.iter().enumerate() {
                    next_cell[axis] += step;
                    inside &= next_cell[axis].abs() <= 2;
                }
                if !inside {
                    clipped = true;
                    continue;
                }
                let slot = *neighbor * total_cells + cell_index(&next_cell);
                if !lifted_seen[slot] {
                    lifted_seen[slot] = true;
                    queue.push((*neighbor, next_cell));
                }
            }
        }
        // Real rank via f64 elimination.
        let mut rows: Vec<Vec<f64>> = offsets
            .iter()
            .map(|vector| vector.iter().map(|&component| component as f64).collect())
            .collect();
        let mut rank = 0;
        for column in 0..dimensions {
            let pivot = (rank..rows.len())
                .filter(|&row| rows[row][column].abs() > 1e-9)
                .max_by(|&a, &b| {
                    rows[a][column]
                        .abs()
                        .partial_cmp(&rows[b][column].abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            let Some(pivot) = pivot else { continue };
            rows.swap(rank, pivot);
            for row in (rank + 1)..rows.len() {
                let factor = rows[row][column] / rows[rank][column];
                let pivot_row = rows[rank].clone();
                for (value, pivot_value) in rows[row].iter_mut().zip(&pivot_row) {
                    *value -= factor * pivot_value;
                }
            }
            rank += 1;
        }
        components.push((seed, size, rank));
    }
    (components, clipped)
}
