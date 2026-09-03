use cmc_rs::{
    build_chain, build_hypercubic, build_square, cluster_stats, Bond, BondType, ClusterStats,
    CsrLattice, OccupancyState,
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

pub fn sample_and_analyze(
    case: &Case,
    occupancy: &mut OccupancyState,
    rng: &mut Xoshiro256PlusPlus,
) -> ClusterStats {
    occupancy.resample(P_SITE, P_BOND, rng);
    cluster_stats(&case.lattice, occupancy, &case.from, &case.to)
}
