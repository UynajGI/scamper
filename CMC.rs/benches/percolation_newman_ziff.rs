use std::hint::black_box;
use std::time::Duration;

use cmc_rs::{
    analyze, reweight, ComponentWorkspace, CsrLattice, EdgeActivity, EdgeId, GraphView,
    NewmanZiffWorkspace, Probability, ReweightWorkspace, ScanObservable, UndirectedGraphView,
    VertexActivity, VertexId,
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

#[allow(dead_code)]
mod percolation_support;
use percolation_support::{cases, Case};

const P_COUNT: usize = 17;
const REPLICAS: usize = 64;
const SEED: u64 = 0x5CA9_2A11;

#[derive(Clone, Copy)]
enum Mode {
    Site,
    Bond,
}

impl Mode {
    const ALL: [Self; 2] = [Self::Site, Self::Bond];

    const fn label(self) -> &'static str {
        match self {
            Self::Site => "site",
            Self::Bond => "bond",
        }
    }
}

fn probabilities() -> Vec<Probability> {
    (1..=P_COUNT)
        .map(|index| Probability::new(index as f64 / (P_COUNT + 1) as f64).unwrap())
        .collect()
}

fn benchmark_cases() -> Vec<Case> {
    cases()
        .into_iter()
        .filter(|case| matches!(case.name, "square" | "cubic" | "sparse-er"))
        .collect()
}

struct DirectOccupancy {
    mode: Mode,
    sites: Vec<bool>,
    bonds: Vec<bool>,
}

impl DirectOccupancy {
    fn new(lattice: &CsrLattice, mode: Mode) -> Self {
        Self {
            mode,
            sites: vec![matches!(mode, Mode::Bond); lattice.vertex_count()],
            bonds: vec![matches!(mode, Mode::Site); lattice.edge_count()],
        }
    }

    fn sample(&mut self, probability: f64, rng: &mut Xoshiro256PlusPlus) {
        match self.mode {
            Mode::Site => {
                for active in &mut self.sites {
                    *active = rng.random::<f64>() < probability;
                }
            }
            Mode::Bond => {
                for active in &mut self.bonds {
                    *active = rng.random::<f64>() < probability;
                }
            }
        }
    }
}

impl VertexActivity for DirectOccupancy {
    fn vertex_count(&self) -> usize {
        self.sites.len()
    }

    fn vertex_active(&self, vertex: VertexId) -> bool {
        self.sites[vertex.index()]
    }
}

impl EdgeActivity for DirectOccupancy {
    fn edge_count(&self) -> usize {
        self.bonds.len()
    }

    fn edge_active(&self, edge: EdgeId) -> bool {
        self.bonds[edge.index()]
    }
}

fn bench_scan(criterion: &mut Criterion) {
    for case in benchmark_cases() {
        let (vertices, edges) = (case.lattice.vertex_count(), case.lattice.edge_count());
        let mut workspace = NewmanZiffWorkspace::new();
        workspace.prepare(vertices, edges).unwrap();

        let mut group =
            criterion.benchmark_group(format!("nz_scan/{}/V={vertices}/E={edges}", case.name));
        for mode in Mode::ALL {
            let trials = match mode {
                Mode::Site => vertices,
                Mode::Bond => edges,
            };
            group.throughput(Throughput::Elements(trials as u64));
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(SEED);
            group.bench_function(BenchmarkId::new("scan", mode.label()), |bencher| {
                bencher.iter(|| {
                    let curve = match mode {
                        Mode::Site => {
                            cmc_rs::scan_site(&case.lattice, &mut workspace, &mut rng).unwrap()
                        }
                        Mode::Bond => {
                            cmc_rs::scan_bond(&case.lattice, &mut workspace, &mut rng).unwrap()
                        }
                    };
                    black_box(curve.points()[curve.trials()]);
                });
            });
        }
        group.finish();
    }
}

fn bench_reweight(criterion: &mut Criterion) {
    let probabilities = probabilities();
    for case in benchmark_cases() {
        let mut site_workspace = NewmanZiffWorkspace::new();
        let mut site_rng = Xoshiro256PlusPlus::seed_from_u64(SEED);
        let site_curve =
            cmc_rs::scan_site(&case.lattice, &mut site_workspace, &mut site_rng).unwrap();
        let mut bond_workspace = NewmanZiffWorkspace::new();
        let mut bond_rng = Xoshiro256PlusPlus::seed_from_u64(SEED);
        let bond_curve =
            cmc_rs::scan_bond(&case.lattice, &mut bond_workspace, &mut bond_rng).unwrap();

        let mut group = criterion.benchmark_group(format!(
            "nz_reweight/{}/V={}/E={}",
            case.name,
            case.lattice.vertex_count(),
            case.lattice.edge_count()
        ));
        group.throughput(Throughput::Elements(P_COUNT as u64));
        for mode in Mode::ALL {
            let curve = match mode {
                Mode::Site => &site_curve,
                Mode::Bond => &bond_curve,
            };
            let mut reweight_workspace = ReweightWorkspace::new();
            group.bench_function(
                BenchmarkId::new(format!("{P_COUNT}p"), mode.label()),
                |bencher| {
                    bencher.iter(|| {
                        for &probability in &probabilities {
                            black_box(
                                reweight(
                                    curve,
                                    probability,
                                    ScanObservable::LargestComponentSize,
                                    &mut reweight_workspace,
                                )
                                .unwrap(),
                            );
                        }
                    });
                },
            );
        }
        group.finish();
    }
}

fn bench_direct(criterion: &mut Criterion) {
    let probabilities = probabilities();
    for case in benchmark_cases() {
        let mut workspace = ComponentWorkspace::new();
        let mut group = criterion.benchmark_group(format!(
            "nz_direct/{}/V={}/E={}/{P_COUNT}p/{REPLICAS}r",
            case.name,
            case.lattice.vertex_count(),
            case.lattice.edge_count()
        ));
        group.throughput(Throughput::Elements((P_COUNT * REPLICAS) as u64));
        for mode in Mode::ALL {
            let mut occupancy = DirectOccupancy::new(&case.lattice, mode);
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(SEED);
            group.bench_function(BenchmarkId::new("direct", mode.label()), |bencher| {
                bencher.iter(|| {
                    for _ in 0..REPLICAS {
                        for &probability in &probabilities {
                            occupancy.sample(probability.get(), &mut rng);
                            let summary =
                                analyze(&case.lattice, &occupancy, &occupancy, &[], &mut workspace)
                                    .unwrap()
                                    .summary();
                            black_box(summary.largest_component_size);
                        }
                    }
                });
            });
        }
        group.finish();
    }
}

fn bench_combined(criterion: &mut Criterion) {
    let probabilities = probabilities();
    for case in benchmark_cases() {
        let mut scan_workspace = NewmanZiffWorkspace::new();
        scan_workspace
            .prepare(case.lattice.vertex_count(), case.lattice.edge_count())
            .unwrap();
        let mut reweight_workspace = ReweightWorkspace::new();
        let mut group = criterion.benchmark_group(format!(
            "nz_combined/{}/V={}/E={}/{P_COUNT}p/{REPLICAS}r",
            case.name,
            case.lattice.vertex_count(),
            case.lattice.edge_count()
        ));
        group.throughput(Throughput::Elements((P_COUNT * REPLICAS) as u64));
        for mode in Mode::ALL {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(SEED);
            group.bench_function(BenchmarkId::new("scan+reweight", mode.label()), |bencher| {
                bencher.iter(|| {
                    for _ in 0..REPLICAS {
                        let curve = match mode {
                            Mode::Site => {
                                cmc_rs::scan_site(&case.lattice, &mut scan_workspace, &mut rng)
                                    .unwrap()
                            }
                            Mode::Bond => {
                                cmc_rs::scan_bond(&case.lattice, &mut scan_workspace, &mut rng)
                                    .unwrap()
                            }
                        };
                        for &probability in &probabilities {
                            black_box(
                                reweight(
                                    &curve,
                                    probability,
                                    ScanObservable::LargestComponentSize,
                                    &mut reweight_workspace,
                                )
                                .unwrap(),
                            );
                        }
                    }
                });
            });
        }
        group.finish();
    }
}

criterion_group! {
    name = newman_ziff;
    config = Criterion::default()
        .sample_size(10)
        .warm_up_time(Duration::from_millis(100))
        .measurement_time(Duration::from_millis(250));
    targets = bench_scan, bench_reweight, bench_direct, bench_combined
}
criterion_main!(newman_ziff);
