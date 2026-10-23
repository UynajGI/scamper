use std::hint::black_box;
use std::time::Duration;

use cmc_rs::{
    analyze, reweight, reweight_heterogeneous, ComponentWorkspace, CsrLattice, EdgeActivity,
    EdgeId, GraphView, HetReweightWorkspace, NewmanZiffWorkspace, Probability, ProbabilityField,
    ReweightWorkspace, ScanObservable, UndirectedGraphView, VertexActivity, VertexId,
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

#[allow(dead_code)]
mod percolation_support;
use percolation_support::{cases, Case};

const P_COUNT: usize = 17;
const REPLICAS: usize = 64;
const HET_FIELD_COUNT: usize = 8;
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

    fn sample_heterogeneous(
        &mut self,
        probabilities: &[Probability],
        rng: &mut Xoshiro256PlusPlus,
    ) {
        match self.mode {
            Mode::Site => {
                for (active, &probability) in self.sites.iter_mut().zip(probabilities) {
                    *active = rng.random::<f64>() < probability.get();
                }
            }
            Mode::Bond => {
                for (active, &probability) in self.bonds.iter_mut().zip(probabilities) {
                    *active = rng.random::<f64>() < probability.get();
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

/// T2 heterogeneous comparison: one scan then reweighting at K graded fields
/// vs K independent direct heterogeneous sampling runs at those fields (the
/// N1-style path: per-element draws plus one component analysis each). The
/// fields stay inside the reweighter's usable envelope (narrow ramps; see
/// the `reweight_heterogeneous` variance caveat).
fn bench_heterogeneous(criterion: &mut Criterion) {
    for case in benchmark_cases() {
        let (vertices, edges) = (case.lattice.vertex_count(), case.lattice.edge_count());
        let mut group = criterion.benchmark_group(format!(
            "nz_het/{}/V={vertices}/E={edges}/{HET_FIELD_COUNT}f/{REPLICAS}r",
            case.name
        ));
        group.throughput(Throughput::Elements((HET_FIELD_COUNT * REPLICAS) as u64));
        for mode in Mode::ALL {
            let trials = match mode {
                Mode::Site => vertices,
                Mode::Bond => edges,
            };
            let fields: Vec<Vec<Probability>> = (0..HET_FIELD_COUNT)
                .map(|field_index| {
                    let center = 0.45 + 0.30 * field_index as f64 / (HET_FIELD_COUNT - 1) as f64;
                    (0..trials)
                        .map(|index| {
                            let position = index as f64 / trials.max(1) as f64 - 0.5;
                            Probability::new((center + 0.008 * position).clamp(0.01, 0.99))
                                .expect("bench probability")
                        })
                        .collect()
                })
                .collect();

            let mut scan_workspace = NewmanZiffWorkspace::new();
            scan_workspace.prepare(vertices, edges).unwrap();
            let mut reweight_workspace = HetReweightWorkspace::new();
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
                        for field in &fields {
                            black_box(
                                reweight_heterogeneous(
                                    &curve,
                                    &ProbabilityField::Borrowed(field),
                                    ScanObservable::LargestComponentSize,
                                    &mut reweight_workspace,
                                )
                                .unwrap(),
                            );
                        }
                    }
                });
            });

            let mut workspace = ComponentWorkspace::new();
            let mut occupancy = DirectOccupancy::new(&case.lattice, mode);
            let mut direct_rng = Xoshiro256PlusPlus::seed_from_u64(SEED);
            group.bench_function(BenchmarkId::new("direct", mode.label()), |bencher| {
                bencher.iter(|| {
                    for field in &fields {
                        for _ in 0..REPLICAS {
                            occupancy.sample_heterogeneous(field, &mut direct_rng);
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

criterion_group! {
    name = newman_ziff;
    config = Criterion::default()
        .sample_size(10)
        .warm_up_time(Duration::from_millis(100))
        .measurement_time(Duration::from_millis(250));
    targets = bench_scan, bench_reweight, bench_direct, bench_combined, bench_heterogeneous
}
criterion_main!(newman_ziff);
