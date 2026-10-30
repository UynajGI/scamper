//! Boolean/continuum percolation benchmark (P2).
//!
//! Records the two per-realization figures of the Boolean pipeline at unit
//! number density (`rho = 1`, so the box side is `N^(1/D)`), with the disk
//! filling factor `eta` as the density knob of the overlap criterion (mean
//! overlap degree `z = 2^D * eta` for monodisperse grains):
//!
//! - `realize`: complete sampling realizations per second — point-process
//!   draw, overlap construction (P1 at `2 r_max` plus the polydisperse
//!   contact filter), and the F3 component pass, through the production
//!   `BooleanPercolationMC` sweep+measure path including observable
//!   recording;
//! - `measure`: frozen-configuration analyses per second — the equilibrium
//!   measurement path (overlap construction plus F3, no resampling, no
//!   recording), the same shape as P1's `analyze` figure.
//!
//! Systems: 2D and 3D, N = 4096 and 16384, eta in {0.7, 1.1281 (the 2D disk
//! threshold bracket anchor, used as an operating point), 1.8}, monodisperse
//! and three-species polydisperse (`r0 * {0.5, 1, 1.5}` scaled to the same
//! eta). The polydisperse `measure` figure carries the per-candidate
//! contact-filter cost on top of P1 enumeration.
//!
//! Run with: cargo bench -p cmc-rs --bench boolean_percolation

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    analyze, overlap_pairs, BooleanPercolationMC, ComponentWorkspace, EdgeActivity, EdgeId,
    GraphView, ObservablePlan, OrthorhombicCell, PointProcess, Radii, StaticObservable,
    VertexActivity, VertexId,
};
use criterion::measurement::WallTime;
use criterion::{
    criterion_group, criterion_main, BenchmarkGroup, BenchmarkId, Criterion, Throughput,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// The 2D disk critical filling factor of Mertens & Moore (PRE 86, 061109
/// (2012)) — used here purely as a benchmark operating point near the
/// percolation transition, never as an exact threshold.
const DISK_ETA_C: f64 = 1.128_087_37;

struct AllGrains(usize);

impl VertexActivity for AllGrains {
    fn vertex_count(&self) -> usize {
        self.0
    }

    fn vertex_active(&self, _vertex: VertexId) -> bool {
        true
    }
}

struct AllOverlaps(usize);

impl EdgeActivity for AllOverlaps {
    fn edge_count(&self) -> usize {
        self.0
    }

    fn edge_active(&self, _edge: EdgeId) -> bool {
        true
    }
}

/// Unit-ball volume in `D` dimensions.
fn unit_ball_volume(dimension: usize) -> f64 {
    match dimension {
        2 => std::f64::consts::PI,
        3 => (4.0 / 3.0) * std::f64::consts::PI,
        _ => panic!("benchmark supports 2D and 3D"),
    }
}

/// Monodisperse radius at unit density for filling factor `eta`:
/// `eta = rho * V_D * r^D`.
fn monodisperse_radius(dimension: usize, eta: f64) -> f64 {
    (eta / unit_ball_volume(dimension)).powf(1.0 / dimension as f64)
}

/// Three-species radii `r0 * {0.5, 1, 1.5}` (cycling) scaled so the mean
/// covered fraction matches `eta`; per-particle contact distances stay
/// heterogeneous.
fn polydisperse_radii<const D: usize>(n: usize, eta: f64) -> Radii {
    let powers = [0.5_f64.powi(D as i32), 1.0, 1.5_f64.powi(D as i32)];
    let mean_power = powers.iter().sum::<f64>() / 3.0;
    let r0 = (eta / (unit_ball_volume(D) * mean_power)).powf(1.0 / D as f64);
    Radii::per_particle((0..n).map(|index| r0 * [0.5, 1.0, 1.5][index % 3]))
        .expect("positive species radii")
}

fn bench_dimension<const D: usize>(group: &mut BenchmarkGroup<'_, WallTime>, label: &str) {
    let plan = || {
        ObservablePlan::new(
            [
                StaticObservable::ComponentCount,
                StaticObservable::LargestSize,
                StaticObservable::GiantFraction,
            ],
            Vec::new(),
        )
    };
    for (n, eta) in [
        (4_096_usize, 0.7_f64),
        (4_096, DISK_ETA_C),
        (4_096, 1.8),
        (16_384, 0.7),
        (16_384, DISK_ETA_C),
        (16_384, 1.8),
    ] {
        let mut lengths = [0.0; D];
        lengths.fill((n as f64).powf(1.0 / D as f64));
        let cell = OrthorhombicCell::new(lengths).expect("benchmark cell is valid");
        for (dispersion, radii) in [
            (
                "monodisperse",
                Radii::monodisperse(monodisperse_radius(D, eta)).unwrap(),
            ),
            ("polydisperse", polydisperse_radii::<D>(n, eta)),
        ] {
            let system = format!("{label}_n{n}_eta{eta:.3}_{dispersion}");

            // One frozen fixture shared by the measure path; the realize
            // path resamples every iteration through the adapter.
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xB_05);
            let mut frozen = Vec::with_capacity(n * D);
            for _ in 0..n {
                for &length in &lengths {
                    frozen.push(rng.random_range(0.0..length));
                }
            }
            let pair_count = overlap_pairs(&frozen, &cell, &radii)
                .expect("valid system")
                .pair_count();

            group.throughput(Throughput::Elements(pair_count as u64));
            group.bench_with_input(
                BenchmarkId::new("measure", &system),
                &(cell, radii.clone(), frozen),
                |b, (cell, radii, frozen)| {
                    let mut workspace = ComponentWorkspace::new();
                    b.iter(|| {
                        let graph = overlap_pairs(frozen, cell, radii).expect("valid system");
                        let summary = analyze(
                            &graph,
                            &AllGrains(graph.vertex_count()),
                            &AllOverlaps(graph.pair_count()),
                            &[],
                            &mut workspace,
                        )
                        .expect("analysis succeeds")
                        .summary();
                        std::hint::black_box(summary);
                    });
                },
            );

            group.throughput(Throughput::Elements(n as u64));
            group.bench_with_input(
                BenchmarkId::new("realize", &system),
                &(cell, radii.clone(), plan()),
                |b, (cell, radii, plan)| {
                    let mut adapter = BooleanPercolationMC::<D>::sampling(
                        *cell,
                        radii.clone(),
                        PointProcess::binomial(n).expect("positive count"),
                        plan.clone(),
                    )
                    .expect("valid adapter");
                    let mut context = Context::new_with_binsize(
                        Xoshiro256PlusPlus::seed_from_u64(0xB_06),
                        0,
                        1 << 40,
                    );
                    b.iter(|| {
                        adapter.sweep(&mut context);
                        adapter.measure(&mut context);
                    });
                },
            );

            eprintln!(
                "BOOLEAN_BENCH dim={label} n={n} eta={eta:.4} {dispersion} pairs={pair_count} \
                 pairs_per_particle={:.2}",
                pair_count as f64 / n as f64
            );
        }
    }
}

fn bench_boolean_percolation(c: &mut Criterion) {
    let mut group = c.benchmark_group("boolean_percolation");
    bench_dimension::<2>(&mut group, "2d");
    bench_dimension::<3>(&mut group, "3d");
    group.finish();
}

criterion_group!(benches, bench_boolean_percolation);
criterion_main!(benches);
