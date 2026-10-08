mod percolation_support;

use std::time::Duration;

use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    analyze, BondBernoulli, BoundaryQuery, ComponentSummary, ComponentWorkspace, EdgeActivity,
    GraphView, MixedBernoulli, ObservablePlan, Probability, SiteBernoulli, StaticConfiguration,
    StaticLaw, StaticPercolationMC, UndirectedGraphView, VertexActivity,
};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use percolation_support::{cases, sample_and_analyze, P_BOND, P_SITE};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const PARITY_SAMPLES: usize = 128;
const TIMED_SAMPLE_SEED: u64 = 0x4633_5449_4d45;
const TIMED_INITIAL_SAMPLES: usize = 1;

#[derive(Clone, Copy)]
enum Law {
    Site,
    Bond,
    Mixed,
}

impl Law {
    const ALL: [Self; 3] = [Self::Site, Self::Bond, Self::Mixed];

    const fn label(self) -> &'static str {
        match self {
            Self::Site => "site",
            Self::Bond => "bond",
            Self::Mixed => "mixed",
        }
    }

    fn sample(
        self,
        graph: &impl UndirectedGraphView,
        configuration: &mut StaticConfiguration,
        rng: &mut Xoshiro256PlusPlus,
    ) {
        let site = Probability::new(P_SITE).unwrap();
        let bond = Probability::new(P_BOND).unwrap();
        match self {
            Self::Site => SiteBernoulli::new(site.into()).sample(graph, configuration, rng),
            Self::Bond => BondBernoulli::new(bond.into()).sample(graph, configuration, rng),
            Self::Mixed => {
                MixedBernoulli::new(site.into(), bond.into()).sample(graph, configuration, rng)
            }
        }
        .expect("benchmark dimensions match");
    }
}

#[derive(Debug, Eq, PartialEq)]
struct ReferenceProjection {
    summary: ComponentSummary,
    query_outcomes: Vec<bool>,
}

/// Bench-only flood-fill projection. It intentionally shares no production UF
/// or analyzer helper and may allocate fresh scratch for every call.
fn flood_fill_projection<G, V, E>(
    graph: &G,
    vertex_activity: &V,
    edge_activity: &E,
    queries: &[BoundaryQuery],
) -> ReferenceProjection
where
    G: UndirectedGraphView,
    V: VertexActivity,
    E: EdgeActivity,
{
    let mut active_edge_count = 0usize;
    for edge_index in 0..graph.edge_count() {
        let edge = graph.edge_id(edge_index).expect("dense edge ID");
        let [left, right] = graph.edge_endpoints(edge);
        active_edge_count += usize::from(
            edge_activity.edge_active(edge)
                && vertex_activity.vertex_active(left)
                && vertex_activity.vertex_active(right),
        );
    }

    let mut seen = vec![false; graph.vertex_count()];
    let mut stack = Vec::new();
    let mut component = Vec::new();
    let mut query_outcomes = vec![false; queries.len()];
    let mut summary = ComponentSummary {
        active_vertex_count: 0,
        active_edge_count,
        component_count: 0,
        largest_component: None,
        largest_component_size: 0,
        raw_second_moment: 0,
    };
    for seed_index in 0..graph.vertex_count() {
        let seed = graph.vertex_id(seed_index).expect("dense vertex ID");
        if seen[seed_index] || !vertex_activity.vertex_active(seed) {
            continue;
        }
        seen[seed_index] = true;
        stack.push(seed);
        component.clear();
        let mut identity = seed;
        while let Some(vertex) = stack.pop() {
            component.push(vertex);
            identity = identity.min(vertex);
            for incidence in graph.incidences(vertex) {
                if !edge_activity.edge_active(incidence.edge)
                    || !vertex_activity.vertex_active(incidence.neighbor)
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

        let size = component.len();
        summary.active_vertex_count += size;
        summary.component_count += 1;
        summary.raw_second_moment += (size as u128) * (size as u128);
        if size > summary.largest_component_size
            || (size == summary.largest_component_size
                && summary
                    .largest_component
                    .is_none_or(|current| identity < current))
        {
            summary.largest_component = Some(identity);
            summary.largest_component_size = size;
        }
        for (query_index, query) in queries.iter().enumerate() {
            if query_outcomes[query_index] {
                continue;
            }
            let touches_from = component.iter().any(|vertex| query.from().contains(vertex));
            let touches_to = component.iter().any(|vertex| query.to().contains(vertex));
            query_outcomes[query_index] = touches_from && touches_to;
        }
    }
    ReferenceProjection {
        summary,
        query_outcomes,
    }
}

fn assert_equivalent_samples(
    case: &percolation_support::Case,
    law: Law,
    queries: &[BoundaryQuery],
) {
    let mut configuration =
        StaticConfiguration::new(case.lattice.vertex_count(), case.lattice.edge_count());
    let mut workspace = ComponentWorkspace::new();
    workspace
        .prepare(case.lattice.vertex_count(), queries.len(), false)
        .unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4633_4551_5549);
    for sample in 0..PARITY_SAMPLES {
        law.sample(&case.lattice, &mut configuration, &mut rng);
        let production = analyze(
            &case.lattice,
            &configuration,
            &configuration,
            queries,
            &mut workspace,
        )
        .unwrap();
        let reference =
            flood_fill_projection(&case.lattice, &configuration, &configuration, queries);
        assert_eq!(
            production.summary(),
            reference.summary,
            "{}/{}/sample={sample} summary",
            case.name,
            law.label()
        );
        assert_eq!(
            production.query_outcomes(),
            reference.query_outcomes,
            "{}/{}/sample={sample} queries",
            case.name,
            law.label()
        );
    }
}

fn assert_timed_configuration_sequences_match(case: &percolation_support::Case, law: Law) {
    let mut production =
        StaticConfiguration::new(case.lattice.vertex_count(), case.lattice.edge_count());
    let mut reference =
        StaticConfiguration::new(case.lattice.vertex_count(), case.lattice.edge_count());
    let mut production_rng = Xoshiro256PlusPlus::seed_from_u64(TIMED_SAMPLE_SEED);
    let mut reference_rng = Xoshiro256PlusPlus::seed_from_u64(TIMED_SAMPLE_SEED);
    for sample in 0..TIMED_INITIAL_SAMPLES + PARITY_SAMPLES {
        law.sample(&case.lattice, &mut production, &mut production_rng);
        law.sample(&case.lattice, &mut reference, &mut reference_rng);
        for vertex in case.lattice.vertex_ids() {
            assert_eq!(
                production.vertex_active(vertex),
                reference.vertex_active(vertex),
                "{}/{}/sample={sample}/vertex={vertex} timed sequence",
                case.name,
                law.label()
            );
        }
        for edge_index in 0..case.lattice.edge_count() {
            let edge = case.lattice.edge_id(edge_index).expect("dense edge ID");
            assert_eq!(
                production.edge_active(edge),
                reference.edge_active(edge),
                "{}/{}/sample={sample}/edge={edge} timed sequence",
                case.name,
                law.label()
            );
        }
    }
}

fn bench_percolation_components(criterion: &mut Criterion) {
    for case in cases() {
        let query = BoundaryQuery::new(&case.lattice, &case.from, &case.to).unwrap();
        let queries = [query];
        for law in Law::ALL {
            assert_equivalent_samples(&case, law, &queries);
            assert_timed_configuration_sequences_match(&case, law);

            let mut configuration =
                StaticConfiguration::new(case.lattice.vertex_count(), case.lattice.edge_count());
            let mut workspace = ComponentWorkspace::new();
            workspace
                .prepare(case.lattice.vertex_count(), queries.len(), false)
                .unwrap();
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(TIMED_SAMPLE_SEED);
            for _ in 0..TIMED_INITIAL_SAMPLES {
                law.sample(&case.lattice, &mut configuration, &mut rng);
            }

            let mut group = criterion.benchmark_group(format!(
                "percolation_components/{}/V={}/E={}/p_site={}/p_bond={}",
                case.name,
                case.lattice.vertex_count(),
                case.lattice.edge_count(),
                P_SITE,
                P_BOND
            ));
            group.throughput(Throughput::Elements(case.lattice.edge_count() as u64));
            group.bench_function(BenchmarkId::new("sample+analyze", law.label()), |bencher| {
                bencher.iter(|| {
                    law.sample(&case.lattice, &mut configuration, &mut rng);
                    black_box(
                        analyze(
                            &case.lattice,
                            &configuration,
                            &configuration,
                            &queries,
                            &mut workspace,
                        )
                        .unwrap(),
                    );
                });
            });
            group.bench_function(BenchmarkId::new("analyze-only", law.label()), |bencher| {
                bencher.iter(|| {
                    black_box(
                        analyze(
                            &case.lattice,
                            &configuration,
                            &configuration,
                            &queries,
                            &mut workspace,
                        )
                        .unwrap(),
                    );
                });
            });

            let mut reference_configuration =
                StaticConfiguration::new(case.lattice.vertex_count(), case.lattice.edge_count());
            let mut reference_rng = Xoshiro256PlusPlus::seed_from_u64(TIMED_SAMPLE_SEED);
            for _ in 0..TIMED_INITIAL_SAMPLES {
                law.sample(
                    &case.lattice,
                    &mut reference_configuration,
                    &mut reference_rng,
                );
            }
            group.bench_function(
                BenchmarkId::new("equivalent-reference-sample+analyze", law.label()),
                |bencher| {
                    bencher.iter(|| {
                        law.sample(
                            &case.lattice,
                            &mut reference_configuration,
                            &mut reference_rng,
                        );
                        black_box(flood_fill_projection(
                            &case.lattice,
                            &reference_configuration,
                            &reference_configuration,
                            &queries,
                        ));
                    });
                },
            );

            if case.name == "square" {
                let mode = match law {
                    Law::Site => percolation_support::ReferenceMode::Site,
                    Law::Bond => percolation_support::ReferenceMode::Bond,
                    Law::Mixed => percolation_support::ReferenceMode::Mixed,
                };
                let mut occupancy = percolation_support::ReferenceConfiguration::new(&case, mode);
                let mut historical_rng = Xoshiro256PlusPlus::seed_from_u64(0x4630_5245_4645);
                group.bench_function(
                    BenchmarkId::new("historical-pr4-sample+analyze", law.label()),
                    |bencher| {
                        bencher.iter(|| {
                            black_box(sample_and_analyze(
                                &case,
                                &mut occupancy,
                                &mut historical_rng,
                            ))
                        });
                    },
                );
            }

            // F4 production adapter: sweep + measure through a Carlo.rs
            // Context with the full observable plan. The context has its own
            // RNG seeded like the core paths, so this is an independent
            // measurement of the adapter composition, not a sequence-locked
            // rerun of `sample+analyze`. The difference to `sample+analyze`
            // is the measure overhead: summary projection, boundary events,
            // and recording nine scalars into the binned accumulators.
            let law_static = match law {
                Law::Site => StaticLaw::site(P_SITE),
                Law::Bond => StaticLaw::bond(P_BOND),
                Law::Mixed => StaticLaw::mixed(P_SITE, P_BOND),
            }
            .expect("benchmark probability");
            let plan = ObservablePlan::all(vec![BoundaryQuery::new(
                &case.lattice,
                &case.from,
                &case.to,
            )
            .unwrap()]);
            let mut adapter =
                StaticPercolationMC::new(case.lattice.clone(), law_static, plan).unwrap();
            let mut mc_context = Context::new_with_binsize(
                Xoshiro256PlusPlus::seed_from_u64(TIMED_SAMPLE_SEED),
                0,
                1024,
            );
            for _ in 0..TIMED_INITIAL_SAMPLES {
                adapter.sweep(&mut mc_context);
                adapter.measure(&mut mc_context);
            }
            group.bench_function(
                BenchmarkId::new("static-mc-sweep+measure", law.label()),
                |bencher| {
                    bencher.iter(|| {
                        adapter.sweep(&mut mc_context);
                        adapter.measure(&mut mc_context);
                        black_box(&mut mc_context);
                    });
                },
            );
            group.finish();
        }
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .sample_size(10)
        .warm_up_time(Duration::from_millis(100))
        .measurement_time(Duration::from_millis(250));
    targets = bench_percolation_components
}
criterion_main!(benches);
