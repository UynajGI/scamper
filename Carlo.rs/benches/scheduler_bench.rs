use carlo_rs::{
    CarloError, Context, FromParams, MonteCarlo, Params, RayonBackend, RunConfig, Scheduler,
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand_xoshiro::Xoshiro256PlusPlus;

/// A simple MC for benchmarking - measures sweep count
struct BenchMC {
    sweep_count: u64,
}

impl MonteCarlo for BenchMC {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.sweep_count += 1;
        if ctx.is_thermalized() {
            ctx.measure("SweepCount", self.sweep_count as f64);
        }
    }
}

impl FromParams for BenchMC {
    fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Ok(Self { sweep_count: 0 })
    }
}

fn bench_scheduler(c: &mut Criterion) {
    let mut group = c.benchmark_group("scheduler");

    for sweeps in [100, 1000, 10000].iter() {
        group.bench_with_input(
            BenchmarkId::new("single_task", sweeps),
            sweeps,
            |b, &sweeps| {
                let config = RunConfig {
                    thermalization_sweeps: sweeps / 10,
                    measurement_sweeps: sweeps,
                    binsize: (sweeps / 10) as usize,
                    base_seed: 42,
                    ..Default::default()
                };
                let backend = RayonBackend::new(1);
                let params = Params::new();

                b.iter(|| {
                    Scheduler::new(backend.clone(), config.clone()).run_one::<BenchMC>(&params)
                });
            },
        );
    }

    group.finish();
}

fn bench_parallel(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallel");

    for n_tasks in [1, 2, 4, 8].iter() {
        group.bench_with_input(
            BenchmarkId::new("parallel_tasks", n_tasks),
            n_tasks,
            |b, &n_tasks| {
                let config = RunConfig {
                    thermalization_sweeps: 10,
                    measurement_sweeps: 100,
                    binsize: 10,
                    base_seed: 42,
                    ..Default::default()
                };
                let backend = RayonBackend::new(4);
                let params = Params::new();

                b.iter(|| {
                    Scheduler::new(backend.clone(), config.clone())
                        .run_parallel::<BenchMC>(n_tasks, &params)
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_scheduler, bench_parallel);
criterion_main!(benches);
