//! End-to-end controller/worker-group MPI test.
//!
//! Run separately from other MPI test binaries because an MPI process may only
//! own one initialization lifetime:
//!
//! ```bash
//! mpirun -np 4 cargo test --features mpi --test mpi_distributed_test -- --nocapture
//! ```

#[cfg(feature = "mpi")]
mod distributed_tests {
    use carlo_rs::{
        run_distributed, CarloError, Context, DistributedConfig, FromParams, MonteCarlo, Params,
        RunConfig, TaskSpec,
    };
    use rand_xoshiro::Xoshiro256PlusPlus;
    use std::path::PathBuf;

    struct CounterMc {
        value: u64,
    }

    impl MonteCarlo for CounterMc {
        type Rng = Xoshiro256PlusPlus;

        fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {
            self.value = self.value.saturating_add(1);
        }

        fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
            ctx.measure("Counter", self.value as f64);
        }
    }

    impl FromParams for CounterMc {
        fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
            Ok(Self { value: 0 })
        }
    }

    #[test]
    #[ignore = "requires mpirun"]
    fn dynamic_scheduler_completes_and_merges_runs() {
        let job_dir = std::env::var_os("CARLO_MPI_TEST_JOB_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|| std::env::temp_dir().join("carlo-rs-mpi-distributed-test"));

        let config = DistributedConfig {
            run_config: RunConfig {
                thermalization_sweeps: 2,
                measurement_sweeps: 12,
                binsize: 2,
                base_seed: 1234,
                progress_interval: 1,
                checkpoint_interval: 0,
            },
            ranks_per_run: 1,
            run_time: None,
            checkpoint_time: None,
            job_dir: job_dir.clone(),
            tasks: vec![TaskSpec {
                id: 11,
                target_sweeps: 12,
                thermalization: 2,
                params: Params::new(),
            }],
        };

        match run_distributed::<CounterMc, Xoshiro256PlusPlus>(config) {
            Ok(results) if results.is_empty() => {
                // Worker ranks intentionally return no aggregate.
            }
            Ok(results) => {
                assert_eq!(results.len(), 1);
                let counter = results[0]
                    .get("Counter")
                    .expect("controller must receive Counter");
                assert!(counter.n_bins > 0);
                assert!(counter.mean.is_finite());
                let _ = std::fs::remove_dir_all(job_dir);
            }
            Err(error) => {
                // A normal non-mpirun `cargo test --features mpi` has world
                // size 1 and should fail validation instead of hanging.
                assert!(error.to_string().contains("at least one worker rank"));
            }
        }
    }
}

#[cfg(not(feature = "mpi"))]
#[test]
fn mpi_distributed_feature_disabled_stub_compiles() {}
