//! MPI Backend Tests
//!
//! These tests require MPI to be installed and must be run with `mpirun`.
//!
//! # Running MPI tests
//!
//! ```bash
//! # Install MPI (Ubuntu/Debian)
//! sudo apt-get install libopenmpi-dev openmpi-bin
//!
//! # Build with MPI feature
//! cargo build --features mpi
//!
//! # Run tests with mpirun
//! mpirun -np 4 cargo test --features mpi --test mpi_test
//! ```

#[cfg(feature = "mpi")]
mod mpi_tests {
    use carlo_rs::{Backend, CarloError, Context, FromParams, MonteCarlo, MpiBackend, Params};
    use rand_xoshiro::Xoshiro256PlusPlus;

    /// Simple test MC for MPI
    #[allow(dead_code)]
    struct TestMC {
        sweep_count: u64,
    }

    impl MonteCarlo for TestMC {
        type Rng = Xoshiro256PlusPlus;

        fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
            self.sweep_count += 1;
            if ctx.is_thermalized() {
                ctx.measure("SweepCount", self.sweep_count as f64);
            }
        }
    }

    impl FromParams for TestMC {
        fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
            Ok(Self { sweep_count: 0 })
        }
    }

    #[test]
    fn test_mpi_backend_creation() {
        let backend = MpiBackend::new();
        // This test requires MPI to be initialized via mpirun
        if let Ok(b) = backend {
            assert!(b.size() >= 2, "MPI requires at least 2 ranks");
            if b.is_controller() {
                assert_eq!(b.rank(), 0);
            }
        }
    }

    #[test]
    fn test_mpi_barrier() {
        if let Ok(backend) = MpiBackend::new() {
            // Test that barrier works
            backend.barrier();
        }
    }

    #[test]
    fn test_mpi_communicator_split() {
        if let Ok(backend) = MpiBackend::with_ranks_per_run(2) {
            // Test ranks per run configuration
            assert!((backend.size() - 1) % 2 == 0);
        }
    }
}

#[cfg(not(feature = "mpi"))]
mod no_mpi_tests {
    #[test]
    fn test_mpi_feature_not_enabled() {
        // When MPI feature is not enabled, we can't test MPI backend
        // This test just verifies the code compiles without MPI
    }
}
