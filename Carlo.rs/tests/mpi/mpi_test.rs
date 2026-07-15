//! MPI backend smoke test.
//!
//! Run this test binary under MPI so every rank executes the same test:
//!
//! ```bash
//! mpirun -np 4 cargo test --features mpi --test mpi_test -- --nocapture
//! ```

#[cfg(feature = "mpi")]
mod mpi_tests {
    use carlo_rs::{Backend, MpiBackend};
    use std::sync::Mutex;

    /// MPI generally cannot be initialized and finalized repeatedly in one
    /// process, so this integration test intentionally contains a single test.
    #[test]
    #[ignore = "requires mpirun"]
    fn mpi_backend_smoke_suite() {
        let backend = MpiBackend::new().expect("MPI must be launched under mpirun/mpiexec");
        let rank = backend.rank();
        let size = backend.size();

        assert!(size >= 1);
        assert_eq!(backend.is_controller(), rank == 0);
        assert_eq!(backend.num_workers(), (size - 1).max(0));
        assert_eq!(backend.run_group(), if rank == 0 { 0 } else { rank });
        assert_eq!(backend.rank_in_run(), 0);

        let n_tasks = size as usize * 3 + 1;
        let seen = Mutex::new(Vec::new());
        backend.spawn_tasks(n_tasks, 0x5eed, |task_id, _rng| {
            seen.lock().expect("task list lock poisoned").push(task_id);
        });

        let mut actual = seen.into_inner().expect("task list lock poisoned");
        actual.sort_unstable();
        let expected: Vec<usize> = (0..n_tasks)
            .filter(|task_id| task_id % size as usize == rank as usize)
            .collect();
        assert_eq!(actual, expected);

        backend.barrier();
    }
}

#[cfg(not(feature = "mpi"))]
#[test]
fn mpi_feature_disabled_stub_compiles() {}
