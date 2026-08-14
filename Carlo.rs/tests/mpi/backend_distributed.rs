//! Exactness guarantees for the MPI backend (`MpiBackend`) across rank
//! counts, going beyond the smoke-level checks in `mpi/mpi_test.rs` and the
//! single-task scheduler run in `mpi/distributed.rs`.
//!
//! MPI can be initialized exactly once per process, so this file contains a
//! single `#[test]`. Constructing `MpiBackend` performs that one
//! initialization (the rsmpi `Universe` stays alive inside the backend
//! runtime for the whole test); the world size then selects the scenarios:
//!
//! - **Degenerate singleton (np 1)**: rank 0 is a controller with zero
//!   workers; the complete task set still executes on the controller, in
//!   task order, exactly once, and `barrier()` returns.
//! - **Controller/worker exactness (np >= 2)**: `3*size + 7` independent
//!   tasks (deliberately not divisible by the rank count) are dispatched;
//!   every rank asserts it executed exactly the tasks routed to it, and each
//!   rank (the controller included) aggregates the per-rank partial results
//!   and asserts multiset equality over task ids (no lost, duplicated, or
//!   misrouted tasks), an aggregate count of exactly N, and an aggregate
//!   value equal to the exact expected sum.
//! - **RNG stream independence**: each task's first draw matches the
//!   documented `RngStreamKey` derivation (base seed, task id, replica =
//!   world rank, `RngPhase::BackendTask`) and first draws are pairwise
//!   distinct across all tasks and ranks.
//! - **Dispatch determinism**: an identical seed reproduces the identical
//!   per-rank task/draw sequence, and a different base seed changes every
//!   task's first draw.
//! - **Topology invariants**: with `ranks_per_run == 1`, the rank ->
//!   run-group mapping is a bijection, every rank leads its own run, and
//!   `num_parallel_runs() == num_workers()`.
//! - **Barrier semantics**: markers written before a barrier are visible to
//!   every rank after it, and barriers can be repeated without deadlocking.
//! - **Single-initialization ownership**: once this backend owns the
//!   process's MPI lifetime, a second `MpiBackend::new()` and the
//!   controller/worker scheduler entry point `run_distributed()` both fail
//!   with the documented "already initialized by another owner" error
//!   instead of hanging or corrupting MPI state.
//!
//! Scenarios that need more ranks than available self-skip with a note on
//! stderr; the commands below all perform meaningful work:
//!
//! ```bash
//! mpirun -np 1 cargo test --features mpi --test suite -- --ignored --exact mpi_backend_distributed::mpi_backend_distributed_suite --nocapture
//! mpirun -np 2 cargo test --features mpi --test suite -- --ignored --exact mpi_backend_distributed::mpi_backend_distributed_suite --nocapture
//! mpirun -np 4 cargo test --features mpi --test suite -- --ignored --exact mpi_backend_distributed::mpi_backend_distributed_suite --nocapture
//! ```
//!
//! Cross-rank aggregation uses per-rank files inside a shared scratch
//! directory (overridable via `CARLO_MPI_BACKEND_TEST_DIR`); the backend's
//! own `barrier()` orders writes against reads, so no sleeps are needed.

#![cfg(feature = "mpi")]

use carlo_rs::{
    run_distributed, Backend, CarloError, Context, DistributedConfig, FromParams, MonteCarlo,
    MpiBackend, Params, RngPhase, RngStreamKey, RunConfig, TaskSpec,
};
use rand_core::Rng as _;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

const BASE_SEED: u64 = 0x5EED_C0DE_1234;
const ALT_SEED: u64 = 0x0BAD_1DEA_9876;

/// Minimal model used only to type the `run_distributed` call in the
/// single-initialization check; it never executes a sweep.
struct NoopMc;

impl MonteCarlo for NoopMc {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {}
}

impl FromParams for NoopMc {
    fn from_params(_params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Ok(Self)
    }
}

fn scratch_dir() -> PathBuf {
    std::env::var_os("CARLO_MPI_BACKEND_TEST_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("carlo-rs-mpi-backend-distributed-test"))
}

/// Tasks the backend must dispatch to `rank` out of `n_tasks` total tasks.
fn expected_dispatch(rank: i32, size: i32, n_tasks: usize) -> Vec<usize> {
    (0..n_tasks)
        .filter(|&task| task % size as usize == rank as usize)
        .collect()
}

/// Run `n_tasks` through the backend, recording `(task_id, first_draw)` per
/// executed task in dispatch order.
fn spawn_and_collect(backend: &MpiBackend, n_tasks: usize, base_seed: u64) -> Vec<(usize, u64)> {
    let records: Mutex<Vec<(usize, u64)>> = Mutex::new(Vec::new());
    backend.spawn_tasks(n_tasks, base_seed, |task_id, rng| {
        let draw = rng.next_u64();
        records
            .lock()
            .expect("task record lock poisoned")
            .push((task_id, draw));
    });
    records.into_inner().expect("task record lock poisoned")
}

/// The documented first draw of a backend task stream (see
/// `MpiBackend::spawn_tasks`): a domain-separated key over the base seed,
/// the task id, replica = world rank, and `RngPhase::BackendTask`.
fn expected_first_draw(base_seed: u64, task_id: usize, replica: i32) -> u64 {
    let mut rng: <MpiBackend as Backend>::Rng = RngStreamKey::new(base_seed)
        .with_task(task_id as u64)
        .with_replica(replica as u64)
        .with_phase(RngPhase::BackendTask)
        .seeded();
    rng.next_u64()
}

fn write_scratch(path: &Path, contents: &str) {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).expect("scratch directory creation failed");
    }
    fs::write(path, contents).expect("scratch file write failed");
}

fn read_scratch(path: &Path) -> String {
    fs::read_to_string(path).expect("scratch file read failed")
}

/// Degenerate singleton: rank 0 is a controller with zero workers, yet the
/// complete workload executes on it, once per task, in task order.
fn singleton_controller_completes_everything(backend: &MpiBackend) {
    assert!(backend.is_controller());
    assert_eq!(backend.num_workers(), 0, "np 1 has no worker ranks");
    assert_eq!(backend.num_parallel_runs(), 0);
    assert!(backend.is_run_leader());
    assert_eq!(backend.run_group(), 0);
    assert_eq!(backend.rank_in_run(), 0);

    let n_tasks = 10;
    let records = spawn_and_collect(backend, n_tasks, BASE_SEED);

    let observed: Vec<usize> = records.iter().map(|(task, _)| *task).collect();
    let expected: Vec<usize> = (0..n_tasks).collect();
    assert_eq!(
        observed, expected,
        "at np 1 every task must execute on the controller in task order"
    );

    let mut draws = HashSet::new();
    for (task_id, draw) in records {
        assert_eq!(draw, expected_first_draw(BASE_SEED, task_id, 0));
        assert!(draws.insert(draw), "task {task_id} stream must be unique");
    }
    assert_eq!(draws.len(), n_tasks);

    // A world barrier over a singleton must return immediately.
    backend.barrier();
}

/// Controller/worker exactness: `N = 3*size + 7` tasks (not divisible by the
/// rank count) are dispatched once each; every rank asserts it executed
/// exactly its dispatched tasks, and each rank (the controller included)
/// aggregates all per-rank partial results and checks the multiset, count,
/// and exact sum.
fn controller_receives_every_task_exactly_once(backend: &MpiBackend, dir: &Path) {
    let rank = backend.rank();
    let size = backend.size();
    let n_tasks = 3 * size as usize + 7;
    let records = spawn_and_collect(backend, n_tasks, BASE_SEED);

    // Worker invariant: only tasks routed to this rank executed, in
    // ascending dispatch order, and this rank did real work.
    let expected = expected_dispatch(rank, size, n_tasks);
    assert!(!expected.is_empty(), "every rank must receive work");
    let observed: Vec<usize> = records.iter().map(|(task, _)| *task).collect();
    assert_eq!(
        observed, expected,
        "rank {rank} must execute exactly the tasks dispatched to it"
    );

    // Per-task result record: task id, first stream draw, and the value the
    // task produces (v(t) = 2t + 1, so the exact aggregate is N^2).
    let mut lines = Vec::new();
    for (task_id, draw) in records {
        assert_eq!(
            draw,
            expected_first_draw(BASE_SEED, task_id, rank),
            "task {task_id} on rank {rank} must use the documented stream key"
        );
        lines.push(format!("{task_id} {draw} {}", 2 * task_id as u64 + 1));
    }
    write_scratch(
        &dir.join("parts").join(format!("r{rank}.part")),
        &lines.join("\n"),
    );
    backend.barrier();

    // Every rank aggregates the complete set of per-rank partial results
    // (all of them are visible once the barrier above returns). Rank 0 is
    // the controller that must receive every task exactly once; having each
    // rank verify the same global multiset keeps a failure symmetric, so a
    // broken invariant fails the test instead of deadlocking MPI_Finalize
    // against a peer that is still waiting inside a collective.
    let mut task_ids = Vec::new();
    let mut draws = HashSet::new();
    let mut total: u64 = 0;
    for source in 0..size {
        let part = read_scratch(&dir.join("parts").join(format!("r{source}.part")));
        for line in part.lines() {
            let mut fields = line.split_whitespace();
            let task_id: usize = fields
                .next()
                .expect("task id field")
                .parse()
                .expect("task id is an integer");
            let draw: u64 = fields
                .next()
                .expect("draw field")
                .parse()
                .expect("draw is an integer");
            let value: u64 = fields
                .next()
                .expect("value field")
                .parse()
                .expect("value is an integer");
            task_ids.push(task_id);
            draws.insert(draw);
            total += value;
        }
    }
    // Multiset equality: no lost, duplicated, or misrouted task results.
    task_ids.sort_unstable();
    let all_tasks: Vec<usize> = (0..n_tasks).collect();
    assert_eq!(task_ids, all_tasks, "every task result exactly once");
    assert_eq!(task_ids.len(), n_tasks, "aggregate count matches N");
    let n = n_tasks as u64;
    assert_eq!(
        total,
        n * n,
        "aggregated value equals the exact expected sum"
    );
    assert_eq!(
        draws.len(),
        n_tasks,
        "no two tasks or ranks may share an RNG stream"
    );
    backend.barrier();
}

/// An identical seed reproduces the identical per-rank dispatch and draws,
/// and a different base seed changes every task's first draw.
fn dispatch_is_reproducible_and_seed_sensitive(backend: &MpiBackend) {
    let n_tasks = 2 * backend.size() as usize;

    let first = spawn_and_collect(backend, n_tasks, BASE_SEED);
    let repeat = spawn_and_collect(backend, n_tasks, BASE_SEED);
    assert_eq!(
        first, repeat,
        "the same base seed must reproduce the identical task/draw sequence"
    );

    let alt = spawn_and_collect(backend, n_tasks, ALT_SEED);
    assert_eq!(first.len(), alt.len());
    for ((task_a, draw_a), (task_b, draw_b)) in first.iter().zip(&alt) {
        assert_eq!(task_a, task_b, "task routing must not depend on the seed");
        assert_ne!(
            draw_a, draw_b,
            "task {task_a} stream must depend on the base seed"
        );
    }
}

/// With `ranks_per_run == 1`, every rank forms its own run group: the rank
/// -> run-group mapping is a bijection and all ranks are run leaders.
fn run_group_partition_is_bijective(backend: &MpiBackend, dir: &Path) {
    let rank = backend.rank();
    let size = backend.size();
    write_scratch(
        &dir.join("topology").join(format!("r{rank}.topo")),
        &format!(
            "{} {} {}\n",
            backend.run_group(),
            backend.rank_in_run(),
            i32::from(backend.is_run_leader())
        ),
    );
    backend.barrier();

    let mut groups = Vec::new();
    let mut leader_count = 0;
    for source in 0..size {
        let line = read_scratch(&dir.join("topology").join(format!("r{source}.topo")));
        let mut fields = line.split_whitespace();
        let group: i32 = fields
            .next()
            .expect("group field")
            .parse()
            .expect("group is an integer");
        let position: i32 = fields
            .next()
            .expect("position field")
            .parse()
            .expect("position is an integer");
        let leader: i32 = fields
            .next()
            .expect("leader field")
            .parse()
            .expect("leader is a flag");
        groups.push(group);
        assert_eq!(
            position, 0,
            "ranks_per_run == 1 puts every rank at position 0"
        );
        if leader == 1 {
            leader_count += 1;
        }
    }
    groups.sort_unstable();
    let all_ranks: Vec<i32> = (0..size).collect();
    assert_eq!(groups, all_ranks, "rank -> run_group must be a bijection");
    assert_eq!(leader_count, size);
    assert_eq!(
        backend.num_parallel_runs(),
        backend.num_workers(),
        "ranks_per_run == 1 must yield one parallel run per worker"
    );
    backend.barrier();
}

/// A barrier orders visibility: markers written before the barrier exist on
/// every rank afterwards, and barriers may be repeated.
fn barrier_orders_visibility_and_repeats(backend: &MpiBackend, dir: &Path) {
    let rank = backend.rank();
    let size = backend.size();
    write_scratch(
        &dir.join("barrier").join(format!("arrived-r{rank}.marker")),
        &rank.to_string(),
    );
    backend.barrier();

    for source in 0..size {
        assert!(
            dir.join("barrier")
                .join(format!("arrived-r{source}.marker"))
                .exists(),
            "barrier must order rank {source}'s writes before every reader"
        );
    }
    backend.barrier();
}

/// rsmpi allows exactly one MPI initialization per process. Once this
/// backend owns it, neither a second backend nor the controller/worker
/// scheduler entry point may initialize again; both must fail with the
/// documented error instead of hanging or corrupting MPI state.
fn single_initialization_is_enforced() {
    let second = match MpiBackend::new() {
        Ok(_) => panic!("a second backend must not re-initialize MPI"),
        Err(error) => error,
    };
    assert!(
        second.to_string().contains("already initialized"),
        "unexpected error for double initialization: {second}"
    );

    let config = DistributedConfig {
        run_config: RunConfig {
            thermalization_sweeps: 0,
            measurement_sweeps: 1,
            binsize: 1,
            base_seed: BASE_SEED,
            progress_interval: 1,
            checkpoint_interval: 0,
        },
        ranks_per_run: 1,
        run_time: None,
        checkpoint_time: None,
        job_dir: PathBuf::from("."),
        tasks: vec![TaskSpec {
            id: 0,
            target_sweeps: 1,
            thermalization: 0,
            params: Params::new(),
        }],
    };
    let error = match run_distributed::<NoopMc, Xoshiro256PlusPlus>(config) {
        Ok(_) => panic!("run_distributed must not re-initialize MPI next to a live backend"),
        Err(error) => error,
    };
    assert!(
        error.to_string().contains("already initialized"),
        "unexpected error from run_distributed: {error}"
    );
}

/// MPI cannot be initialized and finalized repeatedly inside one process, so
/// this file intentionally contains a single test. Constructing the backend
/// performs the process's single MPI initialization; the `Universe` is kept
/// alive by the backend runtime until the test ends.
#[test]
#[ignore = "requires mpirun"]
fn mpi_backend_distributed_suite() {
    let backend = MpiBackend::new().expect("MPI must be launched under mpirun/mpiexec");
    let rank = backend.rank();
    let size = backend.size();

    assert!(size >= 1);
    assert_eq!(backend.is_controller(), rank == 0);
    assert_eq!(backend.num_workers(), (size - 1).max(0));
    assert_eq!(backend.ranks_per_run(), 1);

    // `Backend: Clone` must share the existing MPI runtime, not re-init it.
    let clone = backend.clone();
    assert_eq!((clone.rank(), clone.size()), (rank, size));
    assert_eq!(clone.num_workers(), backend.num_workers());
    drop(clone);

    let dir = scratch_dir();
    if rank == 0 {
        let _ = fs::remove_dir_all(&dir);
        for sub in ["parts", "topology", "barrier"] {
            fs::create_dir_all(dir.join(sub)).expect("scratch directory creation failed");
        }
    }
    backend.barrier();

    if size == 1 {
        singleton_controller_completes_everything(&backend);
    } else {
        eprintln!("singleton scenario requires world size 1; skipping at np {size}");
    }
    if size >= 2 {
        controller_receives_every_task_exactly_once(&backend, &dir);
    } else {
        eprintln!(
            "controller/worker exactness scenario requires at least 2 ranks; skipping at np {size}"
        );
    }

    dispatch_is_reproducible_and_seed_sensitive(&backend);
    run_group_partition_is_bijective(&backend, &dir);
    barrier_orders_visibility_and_repeats(&backend, &dir);
    single_initialization_is_enforced();

    backend.barrier();
    if rank == 0 {
        let _ = fs::remove_dir_all(&dir);
    }
}
