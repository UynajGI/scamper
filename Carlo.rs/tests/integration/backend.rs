use carlo_rs::backend::{Backend, RayonBackend};
use rand_core::Rng;
use std::sync::atomic::{AtomicU64, Ordering};

#[test]
fn test_rayon_backend_spawn_tasks() {
    let backend = RayonBackend::new(4); // 4 threads
    let counter = AtomicU64::new(0);

    backend.spawn_tasks(10, 42, |_task_id, rng| {
        // Each task should have unique RNG seed
        let _seed_offset = rng.next_u64(); // Just read something to prove RNG works
        counter.fetch_add(1, Ordering::SeqCst);
    });

    backend.barrier();

    assert_eq!(counter.load(Ordering::SeqCst), 10);
}

/// Thread-count independence: same base_seed produces identical per-task
/// RNG streams regardless of how many Rayon threads are used. This is
/// because RngStreamKey derives from base_seed + task_id only, never
/// the physical worker index.
#[test]
fn thread_count_does_not_change_rng_streams() {
    use std::sync::{Arc, Mutex};

    fn collect_draws(n_threads: usize) -> Vec<(usize, [u64; 8])> {
        type TaskDraws = Vec<(usize, [u64; 8])>;
        let backend = RayonBackend::new(n_threads);
        let results: Arc<Mutex<TaskDraws>> = Arc::new(Mutex::new(Vec::new()));
        let r = results.clone();
        backend.spawn_tasks(8, 12345, move |task_id, rng| {
            let mut draws = [0u64; 8];
            for slot in &mut draws {
                *slot = rng.next_u64();
            }
            r.lock().unwrap().push((task_id, draws));
        });
        backend.barrier();
        let mut v = results.lock().unwrap().clone();
        v.sort_by_key(|(id, _)| *id);
        v
    }

    let v1 = collect_draws(1);
    let v4 = collect_draws(4);

    assert_eq!(v1.len(), 8);
    assert_eq!(v4.len(), 8);

    for (i, ((id1, draws1), (id4, draws4))) in v1.iter().zip(v4.iter()).enumerate() {
        assert_eq!(*id1, i, "task_id should match index");
        assert_eq!(*id4, i);
        assert_eq!(
            draws1, draws4,
            "RNG stream for task {i} differs between 1-thread and 4-thread"
        );
    }
}
