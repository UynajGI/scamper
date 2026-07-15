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
