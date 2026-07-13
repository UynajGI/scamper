use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

use super::Backend;
use crate::{RngPhase, RngStreamKey};

/// Rayon-based parallel backend (Phase 1).
#[derive(Clone)]
pub struct RayonBackend {
    _n_threads: usize,
}

impl RayonBackend {
    pub fn new(n_threads: usize) -> Self {
        Self {
            _n_threads: n_threads,
        }
    }
}

impl Default for RayonBackend {
    fn default() -> Self {
        Self::new(rayon::current_num_threads())
    }
}

impl Backend for RayonBackend {
    type Rng = Xoshiro256PlusPlus;

    fn spawn_tasks<F>(&self, n_tasks: usize, base_seed: u64, f: F)
    where
        F: Fn(usize, &mut Self::Rng) + Sync,
    {
        (0..n_tasks).into_par_iter().for_each(|task_id| {
            let mut rng: Self::Rng = RngStreamKey::new(base_seed)
                .with_task(task_id as u64)
                // Deliberately omit the physical Rayon worker index: task streams
                // must not change when work stealing assigns a task elsewhere.
                .with_phase(RngPhase::BackendTask)
                .seeded();
            f(task_id, &mut rng);
        });
    }

    fn barrier(&self) {
        // Rayon automatically synchronizes after parallel iteration
        // No explicit barrier needed
    }
}
