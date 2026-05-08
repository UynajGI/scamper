use rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

use super::Backend;

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
        // Simple seed offset strategy (default)
        // For strict-repro, use jump sequence (Phase 2)
        (0..n_tasks).into_par_iter().for_each(|task_id| {
            let seed = base_seed.wrapping_add(task_id as u64);
            let mut rng = Self::Rng::seed_from_u64(seed);
            f(task_id, &mut rng);
        });
    }

    fn barrier(&self) {
        // Rayon automatically synchronizes after parallel iteration
        // No explicit barrier needed
    }
}
