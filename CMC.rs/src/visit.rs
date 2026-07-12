//! Reusable site-visitation schedules.

use rand::{Rng, RngExt};

/// Order in which local-update kernels visit sites.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum VisitSchedule {
    /// Deterministic contiguous traversal; maximizes cache locality.
    Sequential,
    /// One Fisher-Yates random permutation per sweep.
    #[default]
    RandomPermutation,
}

/// Allocation-free-after-warmup visitation workspace.
#[derive(Debug, Clone, Default)]
pub struct SiteOrder {
    sites: Vec<usize>,
    identity_order: bool,
}

impl SiteOrder {
    pub const fn new() -> Self {
        Self {
            sites: Vec::new(),
            identity_order: true,
        }
    }

    pub fn prepare(
        &mut self,
        n_sites: usize,
        schedule: VisitSchedule,
        rng: &mut impl Rng,
    ) -> &[usize] {
        if self.sites.len() != n_sites {
            self.sites.clear();
            self.sites.extend(0..n_sites);
            self.identity_order = true;
        }
        match schedule {
            VisitSchedule::Sequential => {
                if !self.identity_order {
                    for (site, value) in self.sites.iter_mut().enumerate() {
                        *value = site;
                    }
                    self.identity_order = true;
                }
            }
            VisitSchedule::RandomPermutation => {
                for index in (1..n_sites).rev() {
                    let swap_with = rng.random_range(0..=index);
                    self.sites.swap(index, swap_with);
                }
                self.identity_order = false;
            }
        }
        &self.sites
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn sequential_schedule_is_stable() {
        let mut order = SiteOrder::new();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
        let _ = order.prepare(5, VisitSchedule::RandomPermutation, &mut rng);
        assert_eq!(
            order.prepare(5, VisitSchedule::Sequential, &mut rng),
            &[0, 1, 2, 3, 4]
        );
    }
}
