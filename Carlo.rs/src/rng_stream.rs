//! Deterministic, domain-separated random-number stream derivation.
//!
//! A stream is identified by logical simulation coordinates rather than by
//! arithmetic seed offsets or execution order. This keeps independent chains,
//! replicas, backend tasks and phases reproducible when scheduling changes.

use rand_core::SeedableRng;

/// Logical lifecycle domain used when deriving a stream.
///
/// Distinct phases ensure that RNG draws made during thermalization,
/// measurement, exchange, and checkpointing never collide, even when
/// task counts or thread assignments change between runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u64)]
pub enum RngPhase {
    /// Model construction or restored state before the first scheduled sweep.
    #[default]
    Initialization = 0,
    /// Warmup sweeps.
    Thermalization = 1,
    /// Production sweeps with measurement accumulation.
    Measurement = 2,
    /// Parallel-tempering replica exchange.
    Exchange = 3,
    /// Checkpoint write/read.
    Checkpoint = 4,
    /// Per-task streams spawned by a parallel backend.
    BackendTask = 5,
    /// Post-run finalization.
    Finished = 6,
}

/// Complete identity of a deterministic random-number stream.
///
/// Fields are combined through domain-separated SplitMix64 rounds in
/// [`seed()`](RngStreamKey::seed) so that changing any single field
/// produces a statistically independent stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RngStreamKey {
    /// Global seed shared by all tasks in the job.
    pub base_seed: u64,
    /// Parameter-set index within a job.
    pub task_id: u64,
    /// Repetition index within a task.
    pub run_id: u64,
    /// Parallel-tempering chain index.
    pub chain_id: u64,
    /// Replica index for independent replicas within a chain.
    pub replica_id: u64,
    /// Physical thread / worker index.
    pub thread_id: u64,
    /// Lifecycle phase of the simulation.
    pub phase: RngPhase,
    /// Sub-stream counter for fine-grained stream splitting.
    pub substream: u64,
}

impl RngStreamKey {
    /// Create a key with only `base_seed` set; all other fields default to zero.
    #[inline]
    pub const fn new(base_seed: u64) -> Self {
        Self {
            base_seed,
            task_id: 0,
            run_id: 0,
            chain_id: 0,
            replica_id: 0,
            thread_id: 0,
            phase: RngPhase::Initialization,
            substream: 0,
        }
    }

    /// Set the task ID.
    #[inline]
    pub const fn with_task(mut self, value: u64) -> Self {
        self.task_id = value;
        self
    }

    /// Set the run ID.
    #[inline]
    pub const fn with_run(mut self, value: u64) -> Self {
        self.run_id = value;
        self
    }

    /// Set the parallel-tempering chain ID.
    #[inline]
    pub const fn with_chain(mut self, value: u64) -> Self {
        self.chain_id = value;
        self
    }

    /// Set the replica ID.
    #[inline]
    pub const fn with_replica(mut self, value: u64) -> Self {
        self.replica_id = value;
        self
    }

    /// Set the thread ID.
    #[inline]
    pub const fn with_thread(mut self, value: u64) -> Self {
        self.thread_id = value;
        self
    }

    /// Set the lifecycle phase.
    #[inline]
    pub const fn with_phase(mut self, value: RngPhase) -> Self {
        self.phase = value;
        self
    }

    /// Set the sub-stream counter.
    #[inline]
    pub const fn with_substream(mut self, value: u64) -> Self {
        self.substream = value;
        self
    }

    /// Derive a stable 64-bit seed using domain-separated SplitMix64 rounds.
    #[inline]
    pub fn seed(self) -> u64 {
        let mut state = mix64(self.base_seed ^ 0x5343_5554_544C_4552);
        state = fold(state, 0x5441_534B, self.task_id);
        state = fold(state, 0x5255_4E00, self.run_id);
        state = fold(state, 0x4348_4149_4E00, self.chain_id);
        state = fold(state, 0x5245_504C_4943_4100, self.replica_id);
        state = fold(state, 0x5448_5245_4144, self.thread_id);
        state = fold(state, 0x5048_4153_4500, self.phase as u64);
        fold(state, 0x5355_4253_5452, self.substream)
    }

    #[inline]
    pub fn seeded<R: SeedableRng>(self) -> R {
        R::seed_from_u64(self.seed())
    }
}

#[inline]
fn fold(state: u64, domain: u64, value: u64) -> u64 {
    mix64(state ^ domain ^ mix64(value.wrapping_add(domain.rotate_left(17))))
}

#[inline]
fn mix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use rand_xoshiro::Xoshiro256PlusPlus;
    use std::collections::HashSet;

    #[test]
    fn stream_derivation_is_deterministic() {
        let key = RngStreamKey::new(42)
            .with_task(3)
            .with_run(7)
            .with_chain(2)
            .with_replica(5)
            .with_thread(1)
            .with_phase(RngPhase::Measurement)
            .with_substream(11);
        assert_eq!(key.seed(), key.seed());
        let mut left: Xoshiro256PlusPlus = key.seeded();
        let mut right: Xoshiro256PlusPlus = key.seeded();
        assert_eq!(left.next_u64(), right.next_u64());
    }

    #[test]
    fn every_identity_component_domain_separates_streams() {
        let base = RngStreamKey::new(123);
        let seeds = [
            base.seed(),
            base.with_task(1).seed(),
            base.with_run(1).seed(),
            base.with_chain(1).seed(),
            base.with_replica(1).seed(),
            base.with_thread(1).seed(),
            base.with_phase(RngPhase::Measurement).seed(),
            base.with_substream(1).seed(),
        ];
        assert_eq!(seeds.into_iter().collect::<HashSet<_>>().len(), seeds.len());
    }
}
