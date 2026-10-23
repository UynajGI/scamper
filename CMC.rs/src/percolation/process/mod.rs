//! E track: cascade-process primitives and adapters.
//!
//! E1 shares allocation, queue, and counter primitives ONLY — there is
//! deliberately no common "physical transition type" (roadmap §12): each
//! process (E2 bootstrap, E3 k-core) defines its own state machine and
//! borrows these reusable buffers. Everything is prepared once and reused
//! across realizations with zero steady-state allocation.

mod bootstrap;
mod incremental;
mod kcore;

use std::collections::TryReserveError;

pub use bootstrap::{
    BootstrapDynamics, BootstrapError, BootstrapOutcome, BootstrapPercolationMC, BootstrapSeeds,
    BootstrapThreshold, BootstrapWorkspace,
};
pub use incremental::{
    IncrementalError, IncrementalSummary, IncrementalUnionWorkspace, ProcessTrajectory,
};
pub use kcore::{
    kcore_decompose, CoreThresholds, KCoreError, KCoreMC, KCoreSummary, KCoreWorkspace,
};

/// Two-phase FIFO work queue over dense vertex indices.
///
/// `pop` advances a head cursor (elements are not shifted); `reset`
/// returns the cursor to zero in O(1) and is what a new realization calls
/// after the previous queue has fully drained. Capacity is retained
/// across resets. Not a `VecDeque`: the head-cursor form keeps the hot
/// loop over a flat `Vec` cache-friendly and avoids ring-wrap branches.
#[derive(Debug, Default)]
pub struct ProcessQueue {
    items: Vec<usize>,
    head: usize,
}

impl ProcessQueue {
    /// Reserve capacity for at most `capacity` queued indices.
    pub fn with_capacity(capacity: usize) -> Result<Self, ProcessError> {
        let mut items = Vec::new();
        items.try_reserve_exact(capacity)?;
        Ok(Self { items, head: 0 })
    }

    /// Grow the retained capacity to at least `capacity` (no-op when
    /// sufficient).
    pub fn reserve(&mut self, capacity: usize) -> Result<(), ProcessError> {
        if self.items.capacity() < capacity {
            self.items.clear();
            self.items.try_reserve_exact(capacity)?;
        }
        Ok(())
    }

    /// Enqueue one index (amortized O(1); never allocates once prepared).
    #[inline]
    pub fn push(&mut self, item: usize) {
        self.items.push(item);
    }

    /// Dequeue the oldest index, or `None` when empty.
    #[inline]
    pub fn pop(&mut self) -> Option<usize> {
        if self.head < self.items.len() {
            let item = self.items[self.head];
            self.head += 1;
            Some(item)
        } else {
            None
        }
    }

    /// Number of queued items not yet popped.
    #[inline]
    pub fn pending(&self) -> usize {
        self.items.len() - self.head
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.pending() == 0
    }

    /// O(1) reset for the next realization: drops the logical content and
    /// the head cursor while retaining capacity.
    #[inline]
    pub fn reset(&mut self) {
        self.items.clear();
        self.head = 0;
    }

    /// Retained capacity (probe/diagnostics).
    pub fn capacity(&self) -> usize {
        self.items.capacity()
    }
}

/// Integer neighbor/residual-degree counters over dense vertex indices.
///
/// `reset` fills with zero (O(V), unavoidable for exact counters); the
/// array is reused across realizations with retained capacity. Arithmetic
/// is plain `u32` — an underflow would be a process bug and should panic
/// loudly rather than saturate silently.
#[derive(Debug, Default)]
pub struct NeighborCounters {
    counts: Vec<u32>,
}

impl NeighborCounters {
    /// Prepare `size` zeroed counters.
    pub fn prepare(&mut self, size: usize) -> Result<(), ProcessError> {
        if self.counts.capacity() < size {
            self.counts.clear();
            self.counts.try_reserve_exact(size)?;
        }
        self.counts.clear();
        self.counts.resize(size, 0);
        Ok(())
    }

    #[inline]
    pub fn get(&self, index: usize) -> u32 {
        self.counts[index]
    }

    #[inline]
    pub fn set(&mut self, index: usize, value: u32) {
        self.counts[index] = value;
    }

    #[inline]
    pub fn increment(&mut self, index: usize, delta: u32) {
        self.counts[index] += delta;
    }

    #[inline]
    pub fn decrement(&mut self, index: usize, delta: u32) {
        self.counts[index] -= delta;
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// The counters as a flat slice (analysis and tests).
    #[inline]
    pub fn as_slice(&self) -> &[u32] {
        &self.counts
    }
}

/// Stamped per-vertex flags that skip the O(V) clear between realizations.
///
/// A flag is "marked" iff its slot holds the current stamp; a new
/// generation invalidates every old mark in O(1). The stamp wraps at
/// `u32::MAX` by clearing (correct because a fresh generation never
/// compares against cleared-zero... it clears to 0 and restarts at 1).
#[derive(Debug, Default)]
pub struct StampFlags {
    slots: Vec<u32>,
    stamp: u32,
}

impl StampFlags {
    /// Prepare flags for `size` vertices, all unmarked.
    pub fn prepare(&mut self, size: usize) -> Result<(), ProcessError> {
        if self.slots.capacity() < size {
            self.slots.clear();
            self.slots.try_reserve_exact(size)?;
        }
        self.slots.clear();
        self.slots.resize(size, 0);
        self.stamp = 0;
        Ok(())
    }

    /// Start a fresh generation: every mark becomes invalid in O(1).
    pub fn new_generation(&mut self) {
        if self.stamp == u32::MAX {
            self.slots.iter_mut().for_each(|slot| *slot = 0);
            self.stamp = 0;
        }
        self.stamp += 1;
    }

    #[inline]
    pub fn mark(&mut self, index: usize) {
        self.slots[index] = self.stamp;
    }

    #[inline]
    pub fn is_marked(&self, index: usize) -> bool {
        self.slots[index] == self.stamp
    }
}

/// Capacity failure for the E-track primitives.
#[derive(Debug)]
pub struct ProcessError {
    source: TryReserveError,
}

impl From<TryReserveError> for ProcessError {
    fn from(source: TryReserveError) -> Self {
        Self { source }
    }
}

impl core::fmt::Display for ProcessError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "could not reserve process buffers: {}",
            self.source
        )
    }
}

impl std::error::Error for ProcessError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn queue_is_fifo_and_resets_in_place() {
        let mut queue = ProcessQueue::with_capacity(4).unwrap();
        queue.push(7);
        queue.push(9);
        queue.push(11);
        assert_eq!(queue.pop(), Some(7));
        assert_eq!(queue.pop(), Some(9));
        queue.push(13); // pushed after pops, still behind 11
        assert_eq!(queue.pop(), Some(11));
        assert_eq!(queue.pop(), Some(13));
        assert_eq!(queue.pop(), None);
        assert!(queue.is_empty());
        let capacity = queue.capacity();
        queue.reset();
        assert_eq!(queue.capacity(), capacity, "reset must retain capacity");
        assert!(queue.is_empty());
        queue.push(1);
        assert_eq!(queue.pop(), Some(1));
    }

    #[test]
    fn counters_reset_to_zero() {
        let mut counters = NeighborCounters::default();
        counters.prepare(3).unwrap();
        counters.increment(1, 5);
        assert_eq!(counters.get(1), 5);
        counters.decrement(1, 2);
        assert_eq!(counters.get(1), 3);
        counters.prepare(3).unwrap();
        assert_eq!(counters.as_slice(), &[0, 0, 0]);
    }

    #[test]
    fn stamps_invalidate_across_generations_without_clearing() {
        let mut flags = StampFlags::default();
        flags.prepare(2).unwrap();
        flags.new_generation();
        flags.mark(0);
        assert!(flags.is_marked(0));
        assert!(!flags.is_marked(1));
        flags.new_generation();
        assert!(!flags.is_marked(0), "old marks die with the generation");
        flags.mark(1);
        assert!(flags.is_marked(1));
    }
}
