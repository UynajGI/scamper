//! D2: static directed-graph analysis over the D1 capabilities.
//!
//! Three questions, three distinct result types (never one bag of flags):
//!
//! - [`analyze_weak_components`]: connectivity when arc direction is
//!   ignored — union-find over the physical arc table.
//! - [`analyze_strong_components`]: strongly connected components by an
//!   iterative (non-recursive) Tarjan pass; component identity is the
//!   lowest vertex ID inside the SCC, matching the platform's canonical
//!   identity convention.
//! - [`reachable`]: single-pair source-to-target reachability by
//!   early-exit BFS with stamped visitation; `reachable(v, v)` is true by
//!   definition (the empty path).
//!
//! All three reuse caller-owned workspaces: after `prepare`, repeated
//! analyses of same-size graphs allocate nothing.

use core::fmt;
use std::collections::TryReserveError;

mod process;
mod reach;
mod strong;
mod weak;

pub use process::{DpProcessError, DpSeedPolicy, LayeredDirectedPercolationMC};
pub use reach::{reachable, ReachabilityWorkspace};
pub use strong::{
    analyze_strong_components, StrongComponentResult, StrongComponentSummary,
    StrongComponentWorkspace,
};
pub use weak::{
    analyze_weak_components, WeakComponentResult, WeakComponentSummary, WeakComponentWorkspace,
};

/// Capacity failure while growing a workspace.
#[derive(Debug)]
pub struct DirectedAnalysisError {
    source: TryReserveError,
}

impl From<TryReserveError> for DirectedAnalysisError {
    fn from(source: TryReserveError) -> Self {
        Self { source }
    }
}

impl fmt::Display for DirectedAnalysisError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "could not reserve directed-analysis workspace: {}",
            self.source
        )
    }
}

impl std::error::Error for DirectedAnalysisError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.source)
    }
}

pub(crate) fn find_root(parent: &mut [usize], node: usize) -> usize {
    let mut root = node;
    while parent[root] != root {
        root = parent[root];
    }
    let mut current = node;
    while parent[current] != root {
        let next = parent[current];
        parent[current] = root;
        current = next;
    }
    root
}

/// Union by size with canonical-lowest identity — the F3 convention.
pub(crate) fn union_by_size(
    parent: &mut [usize],
    size: &mut [usize],
    canonical: &mut [usize],
    left: usize,
    right: usize,
) {
    let left = find_root(parent, left);
    let right = find_root(parent, right);
    if left == right {
        return;
    }
    let (big, small) = if size[left] >= size[right] {
        (left, right)
    } else {
        (right, left)
    };
    parent[small] = big;
    size[big] += size[small];
    canonical[big] = canonical[big].min(canonical[small]);
}
