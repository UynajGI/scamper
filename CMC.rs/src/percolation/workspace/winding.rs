use std::collections::TryReserveError;

use super::components::reserve_len;
use crate::percolation::analyzer::ComponentWinding;

/// Reusable storage for winding-capable component analysis (W2).
///
/// One weighted union-find over active vertices: beside the F3-style
/// connectivity arrays it carries one `i32` potential per periodic axis per
/// vertex (the displacement from a vertex toward its union-find parent) and
/// one winding-basis block per potential root (at most `D` linearly
/// independent generator vectors of `D` components each, maintained in
/// first-nonzero-pivot echelon form during the scan).
///
/// Call [`Self::prepare`] before a repeated workload to make each subsequent
/// analysis allocation-free. Internal scratch is private.
#[derive(Clone, Debug, Default)]
pub struct WindingWorkspace {
    pub(crate) parent: Vec<usize>,
    pub(crate) size: Vec<usize>,
    pub(crate) canonical: Vec<usize>,
    /// Flat `V * dimensions`: displacement from vertex toward its parent.
    pub(crate) potential: Vec<i32>,
    /// Flat `V * dimensions * dimensions`: per-root winding basis blocks.
    pub(crate) basis: Vec<i32>,
    /// Rows used in each root's basis block.
    pub(crate) basis_counts: Vec<u32>,
    /// Root -> component slot map used to list components in ascending
    /// canonical-identity order without sorting.
    pub(crate) root_slot: Vec<u32>,
    /// Flat scratch of `4 * dimensions`: [left potential | right potential |
    /// combined vector | basis-insert work].
    pub(crate) scratch: Vec<i32>,
    /// Union-find chain scratch for potential-aware path compression.
    pub(crate) chain: Vec<usize>,
    pub(crate) components: Vec<ComponentWinding>,
    pub(crate) dimensions: usize,
}

impl WindingWorkspace {
    pub const fn new() -> Self {
        Self {
            parent: Vec::new(),
            size: Vec::new(),
            canonical: Vec::new(),
            potential: Vec::new(),
            basis: Vec::new(),
            basis_counts: Vec::new(),
            root_slot: Vec::new(),
            scratch: Vec::new(),
            chain: Vec::new(),
            components: Vec::new(),
            dimensions: 0,
        }
    }

    /// Reserve and initialize storage for a repeated workload.
    ///
    /// `dimensions` is the embedding dimension `D` of the analyses this
    /// workspace serves; switching dimensions re-reserves every
    /// dimension-dependent array.
    pub fn prepare(
        &mut self,
        vertex_capacity: usize,
        dimensions: usize,
    ) -> Result<(), TryReserveError> {
        reserve_len(&mut self.parent, vertex_capacity)?;
        reserve_len(&mut self.size, vertex_capacity)?;
        reserve_len(&mut self.canonical, vertex_capacity)?;
        reserve_len(&mut self.basis_counts, vertex_capacity)?;
        reserve_len(&mut self.root_slot, vertex_capacity)?;
        reserve_len(&mut self.chain, vertex_capacity)?;
        // Worst case every active vertex is its own component.
        reserve_len(&mut self.components, vertex_capacity)?;
        if dimensions != self.dimensions {
            self.potential.clear();
            self.basis.clear();
            self.scratch.clear();
        }
        reserve_len(
            &mut self.potential,
            vertex_capacity.saturating_mul(dimensions),
        )?;
        reserve_len(
            &mut self.basis,
            vertex_capacity
                .saturating_mul(dimensions)
                .saturating_mul(dimensions),
        )?;
        reserve_len(&mut self.scratch, 4 * dimensions)?;
        self.dimensions = dimensions;
        Ok(())
    }

    /// Reserved vertex capacity shared by all per-vertex arrays.
    pub fn vertex_capacity(&self) -> usize {
        self.parent
            .capacity()
            .min(self.size.capacity())
            .min(self.canonical.capacity())
            .min(self.basis_counts.capacity())
            .min(self.root_slot.capacity())
            .min(self.chain.capacity())
            .min(self.potential.capacity() / self.dimensions.max(1))
            .min(self.basis.capacity() / (self.dimensions.max(1) * self.dimensions.max(1)))
    }

    /// Embedding dimension this workspace is currently reserved for.
    pub const fn dimensions(&self) -> usize {
        self.dimensions
    }

    /// Audit-only retained-capacity accounting for the allocation probe.
    #[cfg(feature = "allocation-probe")]
    #[doc(hidden)]
    pub fn capacity_audit(&self) -> WorkspaceWindingCapacityAudit {
        let vector_bytes = |storage: &Vec<i32>| storage.capacity() * core::mem::size_of::<i32>();
        WorkspaceWindingCapacityAudit {
            parent: self.parent.capacity() * core::mem::size_of::<usize>(),
            size: self.size.capacity() * core::mem::size_of::<usize>(),
            canonical: self.canonical.capacity() * core::mem::size_of::<usize>(),
            potential: vector_bytes(&self.potential),
            basis: vector_bytes(&self.basis),
            basis_counts: self.basis_counts.capacity() * core::mem::size_of::<u32>(),
            root_slot: self.root_slot.capacity() * core::mem::size_of::<u32>(),
            scratch: vector_bytes(&self.scratch),
            chain: self.chain.capacity() * core::mem::size_of::<usize>(),
            components: self.components.capacity() * core::mem::size_of::<ComponentWinding>(),
        }
    }

    pub(crate) fn resize_for(
        &mut self,
        vertex_count: usize,
        dimensions: usize,
    ) -> Result<(), TryReserveError> {
        self.prepare(vertex_count, dimensions)?;
        // Every resize stays within the capacity reserved by `prepare`, so
        // none of these calls can allocate.
        self.parent.resize(vertex_count, usize::MAX);
        self.size.resize(vertex_count, 0);
        self.canonical.resize(vertex_count, usize::MAX);
        self.basis_counts.resize(vertex_count, 0);
        self.root_slot.resize(vertex_count, u32::MAX);
        self.chain.clear();
        self.potential.resize(vertex_count * dimensions, 0);
        self.basis.resize(vertex_count * dimensions * dimensions, 0);
        self.scratch.resize(4 * dimensions, 0);
        self.components.clear();
        Ok(())
    }
}

/// Capacity accounting for the W2 workspace (allocation-probe audits only).
#[cfg(feature = "allocation-probe")]
#[derive(Clone, Copy, Debug, Default)]
pub struct WorkspaceWindingCapacityAudit {
    pub parent: usize,
    pub size: usize,
    pub canonical: usize,
    pub potential: usize,
    pub basis: usize,
    pub basis_counts: usize,
    pub root_slot: usize,
    pub scratch: usize,
    pub chain: usize,
    pub components: usize,
}

#[cfg(feature = "allocation-probe")]
impl WorkspaceWindingCapacityAudit {
    pub fn total_bytes(&self) -> usize {
        self.parent
            + self.size
            + self.canonical
            + self.potential
            + self.basis
            + self.basis_counts
            + self.root_slot
            + self.scratch
            + self.chain
            + self.components
    }
}
