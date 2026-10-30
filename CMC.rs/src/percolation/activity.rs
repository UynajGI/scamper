use crate::{EdgeId, VertexId};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MaskState {
    All,
    None,
    Dense,
}

/// Replaceable activity storage with constant-time all-active/all-inactive states.
///
/// The byte layout is private. Dense storage is retained when switching to an
/// endpoint state so later sampling can reuse its allocation.
#[derive(Clone, Debug)]
pub(crate) struct ActivityMask {
    len: usize,
    state: MaskState,
    dense: Vec<u8>,
}

#[allow(dead_code)]
impl ActivityMask {
    /// An all-active mask without allocating dense storage.
    pub(crate) const fn all(len: usize) -> Self {
        Self {
            len,
            state: MaskState::All,
            dense: Vec::new(),
        }
    }

    /// An all-inactive mask without allocating dense storage.
    pub(crate) const fn none(len: usize) -> Self {
        Self {
            len,
            state: MaskState::None,
            dense: Vec::new(),
        }
    }

    /// Build a dense mask from byte flags. Nonzero bytes are normalized to one.
    pub(crate) fn dense(flags: &[u8]) -> Self {
        let dense = flags.iter().map(|&flag| u8::from(flag != 0)).collect();
        Self {
            len: flags.len(),
            state: MaskState::Dense,
            dense,
        }
    }

    #[inline]
    pub(crate) const fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub(crate) const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns activity at a dense index, or `None` when out of bounds.
    #[inline]
    pub(crate) fn get(&self, index: usize) -> Option<bool> {
        (index < self.len).then(|| self.get_valid(index))
    }

    /// Set one dense index. Returns `false` without mutation when out of bounds.
    pub(crate) fn set(&mut self, index: usize, active: bool) -> bool {
        if index >= self.len {
            return false;
        }
        self.make_dense();
        self.dense[index] = u8::from(active);
        true
    }

    /// Fill the logical mask. Endpoint states retain any dense allocation.
    #[inline]
    pub(crate) fn fill(&mut self, active: bool) {
        self.state = if active {
            MaskState::All
        } else {
            MaskState::None
        };
    }

    /// Change the logical length and fill value, retaining reusable capacity.
    pub(crate) fn reconfigure(&mut self, len: usize, active: bool) {
        self.len = len;
        self.fill(active);
        if self.dense.len() > len {
            self.dense.truncate(len);
        }
    }

    /// Bytes reserved by the retained dense backend.
    #[inline]
    pub(crate) fn owned_bytes(&self) -> usize {
        self.dense.capacity()
    }

    #[inline]
    pub(crate) fn get_valid(&self, index: usize) -> bool {
        debug_assert!(index < self.len);
        match self.state {
            MaskState::All => true,
            MaskState::None => false,
            MaskState::Dense => self.dense[index] != 0,
        }
    }

    pub(crate) fn dense_mut(&mut self) -> &mut [u8] {
        self.make_dense();
        &mut self.dense
    }

    fn make_dense(&mut self) {
        if self.dense.len() < self.len {
            self.dense.resize(self.len, 0);
        } else if self.dense.len() > self.len {
            self.dense.truncate(self.len);
        }
        match self.state {
            MaskState::All => self.dense.fill(1),
            MaskState::None => self.dense.fill(0),
            MaskState::Dense => {}
        }
        self.state = MaskState::Dense;
    }
}

/// Minimal read-only vertex-activity capability.
///
/// IDs must belong to the unchanged graph whose vertex count matches this
/// activity source. Violating that caller contract may panic.
pub trait VertexActivity {
    /// Number of vertices in this activity domain.
    fn vertex_count(&self) -> usize;

    fn vertex_active(&self, vertex: VertexId) -> bool;
}

/// Minimal read-only physical-edge-activity capability.
///
/// IDs must belong to the unchanged graph whose edge count matches this
/// activity source. Violating that caller contract may panic.
pub trait EdgeActivity {
    /// Number of physical edges in this activity domain.
    fn edge_count(&self) -> usize;

    fn edge_active(&self, edge: EdgeId) -> bool;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{build_chain, GraphView, StaticConfiguration, UndirectedGraphView};

    #[test]
    fn mask_semantics_and_bounds() {
        let mut mask = ActivityMask::all(3);
        assert_eq!(mask.get(0), Some(true));
        assert_eq!(mask.get(3), None);
        assert!(!mask.set(3, false));
        assert!(mask.set(1, false));
        assert_eq!([mask.get(0), mask.get(1)], [Some(true), Some(false)]);
        assert!(mask.owned_bytes() >= 3);

        mask.fill(false);
        assert_eq!(mask.get(0), Some(false));
        let capacity = mask.owned_bytes();
        mask.fill(true);
        assert_eq!(mask.get(2), Some(true));
        assert_eq!(mask.owned_bytes(), capacity);

        let dense = ActivityMask::dense(&[0, 7, 0]);
        assert_eq!(dense.get(1), Some(true));
    }

    #[test]
    fn reconfigure_changes_length_and_reuses_storage() {
        let mut mask = ActivityMask::none(4);
        assert!(mask.set(0, true));
        let capacity = mask.owned_bytes();
        mask.reconfigure(2, true);
        assert_eq!(mask.len(), 2);
        assert_eq!(mask.get(1), Some(true));
        assert_eq!(mask.get(2), None);
        mask.reconfigure(4, false);
        assert_eq!(mask.owned_bytes(), capacity);
        assert_eq!(mask.get(3), Some(false));
    }

    #[test]
    #[should_panic]
    fn vertex_capability_contract_panics_for_wrong_id_domain() {
        let large = build_chain(3, false);
        let small = build_chain(2, false);
        let foreign = large.vertex_id(2).expect("large-graph vertex");
        let configuration = StaticConfiguration::new(small.vertex_count(), small.edge_count());
        let _ = configuration.vertex_active(foreign);
    }

    #[test]
    #[should_panic]
    fn edge_capability_contract_panics_for_wrong_id_domain() {
        let large = build_chain(3, false);
        let small = build_chain(2, false);
        let foreign = large.edge_id(1).expect("large-graph edge");
        let configuration = StaticConfiguration::new(small.vertex_count(), small.edge_count());
        let _ = configuration.edge_active(foreign);
    }
}
