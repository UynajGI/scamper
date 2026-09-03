use core::fmt;

/// Dense vertex index.
///
/// The transparent newtype guarantees only a dense numeric representation. It
/// does not carry graph provenance: using an ID with a different graph is a
/// caller contract violation and may panic at an infallible view access.
#[repr(transparent)]
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct VertexId(usize);

impl VertexId {
    #[inline(always)]
    pub(crate) const fn new(index: usize, vertex_count: usize) -> Option<Self> {
        if index < vertex_count {
            Some(Self(index))
        } else {
            None
        }
    }

    #[inline(always)]
    pub(crate) const fn from_valid_index(index: usize) -> Self {
        Self(index)
    }

    /// Returns the dense zero-based index.
    #[inline(always)]
    pub const fn index(self) -> usize {
        self.0
    }
}

impl fmt::Display for VertexId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

/// Dense physical-edge index.
///
/// The transparent newtype guarantees only a dense numeric representation. It
/// does not carry graph provenance: using an ID with a different graph is a
/// caller contract violation and may panic at an infallible view access.
#[repr(transparent)]
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct EdgeId(usize);

impl EdgeId {
    #[inline(always)]
    pub(crate) const fn new(index: usize, edge_count: usize) -> Option<Self> {
        if index < edge_count {
            Some(Self(index))
        } else {
            None
        }
    }

    #[inline(always)]
    pub(crate) const fn from_valid_index(index: usize) -> Self {
        Self(index)
    }

    /// Returns the dense zero-based index.
    #[inline(always)]
    pub const fn index(self) -> usize {
        self.0
    }
}

impl fmt::Display for EdgeId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}
