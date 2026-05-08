//! Hilbert space abstractions for quantum systems.

mod spin_half;

pub use spin_half::SpinHalfHS;

/// Local state encoding for lattice sites.
/// Spin-1/2: 0 = Up, 1 = Down
/// Hubbard: 0 = empty, 1 = up, 2 = down, 3 = double
pub type LocalState = u8;

/// Operator type in SSE representation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OpType {
    /// Identity operator (empty vertex)
    Identity,
    /// Diagonal operator (e.g., SzSz, n_i n_j)
    Diagonal,
    /// Off-diagonal operator (e.g., S+S-, hopping)
    OffDiagonal,
}

/// Trait defining Hilbert space rules for operator actions.
pub trait HilbertSpace: Clone {
    /// Local Hilbert space dimension per site.
    fn local_dim(&self) -> usize;

    /// Check if operator is allowed given local states.
    /// states: [source_state, target_state] for bond operators
    fn is_allowed(&self, states: &[LocalState], op: &OpType) -> bool;

    /// Apply operator to local states (in-place modification).
    fn apply(&self, states: &mut [LocalState], op: &OpType);

    /// Compute dimensionless diagonal matrix element.
    /// Returns pure numerical part; engine multiplies by coupling constant.
    fn diagonal_element(&self, states: &[LocalState], op: &OpType) -> f64;
}