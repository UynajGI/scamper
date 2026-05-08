//! Spin-1/2 Hilbert space implementation.

use super::{HilbertSpace, LocalState, OpType};

/// Spin-1/2 Hilbert space (Ising, Heisenberg, XXZ models).
#[derive(Clone, Copy, Debug, Default)]
pub struct SpinHalfHS;

impl HilbertSpace for SpinHalfHS {
    fn local_dim(&self) -> usize {
        2
    }

    fn is_allowed(&self, states: &[LocalState], op: &OpType) -> bool {
        match op {
            OpType::Identity => true,
            OpType::Diagonal => true,
            OpType::OffDiagonal => states[0] != states[1],
        }
    }

    fn apply(&self, states: &mut [LocalState], op: &OpType) {
        if *op == OpType::OffDiagonal {
            states[0] ^= 1;
            states[1] ^= 1;
        }
    }

    fn diagonal_element(&self, states: &[LocalState], op: &OpType) -> f64 {
        if *op == OpType::Diagonal {
            // For the shifted Heisenberg H' = H - 1/4:
            // Diagonal matrix element: -(SzSz - 1/4) = 1/4 - SzSz
            // ↑↑: 1/4 - 1/4 = 0
            // ↑↓: 1/4 - (-1/4) = 1/2 → but we use 1.0 to match weight J/2
            // ↓↑: 1/4 - (-1/4) = 1/2
            // ↓↓: 1/4 - 1/4 = 0
            // We return 1.0 for anti-aligned so that weight * diag_elem = J/2 * 1.0 = J/2
            // matches the Julia vertex weight of 0.5 (for J=1).
            if states[0] != states[1] {
                1.0
            } else {
                0.0
            }
        } else {
            0.0
        }
    }
}