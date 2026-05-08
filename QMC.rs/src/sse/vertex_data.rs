//! Vertex data for spin-1/2 Heisenberg SSE.
//!
//! Vertex index encoding:
//! - 0 = Identity
//! - 1 = Diagonal ↑↑→↑↑
//! - 2 = Diagonal ↑↓→↑↓
//! - 3 = Diagonal ↓↑→↓↑
//! - 4 = Diagonal ↓↓→↓↓
//! - 5 = OffDiagonal ↑↓→↓↑
//! - 6 = OffDiagonal ↓↑→↑↓

use crate::hilbert::OpType;

/// Vertex data for spin-1/2 Heisenberg.
/// All weights are J/4 (diagonal) and J/2 (off-diagonal).
/// For the shifted Hamiltonian, the directed loop solution with
/// zero bounce is deterministic: Diagonal ↔ OffDiagonal.
pub struct VertexData;

impl VertexData {
    /// Get leg states for a vertex.
    /// Returns [spin_i_in, spin_j_in, spin_i_out, spin_j_out]
    #[inline]
    pub const fn leg_states(vertex_idx: u8) -> [u8; 4] {
        match vertex_idx {
            0 => [0, 0, 0, 0], // Identity (unused legs)
            1 => [0, 0, 0, 0], // Diagonal ↑↑→↑↑
            2 => [0, 1, 0, 1], // Diagonal ↑↓→↑↓
            3 => [1, 0, 1, 0], // Diagonal ↓↑→↓↑
            4 => [1, 1, 1, 1], // Diagonal ↓↓→↓↓
            5 => [0, 1, 1, 0], // OffDiagonal ↑↓→↓↑
            6 => [1, 0, 0, 1], // OffDiagonal ↓↑→↑↓
            _ => [0, 0, 0, 0],
        }
    }

    /// Deterministic scatter for Heisenberg with zero bounce.
    /// (leg_out, new_vertex_idx)
    ///
    /// For spin-1/2 Heisenberg with the shifted Hamiltonian where
    /// diagonal and off-diagonal vertices both have equal weight,
    /// the directed loop equations with zero bounce give deterministic
    /// conversion: Diagonal ↔ OffDiagonal.
    ///
    /// Exit leg pairs: (0↔1, 2↔3) matching Julia's xor(leg-1,1)+1 (1-indexed).
    ///   leg 0 (site_i input) ↔ leg 1 (site_j input)
    ///   leg 2 (site_i output) ↔ leg 3 (site_j output)
    ///
    /// The worm flips the spin on BOTH the entering and exiting legs' sites.
    /// This converts between diagonal (same input/output) and off-diagonal
    /// (different input/output) configurations.
    #[inline]
    pub const fn scatter(leg_in: usize, vertex_idx: u8) -> (usize, u8) {
        let leg_out = leg_in ^ 1; // 0↔1, 2↔3
        let legs = Self::leg_states(vertex_idx);
        // Flip spins on both entering and exiting legs
        let mut new_legs = legs;
        new_legs[leg_in] ^= 1;
        new_legs[leg_out] ^= 1;
        // Find matching vertex
        let new_idx = if new_legs[0] == new_legs[2] && new_legs[1] == new_legs[3] {
            // Diagonal: input == output
            Self::diag_vertex(new_legs[0], new_legs[1])
        } else {
            // OffDiagonal: input != output
            Self::offdiag_vertex(new_legs[0], new_legs[1])
        };
        (leg_out, new_idx)
    }

    #[inline]
    pub const fn op_type(vertex_idx: u8) -> OpType {
        match vertex_idx {
            0 => OpType::Identity,
            1..=4 => OpType::Diagonal,
            5 | 6 => OpType::OffDiagonal,
            _ => OpType::Identity,
        }
    }

    #[inline]
    pub const fn diag_vertex(spin_i: u8, spin_j: u8) -> u8 {
        match (spin_i, spin_j) {
            (0, 0) => 1,
            (0, 1) => 2,
            (1, 0) => 3,
            (1, 1) => 4,
            _ => 1,
        }
    }

    #[inline]
    pub const fn offdiag_vertex(spin_i: u8, spin_j: u8) -> u8 {
        // Off-diagonal only exists for anti-aligned spins
        match (spin_i, spin_j) {
            (0, 1) => 5, // ↑↓ → ↓↑
            (1, 0) => 6, // ↓↑ → ↑↓
            _ => 5,      // fallback (shouldn't happen for valid calls)
        }
    }

    /// Weight for a given vertex_idx and coupling J.
    /// Diagonal vertices: J/4, OffDiagonal: J/2.
    #[inline]
    pub const fn weight(vertex_idx: u8, j: f64) -> f64 {
        match vertex_idx {
            1..=4 => j * 0.25,
            5 | 6 => j * 0.5,
            _ => 0.0,
        }
    }

    /// XXZ diagonal matrix element.
    /// H' = H + J*Δ/4 per bond. SSE expands in -H' + C.
    ///   Aligned:   0
    ///   Anti-aligned: J*Δ/2
    /// At Δ=1 (Heisenberg): 0 and J/2 ✓
    /// At Δ=0 (XY): both 0 (pure off-diagonal)
    #[inline]
    pub fn xxz_diag_element(spin_i: u8, spin_j: u8, delta: f64, j: f64) -> f64 {
        if spin_i == spin_j {
            0.0
        } else {
            j * delta / 2.0
        }
    }

    /// XXZ scatter for directed loop algorithm.
    ///
    /// For |Δ| ≤ 1: deterministic switching (same as Heisenberg).
    ///   Diagonal → OffDiagonal and vice versa, no bounce.
    ///
    /// For Δ > 1: bounce with probability (Δ-1)/(Δ+1) for anti-aligned spins.
    ///   The worm can bounce (exit same leg, keep vertex) or switch
    ///   (exit opposite leg, flip vertex type).
    ///
    /// Returns (leg_out, new_vertex_idx, bounce_probability).
    /// When bounce_probability = 0, the scatter is deterministic.
    #[inline]
    pub fn xxz_scatter(leg_in: usize, vertex_idx: u8, delta: f64) -> (usize, u8, f64) {
        let legs = Self::leg_states(vertex_idx);
        let spin_i = legs[0];
        let spin_j = legs[1];
        let aligned = spin_i == spin_j;

        if aligned {
            // Aligned spins: only diagonal vertex exists, no off-diagonal.
            // Bounce back (exit same site leg, leg_in ^ 2).
            let leg_out = leg_in ^ 2;
            (leg_out, vertex_idx, 0.0)
        } else if delta <= 1.0 {
            // Anti-aligned, |Δ| ≤ 1: deterministic switching (same as Heisenberg).
            // Flip spins on BOTH entering and exiting legs (matching `scatter` logic).
            let leg_out = leg_in ^ 1; // 0↔1, 2↔3
            let mut new_legs = legs;
            new_legs[leg_in] ^= 1;
            new_legs[leg_out] ^= 1;

            // For Δ=0 (XY model): no diagonal vertices exist, stay off-diagonal.
            if delta.abs() < 1e-6 {
                let new_idx = Self::offdiag_vertex(new_legs[0], new_legs[1]);
                (leg_out, new_idx, 0.0)
            } else {
                let new_idx = if new_legs[0] == new_legs[2] && new_legs[1] == new_legs[3] {
                    // Input == output → diagonal vertex
                    Self::diag_vertex(new_legs[0], new_legs[1])
                } else {
                    // Input != output → off-diagonal vertex
                    Self::offdiag_vertex(new_legs[0], new_legs[1])
                };
                (leg_out, new_idx, 0.0)
            }
        } else {
            // Anti-aligned, Δ > 1: bounce or switch.
            // Bounce probability: (Δ - 1)/(Δ + 1)
            // On bounce: exit same leg, keep vertex type.
            // On switch: exit opposite leg, flip vertex type.
            let bounce_prob = (delta - 1.0) / (delta + 1.0);
            let leg_out = leg_in ^ 1; // opposite leg for switch
            let mut new_legs = legs;
            new_legs[leg_in] ^= 1;
            new_legs[leg_out] ^= 1;
            let new_idx = if new_legs[0] == new_legs[2] && new_legs[1] == new_legs[3] {
                Self::diag_vertex(new_legs[0], new_legs[1])
            } else {
                Self::offdiag_vertex(new_legs[0], new_legs[1])
            };
            (leg_out, new_idx, bounce_prob)
        }
    }
}
