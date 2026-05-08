//! Worm-based loop update for SSE algorithm.

use super::{SSEEngine, vertex_data::VertexData, vertex_list::VertexList};
use crate::hilbert::HilbertSpace;
use rand::Rng;
use rand::RngExt;

impl<H: HilbertSpace> SSEEngine<H> {
    /// Perform worm update on operator sequence.
    ///
    /// 1. Build vertex list (worldline topology)
    /// 2. Execute multiple worm traversals
    /// 3. Reconstruct spin state from vertex list
    ///
    /// For the shifted Heisenberg where W_diag = W_offdiag = J/2:
    /// The directed loop equations with zero bounce give deterministic
    /// conversion between diagonal and off-diagonal operators.
    pub fn loopupdate<R: Rng>(&mut self, rng: &mut R) {
        if self.op_seq.n_operators == 0 {
            return;
        }

        // Build worldline topology
        let mut vertex_list = VertexList::new(self.lattice.n_sites, self.op_seq.max_length);
        vertex_list.build(&self.op_seq, &self.bond_list);

        // Number of worms: proportional to operator count
        let num_worms = std::cmp::max(1, self.op_seq.n_operators / 4);

        for _ in 0..num_worms {
            self.worm_traverse(&mut vertex_list, rng);
        }

        // Reconstruct spin state
        self.reconstruct_state(&vertex_list, rng);
    }

    /// Execute a single worm traversal.
    ///
    /// For Heisenberg (Δ = 1): deterministic switching.
    /// For XXZ (|Δ| ≤ 1): deterministic switching (same as Heisenberg).
    /// For XXZ (Δ > 1): bounce with probability (Δ-1)/(Δ+1).
    fn worm_traverse<R: Rng>(&mut self, vertex_list: &mut VertexList, rng: &mut R) {
        // Find random non-Identity entry point
        let (p0, l0) = loop {
            let p = rng.random_range(0..self.op_seq.max_length);
            let l = rng.random_range(0..4);
            let (_next_leg, next_pos) = vertex_list.link(l, p);
            if next_pos != usize::MAX {
                break (p, l);
            }
        };

        let mut leg_in = l0;
        let mut p = p0;
        let delta = self.loop_param;

        // Safety limit to prevent infinite loops
        let max_steps = self.op_seq.max_length * 4;
        let mut steps = 0;

        loop {
            let vertex = &mut self.op_seq.vertices[p];
            if vertex.vertex_idx == 0 {
                break; // Identity — shouldn't happen with correct entry selection
            }

            // XXZ scatter: (leg_out, new_vertex_idx, bounce_prob)
            let (leg_out, new_vertex_idx, bounce_prob) =
                VertexData::xxz_scatter(leg_in, vertex.vertex_idx, delta);

            if bounce_prob > 0.0 && rng.random::<f64>() < bounce_prob {
                // Bounce: exit on same site leg (leg ^ 2)
                let leg_bounce = leg_in ^ 2;
                vertex.vertex_idx = new_vertex_idx;

                // Check closure after bounce
                if p == p0 && leg_bounce == l0 {
                    break;
                }

                let (next_leg, next_p) = vertex_list.link(leg_bounce, p);
                if next_p == usize::MAX {
                    break;
                }
                if next_p == p0 && next_leg == l0 {
                    break;
                }
                leg_in = next_leg;
                p = next_p;
            } else {
                // Switch: change vertex type, exit on opposite leg
                vertex.vertex_idx = new_vertex_idx;

                // Check closure: returned to entry point from same leg
                if p == p0 && leg_out == l0 {
                    break;
                }

                // Follow worldline to next vertex
                let (next_leg, next_p) = vertex_list.link(leg_out, p);
                if next_p == usize::MAX {
                    break;
                }

                // Check if we've returned to entry from the other direction
                if next_p == p0 && next_leg == l0 {
                    break;
                }

                leg_in = next_leg;
                p = next_p;
            }

            steps += 1;
            if steps > max_steps {
                break; // Safety limit
            }
        }
    }

    /// Reconstruct spin state after loop update.
    ///
    /// After the worm modifies vertices, we read the spin configuration from the
    /// modified vertices. For each site, we find the first vertex on its worldline
    /// and read the input leg spin. This gives the spin on the worldline at the
    /// "bottom" of the imaginary time axis.
    ///
    /// Since the worm conserves total Sz (each scatter flips one spin ↑→↓ and one ↓→↑
    /// for anti-aligned spins), the reconstructed state should have the same total Sz
    /// as the pre-loop state.
    fn reconstruct_state<R: Rng>(&mut self, vertex_list: &VertexList, _rng: &mut R) {
        for site in 0..self.lattice.n_sites {
            let (leg, pos) = vertex_list.v_first(site);
            if pos == usize::MAX {
                // No operators on this worldline — keep previous spin
                continue;
            }
            let vertex = &self.op_seq.vertices[pos];
            let legs = VertexData::leg_states(vertex.vertex_idx);
            self.spins[site] = legs[leg];
        }
    }
}
