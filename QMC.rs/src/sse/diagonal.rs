//! Diagonal update for SSE algorithm.

use super::{SSEEngine, Vertex};
use super::vertex_data::VertexData;
use crate::hilbert::{HilbertSpace, OpType};
use rand::RngExt;

impl<H: HilbertSpace> SSEEngine<H> {
    /// Perform diagonal update on operator sequence.
    ///
    /// For each position:
    /// - Identity: try to insert diagonal (or off-diagonal for pure off-diagonal models)
    /// - Diagonal: try to remove it
    /// - OffDiagonal: try to remove it (for pure off-diagonal models), propagate spin if kept
    pub fn diagonal_update<R: RngExt>(&mut self, rng: &mut R) {
        let mut current_spins = self.spins.clone();
        let n_bonds = self.bond_list.len();
        let delta = self.loop_param;
        let is_xxz = (delta - 1.0).abs() > 1e-6;
        let is_pure_offdiag = is_xxz && delta.abs() < 1e-6; // XY model (Δ=0)

        for p in 0..self.op_seq.max_length {
            let vertex = &self.op_seq.vertices[p];

            match VertexData::op_type(vertex.vertex_idx) {
                OpType::Identity => {
                    if self.op_seq.n_operators >= self.op_seq.max_length {
                        continue;
                    }

                    let bond_idx = rng.random_range(0..n_bonds);
                    let (site_i, site_j, bond_type) = self.bond_list[bond_idx];
                    let states = [current_spins[site_i], current_spins[site_j]];
                    let n = self.op_seq.n_operators;
                    let m = self.op_seq.max_length;

                    if is_pure_offdiag {
                        // XY model (Δ=0): insert off-diagonal operators only
                        let w_offdiag = self.weights[&bond_type];
                        if states[0] != states[1] {
                            let p_insert = self.beta * w_offdiag * n_bonds as f64 / (m - n) as f64;
                            if rng.random::<f64>() < p_insert {
                                let vertex_idx = VertexData::offdiag_vertex(current_spins[site_i], current_spins[site_j]);
                                self.op_seq.vertices[p] = Vertex { bond_idx, op: OpType::OffDiagonal, vertex_idx };
                                self.op_seq.n_operators += 1;
                            }
                        }
                    } else {
                        // Heisenberg or XXZ (Δ≠0): insert diagonal operators
                        let w = if is_xxz {
                            // XXZ: diagonal only on anti-aligned bonds
                            if states[0] == states[1] {
                                0.0
                            } else {
                                self.weights[&bond_type]
                            }
                        } else {
                            // Heisenberg: use HilbertSpace method
                            let weight = self.weights[&bond_type];
                            weight * self.hs.diagonal_element(&states, &OpType::Diagonal)
                        };

                        if w > 1e-10 {
                            let p_insert = self.beta * w * n_bonds as f64 / (m - n) as f64;
                            if rng.random::<f64>() < p_insert {
                                let vertex_idx = VertexData::diag_vertex(current_spins[site_i], current_spins[site_j]);
                                self.op_seq.vertices[p] = Vertex { bond_idx, op: OpType::Diagonal, vertex_idx };
                                self.op_seq.n_operators += 1;
                            }
                        }
                    }
                }
                OpType::Diagonal => {
                    // Try to remove diagonal operator
                    let (site_i, site_j, bond_type) = self.bond_list[vertex.bond_idx];
                    let states = [current_spins[site_i], current_spins[site_j]];

                    let w = if is_xxz {
                        if states[0] == states[1] {
                            0.0
                        } else {
                            self.weights[&bond_type]
                        }
                    } else {
                        let weight = self.weights[&bond_type];
                        weight * self.hs.diagonal_element(&states, &OpType::Diagonal)
                    };
                    let n = self.op_seq.n_operators;
                    let m = self.op_seq.max_length;
                    let p_remove = (m - n + 1) as f64 / (self.beta * w * n_bonds as f64);

                    if rng.random::<f64>() < p_remove.min(1.0) {
                        self.op_seq.vertices[p] = Vertex::default();
                        self.op_seq.n_operators -= 1;
                    }
                }
                OpType::OffDiagonal => {
                    // For pure off-diagonal models (XY, Δ=0): try to remove off-diagonal operator
                    if is_pure_offdiag {
                        let (_, _, bond_type) = self.bond_list[vertex.bond_idx];
                        let w_offdiag = self.weights[&bond_type];
                        let n = self.op_seq.n_operators;
                        let m = self.op_seq.max_length;
                        let p_remove = (m - n + 1) as f64 / (self.beta * w_offdiag * n_bonds as f64);
                        if rng.random::<f64>() < p_remove.min(1.0) {
                            self.op_seq.vertices[p] = Vertex::default();
                            self.op_seq.n_operators -= 1;
                            continue;
                        }
                    }
                    // Propagate state through off-diagonal operator
                    let (site_i, site_j, _) = self.bond_list[vertex.bond_idx];
                    let mut states = [current_spins[site_i], current_spins[site_j]];
                    self.hs.apply(&mut states, &OpType::OffDiagonal);
                    current_spins[site_i] = states[0];
                    current_spins[site_j] = states[1];
                }
            }
        }

        self.spins = current_spins;
    }
}
