//! Worldline topology for SSE loop update.
//!
//! Tracks worldline connections between vertices as doubly-linked chains
//! through imaginary time. Each vertex has 4 legs:
//!   leg 0: site_i input  (spin entering from below on site_i)
//!   leg 1: site_j input  (spin entering from below on site_j)
//!   leg 2: site_i output (spin exiting upward on site_i)
//!   leg 3: site_j output (spin exiting upward on site_j)
//!
//! Following `links[leg][pos]` traces a worldline through imaginary time.

use super::engine::OperatorSequence;
use crate::lattice::BondType;

/// Tracks worldline connections between vertices.
pub struct VertexList {
    /// Flat array: links[leg * max_length + pos] = (next_leg, next_pos)
    /// next_pos = usize::MAX means no connection.
    links: Vec<(usize, usize)>,
    /// First (leg, pos) per site on its worldline.
    v_first: Vec<(usize, usize)>,
    /// Last (leg, pos) per site on its worldline.
    v_last: Vec<(usize, usize)>,
    /// Maximum sequence length (stride for leg indexing).
    max_length: usize,
}

impl VertexList {
    pub fn new(n_sites: usize, max_length: usize) -> Self {
        let sentinel = (usize::MAX, usize::MAX);
        Self {
            links: vec![sentinel; 4 * max_length],
            v_first: vec![sentinel; n_sites],
            v_last: vec![sentinel; n_sites],
            max_length,
        }
    }

    /// Access link at (leg, pos).
    #[inline]
    pub fn link(&self, leg: usize, pos: usize) -> (usize, usize) {
        self.links[leg * self.max_length + pos]
    }

    /// Get first vertex on a site's worldline.
    #[inline]
    pub fn v_first(&self, site: usize) -> (usize, usize) {
        self.v_first[site]
    }

    /// Build worldline connections from operator sequence.
    ///
    /// Mirrors Julia vertex_list.jl:make_vertex_list!
    /// For each non-Identity vertex, connects input legs to previous output
    /// legs on the same site, forming periodic worldline chains.
    pub fn build(
        &mut self,
        op_seq: &OperatorSequence,
        bond_list: &[(usize, usize, BondType)],
    ) {
        let sentinel = (usize::MAX, usize::MAX);
        for slot in &mut self.links {
            *slot = sentinel;
        }
        for slot in &mut self.v_first {
            *slot = sentinel;
        }
        for slot in &mut self.v_last {
            *slot = sentinel;
        }

        for p in 0..op_seq.max_length {
            let v = &op_seq.vertices[p];
            if v.vertex_idx == 0 {
                continue; // Identity — no legs
            }
            let (site_i, site_j, _) = bond_list[v.bond_idx];

            // site_i uses legs 0 (input) and 2 (output)
            // site_j uses legs 1 (input) and 3 (output)
            let sites = [site_i, site_j];
            let in_legs = [0, 1];
            let out_legs = [2, 3];

            for s in 0..2 {
                let site = sites[s];
                let in_leg = in_legs[s];
                let out_leg = out_legs[s];

                let prev = self.v_last[site];
                if prev.1 != usize::MAX {
                    // Link previous output to current input
                    let (prev_out_leg, prev_pos) = prev;
                    self.links[prev_out_leg * self.max_length + prev_pos] = (in_leg, p);
                    self.links[in_leg * self.max_length + p] = (prev_out_leg, prev_pos);
                } else {
                    self.v_first[site] = (in_leg, p);
                }
                self.v_last[site] = (out_leg, p);
            }
        }

        // Close periodic chains: link last back to first
        for site in 0..self.v_first.len() {
            let first = self.v_first[site];
            let last = self.v_last[site];
            if first.1 != usize::MAX && last.1 != usize::MAX {
                let (f_leg, f_pos) = first;
                let (l_leg, l_pos) = last;
                self.links[l_leg * self.max_length + l_pos] = (f_leg, f_pos);
                self.links[f_leg * self.max_length + f_pos] = (l_leg, l_pos);
            }
        }
    }
}
