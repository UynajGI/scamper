//! Operator string estimators for SSE correlation functions and susceptibilities.
//!
//! Propagates spins through the operator string and computes:
//! - Spatial correlation function C(r)
//! - Structure factor S(q)
//! - Uniform susceptibility chi(0) via energy derivative
//! - Staggered susceptibility chi(pi) via operator string

use super::vertex_data::VertexData;
use super::{SSEEngine, OpType};
use crate::hilbert::HilbertSpace;
use std::f64::consts::PI;

/// Results from operator string correlation measurement.
pub struct CorrelationResult {
    /// Spatial correlation function C(r) = <S^z_i S^z_{i+r}>.
    /// Indexed by minimum periodic distance r = 0, 1, ..., N/2.
    /// For spin-1/2: C(0) = 1/4.
    pub correlation: Vec<f64>,
    /// Structure factor S(q) = sum_{ij} e^{iq(r_i-r_j)} <S^z_i S^z_j>.
    /// Computed for q = 2*pi*k/N, k = 0, ..., N-1.
    /// Note: this is the UNnormalized version (no 1/N factor),
    /// so S(q=0) = <(sum_i Sz_i)^2> and S(q=pi) relates to staggered susceptibility.
    pub structure_factor: Vec<f64>,
    /// Staggered susceptibility chi(pi) via imaginary-time integrated correlations:
    ///   chi_s = (beta/[M(M+1)]) * sum_p M_s(p)^2 / N
    /// where M_s(p) = sum_i (-1)^i Sz_i(p).
    pub staggered_susceptibility: f64,
}

/// S^z value from spin encoding: spin=0 -> +1/2, spin=1 -> -1/2.
#[inline]
fn sz_of(spin: u8) -> f64 {
    if spin == 0 { 0.5 } else { -0.5 }
}

impl<H: HilbertSpace> SSEEngine<H> {
    /// Measure correlation functions by propagating spins through the operator string.
    ///
    /// For each position p in the operator string, the spin configuration |alpha(p)>
    /// is tracked. Correlations are averaged over all M positions.
    ///
    /// Structure factor: S(q) = sum_{ij} e^{iq(r_i-r_j)} <Sz_i Sz_j>
    /// (unnormalized, so S(0) = <M^2> total fluctuation).
    ///
    /// Staggered susceptibility follows Sandvik 1991:
    ///   chi_s = (beta/[M(M+1)]) * sum_p M_s(p)^2 / N
    ///
    /// Note: uniform susceptibility chi(0) is NOT measured here because the SSE
    /// algorithm conserves total Sz. In the Sz=0 sector, chi(0) would always be 0.
    /// It can be computed from the energy derivative dE/dH instead.
    ///
    /// Returns None if there are no operators.
    pub fn measure_correlations(&self) -> Option<CorrelationResult> {
        if self.op_seq.n_operators == 0 {
            return None;
        }

        let n = self.lattice.n_sites;
        let m = self.op_seq.max_length;
        let nf = n as f64;
        let mf = m as f64;
        let half_n = n / 2;

        // Displacement-specific correlation sums:
        // corr_sum_d[d] = sum_{p, i} Sz_i(p) * Sz_{i+d}(p)
        // where d = (j - i) mod N, d = 0, ..., N-1.
        // This preserves the phase information needed for the Fourier transform.
        let mut corr_sum_d = vec![0.0_f64; n];

        // Staggered susceptibility accumulators
        let mut ms2_total = 0.0_f64; // sum_p M_s(p)^2

        // Self-consistent spin state after full sweep
        let mut spins = self.spins.clone();

        for p in 0..m {
            // S^z values for current configuration
            let sz: Vec<f64> = spins.iter().map(|&s| sz_of(s)).collect();

            // Staggered magnetization
            let mut msp = 0.0_f64;
            for i in 0..n {
                let stag = if i & 1 == 0 { 1.0 } else { -1.0 };
                msp += stag * sz[i];
            }

            // Displacement-specific correlation: for each site pair (i,j),
            // accumulate into bin d = (j - i) mod N.
            // This preserves the direction of displacement for correct FT.
            for i in 0..n {
                let si = sz[i];
                for d in 0..n {
                    let j = (i + d) % n;
                    corr_sum_d[d] += si * sz[j];
                }
            }

            // Staggered susceptibility: sum_p M_s(p)^2
            ms2_total += msp * msp;

            // Propagate through operator at position p
            let v = &self.op_seq.vertices[p];
            if v.vertex_idx != 0 && VertexData::op_type(v.vertex_idx) == OpType::OffDiagonal {
                let (si, sj, _) = self.bond_list[v.bond_idx];
                spins[si] ^= 1;
                spins[sj] ^= 1;
            }
        }

        // Spatial correlation function: C(r) = <Sz_i Sz_{i+r}>
        // C(r) = corr_sum_d[r] / (M * N), averaged over all N starting sites.
        let correlation: Vec<f64> = (0..=half_n)
            .map(|r| corr_sum_d[r] / (mf * nf))
            .collect();

        // Structure factor: S(q_k) = sum_d e^{i*q_k*d} * [corr_sum_d[d] / M]
        // = corr_sum_d[0]/M + sum_{d=1}^{N-1} cos(q_k*d) * corr_sum_d[d]/M
        // (imaginary parts cancel since C(-d) = C(d))
        let structure_factor: Vec<f64> = (0..n)
            .map(|k| {
                let q = 2.0 * PI * k as f64 / nf;
                let mut sq = corr_sum_d[0] / mf;
                for d in 1..n {
                    sq += (q * d as f64).cos() * corr_sum_d[d] / mf;
                }
                sq
            })
            .collect();

        // Staggered susceptibility: chi_s = (beta/[M(M+1)]) * sum_p M_s(p)^2 / N
        let chi_staggered = self.beta * ms2_total / (mf * (mf + 1.0) * nf);

        Some(CorrelationResult {
            correlation,
            structure_factor,
            staggered_susceptibility: chi_staggered,
        })
    }

    /// Compute specific heat from operator count fluctuations.
    ///
    /// C = (<n^2> - <n>^2 - <n>) / N
    pub fn compute_specific_heat(&self, n_mean: f64, n_sq_mean: f64) -> f64 {
        let ns = self.lattice.n_sites as f64;
        (n_sq_mean - n_mean * n_mean - n_mean) / ns
    }
}
