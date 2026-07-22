//! QMC verification: Rabi model ⟨n⟩ at T=0 (β=1024) for FSS near r=1.
//!
//! Compares QMC ⟨n⟩ to exact ED values to verify QMC correctness.
//! Output: results/qmc_verify/verify.csv

use qmc_rs::impurity::spin_boson::occupation::{
    model::{CavityMode, OccupationSpinBosonModel},
    OccupationWorldlineSampler,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::fs;
use std::time::Instant;

const DELTA: f64 = 1.0;
const OMEGAS: [f64; 4] = [10.0, 50.0, 200.0, 500.0];
const BETA: f64 = 1024.0;
const CUTOFF: usize = 80;
const SLICES: usize = 8;

const THERMAL_SWEEPS: u64 = 30_000;
const MEASURE_SWEEPS: u64 = 80_000;

fn run_mc(omega: f64, g: f64, seed: u64) -> (f64, f64, f64, f64, f64) {
    let model = OccupationSpinBosonModel::rabi(
        DELTA,
        vec![CavityMode::new(omega, g, CUTOFF).unwrap()],
    )
    .unwrap();

    let mut sampler = OccupationWorldlineSampler::new(model, BETA, SLICES, 0).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    let mut n_sum = 0.0_f64;
    let mut n2_sum = 0.0_f64;
    let mut x2_sum = 0.0_f64;
    let mut x4_sum = 0.0_f64;
    let mut u4_sum = 0.0_f64;
    let mut samples = 0_u64;

    for sweep in 0..(THERMAL_SWEEPS + MEASURE_SWEEPS) {
        sampler.sweep(&mut rng).unwrap();
        if sweep >= THERMAL_SWEEPS {
            let obs = sampler.measure().unwrap();
            let x2 = obs.quadrature_variance;
            let x4 = obs.quadrature_fourth_moment;
            n_sum += obs.total_boson_number;
            n2_sum += obs.total_boson_number * obs.total_boson_number;
            x2_sum += x2;
            x4_sum += x4;
            // U4 from each sample: 1 - x4/(3*x2²)
            if x2.abs() > 1e-30 {
                u4_sum += 1.0 - x4 / (3.0 * x2 * x2);
            }
            samples += 1;
        }
    }

    let n_f = samples as f64;
    let n_mean = n_sum / n_f;
    let n_err = ((n2_sum / n_f - n_mean * n_mean) / n_f).sqrt();
    let x2_mean = x2_sum / n_f;
    let x4_mean = x4_sum / n_f;
    let u4_sample = u4_sum / n_f;
    // Ensemble U4: 1 - ⟨x⁴⟩/(3⟨x²⟩²)
    let u4_ens = if x2_mean.abs() > 1e-30 { 1.0 - x4_mean / (3.0 * x2_mean * x2_mean) } else { 0.0 };
    (n_mean, n_err, x2_mean, x4_mean, u4_ens)
}

// Simple ED for ground-state ⟨n⟩
fn ed_n(omega: f64, g: f64) -> f64 {
    let (_, n_exp) = ed_ground_state(omega, g);
    n_exp
}

fn ed_x2(omega: f64, g: f64) -> f64 {
    let (v, _) = ed_ground_state(omega, g);
    let dim = 2 * CUTOFF;
    let cutoff = CUTOFF;
    let mut x2_exp = 0.0_f64;
    for state in 0..dim {
        let n = (state >> 1) as f64;
        // diagonal: 2n + 1
        x2_exp += (2.0 * n + 1.0) * v[state] * v[state];
        // n → n+2
        let np2 = (state >> 1) + 2;
        if np2 < cutoff {
            let target = (np2 << 1) | (state & 1);
            let amp = ((n + 1.0) * (n + 2.0)).sqrt();
            x2_exp += 2.0 * amp * v[state] * v[target];
        }
    }
    x2_exp
}

fn ed_ground_state(omega: f64, g: f64) -> (Vec<f64>, f64) {
    let dim = 2 * CUTOFF;
    let mut h = vec![vec![0.0_f64; dim]; dim];
    for n in 0..CUTOFF {
        h[2*n][2*n] = omega * n as f64 - 0.5 * DELTA;
        h[2*n+1][2*n+1] = omega * n as f64 + 0.5 * DELTA;
        if n + 1 < CUTOFF {
            let amp = g * (n as f64 + 1.0).sqrt();
            h[2*n+1][2*(n+1)] += amp; h[2*(n+1)][2*n+1] += amp;
            h[2*n][2*(n+1)+1] += amp; h[2*(n+1)+1][2*n] += amp;
        }
    }

    // Power method for lowest eigenvalue
    let shift = h.iter().map(|row| row.iter().map(|x| x.abs()).fold(0.0, f64::max)).fold(0.0, f64::max) + 1.0;
    let mut b = vec![vec![0.0_f64; dim]; dim];
    for i in 0..dim { for j in 0..dim { b[i][j] = -h[i][j]; } }
    for i in 0..dim { b[i][i] += shift; }
    let mut v = vec![0.0_f64; dim];
    v[0] = 1.0;
    for _ in 0..1000 {
        let mut bv = vec![0.0_f64; dim];
        for i in 0..dim { for j in 0..dim { bv[i] += b[i][j] * v[j]; } }
        let norm = bv.iter().map(|x| x*x).sum::<f64>().sqrt();
        if norm < 1e-30 { break; }
        for i in 0..dim { v[i] = bv[i] / norm; }
    }

    let mut n_exp = 0.0_f64;
    for n in 0..CUTOFF {
        n_exp += n as f64 * (v[2*n]*v[2*n] + v[2*n+1]*v[2*n+1]);
    }
    (v, n_exp)
}

fn ed_x2x4(omega: f64, g: f64) -> (f64, f64, f64) {
    let (v, _) = ed_ground_state(omega, g);
    let dim = 2 * CUTOFF;
    let cutoff = CUTOFF;

    // Build x2 matrix
    let mut x2 = vec![vec![0.0_f64; dim]; dim];
    for state in 0..dim {
        let n = (state >> 1) as f64;
        x2[state][state] = 2.0 * n + 1.0;
        let np2 = (state >> 1) + 2;
        if np2 < cutoff {
            let target = (np2 << 1) | (state & 1);
            let amp = ((n + 1.0) * (n + 2.0)).sqrt();
            x2[state][target] = amp;
            x2[target][state] = amp;
        }
    }
    // x4 = x2 × x2
    let mut x4 = vec![vec![0.0_f64; dim]; dim];
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim { x4[i][j] += x2[i][k] * x2[k][j]; }
        }
    }

    // ⟨x²⟩
    let mut x2_exp = 0.0_f64;
    for i in 0..dim {
        for j in 0..dim { x2_exp += v[i] * x2[i][j] * v[j]; }
    }
    // ⟨x⁴⟩
    let mut x4_exp = 0.0_f64;
    for i in 0..dim {
        for j in 0..dim { x4_exp += v[i] * x4[i][j] * v[j]; }
    }
    let u4 = if x2_exp > 1e-30 { 1.0 - x4_exp / (3.0 * x2_exp * x2_exp) } else { 0.0 };
    (x2_exp, x4_exp, u4)
}

fn main() {
    let out_dir = "results/qmc_verify";
    fs::create_dir_all(out_dir).unwrap();

    // Uniform grid: Δλ=0.02 in [0.2, 2.0] → 91 points
    let lam_min = 0.20; let lam_max = 2.0; let n_pts = 91;
    let lambdas: Vec<f64> = (0..n_pts).map(|i| lam_min + (lam_max - lam_min) * i as f64 / (n_pts - 1) as f64).collect();

    eprintln!("=== QMC Verify: ⟨n⟩ + U4(x) (β=1024, Δλ=0.02, {n_pts} pts) ===");
    eprintln!("ηs={OMEGAS:?}");

    let start = Instant::now();

    // CSV 1: verify.csv (⟨n⟩ data)
    let mut csv_n = String::from("lambda");
    for &omega in &OMEGAS {
        let eta = omega / DELTA;
        csv_n.push_str(&format!(",n_qmc_eta{eta:.0},n_err_eta{eta:.0},n_ed_eta{eta:.0},ntilde_qmc_eta{eta:.0},ntilde_ed_eta{eta:.0}"));
    }
    csv_n.push('\n');

    // CSV 2: verify_u4.csv (U4 and x² data)
    let mut csv_u4 = String::from("lambda");
    for &omega in &OMEGAS {
        let eta = omega / DELTA;
        csv_u4.push_str(&format!(",u4_qmc_eta{eta:.0},u4_ed_eta{eta:.0},x2_qmc_eta{eta:.0},x2_ed_eta{eta:.0}"));
    }
    csv_u4.push('\n');

    for (li, &lambda) in lambdas.iter().enumerate() {
        csv_n.push_str(&format!("{lambda:.6}"));
        csv_u4.push_str(&format!("{lambda:.6}"));
        let t0 = Instant::now();

        for (oi, &omega) in OMEGAS.iter().enumerate() {
            let eta = omega / DELTA;
            let g = lambda * (omega * DELTA).sqrt();
            let seed = (li as u64) * 1000 + (oi as u64) * 100 + BETA as u64;

            let (n_mean, n_err, x2_qmc, _x4, u4_qmc) = run_mc(omega, g, seed);
            let (x2_ed, _x4_ed, u4_ed) = ed_x2x4(omega, g);
            let n_ed = ed_n(omega, g);
            let ntilde_qmc = eta * n_mean;
            let ntilde_ed = eta * n_ed;

            csv_n.push_str(&format!(",{n_mean:.10e},{n_err:.10e},{n_ed:.10e},{ntilde_qmc:.10e},{ntilde_ed:.10e}"));
            csv_u4.push_str(&format!(",{u4_qmc:.10e},{u4_ed:.10e},{x2_qmc:.10e},{x2_ed:.10e}"));
        }
        csv_n.push('\n');
        csv_u4.push('\n');
        let dt = t0.elapsed().as_secs_f64();
        if li % 10 == 0 { eprintln!("  λ={lambda:.3} ({li}/{n_pts}) [{dt:.0}s]"); }
    }

    fs::write(format!("{out_dir}/verify.csv"), &csv_n).unwrap();
    fs::write(format!("{out_dir}/verify_u4.csv"), &csv_u4).unwrap();
    eprintln!("\nSaved verify.csv + verify_u4.csv ({:.0}s)", start.elapsed().as_secs_f64());
}
