//! Exact diagonalization study of finite-temperature Rabi QPT.
//!
//! Computes multiple observables via exact thermal expectation:
//!   ñ = η·⟨n⟩     — photon number (insensitive to QPT, for reference)
//!   ΔE₁₂          — energy gap (closes at λ_c)
//!   ⟨σz⟩          — spin polarization
//!   ⟨x²⟩, ⟨x⁴⟩   — position moments → U4(x) Binder cumulant
//!   C_V           — specific heat (β² Var(H))
//!
//! H = (Δ/2)σz + Ωa†a + gσx(a+a†), with x ≡ a+a†.

use nalgebra::{DMatrix, SymmetricEigen};
use std::fs;
use std::io::Write;
use std::time::Instant;

const DELTA: f64 = 1.0;
const ETAS: [f64; 3] = [50.0, 200.0, 1000.0];
const BETAS: [f64; 8] = [0.5, 1.0, 2.0, 4.0, 16.0, 64.0, 256.0, 1024.0];
const CUTOFF: usize = 80;

fn build_rabi_h(omega: f64, g: f64, cutoff: usize) -> DMatrix<f64> {
    let dim = 2 * cutoff;
    let mut h = DMatrix::zeros(dim, dim);
    for n in 0..cutoff {
        h[(2 * n, 2 * n)] = omega * n as f64 - 0.5 * DELTA;
        h[(2 * n + 1, 2 * n + 1)] = omega * n as f64 + 0.5 * DELTA;
        if n + 1 < cutoff {
            let amp = g * (n as f64 + 1.0).sqrt();
            h[(2 * n + 1, 2 * (n + 1))] += amp;
            h[(2 * (n + 1), 2 * n + 1)] += amp;
            h[(2 * n, 2 * (n + 1) + 1)] += amp;
            h[(2 * (n + 1) + 1, 2 * n)] += amp;
        }
    }
    h
}

/// Build the position operator x = a + a† in the occupation basis.
/// x is block-diagonal in spin, connects n ↔ n±1.
fn build_x_op(cutoff: usize) -> DMatrix<f64> {
    let dim = 2 * cutoff;
    let mut x = DMatrix::zeros(dim, dim);
    for n in 0..cutoff {
        if n + 1 < cutoff {
            let amp = (n as f64 + 1.0).sqrt();
            for s in 0..2 {
                x[(2 * n + s, 2 * (n + 1) + s)] += amp;
                x[(2 * (n + 1) + s, 2 * n + s)] += amp;
            }
        }
    }
    x
}

struct Obs {
    n_mean: f64,
    sigma_z: f64,
    x2: f64,
    x4: f64,
    energy: f64,
    energy2: f64,
}

/// Compute all thermal observables from the eigensystem.
fn thermal_obs(
    eigvals_sorted: &[f64],
    eigvecs: &DMatrix<f64>,
    x2_eig: &DMatrix<f64>,
    x4_eig: &DMatrix<f64>,
    sz_diag: &[(usize, f64)], // (index, σz eigenvalue) for diagonal σz
    beta: f64,
    cutoff: usize,
) -> Obs {
    let ground = eigvals_sorted[0];
    let mut z = 0.0_f64;
    let mut n_sum = 0.0_f64;
    let mut sz_sum = 0.0_f64;
    let mut x2_sum = 0.0_f64;
    let mut x4_sum = 0.0_f64;
    let mut e_sum = 0.0_f64;
    let mut e2_sum = 0.0_f64;

    for k in 0..eigvals_sorted.len() {
        let exponent = -beta * (eigvals_sorted[k] - ground);
        if exponent < -700.0 {
            continue;
        }
        let weight = exponent.exp();
        z += weight;
        let ek = eigvals_sorted[k];

        // ⟨k|n̂|k⟩
        let mut n_k = 0.0_f64;
        for n in 0..cutoff {
            let vd = eigvecs[(2 * n, k)];
            let vu = eigvecs[(2 * n + 1, k)];
            n_k += n as f64 * (vd * vd + vu * vu);
        }
        n_sum += weight * n_k;

        // ⟨k|σz|k⟩ = Σ_i σz_i |V_{ik}|²
        let mut sz_k = 0.0_f64;
        for &(i, szi) in sz_diag {
            sz_k += szi * eigvecs[(i, k)] * eigvecs[(i, k)];
        }
        sz_sum += weight * sz_k;

        // ⟨k|x²|k⟩, ⟨k|x⁴|k⟩ from pre-transformed operators
        x2_sum += weight * x2_eig[(k, k)];
        x4_sum += weight * x4_eig[(k, k)];

        e_sum += weight * ek;
        e2_sum += weight * ek * ek;
    }

    let inv_z = 1.0 / z;
    let energy = e_sum * inv_z;
    let energy2 = e2_sum * inv_z;
    Obs {
        n_mean: n_sum * inv_z,
        sigma_z: sz_sum * inv_z,
        x2: x2_sum * inv_z,
        x4: x4_sum * inv_z,
        energy,
        energy2,
    }
}

fn main() {
    let out_dir = "results/ed_fss";
    fs::create_dir_all(out_dir).unwrap();

    // Dense λ scan — cover both near-λ_c and deep broken phase
    let lambda_min = 0.50;
    let lambda_max = 50.0;
    let n_pts = 491; // Δλ = 0.1
    let lambdas: Vec<f64> = (0..n_pts)
        .map(|i| lambda_min + (lambda_max - lambda_min) * i as f64 / (n_pts - 1) as f64)
        .collect();

    // Precompute static operators
    let x_op = build_x_op(CUTOFF);
    let x2_op = &x_op * &x_op;
    let x4_op = &x2_op * &x2_op;
    // σz diagonal: index 2n → -1, 2n+1 → +1
    let sz_diag: Vec<(usize, f64)> = (0..CUTOFF)
        .flat_map(|n| [(2 * n, -1.0_f64), (2 * n + 1, 1.0_f64)])
        .collect();

    eprintln!("=== Rabi QPT ED: multi-observable ===");
    eprintln!("Δ={DELTA}, ηs={ETAS:?}, βs={BETAS:?}");
    eprintln!("cutoff={CUTOFF}, λ∈[{lambda_min},{lambda_max}], {n_pts} pts\n");

    let start = Instant::now();
    let mut log = fs::File::create(format!("{out_dir}/run.log")).unwrap();

    for &beta in &BETAS {
        let t0 = Instant::now();
        eprint!("β={beta} ... ");

        let mut csv = String::from("lambda");
        for &eta in &ETAS {
            csv.push_str(&format!(
                ",ntilde_eta{eta:.0},sigmaz_eta{eta:.0},x2_eta{eta:.0},u4x_eta{eta:.0},gap_eta{eta:.0},cv_eta{eta:.0}"
            ));
        }
        csv.push('\n');

        for &lambda in &lambdas {
            csv.push_str(&format!("{lambda:.6}"));

            for &eta in &ETAS {
                let omega = eta * DELTA;
                let g = lambda * (omega * DELTA).sqrt();
                let h = build_rabi_h(omega, g, CUTOFF);
                let eig = SymmetricEigen::new(h);

                // Sort eigenvalues (and track sort order for eigenvectors)
                let n = eig.eigenvalues.len();
                let mut order: Vec<usize> = (0..n).collect();
                order.sort_by(|&a, &b| {
                    eig.eigenvalues[a]
                        .total_cmp(&eig.eigenvalues[b])
                });

                let sorted_eigvals: Vec<f64> = order.iter().map(|&i| eig.eigenvalues[i]).collect();

                // Permute eigenvector columns
                let mut perm_eigvecs = DMatrix::zeros(n, n);
                for (new_col, &old_col) in order.iter().enumerate() {
                    for row in 0..n {
                        perm_eigvecs[(row, new_col)] = eig.eigenvectors[(row, old_col)];
                    }
                }

                // Transform x², x⁴ to eigenbasis (diagonal elements only needed)
                let x2_eig = perm_eigvecs.transpose() * &x2_op * &perm_eigvecs;
                let x4_eig = perm_eigvecs.transpose() * &x4_op * &perm_eigvecs;

                let obs = thermal_obs(
                    &sorted_eigvals,
                    &perm_eigvecs,
                    &x2_eig,
                    &x4_eig,
                    &sz_diag,
                    beta,
                    CUTOFF,
                );

                let ntilde = eta * obs.n_mean;
                let u4x = if obs.x2.abs() > 1e-30 {
                    1.0 - obs.x4 / (3.0 * obs.x2 * obs.x2)
                } else {
                    0.0
                };
                let gap = sorted_eigvals[1] - sorted_eigvals[0];
                let cv = beta * beta * (obs.energy2 - obs.energy * obs.energy);

                csv.push_str(&format!(
                    ",{ntilde:.10e},{:.10e},{:.10e},{:.10e},{gap:.10e},{cv:.10e}",
                    obs.sigma_z, obs.x2, u4x
                ));
            }
            csv.push('\n');
        }

        let path = format!("{out_dir}/beta_{beta}.csv");
        fs::write(&path, &csv).unwrap();
        let dt = t0.elapsed().as_secs_f64();
        writeln!(log, "β={beta}: {path} ({dt:.0}s)").unwrap();
        eprintln!("{dt:.0}s");
    }

    eprintln!("\nTotal: {:.0}s", start.elapsed().as_secs_f64());
}
