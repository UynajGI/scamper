//! Lanczos algorithm for finding the ground state energy of a sparse Hamiltonian.

use super::SparseHamiltonian;
use rand_xoshiro::Xoshiro256PlusPlus;
use rand_xoshiro::rand_core::SeedableRng;
use rand::{RngExt};

/// Find ground state energy using Lanczos iteration.
///
/// Returns the converged ground state energy estimate.
pub fn lanczos_ground_state(ham: &SparseHamiltonian, tol: f64, max_iter: usize) -> f64 {
    let dim = ham.dim();
    let max_iter = max_iter.min(dim);

    // Random initial vector
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(12345);
    let mut v = vec![0.0f64; dim];
    let mut norm = 0.0;
    for v_i in &mut v {
        *v_i = rng.random::<f64>() - 0.5;
        norm += *v_i * *v_i;
    }
    norm = norm.sqrt();
    for v_i in &mut v {
        *v_i /= norm;
    }

    // Lanczos recurrence: build tridiagonal matrix
    let mut alphas = Vec::with_capacity(max_iter);
    let mut betas = Vec::with_capacity(max_iter);
    let mut w = vec![0.0f64; dim];
    let mut v_prev = vec![0.0f64; dim];

    for iter in 0..max_iter {
        // w = H · v
        ham.mat_vec(&v, &mut w);

        // alpha = v† · w
        let alpha: f64 = v.iter().zip(w.iter()).map(|(a, b)| a * b).sum();
        alphas.push(alpha);

        // w = w - alpha · v - beta · v_prev
        for i in 0..dim {
            w[i] -= alpha * v[i];
            if iter > 0 {
                w[i] -= betas[iter - 1] * v_prev[i];
            }
        }

        // beta = ||w||
        let beta: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        betas.push(beta);

        // Check convergence
        if beta < tol || iter >= max_iter - 1 {
            break;
        }

        // v_prev = v, v = w / beta
        v_prev.clone_from(&v);
        for i in 0..dim {
            v[i] = w[i] / beta;
        }
    }

    // Build dense tridiagonal matrix and find minimum eigenvalue
    // using dense symmetric eigenvalue solver (Jacobi-like)
    let n = alphas.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return alphas[0];
    }

    // Extract off-diagonal elements
    let offdiag: Vec<f64> = if betas.len() >= n - 1 {
        betas[..n - 1].to_vec()
    } else {
        // Pad with zeros if Lanczos terminated early
        let mut od = betas.clone();
        od.resize(n - 1, 0.0);
        od
    };

    tridiagonal_min_eigenvalue(&alphas, &offdiag)
}

/// Find minimum eigenvalue of symmetric tridiagonal matrix using Sturm bisection.
///
/// Uses the Sturm sequence property: the number of sign changes in
/// p_0, p_1, ..., p_n equals the number of eigenvalues less than λ.
/// Combined with bisection, this gives the minimum eigenvalue to
/// arbitrary precision. Guaranteed correct, no rotation bugs.
fn tridiagonal_min_eigenvalue(diag: &[f64], offdiag: &[f64]) -> f64 {
    let n = diag.len();

    // Gershgorin bounds for bisection interval
    let mut lower = f64::INFINITY;
    let mut upper = f64::NEG_INFINITY;
    for i in 0..n {
        let radius = (if i > 0 { offdiag[i - 1].abs() } else { 0.0 })
            + (if i + 1 < n { offdiag[i].abs() } else { 0.0 });
        lower = lower.min(diag[i] - radius);
        upper = upper.max(diag[i] + radius);
    }

    // Sturm count using LDL^T factorization of T - xI:
    // q[0] = d[0] - x
    // q[k] = d[k] - x - e[k-1]^2 / q[k-1]
    // Number of eigenvalues < x = number of negative q[k]
    let sturm_count = |x: f64| -> usize {
        let mut count = 0;
        let safe_q = |q: f64| -> f64 {
            if q.abs() < 1e-30 {
                if q >= 0.0 { 1e-30 } else { -1e-30 }
            } else {
                q
            }
        };
        let mut q = diag[0] - x;
        if q < 0.0 { count += 1; }
        for k in 1..n {
            q = diag[k] - x - offdiag[k - 1] * offdiag[k - 1] / safe_q(q);
            if q < 0.0 { count += 1; }
        }
        count
    };

    // Bisection: find the smallest eigenvalue
    // We want λ_min such that sturm_count(λ_min) = 1 and sturm_count(λ_min - ε) = 0
    let mut lo = lower;
    let mut hi = upper;

    for _ in 0..200 {
        let mid = (lo + hi) / 2.0;
        if sturm_count(mid) == 0 {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-12 {
            break;
        }
    }

    (lo + hi) / 2.0
}
