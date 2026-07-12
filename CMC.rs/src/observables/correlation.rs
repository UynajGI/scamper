//! Post-hoc correlation functions.

/// Compute 1D equal-time correlation function G(r) = ⟨s_i · s_{i+r}⟩.
///
/// `spins` is a flat slice with `spin_dim` components per site.
/// For a chain of `n_sites` with PBC, returns G(r) for r = 0..=n_sites/2.
pub fn compute_correlation_1d(spins: &[f64], spin_dim: usize, n_sites: usize) -> Vec<f64> {
    if spin_dim == 0 || n_sites == 0 {
        return Vec::new();
    }
    assert_eq!(
        spins.len(),
        spin_dim * n_sites,
        "spin buffer length does not match spin_dim * n_sites"
    );
    let max_r = n_sites / 2;
    let mut g = vec![0.0; max_r + 1];
    let mut counts = vec![0usize; max_r + 1];

    for i in 0..n_sites {
        for j in 0..n_sites {
            let dist = i.abs_diff(j);
            let r = dist.min(n_sites - dist);
            if r > max_r {
                continue;
            }

            let dot: f64 = (0..spin_dim)
                .map(|k| spins[i * spin_dim + k] * spins[j * spin_dim + k])
                .sum();
            g[r] += dot;
            counts[r] += 1;
        }
    }

    for r in 0..=max_r {
        g[r] /= counts[r] as f64;
    }
    g
}
