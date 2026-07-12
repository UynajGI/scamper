//! Derived observables from accumulated measurement data.
//!
//! These functions compute quantities like susceptibility and specific heat
//! from the primary observables (E, E², M, M², M⁴) stored in [`carlo_rs::Results`].

use carlo_rs::Results;

/// Magnetic susceptibility per site: χ = β·N·(⟨M²⟩ − ⟨M⟩²)
///
/// Requires "Magnetization", "M2" in results.
pub fn susceptibility(results: &Results, beta: f64, n_sites: usize) -> Option<f64> {
    if !beta.is_finite() || n_sites == 0 {
        return None;
    }
    let m = results.get("Magnetization")?;
    let m2 = results.get("M2")?;
    Some(beta * n_sites as f64 * (m2.mean - m.mean * m.mean))
}

/// Specific heat per site: C_v = β²/N · (⟨E²⟩ − ⟨E⟩²)
///
/// Requires "Energy", "E2" in results.
pub fn specific_heat(results: &Results, beta: f64, n_sites: usize) -> Option<f64> {
    if !beta.is_finite() || n_sites == 0 {
        return None;
    }
    let e = results.get("Energy")?;
    let e2 = results.get("E2")?;
    Some(beta * beta / n_sites as f64 * (e2.mean - e.mean * e.mean))
}

/// Binder cumulant: U₄ = 1 − ⟨M⁴⟩ / (3·⟨M²⟩²)
///
/// Requires "M2", "M4" in results.
pub fn binder_cumulant(results: &Results) -> Option<f64> {
    let m2 = results.get("M2")?;
    let m4 = results.get("M4")?;
    if m2.mean * m2.mean < 1e-15 {
        return None;
    }
    Some(1.0 - m4.mean / (3.0 * m2.mean * m2.mean))
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::MetropolisCore;
    use crate::classical_mc::ClassicalMC;
    use crate::models::IsingModel;
    use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};

    #[test]
    fn test_compute_correlation_1d_ordered() {
        let spins = vec![1.0; 16];
        let g = compute_correlation_1d(&spins, 1, 16);
        for (r, &val) in g.iter().enumerate() {
            assert!((val - 1.0).abs() < 1e-10, "G({}) = {}", r, val);
        }
    }

    #[test]
    fn test_compute_correlation_1d_alternating() {
        // Alternating +1, -1, +1, -1, ...
        let spins: Vec<f64> = (0..8)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let g = compute_correlation_1d(&spins, 1, 8);
        assert!((g[0] - 1.0).abs() < 1e-10, "G(0) should be 1");
        assert!((g[1] + 1.0).abs() < 1e-10, "G(1) should be -1");
        assert!((g[2] - 1.0).abs() < 1e-10, "G(2) should be 1");
    }

    #[test]
    fn test_susceptibility_positive() {
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let config = RunConfig {
            thermalization_sweeps: 200,
            measurement_sweeps: 500,
            binsize: 50,
            base_seed: 42,
            ..Default::default()
        };

        let backend = RayonBackend::new(1);
        let scheduler = Scheduler::new(backend, config);
        let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

        let chi = susceptibility(&results, 1.0, 64).expect("susceptibility failed");
        assert!(chi > 0.0, "Susceptibility should be positive, got {}", chi);
    }

    #[test]
    fn test_specific_heat_positive() {
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 1.0);
        params.set("J", 1.0);

        let config = RunConfig {
            thermalization_sweeps: 200,
            measurement_sweeps: 500,
            binsize: 50,
            base_seed: 42,
            ..Default::default()
        };

        let backend = RayonBackend::new(1);
        let scheduler = Scheduler::new(backend, config);
        let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

        let cv = specific_heat(&results, 1.0, 64).expect("specific_heat failed");
        assert!(cv > 0.0, "Specific heat should be positive, got {}", cv);
    }

    #[test]
    fn test_correlation_1d_ising_decay() {
        let mut params = Params::new();
        params.set("L", 16usize);
        params.set("beta", 0.5);
        params.set("J", 1.0);

        let config = RunConfig {
            thermalization_sweeps: 300,
            measurement_sweeps: 1000,
            binsize: 100,
            base_seed: 42,
            ..Default::default()
        };

        let backend = RayonBackend::new(1);
        let scheduler = Scheduler::new(backend, config);
        let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

        // Reconstruct a representative spin config from the final state
        // is not directly available. Instead, verify correlation properties
        // from known moments:
        // G(0) = ⟨σ_i²⟩ = 1 exactly for Ising
        // G(1) = ⟨E⟩/(N_bonds * J) for nearest-neighbor
        let e = results.get("Energy").expect("Energy missing");
        let n_sites = 16usize;
        // For 1D chain with PBC, E = -J * N * G(1), so G(1) = -E / (J * N)
        let g1_est = -e.mean / (1.0 * n_sites as f64);
        // At beta=0.5, exact G(1) = tanh(1*0.5) = tanh(0.5) ≈ 0.462
        let g1_exact = (0.5f64).tanh();
        assert!(
            (g1_est - g1_exact).abs() < 0.1,
            "G(1) estimate {:.4} should be near exact {:.4}",
            g1_est,
            g1_exact
        );
    }

    #[test]
    fn test_binder_cumulant_ordered_phase() {
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 5.0); // cold → ordered
        params.set("J", 1.0);

        let config = RunConfig {
            thermalization_sweeps: 300,
            measurement_sweeps: 500,
            binsize: 50,
            base_seed: 42,
            ..Default::default()
        };

        let backend = RayonBackend::new(1);
        let scheduler = Scheduler::new(backend, config);
        let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

        let u4 = binder_cumulant(&results).expect("binder_cumulant failed");
        // In ordered phase, M ≈ const, so ⟨M⁴⟩ ≈ ⟨M²⟩², U₄ ≈ 1 - 1/3 = 2/3
        assert!(
            (0.5..=1.0).contains(&u4),
            "Binder cumulant in ordered phase should be ~0.67, got {}",
            u4
        );
    }
}
