//! Derived observables from accumulated measurement data.
//!
//! These functions compute quantities like susceptibility and specific heat
//! from the primary observables (E, E², M, M², M⁴) stored in [`carlo_rs::Results`].

use carlo_rs::Results;

/// Magnetic susceptibility per site: χ = β·N·(⟨M²⟩ − ⟨M⟩²)
///
/// Requires "Magnetization", "M2" in results.
pub fn susceptibility(results: &Results, beta: f64, n_sites: usize) -> Option<f64> {
    let m = results.get("Magnetization")?;
    let m2 = results.get("M2")?;
    Some(beta * n_sites as f64 * (m2.mean - m.mean * m.mean))
}

/// Specific heat per site: C_v = β²/N · (⟨E²⟩ − ⟨E⟩²)
///
/// Requires "Energy", "E2" in results.
pub fn specific_heat(results: &Results, beta: f64, n_sites: usize) -> Option<f64> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::MetropolisCore;
    use crate::classical_mc::ClassicalMC;
    use crate::models::IsingModel;
    use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};

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
