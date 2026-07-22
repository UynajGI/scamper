//! Ergodicity tests for continuous-time lattice QMC.
//!
//! Verifies that independent runs from widely separated initial states
//! (ferromagnetic, Néel, random) converge to the same thermal expectation
//! values, confirming that the update schedule is ergodic.

use qmc_rs::lattice::ContinuousLatticeEngine;
use qmc_rs::{
    CsrGraph, EdgeCoupling, LatticeConfiguration, QmcKernel, SpinModelBuilder, SpinSpace,
    UpdateSchedule,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// 4-site PBC Heisenberg chain at β=1.0, J=1.0.
/// Three distinct initial states must converge to the same ⟨E⟩ and ⟨m²⟩.
#[test]
fn lattice_ergodicity_multi_init_convergence() {
    let n_sites = 4;
    let beta = 1.0;
    let j = 1.0;
    let n_thermalization = 10_000;
    let n_measurement = 30_000;

    let inits: Vec<(&str, Vec<u16>)> = vec![
        ("ferromagnetic", vec![0, 0, 0, 0]),
        ("neel", vec![0, 1, 0, 1]),
        ("random", vec![1, 0, 1, 0]),
    ];

    let mut results: Vec<(f64, f64)> = Vec::new();

    for (label, init_state) in &inits {
        let graph = CsrGraph::chain(n_sites, true).expect("graph");
        let space = SpinSpace::uniform(n_sites, 1).expect("space");
        let model = SpinModelBuilder::new(graph, space)
            .uniform_edge(EdgeCoupling::heisenberg(j))
            .build()
            .expect("model");
        let mut configuration =
            LatticeConfiguration::new(beta, init_state.clone(), &model).expect("config");
        let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(2, 2, 16));
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

        let mut energy_sum = 0.0;
        let mut m2_sum = 0.0;
        let mut samples = 0u64;

        for sweep in 0..(n_thermalization + n_measurement) {
            engine.sweep(&mut configuration, &mut rng).expect("sweep");
            if sweep >= n_thermalization {
                let obs = qmc_rs::lattice::measure_observables(&configuration, engine.model())
                    .expect("obs");
                energy_sum += obs.energy_total;
                m2_sum += obs.magnetization_z_squared;
                samples += 1;
            }
        }

        assert!(samples > 0, "no samples collected for init '{label}'");
        let n = samples as f64;
        let avg_e = energy_sum / n;
        let avg_m2 = m2_sum / n;
        results.push((avg_e, avg_m2));
        eprintln!("  init={label:14}  ⟨E⟩={avg_e:.6}  ⟨m²⟩={avg_m2:.6}  (n={samples})");
    }

    // All pairs must agree within tolerance.
    let tol = 0.05;
    for i in 0..results.len() {
        for j_idx in (i + 1)..results.len() {
            let de = (results[i].0 - results[j_idx].0).abs();
            let dm2 = (results[i].1 - results[j_idx].1).abs();
            assert!(
                de < tol,
                "⟨E⟩ mismatch between init {} and {}: |{:.6} - {:.6}| = {:.6} ≥ {tol}",
                inits[i].0,
                inits[j_idx].0,
                results[i].0,
                results[j_idx].0,
                de
            );
            assert!(
                dm2 < tol,
                "⟨m²⟩ mismatch between init {} and {}: |{:.6} - {:.6}| = {:.6} ≥ {tol}",
                inits[i].0,
                inits[j_idx].0,
                results[i].1,
                results[j_idx].1,
                dm2
            );
        }
    }
}
