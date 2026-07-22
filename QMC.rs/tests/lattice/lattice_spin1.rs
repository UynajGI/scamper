//! S>1/2 lattice QMC validation.
//!
//! The S>1/2 bounce fallback is documented as potentially broken in the README.
//! This test verifies that S=1 produces finite results and documents the
//! expected limitations.

use qmc_rs::lattice::ContinuousLatticeEngine;
use qmc_rs::{
    CsrGraph, EdgeCoupling, LatticeConfiguration, QmcKernel, SpinModelBuilder, SpinSpace,
    UpdateSchedule,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// S=1 Heisenberg chain: verify the solver runs and produces finite results.
/// Note: S>1/2 uses a bounce fallback that is known to be approximate.
/// This test does NOT compare against ED — it only verifies the solver
/// doesn't crash and produces physically reasonable output.
/// If S=1 is not supported, this test FAILS (not silently passes).
#[test]
fn s1_heisenberg_chain_produces_finite_results() {
    let n_sites = 3;
    let beta = 2.0;
    let j = 1.0;

    let graph = CsrGraph::chain(n_sites, false).expect("graph"); // open chain (PBC frustrates odd AFM ring)
    let space = SpinSpace::uniform(n_sites, 2).expect("space"); // S=1 → 2S+1=3 states
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::heisenberg(j))
        .build()
        .expect("S=1 model construction should succeed");

    let init: Vec<u16> = vec![1, 1, 1]; // m=0 state for S=1
    let mut configuration = LatticeConfiguration::new(beta, init, &model).expect("config");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(2, 2, 16));
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    let mut energy_sum = 0.0;
    let mut m2_sum = 0.0;
    let mut samples = 0u64;
    for sweep in 0..20_000 {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
        if sweep >= 5_000 {
            let obs =
                qmc_rs::lattice::measure_observables(&configuration, engine.model()).expect("obs");
            energy_sum += obs.energy_total;
            m2_sum += obs.magnetization_z_squared;
            samples += 1;
        }
    }
    let n = samples as f64;
    let energy = energy_sum / n;
    let m2 = m2_sum / n;

    assert!(
        energy.is_finite(),
        "S=1 energy should be finite, got {energy}"
    );
    assert!(
        m2.is_finite() && m2 >= 0.0,
        "S=1 ⟨m²⟩ should be finite and non-negative, got {m2}"
    );

    eprintln!("S=1 Heisenberg: E={energy:.4}, ⟨m²⟩={m2:.4}");
    eprintln!("NOTE: S>1/2 bounce fallback is approximate — results may not match ED");
}
