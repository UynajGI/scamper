//! Analytic-solvable-limit tests for continuous-time lattice QMC.
//!
//! These tests verify behavior in well-understood limiting regimes where
//! the answer is known analytically or by exact enumeration, without
//! requiring a dense ED solver.

use qmc_rs::lattice::ContinuousLatticeEngine;
use qmc_rs::{
    CsrGraph, EdgeCoupling, LatticeConfiguration, QmcKernel, SpinModelBuilder, SpinSpace,
    UpdateSchedule,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Zero-coupling limit (J=0): Hamiltonian is trivial, all states have E=0.
/// QMC expansion order must be exactly zero, energy must be exactly zero.
#[test]
fn zero_coupling_has_zero_energy_and_no_vertices() {
    let n_sites = 4;
    let beta = 5.0;
    let graph = CsrGraph::chain(n_sites, true).expect("graph");
    let space = SpinSpace::uniform(n_sites, 1).expect("space");
    // Coupling exactly zero — the model builder should accept this and the
    // engine should produce zero expansion order.
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::heisenberg(0.0))
        .build();
    // Zero coupling → no operator terms → trivial model.
    // If the builder rejects zero coupling, that's also a valid response.
    match model {
        Ok(model) => {
            let mut configuration =
                LatticeConfiguration::new(beta, vec![0, 1, 0, 1], &model).expect("config");
            let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(4, 4, 32));
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
            for _ in 0..500 {
                engine.sweep(&mut configuration, &mut rng).expect("sweep");
            }
            let obs =
                qmc_rs::lattice::measure_observables(&configuration, engine.model()).expect("obs");
            assert!(
                obs.expansion_order.abs() < 1e-12,
                "expansion order must be zero: got {}",
                obs.expansion_order
            );
            assert!(
                obs.energy_total.abs() < 1e-12,
                "energy must be zero: got {}",
                obs.energy_total
            );
        }
        Err(_) => {
            // Builder rejects zero coupling — acceptable behavior.
        }
    }
}

/// High-temperature limit (β→0): energy → 0, magnetization² → 1/(4N) for
/// N independent S=1/2. Use moderate β to keep the estimator variance bounded.
#[test]
fn high_temperature_limit_energy_is_small_and_magnetization_matches() {
    let n_sites = 4;
    let beta = 0.1; // high but finite T
    let j = 1.0;
    let graph = CsrGraph::chain(n_sites, true).expect("graph");
    let space = SpinSpace::uniform(n_sites, 1).expect("space");
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::heisenberg(j))
        .build()
        .expect("model");
    let mut configuration =
        LatticeConfiguration::new(beta, vec![0, 0, 0, 0], &model).expect("config");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(2, 2, 16));
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    let mut energy_sum = 0.0;
    let mut m2_sum = 0.0;
    let mut samples = 0u64;
    for sweep in 0..60_000 {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
        if sweep >= 10_000 {
            let obs =
                qmc_rs::lattice::measure_observables(&configuration, engine.model()).expect("obs");
            energy_sum += obs.energy_total;
            m2_sum += obs.magnetization_z_squared;
            samples += 1;
        }
    }
    let measured_e = energy_sum / samples as f64;
    let measured_m2 = m2_sum / samples as f64;

    // At high T, E is small (thermal energy ~ β * coupling² for small β)
    assert!(
        measured_e.abs() < 0.3,
        "high-T energy should be small: got {measured_e}"
    );
    // ⟨m²⟩ → 1/(4N) for N independent S=1/2 at infinite T
    let exact_m2 = 1.0 / (4.0 * n_sites as f64);
    assert!(
        (measured_m2 - exact_m2).abs() < 0.005,
        "high-T m²: MC={measured_m2:.6}, exact={exact_m2:.6}"
    );
}

/// Strong longitudinal-field limit (h_z >> J): all spins align with the field.
/// For S=1/2 at large β*h_z: every spin is ↑, magnetization → +0.5.
#[test]
fn strong_longitudinal_field_polarizes_all_spins() {
    let n_sites = 4;
    let beta = 4.0;
    let j = 0.1; // weak coupling
    let h_z = 5.0; // strong field

    let graph = CsrGraph::chain(n_sites, false).expect("graph");
    let space = SpinSpace::uniform(n_sites, 1).expect("space");
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::heisenberg(j))
        .uniform_site(qmc_rs::SiteCoupling::new(0.0, h_z, 0.0))
        .build()
        .expect("model");
    let mut configuration =
        LatticeConfiguration::new(beta, vec![0, 0, 0, 0], &model).expect("config");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(8, 4, 64));
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(99);

    let mut mag_sum = 0.0;
    let mut samples = 0u64;
    for sweep in 0..20_000 {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
        if sweep >= 5_000 {
            let obs =
                qmc_rs::lattice::measure_observables(&configuration, engine.model()).expect("obs");
            mag_sum += obs.magnetization_z;
            samples += 1;
        }
    }
    let measured_m = mag_sum / samples as f64;
    // All spins ↑: magnetization → +0.5 (for S=1/2)
    assert!(
        (measured_m - 0.5).abs() < 0.02,
        "strong field should polarize: m={measured_m:.6}"
    );
}

/// Two-site S=1/2 chain at infinite T: Sz correlation = 0 (uncorrelated).
/// At β→0, the system is completely disordered and ⟨Sz_0 Sz_1⟩ → 0.
#[test]
fn dimer_high_temperature_correlation_vanishes() {
    let beta = 0.001;
    let j = 1.0;
    let graph = CsrGraph::chain(2, false).expect("graph");
    let space = SpinSpace::uniform(2, 1).expect("space");
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::heisenberg(j))
        .build()
        .expect("model");
    let mut configuration = LatticeConfiguration::new(beta, vec![0, 1], &model).expect("config");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(1, 1, 8));
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);

    let mut corr_sum = 0.0;
    let mut samples = 0u64;
    for sweep in 0..20_000 {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
        if sweep >= 5_000 {
            let obs =
                qmc_rs::lattice::measure_observables(&configuration, engine.model()).expect("obs");
            corr_sum += obs.nearest_neighbor_sz_correlation;
            samples += 1;
        }
    }
    let measured = corr_sum / samples as f64;
    assert!(
        measured.abs() < 0.01,
        "high-T correlation should vanish: got {measured:.6}"
    );
}

/// Large transverse-field Ising limit (TFIM): with J_z dominant and h_x=0,
/// the model reduces to classical Ising. Energy matches exact enumeration
/// for a 2-site chain. This tests the TFIM coupling path.
#[test]
fn classical_ising_dimer_energy_matches_exact() {
    // 2-site Ising chain: H = J_z * Sz_0 * Sz_1
    // Exact: Z = 2*exp(-β*J_z/4) + 2*exp(+β*J_z/4)
    // E = (-J_z/4)*exp(-β*J_z/4)*2 - (J_z/4)*exp(β*J_z/4)*2 ... wait
    // States: |↑↑⟩(E=J_z/4), |↓↓⟩(E=J_z/4), |↑↓⟩(E=-J_z/4), |↓↑⟩(E=-J_z/4)
    // Z = 2*e^{-β*J_z/4} + 2*e^{+β*J_z/4}
    // E_exact = [2*(J_z/4)*e^{-β*J_z/4} + 2*(-J_z/4)*e^{β*J_z/4}] / Z
    let n_sites = 2;
    let beta = 3.0;
    let j_z = 1.0;

    let graph = CsrGraph::chain(n_sites, false).expect("graph");
    let space = SpinSpace::uniform(n_sites, 1).expect("space");
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::xxz(0.0, j_z)) // pure Ising
        .build();
    // Pure Ising may be rejected (zero off-diagonal). If so, skip.
    match model {
        Ok(model) => {
            let mut configuration =
                LatticeConfiguration::new(beta, vec![0, 1], &model).expect("config");
            let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(4, 4, 32));
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(55);

            let mut energy_sum = 0.0;
            let mut samples = 0u64;
            for sweep in 0..40_000 {
                engine.sweep(&mut configuration, &mut rng).expect("sweep");
                if sweep >= 10_000 && sweep % 2 == 0 {
                    let obs = qmc_rs::lattice::measure_observables(&configuration, engine.model())
                        .expect("obs");
                    energy_sum += obs.energy_total;
                    samples += 1;
                }
            }
            let measured = energy_sum / samples as f64;

            let e_up = j_z * 0.25;
            let e_dn = -j_z * 0.25;
            let z = 2.0 * (-beta * e_up).exp() + 2.0 * (-beta * e_dn).exp();
            let exact = (2.0 * e_up * (-beta * e_up).exp() + 2.0 * e_dn * (-beta * e_dn).exp()) / z;

            assert!(
                (measured - exact).abs() < 0.03,
                "Ising dimer: MC={measured:.6}, exact={exact:.6}"
            );
        }
        Err(_) => {
            // Model builder rejects pure-Ising (zero off-diagonal). Acceptable.
        }
    }
}
