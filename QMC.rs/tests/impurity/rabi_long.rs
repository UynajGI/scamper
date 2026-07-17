//! Rabi model QMC validation: long stochastic test at large β.
//!
//! Runs the occupation worldline sampler at β=2^10 across a range of
//! coupling ratios r = g/g_c and compares ⟨n⟩ against exact ED.
//! This validates the QMC program produces correct physics across
//! the symmetric-to-broken crossover region.
//!
//! #[ignore] — takes minutes to run.

#![allow(clippy::needless_range_loop)]

use qmc_rs::impurity::spin_boson::occupation::{
    model::{CavityMode, OccupationSpinBosonModel},
    transfer::SymmetricEigensystem,
    OccupationWorldlineSampler,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Run MC at a given parameter point and return (⟨n⟩_MC, num_samples).
fn run_mc(delta: f64, omega: f64, g: f64, cutoff: usize, beta: f64, seed: u64) -> (f64, u64) {
    let model =
        OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
            .unwrap();
    let mut sampler = OccupationWorldlineSampler::new(model, beta, 8, 0).unwrap();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    let thermal = 50_000;
    let measure = 200_000;
    let mut n_sum = 0.0;
    let mut samples = 0u64;
    for sweep in 0..(thermal + measure) {
        sampler.sweep(&mut rng).unwrap();
        if sweep >= thermal {
            let obs = sampler.measure().unwrap();
            n_sum += obs.total_boson_number;
            samples += 1;
        }
    }
    (n_sum / samples as f64, samples)
}

/// Exact ED ⟨n⟩ at β→∞ (ground state).
fn ed_photon_number(delta: f64, omega: f64, g: f64, cutoff: usize) -> f64 {
    let model =
        OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
            .unwrap();
    let dim = model.basis().dimension();
    let eigen = SymmetricEigensystem::diagonalize(model.hamiltonian()).unwrap();
    let mut n = 0.0;
    for state in 0..dim {
        n += eigen.vectors[state][0].powi(2) * model.basis().occupation(state, 0) as f64;
    }
    n
}

/// Exact ED ⟨n⟩ at finite β (thermal).
fn ed_photon_number_finite_beta(delta: f64, omega: f64, g: f64, cutoff: usize, beta: f64) -> f64 {
    let model =
        OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
            .unwrap();
    let dim = model.basis().dimension();
    let eigen = SymmetricEigensystem::diagonalize(model.hamiltonian()).unwrap();
    let ground = eigen.values[0];
    let mut n_weighted = 0.0;
    let mut z = 0.0;
    for k in 0..dim {
        let boltz = (-beta * (eigen.values[k] - ground)).exp();
        let mut n_k = 0.0;
        for state in 0..dim {
            n_k += eigen.vectors[state][k].powi(2) * model.basis().occupation(state, 0) as f64;
        }
        n_weighted += boltz * n_k;
        z += boltz;
    }
    n_weighted / z
}

/// Scan r = g/g_c across [0.5, 3.0] and verify MC matches finite-β ED
/// at every point. 4σ tolerance.
#[test]
#[ignore]
fn mc_matches_ed_across_full_rabi_crossover() {
    let delta = 1.0_f64;
    let eta = 5.0_f64; // Ω/Δ = 5
    let omega = eta * delta;
    let gc = (omega * delta).sqrt() / 2.0;
    let cutoff = 60;
    let beta = (1u64 << 10) as f64; // β = 1024

    // Test points spanning weak → strong coupling
    let r_values = [0.5_f64, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 4.0, 5.0];

    eprintln!("η={eta}, β={beta}, cutoff={cutoff}");
    eprintln!("r      g       ⟨n⟩_ED     ⟨n⟩_MC     rel_err%");

    for (i, &r) in r_values.iter().enumerate() {
        let g = r * gc;
        let n_ed = ed_photon_number_finite_beta(delta, omega, g, cutoff, beta);
        let (n_mc, n_samples) = run_mc(delta, omega, g, cutoff, beta, 1000 + i as u64);

        let rel_err = if n_ed > 1e-10 {
            100.0 * (n_mc - n_ed).abs() / n_ed
        } else {
            if (n_mc - n_ed).abs() < 0.001 {
                0.0
            } else {
                100.0
            }
        };

        eprintln!("{r:.2}   {g:.4}   {n_ed:.6}   {n_mc:.6}   {rel_err:.2}%  ({n_samples} samples)");

        // At β=1024, the occupation solver (exact transfer matrix) should
        // match ED to within a few percent everywhere.
        assert!(
            rel_err < 5.0,
            "r={r:.2}: MC ⟨n⟩={n_mc:.6}, ED ⟨n⟩={n_ed:.6}, rel_err={rel_err:.2}%"
        );
    }
}

/// Verify that at large β and moderate η, the system is approaching
/// the ground state. At β→∞, thermal ⟨n⟩ should equal ground-state ⟨n⟩.
#[test]
#[ignore]
fn mc_converges_to_ground_state_at_large_beta() {
    let delta = 1.0_f64;
    let omega = 5.0_f64;
    let gc = (omega * delta).sqrt() / 2.0;
    let g = gc; // r = 1
    let cutoff = 40;

    let n_gs = ed_photon_number(delta, omega, g, cutoff);

    // β = 2^12 = 4096
    let beta = (1u64 << 12) as f64;
    let (n_mc, samples) = run_mc(delta, omega, g, cutoff, beta, 42);

    eprintln!("β={beta}: ⟨n⟩_GS={n_gs:.6}, ⟨n⟩_MC={n_mc:.6}, samples={samples}");

    let rel_err = 100.0 * (n_mc - n_gs).abs() / n_gs;
    assert!(
        rel_err < 5.0,
        "MC should converge to ground state: ⟨n⟩_MC={n_mc:.6}, ⟨n⟩_GS={n_gs:.6}, rel_err={rel_err:.2}%"
    );
}

/// Verify ⟨σz⟩ matches ED across the crossover.
#[test]
#[ignore]
fn mc_sigma_z_matches_ed_across_crossover() {
    let delta = 1.0_f64;
    let eta = 5.0_f64;
    let omega = eta * delta;
    let gc = (omega * delta).sqrt() / 2.0;
    let cutoff = 60;
    let beta = (1u64 << 10) as f64;

    let r_values = [0.5_f64, 1.0, 2.0, 3.0];

    eprintln!("σ_z validation: η={eta}, β={beta}");
    eprintln!("r      ⟨σz⟩_ED    ⟨σz⟩_MC    rel_err%");

    for (i, &r) in r_values.iter().enumerate() {
        let g = r * gc;

        // ED finite-β ⟨σz⟩
        let model_ed =
            OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
                .unwrap();
        let dim = model_ed.basis().dimension();
        let eigen = SymmetricEigensystem::diagonalize(model_ed.hamiltonian()).unwrap();
        let ground = eigen.values[0];
        let mut sz_weighted = 0.0;
        let mut z = 0.0;
        for k in 0..dim {
            let boltz = (-beta * (eigen.values[k] - ground)).exp();
            let mut sz_k = 0.0;
            for state in 0..dim {
                sz_k += eigen.vectors[state][k].powi(2) * model_ed.basis().spin(state).sigma_z();
            }
            sz_weighted += boltz * sz_k;
            z += boltz;
        }
        let sz_ed = sz_weighted / z;

        // MC
        let model_mc =
            OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
                .unwrap();
        let mut sampler = OccupationWorldlineSampler::new(model_mc, beta, 8, 0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(2000 + i as u64);
        let mut sz_sum = 0.0;
        let mut samples = 0u64;
        for sweep in 0..250_000u64 {
            sampler.sweep(&mut rng).unwrap();
            if sweep >= 50_000 {
                let obs = sampler.measure().unwrap();
                sz_sum += obs.sigma_z;
                samples += 1;
            }
        }
        let sz_mc = sz_sum / samples as f64;

        eprintln!(
            "{r:.2}   {sz_ed:.6}   {sz_mc:.6}   {:.2}%",
            100.0 * (sz_mc - sz_ed).abs() / sz_ed.abs().max(0.01)
        );

        // ⟨σz⟩ should match ED within 5%
        let rel_err_sz = 100.0 * (sz_mc - sz_ed).abs() / sz_ed.abs().max(0.01);
        assert!(
            rel_err_sz < 5.0,
            "r={r:.2}: ⟨σz⟩_MC={sz_mc:.4}, ⟨σz⟩_ED={sz_ed:.4}, rel_err={rel_err_sz:.2}%"
        );
    }
}

/// Verify energy matches ED across the crossover.
#[test]
#[ignore]
fn mc_energy_matches_ed_across_crossover() {
    let delta = 1.0_f64;
    let eta = 5.0_f64;
    let omega = eta * delta;
    let gc = (omega * delta).sqrt() / 2.0;
    let cutoff = 60;
    let beta = (1u64 << 10) as f64;

    let r_values = [0.5_f64, 1.0, 2.0, 3.0, 5.0];

    eprintln!("Energy validation: η={eta}, β={beta}");
    eprintln!("r      ⟨E⟩_ED      ⟨E⟩_MC      rel_err%");

    for (i, &r) in r_values.iter().enumerate() {
        let g = r * gc;

        // ED finite-β ⟨E⟩
        let model_ed =
            OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
                .unwrap();
        let dim = model_ed.basis().dimension();
        let eigen = SymmetricEigensystem::diagonalize(model_ed.hamiltonian()).unwrap();
        let ground = eigen.values[0];
        let mut e_weighted = 0.0;
        let mut z = 0.0;
        for k in 0..dim {
            let boltz = (-beta * (eigen.values[k] - ground)).exp();
            e_weighted += boltz * eigen.values[k];
            z += boltz;
        }
        let e_ed = e_weighted / z;

        // MC
        let model_mc =
            OccupationSpinBosonModel::rabi(delta, vec![CavityMode::new(omega, g, cutoff).unwrap()])
                .unwrap();
        let mut sampler = OccupationWorldlineSampler::new(model_mc, beta, 8, 0).unwrap();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(3000 + i as u64);
        let mut e_sum = 0.0;
        let mut samples = 0u64;
        for sweep in 0..250_000u64 {
            sampler.sweep(&mut rng).unwrap();
            if sweep >= 50_000 {
                let obs = sampler.measure().unwrap();
                e_sum += obs.energy;
                samples += 1;
            }
        }
        let e_mc = e_sum / samples as f64;

        let rel_err = 100.0 * (e_mc - e_ed).abs() / e_ed.abs().max(0.01);
        eprintln!("{r:.2}   {e_ed:.6}   {e_mc:.6}   {rel_err:.2}%");

        assert!(
            rel_err < 5.0,
            "r={r:.2}: ⟨E⟩_MC={e_mc:.6}, ⟨E⟩_ED={e_ed:.6}, rel_err={rel_err:.2}%"
        );
    }
}
