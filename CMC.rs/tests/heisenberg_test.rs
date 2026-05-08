//! Heisenberg model unit tests.

use cmc_rs::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn run_sim<T: MonteCarlo + FromParams>(
    params: &Params,
    thermalization: u64,
    measurements: u64,
    binsize: usize,
) -> Results {
    let backend = RayonBackend::new(1);
    let config = RunConfig {
        thermalization_sweeps: thermalization,
        measurement_sweeps: measurements,
        binsize,
        base_seed: 42,
        progress_interval: 0,
        checkpoint_interval: 0,
    };
    let scheduler = Scheduler::new(backend, config);
    scheduler.run_one::<T>(params)
}

fn make_params(l: usize, beta: f64, j: f64, pbc: bool) -> Params {
    let mut p = Params::new();
    p.set("L", l);
    p.set("beta", beta);
    p.set("J", j);
    p.set("pbc", pbc);
    p
}

#[test]
fn test_heisenberg_ground_state() {
    // All spins +z → each bond contributes -J (S·S = 1)
    let lattice = cmc_rs::lattice::build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);
    assert!((model.total_energy() - (-4.0)).abs() < 1e-10);
}

#[test]
fn test_heisenberg_spin_norm() {
    let lattice = cmc_rs::lattice::build_chain(8, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);
    for i in 0..8 {
        let sx = model.spins()[3 * i];
        let sy = model.spins()[3 * i + 1];
        let sz = model.spins()[3 * i + 2];
        let norm = (sx * sx + sy * sy + sz * sz).sqrt();
        assert!((norm - 1.0).abs() < 1e-10, "Spin {} has norm {}", i, norm);
    }
}

#[test]
fn test_heisenberg_magnetization_all_z() {
    let lattice = cmc_rs::lattice::build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);
    assert!((model.magnetization() - 1.0).abs() < 1e-10);
}

#[test]
fn test_heisenberg_metropolis_sweep() {
    let lattice = cmc_rs::lattice::build_chain(8, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);
    let mut core = MetropolisCore::new(model);
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 100);

    for _ in 0..100 {
        core.sweep(&mut ctx);
    }
    // Spins must remain unit vectors
    for i in 0..8 {
        let sx = core.model().spins()[3 * i];
        let sy = core.model().spins()[3 * i + 1];
        let sz = core.model().spins()[3 * i + 2];
        let norm = (sx * sx + sy * sy + sz * sz).sqrt();
        assert!((norm - 1.0).abs() < 0.01, "Spin {} norm = {}", i, norm);
    }
}

#[test]
fn test_heisenberg_from_params() {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut params = Params::new();
    params.set("L", 8usize);
    params.set("beta", 2.0f64);
    params.set("J", 1.0f64);

    let model = HeisenbergModel::from_params(&params, &mut rng).unwrap();
    assert_eq!(model.n_sites(), 8);
    assert_eq!(model.spin_dim(), 3);
}

#[test]
fn test_heisenberg_energy_extensive() {
    let params = make_params(16, 1.0, 1.0, true);
    let params2 = make_params(64, 1.0, 1.0, true);

    let res1 = run_sim::<MetropolisCore<HeisenbergModel>>(&params, 1000, 5000, 100);
    let res2 = run_sim::<MetropolisCore<HeisenbergModel>>(&params2, 2000, 8000, 100);

    let e1 = res1.get("Energy").unwrap().mean / 16.0;
    let e2 = res2.get("Energy").unwrap().mean / 64.0;

    assert!(
        (e1 - e2).abs() < 0.10 * e1.abs(),
        "Per-site energy should be similar: {} vs {}",
        e1, e2
    );
}

#[test]
fn test_heisenberg_local_energy_change_spin() {
    // 4-site PBC chain, all spins +z. Flip site 0 to -z.
    // Site 0 has 2 outgoing neighbors (sites 1, 3), each with spin +z.
    // Each bond: S_old·S_neighbor = 1, S_new·S_neighbor = -1
    // delta_e = -J * [(-1 - 1) + (-1 - 1)] = -J * (-4) = 4J = 4.0
    let lattice = cmc_rs::lattice::build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);

    let old = vec![0.0, 0.0, 1.0]; // +z
    let new = vec![0.0, 0.0, -1.0]; // -z
    let de = model.local_energy_change_spin(0, &old, &new);
    assert!((de - 4.0).abs() < 1e-10, "delta_e = {}, expected 4.0", de);

    // Verify against total_energy difference:
    // E_before = -4.0 (all +z)
    // E_after: site 0 is -z, others +z
    // Bonds: (0,1): -1, (1,2): 1, (2,3): 1, (3,0): -1 → sum = 0, /2 = 0.0
    // So delta_E = 0.0 - (-4.0) = 4.0
    let e_before = model.total_energy();
    assert!((e_before - (-4.0)).abs() < 1e-10);

    // After flipping site 0: energy = 0.0
    let mut model_flipped = HeisenbergModel::new(
        cmc_rs::lattice::build_chain(4, true), 1.0, 1.0, std::f64::consts::PI / 8.0
    );
    model_flipped.spins_mut()[0] = 0.0;
    model_flipped.spins_mut()[1] = 0.0;
    model_flipped.spins_mut()[2] = -1.0;
    let e_after = model_flipped.total_energy();
    let actual_delta = e_after - e_before;
    assert!(
        (actual_delta - de).abs() < 1e-10,
        "total_energy delta {} doesn't match local_energy_change_spin {}",
        actual_delta, de
    );
}

#[test]
fn test_heisenberg_local_vs_total_energy_consistency() {
    let lattice = cmc_rs::lattice::build_chain(8, true);
    let mut model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);

    // Set non-trivial spin configuration
    for i in 0..8 {
        let theta = i as f64 * std::f64::consts::PI / 4.0;
        model.spins_mut()[3 * i] = theta.sin();    // x
        model.spins_mut()[3 * i + 1] = 0.0;         // y
        model.spins_mut()[3 * i + 2] = theta.cos(); // z
    }
    // Normalize all spins
    for i in 0..8 {
        let norm = (model.spins()[3 * i].powi(2) + model.spins()[3 * i + 1].powi(2) + model.spins()[3 * i + 2].powi(2)).sqrt();
        model.spins_mut()[3 * i] /= norm;
        model.spins_mut()[3 * i + 1] /= norm;
        model.spins_mut()[3 * i + 2] /= norm;
    }

    // Test each site: flip (negate all components) and verify consistency
    for site in 0..8 {
        let old = vec![model.spins()[3 * site], model.spins()[3 * site + 1], model.spins()[3 * site + 2]];
        let new = vec![-old[0], -old[1], -old[2]];

        let de_local = model.local_energy_change_spin(site, &old, &new);

        let e_before = model.total_energy();
        model.spins_mut()[3 * site] = new[0];
        model.spins_mut()[3 * site + 1] = new[1];
        model.spins_mut()[3 * site + 2] = new[2];
        let e_after = model.total_energy();
        let de_global = e_after - e_before;

        assert!(
            (de_local - de_global).abs() < 1e-8,
            "Site {}: local ΔE={} vs global ΔE={}",
            site, de_local, de_global
        );

        // Restore
        model.spins_mut()[3 * site] = old[0];
        model.spins_mut()[3 * site + 1] = old[1];
        model.spins_mut()[3 * site + 2] = old[2];
    }
}

/// Type alias for OPSS Metropolis simulation (HeisenbergModel + OPSSStrategy).
type OPSSMetropolis = MetropolisCore<HeisenbergModel, OPSSStrategy>;

#[test]
fn test_opss_does_not_crash() {
    let params = make_params(8, 1.0, 1.0, true);
    let results = run_sim::<OPSSMetropolis>(&params, 500, 2000, 50);
    let energy = results.get("Energy").unwrap();
    assert!(energy.mean < 0.0, "Energy should be negative for ferromagnet");
}

#[test]
fn test_opss_vs_naive_agree() {
    // OPSS and naive Metropolis should give similar energies at moderate temperature
    let params = make_params(16, 1.0, 1.0, true);

    let naive = run_sim::<MetropolisCore<HeisenbergModel>>(&params, 2000, 5000, 100);
    let opss = run_sim::<OPSSMetropolis>(&params, 2000, 5000, 100);

    let e_naive = naive.get("Energy").unwrap().mean / 16.0;
    let e_opss = opss.get("Energy").unwrap().mean / 16.0;

    assert!(
        (e_naive - e_opss).abs() < 0.15 * e_naive.abs(),
        "OPSS per-site energy {} vs naive {} differ too much",
        e_opss, e_naive
    );
}

#[test]
fn test_opss_high_temperature() {
    let params = make_params(32, 0.01, 1.0, true);
    let results = run_sim::<OPSSMetropolis>(&params, 500, 5000, 100);
    let energy = results.get("Energy").unwrap();
    assert!(
        energy.mean.abs() < 2.0,
        "High-T energy should be near zero, got {}",
        energy.mean
    );
}

#[test]
fn test_opss_sigma_bounds() {
    // Run OPSS for a few sweeps and verify sigma stays bounded.
    // Very high T → near-100% acceptance → tests upper bound (sigma capped at 60.0).
    let params = make_params(8, 0.01, 1.0, true);
    let results = run_sim::<OPSSMetropolis>(&params, 100, 200, 50);
    let energy = results.get("Energy").unwrap();
    assert!(energy.mean.is_finite(), "Energy should be finite (sigma bounded)");
    assert!(!energy.mean.is_nan(), "Energy should not be NaN");
}
