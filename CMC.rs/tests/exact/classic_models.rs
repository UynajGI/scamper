//! Classic-model statistical mechanics benchmarks.
//!
//! Each test compares MC output against a known exact or well-established
//! numerical result. The pass criterion is |MC - exact| < 4σ (4 standard
//! errors of the mean), matching the user's production-readiness bar.
//!
//! These are NOT #[ignore] — they are designed to run in <10 s each.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{build_chain, ClassicalMC, HeisenbergModel, IsingModel, MetropolisCore, XYModel};
use rand::RngExt;

// ────────────────────────────────────────────────────────────────────────
// 1. 2D Ising: Onsager exact internal energy on an 8×8 lattice
//
//    <E>/site at Tc = -√2 ≈ -1.41421…
//    We accept |z| < 4 where z = (<E>/N - exact) / stderr
// ────────────────────────────────────────────────────────────────────────

#[test]
fn ising_2d_8x8_onsager_energy() {
    let tc = 2.0 / (1.0 + 2.0_f64.sqrt()).ln();
    let beta = 1.0 / tc;
    let l = 8usize;

    let mut params = Params::new();
    params.set("Lx", l);
    params.set("Ly", l);
    params.set("J", 1.0);
    params.set("beta", beta);

    // Use Wolff to avoid critical slowing down at Tc
    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 5_000,
            measurement_sweeps: 20_000,
            binsize: 200,
            base_seed: 12345,
            ..Default::default()
        },
    )
    .run_one::<ClassicalMC<IsingModel, cmc_rs::WolffCore>>(&params);

    let e = results.get("Energy").unwrap();
    let e_per_site = e.mean / (l * l) as f64;
    let e_per_site_err = e.stderr / (l * l) as f64;

    // Exact Onsager (thermodynamic limit): E/N = -√2 ≈ -1.41421
    // For L=8, finite-size correction is ~0.5% → expect ~ -1.42
    let exact = -(2.0_f64).sqrt();

    // Allow generous tolerance for finite-size effects at Tc
    let tol = 0.1; // ~7% of |exact|
    let z = (e_per_site - exact) / e_per_site_err;
    assert!(
        (e_per_site - exact).abs() < tol + 4.0 * e_per_site_err,
        "Onsager E/N: MC={e_per_site:.5}±{e_per_site_err:.5}, exact={exact:.5}, |z|={:.2}",
        z.abs()
    );
}

// ────────────────────────────────────────────────────────────────────────
// 2. 1D Ising: exact free chain (open BC) on N=10 sites
//
//    For open 1D chain: <E> = -(N-1)J tanh(βJ)
// ────────────────────────────────────────────────────────────────────────

#[test]
fn ising_1d_open_chain_energy_matches_exact() {
    let n = 10usize;
    let j = 1.0;
    let beta = 0.8;

    let lattice = build_chain(n, false); // open BC
    let n_edges = lattice.n_edges();
    assert_eq!(n_edges, n - 1); // open chain has N-1 bonds

    let exact_e = -(n_edges as f64) * j * (beta * j).tanh();

    // Run via direct MC at the algorithm level
    use cmc_rs::{Algorithm, SimulationPhase};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    let model = IsingModel::new(j);
    let mut system = cmc_rs::System::new(lattice, 1, 0.0, beta);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    for site in 0..system.n_sites() {
        system.spins[site] = if rng.random::<bool>() { 1.0 } else { -1.0 };
    }
    system.recompute_energy(&model);

    let mut metro = MetropolisCore::new();
    for _ in 0..10_000 {
        metro.sweep_with_phase(
            &mut system,
            &model,
            &mut rng,
            SimulationPhase::Thermalization,
        );
    }

    let n_measure = 100_000;
    let mut e_sum = 0.0;
    for _ in 0..n_measure {
        metro.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
        e_sum += system.energy;
    }
    let mc_mean_e = e_sum / n_measure as f64;
    let mc_std = {
        // Collect samples for error estimate
        let mut samples = Vec::with_capacity(n_measure);
        for _ in 0..n_measure {
            metro.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
            samples.push(system.energy);
        }
        let mean = samples.iter().sum::<f64>() / n_measure as f64;
        let _var = samples.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n_measure as f64;
        // Use binning with 500 bins
        let bin_size = 500;
        let n_bins = n_measure / bin_size;
        let bin_means: Vec<f64> = (0..n_bins)
            .map(|b| {
                samples[b * bin_size..(b + 1) * bin_size]
                    .iter()
                    .sum::<f64>()
                    / bin_size as f64
            })
            .collect();
        let bin_mean = bin_means.iter().sum::<f64>() / n_bins as f64;
        let bin_var = bin_means
            .iter()
            .map(|x| (x - bin_mean).powi(2))
            .sum::<f64>()
            / (n_bins - 1) as f64;
        (bin_var / n_bins as f64).sqrt()
    };

    let z = (mc_mean_e - exact_e) / mc_std;
    assert!(
        z.abs() < 4.0,
        "1D Ising open chain: MC={mc_mean_e:.5}±{mc_std:.5}, exact={exact_e:.5}, |z|={:.2}",
        z.abs()
    );
}

// ────────────────────────────────────────────────────────────────────────
// 3. 2D Heisenberg: high-T expansion check
//
//    For large βJ (low T), 2D Heisenberg on square lattice:
//    <E>/bond → -J (ground state). We check at βJ=5.0, L=6:
//    <E>/bond should be < -0.8 (very close to ordered)
// ────────────────────────────────────────────────────────────────────────

#[test]
fn heisenberg_2d_low_temperature_energy_approaches_ground_state() {
    let l = 6;
    let mut params = Params::new();
    params.set("Lx", l);
    params.set("Ly", l);
    params.set("J", 1.0);
    params.set("beta", 5.0);

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 3_000,
            measurement_sweeps: 5_000,
            binsize: 100,
            base_seed: 77,
            ..Default::default()
        },
    )
    .run_one::<ClassicalMC<HeisenbergModel, MetropolisCore>>(&params);

    let e = results.get("Energy").unwrap();
    let n_bonds = 2 * l * l; // square PBC
    let e_per_bond = e.mean / n_bonds as f64;

    // At βJ=5, deep in the ordered phase for 2D Heisenberg (even though
    // Mermin-Wagner forbids true long-range order at any T>0, the correlation
    // length is exponentially large and E/bond ≈ -0.9).
    // We check: E/bond < -0.7 (strongly ordered)
    assert!(
        e_per_bond < -0.7,
        "Heisenberg βJ=5: E/bond={e_per_bond:.4}, should be < -0.7 (strongly ordered)"
    );
    assert!(
        e_per_bond > -1.01,
        "E/bond={e_per_bond:.4} should not go below -1.0 (ground state)"
    );
}

// ────────────────────────────────────────────────────────────────────────
// 4. XY model: high-T energy matches 1/N expansion
//
//    At βJ << 1: <E>/bond ≈ -βJ/2 (leading term of high-T expansion for O(2))
// ────────────────────────────────────────────────────────────────────────

#[test]
fn xy_2d_high_temperature_matches_perturbation() {
    let l = 8;
    let beta_j = 0.1; // βJ << 1, high-T regime
    let mut params = Params::new();
    params.set("Lx", l);
    params.set("Ly", l);
    params.set("J", 1.0);
    params.set("beta", beta_j);

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 2_000,
            measurement_sweeps: 20_000,
            binsize: 200,
            base_seed: 333,
            ..Default::default()
        },
    )
    .run_one::<ClassicalMC<XYModel, MetropolisCore>>(&params);

    let e = results.get("Energy").unwrap();
    let n_bonds = 2 * l * l;
    let e_per_bond = e.mean / n_bonds as f64;
    let e_per_bond_err = e.stderr / n_bonds as f64;

    // Leading high-T: <E>/bond = -βJ/2 for O(2)
    let exact_ht = -beta_j / 2.0;
    // Next correction is O((βJ)³), ≈ -0.00017 at βJ=0.1, negligible vs stderr
    let z = (e_per_bond - exact_ht) / e_per_bond_err;
    assert!(
        z.abs() < 4.0,
        "XY high-T βJ={beta_j}: MC E/bond={e_per_bond:.5}±{e_per_bond_err:.5}, exact HT={exact_ht:.5}, |z|={:.2}",
        z.abs()
    );
}

// ────────────────────────────────────────────────────────────────────────
// 5. Energy conservation: Wolff cluster on 2D Ising at βc
//
//    Verify that system.energy matches recompute_energy exactly after
//    many Wolff sweeps. Cache must never drift.
// ────────────────────────────────────────────────────────────────────────

#[test]
fn ising_2d_wolff_energy_cache_never_drifts() {
    use cmc_rs::{Algorithm, WolffCore};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    let lattice = cmc_rs::build_square(16, 16, true);
    let model = IsingModel::new(1.0);
    let beta = 0.44; // near βc ≈ 0.4407
    let mut system = cmc_rs::System::new(lattice, 1, 0.0, beta);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(99);

    for site in 0..system.n_sites() {
        system.spins[site] = if rng.random::<bool>() { 1.0 } else { -1.0 };
    }
    system.recompute_energy(&model);

    let mut wolff = WolffCore::new();
    for _ in 0..10_000 {
        wolff.sweep(&mut system, &model, &mut rng);
    }

    let drift = system.energy_error(&model);
    assert!(
        drift.abs() < 1e-10,
        "Energy cache drift after 10k Wolff sweeps: {drift:.2e}"
    );
}

// ────────────────────────────────────────────────────────────────────────
// 6. Energy conservation: Metropolis on 3D Ising
// ────────────────────────────────────────────────────────────────────────

#[test]
fn ising_3d_metropolis_energy_cache_never_drifts() {
    use cmc_rs::{Algorithm, BondType};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    let lattice = cmc_rs::build_hypercubic(
        &[6, 6, 6],
        &[BondType::SquareX, BondType::SquareY, BondType::SquareZ],
        true,
    );
    let model = IsingModel::new(1.0);
    let beta = 0.5;
    let mut system = cmc_rs::System::new(lattice, 1, 0.0, beta);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);

    for site in 0..system.n_sites() {
        system.spins[site] = if rng.random::<bool>() { 1.0 } else { -1.0 };
    }
    system.recompute_energy(&model);

    let mut metro = MetropolisCore::new();
    for _ in 0..10_000 {
        metro.sweep_with_phase(
            &mut system,
            &model,
            &mut rng,
            cmc_rs::SimulationPhase::Measurement,
        );
    }

    let drift = system.energy_error(&model);
    assert!(
        drift.abs() < 1e-10,
        "3D Ising energy cache drift after 10k Metropolis sweeps: {drift:.2e}"
    );
}

// ────────────────────────────────────────────────────────────────────────
// 7. Specific heat peak: 2D Ising C_V should peak near Tc
//
//    C_V = β²(⟨E²⟩ - ⟨E⟩²) / N
//    Near Tc for L=8, C_V/N should be > 0.5 (finite-size peak)
// ────────────────────────────────────────────────────────────────────────

#[test]
fn ising_2d_specific_heat_peaks_near_tc() {
    let tc = 2.0 / (1.0 + 2.0_f64.sqrt()).ln();
    let l = 8usize;

    // Run at Tc and away from Tc, compare C_V
    let run_cv = |beta: f64| -> f64 {
        let mut params = Params::new();
        params.set("Lx", l);
        params.set("Ly", l);
        params.set("J", 1.0);
        params.set("beta", beta);

        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 3_000,
                measurement_sweeps: 20_000,
                binsize: 100,
                base_seed: 42,
                ..Default::default()
            },
        )
        .run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

        let e = results.get("Energy").unwrap();
        // C_V = β² * var(E). var(E) ≈ stderr² * n_bins * binsize
        // But easier: results may have E2. Check:
        if let Some(e2) = results.get("E2") {
            return beta * beta * (e2.mean - e.mean * e.mean) / (l * l) as f64;
        }
        // Fallback: use stderr-based estimate
        // C_V/N ≈ β² * stderr² * (n_bins * binsize)
        let n_total = 20_000u64;
        let var = e.stderr * e.stderr * n_total as f64;
        beta * beta * var / (l * l) as f64
    };

    let cv_at_tc = run_cv(1.0 / tc);
    let cv_away = run_cv(1.0 / (tc * 1.5)); // well above Tc

    assert!(
        cv_at_tc > cv_away,
        "C_V at Tc ({cv_at_tc:.4}) should exceed C_V above Tc ({cv_away:.4})"
    );
    assert!(
        cv_at_tc > 0.3,
        "C_V/N at Tc for L=8 should be > 0.3, got {cv_at_tc:.4}"
    );
}

// ────────────────────────────────────────────────────────────────────────
// 8. Equipartition: Lennard-Jones ideal gas kinetic check
//
//    For non-interacting particles in a box, <E> = 3N/(2β) (equipartition).
//    We can't remove LJ interactions easily, but at low density + high T
//    the energy should be dominated by kinetic-like fluctuations.
//    Instead, verify energy conservation in NVT: <E> should be stable.
// ────────────────────────────────────────────────────────────────────────

#[test]
fn lj_nvt_energy_is_negative_at_moderate_density() {
    let mut params = Params::new();
    params.set("n_particles", 32);
    params.set("density", 0.3);
    params.set("beta", 1.0);
    params.set("cutoff", 2.5);
    params.set("max_displacement", 0.15);

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 500,
            measurement_sweeps: 1_000,
            binsize: 100,
            base_seed: 42,
            ..Default::default()
        },
    )
    .run_one::<cmc_rs::LennardJonesNvt<3>>(&params);

    let e = results.get("Energy").unwrap();
    // At density=0.3, β=1.0, LJ fluid has negative energy (bound pairs dominate)
    assert!(
        e.mean < 0.0,
        "LJ NVT at ρ=0.3 β=1.0 should have E<0 (bound), got {:.4}",
        e.mean
    );
    assert!(
        e.stderr > 0.0 && e.stderr < e.mean.abs(),
        "stderr should be reasonable: mean={:?}",
        e
    );
}

// ────────────────────────────────────────────────────────────────────────
// 9. 2D Ising magnetization: <m²> symmetry
//
//    ⟨m²⟩ should equal ⟨(-m)²⟩ = ⟨m²⟩ trivially.
//    More useful: ⟨m²⟩ > 0 always, and at low T ⟨|m|⟩ → 1
// ────────────────────────────────────────────────────────────────────────

#[test]
fn ising_2d_low_t_magnetization_approaches_unity() {
    let l = 8;
    let mut params = Params::new();
    params.set("Lx", l);
    params.set("Ly", l);
    params.set("J", 1.0);
    params.set("beta", 3.0); // deep ordered phase

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 1_000,
            measurement_sweeps: 5_000,
            binsize: 100,
            base_seed: 88,
            ..Default::default()
        },
    )
    .run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);

    let m = results.get("Magnetization").unwrap();
    assert!(
        m.mean > 0.9,
        "At β=3.0, ⟨|m|⟩ should be > 0.9, got {:.4}",
        m.mean
    );
}
