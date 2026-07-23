//! Long stochastic convergence tests — `#[ignore]` by default.
//!
//! Run with: `cargo test -p cmc-rs --test suite -- --ignored`
//!
//! These tests compare MC output to exact analytical results on systems
//! large enough that the calculation takes seconds to minutes. They are
//! the "gold standard" physics validation layer.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::*;

// ════════════════════════════════════════════════════════════════════════
// 1. FINITE-SIZE SCALING: Binder cumulant at Tc for 2D Ising
// ════════════════════════════════════════════════════════════════════════

#[test]
#[ignore = "long: Binder cumulant finite-size scaling (~30 s)"]
fn binder_cumulant_at_tc_matches_universal_value() {
    let tc = 2.0 / (1.0 + 2.0_f64.sqrt()).ln();
    let beta_c = 1.0 / tc;

    let run_binder = |l: usize| -> f64 {
        let mut params = Params::new();
        params.set("Lx", l);
        params.set("Ly", l);
        params.set("J", 1.0);
        params.set("beta", beta_c);

        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 10_000,
                measurement_sweeps: 80_000,
                binsize: 400,
                base_seed: 42,
                ..Default::default()
            },
        )
        .run_one::<ClassicalMC<IsingModel, WolffCore>>(&params);

        let m2 = results
            .get("M2")
            .or_else(|| results.get("Magnetization"))
            .unwrap();
        let m4 = results.get("M4").expect("M4 observable");

        1.0 - m4.mean / (3.0 * m2.mean * m2.mean)
    };

    let u16 = run_binder(16);
    let u32 = run_binder(32);

    assert!(
        (u32 - 0.610).abs() < 0.04,
        "Binder cumulant at Tc for L=32: U={u32:.4}, expected ≈0.610"
    );
    assert!(
        (u16 - u32).abs() < 0.05,
        "U should be size-independent at Tc: U(L=16)={u16:.4}, U(L=32)={u32:.4}"
    );
}

// ════════════════════════════════════════════════════════════════════════
// 2. FINITE-SIZE SCALING: susceptibility peak narrows with larger L
// ════════════════════════════════════════════════════════════════════════

#[test]
#[ignore = "long: susceptibility peak height grows with L (~60 s)"]
fn susceptibility_peak_height_grows_with_system_size() {
    let tc = 2.0 / (1.0 + 2.0_f64.sqrt()).ln();
    let beta_c = 1.0 / tc;

    let measure_chi = |l: usize| -> f64 {
        let mut params = Params::new();
        params.set("Lx", l);
        params.set("Ly", l);
        params.set("J", 1.0);
        params.set("beta", beta_c);

        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 5_000,
                measurement_sweeps: 40_000,
                binsize: 400,
                base_seed: 42,
                ..Default::default()
            },
        )
        .run_one::<ClassicalMC<IsingModel, WolffCore>>(&params);

        let m2 = results
            .get("M2")
            .or_else(|| results.get("Magnetization"))
            .unwrap();
        beta_c * (l * l) as f64 * m2.mean
    };

    // χ at Tc should grow as L^{2-η} with η=1/4 for 2D Ising → χ ∝ L^{7/4}
    let chi8 = measure_chi(8);
    let chi16 = measure_chi(16);

    // χ(16)/χ(8) ≈ (16/8)^{7/4} = 2^{1.75} ≈ 3.36
    let ratio = chi16 / chi8;
    assert!(
        ratio > 2.5 && ratio < 4.5,
        "Susceptibility ratio χ(16)/χ(8)={ratio:.2}, expected ≈3.36 (2D Ising η=1/4)"
    );
}

// ════════════════════════════════════════════════════════════════════════
// 3. WANG-LANDAU: converged DOS matches exact enumeration
// ════════════════════════════════════════════════════════════════════════

#[test] // was #[ignore] — runs in ~11s, fast enough for CI
fn wang_landau_dos_matches_exact_on_4x4_ising() {
    let lattice = build_square(4, 4, true);
    let model = IsingModel::new(1.0);
    let exact_dos = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    let exact_log = exact_dos.log_density().unwrap();

    let mut params = Params::new();
    params.set("Lx", 4);
    params.set("Ly", 4);
    params.set("J", 1.0);
    params.set("wl_final_log_f", 1e-6);
    params.set("wl_flatness", 0.8);
    params.set("wl_flatness_check_interval", 100);
    params.set("wl_discovery_sweeps", 0);
    params.set("wl_minimum_visited_fraction", 0.8);
    params.set("wl_max_adaptation_sweeps", 500000);

    let scheduler = Scheduler::new(RayonBackend::new(1), RunConfig::default());
    let (mc, _results) = scheduler
        .run_controlled_with_state::<IsingWangLandau, WangLandauRunControl>(
            &params,
            WangLandauRunControl::new(0),
        )
        .expect("WL run should succeed");

    let wl_log = mc.estimator().log_density();
    let n_bins = exact_log.bins();
    assert_eq!(wl_log.bins(), n_bins);

    let exact_max = (0..n_bins)
        .filter(|&b| exact_log.is_visited(b))
        .map(|b| exact_log.value(b))
        .fold(f64::NEG_INFINITY, f64::max);
    let wl_max = (0..n_bins)
        .filter(|&b| wl_log.is_visited(b))
        .map(|b| wl_log.value(b))
        .fold(f64::NEG_INFINITY, f64::max);

    let mut max_diff: f64 = 0.0;
    for bin in 0..n_bins {
        if exact_log.is_visited(bin) && wl_log.is_visited(bin) {
            let diff = ((exact_log.value(bin) - exact_max) - (wl_log.value(bin) - wl_max)).abs();
            max_diff = max_diff.max(diff);
        }
    }

    assert!(
        max_diff < 0.15,
        "WL DOS max |Δ ln g(E)| = {max_diff:.4}, expected < 0.15"
    );
}

// ════════════════════════════════════════════════════════════════════════
// 4. KAWASAKI DYNAMICS: high-T equilibration and energy ordering
// ════════════════════════════════════════════════════════════════════════

#[test]
#[ignore = "long: Kawasaki spin diffusion on 1D chain (~10 s)"]
fn kawasaki_1d_high_t_equilibrates_and_energy_is_conserved() {
    let mut params = Params::new();
    params.set("Lx", 32);
    params.set("Ly", 1);
    params.set("beta", 0.01);
    params.set("J", 1.0);

    let results = Scheduler::new(
        RayonBackend::new(1),
        RunConfig {
            thermalization_sweeps: 5_000,
            measurement_sweeps: 20_000,
            binsize: 200,
            base_seed: 42,
            ..Default::default()
        },
    )
    .run_one::<KawasakiIsingMC>(&params);

    let e = results.get("Energy").expect("Energy");
    let e_per_bond = e.mean / 32.0;
    assert!(
        e_per_bond.abs() < 0.1,
        "Kawasaki high-T: E/bond={e_per_bond:.4}, should be near 0"
    );
    assert!(e.stderr > 0.0);
}

#[test]
#[ignore = "long: Kawasaki energy decreases on cooling (~10 s)"]
fn kawasaki_energy_decreases_with_cooling() {
    let run_energy = |beta: f64| -> f64 {
        let mut params = Params::new();
        params.set("Lx", 16);
        params.set("Ly", 1);
        params.set("beta", beta);
        params.set("J", 1.0);

        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 2_000,
                measurement_sweeps: 10_000,
                binsize: 100,
                base_seed: 42,
                ..Default::default()
            },
        )
        .run_one::<KawasakiIsingMC>(&params);

        results.get("Energy").unwrap().mean
    };

    let e_hot = run_energy(0.1);
    let e_cold = run_energy(2.0);

    assert!(
        e_cold < e_hot,
        "Kawasaki E should decrease on cooling: E(β=2)={e_cold:.4} < E(β=0.1)={e_hot:.4}"
    );
}

// ════════════════════════════════════════════════════════════════════════
// 5. NPT: <V> increases when P decreases
// ════════════════════════════════════════════════════════════════════════

#[test]
#[ignore = "long: NPT volume-pressure consistency (~30 s)"]
fn npt_volume_increases_when_pressure_decreases() {
    let run_volume = |pressure: f64| -> f64 {
        let mut params = Params::new();
        params.set("n_particles", 16);
        params.set("density", 0.05);
        params.set("beta", 2.0);
        params.set("pressure", pressure);
        params.set("cutoff", 2.5);
        params.set("max_displacement", 0.2);
        params.set("max_volume_scale", 0.1);

        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 500,
                measurement_sweeps: 2000,
                binsize: 100,
                base_seed: 42,
                ..Default::default()
            },
        )
        .run_one::<LennardJonesNpt<3>>(&params);

        results.get("Volume").expect("Volume observable").mean
    };

    let v_low_p = run_volume(0.05);
    let v_high_p = run_volume(0.5);

    // Lower pressure → larger volume (directional check only;
    // quantitative ratio needs much longer sampling for LJ fluid)
    assert!(
        v_low_p > v_high_p,
        "NPT: V(P=0.05)={v_low_p:.2} should be > V(P=0.5)={v_high_p:.2}"
    );
}

// ════════════════════════════════════════════════════════════════════════
// 5b. NPT: finite-N ideal gas quantitative — ⟨V⟩ = (N+1)kT/P
// ════════════════════════════════════════════════════════════════════════

#[test]
#[ignore = "long: NPT finite-N ideal gas quantitative (~40 s)"]
fn npt_ideal_gas_volume_matches_finite_n_exact() {
    // Finite-N NPT ideal gas: exact ⟨V⟩ = (N+1)kT/P (not NkT/P).
    // The +1 comes from the volume integration measure dV in the partition function.
    // Thermodynamic limit recovers NkT/P as N→∞.
    //
    // Use box_length (not density) to control initial volume directly —
    // density-based box uses reference_particles=108 which creates a huge box
    // for small N.
    let run_volume = |pressure: f64, box_length: f64| -> f64 {
        let mut params = Params::new();
        params.set("n_particles", 4usize);
        params.set("box_length", box_length);
        params.set("beta", 1.0);
        params.set("pressure", pressure);
        params.set("cutoff", box_length * 0.5); // large cutoff → nearly ideal
        params.set("max_displacement", 0.3);
        params.set("max_log_volume_change", 0.5);

        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 10000,
                measurement_sweeps: 40000,
                binsize: 500,
                base_seed: 42,
                ..Default::default()
            },
        )
        .run_one::<LennardJonesNpt<3>>(&params);

        results.get("Volume").expect("Volume observable").mean
    };

    // Finite-N exact: ⟨V⟩ = (N+1)kT/P = 5/P for N=4
    // P=0.1 → 50, P=0.5 → 10.  Ratio = 5.0
    let v_low_p = run_volume(0.1, 3.68); // initial V ≈ 50
    let v_high_p = run_volume(0.5, 2.15); // initial V ≈ 10

    eprintln!("NPT finite-N: V(P=0.1)={v_low_p:.4}, V(P=0.5)={v_high_p:.4}");

    // At low density (large V/N ≈ 12.5), LJ repulsion is negligible,
    // so the system is nearly ideal. Expect ratio close to P2/P1 = 5.0.
    let ratio = v_low_p / v_high_p;
    let exact_ratio = 5.0 / 1.0; // P(high)/P(low) = 0.5/0.1
    assert!(
        (ratio - exact_ratio).abs() / exact_ratio < 0.25,
        "NPT: V(P=0.1)/V(P=0.5) = {ratio:.4}, exact finite-N = {exact_ratio:.1}"
    );
}

// ════════════════════════════════════════════════════════════════════════
// 6b. μVT: ideal gas quantitative — ⟨N⟩ = zV = exp(βμ)Λ⁻³V
// ════════════════════════════════════════════════════════════════════════

#[test]
#[ignore = "long: μVT ideal gas quantitative (~40 s)"]
fn muvt_ideal_gas_particle_number_matches_poisson_most_probable() {
    // Ideal gas μVT: exact ⟨N⟩ = zV where z = exp(βμ)Λ⁻³.
    // N follows Poisson(zV). Ratio N(μ=0)/N(μ=-2) = exp(2β) ≈ 7.39 for β=1.
    //
    // Use box_length (not density) — density-based box uses reference_particles=108
    // which gives a huge box for small N. Use large cutoff to suppress LJ repulsion.
    // thermal_wavelength=1.0 for simplicity (Λ=1).
    let run_n = |mu: f64, box_length: f64| -> f64 {
        let mut params = Params::new();
        params.set("n_particles", 4usize);
        params.set("box_length", box_length);
        params.set("beta", 1.0);
        params.set("sigma", 0.1); // tiny particles → nearly ideal gas
        params.set("cutoff", 0.25); // cutoff = 2.5*sigma
        params.set("max_displacement", 0.3);
        params.set("chemical_potential", mu);
        params.set("thermal_wavelength", 1.0);
        params.set("maximum_particles", 40usize);
        params.set("exchange_attempts", 4u64);

        let results = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                thermalization_sweeps: 4000,
                measurement_sweeps: 16000,
                binsize: 200,
                base_seed: 42,
                ..Default::default()
            },
        )
        .run_one::<LennardJonesMuVt<3>>(&params);

        results
            .get("ParticleNumber")
            .expect("ParticleNumber observable")
            .mean
    };

    // Box V=10: at μ=0, ⟨N⟩ = V × 1.0 × 1.0 = 10. At μ=-2, ⟨N⟩ = V × 0.135 × 1.0 ≈ 1.4.
    // Ratio ≈ 7.4. Use small box to keep particle counts manageable.
    let n_low_mu = run_n(-2.0, 2.154); // V ≈ 10
    let n_high_mu = run_n(0.0, 2.154); // same V

    eprintln!("μVT: N(μ=0)={n_high_mu:.4}, N(μ=-2)={n_low_mu:.4}");

    let ratio = n_high_mu / n_low_mu;
    let exact = (2.0f64).exp(); // exp(β(μ_high - μ_low)) = exp(2) ≈ 7.39
    assert!(
        (ratio - exact).abs() / exact < 0.25,
        "μVT: N(μ=0)/N(μ=-2) = {ratio:.4}, exact = {exact:.2}"
    );
}
