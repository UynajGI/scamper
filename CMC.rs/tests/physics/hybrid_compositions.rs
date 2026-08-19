//! Validation of every hybrid (composed) update actually offered by the API.
//!
//! `HybridCore<A, B>` composes any two `Algorithm<H>` kernels. The offered
//! kernel set per model is fixed by the capability traits:
//!
//! - **Ising** (`MetropolisCore`, `WolffCore`, `SWCore`, `HeatBathCore`): all
//!   six pairwise compositions validated against full 2^N enumeration on the
//!   4-site ring — ⟨E⟩, ⟨m²⟩ and the specific heat C = β²(⟨E²⟩ − ⟨E⟩²), 8
//!   seeds each, |z| < 4 per seed, |z̄| < 2 per combo, pooled one-sided Σz
//!   gate across the whole matrix.
//! - **Composition boundary semantics** (the MCMC.rs-combinator pattern):
//!   same-seed determinism, `Hybrid(A, B)` ≡ manual `A; B` sequencing
//!   bit-for-bit, `repetitions(2, 3)` ≡ `A; A; B; B; B`, and the nested
//!   combinator closure `Hybrid(A, Hybrid(B, C))` ≡ `A; B; C`.
//! - **Continuous (O(2))**: the hybrid-only continuous pairings
//!   (Wolff, SW), (Metropolis, ContinuousHeatBath) and (Metropolis, Wolff)
//!   vs a pure-Wolff reference on 8×8 (pooled cross-solver z on ⟨E⟩, ⟨m²⟩).
//!   (Metropolis, Microcanonical) — the over-relaxation composition — is
//!   validated exhaustively in `over_relaxation.rs`.

use super::common::{exact_ising_moments, zscore_seed_count};
use cmc_rs::{
    build_chain, Algorithm, ClassicalMC, HeatBathCore, HybridCore, IsingModel, MetropolisCore,
    ONModel, OPSSStrategy, SWCore, SimulationPhase, System, WolffCore,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

type Rng = Xoshiro256PlusPlus;

/// Binned point estimate (mean, stderr).
type Estimate = (f64, f64);

// ── Ising: every pairwise composition vs exact enumeration ─────────────────

/// Scheduler run of one hybrid composition on the 4-site PBC Ising ring.
/// Returns ((E, se), (M2, se), (E2, se)).
fn run_ising_hybrid<A, B>(beta: f64, seed: u64) -> ((f64, f64), (f64, f64), (f64, f64))
where
    A: Default + Algorithm<IsingModel>,
    B: Default + Algorithm<IsingModel>,
{
    let mut params = carlo_rs::Params::new();
    params.set("L", 4_usize);
    params.set("J", 1.0_f64);
    params.set("beta", beta);
    let config = carlo_rs::RunConfig {
        thermalization_sweeps: 2_000,
        measurement_sweeps: 24_000,
        binsize: 400,
        base_seed: seed,
        ..Default::default()
    };
    let results = carlo_rs::Scheduler::new(carlo_rs::RayonBackend::new(1), config)
        .run_one::<ClassicalMC<IsingModel, HybridCore<A, B>>>(&params);
    let energy = results.get("Energy").expect("Energy");
    let m2 = results.get("M2").expect("M2");
    let e2 = results.get("E2").expect("E2");
    (
        (energy.mean, energy.stderr),
        (m2.mean, m2.stderr),
        (e2.mean, e2.stderr),
    )
}

/// One combo's multi-seed z-scores for ⟨E⟩, ⟨m²⟩ and C = β²(⟨E²⟩−⟨E⟩²).
fn combo_z_scores<A, B>(
    beta: f64,
    exact: (f64, f64, f64, f64),
    base_seed: u64,
    n_seeds: usize,
) -> Vec<f64>
where
    A: Default + Algorithm<IsingModel>,
    B: Default + Algorithm<IsingModel>,
{
    let (_, exact_e, exact_e2, exact_m2) = exact;
    let exact_c = beta * beta * (exact_e2 - exact_e * exact_e);
    let mut scores = Vec::with_capacity(3 * n_seeds);
    for seed in 0..n_seeds as u64 {
        let ((e, e_se), (m2, m2_se), (e2, e2_se)) =
            run_ising_hybrid::<A, B>(beta, base_seed + seed);
        let c = beta * beta * (e2 - e * e);
        let c_se = beta * beta * (e2_se * e2_se + (2.0 * e * e_se) * (2.0 * e * e_se)).sqrt();
        scores.push((e - exact_e) / e_se.max(1e-10));
        scores.push((m2 - exact_m2) / m2_se.max(1e-10));
        scores.push((c - exact_c) / c_se.max(1e-10));
    }
    scores
}

#[test]
fn hybrid_ising_all_kernel_pairings_match_exact_enumeration() {
    let beta = 0.5;
    let lattice = build_chain(4, true);
    let exact = exact_ising_moments(&lattice, 1.0, beta);
    let n_seeds = zscore_seed_count(8);

    type Metro = MetropolisCore;
    let mut all_scores: Vec<f64> = Vec::new();
    let combos: Vec<(&str, Vec<f64>)> = vec![
        (
            "Metropolis+Wolff",
            combo_z_scores::<Metro, WolffCore>(beta, exact, 0xC011, n_seeds),
        ),
        (
            "Metropolis+SW",
            combo_z_scores::<Metro, SWCore>(beta, exact, 0xC012, n_seeds),
        ),
        (
            "Metropolis+HeatBath",
            combo_z_scores::<Metro, HeatBathCore>(beta, exact, 0xC013, n_seeds),
        ),
        (
            "Wolff+SW",
            combo_z_scores::<WolffCore, SWCore>(beta, exact, 0xC014, n_seeds),
        ),
        (
            "Wolff+HeatBath",
            combo_z_scores::<WolffCore, HeatBathCore>(beta, exact, 0xC015, n_seeds),
        ),
        (
            "SW+HeatBath",
            combo_z_scores::<SWCore, HeatBathCore>(beta, exact, 0xC016, n_seeds),
        ),
    ];
    for (label, scores) in &combos {
        let max_abs_z = scores.iter().fold(0.0f64, |acc, z| acc.max(z.abs()));
        let mean_z = scores.iter().sum::<f64>() / scores.len() as f64;
        eprintln!(
            "[hybrid-ising {label}] n_z = {}, max|z| = {max_abs_z:.2}, z̄ = {mean_z:+.2}",
            scores.len()
        );
        assert!(
            max_abs_z < 4.0,
            "{label}: max |z| = {max_abs_z:.2} vs exact (E={:.5}, m²={:.5})",
            exact.1,
            exact.3
        );
        assert!(mean_z.abs() < 2.0, "{label}: z̄ = {mean_z:.2}");
        all_scores.extend(scores.iter().copied());
    }

    // Pooled one-sided-bias gate over the whole matrix (repo convention:
    // Σz > -2√n flags a systematically low-biased sampler).
    let sum_z: f64 = all_scores.iter().sum();
    let n = all_scores.len() as f64;
    eprintln!(
        "[hybrid-ising pooled] n = {n}, Σz = {sum_z:+.2}, gate = −{:.2}",
        2.0 * n.sqrt()
    );
    assert!(
        sum_z > -2.0 * n.sqrt(),
        "pooled Σz = {sum_z:.2} below the one-sided gate"
    );
}

// ── Composition boundary semantics ─────────────────────────────────────────

fn ising_system(seed: u64, beta: f64) -> (System, IsingModel) {
    let model = IsingModel::new(1.0);
    let mut system = System::new(build_chain(6, true), 1, 0.0, beta);
    let mut rng = Rng::seed_from_u64(seed);
    for spin in &mut system.spins {
        *spin = if rng.random::<bool>() { 1.0 } else { -1.0 };
    }
    system.recompute_energy(&model);
    (system, model)
}

#[test]
fn hybrid_composition_boundaries_match_manual_sequences_exactly() {
    let beta = 0.7;

    // Determinism: identical construction + seed ⇒ identical trajectory.
    let mut left = HybridCore::new(MetropolisCore::new(), WolffCore::new());
    let mut right = HybridCore::new(MetropolisCore::new(), WolffCore::new());
    let (mut left_system, model) = ising_system(0xB0A, beta);
    let (mut right_system, _) = ising_system(0xB0A, beta);
    assert_eq!(left_system.spins, right_system.spins);
    let mut left_rng = Rng::seed_from_u64(0xB0B);
    let mut right_rng = Rng::seed_from_u64(0xB0B);
    for _ in 0..100 {
        left.sweep_with_phase(
            &mut left_system,
            &model,
            &mut left_rng,
            SimulationPhase::Measurement,
        );
        right.sweep_with_phase(
            &mut right_system,
            &model,
            &mut right_rng,
            SimulationPhase::Measurement,
        );
    }
    assert_eq!(left_system.spins, right_system.spins);

    // Hybrid(A, B) ≡ manual "A; B" with the same stream.
    let (mut hybrid_system, model) = ising_system(0xC0A, beta);
    let (mut manual_system, _) = ising_system(0xC0A, beta);
    let mut hybrid = HybridCore::new(MetropolisCore::new(), SWCore::new());
    let mut metro = MetropolisCore::new();
    let mut sw = SWCore::new();
    let mut hybrid_rng = Rng::seed_from_u64(0xC0B);
    let mut manual_rng = Rng::seed_from_u64(0xC0B);
    for _ in 0..50 {
        hybrid.sweep_with_phase(
            &mut hybrid_system,
            &model,
            &mut hybrid_rng,
            SimulationPhase::Measurement,
        );
        metro.sweep_with_phase(
            &mut manual_system,
            &model,
            &mut manual_rng,
            SimulationPhase::Measurement,
        );
        sw.sweep_with_phase(
            &mut manual_system,
            &model,
            &mut manual_rng,
            SimulationPhase::Measurement,
        );
    }
    assert_eq!(
        hybrid_system.spins, manual_system.spins,
        "Hybrid(A, B) sweep must equal manual A; B sequencing"
    );

    // repetitions(2, 3) ≡ A; A; B; B; B per sweep.
    let (mut hybrid_system, model) = ising_system(0xD0A, beta);
    let (mut manual_system, _) = ising_system(0xD0A, beta);
    let mut hybrid = HybridCore::new(WolffCore::new(), HeatBathCore::new()).repetitions(2, 3);
    let mut wolff = WolffCore::new();
    let mut heat_bath = HeatBathCore::new();
    let mut hybrid_rng = Rng::seed_from_u64(0xD0B);
    let mut manual_rng = Rng::seed_from_u64(0xD0B);
    for _ in 0..50 {
        hybrid.sweep_with_phase(
            &mut hybrid_system,
            &model,
            &mut hybrid_rng,
            SimulationPhase::Measurement,
        );
        for _ in 0..2 {
            wolff.sweep_with_phase(
                &mut manual_system,
                &model,
                &mut manual_rng,
                SimulationPhase::Measurement,
            );
        }
        for _ in 0..3 {
            heat_bath.sweep_with_phase(
                &mut manual_system,
                &model,
                &mut manual_rng,
                SimulationPhase::Measurement,
            );
        }
    }
    assert_eq!(
        hybrid_system.spins, manual_system.spins,
        "repetitions(2, 3) must equal A; A; B; B; B sequencing"
    );

    // Nested combinator closure: Hybrid(A, Hybrid(B, C)) ≡ A; B; C.
    let (mut hybrid_system, model) = ising_system(0xE0A, beta);
    let (mut manual_system, _) = ising_system(0xE0A, beta);
    let mut hybrid = HybridCore::new(
        MetropolisCore::new(),
        HybridCore::new(WolffCore::new(), SWCore::new()),
    );
    let mut metro = MetropolisCore::new();
    let mut wolff = WolffCore::new();
    let mut sw = SWCore::new();
    let mut hybrid_rng = Rng::seed_from_u64(0xE0B);
    let mut manual_rng = Rng::seed_from_u64(0xE0B);
    for _ in 0..50 {
        hybrid.sweep_with_phase(
            &mut hybrid_system,
            &model,
            &mut hybrid_rng,
            SimulationPhase::Measurement,
        );
        metro.sweep_with_phase(
            &mut manual_system,
            &model,
            &mut manual_rng,
            SimulationPhase::Measurement,
        );
        wolff.sweep_with_phase(
            &mut manual_system,
            &model,
            &mut manual_rng,
            SimulationPhase::Measurement,
        );
        sw.sweep_with_phase(
            &mut manual_system,
            &model,
            &mut manual_rng,
            SimulationPhase::Measurement,
        );
    }
    assert_eq!(
        hybrid_system.spins, manual_system.spins,
        "nested Hybrid(A, Hybrid(B, C)) must equal A; B; C sequencing"
    );
}

// ── Continuous compositions vs a pure-Wolff reference (8×8 O(2)) ───────────

fn run_8x8_xy<A>(beta: f64, seed: u64) -> (Estimate, Estimate)
where
    A: Default + Algorithm<ONModel<2>>,
    ClassicalMC<ONModel<2>, A>: carlo_rs::FromParams,
{
    let mut params = carlo_rs::Params::new();
    params.set("Lx", 8_usize);
    params.set("Ly", 8_usize);
    params.set("J", 1.0_f64);
    params.set("beta", beta);
    let config = carlo_rs::RunConfig {
        thermalization_sweeps: 3_000,
        measurement_sweeps: 2_000,
        binsize: 100,
        base_seed: seed,
        ..Default::default()
    };
    let results = carlo_rs::Scheduler::new(carlo_rs::RayonBackend::new(1), config)
        .run_one::<ClassicalMC<ONModel<2>, A>>(&params);
    let energy = results.get("Energy").expect("Energy");
    let m2 = results.get("M2").expect("M2");
    ((energy.mean, energy.stderr), (m2.mean, m2.stderr))
}

fn assert_pooled_agreement(
    hybrid: &[(Estimate, Estimate)],
    wolff: &[(Estimate, Estimate)],
    label: &str,
) {
    let pool = |values: &[Estimate]| {
        let n = values.len() as f64;
        let mean = values.iter().map(|(m, _)| m).sum::<f64>() / n;
        let sem = values.iter().map(|(_, s)| s * s).sum::<f64>().sqrt() / n;
        (mean, sem)
    };
    let hybrid_energy: Vec<_> = hybrid.iter().map(|(e, _)| *e).collect();
    let wolff_energy: Vec<_> = wolff.iter().map(|(e, _)| *e).collect();
    let hybrid_m2: Vec<_> = hybrid.iter().map(|(_, m)| *m).collect();
    let wolff_m2: Vec<_> = wolff.iter().map(|(_, m)| *m).collect();
    let (he, he_se) = pool(&hybrid_energy);
    let (we, we_se) = pool(&wolff_energy);
    let (hm, hm_se) = pool(&hybrid_m2);
    let (wm, wm_se) = pool(&wolff_m2);
    let z_e = (he - we) / (he_se * he_se + we_se * we_se).sqrt();
    let z_m = (hm - wm) / (hm_se * hm_se + wm_se * wm_se).sqrt();
    eprintln!(
        "[hybrid-xy {label}] E: {he:.4}±{he_se:.4} vs wolff {we:.4}±{we_se:.4} z={z_e:.2} | \
         m²: {hm:.4}±{hm_se:.4} vs wolff {wm:.4}±{wm_se:.4} z={z_m:.2}"
    );
    assert!(z_e.abs() < 4.0, "{label}: ⟨E⟩ pooled-z = {z_e:.2}");
    assert!(z_m.abs() < 4.0, "{label}: ⟨m²⟩ pooled-z = {z_m:.2}");
}

#[test]
fn hybrid_continuous_compositions_match_wolff_reference() {
    let beta = 0.90;
    let n_seeds = zscore_seed_count(8);
    let wolff: Vec<_> = (0..n_seeds as u64)
        .map(|seed| run_8x8_xy::<WolffCore>(beta, 0xA100 + seed))
        .collect();
    assert_pooled_agreement(
        &(0..n_seeds as u64)
            .map(|seed| run_8x8_xy::<HybridCore<WolffCore, SWCore>>(beta, 0xA200 + seed))
            .collect::<Vec<_>>(),
        &wolff,
        "Wolff+SW",
    );
    assert_pooled_agreement(
        &(0..n_seeds as u64)
            .map(|seed| {
                run_8x8_xy::<HybridCore<MetropolisCore<OPSSStrategy>, WolffCore>>(
                    beta,
                    0xA300 + seed,
                )
            })
            .collect::<Vec<_>>(),
        &wolff,
        "Metropolis+Wolff",
    );
    assert_pooled_agreement(
        &(0..n_seeds as u64)
            .map(|seed| {
                run_8x8_xy::<HybridCore<MetropolisCore<OPSSStrategy>, cmc_rs::ContinuousHeatBathCore>>(
                    beta,
                    0xA400 + seed,
                )
            })
            .collect::<Vec<_>>(),
        &wolff,
        "Metropolis+ContinuousHeatBath",
    );
    // (The Metropolis+Microcanonical continuous composition is validated in
    // over_relaxation.rs against quadrature, analytic limits and Wolff.)
}
