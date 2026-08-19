//! Production validation of `MultiSpinIsing` (64-replica bit-parallel Ising).
//!
//! Closes the two named gaps: **no cross-solver validation** and **no
//! multi-seed z-score coverage**.
//!
//! - **Multi-seed exact z-scores** (criterion B/E statistics): 8-site PBC
//!   chain, β ∈ {0.4, 0.8}, 8 seeds; ⟨E⟩, ⟨m²⟩ and the specific heat
//!   C = β²(⟨E²⟩ − ⟨E⟩²) vs full 2^8 enumeration — |z| < 4 per seed, |z̄| < 2
//!   per (β, observable), pooled one-sided Σz gate.
//! - **All 64 replicas ensemble-consistent**: the per-replica array
//!   observable `Energy_replica` averages over replicas *and* time; its
//!   per-seed mean matches the exact ⟨E⟩ (the replicas are exchangeable
//!   valid chains — replica 0 is not special).
//! - **Cross-solver** (criterion F): MultiSpinIsing vs scalar per-site
//!   `MetropolisCore` on identical physics (same lattice, β, J) — pooled
//!   cross-solver z on ⟨E⟩, ⟨|m|⟩ and ⟨m²⟩ at both temperatures.

use super::common::{exact_ising_moments, zscore_seed_count};
use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};

/// Binned point estimate (mean, stderr).
type Estimate = (f64, f64);
use cmc_rs::{build_chain, ClassicalMC, IsingModel, MetropolisCore, MultiSpinIsing};

const BETAS: [f64; 2] = [0.4, 0.8];
const THERM: u64 = 5_000;
const MEAS: u64 = 40_000;
const BINSIZE: usize = 500;

/// One MultiSpinIsing scheduler run: (E, se), (M2, se), (E2, se),
/// (|m|, se), (replica-averaged energy mean).
fn run_multispin(beta: f64, seed: u64) -> (Estimate, Estimate, Estimate, Estimate, f64) {
    let mut params = Params::new();
    params.set("L", 8_usize);
    params.set("J", 1.0_f64);
    params.set("beta", beta);
    let config = RunConfig {
        thermalization_sweeps: THERM,
        measurement_sweeps: MEAS,
        binsize: BINSIZE,
        base_seed: seed,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config).run_one::<MultiSpinIsing>(&params);
    let energy = results.get("Energy").expect("Energy");
    let m2 = results.get("M2").expect("M2");
    let e2 = results.get("E2").expect("E2");
    let magnetization = results.get("Magnetization").expect("Magnetization");
    let replica_energy = results
        .get("Energy_replica")
        .expect("Energy_replica array observable");
    (
        (energy.mean, energy.stderr),
        (m2.mean, m2.stderr),
        (e2.mean, e2.stderr),
        (magnetization.mean, magnetization.stderr),
        replica_energy.mean,
    )
}

/// One scalar Metropolis run: (E, se), (M2, se), (|m|, se).
fn run_scalar_metropolis(beta: f64, seed: u64) -> (Estimate, Estimate, Estimate) {
    let mut params = Params::new();
    params.set("L", 8_usize);
    params.set("J", 1.0_f64);
    params.set("beta", beta);
    let config = RunConfig {
        thermalization_sweeps: THERM,
        measurement_sweeps: MEAS,
        binsize: BINSIZE,
        base_seed: seed,
        ..Default::default()
    };
    let results = Scheduler::new(RayonBackend::new(1), config)
        .run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);
    let energy = results.get("Energy").expect("Energy");
    let m2 = results.get("M2").expect("M2");
    let magnetization = results.get("Magnetization").expect("Magnetization");
    (
        (energy.mean, energy.stderr),
        (m2.mean, m2.stderr),
        (magnetization.mean, magnetization.stderr),
    )
}

#[test]
fn multispin_multi_seed_z_scores_match_exact_enumeration() {
    let lattice = build_chain(8, true);
    let n_seeds = zscore_seed_count(8);
    let mut all_z: Vec<f64> = Vec::new();
    for &beta in &BETAS {
        let (_, exact_e, exact_e2, exact_m2) = exact_ising_moments(&lattice, 1.0, beta);
        let exact_c = beta * beta * (exact_e2 - exact_e * exact_e);
        let mut z_max = 0.0f64;
        let mut z_sum = [0.0f64; 4];
        let mut replica_bias = 0.0f64;
        for seed in 0..n_seeds as u64 {
            let ((e, e_se), (m2, m2_se), (e2, e2_se), _, replica_mean) =
                run_multispin(beta, 0x475 + seed);
            let c = beta * beta * (e2 - e * e);
            let c_se = beta * beta * (e2_se * e2_se + (2.0 * e * e_se) * (2.0 * e * e_se)).sqrt();
            let zs = [
                (e - exact_e) / e_se.max(1e-10),
                (m2 - exact_m2) / m2_se.max(1e-10),
                (c - exact_c) / c_se.max(1e-10),
            ];
            for (index, z) in zs.iter().enumerate() {
                z_max = z_max.max(z.abs());
                z_sum[index] += z;
            }
            // The 64-replica array observable averages to the same ⟨E⟩.
            replica_bias = replica_bias.max((replica_mean - exact_e).abs());
            z_sum[3] += (replica_mean - exact_e) / e_se.max(1e-10);
            all_z.extend(zs.iter().copied());
        }
        let n = n_seeds as f64;
        let mean_z: Vec<f64> = z_sum.iter().map(|z| z / n).collect();
        eprintln!(
            "[multispin β={beta}] exact E={exact_e:.5} m²={exact_m2:.5} C={exact_c:.5} | \
             max|z|={z_max:.2} z̄=[E {:+.2}, m² {:+.2}, C {:+.2}, replicas {:+.2}] max \
             replica-mean bias {replica_bias:.4}",
            mean_z[0], mean_z[1], mean_z[2], mean_z[3]
        );
        assert!(
            z_max < 4.0,
            "β={beta}: MultiSpinIsing max |z| = {z_max:.2} vs exact enumeration"
        );
        assert!(
            mean_z.iter().all(|z| z.abs() < 2.0),
            "β={beta}: mean z = {mean_z:?}"
        );
        // The replica-averaged energy is within one pooled exact-value sigma
        // of the same reference (all 64 replicas are valid exchangeable
        // chains, not just replica 0).
        assert!(
            mean_z[3].abs() < 1.0,
            "β={beta}: replica-averaged ⟨E⟩ biased, z̄ = {:.2}",
            mean_z[3]
        );
    }
    let sum_z: f64 = all_z.iter().sum();
    let n = all_z.len() as f64;
    eprintln!("[multispin pooled] n = {n}, Σz = {sum_z:+.2}");
    assert!(
        sum_z > -2.0 * n.sqrt(),
        "pooled Σz = {sum_z:.2} below one-sided gate"
    );
}

#[test]
fn multispin_matches_scalar_metropolis_cross_solver() {
    let n_seeds = zscore_seed_count(8);
    for &beta in &BETAS {
        let multispin: Vec<(Estimate, Estimate, Estimate)> = (0..n_seeds as u64)
            .map(|seed| {
                let ((e, e_se), (m2, m2_se), _, (mag, mag_se), _) =
                    run_multispin(beta, 0x5EED + seed);
                ((e, e_se), (m2, m2_se), (mag, mag_se))
            })
            .collect();
        let scalar: Vec<(Estimate, Estimate, Estimate)> = (0..n_seeds as u64)
            .map(|seed| run_scalar_metropolis(beta, 0x5EEE + seed))
            .collect();
        let pool = |values: &[Estimate]| {
            let n = values.len() as f64;
            let mean = values.iter().map(|(m, _)| m).sum::<f64>() / n;
            let sem = values.iter().map(|(_, s)| s * s).sum::<f64>().sqrt() / n;
            (mean, sem)
        };
        let multi_e: Vec<_> = multispin.iter().map(|(e, _, _)| *e).collect();
        let multi_m2: Vec<_> = multispin.iter().map(|(_, m, _)| *m).collect();
        let multi_mag: Vec<_> = multispin.iter().map(|(_, _, m)| *m).collect();
        let scalar_e: Vec<_> = scalar.iter().map(|(e, _, _)| *e).collect();
        let scalar_m2: Vec<_> = scalar.iter().map(|(_, m, _)| *m).collect();
        let scalar_mag: Vec<_> = scalar.iter().map(|(_, _, m)| *m).collect();

        let (me, me_se) = pool(&multi_e);
        let (se, se_se) = pool(&scalar_e);
        let (mm, mm_se) = pool(&multi_m2);
        let (sm, sm_se) = pool(&scalar_m2);
        let (mg, mg_se) = pool(&multi_mag);
        let (sg, sg_se) = pool(&scalar_mag);
        let z_e = (me - se) / (me_se * me_se + se_se * se_se).sqrt();
        let z_m2 = (mm - sm) / (mm_se * mm_se + sm_se * sm_se).sqrt();
        let z_mag = (mg - sg) / (mg_se * mg_se + sg_se * sg_se).sqrt();
        eprintln!(
            "[multispin-cross β={beta}] E: multi {me:.4}±{me_se:.4} scalar {se:.4}±{se_se:.4} \
             z={z_e:+.2} | m²: multi {mm:.4}±{mm_se:.4} scalar {sm:.4}±{sm_se:.4} z={z_m2:+.2} | \
             |m|: multi {mg:.4}±{mg_se:.4} scalar {sg:.4}±{sg_se:.4} z={z_mag:+.2}"
        );
        assert!(
            z_e.abs() < 4.0,
            "β={beta}: MultiSpinIsing vs scalar Metropolis ⟨E⟩ pooled-z = {z_e:.2}"
        );
        assert!(
            z_m2.abs() < 4.0,
            "β={beta}: MultiSpinIsing vs scalar Metropolis ⟨m²⟩ pooled-z = {z_m2:.2}"
        );
        assert!(
            z_mag.abs() < 4.0,
            "β={beta}: MultiSpinIsing vs scalar Metropolis ⟨|m|⟩ pooled-z = {z_mag:.2}"
        );
    }
}
