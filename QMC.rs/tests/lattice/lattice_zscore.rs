//! Multi-seed z-score test for the lattice directed-loop solver.
//!
//! Runs a 3-site AFM Heisenberg chain at β=3.0 from 4 independent seeds
//! (default), computes the z-score of the MC energy against the
//! exact-diagonalization reference, and verifies statistical consistency:
//!   - |z| < 4 for each individual seed
//!   - mean |z| < 2 (no systematic bias)
//!
//! Setting `SCUTTLE_ZSCORE_SEEDS=<n>` raises the seed count for nightly
//! high-power monitoring (unset → the default 4 seeds, unchanged for CI).

use qmc_rs::lattice::ContinuousLatticeEngine;
use qmc_rs::{
    CsrGraph, EdgeCoupling, LatticeConfiguration, QmcKernel, SpinModelBuilder, SpinSpace,
    UpdateSchedule,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::lattice_ed::build_hamiltonian;
use crate::zscore_seeds::zscore_seeds;

/// Batch-means standard error of the mean.
fn batch_means_stderr(samples: &[f64], binsize: usize) -> f64 {
    let n_bins = samples.len() / binsize;
    assert!(n_bins >= 4, "need at least 4 bins for stderr estimate");
    let bin_means: Vec<f64> = (0..n_bins)
        .map(|b| {
            let start = b * binsize;
            samples[start..start + binsize].iter().sum::<f64>() / binsize as f64
        })
        .collect();
    let grand_mean: f64 = bin_means.iter().sum::<f64>() / n_bins as f64;
    let variance: f64 = bin_means
        .iter()
        .map(|&m| (m - grand_mean).powi(2))
        .sum::<f64>()
        / (n_bins - 1) as f64;
    (variance / n_bins as f64).sqrt()
}

/// Run a single MC seed and return (mean_energy, stderr_energy).
fn run_lattice_seed(seed: u64) -> (f64, f64) {
    let n_sites = 3;
    let beta = 3.0;
    let j = 1.0;
    let n_thermalization = 20_000;
    let n_measurement = 80_000;
    let measure_interval = 2;
    let binsize = 100;

    let graph = CsrGraph::chain(n_sites, false).expect("graph");
    let space = SpinSpace::uniform(n_sites, 1).expect("space");
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::heisenberg(j))
        .build()
        .expect("model");
    let mut configuration = LatticeConfiguration::new(beta, vec![0, 1, 0], &model).expect("config");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(8, 4, 64));
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    let mut energy_samples: Vec<f64> = Vec::new();
    for sweep in 0..(n_thermalization + n_measurement) {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
        if sweep >= n_thermalization && sweep % measure_interval == 0 {
            let obs =
                qmc_rs::lattice::measure_observables(&configuration, engine.model()).expect("obs");
            energy_samples.push(obs.energy_total);
        }
    }

    let mean = energy_samples.iter().sum::<f64>() / energy_samples.len() as f64;
    let stderr = batch_means_stderr(&energy_samples, binsize);
    (mean, stderr)
}

#[test]
fn lattice_zscore_energy_4_seeds() {
    // ── ED reference ─────────────────────────────────────────────────────
    let n_sites = 3;
    let beta = 3.0;
    let j = 1.0;
    let edges_pair = [(0, 1), (1, 2)];

    let graph = CsrGraph::chain(n_sites, false).expect("graph");
    let weight = graph.edges().first().unwrap().weight;
    let edges: Vec<(usize, usize, EdgeCoupling)> = edges_pair
        .iter()
        .map(|&(si, sj)| (si, sj, EdgeCoupling::heisenberg(j * weight)))
        .collect();
    let h = build_hamiltonian(n_sites, &edges);
    let rho = h.expm_negative(beta);
    let z = rho.trace();
    let exact_energy = h.multiply(&rho).trace() / z;

    // ── MC runs ──────────────────────────────────────────────────────────
    let seeds = zscore_seeds(&[42u64, 123, 456, 789]);
    let results: Vec<(f64, f64)> = seeds.iter().map(|&s| run_lattice_seed(s)).collect();

    let z_scores: Vec<f64> = results
        .iter()
        .map(|&(mean, stderr)| (mean - exact_energy) / stderr.max(1e-10))
        .collect();

    // Each seed: |z| < 4
    for (i, &z) in z_scores.iter().enumerate() {
        assert!(
            z.abs() < 4.0,
            "Seed {}: z-score = {z:.3} (MC={:.6}, exact={exact_energy:.6}, stderr={:.6})",
            seeds[i],
            results[i].0,
            results[i].1
        );
    }

    // Mean |z| < 2 (no systematic bias)
    let mean_abs_z: f64 = z_scores.iter().map(|z| z.abs()).sum::<f64>() / z_scores.len() as f64;
    assert!(
        mean_abs_z < 2.0,
        "Mean |z| = {mean_abs_z:.3}, should be < 2 (z-scores: {z_scores:?})"
    );
}
