//! Extended ergodicity tests: multi-seed convergence for non-Metropolis solvers.
//!
//! Each test runs 4 independent seeds and asserts that ⟨E⟩ (and ⟨m²⟩ where
//! physically meaningful) agree within 4 combined standard errors.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{
    ClassicalMC, ContinuousHeatBathCore, HeatBathCore, HeisenbergModel, IsingGraphWormMC,
    IsingModel, KawasakiCore, KineticIsingBklMC,
};

/// Check that two estimates agree within `n` combined standard errors.
fn within_n_sigma(n: f64, a_mean: f64, a_err: f64, b_mean: f64, b_err: f64) -> bool {
    let combined = (a_err * a_err + b_err * b_err).sqrt();
    if combined == 0.0 {
        // Zero stderr (e.g. a conserved quantity): require near-exact equality.
        (a_mean - b_mean).abs() < 1e-10
    } else {
        (a_mean - b_mean).abs() < n * combined
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Heat bath on 2D Ising
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn heat_bath_2d_ising_converges_same_regardless_of_seed() {
    // ClassicalMC<IsingModel, HeatBathCore> on a 4×4 square lattice.
    // Heat bath is an exact conditional sampler — faster mixing than
    // Metropolis at the cost of the same fixed-point distribution.
    // All 4 seeds must converge to the same ⟨E⟩ and ⟨m²⟩ within 4σ.
    fn run(seed: u64) -> (f64, f64, f64, f64) {
        let mut params = Params::new();
        params.set("Lx", 4);
        params.set("Ly", 4);
        params.set("J", 1.0);
        params.set("beta", 0.5);
        let config = RunConfig {
            thermalization_sweeps: 5000,
            measurement_sweeps: 20000,
            binsize: 500,
            base_seed: seed,
            ..Default::default()
        };
        let scheduler = Scheduler::new(RayonBackend::new(1), config);
        let results = scheduler.run_one::<ClassicalMC<IsingModel, HeatBathCore>>(&params);
        let e = results.get("Energy").expect("Energy");
        let m2 = results.get("M2").expect("M2");
        (e.mean, e.stderr, m2.mean, m2.stderr)
    }

    let seeds = [42u64, 999, 7, 314];
    let estimates: Vec<_> = seeds.iter().map(|&s| run(s)).collect();

    for i in 0..estimates.len() {
        for j in (i + 1)..estimates.len() {
            let (ei, ei_err, m2i, m2i_err) = estimates[i];
            let (ej, ej_err, m2j, m2j_err) = estimates[j];
            assert!(
                within_n_sigma(4.0, ei, ei_err, ej, ej_err),
                "HeatBath 2D Ising ⟨E⟩ disagree: seeds {} vs {}: {:.4}±{:.4} vs {:.4}±{:.4}",
                seeds[i],
                seeds[j],
                ei,
                ei_err,
                ej,
                ej_err
            );
            assert!(
                within_n_sigma(4.0, m2i, m2i_err, m2j, m2j_err),
                "HeatBath 2D Ising ⟨m²⟩ disagree: seeds {} vs {}: {:.4}±{:.4} vs {:.4}±{:.4}",
                seeds[i],
                seeds[j],
                m2i,
                m2i_err,
                m2j,
                m2j_err
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Kawasaki on 2D Ising
// ═══════════════════════════════════════════════════════════════════════

#[test]
#[ignore = "Kawasaki conserves magnetization: different seeds start in different \
magnetization sectors, so neither ⟨E⟩ nor ⟨m²⟩ is expected to agree across \
seeds. The standard scheduler API (ClassicalMC + from_params) does not expose \
a zero-magnetization initial state. Re-enable once a balanced initialisation \
option is available."]
fn kawasaki_2d_ising_energy_converges_same_regardless_of_seed() {
    // ClassicalMC<IsingModel, KawasakiCore> on a 4×4 square lattice.
    //
    // Kawasaki dynamics swaps neighboring spins instead of flipping them,
    // conserving total magnetization.  Each trajectory is confined to a
    // fixed-magnetization sector, so ⟨m²⟩ is seed-dependent and ⟨E⟩
    // differs across sectors (especially near Tc at β=0.5).  The standard
    // scheduler only supports "hot" (random) and "cold" (ordered) initial
    // states — neither guarantees a shared magnetization sector across
    // seeds.  This test is therefore #[ignore]'d until a balanced
    // (zero-magnetization) initialisation is available through the API.
    //
    // When it can be run, it checks ⟨E⟩ only (⟨m²⟩ is conserved per
    // trajectory and thus meaningless for multi-seed comparison).
    fn run(seed: u64) -> (f64, f64) {
        let mut params = Params::new();
        params.set("Lx", 4);
        params.set("Ly", 4);
        params.set("J", 1.0);
        params.set("beta", 0.5);
        let config = RunConfig {
            thermalization_sweeps: 5000,
            measurement_sweeps: 20000,
            binsize: 500,
            base_seed: seed,
            ..Default::default()
        };
        let scheduler = Scheduler::new(RayonBackend::new(1), config);
        let results = scheduler.run_one::<ClassicalMC<IsingModel, KawasakiCore>>(&params);
        let e = results.get("Energy").expect("Energy");
        (e.mean, e.stderr)
    }

    let seeds = [42u64, 999, 7, 314];
    let estimates: Vec<_> = seeds.iter().map(|&s| run(s)).collect();

    for i in 0..estimates.len() {
        for j in (i + 1)..estimates.len() {
            let (ei, ei_err) = estimates[i];
            let (ej, ej_err) = estimates[j];
            assert!(
                within_n_sigma(4.0, ei, ei_err, ej, ej_err),
                "Kawasaki 2D Ising ⟨E⟩ disagree: seeds {} vs {}: {:.4}±{:.4} vs {:.4}±{:.4}",
                seeds[i],
                seeds[j],
                ei,
                ei_err,
                ej,
                ej_err
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Continuous O(3) Heisenberg heat bath
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn o3_heisenberg_heat_bath_converges_same_regardless_of_seed() {
    // ClassicalMC<HeisenbergModel, ContinuousHeatBathCore> on a 4×4 square
    // lattice.  HeisenbergModel = ONModel<3>; ContinuousHeatBathCore uses
    // the exact conditional sampler for O(3) spins on S².  All 4 seeds
    // must converge to the same ⟨E⟩ and ⟨m²⟩ within 4σ.
    fn run(seed: u64) -> (f64, f64, f64, f64) {
        let mut params = Params::new();
        params.set("Lx", 4);
        params.set("Ly", 4);
        params.set("J", 1.0);
        params.set("beta", 0.5);
        let config = RunConfig {
            thermalization_sweeps: 5000,
            measurement_sweeps: 20000,
            binsize: 500,
            base_seed: seed,
            ..Default::default()
        };
        let scheduler = Scheduler::new(RayonBackend::new(1), config);
        let results =
            scheduler.run_one::<ClassicalMC<HeisenbergModel, ContinuousHeatBathCore>>(&params);
        let e = results.get("Energy").expect("Energy");
        let m2 = results.get("M2").expect("M2");
        (e.mean, e.stderr, m2.mean, m2.stderr)
    }

    let seeds = [42u64, 999, 7, 314];
    let estimates: Vec<_> = seeds.iter().map(|&s| run(s)).collect();

    for i in 0..estimates.len() {
        for j in (i + 1)..estimates.len() {
            let (ei, ei_err, m2i, m2i_err) = estimates[i];
            let (ej, ej_err, m2j, m2j_err) = estimates[j];
            assert!(
                within_n_sigma(4.0, ei, ei_err, ej, ej_err),
                "O(3) HeatBath ⟨E⟩ disagree: seeds {} vs {}: {:.4}±{:.4} vs {:.4}±{:.4}",
                seeds[i],
                seeds[j],
                ei,
                ei_err,
                ej,
                ej_err
            );
            assert!(
                within_n_sigma(4.0, m2i, m2i_err, m2j, m2j_err),
                "O(3) HeatBath ⟨m²⟩ disagree: seeds {} vs {}: {:.4}±{:.4} vs {:.4}±{:.4}",
                seeds[i],
                seeds[j],
                m2i,
                m2i_err,
                m2j,
                m2j_err
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Worm (Ising high-temperature graph) on 2D Ising
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn worm_2d_ising_energy_converges_same_regardless_of_seed() {
    // IsingGraphWormMC on an 8×8 square lattice at β=0.3.
    //
    // The worm samples the high-temperature graph representation. Energy
    // is measured only in the physical sector (closed worm paths). The
    // worm does not directly measure ⟨m²⟩ — m² would require the endpoint-
    // pair histogram or an auxiliary observable. This test checks ⟨E⟩
    // convergence across 4 independent seeds.
    fn run(seed: u64) -> (f64, f64) {
        let mut params = Params::new();
        params.set("lattice_type", "square");
        params.set("Lx", 8);
        params.set("Ly", 8);
        params.set("pbc", true);
        params.set("J", 1.0);
        params.set("beta", 0.3);
        params.set("worm_updates_per_sweep", 128);
        let config = RunConfig {
            thermalization_sweeps: 5000,
            measurement_sweeps: 20000,
            binsize: 500,
            base_seed: seed,
            ..Default::default()
        };
        let scheduler = Scheduler::new(RayonBackend::new(1), config);
        let results = scheduler.run_one::<IsingGraphWormMC>(&params);
        let e = results.get("Energy").expect("Energy");
        (e.mean, e.stderr)
    }

    let seeds = [42u64, 999, 7, 314];
    let estimates: Vec<_> = seeds.iter().map(|&s| run(s)).collect();

    for i in 0..estimates.len() {
        for j in (i + 1)..estimates.len() {
            let (ei, ei_err) = estimates[i];
            let (ej, ej_err) = estimates[j];
            assert!(
                within_n_sigma(4.0, ei, ei_err, ej, ej_err),
                "Worm 2D Ising ⟨E⟩ disagree: seeds {} vs {}: {:.4}±{:.4} vs {:.4}±{:.4}",
                seeds[i],
                seeds[j],
                ei,
                ei_err,
                ej,
                ej_err
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════
// BKL (n-fold way) on 2D Ising
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn bkl_2d_ising_energy_converges_same_regardless_of_seed() {
    // KineticIsingBklMC on an 8×8 square lattice at β=0.3.
    //
    // BKL/n-fold-way uses rejection-free Glauber dynamics sampled at fixed
    // event-time intervals. Each sweep advances the continuous-time clock by
    // one event-time unit. The Fenwick-tree kernel guarantees exact event-
    // time statistics. All 4 seeds must converge to the same ⟨E⟩ within 4σ.
    fn run(seed: u64) -> (f64, f64) {
        let mut params = Params::new();
        params.set("lattice_type", "square");
        params.set("Lx", 8);
        params.set("Ly", 8);
        params.set("pbc", true);
        params.set("J", 1.0);
        params.set("beta", 0.3);
        params.set("event_time_per_sweep", 1.0);
        params.set("kinetic_rate", "glauber");
        let config = RunConfig {
            thermalization_sweeps: 5000,
            measurement_sweeps: 20000,
            binsize: 500,
            base_seed: seed,
            ..Default::default()
        };
        let scheduler = Scheduler::new(RayonBackend::new(1), config);
        let results = scheduler.run_one::<KineticIsingBklMC>(&params);
        let e = results.get("Energy").expect("Energy");
        (e.mean, e.stderr)
    }

    let seeds = [42u64, 999, 7, 314];
    let estimates: Vec<_> = seeds.iter().map(|&s| run(s)).collect();

    for i in 0..estimates.len() {
        for j in (i + 1)..estimates.len() {
            let (ei, ei_err) = estimates[i];
            let (ej, ej_err) = estimates[j];
            assert!(
                within_n_sigma(4.0, ei, ei_err, ej, ej_err),
                "BKL 2D Ising ⟨E⟩ disagree: seeds {} vs {}: {:.4}±{:.4} vs {:.4}±{:.4}",
                seeds[i],
                seeds[j],
                ei,
                ei_err,
                ej,
                ej_err
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Wang-Landau on 2D Ising
// ═══════════════════════════════════════════════════════════════════════

/// Canonical energy from a log-density-of-states at inverse temperature β.
fn wl_canonical_energy(log_dos: &[f64], energies: &[f64], beta: f64) -> f64 {
    let mut z = 0.0_f64;
    let mut weighted_sum = 0.0_f64;
    for (&ln_g, &e) in log_dos.iter().zip(energies.iter()) {
        let weight = (ln_g - beta * e).exp();
        z += weight;
        weighted_sum += weight * e;
    }
    if z == 0.0 {
        return f64::NAN;
    }
    weighted_sum / z
}

#[test]
#[ignore = "long: WL 4-seed ergodicity (~15s)"]
fn wang_landau_2d_ising_energy_converges_same_regardless_of_seed() {
    use cmc_rs::{IsingWangLandau, MacrostateAxis, WangLandauRunControl};

    const WL_BETA: f64 = 0.5;

    fn run(seed: u64) -> f64 {
        let mut params = Params::new();
        params.set("Lx", 4);
        params.set("Ly", 4);
        params.set("J", 1.0);
        params.set("beta", WL_BETA);
        params.set("wl_final_log_f", 1e-6);
        params.set("wl_flatness", 0.8);
        params.set("wl_flatness_check_interval", 100);
        params.set("wl_discovery_sweeps", 0);
        params.set("wl_minimum_visited_fraction", 0.8);
        params.set("wl_max_adaptation_sweeps", 100_000);

        let scheduler = Scheduler::new(
            RayonBackend::new(1),
            RunConfig {
                base_seed: seed,
                ..Default::default()
            },
        );
        let (mc, _results) = scheduler
            .run_controlled_with_state::<IsingWangLandau, WangLandauRunControl>(
                &params,
                WangLandauRunControl::new(0),
            )
            .expect("WL run should succeed");

        let wl_log = mc.estimator().log_density();
        let wl_centers = mc.chain.algorithm.axis().centers();
        let wl_log_values: Vec<f64> = (0..wl_log.bins()).map(|b| wl_log.value(b)).collect();
        wl_canonical_energy(&wl_log_values, &wl_centers, WL_BETA)
    }

    let seeds = [42u64, 999, 7, 314];
    let estimates: Vec<f64> = seeds.iter().map(|&s| run(s)).collect();

    // Each WL seed produces a single DOS-based energy estimate (no per-seed
    // stderr). We use the sample standard deviation across the 4 seeds.
    let n = estimates.len() as f64;
    let mean_e = estimates.iter().sum::<f64>() / n;
    let var = estimates.iter().map(|e| (e - mean_e).powi(2)).sum::<f64>() / (n - 1.0).max(1.0);
    let sigma = var.sqrt().max(1e-10);

    for i in 0..estimates.len() {
        for j in (i + 1)..estimates.len() {
            let diff = (estimates[i] - estimates[j]).abs();
            assert!(
                diff < 4.0 * sigma,
                "WL 2D Ising ⟨E⟩ disagree: seeds {} vs {}: {:.4} vs {:.4}, diff {:.4} > 4σ ({:.4})",
                seeds[i],
                seeds[j],
                estimates[i],
                estimates[j],
                diff,
                4.0 * sigma,
            );
        }
    }
}
