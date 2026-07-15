//! Statistical correctness tests for classical Monte Carlo sampling.
//!
//! Verifies that measured observables match exact enumeration for small systems.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{ClassicalMC, Hamiltonian, IsingModel, MetropolisCore, PottsModel, SWCore, WolffCore};

// ---------------------------------------------------------------------------
// Exact enumeration helpers
// ---------------------------------------------------------------------------

/// Boltzmann weight for a configuration: exp(-beta * energy)
fn boltzmann_weight(energy: f64, beta: f64) -> f64 {
    (-beta * energy).exp()
}

/// Enumerate all Ising configurations (2^N states) and compute exact <E>.
fn exact_ising_mean_energy(n: usize, j: f64, beta: f64, pbc: bool) -> f64 {
    let model = IsingModel::new(j);
    let lattice = cmc_rs::build_chain(n, pbc);
    let mut z = 0.0_f64;
    let mut weighted_e = 0.0_f64;
    for mask in 0..(1u32 << n) {
        let spins: Vec<f64> = (0..n)
            .map(|i| if (mask >> i) & 1 == 1 { 1.0 } else { -1.0 })
            .collect();
        let e = model.compute_total_energy(&spins, &lattice, 1.0);
        let w = boltzmann_weight(e, beta);
        z += w;
        weighted_e += e * w;
    }
    weighted_e / z
}

// ---------------------------------------------------------------------------
// A.1: Ising exact enumeration (N=2,3,4)
// ---------------------------------------------------------------------------

#[test]
fn exact_ising_n2_energy() {
    let exact = exact_ising_mean_energy(2, 1.0, 0.5, true);
    let mut params = Params::new();
    params.set("L", 2);
    params.set("J", 1.0);
    params.set("beta", 0.5);
    let config = RunConfig {
        thermalization_sweeps: 5000,
        measurement_sweeps: 10000,
        binsize: 200,
        base_seed: 42,
        ..Default::default()
    };
    let scheduler = Scheduler::new(RayonBackend::new(1), config);
    let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);
    let e = results.get("Energy").expect("Energy missing");
    assert!(
        (e.mean - exact).abs() < 3.0 * e.stderr.max(1e-4),
        "MC {:.6} ± {:.6}, exact {:.6}",
        e.mean,
        e.stderr,
        exact
    );
}

#[test]
fn exact_ising_n3_energy() {
    let exact = exact_ising_mean_energy(3, 1.0, 0.5, true);
    let mut params = Params::new();
    params.set("L", 3);
    params.set("J", 1.0);
    params.set("beta", 0.5);
    let config = RunConfig {
        thermalization_sweeps: 5000,
        measurement_sweeps: 10000,
        binsize: 200,
        base_seed: 42,
        ..Default::default()
    };
    let scheduler = Scheduler::new(RayonBackend::new(1), config);
    let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);
    let e = results.get("Energy").expect("Energy missing");
    assert!(
        (e.mean - exact).abs() < 3.0 * e.stderr.max(1e-4),
        "MC {:.6} ± {:.6}, exact {:.6}",
        e.mean,
        e.stderr,
        exact
    );
}

#[test]
fn exact_ising_n4_energy() {
    let exact = exact_ising_mean_energy(4, 1.0, 0.5, true);
    let mut params = Params::new();
    params.set("L", 4);
    params.set("J", 1.0);
    params.set("beta", 0.5);
    let config = RunConfig {
        thermalization_sweeps: 5000,
        measurement_sweeps: 10000,
        binsize: 200,
        base_seed: 42,
        ..Default::default()
    };
    let scheduler = Scheduler::new(RayonBackend::new(1), config);
    let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);
    let e = results.get("Energy").expect("Energy missing");
    assert!(
        (e.mean - exact).abs() < 3.0 * e.stderr.max(1e-4),
        "MC {:.6} ± {:.6}, exact {:.6}",
        e.mean,
        e.stderr,
        exact
    );
}

// ---------------------------------------------------------------------------
// A.2: Potts q=3 energy distribution on N=4 chain
// ---------------------------------------------------------------------------

#[test]
fn potts_q3_n4_energy_distribution() {
    let n = 4;
    let q = 3;
    let j = 1.0;
    let beta = 0.25;
    let model = PottsModel::new(j, q);
    let lattice = cmc_rs::build_chain(n, true);

    // Enumerate all q^N = 81 states
    let mut z = 0.0_f64;
    let mut weighted_e = 0.0_f64;
    for state in 0..(q as u32).pow(n as u32) {
        let mut s = state;
        let spins: Vec<f64> = (0..n)
            .map(|_| {
                let val = (s % q as u32) as f64;
                s /= q as u32;
                val
            })
            .collect();
        let e = model.compute_total_energy(&spins, &lattice, 1.0);
        let w = boltzmann_weight(e, beta);
        z += w;
        weighted_e += e * w;
    }
    let exact_mean = weighted_e / z;

    let mut params = Params::new();
    params.set("L", n);
    params.set("J", j);
    params.set("q", q);
    params.set("beta", beta);
    let config = RunConfig {
        thermalization_sweeps: 5000,
        measurement_sweeps: 10000,
        binsize: 200,
        base_seed: 42,
        ..Default::default()
    };
    let scheduler = Scheduler::new(RayonBackend::new(1), config);
    let results = scheduler.run_one::<ClassicalMC<PottsModel, MetropolisCore>>(&params);
    let e = results.get("Energy").expect("Energy missing");
    assert!(
        (e.mean - exact_mean).abs() < 3.0 * e.stderr.max(5e-3),
        "MC {:.6} ± {:.6}, exact {:.6}",
        e.mean,
        e.stderr,
        exact_mean
    );
}

// ---------------------------------------------------------------------------
// A.3: Algorithm consistency — Metropolis, Wolff, SW give same means
// ---------------------------------------------------------------------------

fn run_ising_8x8<A>(seed: u64) -> (f64, f64, f64, f64)
where
    A: cmc_rs::Algorithm<IsingModel> + Default,
{
    let mut params = Params::new();
    params.set("Lx", 8);
    params.set("Ly", 8);
    params.set("J", 1.0);
    params.set("beta", 0.44);
    let config = RunConfig {
        thermalization_sweeps: 2000,
        measurement_sweeps: 5000,
        binsize: 200,
        base_seed: seed,
        ..Default::default()
    };
    let scheduler = Scheduler::new(RayonBackend::new(1), config);
    let results = scheduler.run_one::<ClassicalMC<IsingModel, A>>(&params);
    let e = results.get("Energy").expect("Energy missing");
    let m = results.get("Magnetization").expect("Magnetization missing");
    (e.mean, e.stderr, m.mean, m.stderr)
}

#[test]
fn algorithm_consistency_ising_8x8() {
    let seed = 42;
    let (e_met, se_met, m_met, sm_met) = run_ising_8x8::<MetropolisCore>(seed);
    let (e_wolff, se_wolff, m_wolff, sm_wolff) = run_ising_8x8::<WolffCore>(seed);
    let (e_sw, se_sw, m_sw, sm_sw) = run_ising_8x8::<SWCore>(seed);

    // 3-sigma check: all pairs must be within 3 sigma
    for (label_a, ea, sea, ma, sma) in [
        ("met", e_met, se_met, m_met, sm_met),
        ("wolff", e_wolff, se_wolff, m_wolff, sm_wolff),
        ("sw", e_sw, se_sw, m_sw, sm_sw),
    ] {
        for (label_b, eb, seb, mb, smb) in [
            ("met", e_met, se_met, m_met, sm_met),
            ("wolff", e_wolff, se_wolff, m_wolff, sm_wolff),
            ("sw", e_sw, se_sw, m_sw, sm_sw),
        ] {
            if label_a == label_b {
                continue;
            }
            let e_diff = (ea - eb).abs();
            let e_sigma = (sea.powi(2) + seb.powi(2)).sqrt().max(1e-4);
            assert!(
                e_diff < 3.0 * e_sigma,
                "Energy {label_a} vs {label_b}: diff={e_diff:.6}, 3σ={:.6}",
                3.0 * e_sigma
            );
            let m_diff = (ma - mb).abs();
            let m_sigma = (sma.powi(2) + smb.powi(2)).sqrt().max(1e-4);
            assert!(
                m_diff < 3.0 * m_sigma,
                "Mag {label_a} vs {label_b}: diff={m_diff:.6}, 3σ={:.6}",
                3.0 * m_sigma
            );
        }
    }
}

// ---------------------------------------------------------------------------
// A.5: Fixed seed reproducibility
// ---------------------------------------------------------------------------

#[test]
fn fixed_seed_reproducibility() {
    let seed = 12345;
    let run = || {
        let mut params = Params::new();
        params.set("Lx", 4);
        params.set("Ly", 4);
        params.set("J", 1.0);
        params.set("beta", 0.5);
        let config = RunConfig {
            thermalization_sweeps: 100,
            measurement_sweeps: 200,
            binsize: 10,
            base_seed: seed,
            ..Default::default()
        };
        let scheduler = Scheduler::new(RayonBackend::new(1), config);
        let results = scheduler.run_one::<ClassicalMC<IsingModel, MetropolisCore>>(&params);
        let e = results.get("Energy").expect("Energy missing");
        (e.mean, e.stderr)
    };

    let (e1, se1) = run();
    let (e2, se2) = run();
    assert!(
        (e1 - e2).abs() < 1e-14,
        "Reproducibility failed: {e1} vs {e2}"
    );
    assert!(
        (se1 - se2).abs() < 1e-14,
        "Reproducibility failed: {se1} vs {se2}"
    );
}

// ---------------------------------------------------------------------------
// A.4: PT energy consistency under beta changes
// ---------------------------------------------------------------------------

#[test]
fn pt_energy_consistency_under_parameter_change() {
    use carlo_rs::{FromParams, ParallelTemperingCompatible};
    use rand::SeedableRng;

    let mut params = Params::new();
    params.set("L", 4usize);
    params.set("J", 1.0);
    params.set("beta", 1.0);

    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(12);
    let mc =
        <ClassicalMC<IsingModel, MetropolisCore> as FromParams>::from_params(&params, &mut rng)
            .unwrap();

    // log_weight_ratio for "beta" uses physical energy (not beta*E):
    // log_weight_ratio("beta", new_beta) = (beta_old - beta_new) * energy
    let beta_old = mc.system.beta;
    let beta_new = 2.0;
    let ratio = mc.log_weight_ratio("beta", beta_new);
    assert_eq!(ratio, (beta_old - beta_new) * mc.system.energy);

    // After change_parameter("beta", new_beta), energy must be recomputed
    // (physical energy is beta-independent, but cache is refreshed)
    let mut mc = mc;
    let energy_before = mc.system.energy;
    mc.change_parameter("beta", 2.0);
    assert_eq!(mc.system.beta, 2.0);
    // Physical energy should be unchanged (or recomputed to same value)
    assert!(
        (mc.system.energy - energy_before).abs() < 1e-12,
        "physical energy must be invariant under beta change"
    );
}
