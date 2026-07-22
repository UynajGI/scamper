//! Remaining validation tests: HybridCore correctness.
//!
//! NPT, μVT, event-chain, WL, Binder cumulant, and Kawasaki long
//! stochastic tests are in `long_convergence.rs`. This file covers
//! HybridCore, which needs a Default impl to work with the scheduler.

use carlo_rs::{Params, RayonBackend, RunConfig, Scheduler};
use cmc_rs::{ClassicalMC, Hamiltonian, IsingModel};

// ═══════════════════════════════════════════════════════════════════════
// P2.2: HybridCore correctness — now works with Default
// ═══════════════════════════════════════════════════════════════════════

fn exact_energy(n: usize, j: f64, beta: f64, pbc: bool) -> f64 {
    let lattice = cmc_rs::build_chain(n, pbc);
    let model = IsingModel::new(j);
    let mut z = 0.0;
    let mut we = 0.0;
    for mask in 0..(1u32 << n) {
        let spins: Vec<f64> = (0..n)
            .map(|i| if (mask >> i) & 1 == 1 { 1.0 } else { -1.0 })
            .collect();
        let e = model.compute_total_energy(&spins, &lattice, 1.0);
        let w = (-beta * e).exp();
        z += w;
        we += e * w;
    }
    we / z
}

fn exact_m2(n: usize, j: f64, beta: f64, pbc: bool) -> f64 {
    let lattice = cmc_rs::build_chain(n, pbc);
    let model = IsingModel::new(j);
    let mut z = 0.0;
    let mut wm2 = 0.0;
    for mask in 0..(1u32 << n) {
        let spins: Vec<f64> = (0..n)
            .map(|i| if (mask >> i) & 1 == 1 { 1.0 } else { -1.0 })
            .collect();
        let e = model.compute_total_energy(&spins, &lattice, 1.0);
        let m: f64 = spins.iter().sum::<f64>() / n as f64;
        let w = (-beta * e).exp();
        z += w;
        wm2 += m * m * w;
    }
    wm2 / z
}

#[test]
fn hybrid_metropolis_wolff_matches_exact_energy_and_magnetization() {
    let n = 4;
    let beta = 0.5;
    let j = 1.0;
    let exact_e = exact_energy(n, j, beta, true);
    let exact_m2 = exact_m2(n, j, beta, true);

    let mut params = Params::new();
    params.set("L", n);
    params.set("J", j);
    params.set("beta", beta);
    let config = RunConfig {
        thermalization_sweeps: 5000,
        measurement_sweeps: 20000,
        binsize: 500,
        base_seed: 77,
        ..Default::default()
    };
    let scheduler = Scheduler::new(RayonBackend::new(1), config);
    let results = scheduler.run_one::<
        ClassicalMC<IsingModel, cmc_rs::HybridCore<cmc_rs::MetropolisCore, cmc_rs::WolffCore>>,
    >(&params);

    let e = results.get("Energy").expect("Energy");
    assert!(
        (e.mean - exact_e).abs() < 3.0 * e.stderr.max(0.1),
        "Hybrid E={:.4} ± {:.4}, exact={:.4}",
        e.mean,
        e.stderr,
        exact_e
    );

    let m2 = results
        .get("MagnetizationSquared")
        .or_else(|| results.get("M2"))
        .or_else(|| results.get("MagSquared"))
        .expect("MagnetizationSquared observable");
    assert!(
        (m2.mean - exact_m2).abs() < 3.0 * m2.stderr.max(0.05),
        "Hybrid ⟨m²⟩={:.4} ± {:.4}, exact={:.4}",
        m2.mean,
        m2.stderr,
        exact_m2
    );
}

// Event-chain, NPT, μVT, WL, Binder cumulant physical validation
// tests are in long_convergence.rs (run via --ignored).
// - NPT: V increases when P decreases
// - μVT: N increases with μ
// - WL: DOS matches exact 4×4 enumeration
// - Binder: U4 at Tc matches universal value
// - Kawasaki: energy decreases on cooling

// P2.3/P2.4: NPT and μVT long stochastic tests live in long_convergence.rs
// alongside the existing stochastic physics suite (WL 4×4, Kawasaki, etc.).
// They are run via: cargo test --test suite -- long_convergence --ignored
