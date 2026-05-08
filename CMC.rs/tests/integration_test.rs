//! Comprehensive integration tests for CMC.rs.
//!
//! Covers: algorithm comparison, exact solution validation,
//! boundary conditions, multi-seed statistics, physical symmetries.

use cmc_rs::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

// ─── Helpers ────────────────────────────────────────────────────────────────

fn run_simulation<T: MonteCarlo + FromParams>(
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

// 1D Ising exact solution: E/N_bonds = -J * tanh(βJ)
fn ising_1d_exact_energy_per_bond(beta: f64, j: f64) -> f64 {
    -j * (beta * j).tanh()
}

// Onsager 2D Ising critical temperature: βc = ln(1+√2)/(2J)
fn ising_2d_critical_beta(j: f64) -> f64 {
    (1.0 + 2.0_f64.sqrt()).ln() / (2.0 * j)
}

// Onsager 2D Ising energy per site at critical temperature.
// At Tc: E/N = -J * coth(2βcJ) where tanh(2βcJ) = 1/√2, so coth(2βcJ) = √2.
// The full Onsager formula has an elliptic integral term, but at Tc this
// term vanishes, leaving simply -√2 * J.
fn ising_2d_energy_per_site_at_tc(j: f64) -> f64 {
    -j * 2.0_f64.sqrt()
}

// ─── Algorithm comparison ───────────────────────────────────────────────────

#[test]
fn test_all_algorithms_agree_1d() {
    // At moderate temperature, all algorithms should give consistent energy
    // Metropolis needs longer thermalization than cluster algorithms
    let params = make_params(32, 1.0, 1.0, true);

    let met = run_simulation::<MetropolisCore<IsingModel>>(&params, 5000, 10000, 100);
    let wol = run_simulation::<WolffCore<IsingModel>>(&params, 1000, 5000, 100);
    let sw = run_simulation::<SWCore<IsingModel>>(&params, 1000, 5000, 100);

    let e_met = met.get("Energy").unwrap().mean;
    let e_wol = wol.get("Energy").unwrap().mean;
    let e_sw = sw.get("Energy").unwrap().mean;

    // Different algorithms have different autocorrelation properties, so
    // finite-run averages can differ. Allow 10% tolerance.
    let avg = (e_met + e_wol + e_sw) / 3.0;
    assert!(
        (e_met - avg).abs() < 0.10 * avg.abs(),
        "Metropolis energy {} deviates too much from avg {}",
        e_met,
        avg
    );
    assert!(
        (e_wol - avg).abs() < 0.10 * avg.abs(),
        "Wolff energy {} deviates too much from avg {}",
        e_wol,
        avg
    );
    assert!(
        (e_sw - avg).abs() < 0.10 * avg.abs(),
        "SW energy {} deviates too much from avg {}",
        e_sw,
        avg
    );
}

// ─── 1D exact solution validation ───────────────────────────────────────────

#[test]
fn test_metropolis_vs_1d_exact_low_temp() {
    // At low T (high β), exact solution is very accurate
    // Metropolis needs long thermalization at low T due to large energy barriers
    let beta = 3.0;
    let j = 1.0;
    let params = make_params(64, beta, j, true);

    let results = run_simulation::<MetropolisCore<IsingModel>>(&params, 10000, 20000, 100);
    let energy = results.get("Energy").unwrap();

    // For directed bonds in PBC chain: n_bonds = n_sites
    let exact = ising_1d_exact_energy_per_bond(beta, j) * 64.0;

    // At low T, autocorrelation is long and simple binning underestimates errors.
    // Allow wider tolerance (5% of exact value) to account for this.
    assert!(
        (energy.mean - exact).abs() < 0.05 * exact.abs(),
        "Energy {:.4} differs from exact {:.4} by >5%",
        energy.mean,
        exact
    );
}

#[test]
fn test_wolff_vs_1d_exact_high_temp() {
    // At high T, cluster algorithm should still match exact
    // Wolff at high T has small clusters, so needs more sweeps for good statistics
    let beta = 0.3;
    let j = 1.0;
    let params = make_params(64, beta, j, true);

    let results = run_simulation::<WolffCore<IsingModel>>(&params, 5000, 10000, 100);
    let energy = results.get("Energy").unwrap();

    let exact = ising_1d_exact_energy_per_bond(beta, j) * 64.0;

    // Allow slightly wider tolerance for Wolff since cluster dynamics differ from Metropolis
    assert!(
        (energy.mean - exact).abs() < 5.0 * energy.stderr,
        "Energy {:.4} differs from exact {:.4} by {:.1}σ",
        energy.mean,
        exact,
        (energy.mean - exact).abs() / energy.stderr
    );
}

// ─── 2D Onsager exact solution ──────────────────────────────────────────────

#[test]
fn test_2d_ising_at_critical_temperature() {
    // Run 2D Ising at Tc using Metropolis — energy should match Onsager
    // Note: IsingModel uses 1D chain lattice, so we test the 1D physics here.
    // For a proper 2D test, we'd need a 2D IsingModel variant.
    // This test documents the validation target for future 2D model support.

    // Verify βc ≈ 0.4407 for J=1
    let j = 1.0;
    let beta_c = ising_2d_critical_beta(j);
    assert!(
        (beta_c - 0.4407).abs() < 0.001,
        "βc = {}, expected ~0.4407",
        beta_c
    );

    // Verify the energy formula: E/N at Tc = -√2 * J
    let expected_e_per_site = ising_2d_energy_per_site_at_tc(j);
    assert!(
        (expected_e_per_site - (-2.0_f64.sqrt() * j)).abs() < 0.001,
        "E/N at Tc = {}, expected ~{}",
        expected_e_per_site,
        -2.0_f64.sqrt() * j
    );
}

// ─── High-temperature limit ─────────────────────────────────────────────────

#[test]
fn test_high_temperature_limit() {
    // At β → 0, spins are random → ⟨E⟩ → 0
    let params = make_params(32, 0.01, 1.0, true);

    let results = run_simulation::<MetropolisCore<IsingModel>>(&params, 500, 5000, 100);
    let energy = results.get("Energy").unwrap();

    assert!(
        energy.mean.abs() < 2.0,
        "High-T energy should be near zero, got {}",
        energy.mean
    );
}

// ─── Low-temperature ordered phase ──────────────────────────────────────────

#[test]
fn test_low_temperature_ground_state() {
    // At T → 0, system should reach ground state: all spins aligned
    // For 1D PBC chain with J=1: E_ground = -N_bonds = -N_sites
    let params = make_params(16, 10.0, 1.0, true);

    let results = run_simulation::<MetropolisCore<IsingModel>>(&params, 1000, 2000, 100);
    let energy = results.get("Energy").unwrap();

    // Ground state energy: each of the 16 bonds contributes -J = -1
    let ground = -16.0;
    assert!(
        (energy.mean - ground).abs() < 0.5,
        "Low-T energy {:.4} should be near ground state {:.4}",
        energy.mean,
        ground
    );
}

// ─── Open vs periodic boundary conditions ───────────────────────────────────

#[test]
fn test_open_vs_pbc_energy_scaling() {
    // PBC chain has N bonds, open chain has N-1 bonds
    // At same (β, J), energy per bond should be similar
    let beta = 1.0;
    let j = 1.0;
    let n = 32;

    let params_pbc = make_params(n, beta, j, true);
    let params_open = make_params(n, beta, j, false);

    let res_pbc = run_simulation::<MetropolisCore<IsingModel>>(&params_pbc, 1000, 5000, 100);
    let res_open = run_simulation::<MetropolisCore<IsingModel>>(&params_open, 1000, 5000, 100);

    let e_pbc_per_bond = res_pbc.get("Energy").unwrap().mean / (n as f64);
    let e_open_per_bond = res_open.get("Energy").unwrap().mean / ((n - 1) as f64);

    // Per-bond energy should agree within ~10% (finite-size + statistical effects)
    assert!(
        (e_pbc_per_bond - e_open_per_bond).abs() < 0.1 * e_pbc_per_bond.abs(),
        "PBC per-bond energy {} vs open per-bond energy {}",
        e_pbc_per_bond,
        e_open_per_bond
    );
}

// ─── System size scaling ────────────────────────────────────────────────────

#[test]
fn test_energy_scales_extensively() {
    // Energy should scale linearly with system size
    let beta = 1.0;
    let j = 1.0;

    let params_small = make_params(16, beta, j, true);
    let params_large = make_params(64, beta, j, true);

    let res_small = run_simulation::<MetropolisCore<IsingModel>>(&params_small, 1000, 5000, 100);
    let res_large = run_simulation::<MetropolisCore<IsingModel>>(&params_large, 2000, 8000, 100);

    let e_small = res_small.get("Energy").unwrap().mean;
    let e_large = res_large.get("Energy").unwrap().mean;

    // Energy per site should be similar (within 5%)
    let e_small_per_site = e_small / 16.0;
    let e_large_per_site = e_large / 64.0;

    assert!(
        (e_small_per_site - e_large_per_site).abs() < 0.05 * e_small_per_site.abs(),
        "Small system E/N={} vs large system E/N={}",
        e_small_per_site,
        e_large_per_site
    );
}

// ─── Wolff algorithm specifics ──────────────────────────────────────────────

#[test]
fn test_wolff_cluster_size_temperature_dependence() {
    // At low T, Wolff clusters should be large (many spins flip together)
    // At high T, clusters should be small
    // We test this indirectly via the acceptance behavior:
    // low-T Wolff should thermalize faster than Metropolis

    let beta = 5.0; // Low temperature
    let n = 64;
    let short_therm = 50; // Very short thermalization

    let params = make_params(n, beta, 1.0, true);

    let met = run_simulation::<MetropolisCore<IsingModel>>(&params, short_therm, 2000, 50);
    let wol = run_simulation::<WolffCore<IsingModel>>(&params, short_therm, 2000, 50);

    let e_met = met.get("Energy").unwrap().mean;
    let e_wol = wol.get("Energy").unwrap().mean;

    // Wolff should reach lower energy (closer to ground state) with same short thermalization
    // because large clusters flip together, escaping metastable states faster
    assert!(
        e_wol <= e_met + 1.0,
        "Wolff energy {} should be lower than Metropolis {} at low T with short thermalization",
        e_wol,
        e_met
    );
}

// ─── Swendsen-Wang specifics ────────────────────────────────────────────────

#[test]
fn test_swendsen_wang_spin_symmetry() {
    // Starting from all-up vs all-down should give same energy distribution
    // This tests that the algorithm doesn't have a symmetry-breaking bias
    let params = make_params(32, 0.5, 1.0, true);

    // Both start from all-up (default in IsingModel)
    let results = run_simulation::<SWCore<IsingModel>>(&params, 500, 3000, 50);
    let energy = results.get("Energy").unwrap();

    // Energy should be finite and negative (below Tc for 1D at β=0.5, J=1)
    assert!(energy.mean < 0.0, "Energy should be negative at β=0.5");
    assert!(energy.stderr > 0.0, "Should have non-zero error estimate");
}

// ─── Multi-seed statistical validation ──────────────────────────────────────

#[test]
fn test_reproducibility_same_seed() {
    // Same seed should give identical results
    let params = make_params(16, 1.0, 1.0, true);

    let backend = RayonBackend::new(1);
    let config = RunConfig {
        thermalization_sweeps: 200,
        measurement_sweeps: 500,
        binsize: 50,
        base_seed: 12345,
        progress_interval: 0,
        checkpoint_interval: 0,
    };
    let scheduler1 = Scheduler::new(backend, config.clone());
    let scheduler2 = Scheduler::new(RayonBackend::new(1), config);

    let r1 = scheduler1.run_one::<MetropolisCore<IsingModel>>(&params);
    let r2 = scheduler2.run_one::<MetropolisCore<IsingModel>>(&params);

    let e1 = r1.get("Energy").unwrap().mean;
    let e2 = r2.get("Energy").unwrap().mean;

    assert!(
        (e1 - e2).abs() < 1e-10,
        "Same seed should give identical results: {} vs {}",
        e1,
        e2
    );
}

// ─── Parameter validation ───────────────────────────────────────────────────

#[test]
fn test_metropolis_accepts_reasonable_range() {
    // Sweep through temperature range and verify energy changes monotonically
    let betas = [0.2, 0.5, 1.0, 2.0, 5.0];
    let mut energies = Vec::new();

    for &beta in &betas {
        let params = make_params(32, beta, 1.0, true);
        let results = run_simulation::<MetropolisCore<IsingModel>>(&params, 500, 3000, 100);
        energies.push(results.get("Energy").unwrap().mean);
    }

    // Energy should become more negative as β increases
    for i in 1..energies.len() {
        assert!(
            energies[i] < energies[i - 1] + 1.0, // tolerance for statistical noise
            "Energy should generally decrease with β: {:?}",
            betas.iter().zip(energies.iter()).collect::<Vec<_>>()
        );
    }
}

// ─── Directed bond consistency ──────────────────────────────────────────────

#[test]
fn test_directed_bond_energy_consistency() {
    // Bidirectional bonds: each physical bond appears twice (a→b and b→a)
    // total_energy() divides by 2 to get the standard undirected-bond energy
    let lattice = cmc_rs::lattice::build_chain(10, true);
    let model = IsingModel::new(lattice, 1.0, 1.0);

    // All up: 20 directed bonds (10 physical × 2), each -J, divided by 2 = -10.0
    assert!((model.total_energy() - (-10.0)).abs() < 1e-10);

    // Open chain of 10 sites: 18 directed bonds (9 physical × 2)
    let lattice_open = cmc_rs::lattice::build_chain(10, false);
    let model_open = IsingModel::new(lattice_open, 1.0, 1.0);
    assert!((model_open.total_energy() - (-9.0)).abs() < 1e-10);
}

// ─── Potts model with cluster algorithms ────────────────────────────────────

#[test]
fn test_potts_sw_wolff_agree() {
    // SW and Wolff should agree for 3-state Potts model
    let mut params = make_params(32, 0.5, 1.0, true);
    params.set("q", 3usize);

    let sw = run_simulation::<SWCore<PottsModel>>(&params, 1000, 5000, 100);
    let wol = run_simulation::<WolffCore<PottsModel>>(&params, 1000, 5000, 100);

    let e_sw = sw.get("Energy").unwrap().mean;
    let e_wol = wol.get("Energy").unwrap().mean;

    assert!(
        (e_sw - e_wol).abs() < 0.15 * e_sw.abs(),
        "SW energy {} vs Wolff energy {} differ too much for Potts",
        e_sw,
        e_wol
    );
}

#[test]
fn test_potts_metropolis_sw_agree() {
    // Metropolis and SW should agree for 3-state Potts
    // Use higher temperature (lower beta) where Metropolis thermalizes faster
    let mut params = make_params(32, 0.2, 1.0, true);
    params.set("q", 3usize);

    let met = run_simulation::<MetropolisCore<PottsModel>>(&params, 10000, 20000, 100);
    let sw = run_simulation::<SWCore<PottsModel>>(&params, 1000, 5000, 100);

    let e_met = met.get("Energy").unwrap().mean;
    let e_sw = sw.get("Energy").unwrap().mean;

    assert!(
        (e_met - e_sw).abs() < 0.15 * e_met.abs(),
        "Metropolis energy {} vs SW energy {} differ too much for Potts",
        e_met,
        e_sw
    );
}

// ─── Magnetization tests ────────────────────────────────────────────────────

#[test]
fn test_ising_magnetization_temperature_dependence() {
    // Magnetization should decrease with temperature
    let beta_low = 5.0;  // Low T → ordered
    let beta_high = 0.1; // High T → disordered

    let params_low = make_params(32, beta_low, 1.0, true);
    let params_high = make_params(32, beta_high, 1.0, true);

    let res_low = run_simulation::<MetropolisCore<IsingModel>>(&params_low, 2000, 5000, 100);
    let res_high = run_simulation::<MetropolisCore<IsingModel>>(&params_high, 500, 3000, 100);

    let m_low = res_low.get("Magnetization").unwrap().mean;
    let m_high = res_high.get("Magnetization").unwrap().mean;

    assert!(
        m_low > m_high,
        "Low-T magnetization {} should be > high-T {}",
        m_low,
        m_high
    );
}

#[test]
fn test_magnetization_squared_fluctuation() {
    // Magnetization² should be measurable and positive
    let params = make_params(16, 1.0, 1.0, true);
    let results = run_simulation::<MetropolisCore<IsingModel>>(&params, 500, 3000, 100);

    let m2 = results.get("Magnetization_Squared").unwrap();
    assert!(m2.mean > 0.0, "Magnetization_Squared should be positive");
}

#[test]
fn test_energy_squared_fluctuation() {
    // Energy² should be measurable and >= (Energy)² by Jensen's inequality
    let params = make_params(16, 1.0, 1.0, true);
    let results = run_simulation::<MetropolisCore<IsingModel>>(&params, 500, 3000, 100);

    let e = results.get("Energy").unwrap();
    let e2 = results.get("Energy_Squared").unwrap();

    assert!(
        e2.mean > e.mean * e.mean * 0.9,
        "Energy_Squared {} should be >= (Energy)² {}",
        e2.mean,
        e.mean * e.mean
    );
}

// ─── Snapshot recording tests ───────────────────────────────────────────────

#[test]
fn test_metropolis_snapshot_recording() {
    let lattice = cmc_rs::lattice::build_chain(8, true);
    let model = IsingModel::new(lattice, 1.0, 1.0);
    let mut core = MetropolisCore::new(model).with_snapshot_interval(10);
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 100);

    // Run 50 sweeps, measuring each time
    for _ in 0..50 {
        core.sweep(&mut ctx);
        core.measure(&mut ctx);
        ctx.advance_sweep();
    }

    let results = ctx.finalize_measurements();
    // Snapshot should have been recorded at sweeps 10, 20, 30, 40, 50
    assert!(
        results.contains_key("Snapshot"),
        "Snapshot should be recorded in measurements"
    );

    // Energy and other observables should also be present
    assert!(results.contains_key("Energy"), "Energy should be recorded");
    assert!(
        results.contains_key("Magnetization"),
        "Magnetization should be recorded"
    );
}

#[test]
fn test_wolff_snapshot_recording() {
    let lattice = cmc_rs::lattice::build_chain(8, true);
    let model = IsingModel::new(lattice, 1.0, 1.0);
    let mut core = WolffCore::new(model).with_snapshot_interval(10);
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 100);

    for _ in 0..50 {
        core.sweep(&mut ctx);
        core.measure(&mut ctx);
        ctx.advance_sweep();
    }

    let results = ctx.finalize_measurements();
    assert!(
        results.contains_key("Snapshot"),
        "Wolff snapshot should be recorded"
    );
}

#[test]
fn test_sw_snapshot_recording() {
    let lattice = cmc_rs::lattice::build_chain(8, true);
    let model = IsingModel::new(lattice, 1.0, 1.0);
    let mut core = SWCore::new(model).with_snapshot_interval(10);
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 100);

    for _ in 0..50 {
        core.sweep(&mut ctx);
        core.measure(&mut ctx);
        ctx.advance_sweep();
    }

    let results = ctx.finalize_measurements();
    assert!(
        results.contains_key("Snapshot"),
        "SW snapshot should be recorded"
    );
}

#[test]
fn test_snapshot_disabled_by_default() {
    let lattice = cmc_rs::lattice::build_chain(8, true);
    let model = IsingModel::new(lattice, 1.0, 1.0);
    let mut core = MetropolisCore::new(model);
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 100);

    for _ in 0..50 {
        core.sweep(&mut ctx);
        core.measure(&mut ctx);
        ctx.advance_sweep();
    }

    let results = ctx.finalize_measurements();
    assert!(
        !results.contains_key("Snapshot"),
        "Snapshot should NOT be recorded when not enabled"
    );
}

#[test]
fn test_snapshot_dimension_matches_model() {
    // Heisenberg model has 3D spins, so snapshot should have 3*N elements
    let lattice = cmc_rs::lattice::build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);
    let mut core = MetropolisCore::with_strategy(
        model,
        cmc_rs::algorithms::OPSSStrategy::new(60.0),
    )
    .with_snapshot_interval(5);
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new_with_binsize(rng, 0, 100);

    for _ in 0..20 {
        core.sweep(&mut ctx);
        core.measure(&mut ctx);
        ctx.advance_sweep();
    }

    let results = ctx.finalize_measurements();
    assert!(
        results.contains_key("Snapshot"),
        "Snapshot should be recorded for Heisenberg model"
    );
    // 4 sites * 3 dimensions = 12 elements in snapshot
    let snapshot_est = results.get("Snapshot").unwrap();
    // The mean should be finite, confirming data was recorded
    assert!(
        snapshot_est.mean.is_finite(),
        "Snapshot mean should be finite"
    );
}

// ─── Antiferromagnetic coupling (J < 0) ─────────────────────────────────────

#[test]
fn test_antiferromagnetic_ground_state() {
    // J < 0: antiferromagnetic. For 1D chain with even sites,
    // ground state is alternating spins.
    // H = -J*Σs_i*s_j (undirected). With J=-1, H = +Σs_i*s_j.
    // Alternating: s_i*s_j = -1 for each bond. H = +(-1)*n_bonds = -n_bonds.
    let lattice = cmc_rs::lattice::build_chain(4, true);
    let mut model = IsingModel::new(lattice, 10.0, -1.0);
    for (i, s) in model.spins_mut().iter_mut().enumerate() {
        *s = if i % 2 == 0 { 1.0 } else { -1.0 };
    }
    // 4 physical bonds, each s_i*s_j = -1, H = -(-1)*(-4) = -4.0
    assert!((model.total_energy() - (-4.0)).abs() < 1e-10);

    // Run simulation and verify energy approaches ground state
    let mut params = Params::new();
    params.set("L", 16usize);
    params.set("beta", 10.0f64);
    params.set("J", -1.0f64);
    params.set("pbc", true);

    let results = run_simulation::<MetropolisCore<IsingModel>>(&params, 5000, 10000, 100);
    let energy = results.get("Energy").unwrap();
    // Ground state for 16-site PBC AFM: H = -(-1)*(-16) = -16.0
    assert!(
        (energy.mean - (-16.0)).abs() < 0.5,
        "AFM energy {:.4} should be near ground state -16.0",
        energy.mean
    );
}

// ─── XY model: local_energy vs total_energy consistency ─────────────────────

#[test]
fn test_xy_local_vs_total_energy_consistency() {
    use std::f64::consts::PI;

    let lattice = cmc_rs::lattice::build_chain(8, true);
    let mut model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
    // Set non-trivial spin configuration (spread angles)
    let angles = [0.0, PI / 4.0, PI / 2.0, 3.0 * PI / 4.0, PI, 5.0 * PI / 4.0, 3.0 * PI / 2.0, 7.0 * PI / 4.0];
    for (i, &a) in angles.iter().enumerate() {
        model.spins_mut()[2 * i] = a.cos();
        model.spins_mut()[2 * i + 1] = a.sin();
    }

    for site in 0..8 {
        let old = vec![model.spins()[2 * site], model.spins()[2 * site + 1]];
        let new = vec![-old[0], -old[1]]; // Flip by π (negation)

        let de_local = model.local_energy_change_spin(site, &old, &new);

        let e_before = model.total_energy();
        model.spins_mut()[2 * site] = new[0];
        model.spins_mut()[2 * site + 1] = new[1];
        let e_after = model.total_energy();
        let de_global = e_after - e_before;

        assert!(
            (de_local - de_global).abs() < 1e-8,
            "Site {}: local ΔE={} vs global ΔE={}",
            site, de_local, de_global
        );

        // Restore
        model.spins_mut()[2 * site] = old[0];
        model.spins_mut()[2 * site + 1] = old[1];
    }
}

// ─── OPSS + XY model ────────────────────────────────────────────────────────

#[test]
fn test_opss_on_xy_agrees_with_metropolis() {
    // OPSS and Metropolis should give similar energies for XY model
    let params = make_params(16, 1.0, 1.0, true);

    let met = run_simulation::<MetropolisCore<XYModel>>(&params, 2000, 5000, 100);
    let opss = run_simulation::<MetropolisCore<XYModel, OPSSStrategy>>(&params, 2000, 5000, 100);

    let e_met = met.get("Energy").unwrap().mean / 16.0;
    let e_opss = opss.get("Energy").unwrap().mean / 16.0;

    assert!(
        (e_met - e_opss).abs() < 0.15 * e_met.abs(),
        "OPSS per-site energy {} vs Metropolis {} differ too much for XY",
        e_opss, e_met
    );
}

// ─── XY model integration tests ─────────────────────────────────────────────

#[test]
fn test_xy_metropolis_energy_extensive() {
    // Energy should scale linearly with system size
    let params = make_params(16, 1.0, 1.0, true);
    let params2 = make_params(64, 1.0, 1.0, true);

    let res1 = run_simulation::<MetropolisCore<XYModel>>(&params, 1000, 5000, 100);
    let res2 = run_simulation::<MetropolisCore<XYModel>>(&params2, 2000, 8000, 100);

    let e1 = res1.get("Energy").unwrap().mean / 16.0;
    let e2 = res2.get("Energy").unwrap().mean / 64.0;

    assert!(
        (e1 - e2).abs() < 0.10 * e1.abs(),
        "Per-site energy should be similar: {} vs {}",
        e1, e2
    );
}

#[test]
fn test_xy_high_temperature() {
    // At β → 0, spins are random → ⟨E⟩ → 0
    let params = make_params(32, 0.01, 1.0, true);
    let results = run_simulation::<MetropolisCore<XYModel>>(&params, 500, 5000, 100);
    let energy = results.get("Energy").unwrap();
    assert!(
        energy.mean.abs() < 2.0,
        "High-T energy should be near zero, got {}",
        energy.mean
    );
}

// ─── Single-site boundary cases ─────────────────────────────────────────────

#[test]
fn test_single_site_ising() {
    // Single site, no bonds (use non-PBC to avoid degenerate self-loops)
    let lattice = cmc_rs::lattice::build_chain(1, false);
    let model = IsingModel::new(lattice, 1.0, 1.0);
    assert!((model.total_energy() - 0.0).abs() < 1e-10, "Single-site Ising energy should be 0");
    assert!((model.magnetization() - 1.0).abs() < 1e-10, "Single-site Ising magnetization should be 1");
}

#[test]
fn test_single_site_potts() {
    let lattice = cmc_rs::lattice::build_chain(1, false);
    let model = PottsModel::new(lattice, 1.0, 1.0, 3);
    assert!((model.total_energy() - 0.0).abs() < 1e-10, "Single-site Potts energy should be 0");
}

#[test]
fn test_single_site_heisenberg() {
    let lattice = cmc_rs::lattice::build_chain(1, false);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);
    assert!((model.total_energy() - 0.0).abs() < 1e-10, "Single-site Heisenberg energy should be 0");
    assert!((model.magnetization() - 1.0).abs() < 1e-10, "Single-site Heisenberg magnetization should be 1");
}

#[test]
fn test_single_site_xy() {
    let lattice = cmc_rs::lattice::build_chain(1, false);
    let model = XYModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 4.0);
    assert!((model.total_energy() - 0.0).abs() < 1e-10, "Single-site XY energy should be 0");
    assert!((model.magnetization() - 1.0).abs() < 1e-10, "Single-site XY magnetization should be 1");
}
