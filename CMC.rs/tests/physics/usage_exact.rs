//! Physics absolute-correctness tests for new code paths.
//!
//! These tests verify deterministic properties: exact energy, detailed
//! balance, conservation laws, and analytical reference comparisons.

use super::common::{assert_close, direct_ising_energy, exact_ising_moments};
use cmc_rs::*;
use rand::RngExt;

// ── 3D Ising: exact energy on a small cube ────────────────────────────────

#[test]
fn three_d_cubic_lattice_has_correct_edge_count() {
    use cmc_rs::BondType;
    let lattice = build_hypercubic(
        &[3, 3, 3],
        &[BondType::SquareX, BondType::SquareY, BondType::SquareZ],
        true,
    );
    // 3x3x3 = 27 sites, 3 edges per site per axis with PBC = 27*3 = 81 edges
    assert_eq!(lattice.n_sites, 27);
    assert_eq!(lattice.n_edges(), 81);
}

#[test]
fn three_d_ising_energy_matches_exact_edge_sum() {
    use cmc_rs::{BondType, Hamiltonian};
    let lattice = build_hypercubic(
        &[2, 2, 2],
        &[BondType::SquareX, BondType::SquareY, BondType::SquareZ],
        true,
    );
    let model = IsingModel::new(1.0);

    // 2x2x2 PBC: each of 8 sites has 3 forward neighbors → 24 physical edges
    assert_eq!(lattice.n_edges(), 24);

    // All spins up → each edge contributes -1.0
    let spins = vec![1.0; 8];
    let energy = model.compute_total_energy(&spins, &lattice, 0.0);
    assert_close(energy, -24.0, 1e-15);

    // Checkered: alternating spins → each edge connects opposite spins → +1.0
    let spins = vec![-1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0];
    let energy = model.compute_total_energy(&spins, &lattice, 0.0);
    assert_close(energy, 24.0, 1e-15);
}

#[test]
fn three_d_ising_every_flip_delta_matches_recomputation() {
    use cmc_rs::{BondType, Hamiltonian};
    let lattice = build_hypercubic(
        &[2, 2, 3],
        &[BondType::SquareX, BondType::SquareY, BondType::SquareZ],
        true,
    );
    let model = IsingModel::new(0.7);

    let spins = vec![
        1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0,
    ];
    let base_energy = model.compute_total_energy(&spins, &lattice, 0.0);

    for site in 0..12 {
        let mut flipped = spins.clone();
        flipped[site] = -flipped[site];
        let new_energy = model.compute_total_energy(&flipped, &lattice, 0.0);
        let delta = model.delta_energy(&spins, &lattice, site, &[flipped[site]]);
        assert_close(new_energy - base_energy, delta, 1e-14);
    }
}

// ── Continuous-spin cluster detailed balance ──────────────────────────────

#[test]
fn xy_cluster_activation_matches_fortuin_kasteleyn() {
    use cmc_rs::{Bond, BondType, ClusterAuxiliary, ClusterModel};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    let model = XYModel::new(1.0);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    // For XY model, the FK probability is p = 1 - exp(-2*beta*J*s_i·s_j)
    // The auxiliary determines the reflection axis.
    let seed = [1.0, 0.0];
    let aux = model.wolff_auxiliary(&seed, &mut rng);

    // For XY, the auxiliary should be a ContinuousTarget (a unit vector)
    match &aux {
        ClusterAuxiliary::Reflection(refl) => {
            let norm: f64 = refl.iter().map(|x| x * x).sum::<f64>().sqrt();
            assert_close(norm, 1.0, 1e-15);
        }
        _ => panic!("XY model should produce Reflection auxiliary"),
    }

    // FK bond probability for aligned spins must be in (0, 1)
    let bond = Bond::new(0, 1, BondType::Generic, 1.0);
    let beta = 1.0;
    let p = model.cluster_bond_probability(&[1.0, 0.0], &[1.0, 0.0], &bond, &aux, beta);
    assert!(
        p > 0.0 && p < 1.0,
        "FK bond probability for aligned XY spins should be in (0,1), got {p}"
    );
    // Verify against the FK formula: p = 1 - exp(-2*beta*J*(s_i·r)*(s_j·r))
    if let ClusterAuxiliary::Reflection(refl) = &aux {
        let proj: f64 = [1.0, 0.0].iter().zip(refl.iter()).map(|(a, b)| a * b).sum();
        let expected = 1.0 - (-2.0 * beta * 1.0 * proj * proj).exp();
        assert_close(p, expected, 1e-14);
    }
    // Anti-aligned spins must give zero bond probability
    let p_anti = model.cluster_bond_probability(&[1.0, 0.0], &[-1.0, 0.0], &bond, &aux, beta);
    assert_eq!(p_anti, 0.0);
}

#[test]
fn heisenberg_cluster_activation_uses_continuous_reflection() {
    use cmc_rs::{Bond, BondType, ClusterAuxiliary, ClusterModel};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    let model = HeisenbergModel::new(1.0);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);

    let seed = [0.0, 0.6, 0.8];
    let aux = model.wolff_auxiliary(&seed, &mut rng);

    match &aux {
        ClusterAuxiliary::Reflection(refl) => {
            assert_eq!(refl.len(), 3);
            let norm: f64 = refl.iter().map(|x| x * x).sum::<f64>().sqrt();
            assert_close(norm, 1.0, 1e-15);
        }
        _ => panic!("Heisenberg model should produce Reflection auxiliary"),
    }

    // FK bond probability for aligned spins must be in (0, 1)
    let bond = Bond::new(0, 1, BondType::Generic, 1.0);
    let beta = 1.0;
    let p = model.cluster_bond_probability(&[0.0, 0.6, 0.8], &[0.0, 0.6, 0.8], &bond, &aux, beta);
    assert!(
        p > 0.0 && p < 1.0,
        "FK bond probability for aligned Heisenberg spins should be in (0,1), got {p}"
    );
    // Verify against the FK formula: p = 1 - exp(-2*beta*J*(s_i·r)*(s_j·r))
    if let ClusterAuxiliary::Reflection(refl) = &aux {
        let proj: f64 = [0.0, 0.6, 0.8]
            .iter()
            .zip(refl.iter())
            .map(|(a, b)| a * b)
            .sum();
        let expected = 1.0 - (-2.0 * beta * 1.0 * proj * proj).exp();
        assert_close(p, expected, 1e-14);
    }
    // Anti-aligned spins must give zero bond probability
    let p_anti =
        model.cluster_bond_probability(&[0.0, 0.6, 0.8], &[0.0, -0.6, -0.8], &bond, &aux, beta);
    assert_eq!(p_anti, 0.0);
}

// ── Over-relaxation energy preservation for XY model ──────────────────────

#[test]
fn xy_over_relaxation_preserves_energy_exactly() {
    use cmc_rs::{Algorithm, MicrocanonicalCore, SimulationPhase};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    let lattice = build_square(4, 4, true);
    let model = XYModel::new(1.0);
    let mut system = System::new(lattice, 2, 0.0, 0.5);

    // Initialize with random unit vectors
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);
    for site in 0..system.n_sites() {
        let theta = rng.random::<f64>() * std::f64::consts::TAU;
        system
            .spin_at_mut(site, 2)
            .copy_from_slice(&[theta.cos(), theta.sin()]);
    }
    system.recompute_energy(&model);
    let initial_energy = system.energy;

    let mut kernel = MicrocanonicalCore::new();
    for _ in 0..50 {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
    }

    // Over-relaxation must preserve energy
    assert_close(system.energy, initial_energy, 1e-10);
    assert_close(system.energy_error(&model), 0.0, 1e-10);

    // All spins must remain unit vectors
    for spin in system.spins.as_chunks::<2>().0 {
        let norm: f64 = spin.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert_close(norm, 1.0, 1e-10);
    }
}

// ── Continuous heat bath: energy stays in valid range ─────────────────────

#[test]
fn continuous_heat_bath_heisenberg_energy_is_physical() {
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    let lattice = build_square(6, 6, true);
    let model = HeisenbergModel::new(1.0);
    let beta = 1.0;
    let mut system = System::new(lattice, 3, 0.0, beta);

    // Random initial state
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    for site in 0..system.n_sites() {
        let theta = rng.random::<f64>() * std::f64::consts::TAU;
        let z = rng.random::<f64>() * 2.0 - 1.0;
        let r = (1.0 - z * z).sqrt();
        system
            .spin_at_mut(site, 3)
            .copy_from_slice(&[r * theta.cos(), r * theta.sin(), z]);
    }
    system.recompute_energy(&model);

    let mut kernel = ContinuousHeatBathCore::new();
    for _ in 0..100 {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
    }

    // Energy must be in physical range: [-n_edges * J, +n_edges * J]
    let n_edges = system.lattice.n_edges() as f64;
    assert!(
        system.energy >= -n_edges * 1.0 - 1e-10,
        "energy {} below physical minimum {}",
        system.energy,
        -n_edges
    );
    assert!(
        system.energy <= n_edges * 1.0 + 1e-10,
        "energy {} above physical maximum {}",
        system.energy,
        n_edges
    );
    // At β=1.0 with ferromagnetic J=1, thermal correlations should make E < 0
    assert!(
        system.energy < 0.0,
        "ferromagnetic Heisenberg at β=1.0 should have negative energy, got {}",
        system.energy
    );
    assert_close(system.energy_error(&model), 0.0, 1e-10);

    // All spins must be unit vectors
    for spin in system.spins.as_chunks::<3>().0 {
        let norm: f64 = spin.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert_close(norm, 1.0, 1e-10);
    }
}

// ── Wang-Landau on 2D Ising: energy levels match exact ────────────────────

#[test]
fn wang_landau_on_2d_square_finds_all_energy_levels() {
    let lattice = build_square(2, 2, true);
    let model = IsingModel::new(1.0);
    let exact = enumerate_ising_density_of_states(&lattice, &model).unwrap();

    // 2x2 PBC Ising: 8 edges (each site has 2 neighbors per axis), J=1
    // Energy levels: -8 (all aligned, deg 2), 0 (deg 12), +8 (deg 2)
    assert_eq!(exact.energies().len(), 3);
    assert_close(exact.energies()[0], -8.0, 1e-14);
    assert_close(exact.energies()[1], 0.0, 1e-14);
    assert_close(exact.energies()[2], 8.0, 1e-14);

    assert_eq!(exact.degeneracies()[0], 2);
    assert_eq!(exact.degeneracies()[1], 12);
    assert_eq!(exact.degeneracies()[2], 2);
    assert_eq!(exact.states(), 16);
}

// ── Wang-Landau reweighting matches exact 2D moments ──────────────────────

#[test]
fn two_d_ising_dos_reweighting_matches_exact_thermodynamics() {
    let lattice = build_square(2, 2, true);
    let model = IsingModel::new(1.0);
    let dos = enumerate_ising_density_of_states(&lattice, &model).unwrap();
    let axis = dos.axis().unwrap();
    let log_density = dos.log_density().unwrap();

    for beta in [0.1, 0.5, 1.0, 2.0] {
        let rw = canonical_reweight(&axis, &log_density, beta).unwrap();

        // Exact from enumeration
        let (_, exact_e, exact_e2, _) = exact_ising_moments(&lattice, 1.0, beta);

        assert_close(rw.mean_energy(), exact_e, 1e-12);
        assert_close(rw.mean_energy_squared(), exact_e2, 1e-12);
    }
}

// ── Honeycomb lattice: coordination number is 3 ───────────────────────────

#[test]
fn honeycomb_lattice_has_coordination_number_three() {
    let lattice = build_honeycomb(4, 4);
    // Honeycomb is 3-regular: every site has exactly 3 neighbors
    for site in 0..lattice.n_sites {
        let degree = lattice.offsets[site + 1] - lattice.offsets[site];
        assert_eq!(degree, 3, "site {site} has degree {degree}, expected 3");
    }
    // 4x4 = 16 sites, 3 edges per site / 2 = 24 edges
    assert_eq!(lattice.n_edges(), 24);
}

// ── 3D Ising Metropolis matches exact 2x2x2 ───────────────────────────────

#[test]
fn metropolis_3d_ising_matches_exact_2x2x2_energy() {
    use cmc_rs::BondType;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    let lattice = build_hypercubic(
        &[2, 2, 2],
        &[BondType::SquareX, BondType::SquareY, BondType::SquareZ],
        true,
    );
    let coupling = 1.0;
    let beta = 0.8;

    // Exact energy from enumeration (256 states for 8 sites)
    let mut z = 0.0;
    let mut e_sum = 0.0;
    for mask in 0..1u32 << lattice.n_sites {
        let spins: Vec<f64> = (0..lattice.n_sites)
            .map(|s| if mask & (1 << s) == 0 { -1.0 } else { 1.0 })
            .collect();
        let energy = direct_ising_energy(&spins, &lattice, coupling);
        let weight = (-beta * energy).exp();
        z += weight;
        e_sum += weight * energy;
    }
    let exact_mean_e = e_sum / z;

    // MC run
    let model = IsingModel::new(coupling);
    let mut system = System::new(lattice, 1, 0.0, beta);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    // Random initial state
    for site in 0..system.n_sites() {
        system.spins[site] = if rng.random::<bool>() { 1.0 } else { -1.0 };
    }
    system.recompute_energy(&model);

    // Thermalization
    let mut metro = MetropolisCore::new();
    for _ in 0..5000 {
        metro.sweep_with_phase(
            &mut system,
            &model,
            &mut rng,
            SimulationPhase::Thermalization,
        );
    }

    // Measurement
    let mut e_samples = Vec::new();
    for _ in 0..10_000 {
        metro.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
        e_samples.push(system.energy);
    }

    let mc_mean_e: f64 = e_samples.iter().sum::<f64>() / e_samples.len() as f64;
    let mc_std_e: f64 = {
        let var = e_samples
            .iter()
            .map(|e| (e - mc_mean_e).powi(2))
            .sum::<f64>()
            / e_samples.len() as f64;
        var.sqrt() / (e_samples.len() as f64).sqrt()
    };

    // Should agree within 4 sigma
    let diff = (mc_mean_e - exact_mean_e).abs();
    assert!(
        diff < 4.0 * mc_std_e,
        "3D Ising mean E mismatch: MC={mc_mean_e:.4}±{mc_std_e:.4}, exact={exact_mean_e:.4}"
    );
}

// ── Parallel tempering: log_weight_ratio calls the real implementation ────

#[test]
fn parallel_tempering_log_weight_ratio_matches_implementation() {
    use carlo_rs::ParallelTemperingCompatible;

    let lattice = build_square(4, 4, true);
    let model = IsingModel::new(1.0);
    let beta_old = 1.0;
    let mut system = System::new(lattice, 1, 0.0, beta_old);
    let mut rng = rand::rng();
    for site in 0..system.n_sites() {
        system.spins[site] = if rng.random::<bool>() { 1.0 } else { -1.0 };
    }
    system.recompute_energy(&model);

    let mc = ClassicalMC {
        system,
        model,
        algorithm: MetropolisCore::new(),
        observables: DefaultObservableSet::default(),
    };

    let beta_new = 2.0;
    let energy = mc.system.energy;
    let expected = (beta_old - beta_new) * energy;
    let actual = mc.log_weight_ratio("beta", beta_new);
    assert_close(actual, expected, 1e-15);
}

// ── Continuous-spin Wolff: energy matches analytical solution ─────────────
//
// For a 2-site O(3) chain with PBC (1 edge), H = -J cos θ where θ is the
// angle between the two spins. The exact mean energy is:
//   <E> = -J * L(βJ)  where L(x) = coth(x) - 1/x  (Langevin function)

#[test]
fn heisenberg_wolff_energy_matches_langevin_exact() {
    use cmc_rs::{Algorithm, Hamiltonian, WolffCore};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    let lattice = build_chain(2, true); // 1 edge
                                        // Diagnostic: check edge count and aligned energy
    let n_edges = lattice.n_edges();
    let j = 1.0;
    let beta = 1.5;
    let model = HeisenbergModel::new(j);

    // Aligned state energy = -n_edges * J
    let aligned = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
    let aligned_e = model.compute_total_energy(&aligned, &lattice, beta);
    assert!(
        (aligned_e - (-(n_edges as f64) * j)).abs() < 1e-14,
        "aligned E={aligned_e}, expected={}",
        -(n_edges as f64) * j
    );

    let mut system = System::new(lattice, 3, 0.0, beta);

    // Random initial state
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    for site in 0..system.n_sites() {
        let theta = rng.random::<f64>() * std::f64::consts::TAU;
        let z = rng.random::<f64>() * 2.0 - 1.0;
        let r = (1.0 - z * z).sqrt();
        system
            .spin_at_mut(site, 3)
            .copy_from_slice(&[r * theta.cos(), r * theta.sin(), z]);
    }
    system.recompute_energy(&model);

    let mut wolff = WolffCore::new();

    // Thermalization
    for _ in 0..5_000 {
        wolff.sweep(&mut system, &model, &mut rng);
    }

    // Measurement
    let n_measure = 50_000;
    let mut e_sum = 0.0;
    for _ in 0..n_measure {
        wolff.sweep(&mut system, &model, &mut rng);
        e_sum += system.energy;
    }
    let mc_mean_e = e_sum / n_measure as f64;

    // Exact: H = -α cos θ with α = n_edges * J.
    // For O(3): <E> = -α * L(βα) where L(x) = coth(x) - 1/x
    let alpha = n_edges as f64 * j;
    let x = beta * alpha;
    let coth = (x.exp() + (-x).exp()) / (x.exp() - (-x).exp());
    let langevin = coth - 1.0 / x;
    let exact_e = -alpha * langevin;

    assert!(
        (mc_mean_e - exact_e).abs() < 0.01,
        "Heisenberg Wolff: MC={mc_mean_e:.6}, exact(Langevin,n_edges={n_edges})={exact_e:.6}"
    );
}

// ── XY Wolff: energy matches Bessel ratio exact ───────────────────────────
//
// For a 2-site O(2) chain with PBC (1 edge), H = -J cos θ.
// <E> = -J * I₁(βJ) / I₀(βJ)  where I₀, I₁ are modified Bessel functions.

#[test]
fn xy_wolff_energy_matches_bessel_ratio_exact() {
    use cmc_rs::{Algorithm, Hamiltonian, WolffCore};
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    let lattice = build_chain(2, true);
    let n_edges = lattice.n_edges();
    let j = 1.0;
    let beta = 2.0;
    let model = XYModel::new(j);

    // Aligned state energy = -n_edges * J
    let aligned = vec![1.0, 0.0, 1.0, 0.0];
    let aligned_e = model.compute_total_energy(&aligned, &lattice, beta);
    assert!(
        (aligned_e - (-(n_edges as f64) * j)).abs() < 1e-14,
        "aligned E={aligned_e}, expected={}",
        -(n_edges as f64) * j
    );

    let mut system = System::new(lattice, 2, 0.0, beta);

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
    for site in 0..system.n_sites() {
        let theta = rng.random::<f64>() * std::f64::consts::TAU;
        system
            .spin_at_mut(site, 2)
            .copy_from_slice(&[theta.cos(), theta.sin()]);
    }
    system.recompute_energy(&model);

    let mut wolff = WolffCore::new();
    for _ in 0..5_000 {
        wolff.sweep(&mut system, &model, &mut rng);
    }

    let n_measure = 50_000;
    let mut e_sum = 0.0;
    for _ in 0..n_measure {
        wolff.sweep(&mut system, &model, &mut rng);
        e_sum += system.energy;
    }
    let mc_mean_e = e_sum / n_measure as f64;

    // Exact: H = -α cos θ with α = n_edges * J.
    // For O(2): <E> = -α * I₁(βα) / I₀(βα)
    let alpha = n_edges as f64 * j;
    let x = beta * alpha;
    let i0 = modified_bessel_i0(x);
    let i1 = modified_bessel_i1(x);
    let exact_e = -alpha * i1 / i0;

    assert!(
        (mc_mean_e - exact_e).abs() < 0.01,
        "XY Wolff: MC={mc_mean_e:.6}, exact(Bessel,n_edges={n_edges})={exact_e:.6}"
    );
}

/// Modified Bessel function I₀(x) via series: Σ (x/2)²ᵏ / (k!)².
fn modified_bessel_i0(x: f64) -> f64 {
    let half_x_sq = (x / 2.0).powi(2);
    let mut sum = 1.0;
    let mut term = 1.0;
    for k in 1..100 {
        term *= half_x_sq / (k as f64 * k as f64);
        sum += term;
        if term < 1e-16 * sum {
            break;
        }
    }
    sum
}

/// Modified Bessel function I₁(x) via series: Σ (x/2)²ᵏ⁺¹ / (k!(k+1)!).
fn modified_bessel_i1(x: f64) -> f64 {
    let half_x = x / 2.0;
    let half_x_sq = half_x * half_x;
    let mut term = half_x; // k=0 term
    let mut sum = term;
    for k in 1..100 {
        term *= half_x_sq / (k as f64 * (k as f64 + 1.0));
        sum += term;
        if term < 1e-16 * sum {
            break;
        }
    }
    sum
}
