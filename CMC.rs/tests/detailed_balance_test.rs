//! Detailed balance tests for classical Monte Carlo transition kernels.
//!
//! For small state spaces (N≤4 Ising chains), we enumerate all states,
//! sample transitions directly, and verify π(x)P(x→y) ≈ π(y)P(y→x).

use cmc_rs::{
    build_chain, Algorithm, Bond, BondType, CanonicalEnsemble, CsrLattice, EnergyPatch,
    Hamiltonian, HeatBathable, IsingModel, ProposalStrategy, ProposedSpin, SiteSpinMove, Spin,
    System, TrialEvaluator, WolffCore,
};
use rand::{RngExt, SeedableRng};

type RngType = rand_xoshiro::Xoshiro256PlusPlus;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn enumerate_ising_states(n: usize) -> Vec<Vec<f64>> {
    (0..(1u32 << n))
        .map(|mask| {
            (0..n)
                .map(|i| if (mask >> i) & 1 == 1 { 1.0 } else { -1.0 })
                .collect()
        })
        .collect()
}

fn boltzmann(spins: &[f64], model: &IsingModel, lattice: &CsrLattice, beta: f64) -> f64 {
    let e = model.compute_total_energy(spins, lattice, 1.0);
    (-beta * e).exp()
}

fn build_system(spins: &[f64], lattice: CsrLattice, beta: f64) -> System {
    let mut system = System::new(lattice, 1, 1.0, beta);
    system.spins.copy_from_slice(spins);
    let model = IsingModel::new(1.0);
    system.recompute_energy(&model);
    system
}

/// Find the index of `spins` in the state list.
fn state_index(spins: &[f64], states: &[Vec<f64>]) -> Option<usize> {
    states.iter().position(|s| s == spins)
}

// ---------------------------------------------------------------------------
// B.1: Asymmetric Hastings proposal via custom ProposalStrategy on N=2
// ---------------------------------------------------------------------------

/// A proposal strategy that always proposes +1 with 80% probability
/// regardless of the current spin. The Hastings correction handles the asymmetry.
struct BiasedStrategy;

impl ProposalStrategy<IsingModel> for BiasedStrategy {
    fn propose(
        &mut self,
        _model: &IsingModel,
        system: &System,
        site: usize,
        rng: &mut impl rand::Rng,
    ) -> ProposedSpin {
        let current = system.spin_at(site, 1)[0];
        let proposed = if rng.random_range(0.0..1.0_f64) < 0.8 {
            1.0
        } else {
            -1.0
        };

        // q(new|old): probability of proposing `proposed` given `current`
        let q_new_given_old = if proposed == 1.0 { 0.8 } else { 0.2 };
        // q(old|new): probability of proposing `current` given `proposed` (reverse)
        let q_old_given_new = if current == 1.0 { 0.8 } else { 0.2 };
        // log_reverse_over_forward = ln q(old|new) - ln q(new|old)
        let ratio: f64 = q_old_given_new / q_new_given_old;
        let log_reverse_over_forward = ratio.ln();

        ProposedSpin {
            spin: Spin::from_slice(&[proposed]),
            log_reverse_over_forward,
        }
    }
}

#[test]
fn asymmetric_hastings_detailed_balance_n2() {
    let n = 2;
    let beta = 0.5;
    let model = IsingModel::new(1.0);
    let lattice = build_chain(n, true);
    let ensemble = CanonicalEnsemble::new(beta);
    let states = enumerate_ising_states(n);
    let n_states = states.len();
    let samples_per_state = 30_000;
    let tolerance = 0.03;

    let mut counts = vec![vec![0u64; n_states]; n_states];

    for (x_idx, x_spins) in states.iter().enumerate() {
        let mut system = build_system(x_spins, lattice.clone(), beta);
        let mut rng = RngType::seed_from_u64(42 + x_idx as u64);
        let mut strategy = BiasedStrategy;

        for _ in 0..samples_per_state {
            let site = rng.random_range(0..n);
            let proposed_spin = strategy.propose(&model, &system, site, &mut rng);

            let movement = SiteSpinMove::new(site, proposed_spin.spin);
            let proposal =
                cmc_rs::ProposedMove::new(movement, proposed_spin.log_reverse_over_forward);
            let mut patch = EnergyPatch::default();
            let outcome = cmc_rs::metropolis_hastings_step(
                &mut system,
                &model,
                &proposal,
                &ensemble,
                &mut patch,
                &mut rng,
            );

            let y_idx = if outcome.accepted {
                state_index(&system.spins, &states).unwrap_or(x_idx)
            } else {
                x_idx
            };
            counts[x_idx][y_idx] += 1;

            // Reset to x for next trial
            if outcome.accepted {
                system.spins.copy_from_slice(x_spins);
                system.recompute_energy(&model);
            }
        }
    }

    // π(x) P(x→y) ≈ π(y) P(y→x)
    let pi: Vec<f64> = states
        .iter()
        .map(|s| boltzmann(s, &model, &lattice, beta))
        .collect();
    let z: f64 = pi.iter().sum();
    let pi_norm: Vec<f64> = pi.iter().map(|p| p / z).collect();

    let mut max_violation = 0.0_f64;
    for x in 0..n_states {
        for y in x + 1..n_states {
            let p_xy = counts[x][y] as f64 / samples_per_state as f64;
            let p_yx = counts[y][x] as f64 / samples_per_state as f64;
            let forward = pi_norm[x] * p_xy;
            let reverse = pi_norm[y] * p_yx;
            let violation = (forward - reverse).abs();
            max_violation = max_violation.max(violation);
            assert!(
                violation < tolerance,
                "DB: π({x})P({x}→{y})={forward:.6} vs π({y})P({y}→{x})={reverse:.6}"
            );
        }
    }
    assert!(max_violation < tolerance);
}

// ---------------------------------------------------------------------------
// B.3: Parallel edges detailed balance on N=2
// ---------------------------------------------------------------------------

#[test]
fn parallel_edge_detailed_balance_n2() {
    let n = 2;
    let beta = 0.5;
    let bonds = vec![
        Bond::new(0, 1, BondType::Generic, 1.0),
        Bond::new(0, 1, BondType::Generic, 0.5),
    ];
    let lattice = CsrLattice::from_edges(n, bonds);
    let model = IsingModel::new(1.0);
    let ensemble = CanonicalEnsemble::new(beta);
    let states = enumerate_ising_states(n);
    let n_states = states.len();
    let samples_per_state = 30_000;
    let tolerance = 0.03;

    let mut counts = vec![vec![0u64; n_states]; n_states];

    for (x_idx, x_spins) in states.iter().enumerate() {
        let mut system = build_system(x_spins, lattice.clone(), beta);
        let mut rng = RngType::seed_from_u64(100 + x_idx as u64);

        for _ in 0..samples_per_state {
            let site = rng.random_range(0..n);
            let current = system.spin_at(site, 1)[0];
            let proposed = Spin::from_slice(&[-current]);

            let movement = SiteSpinMove::new(site, proposed);
            let proposal = cmc_rs::ProposedMove::symmetric(movement);
            let mut patch = EnergyPatch::default();
            let outcome = cmc_rs::metropolis_hastings_step(
                &mut system,
                &model,
                &proposal,
                &ensemble,
                &mut patch,
                &mut rng,
            );

            let y_idx = if outcome.accepted {
                state_index(&system.spins, &states).unwrap_or(x_idx)
            } else {
                x_idx
            };
            counts[x_idx][y_idx] += 1;

            if outcome.accepted {
                system.spins.copy_from_slice(x_spins);
                system.recompute_energy(&model);
            }
        }
    }

    let pi: Vec<f64> = states
        .iter()
        .map(|s| boltzmann(s, &model, &lattice, beta))
        .collect();
    let z: f64 = pi.iter().sum();
    let pi_norm: Vec<f64> = pi.iter().map(|p| p / z).collect();

    let mut max_violation = 0.0_f64;
    for x in 0..n_states {
        for y in x + 1..n_states {
            let p_xy = counts[x][y] as f64 / samples_per_state as f64;
            let p_yx = counts[y][x] as f64 / samples_per_state as f64;
            let forward = pi_norm[x] * p_xy;
            let reverse = pi_norm[y] * p_yx;
            let violation = (forward - reverse).abs();
            max_violation = max_violation.max(violation);
            assert!(violation < tolerance, "DB violated at parallel edges");
        }
    }
    assert!(max_violation < tolerance);
}

// ---------------------------------------------------------------------------
// B.5: Heat bath detailed balance on N=2
// ---------------------------------------------------------------------------

#[test]
fn heatbath_detailed_balance_n2() {
    let n = 2;
    let beta = 0.5;
    let model = IsingModel::new(1.0);
    let lattice = build_chain(n, true);
    let states = enumerate_ising_states(n);
    let n_states = states.len();
    let samples_per_state = 30_000;
    let tolerance = 0.03;

    let mut counts = vec![vec![0u64; n_states]; n_states];

    for (x_idx, x_spins) in states.iter().enumerate() {
        let mut system = build_system(x_spins, lattice.clone(), beta);
        let mut rng = RngType::seed_from_u64(200 + x_idx as u64);

        for _ in 0..samples_per_state {
            let site = rng.random_range(0..n);
            let proposed = model.heat_bath_sample_site(
                &system.spins,
                &system.lattice,
                site,
                system.beta,
                &mut rng,
            );

            let movement = SiteSpinMove::new(site, proposed);
            let mut patch = EnergyPatch::default();
            system.evaluate_trial(&model, &movement, &mut patch);
            <System as TrialEvaluator<IsingModel, SiteSpinMove>>::commit_trial(
                &mut system,
                &movement,
                &patch,
            );

            let y_idx = state_index(&system.spins, &states).unwrap_or(x_idx);
            counts[x_idx][y_idx] += 1;

            system.spins.copy_from_slice(x_spins);
            system.recompute_energy(&model);
        }
    }

    let pi: Vec<f64> = states
        .iter()
        .map(|s| boltzmann(s, &model, &lattice, beta))
        .collect();
    let z: f64 = pi.iter().sum();
    let pi_norm: Vec<f64> = pi.iter().map(|p| p / z).collect();

    let mut max_violation = 0.0_f64;
    for x in 0..n_states {
        for y in x + 1..n_states {
            let p_xy = counts[x][y] as f64 / samples_per_state as f64;
            let p_yx = counts[y][x] as f64 / samples_per_state as f64;
            let forward = pi_norm[x] * p_xy;
            let reverse = pi_norm[y] * p_yx;
            let violation = (forward - reverse).abs();
            max_violation = max_violation.max(violation);
            assert!(
                violation < tolerance,
                "Heat bath DB: forward={forward:.6} reverse={reverse:.6}"
            );
        }
    }
    assert!(max_violation < tolerance);
}

// ---------------------------------------------------------------------------
// B.4: Wolff cluster detailed balance on N=3
// ---------------------------------------------------------------------------

#[test]
fn wolff_detailed_balance_n3() {
    let n = 3;
    let beta = 0.5;
    let model = IsingModel::new(1.0);
    let lattice = build_chain(n, true);
    let states = enumerate_ising_states(n);
    let n_states = states.len();
    let samples = 50_000;
    let tolerance = 0.04;

    let mut counts = vec![vec![0u64; n_states]; n_states];

    for (x_idx, x_spins) in states.iter().enumerate() {
        let mut system = build_system(x_spins, lattice.clone(), beta);
        let mut rng = RngType::seed_from_u64(300 + x_idx as u64);
        let mut wolff = WolffCore::new();

        for _ in 0..samples {
            wolff.sweep(&mut system, &model, &mut rng);

            let y_idx = state_index(&system.spins, &states).unwrap_or(x_idx);
            counts[x_idx][y_idx] += 1;

            system.spins.copy_from_slice(x_spins);
            system.recompute_energy(&model);
        }
    }

    let pi: Vec<f64> = states
        .iter()
        .map(|s| boltzmann(s, &model, &lattice, beta))
        .collect();
    let z: f64 = pi.iter().sum();
    let pi_norm: Vec<f64> = pi.iter().map(|p| p / z).collect();

    let mut max_violation = 0.0_f64;
    for x in 0..n_states {
        for y in x + 1..n_states {
            let p_xy = counts[x][y] as f64 / samples as f64;
            let p_yx = counts[y][x] as f64 / samples as f64;
            let forward = pi_norm[x] * p_xy;
            let reverse = pi_norm[y] * p_yx;
            let violation = (forward - reverse).abs();
            max_violation = max_violation.max(violation);
            assert!(
                violation < tolerance,
                "Wolff DB: π({x})P({x}→{y})={forward:.6} vs π({y})P({y}→{x})={reverse:.6}"
            );
        }
    }
    assert!(max_violation < tolerance);
}

// ---------------------------------------------------------------------------
// B.2: Batch move detailed balance on N=3
// ---------------------------------------------------------------------------

#[test]
fn batch_move_detailed_balance_n3() {
    let n = 3;
    let beta = 0.5;
    let model = IsingModel::new(1.0);
    let lattice = build_chain(n, true);
    let states = enumerate_ising_states(n);
    let n_states = states.len();
    let samples_per_state = 20_000;
    let tolerance = 0.04;

    let mut counts = vec![vec![0u64; n_states]; n_states];

    for (x_idx, x_spins) in states.iter().enumerate() {
        let mut system = build_system(x_spins, lattice.clone(), beta);
        let mut rng = RngType::seed_from_u64(400 + x_idx as u64);

        for _ in 0..samples_per_state {
            let do_batch = rng.random_range(0.0..1.0_f64) < 0.3;
            if do_batch {
                let mut movement = cmc_rs::BatchSpinMove::new(1);
                for site in 0..n {
                    let flipped = -system.spin_at(site, 1)[0];
                    movement.push(site, &[flipped]);
                }
                let mut patch = cmc_rs::BatchEnergyPatch::default();
                system.evaluate_trial(&model, &movement, &mut patch);
                <System as TrialEvaluator<IsingModel, cmc_rs::BatchSpinMove>>::commit_trial(
                    &mut system,
                    &movement,
                    &patch,
                );
            }

            let y_idx = state_index(&system.spins, &states).unwrap_or(x_idx);
            counts[x_idx][y_idx] += 1;

            system.spins.copy_from_slice(x_spins);
            system.recompute_energy(&model);
        }
    }

    let pi: Vec<f64> = states
        .iter()
        .map(|s| boltzmann(s, &model, &lattice, beta))
        .collect();
    let z: f64 = pi.iter().sum();
    let pi_norm: Vec<f64> = pi.iter().map(|p| p / z).collect();

    let mut max_violation = 0.0_f64;
    for x in 0..n_states {
        for y in x + 1..n_states {
            let p_xy = counts[x][y] as f64 / samples_per_state as f64;
            let p_yx = counts[y][x] as f64 / samples_per_state as f64;
            let forward = pi_norm[x] * p_xy;
            let reverse = pi_norm[y] * p_yx;
            let violation = (forward - reverse).abs();
            max_violation = max_violation.max(violation);
            assert!(
                violation < tolerance,
                "Batch move DB: π({x})P({x}→{y})={forward:.6} vs π({y})P({y}→{x})={reverse:.6}"
            );
        }
    }
    assert!(max_violation < tolerance);
}

// ---------------------------------------------------------------------------
// B.3b: Self-loop detailed balance on N=2
// ---------------------------------------------------------------------------

#[test]
fn self_loop_detailed_balance_n2() {
    let n = 2;
    let beta = 0.5;
    let bonds = vec![
        Bond::new(0, 1, BondType::Generic, 1.0),
        Bond::new(0, 0, BondType::Generic, 0.5),
    ];
    let lattice = CsrLattice::from_edges(n, bonds);
    let model = IsingModel::new(1.0);
    let ensemble = CanonicalEnsemble::new(beta);
    let states = enumerate_ising_states(n);
    let n_states = states.len();
    let samples_per_state = 30_000;
    let tolerance = 0.03;

    let mut counts = vec![vec![0u64; n_states]; n_states];

    for (x_idx, x_spins) in states.iter().enumerate() {
        let mut system = build_system(x_spins, lattice.clone(), beta);
        let mut rng = RngType::seed_from_u64(500 + x_idx as u64);

        for _ in 0..samples_per_state {
            let site = rng.random_range(0..n);
            let current = system.spin_at(site, 1)[0];
            let proposed = Spin::from_slice(&[-current]);

            let movement = SiteSpinMove::new(site, proposed);
            let proposal = cmc_rs::ProposedMove::symmetric(movement);
            let mut patch = EnergyPatch::default();
            let outcome = cmc_rs::metropolis_hastings_step(
                &mut system,
                &model,
                &proposal,
                &ensemble,
                &mut patch,
                &mut rng,
            );

            let y_idx = if outcome.accepted {
                state_index(&system.spins, &states).unwrap_or(x_idx)
            } else {
                x_idx
            };
            counts[x_idx][y_idx] += 1;

            if outcome.accepted {
                system.spins.copy_from_slice(x_spins);
                system.recompute_energy(&model);
            }
        }
    }

    let pi: Vec<f64> = states
        .iter()
        .map(|s| boltzmann(s, &model, &lattice, beta))
        .collect();
    let z: f64 = pi.iter().sum();
    let pi_norm: Vec<f64> = pi.iter().map(|p| p / z).collect();

    let mut max_violation = 0.0_f64;
    for x in 0..n_states {
        for y in x + 1..n_states {
            let p_xy = counts[x][y] as f64 / samples_per_state as f64;
            let p_yx = counts[y][x] as f64 / samples_per_state as f64;
            let forward = pi_norm[x] * p_xy;
            let reverse = pi_norm[y] * p_yx;
            let violation = (forward - reverse).abs();
            max_violation = max_violation.max(violation);
            assert!(
                violation < tolerance,
                "Self-loop DB: forward={forward:.6} reverse={reverse:.6}"
            );
        }
    }
    assert!(max_violation < tolerance);
}
