//! Gap 2: Explicit Markov chain connectivity / ergodicity check.
//!
//! For a 2-site PBC Ising model (4 states), we enumerate the full
//! state transition graph by running the MC kernel from each state
//! many times and recording all states ever reached. Then verify:
//!   - The graph is strongly connected (one communicating class)
//!   - Every state can reach every other state
//!   - The chain is aperiodic (self-loops exist)

use cmc_rs::{
    build_chain, Algorithm, Hamiltonian, IsingModel, MetropolisCore, SWCore, SimulationPhase,
    System, WolffCore,
};
use rand::{RngExt, SeedableRng};

type Rng = rand_xoshiro::Xoshiro256PlusPlus;

fn all_ising_states(n: usize) -> Vec<Vec<f64>> {
    (0..(1u32 << n))
        .map(|mask| {
            (0..n)
                .map(|i| if (mask >> i) & 1 == 1 { 1.0 } else { -1.0 })
                .collect()
        })
        .collect()
}

fn state_to_idx(spins: &[f64]) -> usize {
    spins
        .iter()
        .enumerate()
        .map(|(i, &s)| if s > 0.0 { 1usize << i } else { 0 })
        .sum()
}

fn build_system(spins: &[f64], lattice: &cmc_rs::CsrLattice, beta: f64) -> System {
    let mut system = System::new(lattice.clone(), 1, 1.0, beta);
    system.spins.copy_from_slice(spins);
    let model = IsingModel::new(1.0);
    system.recompute_energy(&model);
    system
}

/// Run `n_sequences` multi-sweep trajectories of length `trajectory_length`
/// from `start_state`, record all states visited across all sweeps.
fn reachable_states_multi<A: Algorithm<IsingModel>>(
    kernel: &mut A,
    start_spins: &[f64],
    lattice: &cmc_rs::CsrLattice,
    beta: f64,
    n_sequences: usize,
    trajectory_length: usize,
    seed: u64,
) -> std::collections::HashSet<usize> {
    let mut visited = std::collections::HashSet::new();
    let model = IsingModel::new(1.0);
    for seq in 0..n_sequences {
        let mut system = build_system(start_spins, lattice, beta);
        let mut rng = Rng::seed_from_u64(seed + seq as u64);
        for _ in 0..trajectory_length {
            kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
            visited.insert(state_to_idx(&system.spins));
        }
    }
    visited
}

/// Check that the transition graph is strongly connected.
fn assert_strongly_connected(
    transition_graph: &[std::collections::HashSet<usize>],
    n_states: usize,
    kernel_name: &str,
) {
    // BFS from state 0 — should reach all states
    let mut visited = vec![false; n_states];
    let mut queue = vec![0];
    visited[0] = true;
    while let Some(s) = queue.pop() {
        for &t in &transition_graph[s] {
            if !visited[t] {
                visited[t] = true;
                queue.push(t);
            }
        }
    }
    for (i, &v) in visited.iter().enumerate() {
        assert!(
            v,
            "{kernel_name}: state {i} unreachable from state 0 — not strongly connected"
        );
    }

    // Check aperiodicity: at least one self-loop exists
    let has_self_loop = transition_graph
        .iter()
        .enumerate()
        .any(|(i, reachable)| reachable.contains(&i));
    assert!(
        has_self_loop,
        "{kernel_name}: no self-loops found — chain may be periodic"
    );
}

const N: usize = 2;
const BETA: f64 = 0.5;
const N_TRIALS: usize = 5000;

#[test]
fn metropolis_markov_chain_is_strongly_connected() {
    let lattice = build_chain(N, true);
    let states = all_ising_states(N);
    let n_states = states.len();

    let mut kernel = MetropolisCore::new();
    let mut graph: Vec<std::collections::HashSet<usize>> = Vec::new();
    for spins in &states {
        graph.push(reachable_states_multi(
            &mut kernel,
            spins,
            &lattice,
            BETA,
            N_TRIALS,
            10,
            0xA001,
        ));
    }

    assert_strongly_connected(&graph, n_states, "Metropolis");
}

#[test]
fn wolff_markov_chain_is_strongly_connected() {
    let lattice = build_chain(N, true);
    let states = all_ising_states(N);
    let n_states = states.len();

    let mut kernel = WolffCore::new();
    let mut graph: Vec<std::collections::HashSet<usize>> = Vec::new();
    for spins in &states {
        graph.push(reachable_states_multi(
            &mut kernel,
            spins,
            &lattice,
            BETA,
            N_TRIALS,
            20,
            0xA002,
        ));
    }

    assert_strongly_connected(&graph, n_states, "Wolff");
}

#[test]
fn sw_markov_chain_is_strongly_connected() {
    let lattice = build_chain(N, true);
    let states = all_ising_states(N);
    let n_states = states.len();

    let mut kernel = SWCore::new();
    let mut graph: Vec<std::collections::HashSet<usize>> = Vec::new();
    for spins in &states {
        graph.push(reachable_states_multi(
            &mut kernel,
            spins,
            &lattice,
            BETA,
            N_TRIALS,
            20,
            0xA003,
        ));
    }

    assert_strongly_connected(&graph, n_states, "SW");
}

#[test]
fn metropolis_all_states_visited_from_any_start() {
    // Stronger check: from every starting state, all 4 states are reachable
    // within a multi-sweep trajectory.
    let lattice = build_chain(N, true);
    let states = all_ising_states(N);
    let n_states = states.len();
    let mut kernel = MetropolisCore::new();

    for (i, spins) in states.iter().enumerate() {
        let reachable = reachable_states_multi(
            &mut kernel,
            spins,
            &lattice,
            BETA,
            500,
            50,
            0xB000 + i as u64,
        );
        assert_eq!(
            reachable.len(),
            n_states,
            "Metropolis from state {i}: only reached {}/{n_states} states",
            reachable.len()
        );
    }
}
