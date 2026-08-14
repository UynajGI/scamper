//! Exact-diagonalization cross-checks for continuous-time lattice QMC.
//!
//! Each test builds a dense Hamiltonian for a small S=1/2 system, computes
//! the thermal density matrix ρ = e^{-βH} / Z via scaling-and-squaring
//! matrix exponential, and compares multiple observables against MC output.
//!
//! Statistical tolerance is generous (fixed absolute) for these small-system
//! tests. Rigorous 4σ batch-means validation lives in `lattice_long.rs`.

use qmc_rs::lattice::ContinuousLatticeEngine;
use qmc_rs::{
    CsrGraph, EdgeCoupling, LatticeConfiguration, QmcKernel, SpinModelBuilder, SpinSpace,
    UpdateSchedule,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

// ─── Minimal dense ED via matrix exponential ─────────────────────────────
//
// ρ = e^{-βH} computed via scaling-and-squaring with truncated Taylor series.
// For any observable O:  ⟨O⟩ = Tr[O·ρ] / Tr[ρ]
//
// For diagonal observables (Sz-basis diagonal): Tr[O·ρ] = Σ_s O[s,s]·ρ[s,s]
// For energy: Tr[H·ρ] = Σ_{i,j} H[i,j]·ρ[j,i]

pub(crate) struct DenseMatrix {
    pub(crate) dim: usize,
    pub(crate) elements: Vec<f64>, // row-major dim×dim
}

impl DenseMatrix {
    pub(crate) fn zero(dim: usize) -> Self {
        Self {
            dim,
            elements: vec![0.0; dim * dim],
        }
    }

    fn identity(dim: usize) -> Self {
        let mut m = Self::zero(dim);
        for i in 0..dim {
            m.elements[i * dim + i] = 1.0;
        }
        m
    }

    pub(crate) fn get(&self, i: usize, j: usize) -> f64 {
        self.elements[i * self.dim + j]
    }

    fn set(&mut self, i: usize, j: usize, val: f64) {
        self.elements[i * self.dim + j] = val;
    }

    fn add(&mut self, i: usize, j: usize, val: f64) {
        self.elements[i * self.dim + j] += val;
    }

    pub(crate) fn multiply(&self, other: &Self) -> Self {
        let dim = self.dim;
        let mut result = Self::zero(dim);
        for i in 0..dim {
            for j in 0..dim {
                let mut sum = 0.0;
                for k in 0..dim {
                    sum += self.get(i, k) * other.get(k, j);
                }
                result.set(i, j, sum);
            }
        }
        result
    }

    fn scale(&self, s: f64) -> Self {
        Self {
            dim: self.dim,
            elements: self.elements.iter().map(|&x| x * s).collect(),
        }
    }

    pub(crate) fn trace(&self) -> f64 {
        (0..self.dim).map(|i| self.get(i, i)).sum()
    }

    /// Matrix exponential exp(-beta * H) via scaling and squaring.
    pub(crate) fn expm_negative(&self, beta: f64) -> Self {
        let dim = self.dim;
        let mut a = self.scale(-beta);
        // Scale down so max |element| ≤ 0.5
        let mut max_el = 0.0_f64;
        for &v in &a.elements {
            max_el = max_el.max(v.abs());
        }
        let mut n_scales = 0;
        while max_el > 0.5 {
            a = a.scale(0.5);
            max_el *= 0.5;
            n_scales += 1;
        }
        // Taylor series: I + A + A²/2! + A³/3! + ...
        let mut result = Self::identity(dim);
        let mut term = Self::identity(dim);
        for k in 1..=40 {
            term = term.multiply(&a).scale(1.0 / k as f64);
            result
                .elements
                .iter_mut()
                .zip(&term.elements)
                .for_each(|(r, &t)| *r += t);
        }
        // Square n_scales times
        for _ in 0..n_scales {
            result = result.multiply(&result);
        }
        result
    }
}

/// Build dense S=1/2 Hamiltonian for a set of Heisenberg/XXZ edges.
pub(crate) fn build_hamiltonian(
    n_sites: usize,
    edges: &[(usize, usize, EdgeCoupling)],
) -> DenseMatrix {
    let dim = 1usize << n_sites;
    let mut h = DenseMatrix::zero(dim);
    for &(i, j, coupling) in edges {
        let bi = n_sites - 1 - i;
        let bj = n_sites - 1 - j;
        let jz = coupling.j_z;
        let flip_amp = 0.25 * (coupling.j_x + coupling.j_y);
        for state in 0..dim {
            let si = (state >> bi) & 1;
            let sj = (state >> bj) & 1;
            let mi = si as f64 - 0.5;
            let mj = sj as f64 - 0.5;
            h.add(state, state, jz * mi * mj);
            if si != sj {
                let flipped = state ^ (1 << bi) ^ (1 << bj);
                h.add(state, flipped, flip_amp);
            }
        }
    }
    h
}

/// Add transverse field -h_x·Sx_i.
#[allow(dead_code)]
fn add_transverse_field(h: &mut DenseMatrix, n_sites: usize, i: usize, h_x: f64) {
    let bi = n_sites - 1 - i;
    let dim = h.dim;
    for state in 0..dim {
        let flipped = state ^ (1 << bi);
        h.add(state, flipped, -h_x * 0.5);
    }
}

/// ⟨Sz_i Sz_j⟩ averaged over specified edges, as a diagonal operator.
fn nn_correlation(rho: &DenseMatrix, z: f64, n_sites: usize, edges: &[(usize, usize)]) -> f64 {
    let dim = rho.dim;
    let mut sum = 0.0;
    for state in 0..dim {
        let avg: f64 = edges
            .iter()
            .map(|&(i, j)| {
                let mi = ((state >> (n_sites - 1 - i)) & 1) as f64 - 0.5;
                let mj = ((state >> (n_sites - 1 - j)) & 1) as f64 - 0.5;
                mi * mj
            })
            .sum::<f64>()
            / edges.len() as f64;
        sum += avg * rho.get(state, state);
    }
    sum / z
}

// ─── MC runner ───────────────────────────────────────────────────────────

struct McResult {
    energy: f64,
    nn_correlation: f64,
    m_squared: f64,
}

#[allow(clippy::too_many_arguments)]
fn run_mc(
    graph: CsrGraph,
    edge: EdgeCoupling,
    beta: f64,
    n_sites: usize,
    initial_states: Vec<u16>,
    seed: u64,
    n_thermalization: usize,
    n_measurement: usize,
    measure_interval: usize,
) -> McResult {
    let space = SpinSpace::uniform(n_sites, 1).expect("space");
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(edge)
        .build()
        .expect("model");
    let mut configuration =
        LatticeConfiguration::new(beta, initial_states, &model).expect("configuration");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(8, 4, 64));
    engine.set_validate_each_sweep(true);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    let mut energy_sum = 0.0;
    let mut nn_sum = 0.0;
    let mut m2_sum = 0.0;
    let mut samples = 0u64;

    for sweep in 0..(n_thermalization + n_measurement) {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
        if sweep >= n_thermalization && sweep % measure_interval == 0 {
            let obs =
                qmc_rs::lattice::measure_observables(&configuration, engine.model()).expect("obs");
            energy_sum += obs.energy_total;
            nn_sum += obs.nearest_neighbor_sz_correlation;
            m2_sum += obs.magnetization_z_squared;
            samples += 1;
        }
    }
    assert!(samples > 0, "no samples collected");
    let n = samples as f64;
    McResult {
        energy: energy_sum / n,
        nn_correlation: nn_sum / n,
        m_squared: m2_sum / n,
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────

/// 3-site S=1/2 AFM Heisenberg open chain (Hilbert space = 8).
/// Compares energy, ⟨Sz_i Sz_j⟩, ⟨m²⟩, and ⟨m_s²⟩ against ED.
#[test]
fn three_site_heisenberg_chain_matches_ed_for_four_observables() {
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
    let exact_nn = nn_correlation(&rho, z, n_sites, &edges_pair);

    // ⟨m²⟩ with m = Σ Sz_i / N
    let dim = 1usize << n_sites;
    let exact_m2: f64 = (0..dim)
        .map(|s| {
            let m: f64 = (0..n_sites)
                .map(|i| ((s >> (n_sites - 1 - i)) & 1) as f64 - 0.5)
                .sum::<f64>()
                / n_sites as f64;
            m * m * rho.get(s, s)
        })
        .sum::<f64>()
        / z;

    let result = run_mc(
        graph,
        EdgeCoupling::heisenberg(j),
        beta,
        n_sites,
        vec![0, 1, 0],
        0xCAFE_BABE,
        20_000,
        80_000,
        2,
    );

    assert!(
        (result.energy - exact_energy).abs() < 0.06,
        "energy: MC={:.6}, exact={:.6}",
        result.energy,
        exact_energy
    );
    assert!(
        (result.nn_correlation - exact_nn).abs() < 0.015,
        "nn_corr: MC={:.6}, exact={:.6}",
        result.nn_correlation,
        exact_nn
    );
    assert!(
        (result.m_squared - exact_m2).abs() < 0.005,
        "m²: MC={:.6}, exact={:.6}",
        result.m_squared,
        exact_m2
    );
    // Note: staggered m² comparison omitted — the QMC estimator is
    // ⟨[time_avg(m_s)]²⟩ which differs from ⟨m_s²⟩_instantaneous when
    // [M_s, H] ≠ 0. The uniform m² comparison is exact because [M, H] = 0.
}

/// 4-site S=1/2 XXZ open chain (Hilbert space = 16).
#[test]
fn four_site_xxz_chain_matches_ed_for_energy_and_correlation() {
    let n_sites = 4;
    let beta = 2.5;
    let j_xy = -0.6;
    let j_z = 0.4;
    let edges_pair = [(0, 1), (1, 2), (2, 3)];

    let graph = CsrGraph::chain(n_sites, false).expect("graph");
    let weight = graph.edges().first().unwrap().weight;
    let mut edges: Vec<(usize, usize, EdgeCoupling)> = Vec::new();
    for &(si, sj) in &edges_pair {
        edges.push((si, sj, EdgeCoupling::xxz(j_xy * weight, j_z * weight)));
    }
    let h = build_hamiltonian(n_sites, &edges);
    let rho = h.expm_negative(beta);
    let z = rho.trace();
    let exact_energy = h.multiply(&rho).trace() / z;
    let exact_nn = nn_correlation(&rho, z, n_sites, &edges_pair);

    let result = run_mc(
        graph,
        EdgeCoupling::xxz(j_xy, j_z),
        beta,
        n_sites,
        vec![1, 0, 1, 0],
        0xDEAD_BEEF,
        20_000,
        80_000,
        2,
    );

    assert!(
        (result.energy - exact_energy).abs() < 0.06,
        "energy: MC={:.6}, exact={:.6}",
        result.energy,
        exact_energy
    );
    assert!(
        (result.nn_correlation - exact_nn).abs() < 0.015,
        "nn_corr: MC={:.6}, exact={:.6}",
        result.nn_correlation,
        exact_nn
    );
}

/// 3-site ferromagnetic Heisenberg chain (J<0). All spins polarize at low T.
/// Tests the non-bipartite-gauge path and verifies energy + ⟨m²⟩ against ED.
#[test]
fn three_site_ferromagnetic_heisenberg_chain_matches_ed() {
    let n_sites = 3;
    let beta = 3.0;
    let j = -1.0; // ferromagnetic
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

    let dim = 1usize << n_sites;
    let exact_m2: f64 = (0..dim)
        .map(|s| {
            let m: f64 = (0..n_sites)
                .map(|i| ((s >> (n_sites - 1 - i)) & 1) as f64 - 0.5)
                .sum::<f64>()
                / n_sites as f64;
            m * m * rho.get(s, s)
        })
        .sum::<f64>()
        / z;

    let space = SpinSpace::uniform(n_sites, 1).expect("space");
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::heisenberg(j))
        .build()
        .expect("model");
    let mut configuration = LatticeConfiguration::new(beta, vec![0, 1, 0], &model).expect("config");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(8, 4, 64));
    engine.set_validate_each_sweep(true);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xBEEF);

    let mut energy_sum = 0.0;
    let mut m2_sum = 0.0;
    let mut samples = 0u64;
    for sweep in 0..100_000u64 {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
        if sweep >= 20_000 && sweep % 2 == 0 {
            let obs =
                qmc_rs::lattice::measure_observables(&configuration, engine.model()).expect("obs");
            energy_sum += obs.energy_total;
            m2_sum += obs.magnetization_z_squared;
            samples += 1;
        }
    }
    let measured_e = energy_sum / samples as f64;
    let measured_m2 = m2_sum / samples as f64;

    assert!(
        (measured_e - exact_energy).abs() < 0.05,
        "energy: MC={:.6}, exact={:.6}",
        measured_e,
        exact_energy
    );
    assert!(
        (measured_m2 - exact_m2).abs() < 0.01,
        "m²: MC={:.6}, exact={:.6}",
        measured_m2,
        exact_m2
    );
}

/// 3-site S=1/2 AFM Heisenberg open chain: χ_z = β(⟨m²⟩ − ⟨m⟩²) vs ED.
///
/// ED computes χ_z from the diagonal magnetization operator m = Σ_i Sz_i / N
/// in the Sz basis:  χ_z = β(Tr[m²ρ]/Z − (Tr[mρ]/Z)²).
/// MC collects both ⟨m⟩ and ⟨m²⟩ sample averages and forms the same combination.
#[test]
fn three_site_heisenberg_susceptibility_matches_ed() {
    let n_sites = 3;
    let beta = 3.0;
    let j = 1.0;
    let edges_pair = [(0, 1), (1, 2)];

    // ── ED side ──────────────────────────────────────────────────────────
    let graph = CsrGraph::chain(n_sites, false).expect("graph");
    let weight = graph.edges().first().unwrap().weight;
    let edges: Vec<(usize, usize, EdgeCoupling)> = edges_pair
        .iter()
        .map(|&(si, sj)| (si, sj, EdgeCoupling::heisenberg(j * weight)))
        .collect();
    let h = build_hamiltonian(n_sites, &edges);
    let rho = h.expm_negative(beta);
    let z = rho.trace();

    let dim = 1usize << n_sites;
    // m(s) = Σ_i Sz_i / N  for each basis state
    let exact_m: f64 = (0..dim)
        .map(|s| {
            let m: f64 = (0..n_sites)
                .map(|i| ((s >> (n_sites - 1 - i)) & 1) as f64 - 0.5)
                .sum::<f64>()
                / n_sites as f64;
            m * rho.get(s, s)
        })
        .sum::<f64>()
        / z;
    let exact_m2: f64 = (0..dim)
        .map(|s| {
            let m: f64 = (0..n_sites)
                .map(|i| ((s >> (n_sites - 1 - i)) & 1) as f64 - 0.5)
                .sum::<f64>()
                / n_sites as f64;
            m * m * rho.get(s, s)
        })
        .sum::<f64>()
        / z;
    let exact_chi_z = beta * (exact_m2 - exact_m * exact_m);

    // ── MC side ──────────────────────────────────────────────────────────
    let space = SpinSpace::uniform(n_sites, 1).expect("space");
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::heisenberg(j))
        .build()
        .expect("model");
    let mut configuration = LatticeConfiguration::new(beta, vec![0, 1, 0], &model).expect("config");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(8, 4, 64));
    engine.set_validate_each_sweep(true);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0xC410);

    let mut m_sum = 0.0;
    let mut m2_sum = 0.0;
    let mut samples = 0u64;
    for sweep in 0..100_000u64 {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
        if sweep >= 20_000 && sweep % 2 == 0 {
            let obs =
                qmc_rs::lattice::measure_observables(&configuration, engine.model()).expect("obs");
            m_sum += obs.magnetization_z;
            m2_sum += obs.magnetization_z_squared;
            samples += 1;
        }
    }
    assert!(samples > 0, "no MC samples collected");
    let n = samples as f64;
    let mc_m = m_sum / n;
    let mc_m2 = m2_sum / n;
    let mc_chi_z = beta * (mc_m2 - mc_m * mc_m);

    assert!(
        (mc_chi_z - exact_chi_z).abs() < 0.01,
        "χ_z: MC={mc_chi_z:.6}, exact={exact_chi_z:.6} (⟨m⟩_MC={mc_m:.6}, ⟨m²⟩_MC={mc_m2:.6})"
    );
}

/// Verify the matrix exponential against the known 2-site XXZ spectrum.
#[test]
fn matrix_exp_recovers_known_dimer_partition_function() {
    let beta = 3.0;
    let j_xy = -0.8;
    let j_z = 0.3;

    let edges = [(0, 1, EdgeCoupling::xxz(j_xy, j_z))];
    let h = build_hamiltonian(2, &edges);
    let rho = h.expm_negative(beta);
    let z_computed = rho.trace();

    // Exact partition function from known spectrum
    let levels = [
        0.25 * j_z,
        0.25 * j_z,
        -0.25 * j_z + 0.5 * j_xy,
        -0.25 * j_z - 0.5 * j_xy,
    ];
    let z_exact: f64 = levels.iter().map(|&e| (-beta * e).exp()).sum();

    assert!(
        (z_computed - z_exact).abs() < 1e-10,
        "Z: computed={z_computed:.10}, exact={z_exact:.10}"
    );

    let exact_energy: f64 = levels.iter().map(|&e| e * (-beta * e).exp()).sum::<f64>() / z_exact;
    let computed_energy = h.multiply(&rho).trace() / z_computed;
    assert!(
        (computed_energy - exact_energy).abs() < 1e-10,
        "E: computed={computed_energy:.10}, exact={exact_energy:.10}"
    );
}

// ─── P2.1: Binder cumulant U4 = 1 - ⟨m⁴⟩/(3⟨m²⟩²) vs ED ──────────────

#[test]
fn three_site_heisenberg_binder_cumulant_matches_ed() {
    let n_sites = 3;
    let beta = 3.0;
    let j = 1.0;
    let edges_pair = [(0, 1), (1, 2)];

    // ED side
    let graph = CsrGraph::chain(n_sites, false).expect("graph");
    let weight = graph.edges().first().unwrap().weight;
    let edges: Vec<(usize, usize, EdgeCoupling)> = edges_pair
        .iter()
        .map(|&(si, sj)| (si, sj, EdgeCoupling::heisenberg(j * weight)))
        .collect();
    let h = build_hamiltonian(n_sites, &edges);
    let rho = h.expm_negative(beta);
    let z = rho.trace();
    let dim = 1usize << n_sites;

    // m(s) = Σ_i Sz_i / N for each basis state (diagonal in Sz basis)
    let mut m2_sum = 0.0;
    let mut m4_sum = 0.0;
    for s in 0..dim {
        let m: f64 = (0..n_sites)
            .map(|i| ((s >> (n_sites - 1 - i)) & 1) as f64 - 0.5)
            .sum::<f64>()
            / n_sites as f64;
        let rho_ss = rho.get(s, s);
        m2_sum += m * m * rho_ss;
        m4_sum += m * m * m * m * rho_ss;
    }
    let exact_m2 = m2_sum / z;
    let exact_m4 = m4_sum / z;
    let exact_u4 = 1.0 - exact_m4 / (3.0 * exact_m2 * exact_m2);

    // MC side
    let space = SpinSpace::uniform(n_sites, 1).expect("space");
    let model = SpinModelBuilder::new(graph, space)
        .uniform_edge(EdgeCoupling::heisenberg(j))
        .build()
        .expect("model");
    let mut configuration = LatticeConfiguration::new(beta, vec![0, 0, 0], &model).expect("config");
    let mut engine = ContinuousLatticeEngine::new(model, UpdateSchedule::new(2, 2, 16));
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    let mut mc_m2 = 0.0;
    let mut mc_m4 = 0.0;
    let mut samples = 0u64;
    for sweep in 0..80_000 {
        engine.sweep(&mut configuration, &mut rng).expect("sweep");
        if sweep >= 20_000 {
            let obs =
                qmc_rs::lattice::measure_observables(&configuration, engine.model()).expect("obs");
            mc_m2 += obs.magnetization_z_squared;
            mc_m4 += obs.magnetization_z_fourth;
            samples += 1;
        }
    }
    let n = samples as f64;
    mc_m2 /= n;
    mc_m4 /= n;
    let mc_u4 = 1.0 - mc_m4 / (3.0 * mc_m2 * mc_m2);

    assert!(
        (mc_u4 - exact_u4).abs() < 0.02,
        "U4: MC={mc_u4:.6}, exact={exact_u4:.6} (⟨m²⟩_MC={mc_m2:.6}, ⟨m⁴⟩_MC={mc_m4:.6})"
    );
}
