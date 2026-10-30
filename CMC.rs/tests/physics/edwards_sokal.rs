//! R3 Edwards-Sokal FK exposure: exact validation gates (roadmap §15 R3).
//!
//! # The worm separation note (spec, verbatim requirement)
//!
//! The worm module (`src/worm/ising.rs`) samples the Ising high-temperature
//! graph representation: edge subsets constrained to the **even subgraphs**
//! (every vertex has even degree; the worm sector carries two odd defects),
//! weighted by `prod_e tanh(beta J_e)^{n_e}`. An FK configuration is a
//! DIFFERENT object on a different support: an unrestricted edge subset whose
//! same-spin edges open independently with probability
//! `1 - exp(-2 beta J_e)` and whose measure carries the `2^{k(A)}`
//! spin-symmetry factor. The worm's even subgraph is **not** an FK
//! configuration and must never be relabeled as one; nothing here analyzes
//! worm output through the FK path. The two representations are compared
//! only where they must agree — Ising observables: the partition function
//! (`fk_double_sum_matches_exact_zero_field_ising_partition_function`
//! cross-checks the spin-sum and the even-subgraph high-temperature sum in
//! one gate) and two-point correlations
//! (`fk_connection_probability_equals_spin_correlation` and
//! `worm_representation_correlations_agree_with_spin_enumeration`).
//!
//! # Gates
//!
//! 1. **Exact bond-subset enumeration** (frozen spins): the conditional FK
//!    law given a fixed spin configuration is the product Bernoulli measure
//!    over same-spin edges; all `2^E` subsets are enumerated per
//!    configuration and the sampler's empirical frequencies are z-gated
//!    against the exact law, with the R1 `BondSample`/`analyze_bonds` paths
//!    verified to analyze every subset identically.
//! 2. **Partition function**:
//!    `Z_Ising = exp(beta sum_e J_e) * sum_s sum_A prod p^{|A|} (1-p)^{E-|A|}`
//!    over compatible `(s, A)` pairs — exhaustively on graphs of at most 4
//!    edges (including a parallel-edge + self-loop multigraph), with `k(A)`
//!    from a test-local union-find and the Ising reference from independent
//!    brute-force spin enumeration of `exp(-beta H)`. The double sum and the
//!    random-cluster single sum `sum_A p^{|A|} (1-p)^{E-|A|} 2^{k(A)}` are
//!    the same object and are checked to agree.
//! 3. **`P_FK(i ~ j) = <s_i s_j>`**: a proper Edwards-Sokal chain on the
//!    4-ring; the boundary-crossing mean of the exposed FK draw and the spin
//!    correlator of the SAME chain are each z-gated against the exact
//!    brute-force correlator, and to each other.
//!
//! Draw counts are env-tunable through `SCUTTLE_EDWARDS_DRAWS`.

use std::collections::HashMap;

use super::common::assert_close;
use carlo_rs::{Context, MonteCarlo};
use cmc_rs::{
    analyze, analyze_bonds, build_chain, AllActiveVertices, Bond, BondSample, BondType,
    BoundaryQuery, ComponentWorkspace, CsrLattice, EdwardsSokalMC, GraphView, ObservablePlan,
    StaticObservable, UndirectedGraphView,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Environment variable overriding the per-test draw count.
const DRAWS_ENV: &str = "SCUTTLE_EDWARDS_DRAWS";
/// Upper bound accepted for [`DRAWS_ENV`] — typo protection for nightly runs.
const DRAWS_MAX: usize = 10_000_000;

/// Monte Carlo draws for the statistical gates (unset → `default`).
fn draw_count(default: usize) -> usize {
    let raw = match std::env::var(DRAWS_ENV) {
        Ok(raw) => raw,
        Err(std::env::VarError::NotPresent) => return default,
        Err(std::env::VarError::NotUnicode(raw)) => {
            panic!("{DRAWS_ENV} is not valid Unicode: {raw:?}")
        }
    };
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return default;
    }
    match trimmed.parse::<usize>() {
        Ok(count) if (1..=DRAWS_MAX).contains(&count) => {
            eprintln!("[edwards-sokal] {DRAWS_ENV}={count} (default {default})");
            count
        }
        _ => panic!(
            "{DRAWS_ENV}={raw:?}: expected an integer in 1..={DRAWS_MAX} \
             (unset the variable to use the default of {default} draws)"
        ),
    }
}

/// Labeled seed base ("EDSO") for the statistical gates.
const SEED_BASE: u64 = 0x4544_534f_0000;

/// 4 vertices, 4 heterogeneous edges: 0-1 (J=1.0), 1-2 (J=0.7),
/// 2-3 (J=1.3), 0-2 (J=0.5).
fn wedge_graph() -> CsrLattice {
    CsrLattice::from_edges(
        4,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 0.7),
            Bond::new(2, 3, BondType::Generic, 1.3),
            Bond::new(0, 2, BondType::Generic, 0.5),
        ],
    )
}

/// 3-vertex multigraph: parallel 0-1 pair, one 1-2 edge, one 2-2 self-loop.
fn multigraph() -> CsrLattice {
    CsrLattice::from_edges(
        3,
        vec![
            Bond::new(0, 1, BondType::Generic, 0.8),
            Bond::new(0, 1, BondType::Generic, 0.8),
            Bond::new(1, 2, BondType::Generic, 1.1),
            Bond::new(2, 2, BondType::Generic, 0.6),
        ],
    )
}

/// Per-edge FK open probabilities `1 - exp(-2 beta J_e)`, dense order.
fn bond_probabilities(lattice: &CsrLattice, beta: f64) -> Vec<f64> {
    lattice
        .edges
        .iter()
        .map(|edge| 1.0 - (-2.0 * beta * edge.weight).exp())
        .collect()
}

fn spin_vector(vertex_count: usize, mask: usize) -> Vec<i8> {
    (0..vertex_count)
        .map(|site| if mask & (1 << site) != 0 { 1 } else { -1 })
        .collect()
}

// ── Gate 1: frozen-spins exact bond-subset enumeration ────────────────────
#[test]
fn frozen_bond_draw_matches_exact_subset_enumeration() {
    let lattice = wedge_graph();
    let beta = 0.9;
    let probabilities = bond_probabilities(&lattice, beta);
    let draws = draw_count(60_000);
    let queries = [BoundaryQuery::new(&lattice, &[0], &[3]).unwrap()];
    let plan = ObservablePlan::new([StaticObservable::ActiveEdgeCount], Vec::new());

    for config in [[1, 1, 1, 1], [1, 1, -1, -1]] {
        // Eligible (same-spin) edges carry the product Bernoulli law; the
        // remaining edges are closed with probability one in every draw.
        let eligible: Vec<usize> = (0..lattice.edge_count())
            .filter(|&edge| {
                let bond = &lattice.edges[edge];
                config[bond.source] == config[bond.target]
            })
            .collect();
        assert!(
            !eligible.is_empty(),
            "every test config has open candidates"
        );

        let mut adapter =
            EdwardsSokalMC::frozen(lattice.clone(), beta, config.to_vec(), plan.clone()).unwrap();
        let mut ctx = Context::new_with_binsize(
            Xoshiro256PlusPlus::seed_from_u64(SEED_BASE + config[0] as u64 + 2),
            0,
            draws + 1,
        );
        let mut counts = vec![0usize; 1 << eligible.len()];
        for _ in 0..draws {
            adapter.sweep(&mut ctx);
            let sample = adapter.last_bond_sample();
            let mut subset = 0usize;
            for (bit, &edge) in eligible.iter().enumerate() {
                if sample.contains(lattice.edge_id(edge).unwrap()) {
                    subset |= 1 << bit;
                }
            }
            counts[subset] += 1;
        }

        for (subset, &count) in counts.iter().enumerate() {
            let probability: f64 = eligible
                .iter()
                .enumerate()
                .map(|(bit, &edge)| {
                    if subset & (1 << bit) != 0 {
                        probabilities[edge]
                    } else {
                        1.0 - probabilities[edge]
                    }
                })
                .product();
            let expected = draws as f64 * probability;
            let deviation = (draws as f64 * probability * (1.0 - probability)).sqrt();
            let z = (count as f64 - expected) / deviation;
            assert!(
                z.abs() < 5.0,
                "frozen config {config:?}, subset {subset:0b}: count {count}, expected \
                 {expected:.1}, z = {z:.3} against the product Bernoulli law"
            );
        }

        // Exhaustive per-subset analysis parity: the adapter's measure-time
        // path (BondSample as edge activity, all-active vertices) and the R1
        // one-call raw-mask path agree on every subset of THIS configuration.
        let vertices = AllActiveVertices::new(lattice.vertex_count());
        let mut sample_workspace = ComponentWorkspace::new();
        let mut mask_workspace = ComponentWorkspace::new();
        for subset in 0..(1usize << eligible.len()) {
            let mut mask = vec![false; lattice.edge_count()];
            for (bit, &edge) in eligible.iter().enumerate() {
                mask[edge] = subset & (1 << bit) != 0;
            }
            let raw = {
                let result = analyze_bonds(&lattice, &mask, &queries, &mut mask_workspace)
                    .expect("full-length mask");
                (result.summary(), result.query_outcomes().to_vec())
            };
            let sample = BondSample::from_mask(&lattice, &mask).unwrap();
            let through_sample = {
                let result = analyze(
                    &lattice,
                    &vertices,
                    &sample,
                    &queries,
                    &mut sample_workspace,
                )
                .expect("sample dimensions match the graph");
                (result.summary(), result.query_outcomes().to_vec())
            };
            assert_eq!(
                raw, through_sample,
                "subset {subset:0b} of config {config:?}: BondSample path diverged from R1"
            );
        }
    }
}

// ── Gate 2: partition-function agreement ──────────────────────────────────
/// Number of connected components of `(V, A)` including isolates, via a
/// test-local union-find (independent of every production implementation).
fn component_count(vertex_count: usize, open_edges: &[usize], lattice: &CsrLattice) -> usize {
    fn root_of(parent: &mut [usize], site: usize) -> usize {
        let mut root = site;
        while parent[root] != root {
            root = parent[root];
        }
        let mut current = site;
        while parent[current] != root {
            let next = parent[current];
            parent[current] = root;
            current = next;
        }
        root
    }
    let mut parent: Vec<usize> = (0..vertex_count).collect();
    for &edge in open_edges {
        let bond = &lattice.edges[edge];
        let (left, right) = (
            root_of(&mut parent, bond.source),
            root_of(&mut parent, bond.target),
        );
        if left != right {
            parent[right] = left;
        }
    }
    (0..vertex_count)
        .filter(|&site| root_of(&mut parent, site) == site)
        .count()
}

#[test]
fn fk_double_sum_matches_exact_zero_field_ising_partition_function() {
    for beta in [0.3, 0.9] {
        for lattice in [wedge_graph(), build_chain(4, true), multigraph()] {
            let vertices = lattice.vertex_count();
            let edges = lattice.edge_count();
            let probabilities = bond_probabilities(&lattice, beta);
            let coupling_sum: f64 = lattice.edges.iter().map(|edge| edge.weight).sum();

            // The Edwards-Sokal double sum: every spin configuration, every
            // bond subset (incompatible pairs carry zero weight).
            let mut z_double = 0.0;
            for spin_mask in 0..(1usize << vertices) {
                let spins = spin_vector(vertices, spin_mask);
                for bond_mask in 0..(1usize << edges) {
                    let mut weight = 1.0;
                    for (edge, bond) in lattice.edges.iter().enumerate() {
                        let open = bond_mask & (1 << edge) != 0;
                        if open && spins[bond.source] != spins[bond.target] {
                            weight = 0.0;
                            break;
                        }
                        weight *= if open {
                            probabilities[edge]
                        } else {
                            1.0 - probabilities[edge]
                        };
                    }
                    z_double += weight;
                }
            }

            // The q = 2 random-cluster single sum with 2^{k(A)} and the
            // test-local union-find.
            let mut z_random_cluster = 0.0;
            for bond_mask in 0..(1usize << edges) {
                let open: Vec<usize> = (0..edges)
                    .filter(|&edge| bond_mask & (1 << edge) != 0)
                    .collect();
                let mut weight = 2f64.powi(component_count(vertices, &open, &lattice) as i32);
                for (edge, open_edge) in (0..edges).map(|edge| (edge, open.contains(&edge))) {
                    weight *= if open_edge {
                        probabilities[edge]
                    } else {
                        1.0 - probabilities[edge]
                    };
                }
                z_random_cluster += weight;
            }

            // Independent brute-force zero-field Ising reference.
            let mut z_ising = 0.0;
            for spin_mask in 0..(1usize << vertices) {
                let spins = spin_vector(vertices, spin_mask);
                let alignment: f64 = lattice
                    .edges
                    .iter()
                    .map(|bond| {
                        bond.weight * f64::from(spins[bond.source]) * f64::from(spins[bond.target])
                    })
                    .sum();
                z_ising += (beta * alignment).exp();
            }

            assert_close(z_double, z_random_cluster, 1e-12 * z_double.abs().max(1.0));
            assert_close(
                (beta * coupling_sum).exp() * z_double,
                z_ising,
                1e-12 * z_ising.abs().max(1.0),
            );
        }
    }
}

// ── Gate 3: P_FK(i ~ j) = <s_i s_j> on a proper ES chain ──────────────────
/// Exact `<s_i s_j>` by brute-force spin enumeration of `exp(-beta H)`.
fn exact_spin_correlation(lattice: &CsrLattice, beta: f64, left: usize, right: usize) -> f64 {
    let vertices = lattice.vertex_count();
    let (mut numerator, mut partition) = (0.0, 0.0);
    for spin_mask in 0..(1usize << vertices) {
        let spins = spin_vector(vertices, spin_mask);
        let alignment: f64 = lattice
            .edges
            .iter()
            .map(|bond| bond.weight * f64::from(spins[bond.source]) * f64::from(spins[bond.target]))
            .sum();
        let weight = (beta * alignment).exp();
        partition += weight;
        numerator += weight * f64::from(spins[left]) * f64::from(spins[right]);
    }
    numerator / partition
}

#[test]
fn fk_connection_probability_equals_spin_correlation() {
    let lattice = build_chain(4, true);
    let beta = 0.9;
    let pairs = [(0usize, 1usize), (0usize, 2usize)];
    let exact: Vec<f64> = pairs
        .iter()
        .map(|&(left, right)| exact_spin_correlation(&lattice, beta, left, right))
        .collect();
    for &value in &exact {
        assert!(
            (0.1..0.99).contains(&value),
            "exact correlator {value} is a nontrivial reference"
        );
    }

    let draws = draw_count(60_000);
    let queries = pairs
        .iter()
        .map(|&(left, right)| BoundaryQuery::new(&lattice, &[left], &[right]).unwrap())
        .collect();
    let plan = ObservablePlan::new([StaticObservable::ActiveEdgeCount], queries);
    let mut adapter = EdwardsSokalMC::new(lattice.clone(), beta, plan).unwrap();
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(SEED_BASE), 0, 1);

    // Thermalize the ES chain before measuring (the static-family
    // zero-thermalization convention does NOT apply: this is a Markov chain).
    for _ in 0..1_000 {
        adapter.sweep(&mut ctx);
    }
    let mut correlations = vec![0.0f64; pairs.len()];
    for _ in 0..draws {
        adapter.sweep(&mut ctx);
        adapter.measure(&mut ctx);
        let spins = adapter.spins();
        for (index, &(left, right)) in pairs.iter().enumerate() {
            correlations[index] += f64::from(spins[left]) * f64::from(spins[right]);
        }
    }

    // P_FK comes from the boundary-crossing means recorded on the SAME
    // exposed samples; the spin correlator is read from the chain's spins.
    let values: HashMap<String, carlo_rs::Estimate> = ctx.finalize_measurements();
    let standard_error =
        |probability: f64| (probability * (1.0 - probability) / draws as f64).sqrt();
    for (index, &(left, right)) in pairs.iter().enumerate() {
        let p_fk = values[&format!("BoundaryCrossing{index}")].mean;
        let spin_mean = correlations[index] / draws as f64;
        let reference = exact[index];
        let z_fk = (p_fk - reference) / standard_error(reference);
        let z_spin = (spin_mean - reference) / standard_error(reference);
        let z_gap = (p_fk - spin_mean) / ((1.0 - reference) / draws as f64).sqrt();
        assert!(
            z_fk.abs() < 5.0,
            "P_FK({left} ~ {right}) = {p_fk:.5} vs exact <s s> = {reference:.5}, z = {z_fk:.3}"
        );
        assert!(
            z_spin.abs() < 5.0,
            "<s_{left} s_{right}> = {spin_mean:.5} vs exact {reference:.5}, z = {z_spin:.3}"
        );
        assert!(
            z_gap.abs() < 6.0,
            "P_FK = {p_fk:.5} and <s s> = {spin_mean:.5} disagree on pair ({left}, {right}), \
             z = {z_gap:.3} — the Edwards-Sokal identity failed within the chain"
        );
    }
}

// ── Gate 4: worm-representation consistency (exact, independent) ─────────
/// Degree parity of every vertex over the edge subset `bond_mask`.
fn parity_defects(lattice: &CsrLattice, bond_mask: usize) -> Vec<bool> {
    let mut odd = vec![false; lattice.vertex_count()];
    for (edge, bond) in lattice.edges.iter().enumerate() {
        if bond_mask & (1 << edge) != 0 {
            odd[bond.source] = !odd[bond.source];
            odd[bond.target] = !odd[bond.target];
        }
    }
    odd
}

/// The worm's high-temperature representation computes Ising observables from
/// even-subgraph sums with `tanh(beta J_e)` weights (the physical sector) and
/// two-defect odd subgraphs (the correlator sector). This gate checks those
/// sums against brute-force spin enumeration — the level at which the worm
/// and FK representations must agree. Configuration-level agreement is
/// impossible by construction (see the module note above).
#[test]
fn worm_representation_correlations_agree_with_spin_enumeration() {
    let lattice = build_chain(4, true);
    let edges = lattice.edge_count();
    for beta in [0.45, 0.9] {
        let tanh_weights: Vec<f64> = lattice
            .edges
            .iter()
            .map(|edge| (beta * edge.weight).tanh())
            .collect();

        let mut even_sum = 0.0;
        let mut defect_sums = [0.0f64; 2];
        for bond_mask in 0..(1usize << edges) {
            let weight: f64 = (0..edges)
                .map(|edge| {
                    if bond_mask & (1 << edge) != 0 {
                        tanh_weights[edge]
                    } else {
                        1.0
                    }
                })
                .product();
            let odd = parity_defects(&lattice, bond_mask);
            if odd.iter().all(|&defect| !defect) {
                even_sum += weight;
            }
            // The correlator sector is the EXACT defect pattern: odd at the
            // two correlator endpoints, even everywhere else.
            let two_defect = |left: usize, right: usize| {
                odd[left]
                    && odd[right]
                    && odd
                        .iter()
                        .enumerate()
                        .all(|(site, &defect)| defect == (site == left || site == right))
            };
            if two_defect(0, 1) {
                defect_sums[0] += weight;
            }
            if two_defect(0, 2) {
                defect_sums[1] += weight;
            }
        }

        // Z = 2^V prod_e cosh(beta J_e) * sum_even (the worm normalization).
        let cosh_product: f64 = lattice
            .edges
            .iter()
            .map(|edge| (beta * edge.weight).cosh())
            .product();
        let z_worm = 2f64.powi(lattice.vertex_count() as i32) * cosh_product * even_sum;
        let mut z_ising = 0.0;
        for spin_mask in 0..(1usize << lattice.vertex_count()) {
            let spins = spin_vector(lattice.vertex_count(), spin_mask);
            let alignment: f64 = lattice
                .edges
                .iter()
                .map(|bond| {
                    bond.weight * f64::from(spins[bond.source]) * f64::from(spins[bond.target])
                })
                .sum();
            z_ising += (beta * alignment).exp();
        }
        assert_close(z_worm, z_ising, 1e-12 * z_ising);

        for (index, &(left, right)) in [(0usize, 1usize), (0usize, 2usize)].iter().enumerate() {
            let worm_correlation = defect_sums[index] / even_sum;
            let spin_correlation = exact_spin_correlation(&lattice, beta, left, right);
            assert_close(
                worm_correlation,
                spin_correlation,
                1e-12 * spin_correlation.abs().max(1.0),
            );
        }
    }
}
