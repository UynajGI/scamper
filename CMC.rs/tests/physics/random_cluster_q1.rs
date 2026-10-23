//! R2 q = 1 random-cluster reduction validation (roadmap §15 R2).
//!
//! Exact statement under test: **the q = 1 FK/random-cluster measure on a
//! graph `G = (V, E)` equals the product Bernoulli(p) measure on its
//! bonds**, configuration by configuration and statistically.
//!
//! # Hand derivation of the q = 1 collapse
//!
//! The random-cluster measure with edge parameter `p` and cluster weight
//! `q` assigns each bond subset `A ⊆ E` (all vertices active — the F3 bond
//! convention, `k(A)` = component count of `(V, A)`) the weight
//!
//! ```text
//! w(A)      = p^{|A|} (1-p)^{E-|A|} q^{k(A)},
//! phi(A)    = w(A) / Z(p, q),
//! Z(p, q)   = sum_{A ⊆ E} w(A).
//! ```
//!
//! At `q = 1` the cluster-count factor is exactly `1^{k(A)} = 1` for every
//! configuration — `k(A)` drops out of the weight entirely, whatever the
//! graph (simple, parallel edges, self-loops; the self-loop never changes
//! `k(A)` and a parallel pair never changes it either — they only add
//! `p`-factors). Therefore
//!
//! ```text
//! Z(p, 1) = sum_A p^{|A|} (1-p)^{E-|A|} = (p + (1-p))^E = 1,
//! phi_{p,1}(A) = p^{|A|} (1-p)^{E-|A|}
//!              = product over e of [ p   if e ∈ A,
//!                                    (1-p) if e ∉ A ],
//! ```
//!
//! i.e. the independent bond Bernoulli product measure (heterogeneous
//! per-edge `p_e` analogously: the product of per-edge factors). This is
//! exact mathematics — the founding q = 1 observation of Fortuin and
//! Kasteleyn (Physica 57, 536 (1972); pre-arXiv, nothing archived) — and
//! the gates below prove it by computation rather than fit it.
//!
//! Gates:
//!
//! 1. **Configuration-by-configuration (roadmap gate)**: exhaustive
//!    enumeration of ALL `2^E` bond configurations on a 5-edge chain-ish
//!    graph and a 6-edge multigraph (self-loop + parallel pair). For each
//!    configuration the enumerated FK weight `p^{|A|} (1-p)^{E-|A|} q^{k(A)}`
//!    — with the `q = 1` cluster factor computed via `1.0f64.powi(k(A))`
//!    and multiplied in — (i) sums with its siblings to `Z(p,1) = 1` and
//!    (ii) equals the Bernoulli product weight computed edge by edge; the
//!    two full distributions agree term by term, and every F3 observable
//!    mean (component count, largest size) computed under either
//!    distribution is identical. `k(A)` itself comes from an independent
//!    flood fill cross-checked against the production analyzer.
//! 2. **Bitwise sampler identity (delegation contract)**: same seed → the
//!    adapter's `BondSample`-to-mask equals `BondBernoulli`'s mask read
//!    back from a `StaticConfiguration`, uniform and heterogeneous, at
//!    endpoint and interior `p`.
//! 3. **Statistical (roadmap gate)**: multi-seed two-sample z-scores of
//!    GiantFraction and finite-cluster susceptibility, RandomClusterMC
//!    (Bernoulli delegate → R1 external-mask analysis) versus the
//!    `StaticPercolationMC` bond path (`StaticConfiguration` analysis) at
//!    matched `p` on an open square lattice, independent RNG streams per
//!    path: per-seed `|z| < 4` and `|mean z| < 2`.
//! 4. **q ≠ 1**: typed rejection citing the roadmap's separate-plan
//!    requirement.

use std::collections::HashMap;

use super::common::zscore_seed_count;
use carlo_rs::{Context, Estimate, MonteCarlo};
use cmc_rs::{
    analyze_bonds, build_chain, build_square, Bond, BondBernoulli, BondType, ClusterQ,
    ComponentWorkspace, CsrLattice, EdgeActivity, GraphView, ObservablePlan, Probability,
    ProbabilityField, RandomClusterError, RandomClusterMC, StaticConfiguration, StaticLaw,
    StaticObservable, StaticPercolationMC, UndirectedGraphView,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Independent reference for `k(A)` and `S_max(A)`: adjacency flood fill
/// over the mask-open edges with every vertex active. Shares no code with
/// the production union-find analyzer.
fn independent_components_and_largest<G: UndirectedGraphView>(
    graph: &G,
    mask: &[bool],
) -> (usize, usize) {
    let mut seen = vec![false; graph.vertex_count()];
    let mut components = 0usize;
    let mut largest = 0usize;
    for seed in graph.vertex_ids() {
        if seen[seed.index()] {
            continue;
        }
        seen[seed.index()] = true;
        let mut stack = vec![seed];
        let mut size = 0usize;
        while let Some(vertex) = stack.pop() {
            size += 1;
            for incidence in graph.incidences(vertex) {
                // An open self-loop sees its own endpoint already visited;
                // a parallel pair offers the same neighbor twice. Both are
                // no-ops here exactly as in the F3 union-find.
                if mask[incidence.edge.index()] && !seen[incidence.neighbor.index()] {
                    seen[incidence.neighbor.index()] = true;
                    stack.push(incidence.neighbor);
                }
            }
        }
        components += 1;
        largest = largest.max(size);
    }
    (components, largest)
}

/// Bernoulli product weight of one configuration: the per-edge factor
/// product, computed independently of any power formula.
fn bernoulli_product_weight(probabilities: &[f64], mask: &[bool]) -> f64 {
    probabilities
        .iter()
        .zip(mask)
        .map(|(&probability, &open)| if open { probability } else { 1.0 - probability })
        .product()
}

/// Enumerate every one of the `2^E` bond configurations and prove the
/// q = 1 identity at `probabilities` (one per physical edge).
fn assert_q1_full_distribution_identity<G: UndirectedGraphView>(
    graph: &G,
    probabilities: &[f64],
    context: &str,
) {
    let edge_count = graph.edge_count();
    assert_eq!(probabilities.len(), edge_count, "{context}: field length");
    let uniform = probabilities
        .iter()
        .copied()
        .all(|value| value == probabilities[0]);
    let q: f64 = 1.0; // the R2 cluster weight, kept explicit in the formula below
    let mut workspace = ComponentWorkspace::new();

    let mut fk_total = 0.0;
    let mut bernoulli_total = 0.0;
    let mut fk_mean_components = 0.0;
    let mut bernoulli_mean_components = 0.0;
    let mut fk_mean_largest = 0.0;
    let mut bernoulli_mean_largest = 0.0;
    for bits in 0..(1usize << edge_count) {
        let mask: Vec<bool> = (0..edge_count)
            .map(|index| bits & (1 << index) != 0)
            .collect();
        let (components, largest) = independent_components_and_largest(graph, &mask);
        // Production F3 (through the R1 external path) agrees with the
        // independent flood fill on every single enumerated configuration.
        let summary = analyze_bonds(graph, &mask, &[], &mut workspace)
            .unwrap_or_else(|error| panic!("{context}/mask={bits:?}: analysis failed: {error}"))
            .summary();
        assert_eq!(
            summary.component_count, components,
            "{context}/mask={bits:?}: component count"
        );
        assert_eq!(
            summary.largest_component_size, largest,
            "{context}/mask={bits:?}: largest component"
        );

        // The Bernoulli product weight: per-edge factors.
        let bernoulli = bernoulli_product_weight(probabilities, &mask);
        // The FK weight: the canonical random-cluster formula with the
        // cluster factor computed (as exactly 1.0) and multiplied in. For
        // uniform p this takes the power form p^{|A|} (1-p)^{E-|A|}, a
        // genuinely different floating-point path from the product above.
        let cluster_factor = q.powi(components as i32);
        let open = mask.iter().filter(|&&open| open).count();
        let fk = if uniform {
            let p = probabilities[0];
            p.powi(open as i32) * (1.0 - p).powi((edge_count - open) as i32) * cluster_factor
        } else {
            bernoulli * cluster_factor
        };

        // (ii) per-configuration probability equality: FK weight equals the
        // Bernoulli product weight exactly up to power-vs-product floating
        // association (E <= 6 keeps the discrepancy far below 1e-12).
        let deviation = (fk - bernoulli).abs();
        assert!(
            deviation <= 1e-12,
            "{context}/mask={bits:?}: FK weight {fk:.17e} != Bernoulli {bernoulli:.17e} \
             (deviation {deviation:.3e})"
        );
        assert!((0.0..=1.0).contains(&bernoulli), "{context}: weight range");
        fk_total += fk;
        bernoulli_total += bernoulli;
        fk_mean_components += fk * components as f64;
        bernoulli_mean_components += bernoulli * components as f64;
        fk_mean_largest += fk * largest as f64;
        bernoulli_mean_largest += bernoulli * largest as f64;
    }

    // (i) the enumerated FK probabilities are a distribution: Z(p, 1) = 1.
    assert!(
        (fk_total - 1.0).abs() <= 1e-9,
        "{context}: sum of q = 1 FK weights is {fk_total:.17e}, not 1"
    );
    assert!(
        (bernoulli_total - 1.0).abs() <= 1e-9,
        "{context}: sum of Bernoulli product weights is {bernoulli_total:.17e}, not 1"
    );
    // The two full distributions are identical term by term, hence every
    // observable mean coincides — checked, not implied, for the two F3
    // scalars the statistical gate below samples stochastically.
    assert!(
        (fk_mean_components - bernoulli_mean_components).abs() <= 1e-9,
        "{context}: <k> under FK {fk_mean_components:.17e} vs Bernoulli \
         {bernoulli_mean_components:.17e}"
    );
    assert!(
        (fk_mean_largest - bernoulli_mean_largest).abs() <= 1e-9,
        "{context}: <S_max> under FK {fk_mean_largest:.17e} vs Bernoulli \
         {bernoulli_mean_largest:.17e}"
    );
}

/// 5-edge "chain-ish" graph: a square 0-1-2-3 with a pendant tail 0-4
/// (V = 5, E = 5, so `k(A)` is not a function of `|A|` alone — the cycle
/// makes component structure depend on *which* edges are open).
fn chainish_graph() -> CsrLattice {
    CsrLattice::from_edges(
        5,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
            Bond::new(2, 3, BondType::Generic, 1.0),
            Bond::new(3, 0, BondType::Generic, 1.0),
            Bond::new(0, 4, BondType::Generic, 1.0),
        ],
    )
}

/// 6-edge multigraph (V = 4): one self-loop, one parallel pair, and the
/// 4-cycle 0-1-2-3. Parallel edges and the self-loop are distinct physical
/// edges with independent mask entries and independent `p` factors; none of
/// them can change `k(A)`, so the `q = 1` cluster factor stays exactly 1.
fn multigraph() -> CsrLattice {
    CsrLattice::from_edges(
        4,
        vec![
            Bond::new(0, 0, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
            Bond::new(2, 3, BondType::Generic, 1.0),
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(0, 3, BondType::Generic, 1.0),
        ],
    )
}

#[test]
fn q1_fk_distribution_equals_bernoulli_product_on_chainish_graph() {
    let graph = chainish_graph();
    assert_eq!(graph.edge_count(), 5);
    for p in [0.0f64, 0.25, 0.5, 0.75, 1.0] {
        assert_q1_full_distribution_identity(&graph, &[p; 5], &format!("chainish/uniform p = {p}"));
    }
    assert_q1_full_distribution_identity(
        &graph,
        &[0.1, 0.85, 0.3, 0.6, 0.45],
        "chainish/heterogeneous",
    );
}

#[test]
fn q1_fk_distribution_equals_bernoulli_product_on_multigraph() {
    let graph = multigraph();
    assert_eq!(graph.edge_count(), 6);
    for p in [0.0f64, 0.25, 0.5, 0.75, 1.0] {
        assert_q1_full_distribution_identity(
            &graph,
            &[p; 6],
            &format!("multigraph/uniform p = {p}"),
        );
    }
    assert_q1_full_distribution_identity(
        &graph,
        &[0.05, 0.9, 0.25, 0.7, 0.5, 0.35],
        "multigraph/heterogeneous",
    );
}

#[test]
fn q1_adapter_draw_is_bitwise_bond_bernoulli() {
    let lattice = build_square(6, 6, false);
    let edge_count = lattice.edge_count();
    let plan = ObservablePlan::new([StaticObservable::GiantFraction], Vec::new());

    // Uniform law: the adapter delegates to BondBernoulli, so an
    // identically seeded direct draw must produce the identical mask —
    // through the adapter's dense view AND through R1's BondSample record.
    for p in [0.0f64, 0.25, 0.5, 0.73, 1.0] {
        for seed in [1u64, 42, 90210] {
            let mut adapter = RandomClusterMC::new(
                lattice.clone(),
                ClusterQ::q1(),
                p,
                ObservablePlan::new([StaticObservable::GiantFraction], Vec::new()),
            )
            .expect("valid adapter");
            let law = BondBernoulli::new(Probability::new(p).expect("test p").into());
            let mut reference = StaticConfiguration::new(lattice.vertex_count(), edge_count);
            let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, 1);
            adapter.sweep(&mut ctx);
            law.sample(
                &lattice,
                &mut reference,
                &mut Xoshiro256PlusPlus::seed_from_u64(seed),
            )
            .expect("matching dimensions");
            let expected: Vec<bool> = (0..edge_count)
                .map(|index| reference.edge_active(lattice.edge_id(index).expect("dense edge")))
                .collect();
            assert_eq!(
                adapter.mask(),
                expected.as_slice(),
                "p = {p}, seed = {seed}"
            );
            assert_eq!(
                adapter.realization().to_mask(),
                expected,
                "BondSample mask view, p = {p}, seed = {seed}"
            );
        }
    }

    // Heterogeneous law: the owned per-edge field delegates identically to
    // the borrowed-field reference law.
    let values: Vec<Probability> = [
        0.05, 0.95, 0.3, 0.6, 0.45, 0.8, 0.15, 0.55, 0.25, 0.7, 0.4, 0.9,
    ]
    .iter()
    .chain(std::iter::repeat(&0.5))
    .take(edge_count)
    .map(|&value| Probability::new(value).expect("test probability"))
    .collect();
    let mut adapter =
        RandomClusterMC::new_heterogeneous(lattice.clone(), ClusterQ::q1(), values.clone(), plan)
            .expect("field length matches the topology");
    let law = BondBernoulli::new(ProbabilityField::Borrowed(&values));
    let mut reference = StaticConfiguration::new(lattice.vertex_count(), edge_count);
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(7), 0, 1);
    adapter.sweep(&mut ctx);
    law.sample(
        &lattice,
        &mut reference,
        &mut Xoshiro256PlusPlus::seed_from_u64(7),
    )
    .expect("matching dimensions");
    let expected: Vec<bool> = (0..edge_count)
        .map(|index| reference.edge_active(lattice.edge_id(index).expect("dense edge")))
        .collect();
    assert_eq!(adapter.mask(), expected.as_slice(), "heterogeneous");
    assert_eq!(
        adapter.realization().to_mask(),
        expected,
        "heterogeneous sample"
    );
}

/// Side of the open square lattice for the statistical gate.
const SIDE: usize = 8;
/// Matched bond probability: the square-lattice bond threshold, where
/// cluster observables fluctuate hardest and the gate has most power.
const P_MATCHED: f64 = 0.5;
/// Sweeps (one i.i.d. realization each) per seed and path.
const SWEEPS: u64 = 10_000;
/// Samples per statistics bin.
const BINSIZE: usize = 100;
/// Default seed count (override: `SCUTTLE_ZSCORE_SEEDS`).
const DEFAULT_SEEDS: usize = 16;
/// Base for the per-path, per-seed RNG streams (even = RandomClusterMC,
/// odd = StaticPercolationMC — disjoint streams for both paths).
const SEED_BASE: u64 = 20261022;

fn run_realizations<M>(mut mc: M, seed: u64) -> HashMap<String, Estimate>
where
    M: MonteCarlo<Rng = Xoshiro256PlusPlus>,
{
    let mut ctx = Context::new_with_binsize(Xoshiro256PlusPlus::seed_from_u64(seed), 0, BINSIZE);
    for _ in 0..SWEEPS {
        mc.sweep(&mut ctx);
        mc.measure(&mut ctx);
    }
    ctx.finalize_measurements()
}

#[test]
fn q1_matches_static_bond_path_statistically() {
    let n_seeds = zscore_seed_count(DEFAULT_SEEDS);
    let lattice = build_square(SIDE, SIDE, false);
    let plan = || {
        ObservablePlan::new(
            [
                StaticObservable::GiantFraction,
                StaticObservable::FiniteClusterSusceptibility,
            ],
            Vec::new(),
        )
    };

    for observable in ["GiantFraction", "FiniteClusterSusceptibility"] {
        let mut z_scores = Vec::with_capacity(n_seeds);
        for seed_index in 0..n_seeds as u64 {
            // Two independent code paths with independent RNG streams:
            // RandomClusterMC draws through the BondBernoulli delegate and
            // analyzes the R1 external mask; StaticPercolationMC draws its
            // StaticLaw::bond configuration and analyzes the
            // StaticConfiguration. Matched p, matched topology, matched plan.
            let cluster = run_realizations(
                RandomClusterMC::new(lattice.clone(), ClusterQ::q1(), P_MATCHED, plan())
                    .expect("valid adapter"),
                SEED_BASE + 2 * seed_index,
            );
            let reference = run_realizations(
                StaticPercolationMC::new(
                    lattice.clone(),
                    StaticLaw::bond(P_MATCHED).expect("matched bond law"),
                    plan(),
                )
                .expect("valid static adapter"),
                SEED_BASE + 2 * seed_index + 1,
            );
            let estimate = |values: &HashMap<String, Estimate>| {
                let estimate = values
                    .get(observable)
                    .unwrap_or_else(|| panic!("{observable} missing from {values:?}"));
                (estimate.mean, estimate.stderr)
            };
            let (cluster_mean, cluster_stderr) = estimate(&cluster);
            let (reference_mean, reference_stderr) = estimate(&reference);
            // Two-sample z with independent errors.
            let combined_stderr = (cluster_stderr * cluster_stderr
                + reference_stderr * reference_stderr)
                .sqrt()
                .max(1e-12);
            z_scores.push((cluster_mean - reference_mean) / combined_stderr);
        }

        let max_abs_z = z_scores.iter().fold(0.0f64, |acc, z| acc.max(z.abs()));
        let mean_z = z_scores.iter().sum::<f64>() / z_scores.len() as f64;
        assert!(
            max_abs_z < 4.0,
            "{observable}: max |z| = {max_abs_z:.2} over {n_seeds} seeds should be < 4"
        );
        assert!(
            mean_z.abs() < 2.0,
            "{observable}: mean z = {mean_z:.2} over {n_seeds} seeds should be |mean z| < 2"
        );
    }
}

#[test]
fn q_neq_one_is_a_typed_roadmap_cited_rejection() {
    for value in [0.0f64, 0.5, 1.0 + f64::EPSILON, 2.0, -3.0, f64::NAN] {
        let error = ClusterQ::new(value)
            .err()
            .unwrap_or_else(|| panic!("q = {value} must be rejected"));
        assert!(
            matches!(error, RandomClusterError::UnsupportedQ { .. }),
            "wrong error for q = {value}: {error}"
        );
        let message = error.to_string();
        assert!(message.contains("exactly 1"), "q = {value}: {message}");
        assert!(message.contains("roadmap"), "q = {value}: {message}");
        assert!(message.contains("separate plan"), "q = {value}: {message}");
    }
    // The adapter surface accepts the one legal cluster weight.
    let lattice = build_chain(4, false);
    let plan = ObservablePlan::new([StaticObservable::LargestSize], Vec::new());
    assert!(RandomClusterMC::new(lattice, ClusterQ::q1(), 0.5, plan).is_ok());
}
