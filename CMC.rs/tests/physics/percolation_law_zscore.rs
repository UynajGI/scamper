//! Multi-seed aggregate statistical validation of the F2 Bernoulli laws.

use super::common::zscore_seed_count;
use cmc_rs::{
    build_square, BondBernoulli, EdgeActivity, GraphView, MixedBernoulli, Probability,
    ProbabilityField, SiteBernoulli, StaticConfiguration, UndirectedGraphView, VertexActivity,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const N_SEEDS: usize = 16;
const SAMPLES_PER_SEED: usize = 512;
const SIDE: usize = 8;

// Six independent domain-level gates are evaluated by this test family. Under
// the normal approximation, 6 * 2 * Phi(-4.5) is about 4.1e-5 family-wise.
const AGGREGATE_Z_LIMIT: f64 = 4.5;

fn probability(value: f64) -> Probability {
    Probability::new(value).expect("test probability")
}

fn active_counts(
    graph: &impl UndirectedGraphView,
    configuration: &StaticConfiguration,
) -> (usize, usize) {
    let vertices = graph
        .vertex_ids()
        .filter(|&vertex| configuration.vertex_active(vertex))
        .count();
    let edges = (0..graph.edge_count())
        .filter(|&index| configuration.edge_active(graph.edge_id(index).expect("dense edge ID")))
        .count();
    (vertices, edges)
}

fn moments(probabilities: &[Probability]) -> (f64, f64) {
    probabilities
        .iter()
        .fold((0.0, 0.0), |(mean, variance), p| {
            let p = p.get();
            (mean + p, variance + p * (1.0 - p))
        })
}

fn assert_aggregate_gate(
    successes: usize,
    seeds: usize,
    samples_per_seed: usize,
    one_sample_mean: f64,
    one_sample_variance: f64,
    label: &str,
) {
    let repetitions = (seeds * samples_per_seed) as f64;
    let mean = repetitions * one_sample_mean;
    let variance = repetitions * one_sample_variance;
    if variance == 0.0 {
        assert_eq!(
            successes as f64, mean,
            "{label}: zero-variance law must match its exact activity count"
        );
        return;
    }

    let z = (successes as f64 - mean) / variance.sqrt();
    eprintln!(
        "[percolation-law-zscore] {label}: seeds={seeds}, successes={successes}, mean={mean:.3}, variance={variance:.3}, z={z:.3}"
    );
    assert!(
        z.abs() < AGGREGATE_Z_LIMIT,
        "{label}: aggregate |z| = {:.3}, expected < {AGGREGATE_Z_LIMIT}",
        z.abs()
    );
}

#[test]
fn percolation_law_zscore_uniform_site_and_bond() {
    let graph = build_square(SIDE, SIDE, false);
    let seeds = zscore_seed_count(N_SEEDS);
    let site_p = 0.37;
    let bond_p = 0.61;
    let site = SiteBernoulli::new(probability(site_p).into());
    let bond = BondBernoulli::new(probability(bond_p).into());
    let (mut site_successes, mut bond_successes) = (0, 0);

    for seed in 0..seeds as u64 {
        let mut configuration = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x5349_5445_0000 + seed);
        for _ in 0..SAMPLES_PER_SEED {
            site.sample(&graph, &mut configuration, &mut rng)
                .expect("matching domains");
            site_successes += active_counts(&graph, &configuration).0;
        }

        let mut configuration = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x424f_4e44_0000 + seed);
        for _ in 0..SAMPLES_PER_SEED {
            bond.sample(&graph, &mut configuration, &mut rng)
                .expect("matching domains");
            bond_successes += active_counts(&graph, &configuration).1;
        }
    }

    assert_aggregate_gate(
        site_successes,
        seeds,
        SAMPLES_PER_SEED,
        graph.vertex_count() as f64 * site_p,
        graph.vertex_count() as f64 * site_p * (1.0 - site_p),
        "uniform site activity",
    );
    assert_aggregate_gate(
        bond_successes,
        seeds,
        SAMPLES_PER_SEED,
        graph.edge_count() as f64 * bond_p,
        graph.edge_count() as f64 * bond_p * (1.0 - bond_p),
        "uniform bond activity",
    );
}

#[test]
fn percolation_law_zscore_uniform_mixed_tracks_domains_separately() {
    let graph = build_square(SIDE, SIDE, false);
    let seeds = zscore_seed_count(N_SEEDS);
    let vertex_p = 0.43;
    let edge_p = 0.68;
    let law = MixedBernoulli::new(probability(vertex_p).into(), probability(edge_p).into());
    let (mut vertex_successes, mut edge_successes) = (0, 0);

    for seed in 0..seeds as u64 {
        let mut configuration = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4d49_5845_4400 + seed);
        for _ in 0..SAMPLES_PER_SEED {
            law.sample(&graph, &mut configuration, &mut rng)
                .expect("matching domains");
            let (vertices, edges) = active_counts(&graph, &configuration);
            vertex_successes += vertices;
            edge_successes += edges;
        }
    }

    assert_aggregate_gate(
        vertex_successes,
        seeds,
        SAMPLES_PER_SEED,
        graph.vertex_count() as f64 * vertex_p,
        graph.vertex_count() as f64 * vertex_p * (1.0 - vertex_p),
        "uniform mixed vertex activity",
    );
    assert_aggregate_gate(
        edge_successes,
        seeds,
        SAMPLES_PER_SEED,
        graph.edge_count() as f64 * edge_p,
        graph.edge_count() as f64 * edge_p * (1.0 - edge_p),
        "uniform mixed edge activity",
    );
}

#[test]
fn percolation_law_zscore_heterogeneous_mixed_poisson_binomial() {
    let graph = build_square(SIDE, SIDE, false);
    let seeds = zscore_seed_count(N_SEEDS);
    let vertex_probabilities = (0..graph.vertex_count())
        .map(|index| probability([0.15, 0.35, 0.65, 0.85][index % 4]))
        .collect::<Vec<_>>();
    let edge_probabilities = (0..graph.edge_count())
        .map(|index| probability([0.2, 0.45, 0.7][index % 3]))
        .collect::<Vec<_>>();
    let law = MixedBernoulli::new(
        ProbabilityField::Borrowed(&vertex_probabilities),
        ProbabilityField::Owned(edge_probabilities.clone()),
    );
    let (mut vertex_successes, mut edge_successes) = (0, 0);

    for seed in 0..seeds as u64 {
        let mut configuration = StaticConfiguration::new(graph.vertex_count(), graph.edge_count());
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4845_5445_524f + seed);
        for _ in 0..SAMPLES_PER_SEED {
            law.sample(&graph, &mut configuration, &mut rng)
                .expect("matching domains");
            let (vertices, edges) = active_counts(&graph, &configuration);
            vertex_successes += vertices;
            edge_successes += edges;
        }
    }

    let (vertex_mean, vertex_variance) = moments(&vertex_probabilities);
    let (edge_mean, edge_variance) = moments(&edge_probabilities);
    assert_aggregate_gate(
        vertex_successes,
        seeds,
        SAMPLES_PER_SEED,
        vertex_mean,
        vertex_variance,
        "heterogeneous mixed vertex activity",
    );
    assert_aggregate_gate(
        edge_successes,
        seeds,
        SAMPLES_PER_SEED,
        edge_mean,
        edge_variance,
        "heterogeneous mixed edge activity",
    );
}
