//! X4: contact-process threshold bracket on a random regular graph.
//!
//! The contact process is the absorbing-state directed-percolation model
//! (Hinrichsen 2000, archived). On a locally tree-like `k`-regular graph the
//! branching (mean-field) threshold is `lambda_c = 1/(k-1)`, i.e. 0.5 for
//! `k = 3`. This is a finite-size, finite-time bracket: subcritical runs
//! must die out, supercritical runs must retain positive density, and the
//! density crossing must land near 0.5.

use cmc_rs::{
    contact_ensemble, Bond, BondType, ContactProcessModel, CsrLattice, GraphView,
    UndirectedGraphView,
};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

/// Configuration-model random `degree`-regular multigraph.
fn regular_graph(vertices: usize, degree: usize, seed: u64) -> CsrLattice {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut stubs = Vec::with_capacity(vertices * degree);
    for vertex in 0..vertices {
        stubs.extend(std::iter::repeat_n(vertex, degree));
    }
    for index in (1..stubs.len()).rev() {
        let swap = rng.random_range(0..=index);
        stubs.swap(index, swap);
    }
    let bonds = stubs
        .as_chunks::<2>()
        .0
        .iter()
        .map(|pair| Bond::new(pair[0], pair[1], BondType::Generic, 1.0))
        .collect::<Vec<_>>();
    CsrLattice::from_edges(vertices, bonds)
}

#[test]
fn contact_process_brackets_the_branching_threshold() {
    let vertices = 200;
    let degree = 3;
    let duration = 15.0;
    let initial_density = 0.05;
    let replicas = 4;
    let graphs = 2;
    let mut rows = Vec::new();
    let mut lambda: f64 = 0.2;
    while lambda <= 1.0 + 1e-12 {
        let mut density_sum = 0.0;
        let mut survived = 0usize;
        let mut total = 0usize;
        for graph_index in 0..graphs {
            let graph = regular_graph(
                vertices,
                degree,
                0x4138_0001 ^ lambda.to_bits() ^ ((graph_index as u64) << 32),
            );
            let model = ContactProcessModel::new(graph.clone(), lambda, 1.0).expect("model");
            assert_eq!(graph.vertex_count(), vertices);
            assert_eq!(graph.edge_count(), vertices * degree / 2);
            let ensemble = contact_ensemble(
                &model,
                initial_density,
                replicas,
                duration,
                0x4138_0002 ^ lambda.to_bits() ^ ((graph_index as u64) << 48),
            )
            .expect("ensemble");
            density_sum += ensemble.mean_density * replicas as f64;
            survived += ensemble.survived;
            total += replicas;
        }
        let density = density_sum / total as f64;
        let survival = survived as f64 / total as f64;
        eprintln!("X4 lambda={lambda:.2}: survival {survival:.3}, mean density {density:.4}");
        rows.push((lambda, density));
        lambda += 0.1;
    }

    let (_, first) = rows[0];
    assert!(
        first < 0.01,
        "subcritical lambda=0.2 has density {first:.4}"
    );
    let (_, last) = *rows.last().expect("non-empty ladder");
    assert!(last > 0.1, "supercritical lambda=1.0 has density {last:.4}");

    let crossing = rows
        .windows(2)
        .find_map(|pair| {
            let (l0, d0) = pair[0];
            let (l1, d1) = pair[1];
            (d0 < 0.05 && d1 >= 0.05).then(|| l0 + (0.05 - d0) / (d1 - d0) * (l1 - l0))
        })
        .expect("density crosses 0.05 inside the ladder");
    eprintln!("X4 0.05-density crossing: lambda={crossing:.3} (branching threshold 0.5)");
    assert!(
        (crossing - 0.5).abs() < 0.25,
        "crossing {crossing:.3} is outside the 0.25 bracket around 0.5"
    );
    assert!(
        crossing > 0.3 && crossing < 0.7,
        "crossing {crossing:.3} is outside the finite-size window"
    );
}

#[test]
fn absorbing_runs_report_zero_density() {
    let graph = regular_graph(60, 3, 0x4138_0003);
    let model = ContactProcessModel::new(graph, 0.1, 1.0).expect("model");
    let ensemble = contact_ensemble(&model, 0.01, 8, 20.0, 0x4138_0004).expect("ensemble");
    eprintln!(
        "X4 absorbing control: surviving {}/{}",
        ensemble.survived, ensemble.replicas
    );
    assert_eq!(ensemble.survived, 0);
    assert_eq!(ensemble.survival_fraction, 0.0);
    assert_eq!(ensemble.mean_density, 0.0);
}
