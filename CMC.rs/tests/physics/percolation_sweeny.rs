//! Q1/Q2 validation: Sweeny exact-FK marginal gates, q = 1 reduction,
//! and Chayes-Machta agreement with Sweeny at integer q.
//!
//! The exact random-cluster measure on a small graph is enumerated in full
//! (`2^E` configurations, own union-find for `k(A)`) and compared with the
//! Markov-chain marginals of `|A|` and `k(A)`. The `q = 1` chain must reduce
//! bitwise to the F2 `BondBernoulli` draw on one all-closed sweep.

use cmc_rs::{
    build_square, BondBernoulli, ChayesMachtaMC, EdgeActivity, GraphView, ObservablePlan,
    Probability, StaticConfiguration, SweenyAcceptance, SweenyMC, UndirectedGraphView,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn sweep_once(mc: &mut SweenyMC, rng: &mut Xoshiro256PlusPlus) {
    for edge in 0..mc.lattice().edge_count() {
        mc.propose(edge, rng).expect("edge index in range");
    }
}

/// Component count of one bond mask (own union-find).
fn component_count(lattice: &cmc_rs::CsrLattice, mask: usize) -> usize {
    let vertex_count = lattice.vertex_count();
    let mut parent: Vec<usize> = (0..vertex_count).collect();
    fn find(parent: &mut [usize], mut vertex: usize) -> usize {
        while parent[vertex] != vertex {
            parent[vertex] = parent[parent[vertex]];
            vertex = parent[vertex];
        }
        vertex
    }
    let mut components = vertex_count;
    for edge in 0..lattice.edge_count() {
        if (mask >> edge) & 1 == 1 {
            let bond = &lattice.edges[edge];
            let (left, right) = (
                find(&mut parent, bond.source),
                find(&mut parent, bond.target),
            );
            if left != right {
                parent[left] = right;
                components -= 1;
            }
        }
    }
    components
}

/// Exact unnormalized FK weights and the moments of `|A|` and `k(A)`.
fn exact_moments(lattice: &cmc_rs::CsrLattice, q: f64, p: f64) -> (f64, f64, f64, f64) {
    let edge_count = lattice.edge_count();
    let odds = p / (1.0 - p);
    let mut normalization = 0.0;
    let mut mean_bonds = 0.0;
    let mut mean_square_bonds = 0.0;
    let mut mean_components = 0.0;
    let mut mean_square_components = 0.0;
    for mask in 0_usize..(1_usize << edge_count) {
        let components = component_count(lattice, mask);
        let bonds = mask.count_ones() as f64;
        let weight = odds.powf(bonds) * q.powi(components as i32);
        normalization += weight;
        mean_bonds += weight * bonds;
        mean_square_bonds += weight * bonds * bonds;
        mean_components += weight * components as f64;
        mean_square_components += weight * (components * components) as f64;
    }
    let mean_bonds = mean_bonds / normalization;
    let variance_bonds = mean_square_bonds / normalization - mean_bonds * mean_bonds;
    let mean_components = mean_components / normalization;
    let variance_components =
        mean_square_components / normalization - mean_components * mean_components;
    (
        mean_bonds,
        variance_bonds,
        mean_components,
        variance_components,
    )
}

#[test]
fn sweeny_marginals_match_exact_fk_enumeration() {
    let lattice = build_square(3, 3, false);
    for &(q, p) in &[(1.0_f64, 0.4_f64), (2.0, 0.4), (3.0, 0.55)] {
        let mut chain = SweenyMC::new(
            lattice.clone(),
            q,
            p,
            SweenyAcceptance::HeatBath,
            ObservablePlan::all(vec![]),
        )
        .expect("chain");
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x51EE_0001 + (q * 10.0) as u64);
        for _ in 0..100 {
            sweep_once(&mut chain, &mut rng);
        }
        let samples = 4000;
        let mut bond_sum = 0.0;
        let mut bond_square_sum = 0.0;
        let mut component_sum = 0.0;
        for _ in 0..samples {
            sweep_once(&mut chain, &mut rng);
            let bonds = chain.open_edge_count() as f64;
            bond_sum += bonds;
            bond_square_sum += bonds * bonds;
            component_sum += chain.component_count() as f64;
        }
        let n = samples as f64;
        let empirical_bonds = bond_sum / n;
        let empirical_variance = bond_square_sum / n - empirical_bonds * empirical_bonds;
        let empirical_components = component_sum / n;
        let (exact_bonds, exact_variance, exact_components, exact_component_variance) =
            exact_moments(&lattice, q, p);
        let bond_sem = (exact_variance.max(0.0) / n).sqrt();
        let bond_z = (empirical_bonds - exact_bonds) / bond_sem;
        let component_z = (empirical_components - exact_components)
            / (exact_component_variance.max(0.0) / n).sqrt();
        eprintln!(
            "Q1 exact q={q} p={p}: bonds {empirical_bonds:.4} vs {exact_bonds:.4} (z={bond_z:.2}); \
             components {empirical_components:.4} vs {exact_components:.4} (z={component_z:.2}); \
             empirical var {empirical_variance:.4}"
        );
        assert!(
            bond_z.abs() < 5.0,
            "q={q} p={p}: bond marginal z={bond_z:.2}"
        );
    }
}

#[test]
fn q_one_heat_bath_reduces_bitwise_to_bond_bernoulli() {
    let lattice = build_square(4, 4, false);
    let p = 0.37_f64;
    let seed = 0x51EE_1001_u64;
    let mut chain = SweenyMC::new(
        lattice.clone(),
        1.0, // q = 1
        p,
        SweenyAcceptance::HeatBath,
        ObservablePlan::all(vec![]),
    )
    .expect("chain");
    sweep_once(&mut chain, &mut Xoshiro256PlusPlus::seed_from_u64(seed));

    let mut configuration = StaticConfiguration::new(lattice.vertex_count(), lattice.edge_count());
    BondBernoulli::new(Probability::new(p).expect("probability").into())
        .sample(
            &lattice,
            &mut configuration,
            &mut Xoshiro256PlusPlus::seed_from_u64(seed),
        )
        .expect("Bernoulli draw");
    for edge in 0..lattice.edge_count() {
        assert_eq!(
            chain.mask()[edge],
            configuration.edge_active(lattice.edge_id(edge).expect("edge id")),
            "edge {edge}"
        );
    }
    eprintln!("Q1 q=1: one all-closed heat-bath sweep equals BondBernoulli bitwise");
}

#[test]
fn chayes_machta_agrees_with_sweeny_at_q_two() {
    let lattice = build_square(3, 3, true);
    let p = 0.45_f64;
    let mut sweeny = SweenyMC::new(
        lattice.clone(),
        2.0,
        p,
        SweenyAcceptance::HeatBath,
        ObservablePlan::all(vec![]),
    )
    .expect("Sweeny");
    let mut chayes = ChayesMachtaMC::new(lattice.clone(), 2, p, ObservablePlan::all(vec![]))
        .expect("Chayes-Machta");
    let mut sweeny_rng = Xoshiro256PlusPlus::seed_from_u64(0x51EE_2001);
    let mut chayes_rng = Xoshiro256PlusPlus::seed_from_u64(0x51EE_2002);
    for _ in 0..200 {
        sweep_once(&mut sweeny, &mut sweeny_rng);
        chayes.step(&mut chayes_rng);
    }
    let samples = 4000;
    let mut sweeny_sum = 0.0;
    let mut sweeny_square = 0.0;
    let mut chayes_sum = 0.0;
    let mut chayes_square = 0.0;
    for _ in 0..samples {
        sweep_once(&mut sweeny, &mut sweeny_rng);
        chayes.step(&mut chayes_rng);
        let sweeny_bonds = sweeny.open_edge_count() as f64;
        let chayes_bonds = chayes.open_edge_count() as f64;
        sweeny_sum += sweeny_bonds;
        sweeny_square += sweeny_bonds * sweeny_bonds;
        chayes_sum += chayes_bonds;
        chayes_square += chayes_bonds * chayes_bonds;
    }
    let n = samples as f64;
    let sweeny_mean = sweeny_sum / n;
    let chayes_mean = chayes_sum / n;
    let sweeny_variance = (sweeny_square / n - sweeny_mean * sweeny_mean).max(0.0);
    let chayes_variance = (chayes_square / n - chayes_mean * chayes_mean).max(0.0);
    let difference = (sweeny_mean - chayes_mean).abs();
    let noise = 4.0 * (sweeny_variance / n + chayes_variance / n).sqrt();
    eprintln!(
        "Q2 vs Q1 q=2: Sweeny |A|={sweeny_mean:.4} Chayes-Machta |A|={chayes_mean:.4} \
         difference {difference:.4} noise {noise:.4}"
    );
    assert!(
        difference <= noise,
        "q=2 bond density differs beyond 4 sigma: {difference:.4} > {noise:.4}"
    );
}

#[test]
fn heat_bath_rate_matches_exact_conditional_on_every_state() {
    // Detailed-balance audit at the rate level: for every bond mask A and
    // every edge e, the heat-bath insertion probability must equal the exact
    // FK two-state conditional P(e=1 | rest) = q^{k1} v / (q^{k0} + q^{k1} v).
    let lattice = build_square(3, 3, false);
    let edge_count = lattice.edge_count();
    let p = 0.4_f64;
    let odds = p / (1.0 - p);
    for &q in &[1.0_f64, 2.0, 3.0] {
        let rates = cmc_rs::SweenyRates::new(q, p).expect("rates");
        for mask in 0_usize..(1_usize << edge_count) {
            for edge in 0..edge_count {
                let without = mask & !(1_usize << edge);
                let with = mask | (1_usize << edge);
                let delta_n = component_count(&lattice, with) as i32
                    - component_count(&lattice, without) as i32;
                let ratio = q.powi(delta_n) * odds;
                let exact = ratio / (1.0 + ratio);
                let rate = rates.heat_bath(delta_n, 1);
                assert!(
                    (rate - exact).abs() < 1e-12,
                    "q={q} mask={mask:b} edge={edge}: rate {rate} vs exact {exact}"
                );
                let metropolis = rates.metropolis(delta_n, 1);
                assert!((metropolis - ratio.min(1.0)).abs() < 1e-12);
            }
        }
        eprintln!(
            "Q1 detailed-balance rates: q={q} all masks x edges match the exact FK conditional"
        );
    }
}
