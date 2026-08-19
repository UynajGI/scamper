//! Multi-component worm validation: exact enumeration, partition identity,
//! checkpoint round-trip, and cross-solver agreement.
//!
//! Closes the "multi-component worm" residue by validating the real
//! implementation (not a rejection): the Ising high-temperature graph ensemble
//! factorizes over connected components, and [`IsingGraphWormMC::from_lattice`]
//! samples it with one independent two-defect worm per component on
//! domain-separated derived streams.
//!
//! - **Exact enumeration** (full 2^N spin enumeration, independent reference):
//!   total ⟨E⟩, per-component ⟨E⟩, and same-component two-point correlations
//!   ⟨s_i s_j⟩ (worm endpoint-pair ratio estimator) on a 4-ring + 3-chain +
//!   isolated-site lattice, plus the spin order parameter ⟨m²⟩ reconstructed
//!   from the worm pair correlations (cross-component pairs factorize to zero
//!   exactly in the product ensemble).
//! - **Partition identity** (machine precision): Z_spin = 2^N Π_e cosh(βJ_e) ·
//!   Π_c Z_graph,c against the per-component exact graph enumerations.
//! - **Checkpoints**: v2 multi-component snapshots round-trip bit-exact
//!   trajectories; v1 snapshots cannot silently restore a multi-component
//!   ensemble; corrupted component partitions are rejected.
//! - **Cross-solver**: worm vs spin Metropolis on two disjoint 4×4 squares.

use super::common::{assert_close, enumerate_ising, zscore_seed_count};
use carlo_rs::MonteCarlo;
use cmc_rs::{
    Algorithm, Bond, BondType, CsrLattice, Initializable, IsingGraphWormEnsemble, IsingGraphWormMC,
    IsingModel, MetropolisCore, SimulationPhase, System, WormConfig,
};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

type Rng = Xoshiro256PlusPlus;

/// 4-ring (sites 0-3) + open 3-chain (sites 4-6) + isolated site 7.
fn three_component_lattice() -> CsrLattice {
    CsrLattice::from_edges(
        8,
        vec![
            Bond::new(0, 1, BondType::Generic, 1.0),
            Bond::new(1, 2, BondType::Generic, 1.0),
            Bond::new(2, 3, BondType::Generic, 1.0),
            Bond::new(3, 0, BondType::Generic, 1.0),
            Bond::new(4, 5, BondType::Generic, 1.0),
            Bond::new(5, 6, BondType::Generic, 1.0),
        ],
    )
}

fn tracking_config(local_updates_per_sweep: usize) -> WormConfig {
    // Low worm fugacity keeps every component mostly in the physical sector,
    // raising the all-physical sampling fraction (the physical-sector
    // conditional ensemble is η-independent; η only splits time).
    WormConfig {
        local_updates_per_sweep,
        close_probability: 0.25,
        log_worm_fugacity: (1.0f64 / 64.0).ln(),
        track_endpoint_pairs: true,
        cache_audit_interval: 101,
    }
}

fn component_energy(mc: &IsingGraphWormMC, component: usize) -> Option<f64> {
    let component = mc.ensemble().components().get(component)?;
    if !component.kernel().state().is_physical() {
        return None;
    }
    Some(
        component
            .kernel()
            .model()
            .energy_estimator(component.kernel().state().configuration()),
    )
}

fn spin_energy(spins: &[f64], lattice: &CsrLattice, coupling: f64, edges: &[usize]) -> f64 {
    edges
        .iter()
        .map(|&edge_id| {
            let edge = lattice.edges[edge_id];
            -coupling * edge.weight * spins[edge.source] * spins[edge.target]
        })
        .sum()
}

/// Exact (⟨E_total⟩, ⟨E_ring⟩, ⟨E_chain⟩, ⟨s0 s2⟩, ⟨s4 s6⟩, ⟨m²⟩) by 2^8
/// spin enumeration of the three-component lattice.
fn exact_three_component_moments(
    lattice: &CsrLattice,
    coupling: f64,
    beta: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let mut z = 0.0;
    let mut totals = [0.0f64; 6];
    for spins in enumerate_ising(lattice.n_sites) {
        let energy = spin_energy(&spins, lattice, coupling, &[0, 1, 2, 3, 4, 5]);
        let weight = (-beta * energy).exp();
        let magnetization = spins.iter().sum::<f64>() / lattice.n_sites as f64;
        z += weight;
        totals[0] += weight * energy;
        totals[1] += weight * spin_energy(&spins, lattice, coupling, &[0, 1, 2, 3]);
        totals[2] += weight * spin_energy(&spins, lattice, coupling, &[4, 5]);
        totals[3] += weight * spins[0] * spins[2];
        totals[4] += weight * spins[4] * spins[6];
        totals[5] += weight * magnetization * magnetization;
    }
    (
        totals[0] / z,
        totals[1] / z,
        totals[2] / z,
        totals[3] / z,
        totals[4] / z,
        totals[5] / z,
    )
}

/// ⟨m²⟩ reconstructed from worm pair correlations: intra-component pairs from
/// the endpoint-pair ratio estimator (diagonal exactly 1), cross-component
/// pairs exactly 0 (independent components, ⟨s_i⟩ = 0).
fn worm_magnetization_squared(mc: &IsingGraphWormMC) -> f64 {
    let n = mc.ensemble().lattice().n_sites;
    let mut sum = 0.0;
    for tail in 0..n {
        for head in 0..n {
            if tail == head {
                sum += 1.0;
            } else if let Some(correlation) =
                mc.endpoint_correlation(tail.min(head), tail.max(head))
            {
                sum += correlation;
            }
        }
    }
    sum / (n * n) as f64
}

#[test]
fn worm_multi_component_matches_exact_spin_enumeration() {
    let beta = 0.45;
    let coupling = 1.0;
    let lattice = three_component_lattice();
    let (e_total, e_ring, e_chain, c02, c46, m2) =
        exact_three_component_moments(&lattice, coupling, beta);

    let mut mc =
        IsingGraphWormMC::from_lattice(lattice, beta, coupling, tracking_config(8)).unwrap();
    assert_eq!(mc.n_components(), 3);

    // Drive the real adapter path (sweep + measure) through a Context; a
    // single long measurement loop feeds both the adapter observables and the
    // component-local accumulation below.
    let mut context = carlo_rs::Context::new(Rng::seed_from_u64(0x3C0), 0);
    for _ in 0..50_000 {
        mc.sweep(&mut context);
    }

    // Per-component accumulation on the same chain. Component-local
    // conditioning is valid: the joint extended-sector weight factorizes and
    // so does the "component physical" event.
    let mut ring_samples = 0u64;
    let mut ring_energy_sum = 0.0;
    let mut chain_samples = 0u64;
    let mut chain_energy_sum = 0.0;
    for _ in 0..30_000 {
        context.advance_sweep();
        mc.sweep(&mut context);
        mc.measure(&mut context);
        if let Some(energy) = component_energy(&mc, 0) {
            ring_samples += 1;
            ring_energy_sum += energy;
        }
        if let Some(energy) = component_energy(&mc, 1) {
            chain_samples += 1;
            chain_energy_sum += energy;
        }
    }
    mc.ensemble().validate().unwrap();
    assert!(ring_samples > 10_000 && chain_samples > 10_000);

    let measured_ring = ring_energy_sum / ring_samples as f64;
    let measured_chain = chain_energy_sum / chain_samples as f64;
    eprintln!(
        "[worm-multi] ring ⟨E⟩ {measured_ring:.4} vs exact {e_ring:.4} | chain ⟨E⟩ \
         {measured_chain:.4} vs exact {e_chain:.4} | C(0,2) {:.4} vs {c02:.4} | C(4,6) {:.4} vs \
         {c46:.4} | worm m² {:.4} vs exact {m2:.4}",
        mc.endpoint_correlation(0, 2).unwrap(),
        mc.endpoint_correlation(4, 6).unwrap(),
        worm_magnetization_squared(&mc),
    );
    assert!(
        (measured_ring - e_ring).abs() < 0.08,
        "ring component ⟨E⟩ = {measured_ring:.4}, exact {e_ring:.4}"
    );
    assert!(
        (measured_chain - e_chain).abs() < 0.08,
        "chain component ⟨E⟩ = {measured_chain:.4}, exact {e_chain:.4}"
    );

    // Correlations: worm endpoint-pair ratios vs exact spin correlations.
    let measured_c02 = mc.endpoint_correlation(0, 2).expect("ring pair estimator");
    let measured_c46 = mc.endpoint_correlation(4, 6).expect("chain pair estimator");
    assert!(
        (measured_c02 - c02).abs() < 0.06,
        "C(0,2) = {measured_c02:.4}, exact {c02:.4}"
    );
    assert!(
        (measured_c46 - c46).abs() < 0.06,
        "C(4,6) = {measured_c46:.4}, exact {c46:.4}"
    );
    // Cross-component pairs have no worm estimator (they factorize to zero).
    assert!(mc.endpoint_correlation(0, 4).is_none());

    // Order parameter reconstructed from the worm pair correlations.
    let worm_m2 = worm_magnetization_squared(&mc);
    assert!(
        (worm_m2 - m2).abs() < 0.06,
        "worm-reconstructed ⟨m²⟩ = {worm_m2:.4}, exact {m2:.4}"
    );

    // Total energy through the real adapter measure() path (all-physical
    // conditioning preserves the product ensemble).
    let estimates = context.finalize_measurements();
    let energy = estimates
        .get("Energy")
        .expect("Energy observable from the all-physical sector");
    eprintln!(
        "[worm-multi] adapter ⟨E⟩ = {:.4} ± {:.4} vs exact {e_total:.4}",
        energy.mean, energy.stderr
    );
    assert!(
        (energy.mean - e_total).abs() < 0.08,
        "adapter ⟨E⟩ = {:.4}, exact {e_total:.4}",
        energy.mean
    );
    assert!(estimates.contains_key("WormSector"));
}

#[test]
fn worm_multi_component_partition_identity_factorizes() {
    // Z_spin = 2^N · Π_e cosh(βJ_e) · Π_c Z_graph,c, machine precision.
    let beta = 0.47;
    let coupling = 1.2;
    let lattice = three_component_lattice();
    let ensemble =
        IsingGraphWormEnsemble::new(lattice.clone(), beta, coupling, tracking_config(4)).unwrap();

    let log_graph_partition: f64 = ensemble
        .components()
        .iter()
        .map(|component| {
            cmc_rs::enumerate_ising_graph_expansion(component.model())
                .unwrap()
                .log_reduced_partition
        })
        .sum();
    let spin_partition: f64 = enumerate_ising(lattice.n_sites)
        .iter()
        .map(|spins| (-beta * spin_energy(spins, &lattice, coupling, &[0, 1, 2, 3, 4, 5])).exp())
        .sum();
    let prefactor = 2.0f64.powi(lattice.n_sites as i32)
        * lattice
            .edges
            .iter()
            .map(|edge| (beta * coupling * edge.weight).cosh())
            .product::<f64>();
    assert_close(spin_partition, prefactor * log_graph_partition.exp(), 1e-10);
}

#[test]
fn worm_multi_component_snapshot_round_trips_exact_trajectory() {
    let lattice = three_component_lattice();
    let config = tracking_config(6);
    let mut original =
        IsingGraphWormMC::from_lattice(lattice.clone(), 0.4, 1.0, config.clone()).unwrap();
    let mut rng = Rng::seed_from_u64(0x0DD);
    for _ in 0..300 {
        original.ensemble_mut().sweep(&mut rng).unwrap();
    }
    let snapshot = original.save_snapshot();
    assert_eq!(
        snapshot["format"].as_str(),
        Some("cmc-rs-ising-worm-v2"),
        "multi-component ensembles must use the v2 layout"
    );

    let mut restored = IsingGraphWormMC::from_lattice(lattice, 0.4, 1.0, config).unwrap();
    restored.load_snapshot(&snapshot).unwrap();
    assert_eq!(restored.save_snapshot(), snapshot);

    // Future trajectories stay identical (no hidden RNG state).
    let mut replay = rng.clone();
    for _ in 0..1_000 {
        original.ensemble_mut().sweep(&mut rng).unwrap();
        restored.ensemble_mut().sweep(&mut replay).unwrap();
    }
    assert_eq!(original.save_snapshot(), restored.save_snapshot());

    // A v1 (single-component) snapshot cannot restore a multi-component
    // ensemble — it must fail loudly, not silently misassign components.
    let mut v1 = snapshot.clone();
    v1["format"] = serde_json::json!("cmc-rs-ising-worm-v1");
    let mut fresh =
        IsingGraphWormMC::from_lattice(three_component_lattice(), 0.4, 1.0, tracking_config(6))
            .unwrap();
    assert!(fresh.load_snapshot(&v1).is_err());

    // Corrupting the component partition is rejected.
    let mut corrupted = snapshot;
    corrupted["components"][0]["sites"][0] = serde_json::json!(7);
    assert!(fresh.load_snapshot(&corrupted).is_err());
}

// ── Cross-solver: worm vs spin Metropolis on a disconnected geometry ────────

fn two_disjoint_squares() -> CsrLattice {
    let left = cmc_rs::build_square(4, 4, true);
    let right = cmc_rs::build_square(4, 4, true);
    let offset = left.n_sites;
    let mut edges = left.edges.clone();
    for edge in &right.edges {
        edges.push(Bond::new(
            edge.source + offset,
            edge.target + offset,
            BondType::Generic,
            edge.weight,
        ));
    }
    CsrLattice::from_edges(2 * offset, edges)
}

fn binned_stats(samples: &[f64], binsize: usize) -> (f64, f64) {
    let usable = samples.len() / binsize * binsize;
    assert!(usable > 0, "no complete bins");
    let bins: Vec<f64> = samples[..usable]
        .chunks(binsize)
        .map(|chunk| chunk.iter().sum::<f64>() / chunk.len() as f64)
        .collect();
    let n = bins.len() as f64;
    let mean = bins.iter().sum::<f64>() / n;
    let variance = bins.iter().map(|b| (b - mean).powi(2)).sum::<f64>() / (n - 1.0);
    (mean, (variance / n).sqrt())
}

const THERM: usize = 2_000;
const MEAS: usize = 20_000;
const BIN: usize = 200;

fn run_worm_energy(lattice: &CsrLattice, beta: f64, seed: u64) -> (f64, f64) {
    let config = WormConfig {
        local_updates_per_sweep: lattice.n_edges(),
        close_probability: 0.25,
        // Small fugacity keeps both components mostly physical, so the
        // all-physical conditioning yields enough total-energy samples.
        log_worm_fugacity: (1.0f64 / 512.0).ln(),
        track_endpoint_pairs: false,
        cache_audit_interval: 0,
    };
    let mut mc = IsingGraphWormMC::from_lattice(lattice.clone(), beta, 1.0, config).unwrap();
    let mut rng = Rng::seed_from_u64(seed);
    for _ in 0..THERM {
        mc.ensemble_mut().sweep(&mut rng).unwrap();
    }
    let mut samples = Vec::with_capacity(MEAS);
    for _ in 0..MEAS {
        mc.ensemble_mut().sweep(&mut rng).unwrap();
        if let Some((energy, _)) = mc.ensemble().total_energy_and_occupied_edges() {
            samples.push(energy);
        }
    }
    binned_stats(&samples, BIN)
}

fn run_metropolis_energy(lattice: &CsrLattice, beta: f64, seed: u64) -> (f64, f64) {
    let model = IsingModel::new(1.0);
    let mut system = System::new(lattice.clone(), 1, 0.0, beta);
    let mut rng = Rng::seed_from_u64(seed);
    for site in 0..system.n_sites() {
        system.spins[site] = model.random_spin(&mut rng)[0];
    }
    system.recompute_energy(&model);
    let mut kernel = MetropolisCore::new();
    for _ in 0..THERM {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
    }
    let mut samples = Vec::with_capacity(MEAS);
    for _ in 0..MEAS {
        kernel.sweep_with_phase(&mut system, &model, &mut rng, SimulationPhase::Measurement);
        samples.push(system.energy);
    }
    binned_stats(&samples, BIN)
}

#[test]
fn worm_multi_component_cross_solver_matches_spin_metropolis() {
    let beta = 0.44;
    let lattice = two_disjoint_squares();
    let n_seeds = zscore_seed_count(8);
    let worm: Vec<(f64, f64)> = (0..n_seeds as u64)
        .map(|seed| run_worm_energy(&lattice, beta, 0x0E61 + seed))
        .collect();
    let metro: Vec<(f64, f64)> = (0..n_seeds as u64)
        .map(|seed| run_metropolis_energy(&lattice, beta, 0x0E62 + seed))
        .collect();
    let pool = |results: &[(f64, f64)]| {
        let n = results.len() as f64;
        let mean = results.iter().map(|(m, _)| m).sum::<f64>() / n;
        let sem = results.iter().map(|(_, s)| s * s).sum::<f64>().sqrt() / n;
        (mean, sem)
    };
    let (worm_mean, worm_sem) = pool(&worm);
    let (metro_mean, metro_sem) = pool(&metro);
    let z = (worm_mean - metro_mean) / (worm_sem * worm_sem + metro_sem * metro_sem).sqrt();
    eprintln!(
        "[worm-multi-cross] two 4x4 squares beta={beta}: worm ⟨E⟩ = {worm_mean:.4} ± \
         {worm_sem:.4}, metropolis ⟨E⟩ = {metro_mean:.4} ± {metro_sem:.4}, z = {z:+.2}"
    );
    assert!(
        z.abs() < 4.0,
        "worm vs Metropolis ⟨E⟩ disagree on the multi-component geometry: z = {z:.2} \
         ({worm_mean:.4} vs {metro_mean:.4})"
    );
}
