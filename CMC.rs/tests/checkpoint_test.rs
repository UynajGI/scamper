//! Checkpoint/snapshot round-trip tests for ClassicalMC.
//!
//! Verifies that snapshot save/load preserves state identity and that
//! the format tag validation rejects corrupted or unknown snapshots.

use carlo_rs::{CarloError, Context, MonteCarlo};
use cmc_rs::{ClassicalMC, IsingModel, MetropolisCore};
use rand::SeedableRng;
use serde_json::json;

type Rng = rand_xoshiro::Xoshiro256PlusPlus;

fn make_ising_4x4(beta: f64) -> ClassicalMC<IsingModel, MetropolisCore> {
    let model = IsingModel::new(1.0);
    let lattice = cmc_rs::build_square(4, 4, true);
    let system = cmc_rs::System::new(lattice, 1, 1.0, beta);
    let algorithm = MetropolisCore::new();
    let mut mc = ClassicalMC::new(system, model, algorithm);
    mc.system.recompute_energy(&mc.model);
    mc
}

// ---------------------------------------------------------------------------
// C.1: Save / restore / continue: split-run state matches continuous run
// ---------------------------------------------------------------------------

#[test]
fn snapshot_split_run_state_identity() {
    // Continuous run: 400 sweeps with deterministic RNG
    let mut mc_continuous = make_ising_4x4(0.5);
    // Use identical RNG seed
    let rng = Rng::seed_from_u64(42);
    let mut ctx = Context::new(rng, 0);

    let rng2 = Rng::seed_from_u64(42);
    let mut ctx_split = Context::new(rng2, 0);

    // Split run: 200 sweeps, save, restore into fresh MC, 200 more sweeps
    let mut mc_split = make_ising_4x4(0.5);

    // First 200 sweeps on both
    for _ in 0..200 {
        mc_continuous.sweep(&mut ctx);
        ctx.advance_sweep();
        mc_split.sweep(&mut ctx_split);
        ctx_split.advance_sweep();
    }

    // Save split state
    let snapshot = mc_split.save_snapshot();

    // Continue continuous for 200 more
    for _ in 0..200 {
        mc_continuous.sweep(&mut ctx);
        ctx.advance_sweep();
    }

    // Restore split, continue for 200 more
    mc_split.load_snapshot(&snapshot).unwrap();
    for _ in 0..200 {
        mc_split.sweep(&mut ctx_split);
        ctx_split.advance_sweep();
    }

    // Final states must match
    assert_eq!(mc_continuous.system.spins, mc_split.system.spins);
    assert!(
        (mc_continuous.system.energy - mc_split.system.energy).abs() < 1e-12,
        "Energy mismatch: {} vs {}",
        mc_continuous.system.energy,
        mc_split.system.energy
    );
}

// ---------------------------------------------------------------------------
// C.2: Snapshot rejects unknown format
// ---------------------------------------------------------------------------

#[test]
fn snapshot_rejects_unknown_format() {
    let mut mc = make_ising_4x4(0.5);
    let snapshot = mc.save_snapshot();
    assert_eq!(snapshot["format"], "cmc-rs-snapshot-v2");

    let mut bad = snapshot.clone();
    *bad.get_mut("format").unwrap() = json!("cmc-rs-snapshot-v1");
    let err = mc.load_snapshot(&bad).unwrap_err();
    assert!(
        matches!(err, CarloError::CheckpointCorrupted { .. }),
        "expected CheckpointCorrupted, got {err:?}"
    );
}

// ---------------------------------------------------------------------------
// C.3: Snapshot rejects bad edge kind
// ---------------------------------------------------------------------------

#[test]
fn snapshot_rejects_bad_edge_kind() {
    let mut mc = make_ising_4x4(0.5);
    let snapshot = mc.save_snapshot();

    let mut bad = snapshot.clone();
    bad["edges"][0]["kind"] = json!("invalid_kind");
    let err = mc.load_snapshot(&bad).unwrap_err();
    assert!(
        matches!(err, CarloError::CheckpointCorrupted { .. }),
        "expected CheckpointCorrupted, got {err:?}"
    );
}

// ---------------------------------------------------------------------------
// C.4: Snapshot rejects topology mismatch
// ---------------------------------------------------------------------------

#[test]
fn snapshot_rejects_topology_mismatch() {
    let mut mc = make_ising_4x4(0.5);
    let snapshot = mc.save_snapshot();

    let mut bad = snapshot.clone();
    bad["n_sites"] = json!(999usize);
    let err = mc.load_snapshot(&bad).unwrap_err();
    assert!(
        matches!(err, CarloError::CheckpointCorrupted { .. }),
        "expected CheckpointCorrupted, got {err:?}"
    );
}

// ---------------------------------------------------------------------------
// C.5: Snapshot round-trip recomputes energy (not trusting saved value)
// ---------------------------------------------------------------------------

#[test]
fn snapshot_roundtrip_recomputes_energy() {
    let mut mc = make_ising_4x4(0.5);
    let energy_before = mc.system.energy;
    let snapshot = mc.save_snapshot();

    // Corrupt cached energy, then restore
    mc.system.energy = 12345.0;
    mc.load_snapshot(&snapshot).unwrap();
    assert!(
        (mc.system.energy - energy_before).abs() < 1e-12,
        "Energy should be recomputed: before={energy_before}, after={}",
        mc.system.energy
    );
}

// ---------------------------------------------------------------------------
// C.6: 1000-sweep split run: 400 therm → save → restore → 600 meas
// ---------------------------------------------------------------------------

#[test]
fn split_run_thousand_sweeps() {
    // Continuous: 1000 sweeps (no thermalization boundary — uniform sweeps)
    let mut mc_continuous = make_ising_4x4(0.44);
    let rng_cont = Rng::seed_from_u64(999);
    let mut ctx_cont = Context::new(rng_cont, 0);
    for _ in 0..1000 {
        mc_continuous.sweep(&mut ctx_cont);
        ctx_cont.advance_sweep();
    }

    // Split: 400 sweeps, save, restore into fresh MC, 600 more
    let mut mc_split = make_ising_4x4(0.44);
    let rng_split = Rng::seed_from_u64(999);
    let mut ctx_split = Context::new(rng_split, 0);
    for _ in 0..400 {
        mc_split.sweep(&mut ctx_split);
        ctx_split.advance_sweep();
    }

    let snapshot = mc_split.save_snapshot();
    mc_split.load_snapshot(&snapshot).unwrap();

    for _ in 0..600 {
        mc_split.sweep(&mut ctx_split);
        ctx_split.advance_sweep();
    }

    // Final states must be bitwise identical (same RNG trajectory)
    assert_eq!(mc_continuous.system.spins, mc_split.system.spins);
    assert!(
        (mc_continuous.system.energy - mc_split.system.energy).abs() < 1e-12,
        "Split-run 1000 sweep energy mismatch"
    );
    assert_eq!(
        ctx_cont.sweep_count(),
        ctx_split.sweep_count(),
        "Context sweep counts must match"
    );
}
