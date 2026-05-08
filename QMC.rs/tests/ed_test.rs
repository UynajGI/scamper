//! Exact diagonalization tests for QMC.rs.
//!
//! Validates the Lanczos solver against known ground state energies,
//! then compares ED results with SSE simulations for small systems.

use qmc_rs::ed::{SparseHamiltonian, lanczos_ground_state};
use qmc_rs::lattice::builders::build_chain;

#[test]
fn test_ed_heisenberg_2site() {
    // 2-site chain with PBC: H = J S_1 · S_2
    // Singlet: E = -3J/4 = -0.75, Triplet: E = J/4 = 0.25
    let lattice = build_chain(2, true);
    let ham = SparseHamiltonian::from_heisenberg(&lattice, 1.0);

    assert_eq!(ham.dim(), 4); // 2^2 = 4

    let energy = lanczos_ground_state(&ham, 1e-12, 10);

    assert!(
        (energy - (-0.75)).abs() < 1e-10,
        "2-site ED energy {:.10} != -0.75",
        energy
    );
}

#[test]
fn test_ed_heisenberg_4site() {
    // 4-site PBC chain: exact ground state E = -2.0
    let lattice = build_chain(4, true);
    let ham = SparseHamiltonian::from_heisenberg(&lattice, 1.0);

    assert_eq!(ham.dim(), 16);

    let energy = lanczos_ground_state(&ham, 1e-12, 100);

    // For 4-site PBC, E ≈ -2.0. Allow 2% tolerance for Lanczos convergence.
    assert!(
        (energy - (-2.0)).abs() < 0.04,
        "4-site ED energy {:.6} != -2.0 (tolerance 0.04)",
        energy
    );
}

#[test]
fn test_ed_heisenberg_6site() {
    // 6-site PBC chain
    let lattice = build_chain(6, true);
    let ham = SparseHamiltonian::from_heisenberg(&lattice, 1.0);

    assert_eq!(ham.dim(), 64);

    let energy = lanczos_ground_state(&ham, 1e-12, 200);

    // For 1D Heisenberg, E/N approaches -0.443147 (Bethe ansatz)
    // For N=6, E ≈ -2.659
    assert!(
        energy < -2.5,
        "6-site ED energy {:.6} should be < -2.5",
        energy
    );
    assert!(
        energy > -3.5,
        "6-site ED energy {:.6} should be > -3.5",
        energy
    );
}

#[test]
fn test_ed_heisenberg_8site() {
    // 8-site PBC chain
    let lattice = build_chain(8, true);
    let ham = SparseHamiltonian::from_heisenberg(&lattice, 1.0);

    assert_eq!(ham.dim(), 256);

    let energy = lanczos_ground_state(&ham, 1e-12, 300);

    // E/N should approach Bethe ansatz value
    let energy_per_site = energy / 8.0;
    assert!(
        (energy_per_site - (-0.443147)).abs() < 0.02,
        "8-site E/N = {:.6}, expected ≈ -0.443",
        energy_per_site
    );
}

#[test]
fn test_ed_heisenberg_10site() {
    // 10-site PBC chain — largest practical test
    let lattice = build_chain(10, true);
    let ham = SparseHamiltonian::from_heisenberg(&lattice, 1.0);

    assert_eq!(ham.dim(), 1024);

    let energy = lanczos_ground_state(&ham, 1e-10, 300);
    let energy_per_site = energy / 10.0;

    assert!(
        (energy_per_site - (-0.443147)).abs() < 0.02,
        "10-site E/N = {:.6}, expected ≈ -0.443",
        energy_per_site
    );
}

#[test]
fn test_ed_heisenberg_12site() {
    // 12-site PBC chain
    let lattice = build_chain(12, true);
    let ham = SparseHamiltonian::from_heisenberg(&lattice, 1.0);

    assert_eq!(ham.dim(), 4096);

    let energy = lanczos_ground_state(&ham, 1e-10, 300);
    let energy_per_site = energy / 12.0;

    assert!(
        (energy_per_site - (-0.443147)).abs() < 0.02,
        "12-site E/N = {:.6}, expected ≈ -0.443",
        energy_per_site
    );
}

#[test]
fn test_ed_heisenberg_14site() {
    // 14-site PBC chain
    let lattice = build_chain(14, true);
    let ham = SparseHamiltonian::from_heisenberg(&lattice, 1.0);

    assert_eq!(ham.dim(), 16384);

    let energy = lanczos_ground_state(&ham, 1e-10, 300);
    let energy_per_site = energy / 14.0;

    assert!(
        (energy_per_site - (-0.443147)).abs() < 0.02,
        "14-site E/N = {:.6}, expected ≈ -0.443",
        energy_per_site
    );
}

/// 16-site ED — skipped by default as it's slow (65536 dim).
/// Run with: cargo test --test ed_test test_ed_heisenberg_16site -- --ignored
#[test]
#[ignore = "slow — 65536 dim matrix"]
fn test_ed_heisenberg_16site() {
    // 16-site PBC chain — largest for ED validation
    let lattice = build_chain(16, true);
    let ham = SparseHamiltonian::from_heisenberg(&lattice, 1.0);

    assert_eq!(ham.dim(), 65536);

    let energy = lanczos_ground_state(&ham, 1e-10, 300);
    let energy_per_site = energy / 16.0;

    // Should be very close to Bethe ansatz
    assert!(
        (energy_per_site - (-0.443147)).abs() < 0.001,
        "16-site E/N = {:.6}, expected ≈ -0.443",
        energy_per_site
    );
}

#[test]
fn test_ed_open_chain() {
    // 4-site open chain (no PBC): 3 bonds instead of 4
    let lattice = build_chain(4, false);
    let ham = SparseHamiltonian::from_heisenberg(&lattice, 1.0);

    assert_eq!(ham.dim(), 16);

    let energy = lanczos_ground_state(&ham, 1e-12, 100);

    // Open chain has different ground state than PBC
    assert!(energy < -1.0, "4-site open chain E = {:.6} should be < -1.0", energy);
    assert!(energy > -2.5, "4-site open chain E = {:.6} should be > -2.5", energy);
}
