//! SSE algorithm unit tests.

use qmc_rs::hilbert::SpinHalfHS;
use qmc_rs::lattice::builders::{build_chain, build_square};
use qmc_rs::lattice::BondType;
use qmc_rs::sse::{SSEEngine, Vertex, VertexData, VertexList};
use qmc_rs::hilbert::OpType;
use std::collections::HashMap;
use rand_xoshiro::Xoshiro256PlusPlus;
use rand_xoshiro::rand_core::SeedableRng;

// Thresholds for equilibrium operator density in balance test
const MIN_DENSITY: f64 = 0.01;
const MAX_DENSITY: f64 = 0.60;
const MIN_FLUCTUATION: f64 = 0.005;

#[test]
fn test_bond_list_chain() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 1.0);
    let engine = SSEEngine::new(lattice, hs, 100, weights, 1.0, 0.0, 1.0);

    // Chain with 4 sites and PBC has 4 bonds: (0-1, 1-2, 2-3, 3-0)
    assert_eq!(engine.bond_list.len(), 4);

    // Verify bond content: each bond should have correct site pairs
    for (si, sj, bt) in &engine.bond_list {
        assert_eq!(*bt, BondType::ChainX);
        let (a, b) = if *si < *sj { (*si, *sj) } else { (*sj, *si) };
        assert!(matches!((a, b), (0, 1) | (1, 2) | (2, 3) | (0, 3)));
    }
}

#[test]
fn test_bond_list_chain_non_periodic() {
    let lattice = build_chain(4, false); // No PBC
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 1.0);
    let engine = SSEEngine::new(lattice, hs, 100, weights, 1.0, 0.0, 1.0);

    // Chain with 4 sites without PBC has 3 bonds: (0-1, 1-2, 2-3)
    assert_eq!(engine.bond_list.len(), 3);

    // Verify bond content
    let expected_bonds = vec![
        (0, 1, BondType::ChainX),
        (1, 2, BondType::ChainX),
        (2, 3, BondType::ChainX),
    ];
    assert_eq!(engine.bond_list, expected_bonds);
}

#[test]
fn test_bond_list_square() {
    let lattice = build_square(4, 4, true);
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::SquareX, 1.0);
    weights.insert(BondType::SquareY, 1.0);
    let engine = SSEEngine::new(lattice, hs, 100, weights, 1.0, 0.0, 1.0);

    // Square 4x4 with PBC: 2*N bonds for nearest-neighbor (horizontal + vertical)
    // N = 16, so 32 bonds
    assert_eq!(engine.bond_list.len(), 32);

    // Verify that we have the right distribution of bond types
    let mut count_x = 0;
    let mut count_y = 0;
    for (_, _, bond_type) in &engine.bond_list {
        match bond_type {
            BondType::SquareX => count_x += 1,
            BondType::SquareY => count_y += 1,
            _ => panic!("Unexpected bond type in square lattice"),
        }
    }
    assert_eq!(count_x, 16, "Expected 16 horizontal bonds");
    assert_eq!(count_y, 16, "Expected 16 vertical bonds");
}

#[test]
fn test_bond_list_square_non_periodic() {
    let lattice = build_square(4, 4, false); // No PBC
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::SquareX, 1.0);
    weights.insert(BondType::SquareY, 1.0);
    let engine = SSEEngine::new(lattice, hs, 100, weights, 1.0, 0.0, 1.0);

    // Square 4x4 without PBC:
    // Horizontal: (L-1)*L = 3*4 = 12 bonds
    // Vertical: L*(L-1) = 4*3 = 12 bonds
    // Total: 24 bonds
    assert_eq!(engine.bond_list.len(), 24);

    // Verify that we have the right distribution of bond types
    let mut count_x = 0;
    let mut count_y = 0;
    for (_, _, bond_type) in &engine.bond_list {
        match bond_type {
            BondType::SquareX => count_x += 1,
            BondType::SquareY => count_y += 1,
            _ => panic!("Unexpected bond type in square lattice"),
        }
    }
    assert_eq!(count_x, 12, "Expected 12 horizontal bonds");
    assert_eq!(count_y, 12, "Expected 12 vertical bonds");
}

#[test]
fn test_bond_type_validation() {
    let lattice = build_chain(4, true);

    // Test with missing bond type in weights (should panic in debug mode)
    let hs = SpinHalfHS;
    let incomplete_weights = HashMap::new(); // Empty - missing ChainX

    #[cfg(debug_assertions)]
    {
        // This should trigger the debug assertion
        let result = std::panic::catch_unwind(|| {
            SSEEngine::new(lattice.clone(), hs, 100, incomplete_weights.clone(), 1.0, 0.0, 1.0);
        });
        assert!(result.is_err(), "Expected panic for missing bond type");
    }

    // Test with correct weights (should succeed)
    let mut complete_weights = HashMap::new();
    complete_weights.insert(BondType::ChainX, 1.0);
    let engine = SSEEngine::new(lattice, hs, 100, complete_weights, 1.0, 0.0, 1.0);
    assert_eq!(engine.bond_list.len(), 4);
}

#[test]
fn test_beta_storage() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let beta = 2.5;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 1.0);
    let engine = SSEEngine::new(lattice, hs, 100, weights, beta, 0.0, 1.0);

    // Verify beta is stored correctly
    assert_eq!(engine.beta, beta);
}

#[test]
fn test_state_propagation() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 0.25);
    let mut engine = SSEEngine::new(lattice, hs, 100, weights, 1.0, 0.0, 1.0);

    // Insert an off-diagonal operator manually
    // Bond list order for 4-site PBC chain: [(0, 3), (0, 1), (1, 2), (2, 3)]
    // We use bond_idx 1 to connect sites 0 and 1
    engine.op_seq.vertices[0] = Vertex {
        bond_idx: 1,
        op: OpType::OffDiagonal,
        vertex_idx: 5,
    };
    engine.op_seq.n_operators = 1;

    // Initial spins: [0, 0, 0, 0] (all up)
    engine.spins = vec![0, 0, 0, 0];

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    engine.diagonal_update(&mut rng);

    // After propagation, spins at bond 1 sites should be flipped
    // Bond 1 connects sites 0 and 1
    // Original: [0, 0, ...], after flip: [1, 1, ...]
    assert_eq!(engine.spins[0], 1);
    assert_eq!(engine.spins[1], 1);
    // Other sites unchanged
    assert_eq!(engine.spins[2], 0);
    assert_eq!(engine.spins[3], 0);
}

#[test]
fn test_diagonal_insert_remove_balance() {
    // Test that insert and removal probabilities reach equilibrium
    let lattice = build_chain(6, true); // 6-site chain with PBC
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 1.0);

    // Higher beta for stronger signal
    let beta = 2.0;
    let max_length = 200;
    let mut engine = SSEEngine::new(lattice, hs, max_length, weights, beta, 0.0, 1.0);

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(12345);

    // Thermalization sweeps
    let therm_sweeps = 2000;
    for _ in 0..therm_sweeps {
        engine.diagonal_update(&mut rng);
    }

    // Measurement: track operator density over many sweeps
    let measure_sweeps = 5000;
    let mut total_density = 0.0;
    for _ in 0..measure_sweeps {
        engine.diagonal_update(&mut rng);
        let density = engine.op_seq.n_operators as f64 / engine.op_seq.max_length as f64;
        total_density += density;
    }
    let avg_density = total_density / measure_sweeps as f64;

    // At equilibrium, operator density should be non-trivial (between 1% and 60%)
    // The exact value depends on beta, weights, lattice size, and the diagonal
    // matrix elements which vary with spin configuration.
    // With beta=2.0, weight=1.0, 6-site chain (6 bonds), diagonal_element ~ 0.25:
    // P_insert ~ 2.0 * 0.25 * 6 / M ~ 3 / 200 ~ 1.5% per slot per sweep
    // Expected equilibrium density is modest but clearly non-zero.
    assert!(
        avg_density > MIN_DENSITY,
        "Operator density too low ({:.4}), insert not working",
        avg_density
    );
    assert!(
        avg_density < MAX_DENSITY,
        "Operator density too high ({:.4}), removal not working",
        avg_density
    );

    // Also verify density fluctuates (not stuck at fixed value)
    // This confirms both insert and removal are active
    let mut min_density = 1.0;
    let mut max_density = 0.0;
    for _ in 0..100 {
        for _ in 0..10 {
            engine.diagonal_update(&mut rng);
        }
        let density = engine.op_seq.n_operators as f64 / engine.op_seq.max_length as f64;
        if density < min_density {
            min_density = density;
        }
        if density > max_density {
            max_density = density;
        }
    }
    assert!(
        max_density - min_density > MIN_FLUCTUATION,
        "Density not fluctuating (range: {:.4}), both insert/remove should be active",
        max_density - min_density
    );
}

#[test]
fn test_energy_calculation() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 0.25);
    let mut engine = SSEEngine::new(lattice, hs, 100, weights, 1.0, 0.0, 1.0);

    // No operators = zero energy
    assert_eq!(engine.compute_energy(), 0.0);

    // Add some operators
    engine.op_seq.n_operators = 8;

    // E = -8 / (1.0 * 4) = -2.0
    assert!((engine.compute_energy() - (-2.0)).abs() < 1e-10);
}

#[test]
fn test_loop_update_empty() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 0.25);
    let mut engine = SSEEngine::new(lattice, hs, 100, weights, 1.0, 0.0, 1.0);

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    // Should not crash on empty sequence
    engine.loopupdate(&mut rng);

    assert_eq!(engine.op_seq.n_operators, 0);
}

#[test]
fn test_loop_update_with_operators() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 0.25);
    let mut engine = SSEEngine::new(lattice, hs, 100, weights, 1.0, 0.0, 1.0);

    // Insert some operators
    engine.op_seq.vertices[0] = Vertex {
        bond_idx: 0,
        op: OpType::Diagonal,
        vertex_idx: 1, // Diagonal ↑↑→↑↑
    };
    engine.op_seq.vertices[1] = Vertex {
        bond_idx: 1,
        op: OpType::OffDiagonal,
        vertex_idx: 5, // OffDiagonal ↑↓→↓↑
    };
    engine.op_seq.n_operators = 2;

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    // Should not crash
    engine.loopupdate(&mut rng);
}

// ============================================================================
// VertexData tests
// ============================================================================


#[test]
fn test_vertex_data_leg_states() {
    // Identity
    assert_eq!(VertexData::leg_states(0), [0, 0, 0, 0]);
    // Diagonal ↑↑→↑↑
    assert_eq!(VertexData::leg_states(1), [0, 0, 0, 0]);
    // Diagonal ↑↓→↑↓
    assert_eq!(VertexData::leg_states(2), [0, 1, 0, 1]);
    // Diagonal ↓↑→↓↑
    assert_eq!(VertexData::leg_states(3), [1, 0, 1, 0]);
    // Diagonal ↓↓→↓↓
    assert_eq!(VertexData::leg_states(4), [1, 1, 1, 1]);
    // OffDiagonal ↑↓→↓↑
    assert_eq!(VertexData::leg_states(5), [0, 1, 1, 0]);
    // OffDiagonal ↓↑→↑↓
    assert_eq!(VertexData::leg_states(6), [1, 0, 0, 1]);
}

#[test]
fn test_vertex_data_op_type() {
    assert_eq!(VertexData::op_type(0), OpType::Identity);
    for idx in 1..=4 {
        assert_eq!(VertexData::op_type(idx), OpType::Diagonal);
    }
    for idx in 5..=6 {
        assert_eq!(VertexData::op_type(idx), OpType::OffDiagonal);
    }
}

#[test]
fn test_vertex_data_diag_vertex() {
    assert_eq!(VertexData::diag_vertex(0, 0), 1);
    assert_eq!(VertexData::diag_vertex(0, 1), 2);
    assert_eq!(VertexData::diag_vertex(1, 0), 3);
    assert_eq!(VertexData::diag_vertex(1, 1), 4);
}

#[test]
fn test_vertex_data_offdiag_vertex() {
    assert_eq!(VertexData::offdiag_vertex(0, 1), 5);
    assert_eq!(VertexData::offdiag_vertex(1, 0), 6);
}

#[test]
fn test_vertex_data_scatter_exit_leg() {
    // Exit leg pairs: (0↔1, 2↔3) per Julia's xor(leg-1,1)+1 (1-indexed)
    for leg in 0..4 {
        let expected = leg ^ 1;
        let (leg_out, _) = VertexData::scatter(leg, 1);
        assert_eq!(leg_out, expected, "leg {} should exit on {}", leg, expected);
    }
}

#[test]
fn test_vertex_data_scatter_diag_to_offdiag() {
    // Vertex 2 (↑↓→↑↓): flip legs 0,1 → [1,0,0,1] → OffDiag 6 (↓↑→↑↓)
    let (_, new_idx) = VertexData::scatter(0, 2);
    assert_eq!(new_idx, 6);

    // Vertex 3 (↓↑→↓↑): flip legs 0,1 → [0,1,1,0] → OffDiag 5 (↑↓→↓↑)
    let (_, new_idx) = VertexData::scatter(0, 3);
    assert_eq!(new_idx, 5);
}

#[test]
fn test_vertex_data_scatter_offdiag_to_diag() {
    // OffDiag 5 (↑↓→↓↑): flip legs 0,1 → [1,0,1,0] → Diag 3 (↓↑→↓↑)
    let (_, new_idx) = VertexData::scatter(0, 5);
    assert_eq!(new_idx, 3);

    // OffDiag 6 (↓↑→↑↓): flip legs 0,1 → [0,1,0,1] → Diag 2 (↑↓→↑↓)
    let (_, new_idx) = VertexData::scatter(0, 6);
    assert_eq!(new_idx, 2);
}

#[test]
fn test_vertex_data_weight() {
    let j = 1.0;
    // Diagonal: J/4
    for idx in 1..=4 {
        assert!((VertexData::weight(idx, j) - 0.25).abs() < 1e-10);
    }
    // OffDiagonal: J/2
    for idx in 5..=6 {
        assert!((VertexData::weight(idx, j) - 0.5).abs() < 1e-10);
    }
}

// ============================================================================

// VertexList tests
// ============================================================================


#[test]
fn test_vertex_list_chain() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 1.0);
    let mut engine = SSEEngine::new(lattice, hs, 100, weights, 1.0, 0.0, 1.0);

    // Insert diagonal operators on bonds
    // Bond 1: sites 0-1, Bond 2: sites 1-2, Bond 3: sites 2-3
    engine.op_seq.vertices[0] = Vertex { bond_idx: 1, op: OpType::Diagonal, vertex_idx: 2 }; // ↑↓→↑↓
    engine.op_seq.vertices[1] = Vertex { bond_idx: 2, op: OpType::Diagonal, vertex_idx: 2 }; // ↑↓→↑↓
    engine.spins = vec![0, 1, 0, 1]; // ↑↓↑↓
    engine.op_seq.n_operators = 2;

    // Build vertex list
    let mut vertex_list = VertexList::new(engine.lattice.n_sites, engine.op_seq.max_length);
    vertex_list.build(&engine.op_seq, &engine.bond_list);

    // Sites 0, 1, 2 are involved in the operators (bonds 0-1 and 1-2)
    // Site 3 is NOT involved, so it should have sentinel v_first
    for site in 0..3 {
        let (_, pos) = vertex_list.v_first(site);
        assert!(pos != usize::MAX, "Site {} should have v_first", site);
    }
    // Site 3 should have sentinel (no operators on its worldline)
    let (_, pos) = vertex_list.v_first(3);
    assert_eq!(pos, usize::MAX, "Site 3 should NOT have v_first (no operators)");
}

#[test]
fn test_vertex_list_empty() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 1.0);
    let engine = SSEEngine::new(lattice, hs, 100, weights, 1.0, 0.0, 1.0);

    let mut vertex_list = VertexList::new(engine.lattice.n_sites, engine.op_seq.max_length);
    vertex_list.build(&engine.op_seq, &engine.bond_list);

    // All sites should have sentinel v_first (no operators)
    for site in 0..4 {
        let (_, pos) = vertex_list.v_first(site);
        assert_eq!(pos, usize::MAX, "Site {} should have no v_first in empty VL", site);
    }
}
