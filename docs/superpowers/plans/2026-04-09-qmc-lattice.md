# QMC.rs Lattice Module Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement qmc_lattice module with SSE + Directed Loop algorithm for lattice quantum Monte Carlo simulations.

**Architecture:** Layered design following philosophy: topology layer (Lattice adjacency list), physics layer (HilbertSpace trait), algorithm layer (SSEEngine generic over HilbertSpace). Users implement `bond_operators()`, framework handles operator sequence and updates.

**Tech Stack:** Rust, rand_xoshiro, Carlo.rs framework (MonteCarlo trait), HashMap for bond weights

---

## File Structure

```
QMC.rs/
├── Cargo.toml
├── src/
│   ├── lib.rs                    # Module exports, re-export from Carlo.rs
│   ├── lattice/
│   │   ├── mod.rs                # Lattice, Neighbor structs
│   │   ├── bond.rs               # BondType enum
│   │   └── builders.rs           # build_chain, build_square, etc.
│   ├── hilbert/
│   │   ├── mod.rs                # HilbertSpace trait, LocalState, OpType
│   │   └── spin_half.rs          # SpinHalfHS implementation
│   ├── sse/
│   │   ├── mod.rs                # SSEMonteCarlo trait, SSECore
│   │   ├── engine.rs             # SSEEngine, OperatorSequence, Vertex
│   │   ├── diagonal.rs           # diagonal_update implementation
│   │   ├── loop.rs               # loop_update (Directed Loop) implementation
│   │   └── measurements.rs       # compute_energy, compute_magnetization
│   └── models/
│       ├── mod.rs                # Model exports
│       └── heisenberg.rs         # HeisenbergModel example
└── tests/
    ├── lattice_test.rs           # Lattice and builders tests
    ├── hilbert_test.rs           # HilbertSpace tests
    ├── sse_engine_test.rs        # SSEEngine structure tests
    ├── heisenberg_test.rs        # Heisenberg model integration test
    └── physics_test.rs           # Physical correctness tests (Bethe ansatz)
```

---

## Phase A: Project Setup and Core Types

### Task 1: Create QMC.rs Cargo.toml

**Files:**
- Create: `QMC.rs/Cargo.toml`

- [ ] **Step 1: Write Cargo.toml**

```toml
[package]
name = "qmc-rs"
version = "0.1.0"
edition = "2021"
authors = ["Scuttle Team"]
description = "Quantum Monte Carlo algorithm toolbox (SSE, Directed Loop)"

[dependencies]
carlo-rs = { path = "../Carlo.rs" }
rand = "0.8"
rand_xoshiro = "0.8"

[dev-dependencies]
approx = "0.5"

[features]
default = []
hdf5 = ["carlo-rs/hdf5"]
mpi = ["carlo-rs/mpi"]
```

- [ ] **Step 2: Create src directory structure**

```bash
mkdir -p QMC.rs/src/lattice QMC.rs/src/hilbert QMC.rs/src/sse QMC.rs/src/models QMC.rs/tests
touch QMC.rs/src/lib.rs
touch QMC.rs/src/lattice/mod.rs QMC.rs/src/lattice/bond.rs QMC.rs/src/lattice/builders.rs
touch QMC.rs/src/hilbert/mod.rs QMC.rs/src/hilbert/spin_half.rs
touch QMC.rs/src/sse/mod.rs QMC.rs/src/sse/engine.rs QMC.rs/src/sse/diagonal.rs QMC.rs/src/sse/loop.rs QMC.rs/src/sse/measurements.rs
touch QMC.rs/src/models/mod.rs QMC.rs/src/models/heisenberg.rs
touch QMC.rs/tests/lattice_test.rs
```

- [ ] **Step 3: Verify build works**

```bash
cd QMC.rs && cargo check
```
Expected: Compiles successfully (empty modules)

- [ ] **Step 4: Commit**

```bash
git add QMC.rs/
git commit -m "feat(qmc): initialize QMC.rs project structure"
```

---

### Task 2: Define BondType Enum

**Files:**
- Create: `QMC.rs/src/lattice/bond.rs`
- Test: `QMC.rs/tests/lattice_test.rs`

- [ ] **Step 1: Write the failing test**

```rust
// tests/lattice_test.rs
use qmc_rs::lattice::BondType;

#[test]
fn test_bond_type_variants() {
    // 1D chain
    assert_eq!(BondType::ChainX as u8, 0);

    // 2D square
    let bx = BondType::SquareX;
    let by = BondType::SquareY;
    assert_ne!(bx, by);

    // 2D triangular
    assert_eq!(BondType::TriX, BondType::TriX);

    // Custom
    let custom = BondType::Custom(42);
    assert_eq!(custom, BondType::Custom(42));
}

#[test]
fn test_bond_type_hashable() {
    use std::collections::HashMap;
    let mut weights: HashMap<BondType, f64> = HashMap::new();
    weights.insert(BondType::SquareX, 1.0);
    assert_eq!(weights.get(&BondType::SquareX), Some(&1.0));
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test bond_type --no-run
```
Expected: Compilation error "use of undeclared type `BondType`"

- [ ] **Step 3: Write minimal implementation**

```rust
// src/lattice/bond.rs
use std::hash::Hash;

/// Bond type enum for direction-dependent Hamiltonian parameters.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondType {
    /// 1D chain: horizontal bond
    ChainX,

    /// 2D square lattice
    SquareX,
    SquareY,

    /// 2D triangular lattice (0°, 60°, 120°)
    TriX,
    TriY,
    TriZ,

    /// 2D honeycomb lattice
    HoneyX,
    HoneyY,
    HoneyZ,

    /// Custom bond type for arbitrary networks
    Custom(u8),
}
```

- [ ] **Step 4: Export from mod.rs**

```rust
// src/lattice/mod.rs
mod bond;
mod builders;

pub use bond::BondType;
```

- [ ] **Step 5: Export from lib.rs**

```rust
// src/lib.rs
pub mod lattice;
pub mod hilbert;
pub mod sse;
pub mod models;

// Re-export Carlo.rs types for convenience
pub use carlo_rs::{
    MonteCarlo, Context, CarloError, FromParams, Params,
    Scheduler, RunConfig, RayonBackend, Results, Estimate,
};
```

- [ ] **Step 6: Run test to verify it passes**

```bash
cd QMC.rs && cargo test bond_type
```
Expected: 2 tests pass

- [ ] **Step 7: Commit**

```bash
git add QMC.rs/src/lattice/bond.rs QMC.rs/src/lattice/mod.rs QMC.rs/src/lib.rs QMC.rs/tests/lattice_test.rs
git commit -m "feat(qmc): define BondType enum for lattice bonds"
```

---

### Task 3: Define Neighbor and Lattice Structs

**Files:**
- Modify: `QMC.rs/src/lattice/mod.rs`
- Test: `QMC.rs/tests/lattice_test.rs`

- [ ] **Step 1: Write the failing test**

```rust
// tests/lattice_test.rs (append)
use qmc_rs::lattice::{Neighbor, Lattice};

#[test]
fn test_neighbor_struct() {
    let n = Neighbor {
        target: 5,
        bond_type: BondType::SquareX,
    };
    assert_eq!(n.target, 5);
    assert_eq!(n.bond_type, BondType::SquareX);
}

#[test]
fn test_lattice_basic() {
    let lattice = Lattice {
        sites: vec![
            vec![Neighbor { target: 1, bond_type: BondType::ChainX }],
            vec![Neighbor { target: 0, bond_type: BondType::ChainX }],
        ],
        n_sites: 2,
        n_bonds: 2,
    };
    assert_eq!(lattice.n_sites, 2);
    assert_eq!(lattice.n_bonds, 2);
    assert_eq!(lattice.sites[0].len(), 1);
}

#[test]
fn test_lattice_clone() {
    let lattice = Lattice {
        sites: vec![vec![]],
        n_sites: 1,
        n_bonds: 0,
    };
    let cloned = lattice.clone();
    assert_eq!(cloned.n_sites, 1);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test lattice_basic --no-run
```
Expected: Compilation error "use of undeclared type `Neighbor`"

- [ ] **Step 3: Write minimal implementation**

```rust
// src/lattice/mod.rs (replace content)
mod bond;
mod builders;

pub use bond::BondType;

/// Neighbor entry in adjacency list.
#[derive(Clone, Debug)]
pub struct Neighbor {
    /// Target site index
    pub target: usize,
    /// Bond type for direction-dependent weights
    pub bond_type: BondType,
}

/// Lattice topology represented as adjacency list.
#[derive(Clone, Debug)]
pub struct Lattice {
    /// Adjacency list: sites[i] = neighbors of site i
    pub sites: Vec<Vec<Neighbor>>,
    /// Total number of sites
    pub n_sites: usize,
    /// Total number of bonds (counting each bond once)
    pub n_bonds: usize,
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd QMC.rs && cargo test lattice_basic
```
Expected: 3 tests pass

- [ ] **Step 5: Commit**

```bash
git add QMC.rs/src/lattice/mod.rs QMC.rs/tests/lattice_test.rs
git commit -m "feat(qmc): define Neighbor and Lattice structs"
```

---

### Task 4: Define HilbertSpace Trait and OpType

**Files:**
- Create: `QMC.rs/src/hilbert/mod.rs`
- Test: `QMC.rs/tests/hilbert_test.rs`

- [ ] **Step 1: Create test file**

```rust
// tests/hilbert_test.rs
use qmc_rs::hilbert::{HilbertSpace, LocalState, OpType};

#[test]
fn test_op_type_variants() {
    assert_eq!(OpType::Identity, OpType::Identity);
    assert_ne!(OpType::Diagonal, OpType::OffDiagonal);
}

#[test]
fn test_local_state_type() {
    let state: LocalState = 0;
    assert_eq!(state, 0);
    let state: LocalState = 1;
    assert_eq!(state, 1);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test hilbert --no-run
```
Expected: Compilation error

- [ ] **Step 3: Write minimal implementation**

```rust
// src/hilbert/mod.rs
mod spin_half;

pub use spin_half::SpinHalfHS;

/// Local state encoding for lattice sites.
/// Spin-1/2: 0 = Up, 1 = Down
/// Hubbard: 0 = empty, 1 = up, 2 = down, 3 = double
pub type LocalState = u8;

/// Operator type in SSE representation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OpType {
    /// Identity operator (empty vertex)
    Identity,
    /// Diagonal operator (e.g., SzSz, n_i n_j)
    Diagonal,
    /// Off-diagonal operator (e.g., S+S-, hopping)
    OffDiagonal,
}

/// Trait defining Hilbert space rules for operator actions.
pub trait HilbertSpace: Clone {
    /// Local Hilbert space dimension per site.
    fn local_dim(&self) -> usize;

    /// Check if operator is allowed given local states.
    /// states: [source_state, target_state] for bond operators
    fn is_allowed(&self, states: &[LocalState], op: &OpType) -> bool;

    /// Apply operator to local states (in-place modification).
    fn apply(&self, states: &mut [LocalState], op: &OpType);

    /// Compute dimensionless diagonal matrix element.
    /// Returns pure numerical part; engine multiplies by coupling constant.
    fn diagonal_element(&self, states: &[LocalState], op: &OpType) -> f64;
}
```

- [ ] **Step 4: Create stub spin_half.rs**

```rust
// src/hilbert/spin_half.rs
use crate::hilbert::{HilbertSpace, LocalState, OpType};

/// Spin-1/2 Hilbert space (Ising, Heisenberg, XXZ models).
#[derive(Clone, Copy, Debug, Default)]
pub struct SpinHalfHS;

impl HilbertSpace for SpinHalfHS {
    fn local_dim(&self) -> usize {
        2
    }

    fn is_allowed(&self, states: &[LocalState], op: &OpType) -> bool {
        match op {
            OpType::Identity => true,
            OpType::Diagonal => true,
            OpType::OffDiagonal => states[0] != states[1],
        }
    }

    fn apply(&self, states: &mut [LocalState], op: &OpType) {
        if *op == OpType::OffDiagonal {
            states[0] ^= 1;
            states[1] ^= 1;
        }
    }

    fn diagonal_element(&self, states: &[LocalState], op: &OpType) -> f64 {
        if *op == OpType::Diagonal {
            let s1 = if states[0] == 0 { 0.5 } else { -0.5 };
            let s2 = if states[1] == 0 { 0.5 } else { -0.5 };
            s1 * s2
        } else {
            0.0
        }
    }
}
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd QMC.rs && cargo test hilbert
```
Expected: 2 tests pass

- [ ] **Step 6: Add SpinHalfHS tests**

```rust
// tests/hilbert_test.rs (append)
use qmc_rs::hilbert::SpinHalfHS;

#[test]
fn test_spin_half_local_dim() {
    let hs = SpinHalfHS;
    assert_eq!(hs.local_dim(), 2);
}

#[test]
fn test_spin_half_is_allowed() {
    let hs = SpinHalfHS;

    // Diagonal always allowed
    assert!(hs.is_allowed(&[0, 0], &OpType::Diagonal));
    assert!(hs.is_allowed(&[0, 1], &OpType::Diagonal));
    assert!(hs.is_allowed(&[1, 0], &OpType::Diagonal));
    assert!(hs.is_allowed(&[1, 1], &OpType::Diagonal));

    // Off-diagonal only for antiparallel spins
    assert!(!hs.is_allowed(&[0, 0], &OpType::OffDiagonal));  // ↑↑ - not allowed
    assert!(hs.is_allowed(&[0, 1], &OpType::OffDiagonal));   // ↑↓ - allowed
    assert!(hs.is_allowed(&[1, 0], &OpType::OffDiagonal));   // ↓↑ - allowed
    assert!(!hs.is_allowed(&[1, 1], &OpType::OffDiagonal));  // ↓↓ - not allowed
}

#[test]
fn test_spin_half_apply() {
    let hs = SpinHalfHS;

    // Diagonal does nothing
    let mut states = [0, 1];
    hs.apply(&mut states, &OpType::Diagonal);
    assert_eq!(states, [0, 1]);

    // Off-diagonal flips both spins
    let mut states = [0, 1];  // ↑↓
    hs.apply(&mut states, &OpType::OffDiagonal);
    assert_eq!(states, [1, 0]);  // ↓↑

    let mut states = [1, 0];  // ↓↑
    hs.apply(&mut states, &OpType::OffDiagonal);
    assert_eq!(states, [0, 1]);  // ↑↓
}

#[test]
fn test_spin_half_diagonal_element() {
    let hs = SpinHalfHS;

    // ↑↑: (+1/2)(+1/2) = 0.25
    assert_eq!(hs.diagonal_element(&[0, 0], &OpType::Diagonal), 0.25);

    // ↑↓: (+1/2)(-1/2) = -0.25
    assert_eq!(hs.diagonal_element(&[0, 1], &OpType::Diagonal), -0.25);

    // ↓↑: (-1/2)(+1/2) = -0.25
    assert_eq!(hs.diagonal_element(&[1, 0], &OpType::Diagonal), -0.25);

    // ↓↓: (-1/2)(-1/2) = 0.25
    assert_eq!(hs.diagonal_element(&[1, 1], &OpType::Diagonal), 0.25);

    // Off-diagonal returns 0
    assert_eq!(hs.diagonal_element(&[0, 1], &OpType::OffDiagonal), 0.0);
}
```

- [ ] **Step 7: Run all hilbert tests**

```bash
cd QMC.rs && cargo test hilbert
```
Expected: 7 tests pass

- [ ] **Step 8: Commit**

```bash
git add QMC.rs/src/hilbert/ QMC.rs/tests/hilbert_test.rs
git commit -m "feat(qmc): define HilbertSpace trait and SpinHalfHS implementation"
```

---

## Phase B: Geometry Builders

### Task 5: Implement build_chain

**Files:**
- Modify: `QMC.rs/src/lattice/builders.rs`
- Test: `QMC.rs/tests/lattice_test.rs`

- [ ] **Step 1: Write the failing test**

```rust
// tests/lattice_test.rs (append)
use qmc_rs::lattice::builders::build_chain;

#[test]
fn test_build_chain_open() {
    let lattice = build_chain(4, false);  // 4 sites, open boundary

    assert_eq!(lattice.n_sites, 4);
    assert_eq!(lattice.n_bonds, 3);  // N-1 bonds for open chain

    // Site 0 has 1 neighbor (site 1)
    assert_eq!(lattice.sites[0].len(), 1);
    assert_eq!(lattice.sites[0][0].target, 1);

    // Site 1 has 2 neighbors (sites 0 and 2)
    assert_eq!(lattice.sites[1].len(), 2);

    // Site 3 has 1 neighbor (site 2)
    assert_eq!(lattice.sites[3].len(), 1);
}

#[test]
fn test_build_chain_periodic() {
    let lattice = build_chain(4, true);  // 4 sites, periodic boundary

    assert_eq!(lattice.n_sites, 4);
    assert_eq!(lattice.n_bonds, 4);  // N bonds for periodic chain

    // Every site has 2 neighbors
    for i in 0..4 {
        assert_eq!(lattice.sites[i].len(), 2);
    }
}

#[test]
fn test_build_chain_bond_types() {
    let lattice = build_chain(10, true);

    // All bonds should be ChainX
    for site in &lattice.sites {
        for neighbor in site {
            assert_eq!(neighbor.bond_type, BondType::ChainX);
        }
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test build_chain --no-run
```
Expected: Compilation error "use of undeclared function `build_chain`"

- [ ] **Step 3: Write minimal implementation**

```rust
// src/lattice/builders.rs
use crate::lattice::{BondType, Neighbor, Lattice};

/// Build a 1D chain lattice.
///
/// # Arguments
/// * `n_sites` - Number of sites in the chain
/// * `pbc` - Periodic boundary condition (true = ring, false = open chain)
pub fn build_chain(n_sites: usize, pbc: bool) -> Lattice {
    assert!(n_sites >= 2, "Chain must have at least 2 sites");

    let n_bonds = if pbc { n_sites } else { n_sites - 1 };

    let mut sites = Vec::with_capacity(n_sites);

    for i in 0..n_sites {
        let mut neighbors = Vec::new();

        // Left neighbor (i-1)
        if i > 0 || pbc {
            let left = if i > 0 { i - 1 } else { n_sites - 1 };
            neighbors.push(Neighbor {
                target: left,
                bond_type: BondType::ChainX,
            });
        }

        // Right neighbor (i+1)
        if i < n_sites - 1 || pbc {
            let right = if i < n_sites - 1 { i + 1 } else { 0 };
            neighbors.push(Neighbor {
                target: right,
                bond_type: BondType::ChainX,
            });
        }

        sites.push(neighbors);
    }

    Lattice {
        sites,
        n_sites,
        n_bonds,
    }
}
```

- [ ] **Step 4: Export from mod.rs**

```rust
// src/lattice/mod.rs (update)
mod bond;
mod builders;

pub use bond::BondType;
pub use builders::{build_chain};
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd QMC.rs && cargo test build_chain
```
Expected: 3 tests pass

- [ ] **Step 6: Commit**

```bash
git add QMC.rs/src/lattice/builders.rs QMC.rs/src/lattice/mod.rs QMC.rs/tests/lattice_test.rs
git commit -m "feat(qmc): implement build_chain for 1D lattice"
```

---

### Task 6: Implement build_square

**Files:**
- Modify: `QMC.rs/src/lattice/builders.rs`
- Test: `QMC.rs/tests/lattice_test.rs`

- [ ] **Step 1: Write the failing test**

```rust
// tests/lattice_test.rs (append)
use qmc_rs::lattice::builders::build_square;

#[test]
fn test_build_square_basic() {
    let lattice = build_square(4, 4, true);  // 4x4 square, periodic

    assert_eq!(lattice.n_sites, 16);
    // 16 sites × 4 bonds per site / 2 (each bond counted twice) = 32
    // But n_bonds counts each bond once, so 32 for periodic
    assert_eq!(lattice.n_bonds, 32);
}

#[test]
fn test_build_square_open() {
    let lattice = build_square(3, 3, false);  // 3x3, open

    assert_eq!(lattice.n_sites, 9);
    // Open: horizontal bonds = (3-1)*3 = 6, vertical bonds = 3*(3-1) = 6
    assert_eq!(lattice.n_bonds, 12);
}

#[test]
fn test_build_square_bond_types() {
    let lattice = build_square(2, 2, true);

    // Check that X and Y bond types exist
    let mut has_x = false;
    let mut has_y = false;

    for site in &lattice.sites {
        for neighbor in site {
            if neighbor.bond_type == BondType::SquareX {
                has_x = true;
            }
            if neighbor.bond_type == BondType::SquareY {
                has_y = true;
            }
        }
    }

    assert!(has_x, "Should have SquareX bonds");
    assert!(has_y, "Should have SquareY bonds");
}

#[test]
fn test_build_square_neighbors() {
    let lattice = build_square(4, 4, true);

    // Site at (0, 0) should have 4 neighbors
    assert_eq!(lattice.sites[0].len(), 4);

    // Every site in periodic should have 4 neighbors
    for i in 0..16 {
        assert_eq!(lattice.sites[i].len(), 4);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test build_square --no-run
```
Expected: Compilation error

- [ ] **Step 3: Write minimal implementation**

```rust
// src/lattice/builders.rs (append)
/// Build a 2D square lattice.
///
/// # Arguments
/// * `lx` - Number of sites in x direction
/// * `ly` - Number of sites in y direction
/// * `pbc` - Periodic boundary condition
pub fn build_square(lx: usize, ly: usize, pbc: bool) -> Lattice {
    assert!(lx >= 2 && ly >= 2, "Square lattice must have at least 2x2 sites");

    let n_sites = lx * ly;

    // Count bonds: horizontal bonds + vertical bonds
    let h_bonds = if pbc { lx * ly } else { (lx - 1) * ly };
    let v_bonds = if pbc { lx * ly } else { lx * (ly - 1) };
    let n_bonds = h_bonds + v_bonds;

    let mut sites = Vec::with_capacity(n_sites);

    for y in 0..ly {
        for x in 0..lx {
            let i = y * lx + x;
            let mut neighbors = Vec::new();

            // X-direction neighbor (right)
            if x < lx - 1 || pbc {
                let x_neighbor = if x < lx - 1 { i + 1 } else { y * lx };
                neighbors.push(Neighbor {
                    target: x_neighbor,
                    bond_type: BondType::SquareX,
                });
            }

            // Y-direction neighbor (down)
            if y < ly - 1 || pbc {
                let y_neighbor = if y < ly - 1 { i + lx } else { x };
                neighbors.push(Neighbor {
                    target: y_neighbor,
                    bond_type: BondType::SquareY,
                });
            }

            // X-direction neighbor (left) - for periodic
            if pbc && x == 0 {
                neighbors.push(Neighbor {
                    target: y * lx + lx - 1,
                    bond_type: BondType::SquareX,
                });
            }

            // Y-direction neighbor (up) - for periodic
            if pbc && y == 0 {
                neighbors.push(Neighbor {
                    target: (ly - 1) * lx + x,
                    bond_type: BondType::SquareY,
                });
            }

            sites.push(neighbors);
        }
    }

    Lattice {
        sites,
        n_sites,
        n_bonds,
    }
}
```

- [ ] **Step 4: Export from mod.rs**

```rust
// src/lattice/mod.rs (update exports)
pub use builders::{build_chain, build_square};
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd QMC.rs && cargo test build_square
```
Expected: 4 tests pass

- [ ] **Step 6: Commit**

```bash
git add QMC.rs/src/lattice/builders.rs QMC.rs/src/lattice/mod.rs QMC.rs/tests/lattice_test.rs
git commit -m "feat(qmc): implement build_square for 2D lattice"
```

---

## Phase C: SSE Engine Structure

### Task 7: Define Vertex and OperatorSequence

**Files:**
- Create: `QMC.rs/src/sse/engine.rs`
- Test: `QMC.rs/tests/sse_engine_test.rs`

- [ ] **Step 1: Create test file**

```rust
// tests/sse_engine_test.rs
use qmc_rs::sse::{Vertex, OperatorSequence, OpType};

#[test]
fn test_vertex_struct() {
    let v = Vertex {
        bond_idx: 5,
        op: OpType::Diagonal,
    };
    assert_eq!(v.bond_idx, 5);
    assert_eq!(v.op, OpType::Diagonal);
}

#[test]
fn test_operator_sequence_new() {
    let seq = OperatorSequence::new(100);

    assert_eq!(seq.max_length, 100);
    assert_eq!(seq.n_operators, 0);
    assert_eq!(seq.vertices.len(), 100);

    // All vertices should be Identity initially
    for v in &seq.vertices {
        assert_eq!(v.op, OpType::Identity);
    }
}

#[test]
fn test_operator_sequence_clone() {
    let seq = OperatorSequence::new(50);
    let cloned = seq.clone();
    assert_eq!(cloned.max_length, 50);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test sse_engine --no-run
```
Expected: Compilation error

- [ ] **Step 3: Write minimal implementation**

```rust
// src/sse/engine.rs
use crate::hilbert::OpType;

/// Vertex in SSE operator sequence.
#[derive(Clone, Debug)]
pub struct Vertex {
    /// Bond index in lattice
    pub bond_idx: usize,
    /// Operator type at this position
    pub op: OpType,
}

impl Default for Vertex {
    fn default() -> Self {
        Vertex {
            bond_idx: 0,
            op: OpType::Identity,
        }
    }
}

/// SSE operator sequence (operator string in imaginary time).
#[derive(Clone, Debug)]
pub struct OperatorSequence {
    /// Fixed-length array of vertices, filled with Identity
    pub vertices: Vec<Vertex>,
    /// Count of non-Identity operators
    pub n_operators: usize,
    /// Maximum sequence length (M = N_sites × β × factor)
    pub max_length: usize,
}

impl OperatorSequence {
    /// Create new empty operator sequence.
    pub fn new(max_length: usize) -> Self {
        let vertices = vec![Vertex::default(); max_length];
        OperatorSequence {
            vertices,
            n_operators: 0,
            max_length,
        }
    }
}
```

- [ ] **Step 4: Export from mod.rs**

```rust
// src/sse/mod.rs
mod engine;
mod diagonal;
mod loop_update;
mod measurements;

pub use engine::{Vertex, OperatorSequence};
pub use crate::hilbert::OpType;  // Re-export for convenience
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd QMC.rs && cargo test sse_engine
```
Expected: 3 tests pass

- [ ] **Step 6: Commit**

```bash
git add QMC.rs/src/sse/engine.rs QMC.rs/src/sse/mod.rs QMC.rs/tests/sse_engine_test.rs
git commit -m "feat(qmc): define Vertex and OperatorSequence for SSE"
```

---

### Task 8: Define SSEEngine Structure

**Files:**
- Modify: `QMC.rs/src/sse/engine.rs`
- Test: `QMC.rs/tests/sse_engine_test.rs`

- [ ] **Step 1: Write the failing test**

```rust
// tests/sse_engine_test.rs (append)
use qmc_rs::hilbert::SpinHalfHS;
use qmc_rs::lattice::{build_chain, BondType};
use qmc_rs::sse::SSEEngine;
use std::collections::HashMap;

#[test]
fn test_sse_engine_new() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;
    let max_length = 100;

    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 1.0);

    let engine = SSEEngine::new(lattice, hs, max_length, weights);

    assert_eq!(engine.lattice.n_sites, 4);
    assert_eq!(engine.spins.len(), 4);
    assert_eq!(engine.op_seq.max_length, 100);
}

#[test]
fn test_sse_engine_initial_spins() {
    let lattice = build_chain(10, true);
    let hs = SpinHalfHS;
    let engine = SSEEngine::new(lattice.clone(), hs, 100, HashMap::new());

    // All spins should be initialized to 0 (Up)
    for spin in &engine.spins {
        assert_eq!(*spin, 0);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test sse_engine_new --no-run
```
Expected: Compilation error "use of undeclared type `SSEEngine`"

- [ ] **Step 3: Write minimal implementation**

```rust
// src/sse/engine.rs (append)
use crate::hilbert::{HilbertSpace, LocalState};
use crate::lattice::{Lattice, BondType};
use std::collections::HashMap;

/// SSE engine with generic HilbertSpace for zero-cost abstraction.
pub struct SSEEngine<H: HilbertSpace> {
    /// Lattice topology
    pub lattice: Lattice,
    /// Current spin/particle configuration
    pub spins: Vec<LocalState>,
    /// Operator sequence
    pub op_seq: OperatorSequence,
    /// HilbertSpace implementation (static dispatch)
    pub hs: H,
    /// Coupling constants from bond_operators()
    pub weights: HashMap<BondType, f64>,
}

impl<H: HilbertSpace> SSEEngine<H> {
    /// Create new SSE engine.
    pub fn new(lattice: Lattice, hs: H, max_length: usize, weights: HashMap<BondType, f64>) -> Self {
        let n_sites = lattice.n_sites;
        let spins = vec![0; n_sites];  // Initialize all spins to Up (0)
        let op_seq = OperatorSequence::new(max_length);

        SSEEngine {
            lattice,
            spins,
            op_seq,
            hs,
            weights,
        }
    }
}
```

- [ ] **Step 4: Export from mod.rs**

```rust
// src/sse/mod.rs (update)
pub use engine::{Vertex, OperatorSequence, SSEEngine};
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd QMC.rs && cargo test sse_engine_new
```
Expected: 2 tests pass

- [ ] **Step 6: Commit**

```bash
git add QMC.rs/src/sse/engine.rs QMC.rs/src/sse/mod.rs QMC.rs/tests/sse_engine_test.rs
git commit -m "feat(qmc): define SSEEngine with generic HilbertSpace"
```

---

## Phase D: SSEMonteCarlo Trait and SSECore

### Task 9: Define LatticeQMC and SSEMonteCarlo Traits

**Files:**
- Modify: `QMC.rs/src/sse/mod.rs`
- Test: `QMC.rs/tests/heisenberg_test.rs`

- [ ] **Step 1: Create test file**

```rust
// tests/heisenberg_test.rs
use qmc_rs::{LatticeQMC, SSEMonteCarlo, HilbertSpace, OpType, BondType, Lattice};
use qmc_rs::hilbert::SpinHalfHS;

/// Minimal model for trait test
struct TestModel {
    lattice: Lattice,
    beta: f64,
}

impl LatticeQMC for TestModel {
    fn lattice(&self) -> &Lattice {
        &self.lattice
    }
}

impl SSEMonteCarlo for TestModel {
    type HilbertSpace = SpinHalfHS;

    fn bond_operators(&self, _bond_type: BondType) -> Vec<(OpType, f64)> {
        vec![(OpType::Diagonal, 1.0), (OpType::OffDiagonal, 0.5)]
    }

    fn hilbert_space(&self) -> &SpinHalfHS {
        static HS: SpinHalfHS = SpinHalfHS;
        &HS
    }

    fn beta(&self) -> f64 {
        self.beta
    }
}

#[test]
fn test_sse_monte_carlo_trait() {
    let lattice = qmc_rs::lattice::builders::build_chain(4, true);
    let model = TestModel { lattice, beta: 1.0 };

    assert_eq!(model.beta(), 1.0);
    assert_eq!(model.n_sites(), 4);

    let ops = model.bond_operators(BondType::ChainX);
    assert_eq!(ops.len(), 2);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test heisenberg --no-run
```
Expected: Compilation error "use of undeclared trait `LatticeQMC`"

- [ ] **Step 3: Write trait definitions**

```rust
// src/sse/mod.rs (replace content)
mod engine;
mod diagonal;
mod loop_update;
mod measurements;

pub use engine::{Vertex, OperatorSequence, SSEEngine};
pub use crate::hilbert::OpType;

use crate::{MonteCarlo, Context};
use crate::hilbert::HilbertSpace;
use crate::lattice::{Lattice, BondType};
use rand_xoshiro::Xoshiro256PlusPlus;
use rand_core::Rng;

/// Domain trait for lattice QMC methods.
pub trait LatticeQMC: MonteCarlo {
    /// Access lattice topology.
    fn lattice(&self) -> &Lattice;

    /// Total number of sites.
    fn n_sites(&self) -> usize {
        self.lattice().n_sites
    }
}

/// Method trait for SSE Monte Carlo.
pub trait SSEMonteCarlo: LatticeQMC {
    /// Associated HilbertSpace type.
    type HilbertSpace: HilbertSpace;

    /// Define operators on each bond type.
    /// Returns (OpType, coupling_constant) pairs.
    fn bond_operators(&self, bond_type: BondType) -> Vec<(OpType, f64)>;

    /// Access HilbertSpace implementation.
    fn hilbert_space(&self) -> &Self::HilbertSpace;

    /// Simulation inverse temperature.
    fn beta(&self) -> f64;
}
```

- [ ] **Step 4: Export from lib.rs**

```rust
// src/lib.rs (update)
pub mod lattice;
pub mod hilbert;
pub mod sse;
pub mod models;

pub use carlo_rs::{
    MonteCarlo, Context, CarloError, FromParams, Params,
    Scheduler, RunConfig, RayonBackend, Results, Estimate,
};

pub use lattice::{Lattice, Neighbor, BondType};
pub use hilbert::{HilbertSpace, LocalState, OpType, SpinHalfHS};
pub use sse::{LatticeQMC, SSEMonteCarlo, SSEEngine, Vertex, OperatorSequence};
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd QMC.rs && cargo test heisenberg
```
Expected: 1 test pass

- [ ] **Step 6: Commit**

```bash
git add QMC.rs/src/sse/mod.rs QMC.rs/src/lib.rs QMC.rs/tests/heisenberg_test.rs
git commit -m "feat(qmc): define LatticeQMC and SSEMonteCarlo traits"
```

---

### Task 10: Define SSECore Wrapper

**Files:**
- Modify: `QMC.rs/src/sse/mod.rs`
- Test: `QMC.rs/tests/heisenberg_test.rs`

- [ ] **Step 1: Write the failing test**

```rust
// tests/heisenberg_test.rs (append)
use qmc_rs::sse::SSECore;
use qmc_rs::hilbert::SpinHalfHS;

#[test]
fn test_sse_core_new() {
    let lattice = qmc_rs::lattice::builders::build_chain(4, true);
    let model = TestModel { lattice, beta: 1.0 };

    let core = SSECore::new(model);

    assert_eq!(core.engine.lattice.n_sites, 4);
    assert_eq!(core.engine.op_seq.max_length, 4 * 1 * 2 + 100);  // N*beta*factor + margin
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test sse_core_new --no-run
```
Expected: Compilation error

- [ ] **Step 3: Write minimal implementation**

```rust
// src/sse/mod.rs (append SSECore definition)
/// SSE core wrapper providing default MonteCarlo implementation.
pub struct SSECore<MC: SSEMonteCarlo> {
    /// SSE engine
    pub engine: SSEEngine<MC::HilbertSpace>,
    /// User model
    pub mc: MC,
}

impl<MC: SSEMonteCarlo> SSECore<MC> {
    /// Create SSE core from user model.
    pub fn new(mc: MC) -> Self {
        let lattice = mc.lattice().clone();
        let beta = mc.beta();
        let n_sites = lattice.n_sites;

        // Estimate max_length: M ~ N × β × (average operators per site)
        let max_length = (n_sites as f64 * beta * 2.0) as usize + 100;

        let hs = mc.hilbert_space().clone();

        // Build weights from bond_operators
        let weights = Self::build_weights(&mc, &lattice);

        let engine = SSEEngine::new(lattice, hs, max_length, weights);

        SSECore { engine, mc }
    }

    /// Collect bond weights from model.
    fn build_weights(mc: &MC, lattice: &Lattice) -> HashMap<BondType, f64> {
        let mut weights = HashMap::new();

        // Collect unique bond types from lattice
        for site in &lattice.sites {
            for neighbor in site {
                if !weights.contains_key(&neighbor.bond_type) {
                    let ops = mc.bond_operators(neighbor.bond_type);
                    // Sum weights for diagonal operators (used in diagonal update)
                    let diag_weight = ops.iter()
                        .filter(|(op, _)| *op == OpType::Diagonal)
                        .map(|(_, w)| w)
                        .sum();
                    weights.insert(neighbor.bond_type, diag_weight);
                }
            }
        }

        weights
    }
}
```

- [ ] **Step 4: Export from mod.rs**

```rust
// src/sse/mod.rs (update exports)
pub use engine::{Vertex, OperatorSequence, SSEEngine};
pub use crate::hilbert::OpType;

// Add SSECore to exports
pub use self::SSECore;
```

Wait, SSECore is defined in the same file. Let me fix the structure properly.

Actually, let me reorganize the module structure for clarity:

- [ ] **Step 5: Fix module structure - move SSECore to separate file**

```rust
// src/sse/core.rs (new file)
use super::{SSEEngine, SSEMonteCarlo, OpType};
use crate::hilbert::HilbertSpace;
use crate::lattice::{Lattice, BondType};
use std::collections::HashMap;

/// SSE core wrapper providing default MonteCarlo implementation.
pub struct SSECore<MC: SSEMonteCarlo> {
    /// SSE engine
    pub engine: SSEEngine<MC::HilbertSpace>,
    /// User model
    pub mc: MC,
}

impl<MC: SSEMonteCarlo> SSECore<MC> {
    /// Create SSE core from user model.
    pub fn new(mc: MC) -> Self {
        let lattice = mc.lattice().clone();
        let beta = mc.beta();
        let n_sites = lattice.n_sites;

        // Estimate max_length: M ~ N × β × (average operators per site)
        let max_length = (n_sites as f64 * beta * 2.0) as usize + 100;

        let hs = mc.hilbert_space().clone();

        // Build weights from bond_operators
        let weights = Self::build_weights(&mc, &lattice);

        let engine = SSEEngine::new(lattice, hs, max_length, weights);

        SSECore { engine, mc }
    }

    /// Collect bond weights from model.
    fn build_weights(mc: &MC, lattice: &Lattice) -> HashMap<BondType, f64> {
        let mut weights = HashMap::new();

        for site in &lattice.sites {
            for neighbor in site {
                if !weights.contains_key(&neighbor.bond_type) {
                    let ops = mc.bond_operators(neighbor.bond_type);
                    let diag_weight = ops.iter()
                        .filter(|(op, _)| *op == OpType::Diagonal)
                        .map(|(_, w)| w)
                        .sum();
                    weights.insert(neighbor.bond_type, diag_weight);
                }
            }
        }

        weights
    }
}
```

```rust
// src/sse/mod.rs (update)
mod engine;
mod core;
mod diagonal;
mod loop_update;
mod measurements;

pub use engine::{Vertex, OperatorSequence, SSEEngine};
pub use core::SSECore;
pub use crate::hilbert::OpType;

// Traits defined here
use crate::{MonteCarlo, Context};
use crate::hilbert::HilbertSpace;
use crate::lattice::{Lattice, BondType};

pub trait LatticeQMC: MonteCarlo {
    fn lattice(&self) -> &Lattice;
    fn n_sites(&self) -> usize {
        self.lattice().n_sites
    }
}

pub trait SSEMonteCarlo: LatticeQMC {
    type HilbertSpace: HilbertSpace;
    fn bond_operators(&self, bond_type: BondType) -> Vec<(OpType, f64)>;
    fn hilbert_space(&self) -> &Self::HilbertSpace;
    fn beta(&self) -> f64;
}
```

- [ ] **Step 6: Create core.rs file**

```bash
touch QMC.rs/src/sse/core.rs
```

- [ ] **Step 7: Run test to verify it passes**

```bash
cd QMC.rs && cargo test sse_core_new
```
Expected: 1 test pass

- [ ] **Step 8: Commit**

```bash
git add QMC.rs/src/sse/core.rs QMC.rs/src/sse/mod.rs QMC.rs/tests/heisenberg_test.rs
git commit -m "feat(qmc): define SSECore wrapper for MonteCarlo implementation"
```

---

## Phase E: Heisenberg Model Implementation

### Task 11: Implement HeisenbergModel

**Files:**
- Create: `QMC.rs/src/models/heisenberg.rs`
- Test: `QMC.rs/tests/heisenberg_test.rs`

- [ ] **Step 1: Write the failing test**

```rust
// tests/heisenberg_test.rs (append)
use qmc_rs::models::HeisenbergModel;

#[test]
fn test_heisenberg_model_new() {
    let lattice = qmc_rs::lattice::builders::build_chain(10, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0);  // beta=1.0, J=1.0

    assert_eq!(model.beta(), 1.0);
    assert_eq!(model.j(), 1.0);
    assert_eq!(model.n_sites(), 10);
}

#[test]
fn test_heisenberg_bond_operators() {
    let lattice = qmc_rs::lattice::builders::build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 0.5);  // J=0.5

    let ops = model.bond_operators(BondType::ChainX);

    // Should return [(Diagonal, J/4), (OffDiagonal, J/2)]
    assert_eq!(ops.len(), 2);

    // Find diagonal operator
    let diag = ops.iter().find(|(op, _)| *op == OpType::Diagonal);
    assert!(diag.is_some());
    assert_eq!(diag.unwrap().1, 0.5 * 0.25);  // J/4

    // Find off-diagonal operator
    let offdiag = ops.iter().find(|(op, _)| *op == OpType::OffDiagonal);
    assert!(offdiag.is_some());
    assert_eq!(offdiag.unwrap().1, 0.5 * 0.5);  // J/2
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test heisenberg_model --no-run
```
Expected: Compilation error

- [ ] **Step 3: Write implementation**

```rust
// src/models/heisenberg.rs
use crate::lattice::{Lattice, BondType};
use crate::hilbert::{HilbertSpace, OpType, SpinHalfHS};
use crate::sse::{LatticeQMC, SSEMonteCarlo};

/// Heisenberg model H = J Σ S_i · S_j
pub struct HeisenbergModel {
    lattice: Lattice,
    beta: f64,
    j: f64,
}

impl HeisenbergModel {
    /// Create new Heisenberg model.
    pub fn new(lattice: Lattice, beta: f64, j: f64) -> Self {
        HeisenbergModel { lattice, beta, j }
    }

    /// Access coupling constant J.
    pub fn j(&self) -> f64 {
        self.j
    }
}

impl LatticeQMC for HeisenbergModel {
    fn lattice(&self) -> &Lattice {
        &self.lattice
    }
}

impl SSEMonteCarlo for HeisenbergModel {
    type HilbertSpace = SpinHalfHS;

    fn bond_operators(&self, _bond_type: BondType) -> Vec<(OpType, f64)> {
        // H_b = J S_i · S_j = J (S^z_i S^z_j + 1/2 (S^+_i S^-_j + S^-_i S^+_j))
        // SSE decomposition:
        // - Diagonal: J S^z S^z → weight = J/4 (matrix element s1*s2 × J)
        // - Off-diagonal: J/2 (S^+ S^-) → weight = J/2
        vec![
            (OpType::Diagonal, self.j * 0.25),
            (OpType::OffDiagonal, self.j * 0.5),
        ]
    }

    fn hilbert_space(&self) -> &SpinHalfHS {
        static HS: SpinHalfHS = SpinHalfHS;
        &HS
    }

    fn beta(&self) -> f64 {
        self.beta
    }
}
```

- [ ] **Step 4: Export from models/mod.rs**

```rust
// src/models/mod.rs
mod heisenberg;

pub use heisenberg::HeisenbergModel;
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd QMC.rs && cargo test heisenberg_model
```
Expected: 2 tests pass

- [ ] **Step 6: Commit**

```bash
git add QMC.rs/src/models/heisenberg.rs QMC.rs/src/models/mod.rs QMC.rs/tests/heisenberg_test.rs
git commit -m "feat(qmc): implement HeisenbergModel for spin-1/2 chain"
```

---

### Task 12: Implement MonteCarlo for SSECore

**Files:**
- Modify: `QMC.rs/src/sse/core.rs`
- Test: `QMC.rs/tests/heisenberg_test.rs`

- [ ] **Step 1: Write the failing test**

```rust
// tests/heisenberg_test.rs (append)
use qmc_rs::{MonteCarlo, Context, SSECore};
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn test_sse_core_monte_carlo_trait() {
    let lattice = qmc_rs::lattice::builders::build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0);
    let core = SSECore::new(model);

    // SSECore should implement MonteCarlo
    // Check that it has the Rng type
    let _: <SSECore<HeisenbergModel> as MonteCarlo>::Rng = Xoshiro256PlusPlus::seed_from_u64(42);
}

#[test]
fn test_sse_core_sweep_basic() {
    let lattice = qmc_rs::lattice::builders::build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0);
    let mut core = SSECore::new(model);

    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 1000);

    // Sweep should not crash
    core.sweep(&mut ctx);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test sse_core_sweep --no-run
```
Expected: Compilation error (MonteCarlo not implemented for SSECore)

- [ ] **Step 3: Write MonteCarlo implementation**

```rust
// src/sse/core.rs (append)
use crate::{MonteCarlo, Context};
use rand_xoshiro::Xoshiro256PlusPlus;
use rand_core::Rng;

impl<MC: SSEMonteCarlo> MonteCarlo for SSECore<MC> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // Placeholder: will implement diagonal_update and loop_update later
        // For now, just advance sweep count
        ctx.advance_sweep();
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        // Placeholder: will implement measurements later
        // For now, measure nothing
    }

    fn name(&self) -> &'static str {
        "SSECore"
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd QMC.rs && cargo test sse_core_sweep
```
Expected: 2 tests pass

- [ ] **Step 5: Commit**

```bash
git add QMC.rs/src/sse/core.rs QMC.rs/tests/heisenberg_test.rs
git commit -m "feat(qmc): implement MonteCarlo trait for SSECore"
```

---

## Phase F: SSE Algorithms (diagonal_update)

### Task 13: Implement diagonal_update Skeleton

**Files:**
- Create: `QMC.rs/src/sse/diagonal.rs`
- Test: `QMC.rs/tests/sse_engine_test.rs`

- [ ] **Step 1: Write the failing test**

```rust
// tests/sse_engine_test.rs (append)
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn test_diagonal_update_basic() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;

    let mut weights = HashMap::new();
    weights.insert(BondType::ChainX, 1.0);

    let mut engine = SSEEngine::new(lattice, hs, 100, weights);

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

    // diagonal_update should not crash
    engine.diagonal_update(&mut rng);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test diagonal_update --no-run
```
Expected: Compilation error (method diagonal_update not found)

- [ ] **Step 3: Write placeholder implementation**

```rust
// src/sse/diagonal.rs
use crate::hilbert::{HilbertSpace, LocalState, OpType};
use crate::sse::{SSEEngine, OperatorSequence, Vertex};
use crate::lattice::{BondType, Neighbor};
use rand_core::Rng;
use std::collections::HashMap;

impl<H: HilbertSpace> SSEEngine<H> {
    /// Perform diagonal update on operator sequence.
    ///
    /// For each position in the sequence, decide:
    /// 1. If Identity: try to insert diagonal operator
    /// 2. If Diagonal: try to remove it
    /// 3. If OffDiagonal: check if still allowed, update states if needed
    pub fn diagonal_update<R: Rng>(&mut self, rng: &mut R) {
        // Placeholder: iterate through sequence and update operators
        // Full implementation in next task

        for i in 0..self.op_seq.max_length {
            let vertex = &self.op_seq.vertices[i];

            match vertex.op {
                OpType::Identity => {
                    // Try to insert diagonal operator
                    // TODO: implement insertion logic
                }
                OpType::Diagonal => {
                    // Try to remove diagonal operator
                    // TODO: implement removal logic
                }
                OpType::OffDiagonal => {
                    // Check and update states
                    // TODO: implement state propagation
                }
            }
        }
    }
}
```

- [ ] **Step 4: Export from mod.rs**

```rust
// src/sse/mod.rs (ensure diagonal is included)
mod diagonal;
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd QMC.rs && cargo test diagonal_update
```
Expected: 1 test pass (placeholder doesn't crash)

- [ ] **Step 6: Commit**

```bash
git add QMC.rs/src/sse/diagonal.rs QMC.rs/src/sse/mod.rs QMC.rs/tests/sse_engine_test.rs
git commit -m "feat(qmc): add diagonal_update skeleton for SSE"
```

---

## Phase G: Measurements

### Task 14: Implement compute_energy Skeleton

**Files:**
- Create: `QMC.rs/src/sse/measurements.rs`
- Test: `QMC.rs/tests/sse_engine_test.rs`

- [ ] **Step 1: Write the failing test**

```rust
// tests/sse_engine_test.rs (append)
#[test]
fn test_compute_energy_empty() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;

    let weights = HashMap::new();
    let engine = SSEEngine::new(lattice, hs, 100, weights);

    // Empty operator sequence should have zero energy contribution from operators
    let energy = engine.compute_energy();
    // Placeholder returns 0.0
    assert_eq!(energy, 0.0);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd QMC.rs && cargo test compute_energy --no-run
```
Expected: Compilation error

- [ ] **Step 3: Write placeholder implementation**

```rust
// src/sse/measurements.rs
use crate::hilbert::HilbertSpace;
use crate::sse::SSEEngine;

impl<H: HilbertSpace> SSEEngine<H> {
    /// Compute energy from operator sequence.
    ///
    /// Energy = -<n_operators> / beta for SSE representation
    pub fn compute_energy(&self) -> f64 {
        // Placeholder: will implement proper energy calculation
        // Energy in SSE: E = -M / β where M is average operator count
        0.0
    }

    /// Compute magnetization from spin configuration.
    pub fn compute_magnetization(&self) -> f64 {
        let n_sites = self.lattice.n_sites;

        let m: i32 = self.spins.iter()
            .map(|s| if *s == 0 { 1 } else { -1 })
            .sum();

        m.abs() as f64 / n_sites as f64
    }
}
```

- [ ] **Step 4: Export from mod.rs**

```rust
// src/sse/mod.rs (ensure measurements is included)
mod measurements;
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd QMC.rs && cargo test compute_energy
```
Expected: 1 test pass

- [ ] **Step 6: Add magnetization test**

```rust
// tests/sse_engine_test.rs (append)
#[test]
fn test_compute_magnetization() {
    let lattice = build_chain(4, true);
    let hs = SpinHalfHS;

    let mut engine = SSEEngine::new(lattice, hs, 100, HashMap::new());

    // All spins up (0)
    engine.spins = vec![0, 0, 0, 0];
    assert_eq!(engine.compute_magnetization(), 1.0);

    // All spins down (1)
    engine.spins = vec![1, 1, 1, 1];
    assert_eq!(engine.compute_magnetization(), 1.0);

    // Half up, half down
    engine.spins = vec![0, 1, 0, 1];
    assert_eq!(engine.compute_magnetization(), 0.0);
}
```

- [ ] **Step 7: Run all measurement tests**

```bash
cd QMC.rs && cargo test compute_magnetization
```
Expected: 1 test pass

- [ ] **Step 8: Commit**

```bash
git add QMC.rs/src/sse/measurements.rs QMC.rs/src/sse/mod.rs QMC.rs/tests/sse_engine_test.rs
git commit -m "feat(qmc): add compute_energy and compute_magnetization skeletons"
```

---

## Phase H: Integration Test

### Task 15: Full Integration Test with Carlo.rs

**Files:**
- Create: `QMC.rs/tests/integration_test.rs`

- [ ] **Step 1: Write integration test**

```rust
// tests/integration_test.rs
use qmc_rs::{
    MonteCarlo, Context, Scheduler, RunConfig, RayonBackend,
    SSECore, HeisenbergModel, lattice::builders::build_chain,
};
use rand_xoshiro::Xoshiro256PlusPlus;

#[test]
fn test_heisenberg_chain_with_carlo_rs() {
    // Build lattice
    let lattice = build_chain(4, true);

    // Create model
    let model = HeisenbergModel::new(lattice, 1.0, 1.0);  // beta=1.0, J=1.0

    // Create SSE core
    let core = SSECore::new(model);

    // Create scheduler
    let backend = RayonBackend::new(1);  // Single thread
    let config = RunConfig {
        thermalization_sweeps: 100,
        measurement_sweeps: 1000,
        binsize: 100,
        base_seed: 42,
        progress_interval: 0,  // Disable progress for test
        checkpoint_interval: 0,
    };
    let scheduler = Scheduler::new(backend, config);

    // Run simulation
    let params = qmc_rs::Params::new();
    let results = scheduler.run_one::<SSECore<HeisenbergModel>>(&params);

    // Basic checks
    assert!(results.is_ok());
}
```

- [ ] **Step 2: Run test to verify compilation**

```bash
cd QMC.rs && cargo test integration --no-run
```
Expected: Compilation successful (may fail at runtime if algorithms incomplete)

- [ ] **Step 3: Run test**

```bash
cd QMC.rs && cargo test integration
```
Expected: Test passes (sweep does nothing but advances counter)

- [ ] **Step 4: Commit**

```bash
git add QMC.rs/tests/integration_test.rs
git commit -m "test(qmc): add integration test with Carlo.rs framework"
```

---

## Phase I: Full Algorithm Implementation (Future Tasks)

The diagonal_update and loop_update algorithms require detailed physics logic. These are deferred to subsequent implementation cycles after the basic structure is validated.

**Future Task 16: Full diagonal_update Implementation**
- Insert/remove diagonal operators with proper probabilities
- Propagate spin states through off-diagonal operators
- Handle bond weights and HilbertSpace rules

**Future Task 17: Directed Loop Implementation**
- Build vertex graph for loop traversal
- Implement loop entering/exiting logic
- Update spin states during loop traversal

**Future Task 18: Energy Measurement Implementation**
- Calculate energy from operator density
- Implement proper binning for observables

**Future Task 19: Physical Validation Tests**
- Compare with Bethe ansatz for 1D Heisenberg chain
- Validate against Carlo.jl reference implementation

---

## Verification Checklist

After completing Tasks 1-15, verify:

```bash
cd QMC.rs

# All tests pass
cargo test

# No warnings
cargo clippy --all-targets -- -D warnings

# Documentation builds
cargo doc --no-deps

# Integration with Carlo.rs works
cargo test integration
```

---

*Implementation plan for QMC.rs lattice module. Follow TDD approach: write test first, implement minimal code to pass, commit frequently.*