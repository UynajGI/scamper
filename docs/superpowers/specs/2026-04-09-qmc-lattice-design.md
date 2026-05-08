# Lattice QMC Design Specification

> Design for QMC.rs `qmc_lattice` module implementing SSE + Directed Loop algorithm.
>
> Date: 2026-04-09

---

## I. Context and Philosophy Alignment

This design implements Phase 1 of the [QMC.rs Philosophy](../metaphysics/2026-04-09-qmc-rs-philosophy.md):

- **Position**: Algorithm toolbox (component library) between Carlo.rs framework and model implementations
- **Granularity**: Mid-level — users implement `bond_operators()`, framework handles operator sequence management
- **Trait hierarchy**: `MonteCarlo` → `LatticeQMC` → `SSEMonteCarlo`
- **Design principle**: Progressive complexity — simple defaults, advanced customization available

---

## II. Core Data Structures

### II.1 Topology Layer

Lattice is represented as an **adjacency list** — universal topology representation supporting arbitrary networks.

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondType {
    // 1D
    ChainX,
    // 2D Square
    SquareX, SquareY,
    // 2D Triangular
    TriX, TriY, TriZ,  // 0°, 60°, 120° directions
    // 2D Honeycomb
    HoneyX, HoneyY, HoneyZ,
    // Custom (for arbitrary networks)
    Custom(u8),
}

#[derive(Clone, Debug)]
pub struct Neighbor {
    pub target: usize,       // neighbor site index
    pub bond_type: BondType, // type tag for Hamiltonian weights
}

#[derive(Clone, Debug)]
pub struct Lattice {
    pub sites: Vec<Vec<Neighbor>>,  // adjacency list
    pub n_sites: usize,             // total number of sites
    pub n_bonds: usize,             // total number of bonds (for random selection)
}
```

**Design decisions**:
- Adjacency list stores topology only; bond weights are managed separately in engine
- `BondType` tags allow direction-dependent Hamiltonian parameters
- Geometry builders (`build_square()`, `build_triangular()`) construct adjacency from params

### II.2 Physics Layer

Hilbert space abstraction defines operator rules. Internal state uses compact encoding for performance; trait interface allows extensibility.

```rust
pub type LocalState = u8;  // compact encoding: spin 0/1, Hubbard 0-3

pub trait HilbertSpace {
    fn local_dim(&self) -> usize;

    /// Check if operator is allowed given local states
    /// states: [source_state, target_state] for bond operators
    fn is_allowed(&self, states: &[LocalState], op: &OpType) -> bool;

    /// Apply operator to local states (in-place modification)
    fn apply(&self, states: &mut [LocalState], op: &OpType);

    /// Compute dimensionless diagonal matrix element
    /// Returns pure numerical part; engine multiplies by coupling constant
    fn diagonal_element(&self, states: &[LocalState], op: &OpType) -> f64;
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OpType {
    Identity,
    Diagonal,    // e.g., SzSz, Hubbard interaction
    OffDiagonal, // e.g., S+S-, Hubbard hopping
}
```

**Convention**:
- `states[0]` = source site state, `states[1]` = target site state
- `diagonal_element` returns dimensionless value (e.g., SzSz returns s1*s2, not J*s1*s2)
- Coupling constants come from `bond_operators()` weights

### II.3 Algorithm Layer

SSE engine uses generic HilbertSpace parameter for zero-cost abstraction.

```rust
#[derive(Clone, Debug)]
pub struct Vertex {
    pub bond_idx: usize,  // index in bond list
    pub op: OpType,       // operator type at this position
}

pub struct OperatorSequence {
    pub vertices: Vec<Vertex>,  // fixed length = max_length, filled with Identity
    pub n_operators: usize,     // count of non-Identity operators
    pub max_length: usize,      // M = N_sites × β × (typical factor)
}

pub struct SSEEngine<H: HilbertSpace> {
    pub lattice: Lattice,
    pub spins: Vec<LocalState>,          // current spin/particle configuration
    pub op_seq: OperatorSequence,        // operator sequence
    pub hs: H,                           // HilbertSpace implementation (static dispatch)
    pub weights: HashMap<BondType, f64>, // coupling constants from bond_operators()
}
```

**Performance considerations**:
- `H: HilbertSpace` generic parameter eliminates dynamic dispatch overhead
- `vertices` fixed-length array avoids dynamic allocation during updates
- `HashMap<BondType, f64>` for weights; can optimize to `Vec<f64>` later if profiling shows bottleneck

---

## III. Trait Hierarchy

Following philosophy document's layered design:

```rust
// === Carlo.rs Base (existing) ===
pub trait MonteCarlo: Sized {
    type Rng: Rng + SeedableRng + Send;
    fn sweep(&mut self, ctx: &mut Context<Self::Rng>);
    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {}
}

// === QMC.rs Domain Layer ===
pub trait LatticeQMC: MonteCarlo {
    fn lattice(&self) -> &Lattice;
    fn n_sites(&self) -> usize { self.lattice().n_sites }
}

// === QMC.rs Method Layer ===
pub trait SSEMonteCarlo: LatticeQMC {
    type HilbertSpace: HilbertSpace;

    /// Define operators on each bond type
    /// Returns (OpType, coupling_constant) pairs
    fn bond_operators(&self, bond_type: BondType) -> Vec<(OpType, f64)>;

    /// Access HilbertSpace implementation
    fn hilbert_space(&self) -> &Self::HilbertSpace;

    /// Simulation parameters
    fn beta(&self) -> f64;
}
```

---

## IV. SSECore: Default Implementation

`SSECore` provides complete SSE algorithm, users only implement `bond_operators()`.

```rust
pub struct SSECore<MC: SSEMonteCarlo> {
    pub engine: SSEEngine<MC::HilbertSpace>,
    pub mc: MC,
}

impl<MC: SSEMonteCarlo> MonteCarlo for SSECore<MC> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.engine.diagonal_update(&mut ctx.rng);
        self.engine.loop_update(&mut ctx.rng);
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        // Default measurements
        let energy = self.engine.compute_energy();
        ctx.measure("Energy", energy);
    }
}

impl<MC: SSEMonteCarlo> SSECore<MC> {
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

    fn build_weights(mc: &MC, lattice: &Lattice) -> HashMap<BondType, f64> {
        // Collect all bond types and their operator weights
        // ...
    }
}
```

---

## V. Progressive Complexity Levels

### Level 1: Minimal Implementation

User implements `bond_operators()` only:

```rust
struct HeisenbergModel {
    lattice: Lattice,
    beta: f64,
    j: f64,
}

impl LatticeQMC for HeisenbergModel {
    fn lattice(&self) -> &Lattice { &self.lattice }
}

impl SSEMonteCarlo for HeisenbergModel {
    type HilbertSpace = HeisenbergHS;

    fn bond_operators(&self, _bond_type: BondType) -> Vec<(OpType, f64)> {
        vec![
            (OpType::Diagonal, self.j * 0.25),    // SzSz coupling
            (OpType::OffDiagonal, self.j * 0.5),  // S+S- coupling
        ]
    }

    fn hilbert_space(&self) -> &HeisenbergHS {
        static HS: HeisenbergHS = HeisenbergHS;
        &HS
    }

    fn beta(&self) -> f64 { self.beta }
}

// Usage
let model = HeisenbergModel::new(lattice, beta, j);
let core = SSECore::new(model);
scheduler.run_one::<SSECore<HeisenbergModel>>(&params);
```

### Level 2: Custom Observables

User overrides `measure()`:

```rust
impl MonteCarlo for SSECore<MyModel> {
    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        // Use default energy measurement
        let energy = self.engine.compute_energy();
        ctx.measure("Energy", energy);

        // Add custom observables
        let mag = self.engine.compute_magnetization();
        ctx.measure("Magnetization", mag);
    }
}
```

### Level 3: Replace SSECore

Advanced users can replace `SSECore` with custom implementation for special cases.

### Level 4: Full Manual Implementation

Implement complete `SSEMonteCarlo` trait manually for maximum control.

---

## VI. Default HilbertSpace Implementations

### VI.1 SpinHalf (Ising/Heisenberg/XXZ)

```rust
#[derive(Clone, Copy)]
pub struct SpinHalfHS;

impl HilbertSpace for SpinHalfHS {
    fn local_dim(&self) -> usize { 2 }

    fn is_allowed(&self, states: &[LocalState], op: &OpType) -> bool {
        match op {
            OpType::Identity => true,
            OpType::Diagonal => true,  // SzSz always allowed
            OpType::OffDiagonal => states[0] != states[1],  // antiparallel only
        }
    }

    fn apply(&self, states: &mut [LocalState], op: &OpType) {
        if *op == OpType::OffDiagonal {
            states[0] ^= 1;  // flip source spin
            states[1] ^= 1;  // flip target spin
        }
    }

    fn diagonal_element(&self, states: &[LocalState], op: &OpType) -> f64 {
        if *op == OpType::Diagonal {
            // State 0 = Up (+1/2), 1 = Down (-1/2)
            let s1 = if states[0] == 0 { 0.5 } else { -0.5 };
            let s2 = if states[1] == 0 { 0.5 } else { -0.5 };
            s1 * s2
        } else {
            0.0
        }
    }
}
```

### VI.2 Hubbard (Future)

Four local states: empty (0), up (1), down (2), double (3).

```rust
#[derive(Clone, Copy)]
pub struct HubbardHS {
    pub u: f64,  // interaction strength (for diagonal element)
}

impl HilbertSpace for HubbardHS {
    fn local_dim(&self) -> usize { 4 }
    // ... hopping and interaction operator rules
}
```

---

## VII. Geometry Builders

Functions to construct lattice adjacency from parameters:

```rust
pub fn build_chain(n_sites: usize, pbc: bool) -> Lattice;

pub fn build_square(lx: usize, ly: usize, pbc: bool) -> Lattice;

pub fn build_triangular(lx: usize, ly: usize, pbc: bool) -> Lattice;

pub fn build_honeycomb(lx: usize, ly: usize, pbc: bool) -> Lattice;

pub fn build_custom(adjacency: Vec<Vec<(usize, BondType)>>) -> Lattice;
```

---

## VIII. Module Structure

```
QMC.rs/src/
├── lib.rs
├── lattice/
│   ├── mod.rs          # Lattice, Neighbor, BondType
│   ├── builders.rs     # build_square, build_triangular, etc.
│   └── bond.rs         # BondType enum
├── hilbert/
│   ├── mod.rs          # HilbertSpace trait, LocalState, OpType
│   ├── spin_half.rs    # SpinHalfHS implementation
│   └── hubbard.rs      # HubbardHS implementation (future)
├── sse/
│   ├── mod.rs          # SSEMonteCarlo trait, SSECore
│   ├── engine.rs       # SSEEngine, OperatorSequence, Vertex
│   ├── diagonal.rs     # diagonal_update implementation
│   └── loop.rs         # loop_update (Directed Loop) implementation
│   └── measurements.rs # compute_energy, compute_magnetization
├── models/
│   ├── mod.rs
│   ├── heisenberg.rs   # HeisenbergModel
│   ├── xxz.rs          # XXZModel (future)
│   └── hubbard.rs      # HubbardModel (future)
```

---

## IX. Phase 1 Implementation Scope

| Component | Status |
|-----------|--------|
| `Lattice` + adjacency list | ✓ Designed |
| `BondType` enum | ✓ Designed |
| `HilbertSpace` trait | ✓ Designed |
| `SpinHalfHS` | ✓ Designed |
| `SSEEngine` structure | ✓ Designed |
| `SSEMonteCarlo` trait | ✓ Designed |
| `SSECore` wrapper | ✓ Designed |
| `HeisenbergModel` | ✓ Designed |
| Geometry builders (chain, square) | ✓ Designed |
| `diagonal_update` algorithm | ⏳ Implementation pending |
| `loop_update` (Directed Loop) | ⏳ Implementation pending |
| Measurements (energy, magnetization) | ⏳ Implementation pending |

---

## X. Verification Strategy

1. **Correctness**: Compare Heisenberg chain results with exact Bethe ansatz (1D) and Carlo.jl reference
2. **Performance**: Benchmark sweep rate against Fortran SSE implementations
3. **API usability**: New model implementation should take < 1 day for experienced users

---

*Design established through brainstorming session on 2026-04-09.*