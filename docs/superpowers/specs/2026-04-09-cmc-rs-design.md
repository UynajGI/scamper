# CMC.rs Design Specification

> Design for CMC.rs implementing Classical Monte Carlo algorithms with Metropolis, Wolff, and Swendsen-Wang update strategies.
>
> Date: 2026-04-09

---

## I. Philosophy Alignment

This design implements Phase 1 of the [CMC.rs Philosophy](../metaphysics/2026-04-09-cmc-rs-philosophy.md):

- **Position**: Algorithm toolbox (component library) between Carlo.rs framework and model implementations
- **Granularity**: Mid-level — users implement `ModelMC` trait, engines handle update algorithms
- **Relationship**: Equal sibling of QMC.rs — symmetric architecture, independent implementations
- **Design principle**: Progressive complexity (3 layers)

---

## II. Module Structure

```
CMC.rs/
├── Cargo.toml
└── src/
    ├── lib.rs                  # Public exports
    ├── lattice/
    │   ├── mod.rs              # Lattice, Neighbor, LatticeMC trait
    │   ├── bond.rs             # BondType enum
    │   └── builders.rs         # build_chain, build_square
    ├── models/
    │   ├── mod.rs
    │   └── ising.rs            # IsingModel (impl ModelMC + FromParams)
    └── algorithms/
        ├── mod.rs
        ├── metropolis.rs       # MetropolisCore<MC>
        ├── wolff.rs            # WolffCore<MC>
        └── swendsen_wang.rs    # SWCore<MC>
```

**Dependency direction**: `CMC.rs → Carlo.rs` only. No dependency on QMC.rs. All lattice types independently defined.

---

## III. Core Data Structures

### III.1 Topology Layer

Lattice is represented as an **adjacency list** — identical design to QMC.rs but independently defined.

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondType {
    ChainX,
    SquareX, SquareY,
    TriX, TriY, TriZ,
    HoneyX, HoneyY, HoneyZ,
    Custom(u8),
}

#[derive(Clone, Debug)]
pub struct Neighbor {
    pub target: usize,
    pub bond_type: BondType,
}

#[derive(Clone, Debug)]
pub struct Lattice {
    pub sites: Vec<Vec<Neighbor>>,
    pub n_sites: usize,
    pub n_bonds: usize,
}
```

**Design decisions**:
- Adjacency list stores topology only; coupling constants managed by model
- `BondType` tags allow direction-dependent Hamiltonian parameters
- Geometry builders construct adjacency from parameters

### III.2 Physics Layer

```rust
pub trait LatticeMC {
    fn lattice(&self) -> &Lattice;
    fn n_sites(&self) -> usize { self.lattice().n_sites }
}

pub trait ModelMC: LatticeMC {
    fn spin_dim(&self) -> usize;
    fn coupling(&self) -> f64;
    fn beta(&self) -> f64;
    fn local_energy_change(&self, site: usize, old: f64, new: f64) -> f64;
    fn total_energy(&self) -> f64;
    fn spins(&self) -> &[f64];
    fn spins_mut(&mut self) -> &mut [f64];
}
```

**Convention**:
- `local_energy_change` returns ΔE for flipping site from `old` to `new` state
- `spins` uses `f64` for uniform interface across Ising (±1), XY (angle), Heisenberg (vector)
- Coupling constants come from model directly (no HilbertSpace abstraction needed for classical MC)

### III.3 Algorithm Engines

```rust
pub struct MetropolisCore<MC: ModelMC> {
    model: MC,
}

pub struct WolffCore<MC: ModelMC> {
    model: MC,
}

pub struct SWCore<MC: ModelMC> {
    model: MC,
}
```

Each `*Core` implements `MonteCarlo` trait for Carlo.rs scheduler integration.

---

## IV. Trait Hierarchy

```rust
// === Carlo.rs Base (existing) ===
pub trait MonteCarlo: Sized {
    type Rng: Rng + SeedableRng + Send;
    fn sweep(&mut self, ctx: &mut Context<Self::Rng>);
    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {}
}

// === CMC.rs Domain Layer ===
pub trait LatticeMC: MonteCarlo {
    fn lattice(&self) -> &Lattice;
    fn n_sites(&self) -> usize { self.lattice().n_sites }
}

// === CMC.rs Method Layer ===
pub trait ModelMC: LatticeMC {
    fn spin_dim(&self) -> usize;
    fn coupling(&self) -> f64;
    fn beta(&self) -> f64;
    fn local_energy_change(&self, site: usize, old: f64, new: f64) -> f64;
    fn total_energy(&self) -> f64;
    fn spins(&self) -> &[f64];
    fn spins_mut(&mut self) -> &mut [f64];
}
```

**Note**: `LatticeMC: MonteCarlo` bound ensures compatibility with Carlo.rs scheduler. The actual `MonteCarlo` implementation is provided by `*Core` wrappers, not by user models directly.

---

## V. *Core Wrappers: Default Implementations

### V.1 MetropolisCore

```rust
impl<MC: ModelMC> MonteCarlo for MetropolisCore<MC> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        let n = self.model.n_sites();
        for _ in 0..n {
            let site = ctx.rng.random_range(0..n);
            let old = self.model.spins()[site];
            let new = -old; // Ising flip
            let de = self.model.local_energy_change(site, old, new);
            if de < 0.0 || ctx.rng.random::<f64>() < (-self.model.beta() * de).exp() {
                self.model.spins_mut()[site] = new;
            }
        }
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let energy = self.model.total_energy();
        ctx.measure("Energy", energy);
    }
}
```

### V.2 WolffCore

```rust
impl<MC: ModelMC> MonteCarlo for WolffCore<MC> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // Wolff cluster algorithm:
        // 1. Random seed site
        // 2. p_add = 1 - exp(-2*J*beta) grow cluster via neighbors
        // 3. Flip entire cluster
        // 4. Repeat until cluster size saturates
    }
}
```

### V.3 SWCore

```rust
impl<MC: ModelMC> MonteCarlo for SWCore<MC> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // Swendsen-Wang algorithm:
        // 1. Build Fortuin-Kasteleyn bonds with p_add = 1 - exp(-2*J*beta)
        // 2. Identify clusters (Union-Find)
        // 3. Assign random spin to each cluster
    }
}
```

**Each engine requires**: `MC: ModelMC + FromParams` for Carlo.rs scheduler compatibility.

---

## VI. FromParams Integration

```rust
impl<MC: ModelMC + FromParams<Rng = Xoshiro256PlusPlus>> FromParams for MetropolisCore<MC> {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let mc = MC::from_params(params, rng)?;
        Ok(MetropolisCore { model: mc })
    }
}
```

Same pattern for `WolffCore` and `SWCore`. This allows direct use with `Scheduler::run_one::<MetropolisCore<IsingModel>>(&params)`.

---

## VII. Progressive Complexity

### Level 1: Built-in

```rust
let model = IsingModel::new(lattice, beta, j);
let core = MetropolisCore { model };
// core implements MonteCarlo → use with Carlo.rs scheduler
```

### Level 2: Custom Model

```rust
struct MyModel { ... }

impl LatticeMC for MyModel {
    fn lattice(&self) -> &Lattice { &self.lattice }
}

impl ModelMC for MyModel {
    fn spin_dim(&self) -> usize { 1 }
    fn coupling(&self) -> f64 { self.j }
    fn beta(&self) -> f64 { self.beta }
    // ... energy and spin access
}

let core = WolffCore { model: MyModel::new(...) };
scheduler.run_one::<WolffCore<MyModel>>(&params);
```

### Level 3: Custom Engine

Advanced users implement their own `*Core` or new algorithm entirely.

---

## VIII. Geometry Builders

```rust
pub fn build_chain(n_sites: usize, pbc: bool) -> Lattice;

pub fn build_square(lx: usize, ly: usize, pbc: bool) -> Lattice;

pub fn build_triangular(lx: usize, ly: usize, pbc: bool) -> Lattice;

pub fn build_honeycomb(lx: usize, ly: usize, pbc: bool) -> Lattice;
```

Same signatures as QMC.rs lattice builders, independently implemented.

---

## IX. Phase 1 Implementation Scope

| Component | Status |
|-----------|--------|
| `Lattice` + adjacency list | ⏳ Design ready |
| `BondType` enum | ⏳ Design ready |
| `LatticeMC` trait | ⏳ Design ready |
| `ModelMC` trait | ⏳ Design ready |
| `IsingModel` | ⏳ Design ready |
| `MetropolisCore` algorithm | ⏳ Design ready |
| `WolffCore` algorithm | ⏳ Design ready |
| `SWCore` algorithm | ⏳ Design ready |
| Geometry builders (chain, square) | ⏳ Design ready |
| Validation vs Onsager exact solution | ⏳ Design ready |
| Validation vs Carlo.jl reference | ⏳ Design ready |

---

## X. Verification Strategy

1. **Onsager Exact Solution**: 2D Ising critical temperature `Tc = 2/sinh⁻¹(1) ≈ 2.269`, energy at Tc matches exact value within 3σ
2. **Carlo.jl Comparison**: Same seeds → same results within floating-point tolerance
3. **Performance Benchmark**: Sweep rate compared to optimized C++/Fortran implementations

---

*Design established through brainstorming session on 2026-04-09.*
