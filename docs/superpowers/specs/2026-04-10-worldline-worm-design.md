# Worldline + Worm Algorithm Design

> Worldline Monte Carlo for Heisenberg model on QMC.rs, with plain local update and worm update as two update modes sharing the same model layer.

## 1. Motivation

QMC.rs philosophy defines a **progressive complexity** principle: users start simple, then swap in advanced updates. The current SSE implementation provides one algorithmic approach; worldline provides another, with worm being a further improvement over plain local updates.

**Goal:** Implement worldline Monte Carlo for Heisenberg (spin-1/2 chain), with two update modes:
1. **Plain local update** — classical Metropolis single-bond flip
2. **Worm update** — topological defect insertion/movement/closure

Both modes share the same `WorldlineEngine`, `HeisenbergModel`, and Carlo.rs framework. Success is measured by matching SSE results within statistical error.

## 2. Architecture

### 2.1 File Structure

```
QMC.rs/src/
├── sse/              # Existing
├── worldline/        # NEW
│   ├── mod.rs        # Module entry, WorldlineCore, UpdateMode enum
│   ├── engine.rs     # WorldlineEngine, WorldlineConfig
│   ├── local.rs      # Plain Metropolis local update
│   ├── worm.rs       # Worm update (insert/move/close)
│   └── measurements.rs  # Energy, magnetization
QMC.rs/tests/
├── worldline_test.rs     # Unit tests
└── worldline_vs_sse.rs   # Physics validation
```

### 2.2 Core Data Structures

**`WorldlineConfig`** — holds the spacetime spin configuration:

```rust
pub struct WorldlineConfig {
    pub n_sites: usize,
    pub n_trotter: usize,
    /// Flattened: spins[site * n_trotter + time]
    pub spins: Vec<LocalState>,
}
```

- 2D spacetime lattice: N sites × M Trotter slices
- Periodic in both space (via lattice) and imaginary time
- Flattened 1D vector for cache-friendly access

**`WorldlineEngine`** — the computation engine:

```rust
pub struct WorldlineEngine {
    pub config: WorldlineConfig,
    pub lattice: Lattice,
    pub beta: f64,
    pub j: f64,
}
```

No operator sequence — pure classical spin configuration on the spacetime lattice.

**`WorldlineCore<MC>`** — MonteCarlo wrapper (like SSECore):

```rust
pub struct WorldlineCore<MC: WorldlineMonteCarlo> {
    pub engine: WorldlineEngine,
    pub mc: MC,
    pub update_mode: WorldlineUpdateMode,
}

pub enum WorldlineUpdateMode {
    Local,
    Worm,
}
```

### 2.3 Trait Layer

```rust
pub trait WorldlineMonteCarlo: LatticeQMC {
    type HilbertSpace: HilbertSpace;
    fn hilbert_space(&self) -> &Self::HilbertSpace;
    fn beta(&self) -> f64;
    fn trotter_slices(&self) -> usize;
    fn coupling(&self, bond_type: BondType) -> f64;
}
```

The `HeisenbergModel` implements both `SSEMonteCarlo` and `WorldlineMonteCarlo`, demonstrating the progressive complexity principle.

### 2.4 MonteCarlo Implementation

```rust
impl<MC: WorldlineMonteCarlo> MonteCarlo for WorldlineCore<MC> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        match self.update_mode {
            WorldlineUpdateMode::Local => self.engine.local_update(&mut ctx.rng),
            WorldlineUpdateMode::Worm => self.engine.worm_update(&mut ctx.rng),
        }
        ctx.advance_sweep();
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let energy = self.engine.compute_energy();
        ctx.measure("Energy", energy);
    }
}
```

## 3. Plain Local Update

### 3.1 Algorithm

For each sweep, attempt N_bonds × N_trotter local moves:

1. Randomly select a plaquette (bond × time slice)
2. Compute energy change ΔE for flipping the two spins
3. Accept with probability `min(1, exp(-β·ΔE/M))`

The Trotter-discretized Heisenberg Hamiltonian:

```
H ≈ -J/(4M) Σ_{<ij>,τ} (S^z_i(τ)S^z_j(τ) + S^x_i(τ)S^x_j(τ) + S^y_i(τ)S^y_j(τ))
    - (M/β) Σ_{i,τ} S^z_i(τ)S^z_i(τ+1)
```

The local flip changes both the spatial interaction term and the temporal kinetic term.

### 3.2 Implementation (`local.rs`)

```rust
impl WorldlineEngine {
    pub fn local_update<R: Rng>(&mut self, rng: &mut R) {
        let n_bonds = self.lattice.n_bonds;
        for _ in 0..(n_bonds * self.config.n_trotter) {
            let bond = rng.random_range(0..n_bonds);
            let time = rng.random_range(0..self.config.n_trotter);
            let delta_e = self.compute_local_energy_change(bond, time);
            let p = (-self.beta * delta_e / self.config.n_trotter as f64).exp();
            if rng.random::<f64>() < p {
                self.flip_spins_on_bond(bond, time);
            }
        }
    }
}
```

## 4. Worm Update

### 4.1 Algorithm

The worm algorithm introduces a topological defect (worm) that moves through the spacetime lattice:

1. **Insert worm** — select random (site, time), flip spin there. This creates a worm head and breaks the periodic constraint at that point.
2. **Move worm** — the worm head wanders through neighboring sites, flipping spins as it goes. At each step, it chooses a neighbor weighted by the Boltzmann factor of the resulting energy change.
3. **Close worm** — when the worm head returns to its origin position, flip that spin back. The configuration is again periodic. If the worm fails to close after a maximum number of steps, erase it (backtrack all flips).

### 4.2 Implementation (`worm.rs`)

```rust
impl WorldlineEngine {
    pub fn worm_update<R: Rng>(&mut self, rng: &mut R) {
        let start_site = rng.random_range(0..self.config.n_sites);
        let start_time = rng.random_range(0..self.config.n_trotter);
        let start = (start_site, start_time);

        // Insert: flip spin at start
        self.flip_spin(start_site, start_time);

        let mut pos = start;
        let max_steps = self.config.n_sites * self.config.n_trotter * 2;

        for step in 0..max_steps {
            if step > 0 && pos == start {
                // Worm closed successfully - flip back
                self.flip_spin(start_site, start_time);
                return;
            }

            // Choose next move
            let neighbors = self.get_neighbors(pos);
            let (next_site, next_time) = self.choose_neighbor(rng, &neighbors, pos);

            // Flip spin at new position and move worm
            self.flip_spin(next_site, next_time);
            pos = (next_site, next_time);
        }

    // Failed to close: backtrack all flips
    // Each flip is recorded in a Vec<(site, time)> during the walk.
    // On failure, iterate the Vec in reverse and flip each back.
    // The Vec is cleared after each update (success or failure).
    }
}
```

### 4.3 Neighbor Selection

At each step, the worm chooses from valid neighbors:

- **Spatial neighbors** — sites connected by lattice bonds (same time slice)
- **Temporal neighbors** — same site, adjacent time slice

The selection weight for each neighbor is the matrix element of the Hamiltonian between the current and proposed configuration. For Heisenberg, the off-diagonal term gives non-zero weight only for spin-flip moves.

## 5. Measurements

### 5.1 HDF5 Output Format

Worldline measurements produce HDF5 output compatible with Carlo.jl and the SSE implementation, enabling unified merge/analysis:

```
observables/
├── Energy/
│   ├── bin_length: [u64]
│   └── samples: [f64 × N_samples]
└── Magnetization/
    ├── bin_length: [u64]
    └── samples: [f64 × N_samples]
```

This format matches Carlo.rs's merge module expectations, so `carlo merge` works for both SSE and worldline runs.

### 5.2 Energy

The worldline energy estimator:

```
E = -∂ln(Z)/∂β ≈ -(1/N) Σ_{<ij>,τ} <S_i(τ)·S_j(τ)> / M
```

Compared to SSE's `E = -<n>/(β·N)`, this measures the same quantity through a different route.

### 5.3 Magnetization

```
M = (1/N) |Σ_i S^z_i|
```

Measured on a single time slice (all slices are equivalent for diagonal observables).

## 6. Testing Strategy

### 6.1 Unit Tests

- `test_worldline_config_init` — configuration created with correct dimensions
- `test_local_update_does_not_crash` — basic sanity
- `test_worm_closes` — worm update always returns to origin or backtracks cleanly
- `test_energy_zero_empty` — empty configuration has zero energy

### 6.2 Physics Validation

- `test_worldline_vs_sse` — Same parameters (N=8, β=4, J=1), worldline energy matches SSE energy within 3σ
- `test_worldline_bethe_ansatz` — Low-T 1D Heisenberg chain energy per site → -0.443147

### 6.3 Worm Correctness

- `test_worm_detailed_balance` — worm and local update give same energy distribution (within error)
- `test_worm_acceptance_rate` — worm closure rate is reasonable (> 50%)

## 7. Integration with Carlo.rs

The `WorldlineCore` implements `MonteCarlo` and `FromParams`, enabling:

```rust
let scheduler = Scheduler::new(backend, config);
let results = scheduler.run_one::<WorldlineCore<HeisenbergModel>>(&params)?;
```

The CLI can run worldline simulations with the same `carlo run` command.

## 8. Implementation Order

1. **WorldlineConfig + WorldlineEngine** — basic data structures
2. **Plain local update** — Metropolis single-bond flip
3. **WorldlineCore + MonteCarlo impl** — framework integration
4. **Measurements** — energy and magnetization
5. **Worm update** — insert/move/close
6. **Unit tests** — local update, worm closure
7. **Physics validation** — worldline vs SSE comparison
