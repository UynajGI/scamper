# CMC.rs Quality Guidelines

## Required patterns

### Model is stateless

```
// CORRECT: Model is pure parameters + formulas
impl Model for IsingModel {
    fn local_energy(&self, spins: &[f64], lattice: &Lattice, site: usize, proposed: &[f64]) -> f64 {
        let s = proposed[0];
        let mut e = 0.0;
        for nb in &lattice.sites[site] {
            e += -self.j * s * spins[nb.target];
        }
        e
    }
}

// WRONG: Model holds state or mutates itself
struct BadModel {
    j: f64,
    energy: f64,  // state belongs in System, not Model
}
```

### System fields are pub

The `System` struct exposes all fields as `pub`. Algorithms freely read and write `spins` and `energy`. No getters/setters.

```rust
pub struct System {
    pub lattice: Lattice,
    pub spins: Vec<f64>,   // flattened, spin_dim components per site
    pub energy: f64,       // running total, updated incrementally
}
```

### Algorithm directly mutates system.energy

```rust
// CORRECT: energy tracked in system, updated during sweep
fn sweep(&mut self, system: &mut System, model: &M, rng: &mut impl Rng) {
    for &site in &order {
        let old_e = model.local_energy(..., &old_spin);
        let new_e = model.local_energy(..., &proposed);
        let delta = new_e - old_e;
        if accept {
            system.spin_at_mut(site, sd).copy_from_slice(&proposed);
            system.energy += delta;  // ← direct mutation
        }
    }
}

// WRONG: returning delta from sweep, or tracking energy in algorithm struct
```

### ClassicalMC composes, doesn't own physics

`ClassicalMC<M, A>` is a thin wrapper. It does NOT:
- Implement energy formulas (those are in Model)
- Implement update logic (that's in Algorithm)
- Implement lattice building (that's in lattice.rs)

It only wires Carlo.rs traits (MonteCarlo, FromParams) to CMC types.

### Each model has FromModelParams

Every model type gets a `FromModelParams` impl in `classical_mc.rs`, even if it only reads `J` and `beta`. No default/blanket impl.

```rust
impl FromModelParams for XYModel {
    fn from_model_params(params: &Params) -> Result<Self, CarloError> {
        let j = params.get::<f64>("J").unwrap_or(1.0);
        let beta = params.get::<f64>("beta").unwrap_or(1.0);
        Ok(Self::new(j, beta))
    }
}
```

## Forbidden patterns

### Don't put energy state in Model

Model is stateless physics. Energy, spins, lattice topology all belong to System.

### Don't re-implement energy tracking per algorithm

Every algorithm updates `system.energy` the same way: compute `local_energy` before and after, `system.energy += delta`. Don't invent alternative energy bookkeeping.

### Don't put algorithm logic in ClassicalMC

`ClassicalMC` is a pure delegation wrapper. If you find yourself writing sweep logic in `classical_mc.rs`, it belongs in `algorithm.rs`.

### Lattice routing via `lattice_type` param

Use `build_lattice_from_params(params, pbc)` to route lattice construction. Add new lattice types by extending the match in this function.

```rust
// CORRECT: single function routes to all builders
fn build_lattice_from_params(params: &Params, pbc: bool) -> Result<CsrLattice, CarloError> {
    match params.get::<String>("lattice_type").unwrap_or_default().as_str() {
        "triangular" => Ok(build_triangular(lx, ly)),
        "honeycomb" => Ok(build_honeycomb(lx, ly)),
        "kagome" => Ok(build_kagome(lx, ly)),
        _ => { /* hypercubic family */ }
    }
}

// WRONG: duplicating lattice-building logic in each FromParams impl
```

### Don't compute total_energy from scratch per sweep

Use `model.compute_total_energy()` ONLY for initialization. During sweeps, use incremental updates (`system.energy += delta`). Computing `total_energy()` from scratch every sweep is O(N) unnecessary work.

## Testing requirements

### Per model
- `local_energy` for aligned and anti-aligned neighbor pairs
- `total_energy` on small ring (known exact value)
- `magnetization` for ordered and random configurations
- `fk_bond_probability` if model overrides default
- `normalize_spin` for vector models (XY, Heisenberg)

### Per algorithm
- Single-site lattice (no bonds, energy stays 0)
- Convergence to ground state at high beta on small ring

### Integration
- End-to-end `Scheduler.run_one::<ClassicalMC<M, A>>()` at moderate beta
- Onsager 2D Ising validation (energy per site at Tc, high-T magnetization vanish, low-T magnetization appear)

### Derived observables: measure raw moments, post-process

Susceptibility, specific heat, and Binder cumulant require statistics of fluctuations (⟨M²⟩, ⟨E²⟩, ⟨M⁴⟩). Record the raw moments in `measure()` with `ctx.measure()`, then compute derived quantities from `Results` in `postprocess.rs`.

```rust
// CORRECT: record raw moments during measurement
fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
    ctx.measure("E2", e * e);
    ctx.measure("M2", m * m);
    ctx.measure("M4", m * m * m * m);
}

// postprocess.rs — pure functions, no CMC internals
pub fn susceptibility(results: &Results, beta: f64, n: usize) -> Option<f64> {
    let m = results.get("Magnetization")?;
    let m2 = results.get("M2")?;
    Some(beta * n as f64 * (m2.mean - m.mean * m.mean))
}

// WRONG: computing derived quantities inside measure() or sweep()
// WRONG: Observable<H> trait trying to access time-series statistics
```

### JSON checkpoint when HDF5 unavailable

When the `hdf5` feature doesn't compile, use JSON snapshot as a practical alternative. Serialize full state (spins, energy, beta, lattice topology) to serde_json::Value.

```rust
// CORRECT: JSON round-trip captures config + state
impl ClassicalMC {
    pub fn save_snapshot(&self) -> Json { ... }
    pub fn load_snapshot(&mut self, snapshot: &Json) -> Result<(), CarloError> { ... }
}

// WRONG: implementing MonteCarloCheckpoint with HDF5 when the feature is broken
```

### MultiSpinIsing MonteCarlo integration

Non-ClassicalMC types can implement MonteCarlo directly by owning system/model/lattice fields. The sweep delegates to the internal sweep method; measure follows the same moments pattern as ClassicalMC.

```rust
// CORRECT: own system, model; delegate in MonteCarlo::sweep
impl MonteCarlo for MultiSpinIsing {
    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.sweep(&mut ctx.rng);
    }
}

// WRONG: trying to fit bit-packed code into ClassicalMC's generic structure
```

### Don't store duplicate lattice

`System` already owns the lattice. Don't store a second copy.

```rust
// CORRECT
pub struct MultiSpinIsing {
    pub system: System,  // use self.system.lattice
    pub model: IsingModel,
}

// WRONG
pub struct MultiSpinIsing {
    pub system: System,
    pub lattice: CsrLattice,  // duplicate! use self.system.lattice
    pub model: IsingModel,
}
```

### Don't use bit-plane encoding for multi-spin neighbor counting

The original `pack_bit_planes` transposition conflates neighbor index with replica index. Use a direct per-replica anti-aligned counter instead.

```rust
// CORRECT: count anti-aligned neighbors per replica directly
let mut anti_counts = [0u8; 64];
for &nb in self.system.lattice.neighbors(site) {
    let xor = site_word ^ self.packed_spins[nb];
    for (r, count) in anti_counts.iter_mut().enumerate() {
        *count += ((xor >> r) & 1) as u8;
    }
}

// WRONG: pack_bit_planes — shifts bits into wrong positions,
// conflates replica r with neighbor r, overflows for anti > z
```

## Pre-commit gates

```bash
cargo fmt --check    # rustfmt
cargo clippy -- -D warnings   # no warnings allowed
cargo test --all-targets      # all tests pass
```
