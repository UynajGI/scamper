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

## Pre-commit gates

```bash
cargo fmt --check    # rustfmt
cargo clippy -- -D warnings   # no warnings allowed
cargo test --all-targets      # all tests pass
```
