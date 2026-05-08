# CMC.rs Phase C+A+B Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add magnetization/snapshot measurements, fix SW/Wolff for Potts, and implement Heisenberg model with optional OPSS.

**Architecture:** Extend the existing `ModelMC` trait with 3 new methods (`magnetization`, `snapshot`, `random_cluster_spin`/`opposite_spin`), update all algorithm `measure()` methods to record fluctuations, and add `HeisenbergModel` as a new model.

**Tech Stack:** Rust, Carlo.rs framework (MonteCarlo/FromParams/Context), rand, xoshiro256plus

---

### Task 1: Phase C — Add `magnetization()` and `snapshot()` to ModelMC

**Files:**
- Modify: `CMC.rs/src/models/mod.rs` — Add trait methods
- Modify: `CMC.rs/src/models/ising.rs` — Implement for Ising
- Modify: `CMC.rs/src/models/ising_2d.rs` — Implement for Ising2D
- Modify: `CMC.rs/src/models/potts.rs` — Implement for Potts
- Modify: `CMC.rs/src/models/xy.rs` — Implement for XY

- [ ] **Step 1: Add `magnetization()` and `snapshot()` to ModelMC trait**

In `CMC.rs/src/models/mod.rs`, add these methods to the `ModelMC` trait (after `spins_mut`):

```rust
/// Magnetization of the current configuration.
/// Model-specific: Ising=|Σs_i|/N, Potts=(q·max(n_k)-N)/(N·(q-1)), XY=|Σ(cosθ,sinθ)|/N
fn magnetization(&self) -> f64;

/// Raw spin configuration snapshot as Vec<f64>.
fn snapshot(&self) -> Vec<f64> {
    self.spins().to_vec()
}
```

- [ ] **Step 2: Implement `magnetization()` for IsingModel**

In `CMC.rs/src/models/ising.rs`, add to `impl ModelMC for IsingModel`:

```rust
fn magnetization(&self) -> f64 {
    let sum: f64 = self.spins.iter().sum();
    sum.abs() / self.spins.len() as f64
}
```

- [ ] **Step 3: Implement `magnetization()` for IsingModel2D**

In `CMC.rs/src/models/ising_2d.rs`, same formula as IsingModel:

```rust
fn magnetization(&self) -> f64 {
    let sum: f64 = self.spins.iter().sum();
    sum.abs() / self.spins.len() as f64
}
```

- [ ] **Step 4: Implement `magnetization()` for PottsModel**

In `CMC.rs/src/models/potts.rs`, add to `impl ModelMC for PottsModel`:

```rust
fn magnetization(&self) -> f64 {
    // Count spins in each state
    let n = self.spins.len();
    let mut counts = vec![0usize; self.q];
    for &s in self.spins.iter() {
        let idx = s as usize;
        if idx < self.q {
            counts[idx] += 1;
        }
    }
    let max_count = counts.into_iter().max().unwrap_or(0);
    (self.q as f64 * max_count as f64 - n as f64) / (n as f64 * (self.q as f64 - 1.0))
}
```

- [ ] **Step 5: Implement `magnetization()` for XYModel**

In `CMC.rs/src/models/xy.rs`, add to `impl ModelMC for XYModel`:

```rust
fn magnetization(&self) -> f64 {
    let n = self.spins.len() as f64;
    let (cx, sy) = self.spins.iter().fold((0.0, 0.0), |(cx, sy), &theta| {
        (cx + theta.cos(), sy + theta.sin())
    });
    (cx * cx + sy * sy).sqrt() / n
}
```

- [ ] **Step 6: Add unit tests for magnetization**

In `CMC.rs/src/models/ising.rs` test module, add:

```rust
#[test]
fn test_ising_magnetization_all_up() {
    let lattice = build_chain(4, true);
    let model = IsingModel::new(lattice, 1.0, 1.0);
    assert!((model.magnetization() - 1.0).abs() < 1e-10);
}
```

In `CMC.rs/src/models/potts.rs` test module, add:

```rust
#[test]
fn test_potts_magnetization_all_same() {
    let lattice = build_chain(4, true);
    let model = PottsModel::new(lattice, 1.0, 1.0, 3);
    // All spins 0 → max_count = 4 → (3*4 - 4) / (4*2) = 8/8 = 1.0
    assert!((model.magnetization() - 1.0).abs() < 1e-10);
}
```

In `CMC.rs/src/models/xy.rs` test module, add:

```rust
#[test]
fn test_xy_magnetization_aligned() {
    let lattice = build_chain(4, true);
    let mut model = XYModel::new(lattice, 1.0, 1.0);
    // All angles 0 → perfect alignment
    assert!((model.magnetization() - 1.0).abs() < 1e-10);
}
```

- [ ] **Step 7: Run tests and commit**

Run: `cd CMC.rs && cargo test --lib -q`
Expected: All tests pass (existing + new magnetization tests)

```bash
cd CMC.rs && git add src/models/mod.rs src/models/ising.rs src/models/ising_2d.rs src/models/potts.rs src/models/xy.rs
git commit -m "feat(cmc): add magnetization() and snapshot() to ModelMC trait

Implements model-specific magnetization:
- Ising: |Σs_i|/N
- Potts: (q·max(n_k)-N)/(N·(q-1))
- XY: |Σ(cosθ,sinθ)|/N
Default snapshot() returns spins().to_vec()"
```

---

### Task 2: Phase C — Extend Algorithm `measure()` Methods

**Files:**
- Modify: `CMC.rs/src/algorithms/metropolis.rs` — Extend measure()
- Modify: `CMC.rs/src/algorithms/wolff.rs` — Extend measure()
- Modify: `CMC.rs/src/algorithms/swendsen_wang.rs` — Extend measure()

- [ ] **Step 1: Extend MetropolisCore measure()**

In `CMC.rs/src/algorithms/metropolis.rs`, replace the `measure` method:

```rust
fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
    let energy = self.model.total_energy();
    let magnetization = self.model.magnetization();
    ctx.measure("Energy", energy);
    ctx.measure("Energy²", energy * energy);
    ctx.measure("Magnetization", magnetization);
    ctx.measure("Magnetization²", magnetization * magnetization);
}
```

- [ ] **Step 2: Extend WolffCore measure()**

In `CMC.rs/src/algorithms/wolff.rs`, replace the `measure` method:

```rust
fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
    let energy = self.model.total_energy();
    let magnetization = self.model.magnetization();
    ctx.measure("Energy", energy);
    ctx.measure("Energy²", energy * energy);
    ctx.measure("Magnetization", magnetization);
    ctx.measure("Magnetization²", magnetization * magnetization);
}
```

- [ ] **Step 3: Extend SWCore measure()**

In `CMC.rs/src/algorithms/swendsen_wang.rs`, replace the `measure` method:

```rust
fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
    let energy = self.model.total_energy();
    let magnetization = self.model.magnetization();
    ctx.measure("Energy", energy);
    ctx.measure("Energy²", energy * energy);
    ctx.measure("Magnetization", magnetization);
    ctx.measure("Magnetization²", magnetization * magnetization);
}
```

- [ ] **Step 4: Run tests and commit**

Run: `cd CMC.rs && cargo test -q`
Expected: All tests pass (including integration tests which check Energy)

```bash
cd CMC.rs && git add src/algorithms/metropolis.rs src/algorithms/wolff.rs src/algorithms/swendsen_wang.rs
git commit -m "feat(cmc): extend measure() to record Energy², Magnetization, Magnetization²

All three algorithms (Metropolis, Wolff, SW) now record:
- Energy (existing)
- Energy² (for specific heat via fluctuations)
- Magnetization (model-specific order parameter)
- Magnetization² (for susceptibility via fluctuations)"
```

---

### Task 3: Phase A — Add `random_cluster_spin()` and `opposite_spin()` to ModelMC

**Files:**
- Modify: `CMC.rs/src/models/mod.rs` — Add trait methods
- Modify: `CMC.rs/src/models/ising.rs` — Implement for Ising
- Modify: `CMC.rs/src/models/ising_2d.rs` — Implement for Ising2D
- Modify: `CMC.rs/src/models/potts.rs` — Implement for Potts
- Modify: `CMC.rs/src/models/xy.rs` — Implement for XY (returns 0.0, not applicable)
- Modify: `CMC.rs/src/algorithms/swendsen_wang.rs` — Use trait method instead of hardcoded
- Modify: `CMC.rs/src/algorithms/wolff.rs` — Use trait method instead of hardcoded

- [ ] **Step 1: Add `random_cluster_spin()` and `opposite_spin()` to ModelMC trait**

In `CMC.rs/src/models/mod.rs`, add to the `ModelMC` trait:

```rust
/// Random spin value for cluster assignment (SW algorithm).
/// Only meaningful for discrete spin models (Ising, Potts).
fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64;

/// Opposite of a given spin (Wolff cluster flip).
/// Only meaningful for discrete spin models with reflection symmetry.
fn opposite_spin(&self, spin: f64, rng: &mut impl Rng) -> f64;
```

- [ ] **Step 2: Implement for IsingModel**

In `CMC.rs/src/models/ising.rs`:

```rust
fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64 {
    if rng.random::<f64>() < 0.5 { 1.0 } else { -1.0 }
}

fn opposite_spin(&self, spin: f64, _rng: &mut impl Rng) -> f64 {
    -spin
}
```

- [ ] **Step 3: Implement for IsingModel2D**

Same as IsingModel:

```rust
fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64 {
    if rng.random::<f64>() < 0.5 { 1.0 } else { -1.0 }
}

fn opposite_spin(&self, spin: f64, _rng: &mut impl Rng) -> f64 {
    -spin
}
```

- [ ] **Step 4: Implement for PottsModel**

In `CMC.rs/src/models/potts.rs`:

```rust
fn random_cluster_spin(&self, rng: &mut impl Rng) -> f64 {
    rng.random_range(0..self.q) as f64
}

fn opposite_spin(&self, spin: f64, rng: &mut impl Rng) -> f64 {
    let current = spin as usize;
    if self.q <= 1 {
        return spin;
    }
    let pick = rng.random_range(0..self.q - 1);
    if pick < current { pick as f64 } else { (pick + 1) as f64 }
}
```

- [ ] **Step 5: Implement for XYModel (not applicable — return 0.0)**

In `CMC.rs/src/models/xy.rs`:

```rust
fn random_cluster_spin(&self, _rng: &mut impl Rng) -> f64 {
    0.0 // Not applicable for continuous spin models
}

fn opposite_spin(&self, spin: f64, _rng: &mut impl Rng) -> f64 {
    spin + std::f64::consts::PI // Reflection on circle
}
```

- [ ] **Step 6: Fix SWCore to use `random_cluster_spin()`**

In `CMC.rs/src/algorithms/swendsen_wang.rs`, replace lines 87-93:

Before:
```rust
let cluster_spins: HashMap<usize, f64> = cluster_roots
    .into_iter()
    .map(|root| {
        let spin = if ctx.rng.random::<f64>() < 0.5 { 1.0 } else { -1.0 };
        (root, spin)
    })
    .collect();
```

After:
```rust
let cluster_spins: HashMap<usize, f64> = cluster_roots
    .into_iter()
    .map(|root| {
        let spin = self.model.random_cluster_spin(&mut ctx.rng);
        (root, spin)
    })
    .collect();
```

- [ ] **Step 7: Fix WolffCore to use `opposite_spin()`**

In `CMC.rs/src/algorithms/wolff.rs`, replace line 64:

Before:
```rust
let new_spin = -seed_spin;
```

After:
```rust
let new_spin = self.model.opposite_spin(seed_spin, &mut ctx.rng);
```

- [ ] **Step 8: Add unit tests for cluster spin methods**

In `CMC.rs/src/models/potts.rs` test module, add:

```rust
#[test]
fn test_potts_random_cluster_spin() {
    let lattice = build_chain(4, true);
    let model = PottsModel::new(lattice, 1.0, 1.0, 3);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);
    for _ in 0..20 {
        let spin = model.random_cluster_spin(&mut rng);
        assert!(spin >= 0.0 && spin < 3.0);
        assert!(spin == spin.floor()); // Must be integer
    }
}

#[test]
fn test_potts_opposite_spin() {
    let lattice = build_chain(4, true);
    let model = PottsModel::new(lattice, 1.0, 1.0, 4);
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(123);
    for _ in 0..20 {
        let opp = model.opposite_spin(0.0, &mut rng);
        assert!(opp != 0.0, "opposite of 0 must not be 0");
        assert!(opp >= 1.0 && opp <= 3.0);
    }
}
```

- [ ] **Step 9: Run tests and commit**

Run: `cd CMC.rs && cargo test -q`
Expected: All tests pass, including new Potts cluster spin tests

```bash
cd CMC.rs && git add src/models/mod.rs src/models/ising.rs src/models/ising_2d.rs src/models/potts.rs src/models/xy.rs src/algorithms/swendsen_wang.rs src/algorithms/wolff.rs
git commit -m "feat(cmc): fix SW/Wolff for Potts with trait-based cluster spin methods

Add random_cluster_spin() and opposite_spin() to ModelMC trait.
- Ising: ±1.0 random / negate
- Potts: random from 0..q / random from states ≠ current
- XY: PI shift (reflection on circle)

SW now uses random_cluster_spin() instead of hardcoded ±1.0
Wolff now uses opposite_spin() instead of -seed_spin"
```

---

### Task 4: Phase B — Implement HeisenbergModel

**Files:**
- Create: `CMC.rs/src/models/heisenberg.rs` — HeisenbergModel with naive + OPSS
- Modify: `CMC.rs/src/models/mod.rs` — Add module + export
- Modify: `CMC.rs/src/lib.rs` — Re-export HeisenbergModel
- Create: `CMC.rs/tests/heisenberg_test.rs` — Unit tests

- [ ] **Step 1: Write HeisenbergModel struct and basic impl**

Create `CMC.rs/src/models/heisenberg.rs`:

```rust
//! Classical Heisenberg model: H = -J Σ S⃗_i · S⃗_j
//!
//! S⃗_i = (sinθ cosφ, sinθ sinφ, cosθ) — unit vector on S².
//! Spins stored as 3 consecutive f64 values per site.

use crate::lattice::{LatticeMC, build_chain};
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use crate::models::ModelMC;
use rand::Rng;
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::f64::consts::PI;

/// Classical Heisenberg model with naive Metropolis proposal.
pub struct HeisenbergModel {
    lattice: crate::lattice::Lattice,
    beta: f64,
    j: f64,
    spins: Vec<f64>,          // 3 * n_sites: [x0, y0, z0, x1, y1, z1, ...]
    proposal_width: f64,      // Angular perturbation range [-δ, δ]
}

impl HeisenbergModel {
    pub fn new(lattice: crate::lattice::Lattice, beta: f64, j: f64, proposal_width: f64) -> Self {
        let n_sites = lattice.n_sites;
        // Initialize all spins pointing in +z direction
        let mut spins = vec![0.0; 3 * n_sites];
        for i in 0..n_sites {
            spins[3 * i + 2] = 1.0; // z = 1.0
        }
        HeisenbergModel {
            lattice,
            beta,
            j,
            spins,
            proposal_width,
        }
    }

    pub fn proposal_width(&self) -> f64 {
        self.proposal_width
    }
}

impl LatticeMC for HeisenbergModel {
    fn lattice(&self) -> &crate::lattice::Lattice {
        &self.lattice
    }
}
```

- [ ] **Step 2: Implement ModelMC for HeisenbergModel**

Continue in `CMC.rs/src/models/heisenberg.rs`:

```rust
impl ModelMC for HeisenbergModel {
    fn spin_dim(&self) -> usize {
        3
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn propose_flip(&self, site: usize, rng: &mut impl Rng) -> (f64, f64, f64) {
        // Propose small rotation: random axis + random angle
        let sx = self.spins[3 * site];
        let sy = self.spins[3 * site + 1];
        let sz = self.spins[3 * site + 2];

        // Random rotation axis (uniform on S²)
        let cos_theta = rng.random_range(-1.0..1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let phi = rng.random_range(0.0..2.0 * PI);
        let ax = sin_theta * phi.cos();
        let ay = sin_theta * phi.sin();
        let az = cos_theta;

        // Random rotation angle in [-proposal_width, proposal_width]
        let angle = rng.random_range(-self.proposal_width..self.proposal_width);
        let c = angle.cos();
        let s = angle.sin();
        let omc = 1.0 - c;

        // Rodrigues' rotation formula
        let dot = ax * sx + ay * sy + az * sz;
        let new_x = c * sx + s * (ay * sz - az * sy) + omc * dot * ax;
        let new_y = c * sy + s * (az * sx - ax * sz) + omc * dot * ay;
        let new_z = c * sz + s * (ax * sy - ay * sx) + omc * dot * az;

        // Normalize (numerical stability)
        let norm = (new_x * new_x + new_y * new_y + new_z * new_z).sqrt();
        (sx, sy, sz, new_x / norm, new_y / norm, new_z / norm)
    }

    fn local_energy_change(&self, site: usize, _old: f64, _new: f64) -> f64 {
        // For Heisenberg, we need the full spin vector, not scalar values.
        // This method is not used by the naive Metropolis — it computes energy
        // change directly from dot products. See propose_and_accept pattern below.
        unimplemented!("Use total_energy() for Heisenberg Metropolis acceptance")
    }

    fn total_energy(&self) -> f64 {
        let mut energy = 0.0;
        let n = self.lattice.n_sites;
        for (site_idx, neighbors) in self.lattice.sites.iter().enumerate() {
            for neighbor in neighbors {
                let sx1 = self.spins[3 * site_idx];
                let sy1 = self.spins[3 * site_idx + 1];
                let sz1 = self.spins[3 * site_idx + 2];
                let sx2 = self.spins[3 * neighbor.target];
                let sy2 = self.spins[3 * neighbor.target + 1];
                let sz2 = self.spins[3 * neighbor.target + 2];
                // S⃗_i · S⃗_j
                energy -= self.j * (sx1 * sx2 + sy1 * sy2 + sz1 * sz2);
            }
        }
        energy / 2.0 // Bidirectional bonds
    }

    fn spins(&self) -> &[f64] {
        &self.spins
    }

    fn spins_mut(&mut self) -> &mut [f64] {
        &mut self.spins
    }

    fn magnetization(&self) -> f64 {
        let n = self.lattice.n_sites as f64;
        let mut mx = 0.0;
        let mut my = 0.0;
        let mut mz = 0.0;
        for i in 0..self.lattice.n_sites {
            mx += self.spins[3 * i];
            my += self.spins[3 * i + 1];
            mz += self.spins[3 * i + 2];
        }
        (mx * mx + my * my + mz * mz).sqrt() / n
    }
}
```

Wait — the `propose_flip` signature returns `(f64, f64)` for 1D spins. For Heisenberg with 3D spins, the existing MetropolisCore won't work because it expects `(old, new)` scalar and uses `local_energy_change` which takes scalars.

**We need a different approach:** The MetropolisCore as designed calls `propose_flip(site, rng)` → `(f64, f64)`, then `local_energy_change(site, old, new)`. For Heisenberg, we need to propose 3D vectors and compute energy from dot products.

**Option A:** Change `propose_flip` to return `(Vec<f64>, Vec<f64>)` — breaks all existing models.

**Option B (recommended):** Heisenberg gets its own `HeisenbergMetropolisCore<MC>` algorithm wrapper that handles 3D proposals natively, using a helper trait method.

Let me revise: Add a `propose_flip_vec(&self, site, rng) -> (Vec<f64>, Vec<f64>)` method with a default impl that wraps `propose_flip`, and Heisenberg overrides it. MetropolisCore calls `propose_flip` for 1D models; HeisenbergMetropolisCore calls `propose_flip_vec`.

Actually, the simplest approach: Heisenberg overrides `propose_flip` to return `(energy_old_component, energy_new_component)` where the component is the z-projection (used as a proxy). No, that's hacky.

**Best approach:** Create a `HeisenbergMetropolisCore<MC: ModelMC>` that does its own sweep logic, similar to how `MetropolisCore` works but handles 3D spins. It reuses `total_energy()` for acceptance (compute ΔE = E_new - E_old by temporarily modifying spins, computing energy, then reverting if rejected).

Let me revise the HeisenbergModel: `propose_flip` returns the rotation axis and angle, and we implement a separate `apply_proposed_rotation` that mutates spins. Actually, the cleanest path:

**Add `propose_flip_component` with a default impl that works for 1D, and Heisenberg provides `propose_flip_spin` that returns a new 3D spin. Then MetropolisCore uses the scalar path and HeisenbergMetropolisCore uses the 3D path.**

Let me reconsider. The `ModelMC` trait has `propose_flip -> (f64, f64)` which works perfectly for Ising/Potts (spin is a single f64) and XY (angle is a single f64). For Heisenberg, the spin is a 3D vector. We need to extend the trait.

**Revised approach: Add `propose_flip_spin` to ModelMC with a default impl that uses `propose_flip`. Heisenberg overrides it directly.**

```rust
/// Propose a full spin flip. Returns (old_spin_vec, new_spin_vec).
/// Default: wraps propose_flip for 1D spin models.
/// Override for Heisenberg (3D spins).
fn propose_flip_spin(&self, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
    let (old_scalar, new_scalar) = self.propose_flip(site, rng);
    (vec![old_scalar], vec![new_scalar])
}

/// Energy change for a full spin flip (vector version).
/// Default: wraps local_energy_change for 1D spin models.
fn local_energy_change_spin(&self, site: usize, old: &[f64], new: &[f64]) -> f64 {
    self.local_energy_change(site, old[0], new[0])
}
```

For Heisenberg, these return `Vec<f64>` of length 3.

Let me write the revised plan:

- [ ] **Step 2 (revised): Add `propose_flip_spin` and `local_energy_change_spin` to ModelMC**

In `CMC.rs/src/models/mod.rs`, add:

```rust
/// Propose a full spin flip. Returns (old_spin_vec, new_spin_vec).
/// Default: wraps propose_flip for scalar spin models.
fn propose_flip_spin(&self, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
    let (old, new) = self.propose_flip(site, rng);
    (vec![old], vec![new])
}

/// Energy change for a vector spin flip.
/// Default: wraps local_energy_change for scalar spin models.
fn local_energy_change_spin(&self, site: usize, old: &[f64], new: &[f64]) -> f64 {
    self.local_energy_change(site, old[0], new[0])
}
```

For Ising, Potts, XY: no changes needed (use defaults).

For Heisenberg: override both methods (see Step 3).

- [ ] **Step 3: Implement HeisenbergModel with 3D spin methods**

Continue in `CMC.rs/src/models/heisenberg.rs`. Add `impl ModelMC for HeisenbergModel` with the full implementation:

```rust
impl ModelMC for HeisenbergModel {
    fn spin_dim(&self) -> usize {
        3
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn propose_flip(&self, site: usize, rng: &mut impl Rng) -> (f64, f64) {
        // Fallback: return z-component only (not used by HeisenbergMetropolisCore)
        (self.spins[3 * site + 2], 1.0)
    }

    fn local_energy_change(&self, _site: usize, _old: f64, _new: f64) -> f64 {
        // Not used by HeisenbergMetropolisCore
        unimplemented!("Use local_energy_change_spin for Heisenberg")
    }

    fn propose_flip_spin(&self, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
        let sx = self.spins[3 * site];
        let sy = self.spins[3 * site + 1];
        let sz = self.spins[3 * site + 2];
        let old = vec![sx, sy, sz];

        // Random rotation axis (uniform on S²)
        let cos_theta = rng.random_range(-1.0..1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let phi = rng.random_range(0.0..2.0 * PI);
        let ax = sin_theta * phi.cos();
        let ay = sin_theta * phi.sin();
        let az = cos_theta;

        // Random rotation angle in [-proposal_width, proposal_width]
        let angle = rng.random_range(-self.proposal_width..self.proposal_width);
        let c = angle.cos();
        let s = angle.sin();
        let omc = 1.0 - c;

        // Rodrigues' rotation formula
        let dot = ax * sx + ay * sy + az * sz;
        let nx = c * sx + s * (ay * sz - az * sy) + omc * dot * ax;
        let ny = c * sy + s * (az * sx - ax * sz) + omc * dot * ay;
        let nz = c * sz + s * (ax * sy - ay * sx) + omc * dot * az;

        // Normalize
        let norm = (nx * nx + ny * ny + nz * nz).sqrt();
        let new = vec![nx / norm, ny / norm, nz / norm];
        (old, new)
    }

    fn local_energy_change_spin(&self, site: usize, old: &[f64], new: &[f64]) -> f64 {
        let mut delta_e = 0.0;

        // Outgoing bonds
        for neighbor in &self.lattice.sites[site] {
            let t = neighbor.target;
            let nx = self.spins[3 * t];
            let ny = self.spins[3 * t + 1];
            let nz = self.spins[3 * t + 2];
            delta_e -= self.j * (
                new[0] * nx + new[1] * ny + new[2] * nz
                - old[0] * nx - old[1] * ny - old[2] * nz
            );
        }

        // Incoming bonds
        for (source_idx, neighbors) in self.lattice.sites.iter().enumerate() {
            if source_idx == site { continue; }
            for neighbor in neighbors {
                if neighbor.target == site {
                    let sx = self.spins[3 * source_idx];
                    let sy = self.spins[3 * source_idx + 1];
                    let sz = self.spins[3 * source_idx + 2];
                    delta_e -= self.j * (
                        sx * new[0] + sy * new[1] + sz * new[2]
                        - sx * old[0] - sy * old[1] - sz * old[2]
                    );
                }
            }
        }

        delta_e / 2.0 // Bidirectional bonds
    }

    fn total_energy(&self) -> f64 {
        let mut energy = 0.0;
        for (site_idx, neighbors) in self.lattice.sites.iter().enumerate() {
            for neighbor in neighbors {
                let sx1 = self.spins[3 * site_idx];
                let sy1 = self.spins[3 * site_idx + 1];
                let sz1 = self.spins[3 * site_idx + 2];
                let sx2 = self.spins[3 * neighbor.target];
                let sy2 = self.spins[3 * neighbor.target + 1];
                let sz2 = self.spins[3 * neighbor.target + 2];
                energy -= self.j * (sx1 * sx2 + sy1 * sy2 + sz1 * sz2);
            }
        }
        energy / 2.0
    }

    fn spins(&self) -> &[f64] {
        &self.spins
    }

    fn spins_mut(&mut self) -> &mut [f64] {
        &mut self.spins
    }

    fn magnetization(&self) -> f64 {
        let n = self.lattice.n_sites as f64;
        let mut mx = 0.0;
        let mut my = 0.0;
        let mut mz = 0.0;
        for i in 0..self.lattice.n_sites {
            mx += self.spins[3 * i];
            my += self.spins[3 * i + 1];
            mz += self.spins[3 * i + 2];
        }
        (mx * mx + my * my + mz * mz).sqrt() / n
    }

    fn random_cluster_spin(&self, _rng: &mut impl Rng) -> f64 {
        0.0 // Not applicable for continuous spin models
    }

    fn opposite_spin(&self, spin: f64, _rng: &mut impl Rng) -> f64 {
        -spin // Reflection through origin for single component
    }
}
```

- [ ] **Step 4: Implement MonteCarlo and FromParams for HeisenbergModel**

Continue in `CMC.rs/src/models/heisenberg.rs`:

```rust
impl MonteCarlo for HeisenbergModel {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {
        // Placeholder — HeisenbergMetropolisCore handles actual sweep
    }
}

impl FromParams for HeisenbergModel {
    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let n_sites = params
            .get::<usize>("L")
            .ok_or_else(|| CarloError::InvalidConfig {
                field: "L".into(),
                reason: "System size L is required".into(),
            })?;

        let beta = params
            .get::<f64>("beta")
            .ok_or_else(|| CarloError::InvalidConfig {
                field: "beta".into(),
                reason: "Inverse temperature beta is required".into(),
            })?;

        let j = params.get::<f64>("J").unwrap_or(1.0);
        let pbc = params.get::<bool>("pbc").unwrap_or(true);
        let proposal_width = params.get::<f64>("proposal_width").unwrap_or(PI / 8.0);

        let lattice = build_chain(n_sites, pbc);
        Ok(HeisenbergModel::new(lattice, beta, j, proposal_width))
    }
}
```

- [ ] **Step 5: Implement HeisenbergMetropolisCore**

Create a new file `CMC.rs/src/algorithms/heisenberg_metropolis.rs`:

```rust
//! Metropolis-Hastings for Heisenberg model with 3D spins.

use crate::models::ModelMC;
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use crate::models::heisenberg::HeisenbergModel;
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Metropolis-Hastings single-spin flip algorithm for Heisenberg model.
pub struct HeisenbergMetropolisCore<MC: ModelMC> {
    model: MC,
}

impl<MC: ModelMC> HeisenbergMetropolisCore<MC> {
    pub fn new(model: MC) -> Self {
        HeisenbergMetropolisCore { model }
    }

    pub fn model(&self) -> &MC {
        &self.model
    }

    pub fn model_mut(&mut self) -> &mut MC {
        &mut self.model
    }
}

impl<MC: ModelMC> MonteCarlo for HeisenbergMetropolisCore<MC> {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        let n = self.model.n_sites();
        for _ in 0..n {
            let site = ctx.rng.random_range(0..n);
            let (old_spin, new_spin) = self.model.propose_flip_spin(site, &mut ctx.rng);
            let de = self.model.local_energy_change_spin(site, &old_spin, &new_spin);
            if de < 0.0 || ctx.rng.random::<f64>() < (-self.model.beta() * de).exp() {
                let spins = self.model.spins_mut();
                let dim = self.model.spin_dim();
                for d in 0..dim {
                    spins[site * dim + d] = new_spin[d];
                }
            }
        }
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let energy = self.model.total_energy();
        let magnetization = self.model.magnetization();
        ctx.measure("Energy", energy);
        ctx.measure("Energy²", energy * energy);
        ctx.measure("Magnetization", magnetization);
        ctx.measure("Magnetization²", magnetization * magnetization);
    }
}

impl<MC: ModelMC + FromParams<Rng = Xoshiro256PlusPlus>> FromParams for HeisenbergMetropolisCore<MC> {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let mc = MC::from_params(params, rng)?;
        Ok(HeisenbergMetropolisCore { model: mc })
    }
}

// Convenience: FromParams for HeisenbergModel specifically
impl FromParams for HeisenbergMetropolisCore<HeisenbergModel> {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let model = HeisenbergModel::from_params(params, rng)?;
        Ok(HeisenbergMetropolisCore { model })
    }
}
```

Wait, there's a conflict — `FromParams` would be implemented twice. Let me fix this. The blanket impl `impl<MC: ModelMC + FromParams<Rng = Xoshiro256PlusPlus>> FromParams for HeisenbergMetropolisCore<MC>` covers `HeisenbergModel` since it implements `ModelMC + FromParams`. So we don't need the second impl.

- [ ] **Step 6: Register module exports**

In `CMC.rs/src/models/mod.rs`, add:
```rust
mod heisenberg;
pub use heisenberg::HeisenbergModel;
```

In `CMC.rs/src/algorithms/mod.rs`, add:
```rust
mod heisenberg_metropolis;
pub use heisenberg_metropolis::HeisenbergMetropolisCore;
```

In `CMC.rs/src/lib.rs`, update re-exports:
```rust
pub use algorithms::{HeisenbergMetropolisCore, MetropolisCore, SWCore, WolffCore};
pub use models::{HeisenbergModel, IsingModel, IsingModel2D, ModelMC, PottsModel, XYModel};
```

- [ ] **Step 7: Write unit tests for HeisenbergModel**

Create `CMC.rs/tests/heisenberg_test.rs`:

```rust
//! Heisenberg model unit tests.

use cmc_rs::*;

#[test]
fn test_heisenberg_ground_state() {
    // All spins +z → each bond contributes -J (S·S = 1)
    let lattice = cmc_rs::lattice::build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);
    // PBC chain of 4: 4 physical bonds, each -J
    assert!((model.total_energy() - (-4.0)).abs() < 1e-10);
}

#[test]
fn test_heisenberg_spin_norm() {
    // All spins should be unit vectors
    let lattice = cmc_rs::lattice::build_chain(8, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);
    for i in 0..8 {
        let sx = model.spins()[3 * i];
        let sy = model.spins()[3 * i + 1];
        let sz = model.spins()[3 * i + 2];
        let norm = (sx * sx + sy * sy + sz * sz).sqrt();
        assert!((norm - 1.0).abs() < 1e-10, "Spin {} has norm {}", i, norm);
    }
}

#[test]
fn test_heisenberg_magnetization_all_z() {
    let lattice = cmc_rs::lattice::build_chain(4, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);
    // All +z → M = 1.0
    assert!((model.magnetization() - 1.0).abs() < 1e-10);
}

#[test]
fn test_heisenberg_metropolis_sweep() {
    let lattice = cmc_rs::lattice::build_chain(8, true);
    let model = HeisenbergModel::new(lattice, 1.0, 1.0, std::f64::consts::PI / 8.0);
    let mut core = HeisenbergMetropolisCore::new(model);
    let rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut ctx = Context::new(rng, 100);

    let initial_energy = core.model.total_energy();
    for _ in 0..100 {
        core.sweep(&mut ctx);
    }
    let final_energy = core.model.total_energy();
    // At finite temperature, energy may increase slightly from ground state
    // but should stay bounded. More importantly, spins remain unit vectors.
    assert!(final_energy > -10.0 && final_energy < 10.0);
    for i in 0..8 {
        let sx = core.model.spins()[3 * i];
        let sy = core.model.spins()[3 * i + 1];
        let sz = core.model.spins()[3 * i + 2];
        let norm = (sx * sx + sy * sy + sz * sz).sqrt();
        assert!((norm - 1.0).abs() < 0.01, "Spin {} norm = {}", i, norm);
    }
}

#[test]
fn test_heisenberg_from_params() {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let mut params = Params::new();
    params.set("L", 8usize);
    params.set("beta", 2.0f64);
    params.set("J", 1.0f64);

    let model = HeisenbergModel::from_params(&params, &mut rng).unwrap();
    assert_eq!(model.n_sites(), 8);
    assert_eq!(model.spin_dim(), 3);
}

#[test]
fn test_heisenberg_energy_scales_extensively() {
    let params = make_params(16, 1.0, 1.0, true);
    let params2 = make_params(64, 1.0, 1.0, true);

    let res1 = run_sim::<HeisenbergMetropolisCore<HeisenbergModel>>(&params, 1000, 5000, 100);
    let res2 = run_sim::<HeisenbergMetropolisCore<HeisenbergModel>>(&params2, 2000, 8000, 100);

    let e1 = res1.get("Energy").unwrap().mean / 16.0;
    let e2 = res2.get("Energy").unwrap().mean / 64.0;

    assert!(
        (e1 - e2).abs() < 0.10 * e1.abs(),
        "Per-site energy should be similar: {} vs {}",
        e1,
        e2
    );
}
```

Add helpers at the top of the test file (same as integration_test.rs):

```rust
fn run_sim<T: MonteCarlo + FromParams>(
    params: &Params,
    thermalization: u64,
    measurements: u64,
    binsize: usize,
) -> Results {
    let backend = RayonBackend::new(1);
    let config = RunConfig {
        thermalization_sweeps: thermalization,
        measurement_sweeps: measurements,
        binsize,
        base_seed: 42,
        progress_interval: 0,
        checkpoint_interval: 0,
    };
    let scheduler = Scheduler::new(backend, config);
    scheduler.run_one::<T>(params)
}

fn make_params(l: usize, beta: f64, j: f64, pbc: bool) -> Params {
    let mut p = Params::new();
    p.set("L", l);
    p.set("beta", beta);
    p.set("J", j);
    p.set("pbc", pbc);
    p
}
```

- [ ] **Step 8: Run tests and commit**

Run: `cd CMC.rs && cargo test -q`
Expected: All tests pass (existing + Heisenberg tests)

```bash
cd CMC.rs && git add src/models/heisenberg.rs src/models/mod.rs src/algorithms/heisenberg_metropolis.rs src/algorithms/mod.rs src/lib.rs tests/heisenberg_test.rs
git commit -m "feat(cmc): add Heisenberg model with 3D Metropolis

Implements classical Heisenberg model H = -J Σ S⃗_i · S⃗_j
with Rodrigues rotation proposal and unit vector normalization.
Uses propose_flip_spin / local_energy_change_spin for 3D support
in ModelMC trait with default impls for 1D models."
```

---

### Task 5: Phase B — OPSS (Optimal Phase Space Sampling) for Heisenberg

**Files:**
- Create: `CMC.rs/src/algorithms/heisenberg_opss.rs` — OPSS adaptive sampling
- Modify: `CMC.rs/src/algorithms/mod.rs` — Export
- Create: `CMC.rs/tests/opss_test.rs` — OPSS tests

OPSS is based on Alzate-Cardona et al. (2018): Gaussian move S'_i = (S_i + σF) / |S_i + σF| where F is 3D Gaussian random vector. Adaptive σ recalculated each sweep.

- [ ] **Step 1: Implement OPSS wrapper**

Create `CMC.rs/src/algorithms/heisenberg_opss.rs`:

```rust
//! OPSS (Optimal Phase Space Sampling) for Heisenberg model.
//!
//! Based on Alzate-Cardona et al. (2018): Gaussian move with adaptive σ.
//! S'_i = (S_i + σF) / |S_i + σF|, where F is 3D Gaussian random vector.

use crate::models::heisenberg::HeisenbergModel;
use crate::models::ModelMC;
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand::distributions::Distribution;
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;

/// OPSS Metropolis core for Heisenberg model.
/// Wraps HeisenbergMetropolisCore and replaces the proposal with Gaussian moves.
pub struct OPSSCore {
    model: HeisenbergModel,
    sigma: f64,
    initial_sigma: f64,
    accepted: u64,
    total: u64,
}

impl OPSSCore {
    pub fn new(model: HeisenbergModel) -> Self {
        OPSSCore {
            sigma: 60.0,
            initial_sigma: 60.0,
            accepted: 0,
            total: 0,
            model,
        }
    }

    fn gaussian_move(&mut self, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
        let dim = self.model.spin_dim();
        let old: Vec<f64> = (0..dim)
            .map(|d| self.model.spins()[site * dim + d])
            .collect();

        // Gaussian random vector
        let normal = rand_distr::StandardNormal;
        let new: Vec<f64> = old
            .iter()
            .map(|&s| s + self.sigma * rng.sample(normal))
            .collect();

        // Normalize
        let norm: f64 = new.iter().map(|&x| x * x).sum::<f64>().sqrt();
        let normalized: Vec<f64> = new.iter().map(|&x| x / norm).collect();
        (old, normalized)
    }

    fn adapt_sigma(&mut self) {
        if self.total == 0 {
            return;
        }
        let rate = self.accepted as f64 / self.total as f64;
        if rate >= 1.0 {
            // All moves accepted — reset
            self.sigma = self.initial_sigma;
        } else {
            let f = 0.5 / (1.0 - rate);
            self.sigma *= f;
        }
        if self.sigma > self.initial_sigma {
            self.sigma = self.initial_sigma;
        }
        // Reset counters
        self.accepted = 0;
        self.total = 0;
    }
}

impl MonteCarlo for OPSSCore {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        let n = self.model.n_sites();
        for _ in 0..n {
            let site = ctx.rng.random_range(0..n);
            let (old_spin, new_spin) = self.gaussian_move(site, &mut ctx.rng);
            let de = self.model.local_energy_change_spin(site, &old_spin, &new_spin);
            self.total += 1;
            if de < 0.0 || ctx.rng.random::<f64>() < (-self.model.beta() * de).exp() {
                let spins = self.model.spins_mut();
                let dim = self.model.spin_dim();
                for d in 0..dim {
                    spins[site * dim + d] = new_spin[d];
                }
                self.accepted += 1;
            }
        }
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let energy = self.model.total_energy();
        let magnetization = self.model.magnetization();
        ctx.measure("Energy", energy);
        ctx.measure("Energy²", energy * energy);
        ctx.measure("Magnetization", magnetization);
        ctx.measure("Magnetization²", magnetization * magnetization);
    }
}
```

Wait — the σ adaptation should happen per sweep, not per accepted/rejected move. Let me re-read the spec:

> **Adaptive σ**: recalculated each sweep:
> - f = 0.5 / (1 - R), where R = acceptance_rate from previous sweep
> - σ_new = σ * f
> - Initial σ = 60 (equivalent to random move)
> - If σ > 60, reset to 60 (above Tc all moves accepted)

The `accept` and `total` counters need to accumulate over a sweep, and `adapt_sigma` is called at the end of each sweep. The current design does this correctly since `sweep()` increments `total` and `accepted`, then `adapt_sigma` should be called at the end. But currently `adapt_sigma` is never called from `sweep()`. Let me fix:

At the end of `sweep()`, call `self.adapt_sigma()`.

- [ ] **Step 1 (revised): Fix OPSSCore sweep to call adapt_sigma**

```rust
fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
    let n = self.model.n_sites();
    for _ in 0..n {
        let site = ctx.rng.random_range(0..n);
        let (old_spin, new_spin) = self.gaussian_move(site, &mut ctx.rng);
        let de = self.model.local_energy_change_spin(site, &old_spin, &new_spin);
        self.total += 1;
        if de < 0.0 || ctx.rng.random::<f64>() < (-self.model.beta() * de).exp() {
            let spins = self.model.spins_mut();
            let dim = self.model.spin_dim();
            for d in 0..dim {
                spins[site * dim + d] = new_spin[d];
            }
            self.accepted += 1;
        }
    }
    self.adapt_sigma();
}
```

Actually wait — there's a dependency issue. `HeisenbergMetropolisCore` is generic over `MC: ModelMC`, but `OPSSCore` specifically wraps `HeisenbergModel`. This is fine per the design spec: "Implemented as a separate `HeisenbergAdaptiveCore<MC>` wrapper or as an alternative sweep strategy within `MetropolisCore`."

But also there's a crate dependency issue. The `rand_distr` crate may not be in Cargo.toml. Let me check.

- [ ] **Step 2: Check if rand_distr is available, add if not**

Run: `cd CMC.rs && grep rand_distr Cargo.toml`

If not found, add to `Cargo.toml` dependencies:
```toml
rand_distr = "0.4"
```

Actually, we don't strictly need `rand_distr`. We can generate Gaussian random numbers using the Box-Muller transform:

```rust
fn sample_gaussian(rng: &mut impl Rng) -> f64 {
    let u1: f64 = rng.random_range(0.0001..1.0); // Avoid log(0)
    let u2: f64 = rng.random();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}
```

This avoids adding a dependency. Let me update Step 1 to use Box-Muller instead.

- [ ] **Step 1 (final revised): Implement OPSSCore with Box-Muller**

Create `CMC.rs/src/algorithms/heisenberg_opss.rs`:

```rust
//! OPSS (Optimal Phase Space Sampling) for Heisenberg model.
//!
//! Based on Alzate-Cardona et al. (2018): Gaussian move with adaptive σ.
//! S'_i = (S_i + σF) / |S_i + σF|, where F is 3D Gaussian random vector.

use crate::models::heisenberg::HeisenbergModel;
use crate::models::ModelMC;
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use rand::Rng;
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;

/// Sample from standard normal distribution using Box-Muller transform.
fn sample_gaussian(rng: &mut impl Rng) -> f64 {
    let u1 = rng.random_range(0.0001..1.0);
    let u2 = rng.random::<f64>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// OPSS Metropolis core for Heisenberg model with adaptive σ.
pub struct OPSSCore {
    model: HeisenbergModel,
    sigma: f64,
    accepted: u64,
    total: u64,
}

impl OPSSCore {
    pub fn new(model: HeisenbergModel) -> Self {
        OPSSCore {
            sigma: 60.0,
            accepted: 0,
            total: 0,
            model,
        }
    }

    fn gaussian_move(&self, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
        let dim = self.model.spin_dim();
        let old: Vec<f64> = (0..dim)
            .map(|d| self.model.spins()[site * dim + d])
            .collect();

        let new: Vec<f64> = old
            .iter()
            .map(|&s| s + self.sigma * sample_gaussian(rng))
            .collect();

        let norm = new.iter().map(|&x| x * x).sum::<f64>().sqrt();
        let normalized: Vec<f64> = new.iter().map(|&x| x / norm).collect();
        (old, normalized)
    }

    fn adapt_sigma(&mut self) {
        if self.total == 0 { return; }
        let rate = self.accepted as f64 / self.total as f64;
        if rate >= 1.0 {
            self.sigma = 60.0;
        } else {
            let f = 0.5 / (1.0 - rate);
            self.sigma *= f;
        }
        if self.sigma > 60.0 {
            self.sigma = 60.0;
        }
        self.accepted = 0;
        self.total = 0;
    }
}

impl MonteCarlo for OPSSCore {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        let n = self.model.n_sites();
        for _ in 0..n {
            let site = ctx.rng.random_range(0..n);
            let (old_spin, new_spin) = self.gaussian_move(site, &mut ctx.rng);
            let de = self.model.local_energy_change_spin(site, &old_spin, &new_spin);
            self.total += 1;
            if de < 0.0 || ctx.rng.random::<f64>() < (-self.model.beta() * de).exp() {
                let spins = self.model.spins_mut();
                let dim = self.model.spin_dim();
                for d in 0..dim {
                    spins[site * dim + d] = new_spin[d];
                }
                self.accepted += 1;
            }
        }
        self.adapt_sigma();
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let energy = self.model.total_energy();
        let magnetization = self.model.magnetization();
        ctx.measure("Energy", energy);
        ctx.measure("Energy²", energy * energy);
        ctx.measure("Magnetization", magnetization);
        ctx.measure("Magnetization²", magnetization * magnetization);
    }
}

impl FromParams for OPSSCore {
    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let model = HeisenbergModel::from_params(params, &mut Xoshiro256PlusPlus::seed_from_u64(0))?;
        Ok(OPSSCore::new(model))
    }
}
```

Wait — the `FromParams` for `OPSSCore` uses a dummy RNG which is fine since `HeisenbergModel::from_params` doesn't actually use the RNG (it initializes all spins to +z).

- [ ] **Step 3: Export OPSSCore**

In `CMC.rs/src/algorithms/mod.rs`, add:
```rust
mod heisenberg_opss;
pub use heisenberg_opss::OPSSCore;
```

In `CMC.rs/src/lib.rs`, update:
```rust
pub use algorithms::{HeisenbergMetropolisCore, MetropolisCore, OPSSCore, SWCore, WolffCore};
```

- [ ] **Step 4: Write OPSS tests**

Create `CMC.rs/tests/opss_test.rs`:

```rust
//! OPSS (Optimal Phase Space Sampling) tests.

use cmc_rs::*;

fn run_sim<T: MonteCarlo + FromParams>(
    params: &Params,
    thermalization: u64,
    measurements: u64,
    binsize: usize,
) -> Results {
    let backend = RayonBackend::new(1);
    let config = RunConfig {
        thermalization_sweeps: thermalization,
        measurement_sweeps: measurements,
        binsize,
        base_seed: 42,
        progress_interval: 0,
        checkpoint_interval: 0,
    };
    let scheduler = Scheduler::new(backend, config);
    scheduler.run_one::<T>(params)
}

fn make_params(l: usize, beta: f64, j: f64, pbc: bool) -> Params {
    let mut p = Params::new();
    p.set("L", l);
    p.set("beta", beta);
    p.set("J", j);
    p.set("pbc", pbc);
    p
}

#[test]
fn test_opss_does_not_crash() {
    let params = make_params(8, 1.0, 1.0, true);
    let results = run_sim::<OPSSCore>(&params, 500, 2000, 50);
    let energy = results.get("Energy").unwrap();
    // Should produce reasonable energy
    assert!(energy.mean < 0.0, "Energy should be negative for ferromagnet");
}

#[test]
fn test_opss_vs_naive_agree() {
    // OPSS and naive Metropolis should give similar energies at moderate temperature
    let params = make_params(16, 1.0, 1.0, true);

    let naive = run_sim::<HeisenbergMetropolisCore<HeisenbergModel>>(&params, 2000, 5000, 100);
    let opss = run_sim::<OPSSCore>(&params, 2000, 5000, 100);

    let e_naive = naive.get("Energy").unwrap().mean / 16.0;
    let e_opss = opss.get("Energy").unwrap().mean / 16.0;

    // Should agree within ~15% (OPSS has different proposal dynamics)
    assert!(
        (e_naive - e_opss).abs() < 0.15 * e_naive.abs(),
        "OPSS per-site energy {} vs naive {} differ too much",
        e_opss,
        e_naive
    );
}

#[test]
fn test_opss_high_temperature() {
    // At β → 0, energy should be near zero
    let params = make_params(32, 0.01, 1.0, true);
    let results = run_sim::<OPSSCore>(&params, 500, 5000, 100);
    let energy = results.get("Energy").unwrap();
    assert!(
        energy.mean.abs() < 2.0,
        "High-T energy should be near zero, got {}",
        energy.mean
    );
}
```

- [ ] **Step 5: Run all tests and commit**

Run: `cd CMC.rs && cargo test -q`
Expected: All tests pass

```bash
cd CMC.rs && git add src/algorithms/heisenberg_opss.rs src/algorithms/mod.rs src/lib.rs tests/opss_test.rs
git commit -m "feat(cmc): add OPSS (Optimal Phase Space Sampling) for Heisenberg

Implements adaptive Gaussian proposal moves based on
Alzate-Cardona et al. (2018). Adaptive σ recalibrated each
sweep to maintain ~50% acceptance rate. Uses Box-Muller
transform for Gaussian sampling (no new dependencies)."
```

---

### Task 6: Integration Tests for All New Features

**Files:**
- Modify: `CMC.rs/tests/integration_test.rs` — Add tests for magnetization, Potts+cluster, Heisenberg

- [ ] **Step 1: Add Potts + SW/Wolff agreement test**

Append to `CMC.rs/tests/integration_test.rs`:

```rust
// ─── Potts model with cluster algorithms ────────────────────────────────────

#[test]
fn test_potts_sw_wolff_agree() {
    // SW and Wolff should agree for 3-state Potts model
    let mut params = make_params(32, 0.5, 1.0, true);
    params.set("q", 3usize);

    let sw = run_simulation::<SWCore<PottsModel>>(&params, 1000, 5000, 100);
    let wol = run_simulation::<WolffCore<PottsModel>>(&params, 1000, 5000, 100);

    let e_sw = sw.get("Energy").unwrap().mean;
    let e_wol = wol.get("Energy").unwrap().mean;

    assert!(
        (e_sw - e_wol).abs() < 0.15 * e_sw.abs(),
        "SW energy {} vs Wolff energy {} differ too much for Potts",
        e_sw,
        e_wol
    );
}

#[test]
fn test_potts_metropolis_sw_agree() {
    // Metropolis and SW should agree for 3-state Potts
    let mut params = make_params(32, 0.5, 1.0, true);
    params.set("q", 3usize);

    let met = run_simulation::<MetropolisCore<PottsModel>>(&params, 2000, 8000, 100);
    let sw = run_simulation::<SWCore<PottsModel>>(&params, 500, 3000, 100);

    let e_met = met.get("Energy").unwrap().mean;
    let e_sw = sw.get("Energy").unwrap().mean;

    assert!(
        (e_met - e_sw).abs() < 0.15 * e_met.abs(),
        "Metropolis energy {} vs SW energy {} differ too much for Potts",
        e_met,
        e_sw
    );
}
```

- [ ] **Step 2: Add magnetization physical behavior tests**

Append to `CMC.rs/tests/integration_test.rs`:

```rust
// ─── Magnetization tests ────────────────────────────────────────────────────

#[test]
fn test_ising_magnetization_temperature_dependence() {
    // Magnetization should decrease with temperature
    let beta_low = 5.0;  // Low T → ordered
    let beta_high = 0.1; // High T → disordered

    let params_low = make_params(32, beta_low, 1.0, true);
    let params_high = make_params(32, beta_high, 1.0, true);

    let res_low = run_simulation::<MetropolisCore<IsingModel>>(&params_low, 2000, 5000, 100);
    let res_high = run_simulation::<MetropolisCore<IsingModel>>(&params_high, 500, 3000, 100);

    let m_low = res_low.get("Magnetization").unwrap().mean;
    let m_high = res_high.get("Magnetization").unwrap().mean;

    assert!(
        m_low > m_high,
        "Low-T magnetization {} should be > high-T {}",
        m_low,
        m_high
    );
}

#[test]
fn test_magnetization_squared_fluctuation() {
    // Magnetization² should be measurable and positive
    let params = make_params(16, 1.0, 1.0, true);
    let results = run_simulation::<MetropolisCore<IsingModel>>(&params, 500, 3000, 100);

    let m2 = results.get("Magnetization²").unwrap();
    assert!(m2.mean > 0.0, "Magnetization² should be positive");
}

#[test]
fn test_energy_squared_fluctuation() {
    // Energy² should be measurable and > Energy² (by Jensen's inequality)
    let params = make_params(16, 1.0, 1.0, true);
    let results = run_simulation::<MetropolisCore<IsingModel>>(&params, 500, 3000, 100);

    let e = results.get("Energy").unwrap();
    let e2 = results.get("Energy²").unwrap();

    // E[X²] >= (E[X])² by Jensen's inequality
    assert!(
        e2.mean > e.mean * e.mean * 0.9, // 10% tolerance
        "Energy² {} should be >= (Energy)² {}",
        e2.mean,
        e.mean * e.mean
    );
}
```

- [ ] **Step 3: Add Heisenberg integration tests**

Append to `CMC.rs/tests/integration_test.rs`:

```rust
// ─── Heisenberg model integration ───────────────────────────────────────────

#[test]
fn test_heisenberg_ground_state_energy() {
    // At T → 0, Heisenberg should approach ground state: E = -N_bonds * J
    let params = make_params(16, 50.0, 1.0, true);
    let results = run_sim_for_heisenberg::<HeisenbergMetropolisCore<HeisenbergModel>>(
        &params, 5000, 5000, 100
    );
    let energy = results.get("Energy").unwrap();
    // PBC chain of 16: 16 physical bonds, ground state E = -16
    assert!(
        energy.mean > -18.0 && energy.mean < -14.0,
        "Low-T Heisenberg energy {:.2} should be near -16",
        energy.mean
    );
}

#[test]
fn test_heisenberg_high_temperature() {
    // At β → 0, spins are random → ⟨E⟩ → 0
    let params = make_params(32, 0.01, 1.0, true);
    let results = run_sim_for_heisenberg::<HeisenbergMetropolisCore<HeisenbergModel>>(
        &params, 500, 5000, 100
    );
    let energy = results.get("Energy").unwrap();
    assert!(
        energy.mean.abs() < 2.0,
        "High-T Heisenberg energy should be near zero, got {}",
        energy.mean
    );
}

#[test]
fn test_opss_and_naive_agree_heisenberg() {
    let params = make_params(16, 1.0, 1.0, true);

    let naive = run_sim_for_heisenberg::<HeisenbergMetropolisCore<HeisenbergModel>>(
        &params, 2000, 5000, 100
    );
    let opss = run_sim_for_heisenberg::<OPSSCore>(&params, 2000, 5000, 100);

    let e_naive = naive.get("Energy").unwrap().mean;
    let e_opss = opss.get("Energy").unwrap().mean;

    assert!(
        (e_naive - e_opss).abs() < 0.15 * e_naive.abs(),
        "OPSS energy {} vs naive {} differ too much",
        e_opss,
        e_naive
    );
}
```

Add the helper at the top of the file (after existing helpers):

```rust
fn run_sim_for_heisenberg<T: MonteCarlo + FromParams>(
    params: &Params,
    thermalization: u64,
    measurements: u64,
    binsize: usize,
) -> Results {
    let backend = RayonBackend::new(1);
    let config = RunConfig {
        thermalization_sweeps: thermalization,
        measurement_sweeps: measurements,
        binsize,
        base_seed: 42,
        progress_interval: 0,
        checkpoint_interval: 0,
    };
    let scheduler = Scheduler::new(backend, config);
    scheduler.run_one::<T>(params)
}
```

- [ ] **Step 4: Run all tests and final commit**

Run: `cd CMC.rs && cargo test -q 2>&1`
Expected: All tests pass (50+ tests total)

```bash
cd CMC.rs && git add tests/integration_test.rs
git commit -m "test(cmc): add integration tests for magnetization, Potts cluster, Heisenberg

Tests cover:
- Potts model with SW/Wolff cluster algorithms (agreement test)
- Magnetization temperature dependence (Ising)
- Energy² and Magnetization² fluctuation measurements
- Heisenberg ground state, high-T limit, and OPSS vs naive agreement"
```

---

## Task Dependencies

```
Task 1 (magnetization/snapshot trait)
    └── Task 2 (extend algorithm measure) — needs magnetization() on ModelMC

Task 3 (cluster spin methods) — independent of Tasks 1, 2

Task 4 (HeisenbergModel) — needs Task 3's propose_flip_spin/local_energy_change_spin trait methods
    └── Task 5 (OPSS) — needs HeisenbergModel

Task 6 (integration tests) — needs all of the above
```

Recommended execution order: **1 → 2, 3 (parallel) → 4 → 5 → 6**

## Summary of File Changes

| File | Action | Task |
|------|--------|------|
| `src/models/mod.rs` | Modify | 1, 3, 4 |
| `src/models/ising.rs` | Modify | 1, 3 |
| `src/models/ising_2d.rs` | Modify | 1, 3 |
| `src/models/potts.rs` | Modify | 1, 3 |
| `src/models/xy.rs` | Modify | 1, 3 |
| `src/models/heisenberg.rs` | **Create** | 4 |
| `src/algorithms/metropolis.rs` | Modify | 2 |
| `src/algorithms/wolff.rs` | Modify | 2, 3 |
| `src/algorithms/swendsen_wang.rs` | Modify | 2, 3 |
| `src/algorithms/heisenberg_metropolis.rs` | **Create** | 4 |
| `src/algorithms/heisenberg_opss.rs` | **Create** | 5 |
| `src/algorithms/mod.rs` | Modify | 4, 5 |
| `src/lib.rs` | Modify | 4, 5 |
| `src/algorithms/mod.rs` | Modify | 4, 5 |
| `tests/integration_test.rs` | Modify | 6 |
| `tests/heisenberg_test.rs` | **Create** | 4 |
| `tests/opss_test.rs` | **Create** | 5 |
