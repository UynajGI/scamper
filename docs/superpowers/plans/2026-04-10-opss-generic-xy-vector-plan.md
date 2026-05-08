# OPSS Generic + XY Vector Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make OPSS applicable to all continuous spin models (O(N) ≥ 2) by genericizing OPSSStrategy, and unify XY model to vector (cos θ, sin θ) storage like Heisenberg.

**Architecture:** XYModel changes from angle storage to 2D vector storage. OPSSStrategy becomes generic over any ModelMC, using the model's spin_dim() for Gaussian perturbation + normalization. The same code path works for XY (dim=2), Heisenberg (dim=3), and any future O(N) model.

**Tech Stack:** Rust, Carlo.rs framework traits (ModelMC, ProposalStrategy, MetropolisCore)

---

### Task 1: Rewrite XYModel to vector (x, y) storage

**Files:**
- Modify: `CMC.rs/src/models/xy.rs`
- Test: `CMC.rs/tests/integration_test.rs` (update XY test)

- [ ] **Step 1: Write failing tests for XY vector model**

Add to `CMC.rs/tests/integration_test.rs`:

```rust
#[test]
fn test_xy_vector_ground_state() {
    // All spins (1, 0) → dot product = 1 for each bond
    // 4-site PBC: 4 physical bonds, energy = -J * 4 = -4.0
    let lattice = cmc_rs::lattice::build_chain(4, true);
    let model = XYModel::new(lattice, 1.0, 1.0);
    assert!((model.spin_dim() - 2).abs() < 1e-10, "XY spin_dim should be 2");
    // Initial state: all spins (1, 0)
    assert!((model.spins()[0] - 1.0).abs() < 1e-10);
    assert!((model.spins()[1] - 0.0).abs() < 1e-10);
    assert!((model.total_energy() - (-4.0)).abs() < 1e-10);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd CMC.rs && cargo test test_xy_vector_ground_state`
Expected: FAIL — XYModel currently stores angles (spin_dim=1), not vectors

- [ ] **Step 3: Rewrite XYModel to vector storage**

Replace the entire content of `CMC.rs/src/models/xy.rs`:

```rust
//! Classical XY model: H = -J Σ S⃗_i · S⃗_j
//!
//! S⃗_i = (cos θ_i, sin θ_i) — unit vector on S¹.
//! Spins stored as 2 consecutive f64 values per site: [x0, y0, x1, y1, ...]

use crate::lattice::{LatticeMC, build_chain};
use crate::{CarloError, Context, FromParams, MonteCarlo, Params};
use crate::models::ModelMC;
use rand::Rng;
use rand::RngExt;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::f64::consts::PI;

/// Classical XY model: H = -J Σ S⃗_i · S⃗_j
///
/// Each spin is a unit 2D vector (cos θ, sin θ).
pub struct XYModel {
    lattice: crate::lattice::Lattice,
    beta: f64,
    j: f64,
    spins: Vec<f64>,          // 2 * n_sites: [x0, y0, x1, y1, ...]
    proposal_width: f64,      // Angular perturbation range [-δ, δ]
}

impl XYModel {
    pub fn new(lattice: crate::lattice::Lattice, beta: f64, j: f64, proposal_width: f64) -> Self {
        let n_sites = lattice.n_sites;
        let mut spins = vec![0.0; 2 * n_sites];
        for i in 0..n_sites {
            spins[2 * i] = 1.0; // x = 1.0 (cos 0)
            spins[2 * i + 1] = 0.0; // y = 0.0 (sin 0)
        }
        XYModel {
            lattice, beta, j, spins, proposal_width,
        }
    }

    pub fn proposal_width(&self) -> f64 {
        self.proposal_width
    }
}

impl LatticeMC for XYModel {
    fn lattice(&self) -> &crate::lattice::Lattice {
        &self.lattice
    }
}

impl ModelMC for XYModel {
    fn spin_dim(&self) -> usize {
        2
    }

    fn coupling(&self) -> f64 {
        self.j
    }

    fn beta(&self) -> f64 {
        self.beta
    }

    fn propose_flip(&self, site: usize, rng: &mut impl Rng) -> (f64, f64) {
        let (old, new) = self.propose_flip_spin(site, rng);
        (old[0], new[0]) // Fallback: return x-component
    }

    fn local_energy_change(&self, site: usize, _old: f64, _new: f64) -> f64 {
        // Not applicable for vector model; use local_energy_change_spin
        unimplemented!("Use local_energy_change_spin for XYModel")
    }

    fn propose_flip_spin(&self, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
        let x = self.spins[2 * site];
        let y = self.spins[2 * site + 1];

        // 2D rotation by random angle in [-proposal_width, proposal_width]
        let angle = rng.random_range(-self.proposal_width..self.proposal_width);
        let c = angle.cos();
        let s = angle.sin();

        let nx = c * x - s * y;
        let ny = s * x + c * y;

        // Normalize for numerical stability
        let norm = (nx * nx + ny * ny).sqrt();
        (vec![x, y], vec![nx / norm, ny / norm])
    }

    fn local_energy_change_spin(&self, site: usize, old: &[f64], new: &[f64]) -> f64 {
        let mut delta_e = 0.0;
        for neighbor in &self.lattice.sites[site] {
            let t = neighbor.target;
            let nx = self.spins[2 * t];
            let ny = self.spins[2 * t + 1];
            delta_e -= self.j * (
                new[0] * nx + new[1] * ny
                - old[0] * nx - old[1] * ny
            );
        }
        delta_e
    }

    fn total_energy(&self) -> f64 {
        let mut energy = 0.0;
        for (site_idx, neighbors) in self.lattice.sites.iter().enumerate() {
            for neighbor in neighbors {
                let sx1 = self.spins[2 * site_idx];
                let sy1 = self.spins[2 * site_idx + 1];
                let sx2 = self.spins[2 * neighbor.target];
                let sy2 = self.spins[2 * neighbor.target + 1];
                energy -= self.j * (sx1 * sx2 + sy1 * sy2);
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
        for i in 0..self.lattice.n_sites {
            mx += self.spins[2 * i];
            my += self.spins[2 * i + 1];
        }
        (mx * mx + my * my).sqrt() / n
    }

    fn random_cluster_spin(&self, _rng: &mut impl Rng) -> f64 {
        0.0 // Not applicable for continuous spin models
    }

    fn opposite_spin(&self, spin: f64, _rng: &mut impl Rng) -> f64 {
        -spin // Negation: (x,y) → (-x,-y)
    }
}

impl MonteCarlo for XYModel {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, _ctx: &mut Context<Self::Rng>) {
        // Placeholder — engines handle the actual sweep
    }
}

impl FromParams for XYModel {
    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let n_sites = params.get::<usize>("L").ok_or_else(|| CarloError::InvalidConfig {
            field: "L".into(),
            reason: "System size L is required".into(),
        })?;
        let beta = params.get::<f64>("beta").ok_or_else(|| CarloError::InvalidConfig {
            field: "beta".into(),
            reason: "Inverse temperature beta is required".into(),
        })?;
        let j = params.get::<f64>("J").unwrap_or(1.0);
        let pbc = params.get::<bool>("pbc").unwrap_or(true);
        let proposal_width = params.get::<f64>("proposal_width").unwrap_or(PI / 4.0);
        let lattice = build_chain(n_sites, pbc);
        Ok(XYModel::new(lattice, beta, j, proposal_width))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::build_chain;
    use rand_xoshiro::rand_core::SeedableRng;

    #[test]
    fn test_xy_ground_state() {
        // All spins (1, 0) → each bond dot product = 1
        // PBC chain of 4: 4 physical bonds, energy = -4.0
        let lattice = build_chain(4, true);
        let model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
        assert!((model.total_energy() - (-4.0)).abs() < 1e-10);
    }

    #[test]
    fn test_xy_spin_norm() {
        let lattice = build_chain(8, true);
        let model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
        for i in 0..8 {
            let sx = model.spins()[2 * i];
            let sy = model.spins()[2 * i + 1];
            let norm = (sx * sx + sy * sy).sqrt();
            assert!((norm - 1.0).abs() < 1e-10, "Spin {} has norm {}", i, norm);
        }
    }

    #[test]
    fn test_xy_magnetization_all_aligned() {
        let lattice = build_chain(4, true);
        let model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
        assert!((model.magnetization() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_xy_from_params() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let mut params = Params::new();
        params.set("L", 8usize);
        params.set("beta", 2.0f64);
        params.set("J", 1.0f64);

        let model = XYModel::from_params(&params, &mut rng).unwrap();
        assert_eq!(model.n_sites(), 8);
        assert_eq!(model.spin_dim(), 2);
    }
}
```

- [ ] **Step 4: Run XY unit tests to verify they pass**

Run: `cd CMC.rs && cargo test --lib xy`
Expected: 4 XY unit tests pass

- [ ] **Step 5: Update the XY energy consistency test in integration_test.rs**

Replace the existing `test_xy_local_vs_total_energy_consistency` in `CMC.rs/tests/integration_test.rs`:

```rust
#[test]
fn test_xy_local_vs_total_energy_consistency() {
    use std::f64::consts::PI;

    let lattice = cmc_rs::lattice::build_chain(8, true);
    let mut model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
    // Set non-trivial spin configuration (spread angles)
    let angles = [0.0, PI / 4.0, PI / 2.0, 3.0 * PI / 4.0, PI, 5.0 * PI / 4.0, 3.0 * PI / 2.0, 7.0 * PI / 4.0];
    for (i, &a) in angles.iter().enumerate() {
        model.spins_mut()[2 * i] = a.cos();
        model.spins_mut()[2 * i + 1] = a.sin();
    }

    for site in 0..8 {
        let old = vec![model.spins()[2 * site], model.spins()[2 * site + 1]];
        let new = vec![-old[0], -old[1]]; // Flip by π (negation)

        let de_local = model.local_energy_change_spin(site, &old, &new);

        let e_before = model.total_energy();
        model.spins_mut()[2 * site] = new[0];
        model.spins_mut()[2 * site + 1] = new[1];
        let e_after = model.total_energy();
        let de_global = e_after - e_before;

        assert!(
            (de_local - de_global).abs() < 1e-8,
            "Site {}: local ΔE={} vs global ΔE={}",
            site, de_local, de_global
        );

        // Restore
        model.spins_mut()[2 * site] = old[0];
        model.spins_mut()[2 * site + 1] = old[1];
    }
}
```

Also add a single-site XY boundary test in the single-site section:

```rust
#[test]
fn test_single_site_xy() {
    let lattice = cmc_rs::lattice::build_chain(1, false);
    let model = XYModel::new(lattice, 1.0, 1.0, PI / 4.0);
    assert!((model.total_energy() - 0.0).abs() < 1e-10, "Single-site XY energy should be 0");
    assert!((model.magnetization() - 1.0).abs() < 1e-10, "Single-site XY magnetization should be 1");
}
```

- [ ] **Step 6: Run all tests**

Run: `cd CMC.rs && cargo test`
Expected: All tests pass (90+)

- [ ] **Step 7: Commit**

Run:
```bash
jj commit -m "refactor(cmc): rewrite XYModel to vector (x,y) storage

- Store spins as (cos θ, sin θ) 2D vectors instead of scalar angles
- spin_dim() = 2, consistent with Heisenberg (dim=3)
- Energy uses dot product: -J*(x_i*x_j + y_i*y_j) = -J*cos(θ_i-θ_j)
- propose_flip_spin uses 2D rotation matrix
- local_energy_change_spin uses vector dot product
- Update integration test for vector storage
- Add single-site XY boundary test"
```

---

### Task 2: Genericize OPSSStrategy

**Files:**
- Modify: `CMC.rs/src/algorithms/opss_strategy.rs`
- Test: `CMC.rs/tests/integration_test.rs` (add XY+OPSS test)

- [ ] **Step 1: Write failing test for OPSS on XY**

Add to `CMC.rs/tests/integration_test.rs`:

```rust
#[test]
fn test_opss_on_xy_agrees_with_metropolis() {
    // OPSS and Metropolis should give similar energies for XY model
    let params = make_params(16, 1.0, 1.0, true);

    let met = run_simulation::<MetropolisCore<XYModel>>(&params, 2000, 5000, 100);
    let opss = run_simulation::<MetropolisCore<XYModel, OPSSStrategy>>(&params, 2000, 5000, 100);

    let e_met = met.get("Energy").unwrap().mean / 16.0;
    let e_opss = opss.get("Energy").unwrap().mean / 16.0;

    assert!(
        (e_met - e_opss).abs() < 0.15 * e_met.abs(),
        "OPSS per-site energy {} vs Metropolis {} differ too much for XY",
        e_opss, e_met
    );
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd CMC.rs && cargo test test_opss_on_xy_agrees_with_metropolis`
Expected: FAIL — OPSSStrategy is hardcoded to `ProposalStrategy<HeisenbergModel>`

- [ ] **Step 3: Make OPSSStrategy generic over ModelMC**

Replace the entire content of `CMC.rs/src/algorithms/opss_strategy.rs`:

```rust
//! OPSS (Optimal Phase Space Sampling) strategy for continuous spin models.
//!
//! Gaussian move with adaptive sigma. S'_i = (S_i + sigma * F) / |S_i + sigma * F|
//! where F is a Gaussian random vector. Based on Alzate-Cardona et al. (2018).
//! Applicable to any O(N) model with spin_dim() >= 2.

use crate::algorithms::proposal_strategy::ProposalStrategy;
use crate::models::ModelMC;
use rand::Rng;
use rand::RngExt;

/// Sample from standard normal distribution using Box-Muller transform.
fn sample_gaussian(rng: &mut impl Rng) -> f64 {
    let u1: f64 = rng.random_range(0.0001..1.0);
    let u2 = rng.random::<f64>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// OPSS proposal strategy for continuous spin models (O(N), N >= 2).
/// Uses Gaussian perturbation + normalization to unit sphere, with adaptive sigma.
pub struct OPSSStrategy {
    sigma: f64,
    accepted: u64,
    total: u64,
    sigma_min: f64,
    sigma_max: f64,
}

impl OPSSStrategy {
    pub fn new(initial_sigma: f64) -> Self {
        OPSSStrategy {
            sigma: initial_sigma,
            accepted: 0,
            total: 0,
            sigma_min: 1e-6,
            sigma_max: 60.0,
        }
    }

    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    fn adapt_sigma(&mut self) {
        if self.total == 0 {
            return;
        }
        let rate = self.accepted as f64 / self.total as f64;
        if rate >= 1.0 {
            self.sigma = self.sigma_max;
        } else {
            let f = 0.5 / (1.0 - rate);
            self.sigma *= f;
        }
        self.sigma = self.sigma.clamp(self.sigma_min, self.sigma_max);
        self.accepted = 0;
        self.total = 0;
    }
}

impl<MC: ModelMC> ProposalStrategy<MC> for OPSSStrategy {
    fn propose_flip(
        &mut self,
        model: &MC,
        site: usize,
        rng: &mut impl Rng,
    ) -> (Vec<f64>, Vec<f64>) {
        let dim = model.spin_dim();
        let old: Vec<f64> = (0..dim)
            .map(|d| model.spins()[site * dim + d])
            .collect();

        let new: Vec<f64> = old
            .iter()
            .map(|&s| s + self.sigma * sample_gaussian(rng))
            .collect();

        let norm = new.iter().map(|&x| x * x).sum::<f64>().sqrt();
        let normalized: Vec<f64> = new.iter().map(|&x| x / norm).collect();
        (old, normalized)
    }

    fn compute_delta_e(
        &self,
        model: &MC,
        site: usize,
        old_spin: &[f64],
        new_spin: &[f64],
    ) -> f64 {
        model.local_energy_change_spin(site, old_spin, new_spin)
    }

    fn adapt_after_sweep(&mut self, _model: &mut MC) {
        self.adapt_sigma();
    }
}
```

- [ ] **Step 4: Update module exports**

Update `CMC.rs/src/algorithms/mod.rs` — no changes needed (OPSSStrategy already exported).

Update `CMC.rs/src/lib.rs` — check that `XYModel` already has `proposal_width` parameter in its constructor call if needed.

Actually, `XYModel::new()` now takes 4 params (lattice, beta, j, proposal_width). Check `FromParams` — it already calls `XYModel::new(lattice, beta, j, proposal_width)` in the new code.

No export changes needed.

- [ ] **Step 5: Run all tests**

Run: `cd CMC.rs && cargo test`
Expected: All 91+ tests pass (including new `test_opss_on_xy_agrees_with_metropolis`)

- [ ] **Step 6: Commit**

Run:
```bash
jj commit -m "refactor(cmc): genericize OPSSStrategy for all O(N) models

- Change from ProposalStrategy<HeisenbergModel> to ProposalStrategy<MC: ModelMC>
- Gaussian perturbation + normalization works for any spin_dim() >= 2
- OPSS now works for XY (dim=2), Heisenberg (dim=3), and future O(N) models
- Add test: OPSS agrees with Metropolis for XY model
- Update docstrings to reflect generic applicability"
```

---

### Task 3: XY model integration tests (end-to-end)

**Files:**
- Test: `CMC.rs/tests/integration_test.rs`

- [ ] **Step 1: Add XY Metropolis simulation test**

Add to `CMC.rs/tests/integration_test.rs`:

```rust
#[test]
fn test_xy_metropolis_energy_extensive() {
    let params = make_params(16, 1.0, 1.0, true);
    let params2 = make_params(64, 1.0, 1.0, true);

    let res1 = run_simulation::<MetropolisCore<XYModel>>(&params, 1000, 5000, 100);
    let res2 = run_simulation::<MetropolisCore<XYModel>>(&params2, 2000, 8000, 100);

    let e1 = res1.get("Energy").unwrap().mean / 16.0;
    let e2 = res2.get("Energy").unwrap().mean / 64.0;

    assert!(
        (e1 - e2).abs() < 0.10 * e1.abs(),
        "Per-site energy should be similar: {} vs {}",
        e1, e2
    );
}

#[test]
fn test_xy_high_temperature() {
    let params = make_params(32, 0.01, 1.0, true);
    let results = run_simulation::<MetropolisCore<XYModel>>(&params, 500, 5000, 100);
    let energy = results.get("Energy").unwrap();
    assert!(
        energy.mean.abs() < 2.0,
        "High-T energy should be near zero, got {}",
        energy.mean
    );
}
```

- [ ] **Step 2: Run all tests**

Run: `cd CMC.rs && cargo test`
Expected: All 93+ tests pass

- [ ] **Step 3: Commit**

Run:
```bash
jj commit -m "test(cmc): add XY model integration tests

- Add XY Metropolis energy extensivity test (16 vs 64 sites)
- Add XY high-temperature limit test (energy → 0)
- All XY tests verify vector storage works end-to-end"
```

---
