# OPSS Generic + XY Vector Model Design

> **Goal:** Make OPSS applicable to all continuous spin models (O(N) ≥ 2), and unify XY model to vector storage like Heisenberg.

---

## I. Problem

### XY stores angles, Heisenberg stores vectors — inconsistency

- `XYModel`: `spins[i] = θ_i` (scalar, dim=1)
- `HeisenbergModel`: `spins[3*i..3*i+3] = (x_i, y_i, z_i)` (vector, dim=3)

This forces two separate Metropolis implementations (before the recent refactor), and OPSS cannot apply to XY because Gaussian perturbation + normalization requires vector representation.

### OPSSStrategy is hardcoded to HeisenbergModel

`OPSSStrategy` implements `ProposalStrategy<HeisenbergModel>` — not reusable for XY or any future O(N) model.

---

## II. Design

### A. XY Model: Vector Storage (cos θ, sin θ)

**Before:**
```rust
// spins: [θ₀, θ₁, θ₂, ...]  — scalar angles
spin_dim() = 1
total_energy: -J * cos(θ_i - θ_j)
```

**After:**
```rust
// spins: [x₀, y₀, x₁, y₁, ...]  — unit 2D vectors
// where (x_i, y_i) = (cos θ_i, sin θ_i)
spin_dim() = 2
total_energy: -J * (x_i * x_j + y_i * y_j)  // dot product = cos(θ_i - θ_j)
```

**Energy equivalence:** `x_i * x_j + y_i * y_j = cos(θ_i) * cos(θ_j) + sin(θ_i) * sin(θ_j) = cos(θ_i - θ_j)`. Same Hamiltonian, different representation.

### B. OPSSStrategy: Generic over Model

**Before:**
```rust
impl ProposalStrategy<HeisenbergModel> for OPSSStrategy { ... }
```

**After:**
```rust
impl<MC: ModelMC> ProposalStrategy<MC> for OPSSStrategy {
    fn propose_flip(&self, model: &MC, site: usize, rng: &mut impl Rng) -> (Vec<f64>, Vec<f64>) {
        let dim = model.spin_dim();
        // ... Gaussian perturbation + normalization (already dimension-agnostic)
    }
}
```

OPSS works for any model with `spin_dim() >= 2`. For dim=1 (Ising), normalization degenerates to ±1, making OPSS equivalent to a random sign flip — technically works but physically useless. No compile-time guard needed; users should know what they're doing.

### C. XY ModelMC Implementation Changes

| Method | Old (angle) | New (vector) |
|--------|-------------|--------------|
| `spin_dim()` | 1 | 2 |
| `propose_flip_spin` | Rodrigues on 1D (degenerate) | Rodrigues 2D rotation matrix |
| `local_energy_change_spin` | wraps `local_energy_change` (angles) | Direct dot-product computation |
| `total_energy` | `-J * cos(θ_i - θ_j)` | `-J * (x_i*x_j + y_i*y_j)` |
| `magnetization` | `\|Σ(cos θ, sin θ)\| / N` | `\|Σ(x_i, y_i)\| / N` |
| `propose_flip` | angle perturbation | wraps `propose_flip_spin` |
| `local_energy_change` | angle-based | wraps `local_energy_change_spin` |
| `random_cluster_spin` | `0.0` (N/A) | `0.0` (N/A) |
| `opposite_spin` | θ + π | (-x, -y) negation |

**Rodrigues 2D rotation:** For a unit vector `(x, y)` rotated by angle δ:
```
x' = x*cos(δ) - y*sin(δ)
y' = x*sin(δ) + y*cos(δ)
```
This is just the standard 2D rotation matrix. No normalization needed (rotation preserves unit length), but we normalize anyway for numerical stability.

---

## III. Files to Change

### 1. `CMC.rs/src/models/xy.rs` — Complete rewrite to vector storage

- `spins: Vec<f64>` stores `[x₀, y₀, x₁, y₁, ...]`
- `spin_dim()` returns 2
- All energy, magnetization, proposal methods use vector form
- Scalar `propose_flip` / `local_energy_change` wrap vector methods
- `FromParams` still accepts `proposal_width` parameter (same interface)

### 2. `CMC.rs/src/algorithms/opss_strategy.rs` — Generic implementation

- Change `impl ProposalStrategy<HeisenbergModel>` to `impl<MC: ModelMC> ProposalStrategy<MC>`
- The Gaussian perturbation + normalization code already works for any dimension
- No structural changes needed

### 3. `CMC.rs/tests/integration_test.rs` — Update XY tests

- Ground state test: all spins (1, 0) instead of angle 0
- Energy consistency test: already added in last session, verify still works
- Any other XY-specific tests

### 4. `CMC.rs/tests/heisenberg_test.rs` — No changes needed

### 5. No changes to `lib.rs`, `mod.rs`, or algorithm files

---

## IV. Backwards Compatibility

- `XYModel::new()` constructor changes signature: takes same lattice/beta/j, no proposal_width stored as angle anymore
- `FromParams` interface unchanged (still accepts `proposal_width` for Metropolis rotation width)
- Tests that construct XY models with angles need updating
- **No changes to the public API of algorithms or traits**

---

## V. Testing Plan

1. XY ground state: all spins (1,0) → energy = -J * n_bonds (same as before)
2. XY spin norm: verify all (x_i, y_i) are unit vectors after sweeps
3. XY magnetization: all aligned → M=1, uniform spread → M≈0
4. XY energy extensivity: per-site energy consistent across system sizes
5. XY local vs total energy consistency: already tested in last session
6. XY + OPSS: run OPSS on XY model, verify energy matches Metropolis
7. XY + Wolff/SW: `opposite_spin` = negation, `random_cluster_spin` unchanged

---

## VI. Risks

- **XY tests may fail** if angle-based assumptions are hardcoded anywhere (e.g., in integration tests that create XY models with specific angles)
- **OPSS on XY (dim=2)** — the algorithm is physically correct (Gaussian perturbation on S^1 + normalization), but may have different autocorrelation properties than Metropolis rotation. Need a test to verify agreement.

---

*Design established 2026-04-10.*
