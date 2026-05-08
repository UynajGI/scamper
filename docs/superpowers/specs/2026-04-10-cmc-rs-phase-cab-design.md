# CMC.rs Phase C+A+B Design

> Implementation plan for measurements, algorithm fixes, and Heisenberg model.
>
> Date: 2026-04-10

---

## Phase C: Basic Measurements + Configuration Snapshots

### C.1 Measurement Extension

Add `magnetization()` to `ModelMC` trait. Each model implements its own definition:

| Model | Magnetization |
|-------|---------------|
| Ising | \|Σ s_i\| / N, s_i = ±1 |
| Potts | (q · max(n_0, ..., n_{q-1}) - N) / (N · (q-1)) |
| XY | \|Σ (cos θ_i, sin θ_i)\| / N |
| Heisenberg | \|Σ S⃗_i\| / N |

Extend each `*Core::measure()` to also record:
- `Magnetization` — mean magnetization
- `Magnetization²` — for susceptibility via fluctuations
- `Energy²` — for specific heat via fluctuations

### C.2 Configuration Snapshots

Add `snapshot()` method to `ModelMC` returning `Vec<f64>` (raw spin configuration).
Recorded via `ctx.measure("Snapshot", ...)` at configurable intervals.
No new file format — uses existing Carlo.rs measurement pipeline.

---

## Phase A: Fix SW/Wolff for Non-Ising Models

### A.1 SWCore Cluster Spin Assignment

Current code hardcodes `±1.0` for cluster spin assignment (line 90 of swendsen_wang.rs).

Add `random_cluster_spin(rng: &mut impl Rng) -> f64` to `ModelMC`:
- Ising: returns `+1.0` or `-1.0` with 50% probability
- Potts: returns uniform random integer in `0..q` as f64
- XY/Heisenberg: not applicable (SW only valid for discrete spin models with Z_q symmetry)

### A.2 WolffCore Cluster Flip

Current code flips cluster by setting `new_spin = -seed_spin` (line 59 of wolff.rs).

Add `opposite_spin(spin: f64) -> f64` to `ModelMC`:
- Ising: returns `-spin`
- Potts: returns uniform random from states ≠ current
- XY/Heisenberg: not applicable (Wolff only valid for O(n) models with reflection symmetry — future extension)

### A.3 Algorithm-Model Compatibility

SW and Wolff only work for models with discrete symmetry (Ising, Potts).
XY and Heisenberg require continuous cluster algorithms (not implemented).
No runtime guard needed — if a model's `random_cluster_spin` is nonsensical, the simulation produces wrong results (user responsibility).

---

## Phase B: Heisenberg Model

### B.1 Model Definition

```
H = -J Σ S⃗_i · S⃗_j
```

S⃗_i = (sin θ cos φ, sin θ sin φ, cos θ) — unit vector on S².

### B.2 Storage

Spin stored as three consecutive f64 values: `spins[3*i]`, `spins[3*i+1]`, `spins[3*i+2]`.
`spins()` returns flat `Vec<f64>` of length `3 * n_sites`.
`spin_dim() -> 3`.

### B.3 Propose Flip (朴素 Metropolis)

Small angular perturbation: generate random rotation axis on S², rotate by angle
uniform in [-δ, δ]. Re-normalize after perturbation. Consistent with XYModel's
approach but on the sphere.

### B.4 OPSS (Optimal Phase Space Sampling) — Optional

Based on Alzate-Cardona et al. (2018) [J. Phys.: Condens. Matter, doi:10.1088/1361-648X/aaf852].

**Gaussian move**: S'_i = (S_i + σF) / |S_i + σF|, where F is a 3D Gaussian random vector.

**Adaptive σ**: recalculated each sweep:
- f = 0.5 / (1 - R), where R = acceptance_rate from previous sweep
- σ_new = σ * f
- Initial σ = 60 (equivalent to random move)
- If σ > 60, reset to 60 (above Tc all moves accepted)

Implemented as a separate `HeisenbergAdaptiveCore<MC>` wrapper or as an alternative
sweep strategy within `MetropolisCore`.

### B.5 FromParams

Reads "L", "beta", "J", "pbc", "proposal_width" (default π/8 for朴素).
For OPSS: "opss" (default false) to enable adaptive sampling.
Uses `build_chain` lattice (same as other models).

---

## Implementation Order

1. **Phase C**: Add `magnetization()`, `snapshot()` to trait + implementations, extend `measure()`
2. **Phase A**: Add `random_cluster_spin()`, `opposite_spin()` to trait, fix SW/Wolff
3. **Phase B**: Add `HeisenbergModel` with朴素 Metropolis, then OPSS adaptive sampling
