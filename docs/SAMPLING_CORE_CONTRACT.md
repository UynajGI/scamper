# Sampling Core Contract — Version 2

This document defines the invariants and contracts for the CMC.rs + Carlo.rs
sampling core. Any code touching the trial evaluation pipeline, the energy
cache, or the snapshot format must respect these rules.

---

## 1. Trial Evaluation Protocol

### 1.1 evaluate_trial MUST NOT mutate formal state

- `TrialEvaluator::evaluate_trial(&self, ...)` takes `&self` (immutable).
- The state (spins, energy, any cached fields) must be bit-identical before and
  after the call.
- Only the `Patch` buffer may be modified.
- **Verifiable property:** calling `evaluate_trial` twice with the same movement
  and a fresh `Patch` each time MUST produce the same `Delta`.

### 1.2 commit_trial can only be called once per trial

- Each `evaluate_trial → commit_trial` pair forms a transactional boundary.
- Calling `commit_trial` without a preceding `evaluate_trial` is undefined
  behavior (the Patch contents are unspecified).
- Calling `commit_trial` twice for the same patch is undefined behavior
  (double-counts the delta energy).

### 1.3 Energy is physical energy, no β factor

- `Hamiltonian::local_energy` and `delta_energy` return bare physical energy
  (not β·E).
- `batch_delta_energy` likewise returns physical energy.
- The `Ensemble` trait applies β exactly once in `log_weight_ratio`.
- `CanonicalEnsemble` computes `-beta * delta.energy + delta.log_jacobian`.

---

## 2. Proposal Ratio Convention

### 2.1 Direction

```
log_reverse_over_forward = ln q(old | new) − ln q(new | old)
```

- Symmetric proposals set `log_reverse_over_forward = 0.0`.
- The Metropolis-Hastings acceptance log-probability is:

```
log_acceptance = log_weight_ratio + log_reverse_over_forward
```

- Acceptance criterion:

```
accepted = log_acceptance >= 0.0 || ln(U(0,1]) < log_acceptance
```

### 2.2 Hastings correction is bounded

- `log_reverse_over_forward` must be finite (NaN is rejected by `assert!` in
  `metropolis_hastings_step`).
- The random draw uses `U.max(f64::MIN_POSITIVE)` to avoid `ln(0)`.

---

## 3. Phase Lifecycle

### 3.1 Phase transitions

```
Initialization → Thermalization → Measurement → Finished
```

- Transitions are monotonic (never go backward).
- The Carlo.rs `Scheduler` owns `Context::enter_phase()` calls.
- `MonteCarlo::on_phase_start(phase)` fires after the context enters the phase.
- `MonteCarlo::on_phase_end(phase)` fires before the context leaves the phase.

### 3.2 Phase semantics

| Phase | Adaptation | Measurements | Sweeps |
|-------|-----------|-------------|--------|
| Initialization | N/A | No | None |
| Thermalization | **Allowed** | Not accumulated | Yes |
| Measurement | **Frozen** | Accumulated | Yes |
| Finished | N/A | No | None |

### 3.3 Two execution paths

- `run_one`: fixed sweep counts (RunConfig), scheduler owns phase transitions.
- `run_controlled`: `AdaptiveRunControl` trait decides when to transition.
  Both paths emit the same hook call sequence.

---

## 4. Cache Invariants

### 4.1 Energy cache

- `System::energy` MUST remain finite at all times.
- After `commit_trial`: `energy += patch.delta_energy`.
- `assert!(energy.is_finite())` guards this in `commit_trial`.
- `recompute_energy()` performs exact recomputation and can repair accumulated
  floating-point drift.
- Energy audit: `MetropolisCore` periodically recomputes exact energy when
  `energy_check_interval > 0`.

### 4.2 Spin cache

- All spin components MUST remain finite at all times.
- Spins are validated in `System::validate()`.
- Proposals must not introduce NaN or Inf components (asserted in
  `evaluate_trial`).

### 4.3 ThermodynamicDelta

- `delta_energy` must be finite (asserted in `evaluate_trial`).
- `log_jacobian` and `volume` must be finite (checked via `is_finite()`).
- Users of `log_weight_ratio` assert `!is_nan` on the result.

### 4.4 BatchEnergyWorkspace

- Generation-stamped to avoid clearing on every use.
- `prepare()` asserts no duplicate sites in a batch move.
- `mark_edge_once()` prevents double-counting physical edges.

---

## 5. Ensemble Independence

### 5.1 Ensemble trait

- Converts `ThermodynamicDelta → ln π(new) − ln π(old)`.
- Stateless: the same delta always produces the same ratio.
- `Send + Sync`: safe for Carlo.rs multi-threaded scheduling.

### 5.2 Current implementations

| Ensemble | Formula |
|----------|---------|
| `CanonicalEnsemble` | `-β·ΔE + Δlog_jacobian` |

- Future: grand canonical would add `μ·ΔN`.

---

## 6. Snapshot Format

### 6.1 Version tag

- Format identifier: `"cmc-rs-snapshot-v2"`.
- Validated on `load_snapshot`; unknown formats are rejected with
  `CarloError::CheckpointCorrupted`.

### 6.2 Bond type encoding

- `BondType` serialized via stable `as_label()` method, NOT `Debug` display.
- Labels use snake_case: `"generic"`, `"chain_x"`, `"square_x"`, `"square_y"`,
  `"square_z"`, `"cubic_x"`, `"cubic_y"`, `"cubic_z"`, `"tri_x"`, `"tri_y"`,
  `"tri_diag"`, `"honey_x"`, `"honey_y"`, `"kagome"`.
- `BondType::from_label()` is the stable inverse; unknown labels are rejected.

### 6.3 JSON schema

```json
{
  "format": "cmc-rs-snapshot-v2",
  "spins": [f64, ...],          // flat site-major: n_sites * spin_dim
  "beta": f64,
  "spin_dim": usize,
  "n_sites": usize,
  "n_edges": usize,
  "offsets": [usize, ...],
  "neighbors": [usize, ...],
  "edge_ids": [usize, ...],
  "edges": [
    { "source": usize, "target": usize, "kind": "generic", "weight": f64 }
  ]
}
```

### 6.4 Energy handling

- The snapshot does NOT store cached energy.
- On load, energy is always recomputed via `recompute_energy()`.
- This guarantees correctness even if the snapshot was produced by a different
  compiler version or optimization level.

---

## 7. Error Types

### 7.1 Checkpoint errors

Use `CarloError::CheckpointCorrupted` for:
- Snapshot format validation failures
- Topology mismatches
- Corrupted field values

### 7.2 Configuration errors

Use `CarloError::InvalidConfig` for:
- Invalid model/lattice parameters
- Parameter parsing failures
- Run configuration errors

### 7.3 Assertion panics

Assertions guard programming errors that should never occur in production:
- NaN energy or log-probability
- Dimension mismatches
- Out-of-range site indices
- Duplicate batch sites

---

## 8. Testing Contract

### 8.1 Statistical correctness

- Small systems (N ≤ 4) must match exact enumeration within 3σ.
- Different algorithms (Metropolis, Wolff, SW) must agree within statistical
  error at identical parameters.
- Fixed seeds must produce bitwise-identical results.

### 8.2 Detailed balance

- For N ≤ 4, the transition matrix T(x→y) must satisfy:

```
π(x) · T(x→y) ≈ π(y) · T(y→x)
```

- Tolerance accounts for Monte Carlo sampling noise (~3%).

### 8.3 Snapshot persistence

- `save_snapshot → load_snapshot` must round-trip spin state and energy.
- Split runs (therm → save → restore → meas) must produce identical final
  state as continuous runs (therm → meas).
- Format corruption must be rejected with `CheckpointCorrupted`.
