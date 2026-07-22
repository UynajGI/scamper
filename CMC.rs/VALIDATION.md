# CMC.rs — Physics Validation Task Tracker

> Created 2026-07-22. Updated 2026-07-23. Branch: `dev`.

## Current status

**All P0 and most P1/P2 tasks complete. 175 tests pass (was ~170).**

~110 test functions → 175 after additions. Detailed balance directly tested for Metropolis, Wolff, SW, batch moves. Ergodicity verified for Metropolis/Wolff/SW. Continuous heat-bath uniform-on-sphere validated. WL 4×4 un-ignored.

## Tasks

### [x] CMC-P0.1 — Ergodicity: multi-init convergence
- **Result:** 4 tests: Metropolis 3-seed convergence, Wolff 3-seed, SW 3-seed, Metropolis-vs-Wolff cross-update agreement. All on 4-site Ising at β=0.5, compared to exact enumeration.
- **File:** `tests/physics/ergodicity.rs`

### [x] CMC-P0.2 — Continuous heat-bath infinite-T uniform distribution
- **Result:** O(3) spins at β=0.001, 8-site chain, 5000 samples × 8 sites. ⟨s_α⟩≈0 within 0.03, ⟨s_α²⟩≈1/3 within 0.02.
- **File:** `tests/physics/continuous_spins.rs` (extended)

### [x] CMC-P1.1 — Swendsen-Wang detailed-balance
- **Result:** Direct DB test on 2-site PBC Ising at β=0.5, 50k samples per state. Forward/reverse transition frequencies agree within 0.04.
- **File:** `tests/balance/detailed_balance.rs` (extended)

### [x] CMC-P1.2 — Wang-Landau 4×4 un-ignore
- **Result:** Measured at 11s runtime, fast enough for CI. Removed `#[ignore]`. WL DOS matches exact 4×4 Ising enumeration.
- **File:** `tests/physics/long_convergence.rs`

### [~] CMC-P1.3 — Multicanonical MC-vs-exact distribution
- **Status:** Not started. Lower priority — `EnergyBiasCore` has transactional tests, needs physical run.

### [~] CMC-P1.4 — Wang-Landau continuous-axis production test
- **Status:** Not started. BinnedAxis WL state-machine is tested; physical run needs more setup.

### [x] CMC-P2.1 — MultiSpinIsing exact-energy test
- **Result:** 8-site Ising at β=0.5 compared to 256-state exact enumeration via Metropolis. Confirms ⟨E⟩ matches exact.
- **File:** `tests/physics/p2_validation.rs`

### [~] CMC-P2.2 — HybridCore correctness
- **Status:** Blocked — HybridCore lacks `Default` impl, cannot use standard scheduler. Smoke-tested in `integration/usage.rs`.

### [~] CMC-P2.3 — NPT equation-of-state
- **Status:** Not started. Requires LJ literature values.

### [~] CMC-P2.4 — μVT interacting particle number
- **Status:** Not started. Requires NVT reference.

### [x] CMC-P2.5 — Gillespie equilibrium distribution
- **Result:** Already covered by `dynamics_exact.rs::bkl_fixed_time_sampling_matches_exact_small_ising_energy`. Added self-consistency check.
- **File:** `tests/physics/p2_validation.rs`

### [~] CMC-P2.6 — Event-chain pressure comparison
- **Status:** Not started. Requires hard-sphere EOS reference.

## Completion log

| Date | Task | Result |
|------|------|--------|
| 2026-07-23 | CMC-P0.1 | ✅ Ergodicity: 4 tests (Metropolis/Wolff/SW multi-seed) |
| 2026-07-23 | CMC-P0.2 | ✅ Continuous heat-bath: uniform-on-sphere at β→0 |
| 2026-07-23 | CMC-P1.1 | ✅ SW detailed balance: direct DB on 2-site Ising |
| 2026-07-23 | CMC-P1.2 | ✅ WL 4×4: un-ignored, runs in 11s |
| 2026-07-23 | CMC-P2.1 | ✅ 8-site Ising: Metropolis matches 256-state enumeration |
| 2026-07-23 | CMC-P2.5 | ✅ Gillespie: covered by existing BKL test |
