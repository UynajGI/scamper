# QMC.rs — Physics Validation Task Tracker

> Created 2026-07-22. Branch: `dev`.
> Baseline: `../MATURITY_ASSESSMENT.md`.

## Current status

**Occupation & Cluster: research-grade with genuine MC-vs-ED.**
**Wormhole: research-grade but only free-limit validated.**
**Lattice: research-grade for S=1/2 but missing analytic limits.**

~150 tests across 4 solvers. Scattering table detailed balance verified at machine precision (both policies). Worldline invariants checked each sweep. Cross-solver validation completely absent.

## Tasks

### [ ] QMC-P0.1 — Wormhole interacting MC-vs-ED (≥3 observables)
- **Problem:** Wormhole's interacting regime is completely unvalidated. Only free limits (Poisson expansion, 2-state tanh, spin-inversion symmetry) are tested. An error in vertex weights or retarded interaction could produce correct free limits but wrong interacting results.
- **Plan:** Build single-mode Rabi: Δ=1, Ω=1, g=0.3, β=20. ED via nalgebra (dim 2×N_b=40). Run WormholeEngine: 30k thermal + 100k measure. Compare: (1) ⟨σz⟩_sampled (=physical ⟨σx⟩), (2) expansion order (=−β⟨E⟩), (3) C(β/2). Use 10 seeds, |z|<4.
- **File:** `QMC.rs/tests/impurity/wormhole_interacting_ed.rs` (new)
- **Status:** not started

### [ ] QMC-P0.2 — Lattice QMC analytic limits
- **Problem:** Lattice solver has ED cross-checks but no analytic-limit tests. Zero-coupling, high-T, strong-field are basic sanity checks that could catch gross bugs.
- **Plan:** 3 tests: (1) J=0 → E=0, ⟨m²⟩=S(S+1)/3, zero vertices; (2) β=0.1 → ⟨m²⟩≈infinite-T value; (3) h_z≫J → all spins align, ⟨m_z⟩→S.
- **File:** `QMC.rs/tests/lattice/lattice_limits.rs` (new)
- **Status:** not started

### [ ] QMC-P0.3 — Lattice QMC susceptibility χ_z vs ED
- **Crate:** QMC.rs
- **Problem:** Susceptibility χ_z = β(⟨m²⟩−⟨m⟩²) measured but never compared to ED. The connected subtraction could have a sign error or normalization bug.
- **Plan:** Extend the existing 3-site Heisenberg ED test to also compare χ_z. Compute ⟨m²⟩_connected from MC (⟨m²⟩−⟨m⟩²), from ED (Tr(m²e^{-βH})/Z − (Tr(me^{-βH})/Z)²), compare.
- **File:** `QMC.rs/tests/lattice/lattice_ed.rs` (extend)
- **Status:** not started

### [ ] QMC-P1.1 — Cross-solver: wormhole↔occupation (free limit)
- **Problem:** Two solvers handle overlapping physics domains but are never compared. Convention differences (basis rotation, bath representation) documented but not numerically verified.
- **Plan:** Free two-level system (g→0): both solvers reduce to isolated spin. Compare ⟨σz⟩, ⟨E⟩, partition function. Then weakly interacting (small g): compare ⟨σz⟩ within convention-adjusted tolerance.
- **File:** `QMC.rs/tests/impurity/cross_solver.rs` (extend)
- **Status:** not started

### [ ] QMC-P1.2 — Cross-solver: wormhole↔cluster (longitudinal coupling)
- **Problem:** Both solvers handle longitudinal spin-boson. Never compared numerically.
- **Plan:** Single-mode bath, Δ=0 (no transverse). Both solvers sign-free in same basis. Compare ⟨σz⟩, kink count / expansion order, C(β/2).
- **File:** `QMC.rs/tests/impurity/cross_solver.rs` (extend)
- **Status:** not started

### [ ] QMC-P1.3 — Lattice QMC ergodicity (multi-init)
- **Problem:** No test for convergence from different initial states.
- **Plan:** 4-site Heisenberg chain. Start from ferromagnetic, Néel, random states. Run 20k sweeps. Compare ⟨E⟩ and ⟨m²⟩ across inits.
- **File:** `QMC.rs/tests/lattice/lattice_limits.rs` (extend)
- **Status:** not started

### [ ] QMC-P1.4 — Impurity ergodicity (multi-init)
- **Problem:** All impurity solvers start from empty worldline. No test for sector accessibility.
- **Plan:** Wormhole: start from empty, saturated (β×many vertices), and random configurations. Run sweeps. Compare observables.
- **File:** `QMC.rs/tests/impurity/ergodicity.rs` (new)
- **Status:** not started

### [ ] QMC-P1.5 — Cluster solver multi-mode interacting ED
- **Problem:** Cluster solver validated against interacting ED only for single-mode. Multi-mode/power-law bath interacting case is untested.
- **Plan:** 2-mode bath, finite coupling. ED (larger matrix). Compare ⟨σz⟩, ⟨σx⟩, C(β/2).
- **File:** `QMC.rs/tests/impurity_cluster_test.rs` (extend)
- **Status:** not started

### [ ] QMC-P2.1 — Binder cumulant M⁴ vs ED
- **Problem:** ⟨M⁴⟩ measured everywhere, validated nowhere. U4 could have normalization errors.
- **Plan:** 3-site Heisenberg. Compare MC U4 to ED U4 (requires computing ⟨m⁴⟩ from density matrix).
- **File:** `QMC.rs/tests/lattice/lattice_ed.rs` (extend)
- **Status:** not started

### [ ] QMC-P2.2 — Full imaginary-time C(τ) profile vs ED
- **Problem:** Only C(β/2) tested for cluster. Full C(τ) at multiple τ points untested.
- **Plan:** Compare MC C(τ) at τ = β/4, β/2, 3β/4 to ED values for 3-site Heisenberg.
- **File:** `QMC.rs/tests/lattice/lattice_ed.rs` (extend)
- **Status:** not started

### [ ] QMC-P2.3 — Lattice S>1/2 ED validation
- **Problem:** S>1/2 only has worldline invariant tests (bounce fallback). No ED comparison.
- **Plan:** S=1 chain, 3 sites. ED with 3²=9 states per site. Compare ⟨E⟩ and ⟨m²⟩. Document caveat if results diverge (expected: bounce fallback is known broken).
- **File:** `QMC.rs/tests/lattice/lattice_ed.rs` (extend)
- **Status:** not started

### [ ] QMC-P2.4 — Thread-count independence
- **Problem:** No test for multi-thread statistical equivalence.
- **Plan:** 8-task QMC run with 1 vs 4 threads via Carlo.rs RayonBackend.
- **File:** `QMC.rs/tests/integration/thread_count.rs` (new)
- **Status:** not started

## Completion log

| Date | Task | Result |
|------|------|--------|
