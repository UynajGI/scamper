# Maturity Assessment — Physics Validation Coverage

> Updated 2026-08-19. Covers all 4 crates: Carlo.rs, CMC.rs, QMC.rs, MCMC.rs.
> Per-solver assessment against 8 production-readiness criteria.
> All four crates' physics validation complete (tracker 22/22).
> CMC.rs production hardening complete 2026-08-19: named residues closed,
> input validation audited for all 19 solvers (criterion G evidence in
> `CMC.rs/tests/physics/input_validation.rs`, per-solver domains in
> `CMC.rs/VALIDATION.md`).

## Maturity levels

| Level | Meaning |
|-------|---------|
| **Experimental** | Algorithm implemented, minimal validation |
| **Research-grade** | Validated vs ED/analytic on small systems |
| **Production-ready** | Physics + statistics + software sustained guarantees |
| **Stable** | Production + backward-compat commitment |

## Criteria

| # | Criterion | Standard |
|---|-----------|----------|
| A | Deterministic identities (machine precision) | `abs_error < 1e-12` |
| B | ED/exact comparison (≥3 observables) | E, ⟨M²⟩, C(τ) etc. |
| C | Analytic solvable limits | g=0, high-T, strong field |
| D | Per-update balance & invariants | Each update independently |
| E | Ergodicity evidence | Multi-init convergence, sector access |
| F | Cross-solver validation | Two independent solvers agree |
| G | Input validation | Invalid → error, never silent |
| H | Documented limitations | User-facing validated domain |

---

## 1. Summary status table

### Carlo.rs (framework)

| Component | Status | Notes |
|-----------|--------|-------|
| MC loop / Scheduler / Context | **stable** | 302 suite tests (297 + 5 MPI-ignored), lifecycle, reproducibility |
| Measurements / Accumulators | **stable** | Binning, covariance, jackknife, real autocorr_time |
| Checkpoint (in-memory) | **stable** | JSON serde, clock round-trip |
| Checkpoint (HDF5) | **production-ready** | 7 tests: sweep/clocks/measurements/RNG round-trip, legacy back-compat + loud `CheckpointLoadReport` fallback |
| Error analysis / merge | **research-grade** | Regular + decorrelated autocorr tested; no error-bar coverage calibration |
| Parallel tempering (MPI) | **production-ready** | Exchange acceptance vs analytic (5σ, 50k attempts), np 4 replica round-trip, clean mismatch errors |
| MPI backend | **production-ready** | np 1/2/4 exact fan-out (multiset + exact sum), RNG stream pinning, single-init exclusivity |
| RNG stream derivation | **stable** | Domain separation, reproducibility, thread-count independence |

### CMC.rs (19 solvers)

Production status since 2026-08-19: every solver has criterion-G evidence
(input-validation audit; three source-side holes found and fixed) and an
up-to-date validated domain (criterion H, `CMC.rs/VALIDATION.md`). 15 of 19
are **production-ready**; 4 remain research-grade for physics/evidence
reasons documented below.

| Solver | Status | Key gap / evidence beyond research-grade |
|--------|--------|---------|
| Local Metropolis | **production-ready** | z-score + connectivity + direct DB; continuous-spin and frustrated-triangle cross-solver; Potts q=3/4 cross-solver + analytic limits; input validation audited |
| Wolff cluster | **production-ready** | z-score; direct DB N=3; Potts q=3/4 exact enumeration (E, m², C) |
| Swendsen-Wang | **production-ready** | Direct DB on 2-site Ising; Potts q=3/4 exact enumeration (E, m², C); β=0 uniform limit; 4-solver cross-solver at q=3 βc |
| Heat bath (discrete) | **production-ready** | Exact energies N=4/8 z-scores + direct DB; Potts q=3/4 exact enumeration |
| Heat bath (continuous O(N)) | **production-ready** | Finite-T conditional vs Bessel/Langevin (O(2)+O(3)); infinite-T uniform; DB by construction (exact conditional sampling) |
| Microcanonical over-relaxation | **research-grade** | Energy conservation tested; not ergodic alone (physical, documented) — use composed with ergodic kernels |
| Hybrid (composed) | **research-grade** | Hybrid(Metropolis, Wolff) ⟨E⟩/⟨m²⟩ vs exact enumeration; other compositions and hybrid-specific ergodicity/cross-solver evidence untested (documented) |
| MultiSpinIsing (64-replica) | **research-grade** | 8-site exact enumeration cross-check; no cross-solver vs scalar Metropolis, no multi-seed z coverage (documented) |
| Wang-Landau DOS | **production-ready** | 2/4-site exact + 4×4 (11s, CI-ready) + BinnedAxis production run; unattainable `minimum_visited_fraction` now terminates loudly (`UnreachableBins`) instead of silently non-converging |
| Multicanonical / umbrella | **production-ready** | Full MC-vs-exact distribution (moments + P(E) + flatness); transactional bias algebra |
| Worm (Ising HT graph) | **production-ready** | Exact energy + endpoint correlation; cross-solver vs spin Metropolis; multi-component lattices rejected loudly (silent-freeze failure mode closed); multi-defect worms documented as not implemented |
| Kawasaki dynamics | **production-ready** | Sector-restricted exact equilibrium (BFS-reachable sector); β/J input validation fixed (was a panic path) |
| BKL / n-fold-way | **production-ready** | Residence-time-weighted equilibrium vs exact (N=4/8, 3 β); input validation fixed (β/J/event window) |
| Gillespie | **production-ready** | 3-state CTMC exact stationary π + Ising equilibrium |
| Event chain (hard sphere) | **production-ready** | EOS via exact virial B₂/B₃ (contact-value + Richardson); cross-solver vs Metropolis NVT |
| Particle NVT | **production-ready** | Energy distribution vs quadrature; input validation audited |
| Particle NPT | **production-ready** | Finite-N ideal gas ⟨V⟩ = (N+1)kT/P exact (long test) + directional response |
| Particle μVT | **production-ready** | Ideal gas Poisson ⟨N⟩ exact (long test) + directional response |
| Rigid molecule | **production-ready** | Equilibrium vs quadrature references; one-body dipolar external field (`DipolarExternalField`) validated vs 2D von Mises + 3D Langevin free-rotor answers with a machine-precision −E·μ identity |

### QMC.rs (4 solvers)

| Solver | Status | Key gap |
|--------|--------|---------|
| Lattice directed-loop | **research-grade** | E+⟨m²⟩+NN Sz+χ_z+Binder U4 vs ED; analytic limits; ergodicity |
| Wormhole (spin-boson) | **research-grade** | Interacting MC-vs-ED (3 obs, z-scores < 4); cross-solver ×2; ergodicity |
| Occupation (cavity-QED) | **research-grade** | ⟨σz⟩,⟨σx⟩,E,⟨n⟩ vs ED (Rabi+JC); cross-solver vs wormhole |
| Cluster (longitudinal SB) | **research-grade** | ⟨σz⟩,⟨σx⟩,C(τ) vs ED; cross-solver vs wormhole; single-mode only |

### MCMC.rs (6 kernels + 3 combinators)

| Solver | Status | Key gap |
|--------|--------|---------|
| NUTS | **research-grade** | Accept-step DB not directly tested (leapfrog reversibility + 6-solver posterior agreement instead) |
| Static HMC | **research-grade** | Same as NUTS |
| Random-walk Metropolis | **research-grade** | Machine-precision DB rule + binned flow balance + cross-solver + non-Gaussian |
| Component-wise Metropolis | **research-grade** | Flow balance + cross-solver + non-Gaussian recovery |
| Slice | **research-grade** | Moment recovery (Gaussian/t/bimodal) + cross-solver agreement |
| Gibbs | **research-grade** | Exact conditional + atomicity |
| Combinators (Then/Repeat/Mixture) | **research-grade** | Determinism at boundaries |

---

## 2. Per-solver 8-criteria breakdown

### CMC.rs — Local Metropolis

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | Deterministic identities | ✅ PASS | `metropolis_cache_matches_exact_energy`; energy-cache audit |
| B | ED/exact (≥3 obs) | ✅ PASS | Onsager E, exact finite-Ising N=2–4, 2×2×2 within 4σ |
| C | Analytic limits | ✅ PASS | g=0, high-T ⟨m²⟩→1/4N, strong field polarizes, dimer exact |
| D | Per-update balance | ✅ PASS | `asymmetric_hastings_detailed_balance_n2` (direct DB) |
| E | Ergodicity | ✅ PASS | Multi-seed convergence + BFS strong connectivity (N=2) + aperiodicity |
| F | Cross-solver | ✅ PASS | Metropolis vs Wolff vs SW on 8×8 (3σ); Potts q=3 4-solver at βc (2026-08-19) |
| G | Input validation | ✅ PASS | Topology mismatch, malformed params, unknown lattice; `input_validation.rs` audit (2026-08-19) |
| H | Documented limits | ✅ PASS | VALIDATION.md: 19 solvers with validated domain |

### CMC.rs — Wolff cluster

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | `wolff_batch_cache_matches_exact_energy` |
| B | ✅ PASS | Langevin O(3), Bessel-ratio O(2) exact; 8×8 cross-solver; Potts q=3/4 full enumeration E+m²+C (2026-08-19) |
| C | ✅ PASS | XY/Heisenberg cluster activation, low-T M→1; Potts β=0 uniform, β=8 frozen |
| D | ✅ PASS | `wolff_detailed_balance_n3` (direct DB) |
| E | ✅ PASS | Multi-seed convergence test |
| F | ✅ PASS | Cross-solver 8×8; Potts q=3 4-solver at βc |
| G | ✅ PASS | Sign-problem rejection, model manifold; input-validation audit |
| H | ✅ PASS | VALIDATION.md documents validated domain |

### CMC.rs — Swendsen-Wang

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Batch delta cache |
| B | ✅ PASS | Potts q=3/4 full q^N enumeration: ⟨E⟩, ⟨m²⟩, C at β=0.3/0.8 on N=4/N=8 (2026-08-19; q=2 mapped onto Ising by the enumerator anchor) |
| C | ✅ PASS | β=0 exactly-uniform states (⟨E⟩ = −JΣw/q closed form); strong-coupling frozen ground state |
| D | ✅ PASS | Direct DB on 2-site Ising |
| E | ✅ PASS | Multi-seed convergence test; strongly connected N=2 |
| F | ✅ PASS | 8×8 cross-solver; Potts q=3 4-solver at βc = ln(1+√3) |
| G | ✅ PASS | Input-validation audit |
| H | ✅ PASS | VALIDATION.md documents validated domain |

### QMC.rs — Lattice directed-loop

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ⚠️ PARTIAL | Scattering table balance (S=1, 1e-12); generic S missing |
| B | ✅ PASS | 3-site Heisenberg: E, ⟨Sz_iSz_j⟩, ⟨m²⟩, χ_z, Binder U4 vs ED |
| C | ✅ PASS | Zero-coupling, high-T, strong-field, Ising dimer, dimer correlation |
| D | ✅ PASS | Scattering table detailed balance (both policies) |
| E | ✅ PASS | 4-site Heisenberg from 3 initial states converge |
| F | N/A | Single lattice solver |
| G | ✅ PASS | Sign-problem rejection, frustration detection |
| H | ✅ PASS | README documents S>1/2 caveat, validated domain |

### QMC.rs — Wormhole (spin-boson)

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Scattering DB both policies (1e-12); vertex weights |
| B | ✅ PASS | Interacting MC-vs-ED: ⟨σx⟩, expansion order, C(β/2); z-scores < 4 |
| C | ✅ PASS | Poisson expansion, 2-state partition, spin-inversion symmetry |
| D | ✅ PASS | Table-level DB; loop-abort rollback; worldline invariants |
| E | ✅ PASS | 4-seed convergence + z-score framework |
| F | ✅ PASS | wormhole↔occupation (2 tests), wormhole↔cluster (1 test) |
| G | ✅ PASS | Channel validation, non-stoquastic rejection |
| H | ✅ PASS | README documents basis rotation, convention differences |

### QMC.rs — Occupation (cavity-QED)

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Transfer matrix eigensystem; sqrt(n) matrix elements |
| B | ✅ PASS | ⟨σz⟩,⟨σx⟩,E,⟨n⟩ vs ED (Rabi+JC); long crossover #[ignore] |
| C | ✅ PASS | Free spin tanh, uncoupled Bose distribution, exact Z |
| D | ⚠️ PARTIAL | Non-stoquastic rejection; no per-update DB (transfer matrix) |
| E | ⚠️ PARTIAL | `sampler_only_visits_states_within_basis` (bounds, not connectivity) |
| F | ✅ PASS | Cross-solver vs wormhole (free two-level system) |
| G | ✅ PASS | Beta/slices/cutoff validation |
| H | ✅ PASS | README + Rabi QPT tests document domain |

### QMC.rs — Cluster (longitudinal SB)

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Retarded kernel KMS symmetry; quadrature normalization |
| B | ✅ PASS | ⟨σz⟩,⟨σx⟩ (via kinks), C(β/2) vs exact cosh; interacting ED |
| C | ✅ PASS | Free 2-level: kink count vs exact mean; ⟨σz⟩ vs tanh |
| D | ✅ PASS | Worldline invariants (even kink); 10k updates validated |
| E | ❌ MISSING | |
| F | ✅ PASS | Cross-solver vs wormhole (longitudinal model) |
| G | ✅ PASS | |
| H | ⚠️ PARTIAL | Single-mode only; multi-mode interacting untested |

### MCMC.rs — NUTS

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Leapfrog reversibility (1e-12); metric velocity exact |
| B | ✅ PASS | Correlated Gaussian (5 moments) + Student-t cov + bimodal moments vs analytic |
| C | ✅ PASS | Non-Gaussian recovery (bimodal mixture occupancy + Student-t ν=5) |
| D | ⚠️ PARTIAL | Leapfrog reversibility; accept-step DB not directly tested for HMC family |
| E | ⚠️ PARTIAL | Step-size search rescues bad scale; no sector test |
| F | ✅ PASS | 6-solver posterior agreement (15 pairs × 5 moments, \|z\| < 4) |
| G | ✅ PASS | NaN, dim mismatch, non-PD mass matrix, bad configs |
| H | ✅ PASS | VALIDATION.md documents validated targets (Gaussian, Student-t, bimodal) |

### MCMC.rs — Static HMC

Same as NUTS minus U-turn tests.

### MCMC.rs — Random-walk Metropolis

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Acceptance statistic bit-for-bit `exp(min(0, Δlog p))`; scale adaptation freeze |
| B | ✅ PASS | Correlated Gaussian (5 moments) + Student-t cov + bimodal moments vs analytic |
| C | ✅ PASS | Bimodal mixture (occupancy symmetry) + Student-t ν=5 |
| D | ✅ PASS | Machine-precision DB identity + binned empirical flow balance |
| E | ✅ PASS | Alternating mode inits on bimodal target (16 seeds); multi-seed convergence |
| F | ✅ PASS | 6-solver posterior agreement (\|z\| < 4) |
| G | ✅ PASS | Zero-dim, bad covariance rejected |
| H | ✅ PASS | VALIDATION.md documents validated domain |

### Carlo.rs — Checkpoint (HDF5)

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| A | ✅ PASS | Sweep count, clocks round-trip |
| B | ✅ PASS | Measurements preserved across checkpoint |
| C | ⚠️ PARTIAL | Legacy fallback (unwrap_or(0)) tested |
| D | ✅ PASS | Clocks (attempted/accepted/event_time) now written+read |
| E | N/A | |
| F | N/A | |
| G | ⚠️ PARTIAL | Missing datasets fall back gracefully |
| H | ✅ PASS | VALIDATION.md documents test coverage |

---

## 3. Cross-cutting gaps (universal blockers)

### Resolved (2026-07-23)

1. ~~Ergodicity entirely untested.~~ → CMC: 4 multi-seed tests. QMC: lattice 3-init + impurity 4-seed. Carlo: thread-count independence.
2. ~~No multi-seed statistical monitoring.~~ → CMC: 16-seed z-score framework. QMC: 3–4 seed z-score per test.
3. ~~QMC cross-solver validation absent.~~ → wormhole↔occupation (2 tests), wormhole↔cluster (1 test).
4. ~~QMC susceptibility/Binder unvalidated.~~ → χ_z vs ED, Binder U4 vs ED.
5. ~~Carlo HDF5 checkpoint untested.~~ → 5 tests (clocks fixed).
6. ~~Carlo autocorr_time hardcoded.~~ → Real estimation + 6 AR(1) tests.
7. ~~Carlo MPI untested.~~ → PT exchange verified under mpirun -np 2.
8. ~~strict-repro dead feature.~~ → Removed.

### Still open

9. ~~MCMC: No detailed-balance test~~ → 2026-08-14: machine-precision log-Metropolis rule + binned empirical flow balance (RW + ComponentWise).
10. ~~MCMC: No AR(1) reference test for ESS~~ → 2026-08-14: AR(1) ρ ∈ {0, 0.5, 0.9, 0.99} vs closed-form ESS.
11. ~~MCMC: No cross-solver moment agreement~~ → 2026-08-14: 6 solvers × 15 pairs × 5 moments, |z| < 4.
12. ~~MCMC: No non-Gaussian recovery~~ → 2026-08-14: bimodal mixture + Student-t ν=5, 4 solvers.
13. ~~CMC NPT/μVT equilibrium values wrong~~ → Resolved by e1a07e4: finite-N reference (⟨V⟩=(N+1)kT/P, Poisson ⟨N⟩) matches exactly; earlier mismatch was a test-formula error.
14. ~~Nightly z-score infrastructure~~ → 2026-08-14: `zscore-monitor` job @64 seeds via `SCUTTLE_ZSCORE_SEEDS`; `just nightly-zscore` locally.
15. ~~Carlo PT exchange dynamics untested; MPI only at -np 2; HDF5 legacy fallback silent~~ → 2026-08-14: analytic acceptance-rate validation (5σ), np 1/2/4 backend exactness + replica round-trips, loud `CheckpointLoadReport`; nightly `carlo-framework` job (`just mpi-test`).

### Still open (2026-08-14, Carlo production hardening residue)

16. Carlo error analysis / merge remains research-grade: no error-bar coverage-calibration test (fraction of independent replicas whose truth falls within 1σ ≈ 68%).
17. `Measurements::read_checkpoint_hdf5` still contains 4 silent fallbacks (shape→zeros, flag→false, member_names→empty ×2) — potential silent data loss on corrupt files; report-only for now.
18. MPI failure modes: a rank-local test panic under mpirun deadlocks peers in collectives (finalize semantics); `with_ranks_per_run(k>1)` topology is untestable in-process (single-init exclusivity). The `pt_exchange` end-to-end entry-point test is np 2 only (its entry point owns the MPI init and cannot probe world size).
19. ~~CMC named gaps residue~~ → Closed 2026-08-19:
    - **Potts q>2**: heat bath + SW + Wolff validated against full q^N enumeration (⟨E⟩, ⟨m²⟩, C; N=4/N=8; q=3/4; β=0.3/0.8), plus β=0/strong-coupling analytic limits, 4-solver cross-solver at q=3 βc, q=4 βc = ln 3 anchor (`CMC.rs/tests/physics/potts_exact.rs`).
    - **Molecule one-body external field**: `DipolarExternalField` API added (additive, backward-compatible; neutral molecules enforced loudly); validated against 2D von Mises and 3D Langevin free-rotor moments with a machine-precision −E·μ identity (`CMC.rs/tests/physics/molecule_external_field.rs`).
    - **Wang-Landau `minimum_visited_fraction`**: unattainable fractions now terminate loudly (`WangLandauTermination::UnreachableBins`) once the discovery plateau establishes the reachable set, instead of silently running to the sweep guard; checkpoints carry the evidence (`CMC.rs/tests/physics/wang_landau_convergence.rs`).
    - **Multi-component worm**: honest-limitation route — multi-component (disconnected/isolated-site) lattices, whose components beyond the defect pair would silently freeze, are rejected loudly at construction; multi-defect/multi-leg worms are documented as not implemented; worm-vs-Metropolis cross-solver added.
    - **Input validation (criterion G) for all 19 solvers**: audited with rejection tests (`CMC.rs/tests/physics/input_validation.rs`); three source holes fixed (Kawasaki/BKL β+J panic path, MultiSpinIsing lattice validation, BKL event-window validation).

### Still open (2026-08-19, CMC production residue)

20. CMC solvers remaining research-grade and why: **microcanonical over-relaxation** (not ergodic alone — physical; use in composition), **hybrid compositions** (only Metropolis+Wolff validated), **MultiSpinIsing** (no cross-solver vs scalar Metropolis, no multi-seed z coverage). Their criterion-G input validation and criterion-H domains are covered.

---

## 4. Prioritized roadmap

### P0 — Correctness risks (must fix)

| # | Task | Crate | Why |
|---|------|-------|-----|
| P0.1 | Wormhole interacting MC-vs-ED (≥3 obs) | QMC | Weakest solver; only free-limit validated |
| P0.2 | Lattice QMC analytic limits (g=0, high-T, strong field) | QMC | Basic limits missing |
| P0.3 | Lattice QMC susceptibility χ_z vs ED | QMC | χ measured but unvalidated |
| P0.4 | MCMC detailed-balance test for accept/reject | MCMC | Core MCMC correctness unverified |
| P0.5 | MCMC AR(1) reference test for ESS | MCMC | ESS estimator uncalibrated |
| P0.6 | Carlo.rs HDF5 checkpoint round-trip | Carlo | Production restart path untested |

### P1 — Production blockers

| # | Task | Crate |
|---|------|-------|
| P1.1 | Ergodicity: multi-init convergence test (all solvers) | ALL |
| P1.2 | QMC cross-solver: wormhole↔occupation numerical agreement | QMC |
| P1.3 | QMC cross-solver: wormhole↔cluster numerical agreement | QMC |
| P1.4 | SW detailed-balance test (direct, not just cross-solver) | CMC |
| P1.5 | MCMC cross-solver: NUTS vs RWM vs Slice posterior agreement | MCMC |
| P1.6 | Carlo.rs MPI PT exchange protocol test | Carlo |
| P1.7 | CMC Wang-Landau 4×4 un-#[ignore] | CMC |
| P1.8 | MCMC non-Gaussian recovery (mixture, Student-t) | MCMC |

### P2 — Robustness

| # | Task | Crate |
|---|------|-------|
| P2.1 | Acceptance formula machine-precision validation | CMC, QMC |
| P2.2 | Continuous heat-bath analytic distribution match | CMC |
| P2.3 | NPT/μVT interacting equation-of-state | CMC |
| P2.4 | MultiSpinIsing exact-energy test | CMC |
| P2.5 | Thread-count independence test | ALL |
| P2.6 | Carlo.rs Estimate.autocorr_time real estimation | Carlo |
| P2.7 | strict-repro feature test | Carlo |
| P2.8 | Multi-seed nightly z-score monitoring infrastructure | ALL |

---

## 5. Status statement

**Repository: research-grade overall (QMC/MCMC physics complete); CMC.rs and the Carlo.rs framework core reached their target maturity on 2026-08-19.**

- **Carlo.rs**: stable framework core; HDF5 checkpoint, MPI backend, and PT exchange now **production-ready** (analytic exchange-acceptance validation, np 1/2/4 exact fan-out, loud legacy fallback; nightly `carlo-framework` regression job) — 302 suite tests (297 + 5 MPI-ignored); error analysis/merge still research-grade
- **CMC.rs**: 273 suite tests (258 + 15 long) + 69 lib. **Production hardening complete (2026-08-19)**: 15 of 19 solvers **production-ready** (physics A–F + audited input validation G + documented domains H); Potts q=3/4 validated against full enumeration for heat bath/SW/Wolff with cross-solver and analytic limits; molecule one-body dipolar external field validated against the Langevin/von Mises free-rotor answers; Wang-Landau unattainable-visited-fraction runs fail loudly (`UnreachableBins`); the worm rejects multi-component lattices loudly and carries a cross-solver check. Remaining research-grade with documented reasons: microcanonical over-relaxation (not ergodic alone), hybrid compositions (one combination validated), MultiSpinIsing (no cross-solver/multi-seed coverage)
- **QMC.rs**: 178 tests (171 + 7 long). All 4 solvers research-grade with ED cross-checks, cross-solver validation, ergodicity, and z-score framework
- **MCMC.rs**: 72 tests (69 + 3 long). All 6 kernels research-grade: detailed balance (machine-precision + statistical), ESS calibrated on AR(1), 6-solver posterior agreement, non-Gaussian recovery; nightly z-score monitoring covers CMC/QMC at 64 seeds
