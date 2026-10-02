# CMC.rs — Physics Validation & Validated Domain

> Updated 2026-10-02. Validation snapshot; current implementation status is
> tracked separately in [PERCOLATION.md](PERCOLATION.md).

## Test suite summary

| Layer | Tests | Runtime |
|-------|-------|---------|
| Default (`cargo test`) | 315 | ~75s |
| Long stochastic (`--ignored`) | 17 | ~40s |
| **Suite total** | **332** | (+104 lib unit tests) |

## Per-solver validated domain

### Local Metropolis (`MetropolisCore`)
- **Validated:** Ising S=1/2 on chain/square/hypercubic graphs, PBC and OBC; q-state Potts (q=3,4) via cross-solver and β=0/strong-coupling limits
- **Observables:** ⟨E⟩, ⟨m²⟩ vs exact enumeration (N=2,3,4,8)
- **Detailed balance:** Direct DB verified on N=2 (1e-15 precision)
- **Ergodicity:** 16-seed z-score framework (|z|<4, no systematic bias)
- **Connectivity:** Explicit Markov chain strongly connected on N=2
- **Cross-solver:** Agrees with Wolff and SW within 3σ on 8×8 Ising; continuous-spin cross-solver XY (O(2)) and Heisenberg (O(3)) vs Wolff on 8×8 (pooled |z|<4 at β=0.7/0.9 and 0.5/0.9); 8×8 q=3 Potts 4-solver agreement (see Potts section)
- **Frustrated:** AFM XY triangle (J<0): full equilibrium vs spectral quadrature reference — ⟨E⟩ and chirality ⟨κ²⟩ at β=1,3 (per-seed |z|<3, ground-state algebra exact)
- **NOT validated:** —

### Wolff cluster (`WolffCore`)
- **Validated:** Ising S=1/2, XY (O(2)), Heisenberg (O(3)) on chain/square; q-state Potts q=3,4 (2026-08-19)
- **Observables:** ⟨E⟩ vs exact enumeration, Langevin/Bessel-ratio analytic limits; Potts ⟨E⟩, ⟨m²⟩, C vs full q^N enumeration (N=4, 8; β=0.3/0.8; per-seed |z|<4, pooled Σz gate)
- **Detailed balance:** Direct DB verified on N=3
- **Ergodicity:** 16-seed z-score, strongly connected on N=2
- **NOT validated:** —

### Swendsen-Wang (`SWCore`)
- **Validated:** Ising S=1/2 on chain/square; q-state Potts q=3,4 (2026-08-19); continuous O(2)/O(3) cluster updates (2026-08-19, see below)
- **Detailed balance:** Direct DB verified on N=2 (50k samples/state)
- **Exact equilibrium:** Potts ⟨E⟩, ⟨m²⟩, C vs full q^N enumeration on N=4/N=8 at β=0.3/0.8 (16-seed z-scores); β=0 exactly-uniform state distribution (⟨E⟩ = −JΣw/q closed form)
- **Continuous spins (2026-08-19):** XY 4-ring vs exact spectral quadrature — ⟨E⟩, ⟨m²⟩, ⟨cos Δθ⟩ at β ∈ {0.6, 1.2}, 8 seeds, |z| < 4 per seed, |z̄| < 2; O(3) 4-ring analytic limits (β→0: ⟨E⟩ = 0 and ⟨m²⟩ = 1/N exactly; β=8 spin-wave ⟨E⟩ = −J·N + (N−1)/β); cross-solver vs Wolff on 8×8 O(2)/O(3) at β=0.9 (pooled |z| < 4 on ⟨E⟩ and ⟨m²⟩) — `sw_continuous.rs`
- **Ergodicity:** Strongly connected on N=2
- **Cross-solver:** Agrees with Metropolis within 3σ on 8×8; 8×8 q=3 near βc = ln(1+√3) 4-solver agreement
- **NOT validated:** —

### Heat bath — discrete (`HeatBathCore`)
- **Validated:** Ising S=1/2; q-state Potts q=3,4 (2026-08-19)
- **Detailed balance:** Direct DB verified on N=2
- **Exact energies:** ⟨E⟩ and ⟨m²⟩ vs exact enumeration on N=4 and N=8 at β=0.3/0.8 (16-seed z-scores, max|z|<2.5); Potts adds ⟨E⟩, ⟨m²⟩, C on N=4/N=8 (2026-08-19)
- **NOT validated:** —

### Heat bath — continuous O(N) (`ContinuousHeatBathCore`)
- **Validated:** O(3) uniform-on-sphere at β→0 (⟨s_α⟩≈0, ⟨s_α²⟩≈1/3)
- **Finite-T conditional distribution:** single-spin ⟨cosθ⟩/⟨cos²θ⟩ vs exact O(2) Bessel ratio I₁(x)/I₀(x) and O(3) Langevin coth(x)−1/x at x ∈ [0.2, 5]; analytic anchors (limits at small/large x, Bessel recurrence) to 1e-4–1e-12; two-spin XY pair equilibrium matches the conditional moment (rotational invariance)
- **NOT validated:** — (finite-T distribution and XY validated 2026-08-14; lattice-scale O(N) cross-solver via Metropolis-vs-Wolff above)

### Microcanonical over-relaxation (`MicrocanonicalCore`)
- **Validated:** Energy conservation to 2e-10, unit norm preservation to 3e-12 (long sweeps)
- **Reflection-map identities (2026-08-19, machine precision):** R(s) = 2(s·ĥ)ĥ − s is an exact involution (R∘R = id, 1e-15), norm-preserving (1e-15), preserves the local-field projection (bond energies unchanged, 1e-15), and is an isometry with |Jacobian| = 1 — O(2): θ' = 2φ − θ exactly (residual < 1e-14); O(3): the map matrix 2ħħᵀ − I is orthogonal with det = +1 (π rotation about the field axis; the field direction is fixed, the orthogonal plane reversed). A deterministic involutive isometry satisfies detailed balance against any rotation-invariant measure: T(s→s') = T(s'→s) = 1
- **Kernel ≡ the reflection map (2026-08-19):** with the sequential visit schedule the kernel sweep is bit-identical to the manual per-site reflection, and every transactional update reports |ΔE| < 1e-12 (per-update machine-precision conservation; `energy_error` after full sweeps ≤ 1e-12)
- **Composition equilibrium (2026-08-19):** Hybrid(Metropolis, Microcanonical) on the XY 4-ring vs exact spectral quadrature (zero mode factored, grid-doubling anchored) — ⟨E⟩, ⟨m²⟩, ⟨cos(θ1−θ3)⟩ at β ∈ {0.6, 1.2}, from hot and cold initializations, 8 seeds, |z| < 4 per seed, |z̄| < 2 per observable
- **Cross-solver (2026-08-19):** the composition vs Wolff on 8×8 O(2) and O(3) at β=0.9 (pooled |z| < 4 on ⟨E⟩ and ⟨m²⟩)
- **Analytic limits in composition (2026-08-19):** β→0 gives exactly ⟨E⟩ = 0 and ⟨m²⟩ = 1/N; β=8 approaches the spin-wave result ⟨E⟩ = −J·N + (N−1)/(2β) (tolerance 0.05; anharmonic remainder O(β⁻²)) with ⟨m²⟩ → 1
- **Physics note:** over-relaxation alone is deterministic and energy-conserving, hence **not ergodic** — this is physics, not a defect; the validated production mode is composition with an ergodic kernel (`over_relaxation.rs`)
- **NOT validated:** —

### Hybrid (`HybridCore<A, B>`)
- **Validated:** Hybrid(Metropolis, Wolff) ⟨E⟩ and ⟨m²⟩ vs exact enumeration on 4-site Ising (p2_remaining)
- **All Ising pairwise compositions (2026-08-19):** every pairing of the four Ising-capable kernels — Metropolis+Wolff, Metropolis+SW, Metropolis+HeatBath, Wolff+SW, Wolff+HeatBath, SW+HeatBath — vs full 2^4 enumeration on the 4-site ring at β=0.5: ⟨E⟩, ⟨m²⟩ and C = β²(⟨E²⟩−⟨E⟩²), 8 seeds each, |z| < 4 per seed, |z̄| < 2 per combo, pooled one-sided Σz gate across the matrix
- **Composition boundary semantics (2026-08-19):** same-seed determinism; Hybrid(A, B) ≡ manual `A; B` sequencing bit-for-bit; `repetitions(2, 3)` ≡ `A; A; B; B; B`; the nested combinator closure Hybrid(A, Hybrid(B, C)) ≡ `A; B; C`
- **Continuous compositions (2026-08-19):** O(2) Wolff+SW, Metropolis+Wolff and Metropolis+ContinuousHeatBath vs a pure-Wolff reference on 8×8 at β=0.9 (pooled |z| < 4 on ⟨E⟩ and ⟨m²⟩); Metropolis+Microcanonical validated in the microcanonical section (quadrature, limits, Wolff cross-solver)
- **NOT validated:** —

### MultiSpinIsing (64-replica bit-parallel, `MultiSpinIsing`)
- **Validated:** 8-site exact enumeration cross-check (`p2_validation.rs`); acceptance LUT anchors, PT weight ratio, per-replica array observables (lib tests)
- **Multi-seed z coverage (2026-08-19):** 8-site PBC chain, β ∈ {0.4, 0.8}, 8 seeds — ⟨E⟩, ⟨m²⟩ and C = β²(⟨E²⟩−⟨E⟩²) vs full 2^8 enumeration, |z| < 4 per seed, |z̄| < 2 per (β, observable), pooled one-sided Σz gate
- **All 64 replicas ensemble-consistent (2026-08-19):** the `Energy_replica` array observable (averaged over replicas and time) matches the same exact ⟨E⟩ (|z̄| < 1) — replica 0 is not special; the replicas are exchangeable valid Metropolis chains
- **Cross-solver (2026-08-19):** vs scalar per-site `MetropolisCore` on identical physics (same lattice, β, J) at both temperatures — pooled cross-solver z on ⟨E⟩, ⟨m²⟩ and ⟨|m|⟩, all |z| < 4 (`multispin_cross_solver.rs`)
- **NOT validated:** —

### Wang-Landau (`WangLandauCore`)
- **Validated:** DOS matches exact 4×4 Ising enumeration (un-ignored, ~11s)
- **Validated:** 2-site exact DOS levels/degeneracies
- **Validated:** Canonical reweighting recovers exact ⟨E⟩
- **BinnedAxis production run:** binned DOS vs exact binned degeneracies (8-level weighted 6-ring, 14-bin axis; per-bin |Δln g| RMS 0.013, gate 0.05); canonical reweighting ⟨E⟩(T) at 3 temperatures vs exact (|z|≤1.2); flat-histogram-only route error floor documented (RMS ≤0.1)
- **Convergence robustness (2026-08-19):** an unattainable `minimum_visited_fraction` (more visited bins demanded than physically reachable) terminates loudly with `WangLandauTermination::UnreachableBins` after the discovery plateau is established (500 stalled flatness checks), instead of silently burning sweeps to the maximum-sweep guard; checkpoint round-trips the failure and its evidence; version-1 checkpoints without the stall fields still load
- **NOT validated:** —

### Multicanonical / umbrella (`EnergyBiasCore`)
- **Validated:** Transactional rejection, bias algebra, axis boundaries
- **Validated:** Discrete DOS reweighting matches enumeration
- **Full MC-vs-exact:** real EnergyBiasCore sampling with exact-DOS bias on the 6-site ring — ⟨E⟩(T), ⟨m²⟩(T) at β ∈ {0.2, 0.5, 1.0} (max|z|=3.4 over 16 seeds, |z̄|≤0.65), full P(E) at β=0.5 within 4σ binomial band (TV distance 0.005 vs gate 0.02), biased histogram flatness ≥0.90 (canonical would give 0.003)
- **NOT validated:** — (full MC distribution validated 2026-08-14; experimental status lifted)

### Worm (Ising HT graph) (`WormKernel`)
- **Validated:** HT graph partition identity matches spin enumeration
- **Validated:** Graph energy estimator matches exact spin energy
- **Validated:** Hastings reciprocity, endpoint correlation
- **Cross-solver (2026-08-19):** 4×4 square ⟨E⟩ at β=0.44 agrees with spin Metropolis within pooled 4σ
- **Multi-component lattices supported (2026-08-19, second pass):** the HT-graph ensemble factorizes over connected components, so `IsingGraphWormMC::from_lattice` / `IsingGraphWormEnsemble` (`src/worm/ensemble.rs`) run one independent two-defect worm per component on domain-separated derived streams (`carlo_rs::RngStreamKey`, component index in the replica field; one salt per component per sweep from the shared context stream, so no RNG state is hidden from a checkpoint). Observables combine additively; total energy is measured under the all-physical conditioning (the product ensemble is preserved); endpoint correlations are per component (cross-component pairs have no worm estimator — they factorize to zero); isolated sites form trivial components sampled exactly by the empty graph. Validated vs full 2^8 spin enumeration on a 4-ring + 3-chain + isolated-site lattice (total and per-component ⟨E⟩, two-point correlations, and ⟨m²⟩ reconstructed from the worm pair correlations); partition identity Z_spin = 2^N Π cosh(βJ_e) Π_c Z_graph,c at 1e-10; cross-solver vs spin Metropolis on two disjoint 4×4 squares (pooled |z| < 4); v2 multi-component snapshots round-trip bit-exact trajectories, v1 single-component snapshots still load (and are rejected loudly for multi-component ensembles)
- **Input rejection:** genuinely invalid input is still rejected loudly at construction — empty lattice, non-finite/negative β, non-finite coupling, `J · weight < 0` on any edge, self-loops. The raw `IsingGraphWormModel` + `WormKernel` pair additionally requires a **connected** lattice: its single defect pair would silently freeze the other components, so `IsingGraphWormModel::new` rejects disconnected input for direct users while the ensemble adapter handles it by decomposition
- **NOT validated:** multi-defect / multi-leg (multi-component) worm algorithms — not implemented; the kernel is a two-defect kernel by design (documented in `worm` module docs)

### Kawasaki dynamics (`KawasakiCore`)
- **Validated:** Signed magnetization conservation (exact)
- **Validated:** Energy decreases on cooling (directional)
- **Validated:** High-T equilibration
- **Quantitative equilibrium:** sector-restricted exact validation — BFS over the exchange graph from the initial state pins the reachable fixed-M sector (full sector on 4-ring M=0 and 8-ring M=+2; no hidden invariants), exact sector Boltzmann reference, ⟨E⟩ at β=0.3/0.8 with multi-seed |z|<2.4; sector ⟨E⟩ demonstrably ≠ canonical ⟨E⟩ (0.08 vs −3.20 at β=0.8)
- **NOT validated:** — (quantitative equilibrium validated 2026-08-14; superseded a physically-wrong cross-seed zombie test, see Known issues 3)

### BKL / n-fold-way (`BklIsingKernel`)
- **Validated:** Exact-trajectory reproducibility (bit-exact checkpoint)
- **Validated:** Fixed-time sampling matches exact small-Ising energy
- **Long-time equilibrium:** residence-time-weighted ⟨E⟩ and ⟨m²⟩ vs exact enumeration on N=4 and N=8 at β ∈ {0.2, 0.6, 1.0} (8 seeds; per-seed stderr from time-blocked means — per-visit averaging is length-biased and would be wrong, the estimator is pinned explicitly)
- **NOT validated:** — (long-time equilibrium validated 2026-08-14)

### Gillespie (`GillespieKernel`)
- **Validated:** Rate selection + exponential wait-time mean
- **Validated:** Absorbing-state clock advance
- **Multi-state equilibrium:** 3-state asymmetric CTMC — occupancy fractions vs exact stationary π (solved in-code via Cramer, πQ=0 verified to 1e-12); π = (0.1955, 0.4965, 0.3080); plus Ising-via-Gillespie ⟨E⟩ vs exact at β=0.6 (|z|≤2.1)
- **NOT validated:** — (multi-state equilibrium validated 2026-08-14)

### Event chain (`HardSphereEventChain`)
- **Validated:** Collision geometry, lifting at exact collision, PBC wrapping
- **Validated:** Snapshot restore
- **Equation of state:** contact-value identity Z = 1 + 2η(1−1/N)g(σ⁺) — event-chain collision rate (sausage average + Richardson extrapolation over chain lengths 1σ/2σ) reproduces the exact hard-disk virial series through B₃ (B₂ = πσ²/2, B₃ = 1.9295σ⁴) at η = 0.04/0.07 within 4σ; B₃ constant independently cross-checked by Mayer triple-integral MC (4σ)
- **Cross-solver:** g(σ⁺) at η=0.2 agrees with particle Metropolis NVT shell estimator within pooled 4σ (both in literature window 1.15–1.5); long variant adds η=0.12 and η=0.3
- **NOT validated:** — (EOS and pressure validated 2026-08-14; experimental status lifted)

### Particle NVT (`ParticleMetropolisCore`)
- **Validated:** Two-particle analytic pair energy
- **Validated:** Energy distribution vs quadrature
- **Validated:** Cache integrity, fixed-seed reproducibility
- **NOT validated:** Multi-particle EOS

### Particle NPT (`ParticleNptMetropolisCore`)
- **Validated:** V increases when P decreases (directional)
- **Validated:** V(P1)/V(P2) > 1.02 (non-trivial response)
- **Validated:** Finite-N ideal gas exact: ⟨V⟩ = (N+1)kT/P (long test, `npt_ideal_gas_volume_matches_finite_n_exact`)
- **Resolved:** Earlier "equilibrium volume mismatch" was a missing finite-N correction in the reference formula, not a solver bug

### Particle μVT (`ParticleGrandCanonicalCore`)
- **Validated:** N increases with μ (directional)
- **Validated:** N(μ1)/N(μ2) > 1.02 (non-trivial response)
- **Validated:** Ideal gas Poisson ⟨N⟩ exact (long test, `muvt_ideal_gas_particle_number_matches_poisson_most_probable`; plus `ideal_gas_grand_canonical_number_mean_is_poisson`)
- **Resolved:** Same finite-N reference correction as NPT

### Rigid molecule (`MolecularMetropolisCore`)
- **Validated:** Bond-length preservation, geometry preservation
- **Equilibrium distribution:** three analytic cases against in-code quadrature references through the real solver (translation + plane-rotation moves): two-molecule pair ⟨U⟩ and bound fraction (1D Simpson); dumbbell+atom probe nematic ⟨cos 2α⟩ and ⟨U⟩ (2D midpoint); rotor-pair alignment ⟨cos 2Δθ⟩ across a coupling sweep ε=1→3 (linear response → saturation, the Langevin-x analog; max|z|=1.26 default, 7-coupling long variant max|z|=1.7). Thermalization-length pitfall documented (20k+ sweeps needed at strong coupling)
- **External field (2026-08-19):** one-body dipolar term `DipolarExternalField` (per-atom charges, wrap-safe minimum-image dipoles, non-neutral molecules rejected loudly); free-rotor equilibrium vs the analytic Langevin-dipole answers through the real kernel — 2D: ⟨cosθ⟩=I₁(x)/I₀(x), ⟨cos²θ⟩=(1+I₂/I₀)/2; 3D: ⟨cosθ⟩=L(x), ⟨cos²θ⟩=1−2L(x)/x; x=βpE grid 0.5–5, per-seed |z|<4; machine-precision identity `external_field_energy = −E·μ` (1e-12) every sweep
- **NOT validated:** — (external-field Langevin case validated 2026-08-19)

### Read-only topology capability (`GraphView`, F1, Experimental, 2026-09-03)
- **Owned/borrowed parity:** `CsrLattice` and `BorrowedUndirectedCsr` agree on vertex/edge counts, physical endpoints, and every incidence for a fixture combining disconnected components, isolates, parallel edges, and a self-loop.
- **Representation invariants:** constructor tests reject missing/nonzero/decreasing/wrong-end offsets, incidence length mismatch, out-of-range neighbors/edge IDs/endpoints (including `usize::MAX` edge ID), incidence-to-endpoint mismatch, incorrect physical-edge multiplicity, two same-direction incidences, a third self-loop incidence, and a large unreferenced edge table. Empty borrowed graphs are explicitly valid as `offsets = [0]`; private helper coverage checks reportable validation-capacity overflow.
- **API invariants:** IDs are publicly produced only by the current view, have checked dense bounds and transparent `usize` layout, but carry no graph provenance. The public incidence contract requires only `Iterator`; directed out/in signatures and their compile fixtures remain crate-private pending D1.
- **Validation memory:** temporary scratch is one byte per physical edge via `try_reserve_exact` and is dropped before returning the zero-payload borrowed view. Reportable reserve failures return `ValidationCapacity`; abort-on-OOM allocators remain outside the recoverability claim.
- **Performance:** the repeated three-path `topology_view` benchmark records no >5% regression for generic owned or generic borrowed edge/incidence scans. The allocation probe measures constructor scratch separately and records zero allocations/bytes across 128 post-construction scans; machine-specific results are in [PERCOLATION.md](PERCOLATION.md). No CI timing threshold.
- **Not implemented:** public owned/borrowed directed CSR (D1). The independently reviewed topology-generic F3 analyzer and the F4 `StaticPercolationMC` composition remain Experimental until the platform merges into `dev`.

### Activity configuration and Bernoulli laws (F2, Experimental, 2026-09-03)
- **Ownership and API:** `StaticConfiguration` owns only private vertex/edge masks. The mask type, `law` module, concrete variants, and storage counters are not public; reviewed configuration/activity/probability/law entry points are selectively re-exported.
- **Validation and atomicity:** `Probability` rejects NaN, infinities, and values outside `[0, 1]`. Sampling checks configuration and every heterogeneous-field length before mutation or RNG consumption. Unit tests preserve the complete logical activity
configuration after a mixed edge-length error and compare the next 16 RNG words against an untouched clone. Uniform site, bond, and all four mixed `p=0/1` combinations consume no RNG; heterogeneous fields claim only per-element endpoint skipping.
- **Correctness:** unit tests cover exact endpoints, fixed-seed reproduction, pure/mixed reductions, borrowed/owned fields, a fixed-RNG heterogeneous reference, zero-size borrowed graphs, ID bounds, stale clearing, and resize/capacity reuse.
- **Statistical:** `tests/physics/percolation_law_zscore.rs` runs 16 seeds by default, scaled through the full `SCUTTLE_ZSCORE_SEEDS=1..=4096` range, with 512 samples/seed. Each of six domains (uniform site, uniform bond, uniform mixed vertex/edge, heterogeneous mixed vertex/edge) has one aggregate gate over all seeds. For probabilities `p_i`, the exact total moments are `mean = seeds * samples * sum(p_i)` and `variance = seeds * samples * sum(p_i * (1-p_i))`; uniform fields reduce to the binomial formula and heterogeneous fields are Poisson-binomial. Each aggregate requires `|z| < 4.5`; across six approximately normal gates this gives a family-wise false-alarm probability of about `4.1e-5`. A zero-variance field instead requires exact equality. No per-seed maximum, sign-fraction, or pseudo-exact dispersion gate is used.
- **Performance:** uniform endpoint probabilities use logical all/none states without RNG loops. The 14-ID `percolation_laws` Criterion target isolates sampling with correct `V`, `E`, and `V+E` throughput denominators. Across 25 law/topology allocation records, the probe
reports zero allocations/bytes for 128 post-warmup iterations of site, bond,
mixed, heterogeneous mixed, and two-sample endpoint switching. Recorded measurements are in [PERCOLATION.md](PERCOLATION.md).
- **Limitations:** the API passed F2 review but remains provisional and
  Experimental until the platform merges into `dev`. Heterogeneous fields are only an Experimental F2 substrate; N1 validation/facades remain Not implemented. The independently reviewed F3 analyzer/workspace and the F4 adapter compose into the Validated static family (see below).

### Undirected component analyzer (F3, Experimental, 2026-09-03)
- **Semantics and API:** `analyze` composes `UndirectedGraphView`, `VertexActivity`, and `EdgeActivity` without a mode enum. An active physical edge is counted and joined only when its edge and both endpoints are active; self-loops count once. `ComponentSummary` reports active vertices/edges, component count, canonical-minimum largest identity with deterministic lowest-ID tie breaking, largest size, and raw `M2`. `AnalysisResult<'_>` borrows reusable query outcomes and optional labels from `ComponentWorkspace`, preventing another analysis while results are live.
- **Boundary queries:** construction from raw indices validates range before allocation, sorts and removes duplicates, accepts empty sets (always false), and defines overlap as crossing when an active common vertex exists. Several queries are answered from one component result without cloning query storage per sample.
- **Reference and edge cases:** eight F3 integration-test entries include exhaustive 2x2 configurations for site/bond/mixed activity; 7 topologies x 3 laws x 128 seeded configurations against an independent incidence/stack flood fill; owned/borrowed CSR parity; self-loops, parallel edges, disconnected graphs, isolates, empty borrowed CSR, overlap, empty/duplicate/multiple queries, all summary fields, labels partition, edge-order-independent identity, typed vertex/edge/query errors, and grow/shrink/labels/query workspace reuse. An internal unit test forces the query generation stamp through `u64::MAX` and verifies clearing before reuse.
- **Atomicity:** activity and query dimensions are validated before workspace mutation. Capacity errors are typed; successful earlier `try_reserve` calls may retain capacity if a later reserve fails, but no partial scientific result is returned. Abort-on-OOM allocators remain outside the recoverability claim.
- **Performance:** 15 allocation records (5 topologies x 3 laws, 128 post-warmup sample+analyze calls each) report zero allocations and zero bytes per sample after `ComponentWorkspace::prepare`, versus the retained F0 reference's 7 allocations / 86,016 bytes per sample. Feature-gated audit accounting measures all six workspace Vec capacities as 196,609 bytes at V=4096, Q=1, labels enabled, exactly matching allocator-observed prepare bytes. The 48-ID Criterion run reports samples/s and edges/s and compares every topology/law against an equivalent independent flood fill after 128-sample parity. Fresh-process Linux `/proc` probes record square and sparse-ER process `VmHWM`; full data are in [PERCOLATION.md](PERCOLATION.md).
- **Maturity and limits:** passed independent review and engineering gates; still Experimental. F4 now composes this analyzer into `ObservablePlan` observables and the `StaticPercolationMC` adapter (see below).

### Static observables and production adapter (F4, `StaticPercolationMC`, Validated — review passed, 2026-10-02)
- **Observable semantics:** `ObservablePlan` explicitly requests `ActiveVertexCount`, `ActiveEdgeCount`, `ComponentCount`, `LargestSize`, `GiantFraction` (denominator fixed to total topology `V`), `RawSecondMoment` (`M2 = sum s_i^2`, unnormalized, largest included), and `FiniteClusterSusceptibility`. The susceptibility removes exactly one largest component per sample — equal sizes choose the canonical lowest identity, matching the F3 tie policy — before averaging, so it cannot be reconstructed from separately aggregated means (unit-tested counterexample against the ratio of means). `chi = (M2 - S_max^2)/(V_active - S_max)`; samples with a zero denominator are undefined: the chi scalar is not written for them (its `n_bins` counts only defined samples), while the always-recorded `FiniteClusterSusceptibilityDefined` indicator (0/1 per sample) keeps the defined fraction explicit in `results.json`. Boundary queries record `BoundaryCrossing0`, `BoundaryCrossing1`, ... in plan order. Labels are intentionally not part of the plan; `analyze_with_labels` serves direct label consumers.
- **Conversion policy:** moment observables convert u128 sums of squares up to `V^2` into f64. Plans whose topology could exceed the exact 2^53 f64 integer range (V > 94,906,265) are rejected at construction (`ObservablePlanError::InexactMomentConversion`) for both `RawSecondMoment` and the susceptibility numerator, so no silent precision loss is possible.
- **Adapter composition:** `StaticPercolationMC` composes an owned `CsrLattice`, the uniform `StaticLaw` enum (site/bond/mixed, only at this adapter boundary), `StaticConfiguration`, a prepared `ComponentWorkspace`, and the `ObservablePlan`. `sweep()` only resamples occupancy; `measure()` analyzes with the preallocated workspace and records exactly the planned scalars. Borrowed graphs and custom laws bypass the adapter through `analyze` or `Run::from_parts()`.
- **Parameter schema:** `mode` site/bond/site-bond with strict `p` vs `p_site`/`p_bond` mutual exclusion; lattice parameters reuse the standard builders; `pbc` defaults to `false`; spanning defaults are square left/right columns and chain endpoints, everything else requires explicit `spanning_from`/`spanning_to` comma lists (always both, or the adapter is constructed directly with a `BoundaryQuery`); `observables` defaults to the full plan. Invalid modes, probabilities, observable names, site lists, mismatched dimensions, and missing probabilities are typed rejections (unit-tested error paths).
- **Migrated validation:** every PR #4 asset now runs against the final API in `tests/physics/percolation.rs` (14 tests, 2 `#[ignore]` long) — 2x2 exhaustive enumeration vs hand-derived moments, site spanning polynomial, mixed closed form `2 p_s^2 p_b - p_s^4 p_b^2`, mixed-to-pure reduction identities, p = 0 / p = 1 occupation extremes, chain closed forms, union-find vs flood-fill parity across seven topologies x three modes, crossing monotonicity (pure and per-probability mixed), scheduler 2x2 exact moments, scheduler chain closed forms, fixed-seed bitwise reproduction, and the ignored 32x32 bond p_c = 1/2 and 16^3 cubic critical-bracket long tests. `tests/physics/percolation_zscore.rs` gates `BoundaryCrossing0` and `LargestSize` for site, bond, and mixed modes over 16 default seeds (SCUTTLE_ZSCORE_SEEDS-scalable, |z| < 4, |zbar| < 1.5, no one-sided bias).
- **New F4 coverage:** per-sample tie exclusion (canonical largest), no-active-vertices and V=0, single-giant-component undefined chi, multi-finite-cluster values, per-sample exclusion vs ratio-of-means counterexample, undefined-chi Results behavior and indicator counts under the scheduler, JSON schema and observable selection, and the full FromParams rejection matrix. Adapter unit tests (observable.rs + carlo.rs) add 10 cases; the whole percolation unit+integration set is green with `SCUTTLE_ZSCORE_SEEDS=64`.
- **Performance:** 15 adapter allocation records (5 topologies x 3 laws, 128 post-warmup sweep+measure cycles each) report 0.000 allocations and 0.000 bytes per cycle in steady state; the Criterion target adds 15 `static-mc-sweep+measure` IDs composing sampling, analysis, and recording. F3 core throughput paths are unchanged code whose parity and sequence assertions still pass; same-run adapter-vs-core comparison shows measure overhead within contention noise, and cross-day point estimates carry documented host-load inflation (see [PERCOLATION.md](PERCOLATION.md)).
- **NOT validated (unchanged scope):** heterogeneous probabilities (N1), wrapping (W2), parameter scans (Z1), process families; chi tie conventions other than single canonical-largest exclusion are separate named policies and do not exist.

### Percolation, site / bond / mixed (`PercolationMC`, 2026-09-02 — superseded)
- **Status:** Experimental PR #4 reference implementation, now **superseded**:
  the API (`PercolationMC`, `PercolationMode`, `OccupancyState`, `cluster_stats`,
  `UnionFind`) has been deleted and every validation asset below is migrated to
  the F4 `StaticPercolationMC` API (observable names updated: `MaxCluster` ->
  `LargestSize`, `SecondMoment` -> `RawSecondMoment`, `NClusters` ->
  `ComponentCount`, `Spanning` -> `BoundaryCrossing0`, `Occupied` ->
  `ActiveVertexCount`/`ActiveEdgeCount`). See
  [PERCOLATION.md](PERCOLATION.md) for the four-state support matrix, scientific
  definitions, limitations, and benchmark baselines.
- **Validation assets:** Ordinary site, bond and mixed site-bond percolation on arbitrary `CsrLattice` graphs (i.i.d. occupancy resampling, union-find cluster analysis; mixed connects a bond only when it and both endpoint sites are open). 2×2 open square: full 16-configuration enumeration vs hand-derived closed-form moments for site and bond — ⟨MaxCluster⟩ = 30/16 and 45/16, ⟨sum(s_i²)⟩ = 76/16 and 164/16, ⟨NClusters⟩ = 17/16 and 33/16, P(spanning) = 7/16 and 12/16 at p = 1/2; site spanning matches the polynomial 2p²(1−p)² + 4p³(1−p) + p⁴ across p ∈ {0.2, 0.44, 0.5927, 0.8}
- **Mixed closed form (hand-derived):** 2×2 site-bond P(span) = 2·p_s²·p_b − p_s⁴·p_b² (only an active horizontal bond crosses; the two rows coincide only when everything is open) — exact at (p_s, p_b) ∈ {(0.6,0.7), (0.9,0.4), (1,0.5), (0.5,1), (0.3,0.3)}; reduces exactly to the pure-mode values in both limits
- **Reduction identities (strict):** mixed at p_site = 1 reproduces pure bond moments, at p_bond = 1 pure site moments (all four, 1e-12, at p ∈ {0.2, 0.5, 0.8})
- **Independent algorithm cross-check:** `cluster_stats` (union find) vs an in-test flood-fill reference sharing no algorithmic path, configuration-by-configuration, all three modes — exhaustive on chain-8, square-3x3, cubic-2x2x2, triangular-2x2, honeycomb-2x2, kagome-2x2 site, random-graph site (≈22k configurations); seeded random configurations beyond
- **1D exact solution:** open chain P(span) = p^L (site), p^(L−1) (bond), p_s^L·p_b^(L−1) (mixed) — exact enumeration at three (p_s, p_b) points and through the full scheduler stack at L = 6 for all three modes (100k sweeps, |z| < 4)
- **Scheduler end-to-end:** 200k i.i.d. sweeps on the 2×2 square reproduce all four enumerated moments within |z| < 4; p = 0 and p = 1 boundary behavior exact in unit tests (no span / single spanning cluster); fixed seed reproduces bitwise for pure and mixed modes, different seed shifts the stream
- **Statistical:** 4×4 site at p = 0.6 (2¹⁶ = 65536 configurations), 3×3 bond at p = 0.55 near p_c = 1/2 (2¹²) and 3×3 mixed at (0.6, 0.7) (2²¹ = 2097152 configurations, shared across its two tests) fully enumerated as references; 16-seed z-scores on `Spanning` and `MaxCluster` for all three modes (|z| < 4, |z̄| < 1.5, no one-sided bias) — `tests/physics/percolation_zscore.rs`
- **Physics sanity:** crossing probability monotone non-decreasing in p (8×8, site and bond, p = 0.1…0.9) and separately in p_site / p_bond for mixed (fixed coordinate at 0.9); overlapping spanning-set degeneracy pinned by unit test; probability-parameter mismatches (p in mixed mode, p_site/p_bond in pure modes) rejected loudly
- **Critical-point checks (long, `#[ignore]`, nightly):** 32×32 bond at p_c = 1/2 (exact self-duality) → crossing within 0.06 of 1/2; 16³ cubic bond brackets the critical region — P(cross) < 0.05 at p = 0.12 and > 0.95 at p = 0.40 around p_c ≈ 0.2488 (no unproven 3D duality assumed)
- **NOT validated:** invasion/kinetic percolation variants (not implemented); spanning defaults limited to square/chain (other graphs require explicit site sets, rejected loudly otherwise); critical crossing on pbc triangular/honeycomb/kagome (builders are periodic-only; no clean crossing convention); no exact mixed critical line exists for the square lattice (validated via reductions and closed forms instead)

## Input-validation coverage (criterion G)

`tests/physics/input_validation.rs` (2026-08-19) sweeps the parameter and
constructor surface of every scheduler-ready solver; invalid input must
return `InvalidConfig` (or the kernel error type), never a silent accept and
never an unintended panic:

- Lattice kernels (Metropolis/Wolff/SW/heat baths/continuous heat
  bath/microcanonical/hybrid): negative/NaN β, unknown lattice names, zero
  or malformed dimensions, OBC triangles, Potts q<2, non-finite J, unknown
  `initial_state`.
- MultiSpinIsing, Wang-Landau (fraction/flatness/log_f/interval ranges,
  24-site exact-axis limit), worm (β, coupling sign, close probability,
  fugacity), Kawasaki/BKL (incl. the fixed β/J validation — previously
  `J = NaN` panicked in an assert-backed constructor), BKL event windows,
  event chain (box/particles/diameter/chain lengths), particle NVT/NPT/μVT
  (density, cutoff, displacement, pressure, activity) and the molecule
  kernel (scales, D≥2, topology corruption).

## Statistical validation framework

- **z-score tests:** 16 independent seeds per solver (8 where chains are more expensive), |z|<4 per seed, |z̄|<2 mean, no one-sided bias. Seed counts scale via `SCUTTLE_ZSCORE_SEEDS` (unset → default unchanged; nightly `zscore-monitor` uses 64; `just nightly-zscore` reproduces locally). Σz thresholds are scale-invariant (−2√n); in test files with many configuration cells (Potts, external field) the Σz gate is pooled over all scores of one solver/test to control the multiple-comparisons false-alarm rate.
- **Cross-solver:** Metropolis vs Wolff pooled z-scores agree (|Δz̄|<2); Potts 4-solver pairwise |z|<4 at βc; worm vs Metropolis pooled |z|<4; MultiSpinIsing vs scalar Metropolis pooled |z|<4; SW-continuous and over-relaxation-composition vs Wolff pooled |z|<4
- **Connectivity:** Explicit Markov chain enumeration on N=2 Ising (4 states), BFS strong connectivity + aperiodicity check

## Known issues

1. ~~NPT/μVT equilibrium values~~ → Resolved 2026-08-14 (fixed in e1a07e4): the finite-N reference formulas (⟨V⟩=(N+1)kT/P, Poisson ⟨N⟩) match exactly; the old "mismatch" was a test-side formula error.
2. **strict-repro feature:** Was defined in Cargo.toml but had zero implementation. Removed.
3. ~~Kawasaki cross-seed zombie test~~ → Removed 2026-08-14: `kawasaki_2d_ising_energy_converges_same_regardless_of_seed` was #[ignore]'d with a reason stating it cannot pass (random starts land in different fixed-M sectors whose ⟨E⟩ genuinely differ), yet `--ignored` runs executed it anyway and it failed deterministically at HEAD. Physically-wrong criterion; superseded by the sector-restricted exact validation in `kawasaki_exact.rs`.
4. ~~Kawasaki/BKL β/J validation~~ → Fixed 2026-08-19: `from_params` previously forwarded user `beta`/`J` into assert-backed constructors, so `J = "NaN"` panicked instead of returning `InvalidConfig`. Both adapters now validate through a shared `validate_kinetic_ising_params`.

## Completion log

| Date | Task | Result |
|------|------|--------|
| 2026-07-23 | Gap 1: z-score framework | ✅ 5 tests (16 seeds × 3 solvers + cross-solver) |
| 2026-07-23 | Gap 2: Markov chain connectivity | ✅ 4 tests (Metropolis/Wolff/SW strong connectivity) |
| 2026-07-23 | Gap 3: Quantitative EOS | ✅ NPT/μVT directional + non-trivial response; documented equilibrium issue |
| 2026-07-23 | Gap 4: Validated domain docs | ✅ This document |
| 2026-08-14 | Production hardening: 8 solvers | ✅ 31 new tests: multicanonical full-MC (4), WL binned (4), event-chain EOS (4+1 long), molecule equilibrium (4+1 long), Kawasaki sector-exact (2), heat-bath discrete+O(N) (5), BKL/Gillespie equilibrium (4), continuous-spin cross-solver + frustrated triangle (6) |
| 2026-08-14 | Zombie test removal | ✅ Physically-wrong Kawasaki cross-seed test deleted; superseded by sector-exact validation |
| 2026-08-19 | Production: Potts q>2 | ✅ 8 tests (`potts_exact.rs`): full q^N enumeration on 2×2 square and N=8 chain, q ∈ {3,4}, β ∈ {0.3, 0.8}, observables ⟨E⟩/⟨m²⟩/C for heat bath + SW + Wolff; enumeration anchors (q=1 degenerate, q=2 ≡ Ising with J/2); β=0 uniform and β=8 frozen analytic limits; 8×8 q=3 4-solver cross-solver at βc; q=4 βc = ln 3 directional anchor |
| 2026-08-19 | Production: molecule external field | ✅ 4 tests (`molecule_external_field.rs`): `DipolarExternalField` API (additive, backward-compatible), 2D von Mises + 3D Langevin free-rotor moments, −E·μ machine-precision identity, loud rejection of non-neutral/short/non-finite charge tables |
| 2026-08-19 | Production: WL convergence robustness | ✅ 4 tests (`wang_landau_convergence.rs`): unattainable `minimum_visited_fraction` → loud `UnreachableBins` termination (not Converged/MaximumSweeps), auto-derived reachable set matches enumeration, checkpoint evidence guards, version-1 back-compat |
| 2026-08-19 | Production: multi-component worm | ✅ Honest limitation route: multi-component (disconnected/isolated-site) lattices now rejected loudly at `IsingGraphWormModel::new` (defect pair confined to one component would silently freeze the rest); multi-defect/multi-leg worms documented as not implemented; + worm-vs-Metropolis cross-solver test |
| 2026-08-19 | Production: input-validation audit | ✅ 10 tests (`input_validation.rs`) across all 19 solvers; fixed three source holes (Kawasaki/BKL β+J panic path, MultiSpinIsing lattice validation, BKL event-window validation) |
| 2026-08-19 | Production: multi-component worm implemented | ✅ `src/worm/ensemble.rs` + `IsingGraphWormMC::from_lattice` (per-component two-defect worms, `RngStreamKey` domain-separated streams, additive observables, v1/v2 checkpoints) + `CsrLattice::connected_components`; 4 suite tests (`worm_multi_component.rs`: exact 2^8 enumeration incl. per-component energies/correlations/worm-reconstructed m², partition identity at 1e-10, snapshot round-trip + loud rejections, cross-solver vs spin Metropolis on two disjoint 4×4 squares) + 3 lib tests; multi-component input no longer rejected at the scheduler surface — genuinely invalid input still is |
| 2026-08-19 | Production: over-relaxation in composition | ✅ 6 tests (`over_relaxation.rs`): reflection-map machine-precision identities + deterministic DB (involution/isometry, O(2) angle form, O(3) orthogonal π-rotation); kernel ≡ manual reflection bit-exact with per-update \|ΔE\|<1e-12; Hybrid(Metropolis, Microcanonical) vs exact XY-ring quadrature from 2 inits; analytic limits (β→0 exact, β=8 spin-wave); cross-solver vs Wolff 8×8 O(2)/O(3) |
| 2026-08-19 | Production: all hybrid compositions | ✅ 3 tests (`hybrid_compositions.rs`): all six Ising pairwise combos vs exact enumeration (E, m², C; multi-seed + pooled Σz); boundary semantics bit-exact (sequencing/repetitions/nesting/determinism); continuous O(2) combos vs Wolff 8×8 |
| 2026-08-19 | Production: MultiSpinIsing | ✅ 2 tests (`multispin_cross_solver.rs`): multi-seed z vs exact enumeration (E, m², C at β=0.4/0.8) + all-64-replica array-observable consistency; cross-solver vs scalar Metropolis (E, m², \|m\| pooled z) |
| 2026-08-19 | Production: SW continuous spins | ✅ 4 tests (`sw_continuous.rs`): XY-ring exact quadrature (E, m², cos Δθ), O(3) analytic limits (β→0, β=8 spin-wave with the zero-mode-counted formula), cross-solver vs Wolff 8×8 O(2)/O(3) |
