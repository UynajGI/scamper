# CMC.rs Percolation Status and Scientific Definitions

> F0 baseline established 2026-09-03. This document describes the current
> `feat/cmc-percolation` reference implementation and freezes vocabulary for the
> staged platform roadmap. It does not make the provisional API stable.

The normative design and phase order are the [platform philosophy](../docs/metaphysics/2026-09-03-percolation-platform-philosophy.md),
the [v1 implementation roadmap](../docs/plans/2026-09-03-percolation-platform-roadmap.md),
the [v2 research roadmap](../docs/plans/2026-10-23-percolation-roadmap-v2.md),
and the [v3 pluggable-framework roadmap](../docs/plans/2026-10-25-percolation-roadmap-v3.md).

## Status vocabulary

- **Supported:** an API with explicit semantics exists, but not every production
  validation gate necessarily passes.
- **Validated:** the family-specific semantic, exact/reference, statistical,
  invalid-domain, reproducibility, and performance gates pass.
- **Experimental:** usable reference behavior exists, but it is not the intended
  production API or has not passed all production gates.
- **Not implemented:** no supported claim; a generic graph routine producing a
  number does not change this status.

## Support matrix

| Family or capability | Status | Current scope or reason |
|---|---|---|
| Read-only undirected topology capability | Experimental | Independently reviewed F1 `GraphView`/`UndirectedGraphView` contracts use dense IDs and static dispatch; `CsrLattice` and validated zero-copy `BorrowedUndirectedCsr` implement them. The API remains provisional until the platform merges into `dev`. |
| Directed topology capability (D1) | Validated | Public `DirectedGraphView` (dense `ArcId`s, `arc_endpoints` = `[source, target]`) with separate `DirectedOutGraphView`/`DirectedInGraphView` capabilities — directedness is structural, never a storage flag. Owned `DirectedCsr` (public fields + `validate`, parallel arcs and self-loops first-class: a self-loop is one out- and one in-incidence) and zero-copy `BorrowedDirectedOutCsr`/`BorrowedDirectedCsr` (the latter requires the transpose arrays). Construction validates offset shape, incidence lengths, ranges, per-arc sense on both sides, and exact arc multiplicity (missing or duplicated arcs are typed errors). Independent review passed 2026-10-22 (builder hand-traced, validator attacks predicted and confirmed empirically, SIMD explanation verified from layout evidence); see the D1 section below and `docs/reviews/`. |
| Static configuration and Bernoulli-law substrate | Experimental | F2 provisional API separates private activity storage from independent site/bond/mixed laws, validates probability domains and topology lengths, and samples allocation-free after preparation. The independently reviewed F3 analyzer/workspace and F4 observables build on this substrate; the layer remains Experimental until the platform merges into `dev`. |
| Open observable extension (v3 A1) | Experimental | `ObservablePlan::with_custom` accepts `Arc<dyn ScalarObservable>` (name + evaluate from a borrowed `SampleView`: total vertices, component summary, boundary outcomes, optional winding/labels); names are validated against built-ins, indicator columns, boundary columns, and other customs. Capability declarations (`requires_labels`, `requires_winding`, `requires_exact_moments`) are aggregated into `requests_labels`/`requests_winding`/`validate`, so every adapter materializes labels (including the bond-mask path via `analyze_bonds_with_labels`) and winding only on request, and rejects exact-moment plans above `V^2 = 2^53`. Built-ins are zero-cost `ObservableSpec` impls dispatched statically through the `StaticObservable` compatibility shim (deprecation attribute deferred to the v4 release boundary). Gates: 6 external-crate integration tests (recording, label request/omission, bond-adapter labels, winding + typed no-embedding rejection, name rejection matrix, byte-identical no-custom plans), 9 unit tests, rustdoc example, a 0-allocation steady-state probe, and the full 557-test cmc regression. |
| Open occupation/weight laws and sample views (v3 A2) | Experimental | `OccupationLaw<G>` is the static-dispatch law extension point: `validate_domain` runs once at adapter construction, `sample(&mut self, graph, configuration, &mut Rng)` draws one realization, and an associated `Error` keeps each family's typed rejection. Implemented by `SiteBernoulli`/`BondBernoulli`/`MixedBernoulli`, `StaticLaw` (which now delegates its whole draw through the trait), and `CorrelatedSiteLaw`; continuum and long-range draws remain processes (A4) because their edge universe is discovered while sampling. `WeightLaw` exposes the FPP `EdgeWeightLaw` transform (`validate` + `draw(u)`). `SampleView::new`/`vacuum`/`from_analysis` are public, so a third-party process or adapter can feed plans; `AnalysisResult` query/label accessors now borrow the workspace lifetime, and `ComponentSummary::EMPTY` is public for empty views. `StaticConfiguration::set_vertex_active`/`set_edge_active` are the checked public mutation surface for external laws. Gates: an external-crate `OccupationLaw` implementation analyzed through the public core; bitwise replay of every trait path against the corresponding inherent draw (Bernoulli, static, correlated), `WeightLaw` bitwise over 1001 uniforms, typed length/rate rejections before any RNG draw, and the full adapter regression. |
| Shared moment estimator (v3 A3a) | Experimental | `Moments` owns the count/sum/sum_sq arithmetic, biased-mean and unbiased-variance (clamped at zero) formulas, and the standard error once. It is allocation-free, keeps the exact floating-point operation order of the v2 expressions, and returns `Option` for undefined statistics. Four previously duplicated estimators now delegate: `PassageEstimator` (scalar passage time plus a per-target profile), `ShapeEstimator` (per-target profile with non-negative validation), `ChemicalDistanceEstimator` (per-entry reachability-conditioned moments plus an unconditional diameter), and `ScanEstimator` (per-scan reweighted values). Reachability and defined-sample policies stay explicit at the call site (`add` only finite entries; `Option` propagates undefined reports). Gates: 3 accumulator unit tests (empty/single/known/clamped), 156 percolation unit tests, 26 FPP + 25 L/C0 + 19 scan integration tests including the statistical brackets, and the full regression. |
| Generic occupation-law adapter (v3 A3b) | Experimental | `PercolationAdapter<L = StaticLaw>` owns the lattice, configuration, component/winding workspaces, and plan once; one `MonteCarlo` impl samples through `OccupationLaw` and records built-in + custom observables. `StaticPercolationMC` remains a compatibility alias for the Bernoulli preset, so its parameters, Results columns, winding path, and error enum are unchanged. Custom laws must use `Error = SamplingError`; `validate_domain` rejections surface as the existing typed `LawDomain` variant, and a post-validation sampling error is a contract panic (the preset behavior). Gates: external custom law + custom observable composed through `Run::from_parts` (2 custom tests), the full 558-test regression, clippy, and rustdoc. The bond-dynamics, continuum, and scan adapters still own their structs; the C0 writer is pending. |
| Versioned network I/O (v3 X1) | Experimental | Feature `io` (optional `serde`) adds exactly two explicitly versioned formats. Edge-list v1 is text with a mandatory `# cmc-edge-list v1` header plus required `# nodes`/`# layer` directives and `source target weight` lines; network JSON v1 is `format: "cmc.network"`, `version: 1` with named layers and per-edge strictly-increasing layer indices. Both parse into `NetworkDocument`, validate node/layer/weight/endpoint invariants, round-trip exactly, and project one layer to an owned CSR topology through `EdgeListSpec` policies while carrying bond weights. Unknown JSON fields, wrong versions/format strings, unsorted or empty layer lists, out-of-range endpoints/layers, non-finite or negative weights, malformed directives, and spec/document vertex-count mismatches are typed rejections; nothing is coerced. A two-layer document refuses edge-list serialization rather than dropping a layer. Gates: 7 unit tests (JSON and text round-trip, reference-graph equality, rejection matrix, multi-layer refusal), clippy with `io`, and a Criterion benchmark (16,384 edges: edge-list parse 2.56 ms / write 2.81 ms, JSON parse 4.17 ms / write 3.05 ms, layer-to-CSR 1.76 ms in a quick local run). Independent review pending. |
| Generic FSS core (v3 A5a) | Experimental | `BinderEstimator` accumulates independent `m = S_max/V` realizations into `U = 1 - <m^4>/(3<m^2>^2)` with a first-order delta-method error (including the `cov(m^2, m^4)` term from a third moment accumulator); it rejects empty topologies and `S_max > V`, and returns `None` before two samples or when every realization is empty. `FssCurve` validates one size on a strictly ascending probability grid; `crossing` finds the first sign change of the upper-minus-lower difference and linearly interpolates the crossing probability with a first-order interval, rejecting size-order violations, grid mismatches, and curves that never cross. Conventions are definitional (Binder's normalization; raw finite-size crossing), not literature values. A5b adds `bootstrap`, a percentile bootstrap of an arbitrary `FnMut(&[f64]) -> f64` statistic: a fresh seeded Xoshiro stream draws `n` indices with `random_range(0..n)` per resample into one reused buffer, so identical inputs replay bitwise; resample mean/standard error use `Moments`, intervals are nearest-rank 2.5/97.5 percentiles, and non-finite samples or statistic values are typed rejections. A5c adds `second_moment_correlation_length`: `xi^2 = sum_r r^2 C(r) / (2 d sum_r C(r))` (Amit–Martin-Mayor normalization) over non-negative finite `CorrelationPoint`s, rejecting zero dimensions, empty input, bad distances/values, and zero total weight. Gates: 11 unit tests (exact binary `U=1/3`, constant `U=2/3`, undefined cases, exact crossing at 0.5, exact flat-correlation length, exponential-decay recovery within 0.01, crossing/Binder/bootstrap/correlation rejection matrices, degenerate constant bootstrap, bitwise bootstrap replay + centering), a rustdoc example, clippy, and the full regression. The scaling-collapse fit and v2 bracket reproduction remain. |
| Two-layer mutual-component cascade (v3 X2) | Experimental | `MultiplexNetwork` pairs two `CsrLattice` layers over one dense vertex set (equal counts, non-empty). `mutual_component` alternates the largest component of layer A and layer B among the surviving vertices until the alive set stops shrinking, returning per-vertex membership, the surviving fraction, and the contraction count. Component extraction uses the shared analyzer with labels and the canonical lowest-ID tie rule, so the cascade is deterministic. Definitional anchor: Saberi 2015 review §V (archived): the order parameter is the fraction in the mutually connected largest component, with the cascade mechanism; the original model is Buldyrev et al., Nature 464, 1025 (2010). Scope: one-to-one dependency at full initial occupancy; random initial failures and >2 layers are not implemented. Gates: exact hand cascades (identical chains, split layer, star/path collapse to one vertex), an independent union-find reference map on four six-vertex graph pairs, typed vertex-count/endpoint/empty rejections, a rustdoc example, clippy, and the full regression. A quantitative bracket against the literature and preset adapters are pending. |
| Static site, bond, and mixed Bernoulli percolation on arbitrary owned or borrowed undirected graphs | Validated | All family-specific gates pass: F3 component analysis with reusable workspace, canonical identities, and multi-query crossing; F4 `ObservablePlan` observables and owned `StaticPercolationMC` with every PR #4 validation asset migrated, the zero-allocation adapter gate, and 16-seed z-scores; independent review passed. Production-ready once the branch merges into `dev`. |
| Heterogeneous occupation probabilities | Validated | N1 promotes the F2 `ProbabilityField` substrate to the `StaticPercolationMC` adapter: `StaticLaw::site_heterogeneous`/`bond_heterogeneous`/`mixed_heterogeneous` own one validated `Probability` per vertex and/or per physical edge, with adapter-constructor domain validation and scheduler scheduling through `Run::from_parts()`. Gates passed: exact non-identically-distributed small-graph enumeration (independently reproduced), bitwise uniform and mixed-pure reductions, endpoint determinism, typed length/invalid-probability rejection, 16-seed scheduler z-scores for site/bond/mixed (enumeration and Poisson-binomial references), uniform-vs-heterogeneous Criterion comparison, and a 0-allocation heterogeneous adapter probe. Independent review passed; production on the `dev` merge. The N2 random-failure facade now builds on this substrate. |
| Network random failure and degree-conditioned observables | Validated | N2 `RandomFailureMC` facade over the F2 uniform site law with an original-degree `NetworkProfile` and six C0 element-preserving degree arrays (`DegreeVertexCount`, `DegreeOccupied`, `ActiveFractionByDegree`, `GiantMembershipByDegree`, `DegreeFirstMoment`, `DegreeSecondMoment`); empty bins exact-zero with explicit occupancy, `V^3` exactness gate, 0-allocation probe. Gates: bitwise site-law reduction, independent-BFS enumeration z-scores on path/cycle/star/complete, relabeling invariance, empty-bin semantics, and literature anchors (Callaway et al. PRL 85, 5468 and Cohen et al. PRL 85, 4626; PDFs in `docs/ref/`). Independent review passed 2026-10-22 (kappa/q_c and the Callaway fixed points independently re-derived; star-graph observables hand-enumerated to exact fractions; V^3 boundary recomputed); see the N2 section below and `docs/reviews/`. |
| Targeted degree attack and robustness curves (N3) | Validated | `DegreeAttackMC` removes vertices in descending original degree with an explicit tie policy (`DegreeDescendingStableId`, seed-independent, or `DegreeDescendingSeededShuffle`, one realization per seed); one sweep is one complete removal trajectory recorded as ONE C0 vector sample (`RobustnessCurve`, `RobustnessSusceptibility` with the F4 defined-indicator convention, `RobustnessRemovedFraction` coordinates), never independent scalar samples. Grids: every removal, every k-th, or explicit checkpoints. Gates: per-step independent-BFS parity, strict survivor monotonicity, tie reproduction, coordinates/provenance through `Results` attributes, Callaway et al. Eq. (14) degree-cutoff anchors, Cohen et al. breakdown bracket; PDFs in `docs/ref/`. Independent review passed 2026-10-22 (both anchor graphs rebuilt independently; all eight Callaway values and both Cohen f_c values reproduced digit-for-digit); see the N3 section below and `docs/reviews/`. |
| Adaptive residual-degree attack (T1) | Validated | Roadmap v2 §2 T1 tail of N3: `AttackPolicy::ResidualDegreeDescendingStableId` / `ResidualDegreeDescendingSeededKey` recompute the removal order from residual degree (degree inside the still-active subgraph) after every removal — Holme et al. PRE 65, 056109 (2002) "RD" — on the same adapter, grids, C0 arrays, and provenance as N3. Order generated online by a lazy max-heap keyed (residual desc, tie key asc) with stale-entry skipping, `O((V+E) log V)` per trajectory; tie policies are explicit — ascending dense ID, or ONE fixed random `u64` tie key per vertex drawn at trajectory start (per-run shuffles are not well-defined for an online order). Gates: exact full-order parity vs a naive full-rescan reference at every step (28 cases: chain/cycle/star/dense-ER/multigraphs with self-loops and parallel edges, both policies, several seeds), per-step invariants, tie reproduction, Holme Fig. 7(e) statistical bracket (BA N=10^4, m0=5, m=3, R=24: RD below ID with disjoint empirical-SE 95% CIs for f >= 0.15, RD ~ ID with overlapping CIs at every grid point f <= 0.125, both attacks collapsed by f=0.30); the random-failure contrast at matched f is standard context, not part of Holme's figure. Also: asserted 0-allocation probe, heap-vs-naive Criterion ratio ~92x (independent re-timing 128x). Independent review passed 2026-10-22 (1492-case adversarial fuzz vs a from-scratch reference — zero mismatches; soundness argument re-derived; tie draw contract replayed bit-exactly; N3 orders bit-identical to 53b2019 via snapshot diff; Holme anchor wording and the BA-builder attribution corrected as follow-ups); see the T1 section below and `docs/reviews/`. |
| Directed static analysis (D2) | Validated | Three distinct questions, three distinct result types over the D1 capabilities: `analyze_weak_components` (union-find over the arc table, canonical lowest-ID labels), `analyze_strong_components` (iterative non-recursive Tarjan; `acyclic` tracks singleton-SCC self-loops separately), and `reachable` (early-exit BFS with stamped visitation; `reachable(v, v)` true by definition). Anchors: exhaustive enumeration of every digraph on 3 and 4 vertices vs transitive-closure ground truth (SCC, weak, and reachability all agree), Tarjan vs an independently written Kosaraju, weak vs undirected flood-fill, SCC transpose invariance, and closed forms for chains/trees/stars/rings/DAGs. Independent review passed 2026-10-22 (Tarjan hand-traced frame by frame; a third independent SCC implementation agreed on 300 fresh multigraphs; stamp rollover empirically tested beyond 2^32 queries); see the D2 section below and `docs/reviews/`. |
| Layered space-time directed percolation (D3) | Validated | `LayeredDirectedPercolationMC`: one sweep = one complete space-time realization over T explicit causal layers (double frontier bitset — memory scales with layer width, not space-time volume), recorded as ONE C0 vector sample (`DpLayerActive`, `DpSurvival` per layer; `DpExtinctionTime` scalar) — never T independent measurements. Directed bond DP with documented bitwise draw contract; single-seed and full-occupation initial conditions. Anchored to Hinrichsen Adv. Phys. 49, 815 (2000) (archived): frontier DP equals explicit-DAG path enumeration (D2 reachability engine); P(t) ~ t^-delta at p_c = 0.6447001 with delta = 0.159464 for both survival and full-seed density decay; sub/super-critical brackets. Independent review passed 2026-10-22 (a from-scratch RNG-identical model reproduced production traces bit-exactly; delta_hat converged at 8000 realizations); see the D3 section below and `docs/reviews/`. |
| Bootstrap percolation | Validated | E2 `BootstrapWorkspace` + `BootstrapPercolationMC` (roadmap §12): monotone activation-threshold closure over the E1 `ProcessQueue`/`NeighborCounters`/`StampFlags` primitives, generic over `UndirectedGraphView` with seeds as any `VertexActivity`, homogeneous `r` or per-vertex `u32` thresholds, and TWO explicit clocks — `SynchronousRounds` (Holroyd's parallel rounds) and `QueueClosure` (order-independent FIFO monotone closure) — recorded as `BootstrapActiveCount` (ONE C0 vector sample, fixed shape, truncate/pad convention) + `BootstrapFinalActive`/`BootstrapActivatedTotal`/`BootstrapRounds`|`BootstrapEvents` scalars. Gates: exhaustive 2^V seed masks on path/cycle/star/complete for every r in 0..=max_degree+1 AND seeded heterogeneous vectors vs an independent repeated-full-scan reference (closure sets, per-round trajectories, scalar invariants; ~2600 production runs); r=0 activates all / r>max_degree preserves seeds verbatim; FIFO queue closure = one-random-vertex-at-a-time and shuffled-sweep closures (test-RNG orders) = synchronous closure on random multigraphs; hand-worked heterogeneous cascade; bitwise replay of the one-f64-per-vertex density draw contract; Holroyd (PTRF 125, 195 (2003), arXiv:math/0206132, archived PDF) finite-size metastability brackets — his own full-activation observable I(L,p) is 1.0000 +/- 0.0013 at p+ = pi^2/(18 ln L) and 0.0000-0.0025 at p+/4 for L = 32/64/128, the measured half-majority density p_hat(L) = 0.0692/0.0577/0.0480 decays monotonically (the 1/ln L law) with p_hat ln L = 0.240/0.240/0.233 strictly inside (lambda/4, lambda), and the metastable window is sharper than a factor 1.5 (P(p_hat/1.5) <= 0.018, P(1.5 p_hat) >= 0.992); 0-allocation realization probe on five topologies; Criterion data. See the E2 section. |
| k-core percolation | Validated | E3 `kcore_decompose` + `KCoreMC` (roadmap §12): residual-degree peeling over the E1 `ProcessQueue`/`NeighborCounters`/`StampFlags` primitives, generic over `UndirectedGraphView`, uniform or per-vertex `u32` thresholds, arbitrary initial activity, synchronous-round survivor accounting (`KCoreSurvivorCount` as ONE C0 vector sample + `KCoreFinalSize`/`KCoreRounds`/`KCoreDeletedTotal` scalars), deterministic deletion order and waves exposed post-run. Gates: exhaustive 2^V initial subsets on path/cycle/star/complete vs an independent repeated-full-scan reference (1456 runs), k=0 removes nothing / k>max_degree removes everything in one round, queue-order independence across 90 shuffled relabeling orders, hand-worked heterogeneous thresholds, 1-core/tree closed forms, bitwise site-law reduction of the activity draw, and the Janson-Luczak (arXiv:math/0508453, archived) 3-core emergence bracket on G(n=5e4, m) — own solver reproduces `lambda_3 = 3.350919`, below c=lambda_3-0.15 the 3-core is extinguished (0.00000 ± 0.00000), above c=lambda_3+0.15 it measures 0.45309 ± 0.00134 vs theory 0.45167. 0 allocations/realization probe on five topologies; Criterion data below. See the E3 section. |
| Incremental component process (G1) | Validated | `IncrementalUnionWorkspace` (`percolation/process/incremental.rs`): generic incremental union substrate over any `UndirectedGraphView` — `activate_vertex`/`activate_edge` (union when both endpoints active) each return an online `IncrementalSummary` (component count, largest size with canonical lowest-index identity, raw second moment, active vertex/edge counts) and append `(step, summary)` to a `ProcessTrajectory` whose budget `V + E` is exact. The step coordinate is caller-defined and never interpreted — the roadmap's "record process coordinates without forcing a specific law"; Achlioptas (G2) and invasion (G3) build on this substrate without a common law trait. Gates: random-activation prefix parity vs fresh F3 `analyze` on chain/square/random/multigraph; exhaustive activation orderings vs an independent flood-fill reference; canonical identity vs F3 labels; random vertex replay reproduces `scan_site` prefix summaries (test-only Z1 demonstration); reuse/reset/re-prepare cycles; asserted 0-allocation probe on five topologies; Criterion unions/s data below. |
| Achlioptas processes | Validated | G2 `AchlioptasProcessMC` (`percolation/process/achlioptas.rs`, roadmap §13): rule-driven competitive edge addition on the G1 `IncrementalUnionWorkspace` — every step proposes `m` candidate edges, scores them by the rule through `component_size` (the substrate's designated scoring primitive), and adds the minimum-score winner. ALL physical parameters explicit: `AchlioptasRule` (Product `s_u*s_v` / Sum `s_u+s_v`, evaluated as u128 — a `score(size_a, size_b)` closure is deliberately insufficient), `CandidateSampling` (WithoutReplacement partial Fisher-Yates pool / WithReplacement iid draws), `EdgeUniverse` (all physical edges or a validated subset), `TiePolicy` (tie = equal score in one step → LowestEdgeIndex default / EarliestDraw), `AchlioptasGrid` (Full/Every/At checkpoints over the occupation-count axis). Null steps under with-replacement are recorded, never filtered (occupied candidates compete and can win: draws consumed, nothing added). Bitwise RNG draw contract documented incl. the winner compaction; `step_with_candidates` replays a fixed candidate trace with no RNG. Runs to complete occupation (terminal = every permitted edge). Records ONE C0 vector sample (`AchlioptasLargestSize` + `AchlioptasOccupiedEdges` coordinates) + scalars (`AchlioptasGiantReached`/`AchlioptasStepsToGiant` — explicit criterion `2*S_max > V` after an acceptance — `AchlioptasSteps`, `AchlioptasNullSteps`). Gates: contract-driver replay reproduces sweeps bit-for-bit + pure-trace re-replay; EVERY step vs independent flood-fill on chain/square/random/multigraph across both rules, both sampling modes, 3 seeds; m=1 without replacement == test-side uniform permutation bitwise AND statistically == `scan_bond` (|z| < 1.2 at four occupation counts); m=1 with replacement == coupon-collector statistics (z = -0.87) and permutation-equivalence; terminal-state exactness (full + subset universes); hand-worked product-vs-sum split; asserted 0-allocation probe on five topologies; Riordan-Warnke (Science 335, 1189 (2012), arXiv:1102.5306, archived) continuity anchor. See the G2 section. |
| Invasion percolation | Integrator-verified | G3 `invasion_run` + `InvasionWorkspace` + `InvasionPercolationMC` (`percolation/process/invasion.rs`, roadmap §13): frozen per-realization threshold fields (site: one `f64` per vertex; bond: one per physical edge, NaN/infinity typed rejections) grown from explicit duplicate-free seeds through a **binary-heap minimum-threshold frontier** (`total_cmp` then lowest dense ID — the explicit tie policy; `-0.0 < +0.0` are distinct keys, not a tie). Explicit stop policies (`GrowAll`, `FirstSpanning{from,to}` checked after planting and after every accepted step, `Steps(n)` with zero legal) and an optional Wilkinson-Willemsen trapping rule (explicit exit set; invadable iff connected to an un-invaded exit through the defending phase; trapped sites retire permanently — monotone, DFS with early exit). Site mode queues a site once carrying its first invader's edge; bond mode queues boundary edges and discards stale pops; parallel edges are independent entries, self-loops never enter. Every accepted operation drives the G1 substrate with invasion times as the step coordinate. The invasion graph is a forest (one tree per seed) by construction — verified. With UNIQUE thresholds, non-trapping bond invasion accepts exactly the Prim growth order (single seed: Prim from the seed, terminal tree = the unique MST of the seed's component; `k` seeds: Prim from a virtual super-root joined to the seeds; `Steps(n)` = the exact prefix) — verified against independent test-side Prim AND Kruskal. Records ONE C0 vector sample (`InvasionThresholdOrder`, fixed shape, final-threshold padding) + scalars (`InvasionAcceptedTotal`/`InvasionInvadedVertices`/`InvasionMaxAcceptedThreshold`). Gates: hand-worked multigraph (parallel + self-loop) sequence/frontier/tree in both modes against a full-frontier-scan reference; 576+-realization heap-vs-scan parity (continuous and tie-dense thresholds, multiseed, all stop policies, both components graphs, trapping); forest invariants via independent union-find replay; BFS-verified spanning steps; bitwise adapter draw contract; asserted 0-allocation probe on five topologies; Criterion steps/s below. Wilkinson-Willemsen (1983) predates arXiv (Tarjan precedent — no archivable PDF; the trapping rule and the Prim/MST identity are verified against independent in-repo implementations, anchored by the archived Sweeney & Middleton PRE 88, 032129 (2013) statement of the correspondence). See the G3 section. |
| Particle candidate-pair adapter (P1) | Validated | `percolation/continuum.rs` borrows the particle module's coordinates, `OrthorhombicCell` minimum-image geometry, and `CellList` stencil search (nothing duplicated, nothing modified there) to enumerate candidate pairs `(i, j)`, `i < j`, boundary-inclusive `|displacement| <= cutoff`, in a documented cell-major canonical order. The output is an ephemeral owned `PairGraph` (vertex count + pair `Vec`, `UndirectedGraphView` with derived incidences) that exists for one F3 analysis pass and is dropped — the deliberate transient counterpart to the roadmap's ban on a persistent per-snapshot CSR. Zero-copy `candidate_pairs_from_configuration` and validated prepared-list `candidate_pairs_with_cell_list` entry points serve equilibrium-simulation callers. Gates: exhaustive all-pairs parity in 2D/3D, exact-boundary inclusivity, periodic-image pairs, determinism, lattice and generic translation invariance, F3-vs-independent-BFS summary and label parity, typed rejections, throughput record; Criterion data below. |
| Continuum/Boolean percolation | Validated | `percolation/boolean.rs` (P2): overlapping disks/spheres in orthorhombic periodic cells with `Radii` (monodisperse or finite-polydisperse, typed validation) and a `PointProcess` sampler (binomial fixed count; Poisson with a documented chunked-Knuth draw contract) plus the equilibrium **measurement** path — caller-supplied coordinates analyzed by `overlap_pairs`/`overlap_pairs_from_configuration` or a frozen `BooleanPercolationMC::measuring` adapter that never resamples. Overlap pairs compose P1 at the maximum contact cutoff `2 r_max` with an in-place per-pair filter `d^2 <= (r_i + r_j)^2` (inclusive, squared, no epsilon; monodisperse radii pass through unfiltered — bit-identical criterion), feeding the F3 analyzer through an ephemeral `OverlapGraph`. Gates: all-pairs overlap parity 2D/3D polydisperse, translation/permutation/periodic-representative invariance, zero-radius and fully-overlapping limits, Poisson count mean/variance CI gates (direct and through the adapter), frozen-configuration zero-variance remeasurement, bitwise per-seed determinism, and the 2D disk covered-fraction threshold bracket `eta_c = 1.12808737(6)` (Mertens & Moore 2012, archived) crossing 1/2 inside `(eta_c/1.15, 1.15 eta_c)` — a statistical bracket, never exact. |
| Continuum winding analysis (T4) | Validated | `percolation/continuum_winding.rs`: `ContinuumEmbedding::from_pairs` derives one integer `Z^D` cell vector per P1 `PairGraph` / P2 `OverlapGraph` physical edge — `n_c = round((d_c - r_c)/L_c)` from the particle module's minimum-image displacement `d` and the raw in-buffer delta `r`, under a documented tolerance (dimensionless `1e-9`; fractional, non-finite, or `> 2e9` ratios are the typed `NonIntegerCellOffset`, edge and axis named — never silent rounding or saturation) — and implements `PeriodicEmbedding`, so the unchanged W2 `analyze_winding` accepts continuum configurations; reverse queries negate exactly through the lazy `DirectedDisplacement`. Load-bearing typed precondition, the offending axis named: `cutoff < L_c/2` strictly on every axis (below it the minimum image is the unique physical connection and `n` is its boundary count; beyond it a pair can be within cutoff through two images, making integer winding well-defined but physically ambiguous — the review demonstrated verdict flips under 1e-15 noise at the tie). Gates: exact negation + independent recomputation over 515,079 edges / 1,222,184 axis entries in 2D and 3D (0 disagreements); hand-derived 1D/2D tables; the P2 physics gate — 80 Poisson Boolean-model configurations (20 seeds x spanning/sub x {2D eta 1.5/0.7, 3D eta 0.6/0.15}), 259,068 particles, 50,484 components, winding rank and wraps agreeing EXACTLY component-by-component with a test-local BFS unfolding; winding and extent-spanning census-counted as independent in both directions (0 events at these densities); precondition boundary (P1 accepts the `L/2` tie, T4 rejects it, just-under works end to end); public-construction tolerance rejections; 0-allocation steady-state probe with per-realization construction accounting. Independent review passed 2026-10-23 (a third arithmetic convention over 1,156 fresh configs / 1.8M edges: zero disagreements; gauge invariance over 241k checks; 47 fresh-seed configs / 10,586 components: zero verdict disagreements; the review's Major — an overstrong extent-implication assertion in a gate — corrected as follow-up by census-counting both independence directions and pinning the hand-ring counterexample in a regression test). Discharges the v1 W-track note; see `docs/reviews/`. |
| External bond-configuration analysis (R1) | Validated | Decoupled analysis contract of the R track: an `ExternalBondConfiguration` borrows a caller-supplied bond mask (one `bool` per physical edge) and feeds the shared F3 `analyze`/`analyze_with_labels` and W2 `analyze_winding` as both activity sources (all vertices active, edges per mask), with the `analyze_bonds` one-call convenience and an owned `BondSample` open-edge record for future sampler exposure. Analysis is deliberately sampler-agnostic; gates and bench below. |
| FK/random-cluster q=1 | Validated | R2 validates the exact q = 1 collapse (the FK measure `p^{\|A\|} (1-p)^{E-\|A\|} q^{k(A)}/Z` becomes the product Bernoulli measure because `1^{k(A)} = 1` and `Z(p,1) = 1`) configuration by configuration: exhaustive 2^E enumeration on a 5-edge chain-ish graph and a 6-edge self-loop/parallel multigraph (uniform and heterogeneous p, per-configuration weight equality, Z = 1, observable-mean identity), bitwise delegation identity to `BondBernoulli`, and multi-seed z-scores of `RandomClusterMC` (Bernoulli delegate → R1 external-mask analysis) versus the `StaticPercolationMC` bond path at matched p on a square lattice. `q ≠ 1` is a typed roadmap-citing rejection. Delivered in-branch on `feat/cmc-percolation`. |
| FK q=2 Edwards-Sokal/Swendsen-Wang | Validated | R3 `EdwardsSokalMC` (roadmap §15): an owned-`CsrLattice` Swendsen-Wang/Edwards-Sokal adapter whose `sweep()` performs the conditional FK bond draw `p_e = 1 - exp(-2 beta J_e)` on same-spin edges (ferromagnetic `J_e > 0` typed-rejected otherwise) and exposes that ACTUAL draw as an R1 `BondSample`; the following `measure()` runs the shared F3 analyzer on that same sample through the F4 `ObservablePlan`. Frozen-spins mode reduces the sweep to the pure conditional draw for exact enumeration. Gates: frozen-spin bond-subset enumeration vs the product Bernoulli law with bitwise RNG-contract replay; partition-function agreement `Z_Ising = exp(beta sum J_e) * sum_s sum_A p^{|A|}(1-p)^{E-|A|}` (exhaustive over spins and bond subsets, incl. a parallel+self-loop multigraph, own union-find `k(A)`, independent brute-force Ising reference); `P_FK(i ~ j) = <s_i s_j>` on a proper ES chain (z-score gates, both sides vs exact and vs each other); worm-representation (even/two-defect subgraph, tanh weights) cross-checks at the observable level only — the worm's even subgraph is not FK and is never relabeled. Edwards-Sokal (1988)/Swendsen-Wang (1987) pre-date arXiv (no archivable PDFs; Tarjan precedent) — the exact double-sum identities are the anchors. Delivered in-branch. |
| Generic real-q random-cluster model | Validated (Q1/Q2); Q3 Experimental | Q1 Sweeny samples arbitrary real q>0 (exact-FK marginals, q=1 bitwise, rate-level detailed balance); Q2 Chayes-Machta covers integer q with exact joint stationarity; the q phase diagram is Q3 (exploratory). See the Q rows above. |
| Network edge-list input | Validated | `csr_from_edge_list` builds an owned `CsrLattice` from text edge lists (`left right` per line, `#` comments and blank lines skipped) under an explicit `EdgeListSpec` (declared vertex count, duplicate policy Reject/AllowParallel, self-loop policy Reject/Allow); every rejection is typed and carries the 1-based line (parse, endpoint range, duplicate, self-loop). Zero-copy external arrays remain the F1 `BorrowedUndirectedCsr`. No binary/memory-mapped format: million-edge ingestion evidence recorded below. Independent review passed 2026-10-22 (see `docs/reviews/`). |
| Boundary crossing between explicit vertex sets | Validated | F3 `BoundaryQuery` is validated for owned/borrowed generic undirected views, supports several queries in one component pass, and defines empty/overlap/dedup semantics; F4 migrates the full crossing test assets. The square left/right-column and chain endpoint defaults live in the `StaticPercolationMC` parameter schema; all other parameter-built topologies require explicit site sets. |
| Periodic embedding capability (W1) | Experimental | Directed-incidence integer cell displacements exist as explicit recorded data: the `PeriodicEmbedding` capability (`topology/embedding.rs`) and the `LatticeEmbedding` table serve the pbc chain/square/hypercubic builders and the triangular builder, which record each displacement at bond emission and never infer it from `BondType` or vertex IDs. Independently reviewed; see the W1 section below. |
| Periodic wrapping/winding analysis (W2) | Validated | Winding-capable component analysis over the W1 displacements: `analyze_winding` (percolation/analyzer/winding.rs) runs one physical-edge pass of union-find with checked `i32` `Z^D` potentials, reports per-component winding rank and an echelon generating set (never only a Boolean), and feeds `Wrapping`/`WindingRank` plan observables on `StaticPercolationMC::with_embedding`; open/non-embedded wrapping requests are typed errors. Engineering gates green and independent review passed 2026-10-22 (two from-scratch brute forces reproduced 187/256 and 186684/262144; potential and echelon algebra hand-verified; exact anchors now asserted); see the W2 section below and `docs/reviews/`. Crossing is not wrapping. |
| Honeycomb/kagome periodic embeddings (T3) | Validated | Roadmap v2 §2 T3 W-track tail: `build_honeycomb_with_embedding` / `build_kagome_with_embedding` (pbc only) restructure both builders to a shared bond-emission pass (`honeycomb_bonds`/`kagome_bonds`) so plain and embedding variants share one construction and cannot drift — plain-builder parity pinned bit-identical in-repo across sizes per family, and pre-refactor parity verified externally by the review — recording each bond's fundamental-cell displacement at emission, including the honeycomb `Lx = 2` doubled horizontal pairs carrying `(0,0)` and `(-1,0)`; kagome never has parallel bonds (algebraic: a repeated neighbour needs `lx|1` or `ly|1`; its `Lx=2` obstruction was a covering-image attribution ambiguity, not multiplicity). Gates: hand-derived 2x2 tables (6 + 24 bonds, asserted in both directions), spot checks at 4x2/3x2/4x3, independent neighbour-family re-derivation over size sweeps (exact negation, `is_zero`, per-pair displacement multiplicity), bond-content closed forms (doubled pairs, wrapping-bond counts), and W2 winding integration — fully occupied tori rank 2 with determinant-±1 generators matching the lift-BFS reference; wrapping-bond cuts never wind. Site numbering documented as a public contract in the builder rustdoc. The delivered gate set subsumes the roadmap's W1-style plaquette/cocycle gates (the winding-rank and generator checks over the hand tables are strictly stronger). Independent review passed 2026-10-23 (hand tables re-derived from the rustdoc alone before reading the tests; an own generator exact at all 25 swept sizes; the pre-T3 builders run verbatim via a path dependency — bit-identical; adversarial winding configurations); see the T3 section and `docs/reviews/`. |
| Newman-Ziff parameter scans (Z1) | Validated | `scan_site`/`scan_bond` activate a uniformly random permutation and record the full occupation-count curve (component count, largest size, raw second moment, active edges); `reweight` turns one curve into scalar estimates at uniform `p` via binomial weights, and `ScanEstimator` reports scan-to-scan uncertainty. Zero allocations after `prepare`. Uniform site/bond only; no heterogeneous/mixed/deletion scans, no scheduler adapter. Independent review passed 2026-10-22 (reweight recurrence verified against the archived paper to 9.5e-14 at n=10^6; from-scratch flood-fill parity; speedup reproduced); see the Z1 section below and `docs/reviews/`. |
| Heterogeneous Newman-Ziff reweighting (T2) | Validated | `reweight_heterogeneous` reweights one unchanged Z1 scan curve to a per-element probability field (Poisson-binomial weights, all log-space with a running-max rescale): the uniform-permutation prefix identity makes `sum_n C(N,n) w(A_n) f_n` an unbiased per-scan estimator and `sum_n C(N,n) w(A_n)` a unit-mass diagnostic with expectation exactly 1; uniform fields delegate to the Z1 binomial `reweight`. Gates: exact 2^N enumeration vs a test-local BFS reference (site/bond, U(0.3,0.8)-spread and bimodal 0.2/0.9 fields, max z 1.35 at R=3200); unit-mass and finiteness at N=10^4; L=64/128 graded-ramp cross-validation vs direct `SiteBernoulli` (max z 1.21); 1e-9 uniform reduction (pinned to <1e-8 at N=10^4); typed boundary/length rejections; asserted 0-allocation probe on five topologies. The estimator is heavy-tailed outside a documented envelope (`sqrt((N/4) Var(logit p))` of order one or below) and the 8-field Criterion speedup is below 1 on dense topologies (0.38-0.80x; 1.40x on sparse-er site) — spread fields belong to the N1 direct path; see the T2 section under Z1. Independent review passed 2026-10-22 (identity brute-forced over all permutations in exact rational arithmetic; the normalized-ratio bias independently confirmed and v1 shown structurally immune; a fully independent sampler/union-find cross-check at L=64 gave z=-0.18/-0.26; saturation documentation corrected as follow-up); see `docs/reviews/`. |
| Weighted shortest-path layer (v2 F1) | Validated | Roadmap v2 §3 F1: `first_passage_run` + `FirstPassageWorkspace` (`percolation/first_passage.rs`) — single-source single-pass Dijkstra with a binary heap over any `UndirectedGraphView`, edge weights a caller-owned slice of one `f64` per physical edge in dense edge order (`lattice.edges[e].weight` for a `CsrLattice`; weights are data, redrawn per call for F2's passage-time fields). Explicit policies: heap keyed `(distance total_cmp, vertex id ascending)` — the id tie key makes the pop order total under ties (equal-distance settles and equal-cost predecessor choices are deterministic, not heap-implementation detail; the same tie discipline as the invasion frontier, and a 16-byte heap entry); zero weights ALLOWED (non-negativity is Dijkstra's requirement; zero-weight cycles correct); parallel edges are independent relaxations with the lighter bond winning (strict-improvement updates); unreachable = `f64::INFINITY` with reachability the derived `is_finite()`; optional predecessor vertex+edge arrays (`GeodesicTracking`, for F2/F3 geodesics) gated off by default. Typed rejections at validation time, before the workspace is touched: negative weight (edge named), NaN/±inf weight, wrong weight length, source domain, empty graph. Gates: parity vs a test-local independent Bellman-Ford over chains/cycles/trees/dense-ER multigraphs (parallel + zero-weight edges, U[0,1) and exact-zero fields) — 249 (graph, source) cases / 3440 vertex distances, bit-exact on unique-path and dyadic-weight families, 1e-12 relative elsewhere, 0 mismatches; hand-derived 9-vertex 11-edge multigraph (parallel pair lighter wins, zero-weight triangle, equal-cost tie observable through the predecessor, self-loop, disconnected sentinel) plus zero-weight ring and single-vertex cases; rejection-with-untouched-workspace; reused-vs-fresh workspace bitwise determinism (48 comparisons) with growth/shrink and same-size graph alternation; asserted 0-allocation 0-byte probe per single-source pass on six topologies incl. E = 64 V. Criterion (V=1e5/E=3e5 ER): early-exit Bellman-Ford 2.27x FASTER than the lazy-deletion Dijkstra on the small-diameter ER (recorded honestly, T2-style) while Dijkstra is 5.81x faster on a diameter-398 open 200x200 grid with random weights — the F2-relevant regime; doubling V at fixed degree 6 costs 2.6x while the working set fits cache, 3.2-3.3x at V=4e5 where the random-graph heap leaves L3 (throughput 3.35 M vertices/s at V=2.5e4 down to 0.78 M at V=4e5). Algorithm layer only — no MonteCarlo adapter, no scheduler integration (F2 builds observables on top); the Auffinger-Damron-Hanson FPP anchor applies from F2 on. Independent review passed 2026-10-23 (docs/reviews/2026-10-23-f1-review.md; 712 fresh randomized + 55 adversarial cases vs a deliberately opposite-tie (LIFO insertion-counter) own Dijkstra plus own SPFA and rotated-rounds Bellman-Ford: bitwise agreement everywhere — stronger than the committed 1e-12 policy; the tie test non-vacuous: 2,545 vertices took different-but-valid predecessors at identical distances; E+1 heap budget tight and never exceeded incl. adversarial chains; saturation, rejection isolation, and the bench trade-off independently reproduced — grid ratio re-measured 3.68x vs the committed 5.81x, same direction, both recorded). Follow-ups applied: the GeodesicTracking::Off rustdoc now states prepare clears predecessors unconditionally (the 'stale' wording was wrong); the v1/v2 F1 name collision disambiguated. |
| First-passage observables (v2 F2) | Validated | `percolation/first_passage_observables.rs` (roadmap v2 §3 F2) over the v2 F1 pass: `PassageGeometry::corridor(length, width)` (open (r+1) x w Z^2 box, source set = left column, target set = right column) and `PassageGeometry::segment(steps, margin, half_height)` (single origin, explicit target points along +x inside a reflection-symmetric box), each carrying one super-source vertex joined to the source set by structural zero-weight edges — ONE `first_passage_run` pass realizes the min over the source set exactly. `EdgeWeightLaw` (Exponential `w = -ln(1-u)/rate`, written `0.0 - ln(1-u)/rate` so `u = 0` yields `+0.0` never `-0.0`; Uniform `w = u`) under a bitwise draw contract: one `random::<f64>()` per box edge in dense `build_square` order, structural zeros consume no draws. `PassageSampler` (reusable weight slice + F1 workspace + profile buffer; `realize_with_weights` runs production on caller-owned fields — the exact-enumeration path), `PassageEstimator`/`PassageEstimate` (scalar mean/unbiased variance/SEM of `T_r` = min over the profile plus per-target profile mean/variance). One realization = ONE C0 vector sample: the full passage-time profile along the target boundary (w entries corridor, `steps.len()` segment); the scalar min is a derived accessor. Gates: exhaustive 2^E two-level dyadic {0.5, 1.0} enumeration vs a test-local Bellman-Ford over test-built topologies — 263,168 corridor assignments (919,552 profile entries bit-exact, incl. 128,272 tied-minimum assignments where the optimal path is non-unique: the profile entries are unique minima over path weights, degeneracy never changes content) plus 4,198,432 segment assignments, 0 mismatches; bitwise draw-contract replay (474 weight entries, both laws, both geometries); same-seed bitwise determinism incl. reseed-replays-chain-one on a history-laden sampler; SEM ~ 1/sqrt(N) (ratio 0.531, N = 128 vs 512). Physics brackets (ADH 2017 anchor, bracket language only): Var(T_r) growth on aspect-1 corridors (w = r, r in [8, 128], R = 600 per law) gives fitted chi = 0.2476 CI (0.2028, 0.2923) for Exp(1) and 0.2658 CI (0.2198, 0.3118) for U[0,1), both inside the rigorous open interval (0, 1/2) — never claimed against 1/3 (the effective exponent below 1/3 at these sizes is the documented small-r crossover); subadditivity `E[T_2r] <= 2 E[T_r]`, exact in expectation on reflection-symmetric segment boxes, holds with z = -10.41 / -13.57 / -15.10 / -10.35 at r = 8/16/32/64; the time-constant sequence mu(r) = E[T_r]/r = 0.5621(47) -> 0.5108(27) -> 0.4743(16) -> 0.4488(10) -> 0.4341(6) decreases and stabilizes (convergence from above, by subadditivity). Asserted 0-allocation 0-byte probe per realization (law draw + F1 pass + estimator add) on 64x16 and 128x32 corridors; 233 and 55 realizations/s at 64x64/128x128 aspect-1 corridors (debug test profile). No MonteCarlo adapter — F3's layer. Independent review passed 2026-10-24 (docs/reviews/2026-10-24-f2-review.md): the reviewer's own Floyd-Warshall over their own geometry and edge order ran 21.3M assignments INCLUDING the literal 4x4 corridor (2^24, 67.1M profile comparisons, 8.15M tied-minima, 0 mismatches) and non-dyadic weight families; the draw contract reproduced against a from-scratch xoshiro256++ (10,000 draws bit-equal, -0.0 bits confirmed); chi reproduced within 0.13sigma with an independent pipeline; the reflection subadditivity argument given a coupling proof and pinned pointwise (0 violations in 40,000 fields, max rel diff 6.3e-16); statistics machinery verified (t-CI df=3, SEM slope -0.5008). Follow-up: construction-cost figures marked machine-dependent (both measurements recorded). |
| First-passage shape and geodesics (v2 F3) | Validated | `percolation/first_passage_shape.rs` (roadmap v2 §3 F3) over the v2 F1 pass; independent review pending. `ShapeGeometry::new(directions, distances, margin)` builds one open centered `side x side` Z^2 box (`side = 2*(max n + margin) + 1`) with targets `n*(1,0), n*(0,1), n*(1,1), n*(1,-1)` at every requested `n`; ONE F1 Dijkstra per realization yields the whole C0 passage-time profile (`ShapeSampler`; F2 bitwise draw contract one uniform per box edge in `build_square` dense order; `realize_tracked` adds the F1 predecessor tree), and `ShapeEstimator` accumulates per-target mean/unbiased variance/SEM. `geodesic_count_run` counts geodesics exactly (u128, explicit saturation) by DP over the tight-edge DAG (`dist[u] + w == dist[v]`) in increasing distance order; parallel edges count separately, zero weights are a typed rejection (a zero-cost tight cycle makes the count infinite), unreachable counts 0, `source == target` counts 1. `geodesic_path` reconstructs the deterministic F1 predecessor path into caller buffers; `geodesic_aging` is an EXPLORATORY finite-horizon persistence statistic (shared leading prefix of the geodesics to `n u` and `2 n u` in one environment), anchored loosely in the survey's Question 30 and explicitly NOT the Question 16 independent-overlap object. Anchor (archived ADH 2017): Theorem 2.1 time constant `T(0, n e_1)/n -> mu(e_1)` a.s./L1 with `lim E[T(0,n e_1)]/n = inf_n E[T(0,n e_1)]/n` (convergence from above; ~lines 392-433), `0 <= mu(e_1) <= E tau_e` with strict inequality for non-trivial `F` (Theorem 2.4; ~lines 518-526), Cox-Durrett Theorem 2.16 convex compact symmetric `B_nu` with `(1-eps)B_nu subset B(t)/t subset (1+eps)B_nu` eventually a.s. (~lines 745-755), Definition 4.1 geodesics and the uniquely-geodesic spanning tree (~lines 3904-4117), Question 16 independent-geodesic overlap (~lines 2280-2283), Question 30 edges in infinitely many geodesics (~lines 5044-5046). Gates: uniform weights reproduce the l1 ball exactly (36 vertex distances on 6x6 plus 12 directional targets = n*(abs(dx)+abs(dy))); uniform geodesic counts exact at corner and interior targets on 6x6/10x10/12x12 = C(10,5)=252, C(18,9)=48,620, C(22,11)=705,432 (interior 126/24,310/352,716); two-level dyadic {0.5,1.0} counts match brute-force simple-path enumeration on 16 graphs (3x3, 4x3) x 8 seeds (7 with non-unique geodesics, 34 total geodesics), all-0.5 counts C(4,2)=6 and C(5,3)=10, and every reconstructed path is simple with weight sum bitwise equal to the F1 distance; directional time constants at exponential rate 1, centered 97x97 box, margin 16, R=150: `mu(east) = 0.5653(108), 0.4978(71), 0.4750(39)`, `mu(north) = 0.5513(107), 0.5077(65), 0.4723(39)`, `mu(northeast) = 0.7637(107), 0.6995(61), 0.6551(44)`, `mu(southeast) = 0.7715(111), 0.6997(62), 0.6539(38)` for n = 8/16/32 (SEM in parentheses); every CI is positive and below `E tau_e = 1`, every consecutive ratio change is negative within 3 combined SEM, the 8->16 majorant exceeds the 16->32 decrement (flattening) and the latter stays inside 0.06, axis symmetry `abs(mu(east)-mu(north))` is within 3 noise, and the paired norm inequality gives `mu(east)+mu(north)-mu(northeast) = 0.29221 ± 0.00586` and `- mu(southeast) = 0.29345 ± 0.00599`, both CI-positive; NO closed-form shape is claimed. Exploratory aging at exponential rate 1, east, R=60: persistence 0.6551(427) / 0.7294(405) / 0.6368(445) at r = 8/16/24 (SEM), 360 exact geodesic path checks, drift between the two largest radii 0.0925 vs the generous 0.10 + 3*noise band — labelled exploratory, no asymptotic claim. Typed rejections (geometry 5, sampler 7, estimator 3, counting 6, path 3) with state isolation. No throughput/allocation probe added at this stage. |
| Long-range power-law edge law (L1) | Validated | `percolation/long_range.rs` (roadmap v2 §4 L1): Biskup's long-range edge law `p = min(1, beta * r^{-s})` (archived Biskup 2004 md; the explicit power-law member of `p_xy = |x-y|^{-s+o(1)}`) on a finite box of Z^d — `LongRangeBox` (per-dimension sides, row-major dense indices, ℓ1 distance: integer-exact, ℓ2-equivalent so Biskup's s-vs-d-vs-2d regime statements transfer, and equal to ℓ2 in d=1), the finite-box restriction of the infinite-measure law (pairs with BOTH endpoints in the box; boundary effects confined to the box surface). Consideration order lexicographic (x,y), x<y dense; RNG draw contract: ONE uniform f64 per pair in ascending order (saturated and near-zero pairs included), occupied iff u < min(1, beta r^-s) — `pair_probability` is the pinned bitwise expression; `min(1,.)` saturates beta>=1 nearest neighbors surely (u<1 always) and beta>=r^s shells; beta*r^-s never overflows past beta (r^-s in (0,1]). One realization = one G1 `IncrementalUnionWorkspace` construction over a zero-materialization complete pair view: all vertices activated at step coordinate 0 (G2 convention), each occupied pair activated at step = 1-based ADDED-EDGE COUNT (the process coordinate; realized count random); trajectory steps exactly [0;n]+[1..=m]. `long_range_run` RNG-free core + `LongRangePercolationMC` adapter recording ONE C0 vector sample (`LongRangeLargestSize`, fixed shape, G3 final-value padding) + scalars (AddedEdges/ComponentCount/SecondMoment). Typed rejections: empty box, zero side, negative/non-finite beta or s, wrong-length or out-of-[0,1] uniform fields, zero steps. Gates: bit-adjacent threshold replay from both sides (28 pairs x 3 combos) pinning law + pair order; per-pair frequency z-scores + pooled Pearson chi-square on 1x8 (3 laws) and 4x4 boxes, 600 realizations (max |z| 3.54; chi2 33.3/36.0/23.2/110.2 vs criticals 57.0/57.0/46.9/173.7) + joint-occupancy independence z in [-2.26, 1.96]; closed-form E[deg(center)] z-gates at s<1, s=1, s>2 (z = -0.41/-0.59/0.66, R=4000); bitwise same-seed replay + adapter-vs-manual-order parity; per-step parity vs independent union-find (773 steps) + final structure vs flood fill (30 configs); saturation regimes (r=1 always at beta>=1 across 3 laws x 40 seeds, 38 saturated pairs at beta=4/s=1, complete graph at beta=1e6, beta=0 empty, singleton box). Asserted 0-allocation probe; Criterion-less eprintln throughput below. See the L1 section. Independent review passed 2026-10-24 (docs/reviews/2026-10-24-l1-review.md): the law replayed bit-adjacent on 309 (pair, law) combos with the reviewer's own sampler matching every one of 100,000 realizations; corner/edge degree moments vs exact rationals (s=1 harmonic closed forms) within |z| <= 1.81 at R=5,000 — extending the committed center-only gate; 21 fresh independence pair-pairs exact-q z <= 2.45; overflow-safe rejection boundaries proven by induction; Biskup anchor quotes verified against the archive. The review's Major — LongRangeBox::coordinates() axis-reversed for d >= 2 (zero in-repo callers, off the realization path) — FIXED as a follow-up with the mixed-radix digit order corrected and a round-trip inverse gate added across 1D/2D/3D asymmetric boxes; the PairUniverse lookup comment corrected to O(log V). |
| Chemical-distance regimes (v2 L2) | Validated | `percolation/long_range_distances.rs` (roadmap v2 §4 L2) over one L1 realization: the chemical distance `D(x,y)` is the unweighted graph distance on the occupied edge set. `ChemicalDistanceWorkspace::profile()` holds one `f64` per box vertex, entry `y` = `min_{x in sources} D(x,y)` (multi-source BFS; `f64::INFINITY` sentinel, reachability is `is_finite`), and one realization is one fixed-length C0-shaped profile vector — never a scalar alone. Optional `DiameterTracking::AllVertexBfs` adds the exact box diameter = max finite distance over all pairs including the diagonal (isolates 0; on a disconnected realization it is the largest within-component distance), with `DIAMETER_VERTEX_CAP = 4096` as a hard guard-rail (`O(V^3)` complete-graph worst case; the profile path is uncapped). The occupied set is the L1 draw contract replayed in-layer (one uniform per pair ascending, `pair_probability` evaluated once per ℓ1 distance — the same pinned bits), materialized as CSR; `ChemicalDistanceSampler` owns the uniform field and workspace (zero steady-state allocations; adjacency capacity is reserved WITHOUT resizing, so pre-run accessors report an empty adjacency), and `ChemicalDistanceEstimator` accumulates per-entry mean/unbiased SEM/reachable fraction (per-entry conditioning on reachability — Biskup's `x,y in C_infinity` conditioning) plus diameter mean/SEM. Typed rejections: non-finite/negative beta or s, empty source set, source domain, wrong-length or out-of-[0,1] uniforms, diameter cap, profile length, zero-vertex estimator, and dimension overflow (`vertices+1`, `distance_range+1`, `2*pairs`). Regime gates (Biskup anchor, regime-correct): (1,2) typical distance at beta=1, n=2048, R=150 for s=1.3/1.5/1.8 — monotone growth, finite-window D-vs-r log-log slope 0.172/0.219/0.329 (CI below 0.5), fitted polylog exponents 0.658/0.838/1.258 all below Delta = 1.609/2.409/6.579, CI-disjoint and increasing in s; s=2.5 linear regime at n=8192, R=24 — the gate is PAIRED per realization (24 slopes and 24 `D/r = A + c/r` fits, 95% CIs): D = 105.6/202.1/404.7/747.2, slope 0.949 +/- 0.043 CI (0.907, 0.992), A = 0.3713 +/- 0.0199 CI (0.351, 0.391) inside (0,1), c = 11.18 +/- 8.55 CI (2.63, 19.73) strictly positive, ratio 0.413 -> 0.365 non-increasing; s=0.5 diameter is a BOUNDEDNESS gate (R=8: diameter exactly 3 at n=256..2048, explicitly reported as a constant sequence with undefined R²; means in [2,4], no positive power growth, and the a.s. sup `d/(d-s)=2` is a reference scale, NOT a measured equality); s=1 diameter 5.0/5.33/5.33/6.0 over n=512..4096 against the `ln N` = 6.2..8.3 upper envelope (slope 0.079 +/- 0.096, CI (-0.017, 0.175)), upper-envelope only, no lower claim; s=2 deliberately not run. Correctness gates: 35 profiles + 35 all-vertex diameters on 1D/2D boxes (incl. beta<1 disconnections) bitwise vs a test-local BFS/all-pairs reference, with the occupied set also equal to a test-local law replay and to the L1 core's `occupied_pairs`; estimator mean/SEM/reachable parity on 40 beta<1 realizations (min reachable fraction 0.40) plus the beta=0 never-reached conventions; same-seed bitwise profiles/diameter/adjacency and L1 draw-stream parity; pre-run empty-adjacency accessors. Asserted 0-allocation probe on three topologies plus diameter mode; debug-profile throughput 192 / 12.1 realizations/s at 1x256 / 1x1024. See the L2 section. Independent review passed 2026-10-24 (docs/reviews/2026-10-24-l2-review.md): the reviewer's own mixed-radix geometry/BFS/all-pairs machinery ran 1,188 exact realizations (21,549 checks, 0 failures) incl. beta<1 disconnections, duplicate sources, odd boxes and caps; the draw stream matched a fresh xoshiro stream and the L1 adapter bitwise; the estimator and rejection/precedence matrix were re-derived; the zero-allocation closure was inspected, not just asserted. The review's Majors — seed-fragile s=2.5 gate at R=12 (now paired at R=24 with CIs) and pre-run `occupied_edge_count()` returning `pairs` (now capacity-only reservation, with a pre-run accessor gate) — are FIXED as follow-ups; the never-reached-SE doc mismatch, public `prepare`/`add_diameter` overflow, capacity-failure wording, s=0.5 boundedness wording, "slower than any power"/"from below" overclaims, diameter-probe bytes assertion, and the s>2 `D_N` citation nuance are also applied. |
| Correlated activity fields (C1) | Validated | `percolation/correlated.rs` (roadmap v2 §5 C1): the Schrenk/Pruessner Fourier-filtering generator, quoted exactly (spectrum `S(f) ~ |f|^{-beta_c}`, `beta_c = 2(H+1)`, `c(r) ~ |r|^{2H}`, i.e. `C(r) ~ |r|^{-a}` with `a = -2H`; target `C(r) = (1+r^2)^{-a/2}`). A radix-2 Cooley-Tukey transform per power-of-two axis writes the fixed draw contract (mode order, Hermitian pairing, normalization, one or two uniforms per mode); the field has unit site variance, and occupation is the `Phi^{-1}(1-p)` threshold so the marginal is exactly p. `CorrelatedPercolationMC` reuses the unchanged F3 analyzer/ObservablePlan; `a = +infinity` is the white-noise/iid delegate (bitwise `SiteBernoulli`, the `a=inf` table row), while `a = 0` is the constant-field limit (all-or-nothing, no iid claim). Gates (unit): design correlation and unit variance in the ensemble, marginal p with super-binomial fraction fluctuations, adapter-vs-direct parity, bitwise determinism, `a=0` constant field, `a=inf` iid delegate, typed rejections. A binary (thresholded) correlation gate measures a=1.5 r=1/r=4 = 0.0993/0.0212, a=2.5 r=4 = 0.0049, and white noise r=1/r=4 = -0.0004/-0.0002: positive, decaying, exponent-ordered, and vanishing in the iid limit. Allocation probe: 0.000 allocations/bytes per sweep+measure on 16x16 and 1D-64. |
| Shifted thresholds (C2) | Validated | `percolation/correlated_thresholds.rs`: periodic `L x L` ladders, one threshold-sorted incremental union-find pass per replica, order-parameter susceptibility peak, bootstrap intervals, `p_c(inf)` extrapolation and log-derivative `nu` fit. The quantitative gate is an ignored long test using L=128/256/512 and 256 replicas (Pruessner's table uses L=64-8192, 10^5 replicas, and a WRAPPING threshold; our estimator is the order-parameter susceptibility peak `N Var(S_max/N)`, not their wrapping estimator): the peaks are 0.5139/0.5268/0.5344 and the `L^{-1/nu}` extrapolation with nu=4/3 gives `p_c(inf)=0.5457`, CI (0.5355, 0.5558), chi2/dof=0.00. The 0.0027 offset from 0.54299 is ~0.65 bootstrap sigma and the reference value lies inside the CI; the gate also keeps a +/-0.005 tolerance band. The largest-fraction crossing estimator is independently calibrated on white noise (0.5987 vs the exact 0.592746). C2 remains Experimental pending independent review. |
| Weinrib–Halperin map (C3) | Experimental — calibrated nu estimator; flip not certified | `percolation/weinrib_halperin.rs`: the extended Harris criterion `a nu_iid = 2`, `a_H = 2/nu_iid = 1.5` for `nu_iid = 4/3`, and the two branches `nu = 4/3` (a >= 1.5) and `nu = 2/a` (a < 1.5), with bracketing helpers. Gates: the analytic branch/nu values at a = 2, 1.5, 1.25, 1, 0.75, 0.5 and typed rejections. The old log-derivative nu estimator is diagnosed as non-scale-covariant (its max is taken on a noisy flank with a one-grid-step difference): on the exact iid control it returns 1.215 (0.739, 1.690) with chi2/dof 11.7, so the old 0.934/1.020/0.804 values were estimator failure, not physics. The replacement (`ThresholdLadder::nu_from_pc_shift` on `ThresholdLadder::pool`) refines each order-parameter susceptibility peak with a local log-quadratic fit over `PEAK_FIT_RADIUS = 0.004` and fits `p_c(inf)-p_c(L) = A L^{-1/nu}` with `p_c(inf)` fixed. Pooled 8x512 replicas at L=128/256/512/1024: iid control nu = 1.379 (1.325, 1.433), chi2/dof 0.11, contains 4/3; a=1.5 1.592 (1.413, 1.771); a=1.25 1.819 (1.542, 2.095); a=1.0 2.529 (1.741, 3.318). The ordering is right but the values are prior-dominated (`d nu/d p_c(inf) ~ 0.12` per 0.001) and the 1.5-vs-1.25 gap is 0.226 +/- 0.208 (1.1 sigma), so the branch flip is NOT certified. A one-off larger-L cross-check (L=256/512/1024/2048, 4 pooled ladders x 128 replicas) gives a=1.5 nu=1.620 (CI 1.356-1.884, chi2/dof 1.76) and iid control nu=1.334-1.512 depending on the fit radius; the flip is still not certified. The nightly control gate asserts only the iid calibration. |
| Sweeny sampler (Q1) | Validated | `percolation/sweeny.rs` (roadmap v2 §6 Q1): single-bond heat-bath/Metropolis dynamics for arbitrary real q>0, with the Geyer Eq. 4/5 acceptance, the `1/2 <= pHB/pMH < 1` bound, and the unconditional accept/reject shortcut. `Delta n` queries use a spanning forest with union-by-size and an interleaved smaller-side replacement search; the worst-case bridge cost is O(V) and the docs say it is not the asymptotically optimal dynamic-connectivity structure. Gates: exact 2^E FK enumeration on 3x3 for q=0.5, 0.75, 1, 2, 3 with own union-find and full Z — empirical |A| and k(A) marginals within |z| <= 1.36 (the q<1 cases caught and fixed the unconditional-shortcut ordering bug, whose two acceptance values are now explicitly sorted); q=1 one all-closed heat-bath sweep is bitwise `BondBernoulli`; update fuzz vs full recomputation; the pHB/pMH bound; typed rejections; and a rate-level detailed-balance audit over all 4096 masks x 12 edges x q=1,2,3, where the heat-bath and Metropolis rates equal the exact FK two-state conditional to 1e-12. |
| Chayes–Machta for integer q (Q2) | Validated | `percolation/chayes_machta.rs`: integer-q activate/sleep with per-cluster colour resampling and awake-color bond updates, one uniform per awake edge, deterministic sweep order. Gates: on a periodic 3x3 at p=0.45, q=2/3 agree with Q1's Sweeny chain within 4 sigma (q=2 diff 0.075 vs noise 0.173; q=3 diff 0.0007 vs noise 0.164), the |A| histogram against the exact FK distribution has max |z| 2.96/3.74, colour marginals have max |z| 1.48/2.93, and at q=2 the CM chain agrees with R3's Edwards–Sokal chain (diff 0.041 vs noise 0.191). Exact joint stationarity is now a gate: on a periodic 2x2 the full (bond-mask, colouring) space has 706 states at q=2 and 1458 at q=3, and one CM step preserves the Edwards-Sokal measure to max |mu P - mu| <= 1.3e-17 for every active colour. Allocation probe: 0.000 allocations/bytes per sweep+measure at q=2 and q=0.5. |
| q phase diagram (Q3) | Experimental — exploratory | `percolation/potts_phase.rs`: the exact 2D FK/Potts self-dual curve `p_c(q) = sqrt(q)/(1+sqrt(q))` (Fortuin-Kasteleyn; Wu, Rev. Mod. Phys. 54, 235 (1982)) is the anchor and is never fitted. `phase_scan` runs Q1 at p_c(q) on a periodic square and reports bond-density mean/variance plus the descriptive bimodality coefficient `(skew^2+1)/kurt`; q=4 is the marginal case, q>4 is first order. Gates: exact p_c at q=1,2,4,8; the scan machinery at q=2/8 reports finite density statistics (L=8: bimodality 0.362 vs 0.505); and an ignored nightly L=16 scan gives bimodality 0.363 (q=2), 0.400 (q=4), 0.519 (q=6), 0.654 (q=8) — below 5/9 for q<=4, above for q=8, with q=6 unresolved at L=16 and monotone ordering in q. No latent-heat or finite-size-scaling claim at this stage. |

| Diameter transitions (v2 L3) | Validated | `percolation/long_range_diameter.rs` (roadmap v2 §4 L3): a thin ladder layer over L2's exact `ChemicalDistanceSampler::realize_with_diameter` — no second diameter implementation. `DiameterRegime::classify(s)` encodes the d=1 Biskup boundaries (`Bounded` s<1, `Logarithmic` s=1, `Polylogarithmic` 1<s<2, `Critical` s=2, `Linear` s>2) and `supports_diameter_gate()` is FALSE for the critical case; `measure_diameter_ladder(sizes, exponents, beta, realizations, seed)` records every realization's exact diameter plus mean/SEM, with the documented row seed contract (`seed + row + i`). Gates: 7 tests — ladder replay bitwise vs direct L2 draws; s=2.5 over N=128..2048, R=24: D_N/N ratios 0.445 -> 0.400 -> 0.429 -> 0.387 -> 0.367 (non-increasing within 2 combined SEM; the middle step is noise), log-log slope 0.939 +/- 0.021 CI (0.871, 1.007), `D_N/N = A + c/N` extrapolation A = 0.3805 +/- 0.0153 CI (0.332, 0.429) strictly inside (0,1); s=1.5 over N=128..1024, R=12: D_N = 9.92/10.92/12.08/13.83, log-log slope 0.159 +/- 0.009 with upper CI 0.199 — sub-polynomial growth only, no sharp exponent (Biskup states the D_N bounds are not sharp in this regime); s=1 over N=512..2048, R=4: D_N = 5.0/5.0/5.5 against ln N = 6.24/6.93/7.63 with ratio to ln N/ln ln N = 1.47/1.40/1.47, upper envelope and non-decreasing within one integer step; s=0.5: D_N = 3.0 at every size, bounded envelope [2,4] and no growth; s=2 deliberately not run. Typed rejections: empty ladders, zero size, size above L2's cap, fewer than two realizations, invalid exponent, geometry and L2 validation errors. See the L3 section. |

At F0 no family was **Validated**. Since F4 passed independent review, the
static family and boundary crossing were the first Validated rows; N1, W2,
Z1, and N2 have since passed independent review on the same
production-conditioned basis.
Their production availability is conditioned on the branch merging into `dev`.

## Scientific definitions

Let one configuration induce tracked connected components
`C = {C_1, ..., C_K}` with sizes `s_i = |C_i|`. In site and mixed percolation,
tracked vertices are occupied sites. In bond percolation every vertex is
tracked, including isolated vertices. An edge connects in site mode when both
endpoints are occupied, in bond mode when the edge is occupied, and in mixed
mode only when the edge and both endpoints are occupied.

### Components and moments

- **Component size:** `s_i` is the number of tracked vertices in component
  `C_i`; edges are never counted as component size.
- **Largest component:** `S_max = max_i s_i`, with `S_max = 0` when `K = 0`.
  The current reference reports only its size, not identity.
- **Giant-component fraction:** `G = S_max / V_norm`. The normalization must be
  named. The platform default is total topology vertices `V_norm = |V|`;
  normalization by active/tracked vertices is a distinct observable. This is
  defined vocabulary, not a current PR #4 output.
- **Raw cluster moment:** `M_k = sum_i s_i^k` over the stated component set,
  with no normalization. PR #4 reports `M_2` as `SecondMoment` and includes the
  largest component.
- **Finite-cluster susceptibility:** choose one largest component per sample,
  breaking equal-size ties by lowest stable component/root ID, and remove it.
  Then `chi_f = (sum_{i != max} s_i^2) / (sum_{i != max} s_i)`. The exclusion
  happens before ensemble averaging. If the denominator is zero, the observable
  is undefined and must be represented explicitly as missing/invalid, not
  silently as zero or NaN in an advertised scalar result. Removing every tied
  largest component or dividing by total `|V|` are different named policies.
  PR #4 does not compute `chi_f` and its aggregate outputs cannot reconstruct it.

### Connectivity events

- **Boundary crossing:** for explicit vertex sets `A` and `B`, crossing is true
  iff one tracked connected component intersects both sets. Sets may overlap;
  an active vertex in `A intersect B` therefore satisfies the event. The event
  is defined by the supplied sets, not inferred geometry. Its ensemble mean is
  the crossing probability.
- **Wrapping/winding:** on a periodic embedded topology, lift each occupied
  path using explicit integer cell displacements. A component wraps when a
  cycle has nonzero accumulated displacement; winding records generators,
  axes, or rank rather than merely boundary-set connectivity. Vertex IDs,
  `BondType`, and periodic adjacency alone do not define winding. W2
  implements this analysis (Experimental, awaiting review; see the W2
  section).

### Network conditioning

- **Original degree** is degree in the immutable input graph before occupation
  or removal. An observable conditioned on original degree keeps this label for
  the complete realization or attack trajectory.
- **Residual degree** is degree in the currently active induced graph and may
  change during a process. Residual-degree conditioning or adaptive attack is a
  different process and must be named explicitly.
- Empty degree bins require an explicit missing/undefined representation; they
  are not zero-valued observations.

### Samples and lifecycle

A **realization** is one complete draw or process trajectory under a stated law.
An **ensemble sample** is one realization contributed to an estimator. For
static Bernoulli percolation these coincide: one `sweep()` redraws all relevant
occupancies independently and one following `measure()` analyzes that draw.
Samples are i.i.d.; set `thermalization_sweeps = 0`. For attacks, cascades,
growth, or directed space-time processes, all steps of one trajectory remain
one realization and must not be misreported as independent scalar samples.

## PR #4 behavioral baseline

PR [#4](https://github.com/UynajGI/scamper/pull/4) is open against `dev`; the
recorded implementation head is `f6031bf88bea05f61ecbc2eb669268712c65c866`.
It provided `PercolationMode`, public `Vec<bool>` fields in `OccupancyState`,
`cluster_stats(&CsrLattice, ...)`, `UnionFind`, and scheduler-ready
`PercolationMC`. These names, storage choices, constructor/parameter schema,
and observable names are provisional and are not compatibility commitments.
**F4 has since deleted this entire public API without a compatibility
wrapper** (per the roadmap section 7 migration rule for never-merged APIs) and
migrated every applicable validation asset to the final API; the description
below is the frozen behavioral record the migration was checked against.

Reference behavior:

- each sweep independently redraws uniform site/bond occupancy;
- `cluster_stats` allocates fresh union-find and boundary scratch on every call;
- physical undirected edges are scanned once, then tracked vertices are scanned;
- site/mixed omit closed sites; bond mode counts every isolated site;
- outputs are `MaxCluster`, raw `SecondMoment`, `NClusters`, `Spanning`, plus
  `Occupied` for pure modes or `OccupiedSites`/`OccupiedBonds` for mixed mode;
- crossing defaults are open square left/right columns and open chain endpoints;
  all other parameter-built topologies require explicit site sets.

The validation assets documented in [VALIDATION.md](VALIDATION.md) — exact
full enumeration and hand-derived 1D/2x2 forms, an independent flood-fill
reference, mixed-to-pure reductions, cross-family and random-graph parity,
scheduler tests, fixed-seed reproducibility, 16-seed z-score gates, and
ignored square/cubic critical-region checks — have all been migrated to the
F4 API (see the F4 section below).

The former PR #4 limitations — public storage layout, fresh per-analysis
allocations, no component identity, no giant fraction or finite-cluster
susceptibility — are resolved by F3/F4. Still true for the whole platform:
a uniform-only `FromParams` parameter path (heterogeneous laws enter via direct construction, N1 Validated),
owned `CsrLattice` coupling in the parameter path (borrowed graphs use the
core APIs directly), scalar results, and no wrapping, directed, process,
continuum, FK, attack, or Newman-Ziff implementation.

## Performance baseline procedure

`percolation_baseline` measures the F0 reference sample+analyze path, rebuilt on
the final API after F4: an i.i.d. redraw into a `Vec<bool>` reference
configuration plus the production `analyze` with a workspace allocated per call
(the historical PR #4 allocation profile). The 2026-09-03 F0 table below was
recorded with the original PR #4 implementation at commit `f6031bf`. It does not
enable, depend on, or link
`allocation-counter`; the allocation instrumentation lives in the separate
`percolation_allocations` target behind the `allocation-probe` feature. Cases
are named with topology, `V`, `E`, `p_site`, and `p_bond`; all modes use
`p_site = 0.5927`, `p_bond = 0.5`. Topologies are an open chain (`V=4096`),
open square (`64x64`), open cubic (`16^3`), seeded sparse Erdos-Renyi-like graph
(`V=4096`, `E=16384`), and seeded degree-proportional preferential-attachment
graph (`V=4096`, four distinct attachments per added vertex). The latter starts
with edge `(0,1)` and attachment pool `[0,1]`; setup asserts both `pool.len() =
2E` and, vertex by vertex, pool multiplicity equals graph degree.

```bash
# Criterion throughput only; allocation-counter is disabled and unlinked.
cargo bench -p cmc-rs --bench percolation_baseline --no-default-features -- --test
cargo bench -p cmc-rs --bench percolation_baseline --no-default-features

# Standalone allocation/storage probe; no Criterion timing occurs here.
cargo bench -p cmc-rs --bench percolation_allocations \
  --features allocation-probe -- --test
```

The standalone probe prints `VEC_BOOL_PROBE`, `PERCOLATION_ALLOC`, and
`TOPOLOGY_VIEW_ALLOC` records. The topology construction/scan probe uses a
separate measurement region from the F0 sample/analyze allocation baseline.
`site_logical_capacity` and `bond_logical_capacity` are the capacities reported
by `Vec<bool>` in logical elements; they are not byte counts by API contract.
`estimated_owned_storage_bytes` sums the owned CSR vector capacities plus
occupancy storage using the byte-per-element behavior measured for this exact
compiler. It excludes struct headers, allocator metadata, boundary query
vectors, RNG state, and analyzer scratch, so it is a documented estimate rather
than RSS or an allocator measurement of the complete case. Its per-vertex and
per-edge columns are two normalizations of the same estimate, not additive.

`allocation-counter` separately measures 128 post-warmup sample+analyze calls
in one thread and reports allocator-observed calls and bytes per sample plus the
peak live allocation count and bytes in the measurement region. There is no CI
performance threshold because allocator and timing measurements are
environment-sensitive.

### Recorded F0 run

Recorded on 2026-09-03 with:

- CPU: Intel Xeon Gold 6148 at 2.40 GHz; Linux 6.8.0-134-generic x86_64.
- Compiler: `rustc 1.98.0 (88d9e12ae 2026-08-18)`, LLVM 22.1.8;
  `cargo 1.98.0 (797e8a9bc 2026-08-05)`.
- Throughput command: `cargo bench -p cmc-rs --bench percolation_baseline
  --no-default-features`.
- Allocation command: `cargo bench -p cmc-rs --bench percolation_allocations
  --features allocation-probe -- --test`.
- Criterion short-baseline configuration: 10 samples, 100 ms warm-up, 250 ms
  measurement per case. Results below are point estimates; timing is not a CI
  gate.

The throughput target's normal dependency tree contained no `allocation-counter`,
and its final executable contained no `allocation_counter` marker. Thus the
following measurements do not include that crate's global allocator wrapper.

| Topology | V | E | Estimated owned storage bytes | Est. bytes/V | Est. bytes/E | Site samples/s | Bond samples/s | Mixed samples/s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| chain | 4096 | 4095 | 303,079 | 73.994 | 74.012 | 14,790 | 14,823 | 11,653 |
| square | 4096 | 8064 | 565,128 | 137.971 | 70.080 | 9,027 | 7,544 | 7,190 |
| cubic | 4096 | 11,520 | 941,320 | 229.814 | 81.712 | 7,469 | 6,392 | 5,308 |
| sparse ER | 4096 | 16,384 | 1,101,832 | 269.002 | 67.250 | 4,243 | 4,203 | 3,852 |
| preferential attachment | 4096 | 16,374 | 1,101,374 | 268.890 | 67.264 | 6,515 | 5,553 | 4,552 |

## F1 read-only topology capability

The F1 API remains **Experimental and provisional** after independent review.
`GraphView` and
`UndirectedGraphView` expose dense `VertexId`/`EdgeId`, physical edge endpoints,
and incidence iteration through static dispatch. The public iterator contract is
only `Iterator`; current CSR implementations may provide stronger iterator
properties internally without imposing them on future mmap, generated, or
filtered views. `GraphView::vertex_ids()` is a default convenience, not a
performance requirement for every implementation.

IDs are created publicly only by `GraphView::vertex_id` and
`UndirectedGraphView::edge_id`. Their transparent newtypes guarantee a dense
numeric representation but carry no graph provenance. An ID must be used only
with the unchanged validated view that produced it. Cross-graph IDs or IDs kept
across structural mutation violate the caller contract and may panic at
infallible accessors, but cannot cause memory unsafety. This matters for
`CsrLattice`: its existing public fields remain mutable, so construction or an
explicit successful `CsrLattice::validate()` is the precondition for trait use.

`BorrowedUndirectedCsr` borrows offsets, neighbors, physical edge IDs, and AoS
`[[usize; 2]]` physical endpoints without copying. Its constructor accepts the
empty graph (`offsets = [0]`) and validates offsets, ranges, endpoint agreement,
and exactly two correctly oriented incidences per physical edge. Parallel edges
and self-loops are preserved; adjacency ordering is unrestricted.
Offsets/neighbors alone are insufficient because stable physical edge IDs cannot
be inferred for parallel edges or self-loops.

Construction uses temporary validation state of one byte per physical edge and
drops it before returning the zero-payload view. `try_reserve_exact` capacity or
allocation failures return `TopologyError::ValidationCapacity`; allocator
configurations that abort on OOM cannot be promised recoverable.

Directed out/in signatures were a crate-private compile sketch at F1 time;
D1 has since promoted them to public traits with owned and borrowed CSR
storage (see the D1 section below).

`topology_view` compares direct owned access, generic owned access, and generic
borrowed access on an open 128x128 square (`V=16,384`, `E=32,512`, `I=65,024`).
Run with:

```bash
cargo bench -p cmc-rs --bench topology_view --no-default-features -- --test
cargo bench -p cmc-rs --bench topology_view --no-default-features
```

Recorded on 2026-09-03 on the same Xeon Gold 6148 and rustc 1.98.0 toolchain as
F0, using 20 samples, 250 ms warm-up, and 1 s measurement per case. These are the
second complete run after the three-path benchmark was added:

| Scan | Direct owned | Generic owned | Owned delta | Generic borrowed | Borrowed delta |
|---|---:|---:|---:|---:|---:|
| Physical edges | 24.211 us | 23.472 us | -3.05% | 13.768 us | -43.13% |
| CSR incidences | 58.033 us | 56.896 us | -1.96% | 57.978 us | -0.09% |

The first complete run also stayed below the 5% regression gate: generic owned
was -0.08% for edge scan and +0.32% for incidence scan; generic borrowed was
-42.93% and +0.89%, respectively. Borrowed edge scans benefit from compact AoS
endpoint input compared with the wider `Bond` records; this is a data-layout
effect, not evidence that trait dispatch accelerates work. There is no CI timing
threshold.

The separate allocation probe constructs and validates the borrowed view outside
the scan measurement region, then repeats a generic edge+incidence scan 128
times. For its open 64x64 square (`V=4,096`, `E=8,064`, `I=16,128`), the 64-byte
view header is an x86_64 rustc 1.98 measurement, not a cross-platform guarantee.
Construction performed one 8,064-byte scratch allocation, exactly one byte per
physical edge. The view owns zero payload bytes. The four borrowed input arrays
total 419,848 bytes (`102.502 bytes/V`, `52.064 bytes/E`); each steady-state scan
performed **0 allocations**, allocated **0 bytes**, and had **0 peak live
allocations/bytes**. Run the probe with:

```bash
cargo bench -p cmc-rs --bench percolation_allocations \
  --features allocation-probe -- --test
```

## F2 activity configuration and Bernoulli laws

The provisional F2 API separates generation from storage. `StaticConfiguration`
owns private vertex and physical-edge masks; the mask type and storage accounting
are crate-private and are not returned by any public method. `VertexActivity`
and `EdgeActivity` accept dense `VertexId`/`EdgeId` values under the same
unchanged-view provenance contract as F1: IDs from another graph or after
structural mutation violate the caller contract and may panic, but cannot cause
memory unsafety. The `law` module is private; reviewed law types are selectively
re-exported. The intended crate-root F2 surface comprises `StaticConfiguration`,
`VertexActivity`, `EdgeActivity`, `Probability`, `ProbabilityField`,
`ProbabilityError`, `SamplingError`, `SiteBernoulli`, `BondBernoulli`, and
`MixedBernoulli`.

The internal mask has logical all, none, and dense states while retaining any
dense byte allocation across endpoint transitions. A new pure site law
materializes only its vertex mask and makes edges logically all-active; a new
pure bond law materializes only its edge mask and makes vertices logically
all-active. Mixed sampling materializes both. Every successful sample writes
complete semantics, so switching law types cannot expose stale state.
`StaticConfiguration::resize` changes both domains and clears logical activity
while retaining usable capacity. Sampling checks graph/configuration and all
heterogeneous-field lengths before any mutation or RNG consumption.

`Probability` validates finite `[0, 1]` values once. `ProbabilityField` supports
uniform, borrowed, and owned validated values without per-sample revalidation.
Uniform exact endpoint fields use logical all/none states and consume no RNG for
site, bond, or any mixed endpoint combination. Heterogeneous fields are sampled
through the dense element loop; exact endpoint elements skip RNG individually,
but no O(1) or whole-field no-RNG contract is claimed for heterogeneous data.
N1 has since promoted heterogeneous sampling through the `StaticPercolationMC`
adapter to the Validated heterogeneous family (see the N1 section); the law
module itself remains the Experimental F2 substrate it always was.
Network/domain facades over heterogeneous fields are provided by N2
(see the N2 section); long-range/correlated combinations remain future
work.

The F2 statistical gate aggregates activity across every requested seed rather
than taking a maximum over a growing set of per-seed scores. For each of six
separate domains (uniform site, uniform bond, uniform mixed vertex/edge, and
heterogeneous mixed vertex/edge), probabilities `p_i` give exact total moments

```text
mean     = n_seeds * samples_per_seed * sum_i p_i
variance = n_seeds * samples_per_seed * sum_i p_i * (1 - p_i).
```

Uniform fields reduce to the binomial formula; heterogeneous fields use the
exact Poisson-binomial mean and variance. Each aggregate requires `|z| < 4.5`.
For six approximately normal two-sided gates, the resulting family-wise false
alarm probability is about `6 * 2 * Phi(-4.5) = 4.1e-5`. Zero variance requires
an exact success count. There is deliberately no fixed per-seed maximum,
sign-fraction, or dispersion gate whose rejection probability changes with
`SCUTTLE_ZSCORE_SEEDS`. The targeted test passed at 1, 16, 64, and 4096 seeds;
the observed six-gate z ranges were `[-1.450, 0.935]`, `[-0.759, 1.417]`,
`[-0.568, 1.049]`, and `[-1.782, 0.815]`, respectively. The 4096-seed run took
37.00 seconds on the recorded Xeon host.

`percolation_laws` times only `sample` on an open 64x64 square (`V=4,096`,
`E=8,064`); graph, laws, probability fields, configuration, and RNG are prepared
outside Criterion iteration. The short run uses 10 samples, 100 ms warm-up, and
250 ms measurement per case:

```bash
cargo bench -p cmc-rs --bench percolation_laws --no-default-features -- --test
cargo bench -p cmc-rs --bench percolation_laws --no-default-features
cargo bench -p cmc-rs --bench percolation_allocations \
  --features allocation-probe -- --test
```

Recorded 2026-09-03 on the same Xeon Gold 6148 and rustc 1.98.0 toolchain as F0.
The 2026-09-03 target contained **14 Criterion IDs** (three site, three bond,
seven mixed uniform, one heterogeneous mixed); N1 later added heterogeneous
site/bond for 16 total — see the N1 section. Throughput denominators are `V` for site,
`E` for bond, and `V + E = 12,160` for every mixed case.

| Law | `p_vertex` | `p_edge` | Point time | Point throughput |
|---|---:|---:|---:|---:|
| Site uniform | 0.01 | all | 6.974 us | 587.30 M vertices/s |
| Site uniform | 0.50 | all | 7.053 us | 580.72 M vertices/s |
| Site uniform | 0.99 | all | 7.200 us | 568.87 M vertices/s |
| Bond uniform | all | 0.01 | 13.708 us | 588.27 M edges/s |
| Bond uniform | all | 0.50 | 13.637 us | 591.32 M edges/s |
| Bond uniform | all | 0.99 | 13.869 us | 581.44 M edges/s |
| Mixed uniform | 0.01 | 0.01 | 20.643 us | 589.07 M entities/s |
| Mixed uniform | 0.50 | 0.50 | 21.097 us | 576.39 M entities/s |
| Mixed uniform | 0.99 | 0.99 | 20.685 us | 587.86 M entities/s |
| Mixed uniform | 0.01 | 0.50 | 20.917 us | 581.34 M entities/s |
| Mixed uniform | 0.99 | 0.50 | 20.592 us | 590.53 M entities/s |
| Mixed uniform | 0.50 | 0.01 | 20.543 us | 591.94 M entities/s |
| Mixed uniform | 0.50 | 0.99 | 21.084 us | 576.73 M entities/s |
| Mixed heterogeneous | alternating 0.25/0.75 | repeating 0.2/0.6/0.6 | 32.144 us | 378.30 M entities/s |

The allocation probe warms a mixed dense configuration, switches through exact
`p=0/1`, then measures 128 iterations for five laws across five topologies (25
law/topology records). Site uniform, bond uniform, mixed uniform, mixed
heterogeneous, and repeated mixed endpoint switching all reported
**0 allocations/iteration**, **0 allocated bytes/iteration**, and zero peak live
allocations for chain, square, cubic, sparse-ER, and preferential-attachment
cases. An endpoint-switch iteration performs two samples, `(0,1)` then `(1,0)`,
after both dense buffers have been materialized. The logical worst-case dense
payload is `V + E` bytes: 8,191 (chain), 12,160 (square), 15,616 (cubic),
20,480 (sparse ER), and 20,470 (preferential attachment). These figures are the
expected byte payload for the current dense backend, not a measurement of
private `Vec` capacity or a public runtime introspection API. They describe the
worst-case reusable two-mask state; a fresh pure law does not materialize the
irrelevant mask. They exclude struct headers, RNG state, probability-field
storage, topology, and allocator metadata.

F2 passed independent review and its engineering gates, but its public API
remains **Experimental** until the foundation milestone stabilizes. The
independently reviewed F3 analyzer/workspace is described below; the F4
composition that turns this substrate into scheduler observables is described
in the F4 section.

## F3 undirected component analyzer

The F3 production path is generic over borrowed `UndirectedGraphView`,
`VertexActivity`, and `EdgeActivity` capabilities. It has no site/bond/mixed mode:
a tracked vertex is exactly an active vertex, and a physical edge is active only
when the edge and both endpoints are active. Thus site laws use all-active edges,
bond laws use all-active vertices, and mixed laws compose naturally. Physical
edges are scanned once in edge-ID order; parallel edges remain distinct and a
self-loop contributes one active edge.

`ComponentSummary` records active vertex/edge counts, component count, largest
component identity and size, and raw `M2`. Scientific identity is the canonical
minimum `VertexId` in a component, never the union-find implementation root.
Equal-size largest components choose the lowest canonical identity. Empty
activity and empty borrowed graphs return zero counts, no largest identity, and
zero `M2`.

`BoundaryQuery::new` checks raw indices against the graph, sorts and deduplicates
both sets, and retains their validated vertex domain. Empty sets are valid and
never cross. Overlap means an active vertex in the intersection crosses by
itself. `analyze` answers several queries without cloning their vectors;
`analyze_with_labels` explicitly requests per-vertex canonical labels, with
inactive vertices labelled `None`. `AnalysisResult<'_>` borrows outcomes and
optional labels from `ComponentWorkspace`, so the borrow checker prevents a
subsequent analysis from invalidating a live result.

Activity and query dimensions are checked before workspace mutation. Workspace
capacity failures are typed. Because preparation is a sequence of fallible
`Vec::try_reserve_exact` calls, successful earlier reservations may remain after
a later failure, but no partial scientific output is returned. Allocators that
abort on OOM remain outside the recoverability guarantee. Equal-size graph ID
domains have no runtime provenance token; using IDs/query semantics from another
graph of the same size remains the documented F1 caller-contract violation.

Validation retains the PR #4 analyzer and its existing tests unchanged. New F3
coverage uses a separate incidence/stack flood fill rather than production UF:
full 2x2 site/bond/mixed enumeration; 7 topology families x 3 laws x 128 seeded
configurations; owned/borrowed CSR parity; and self-loop, parallel-edge,
disconnected, isolate, empty, overlap, empty/duplicate/multiple-query, labels,
stable-identity, typed vertex/edge/query errors, and grow/shrink/labels/query
workspace reuse cases. Eight targeted integration tests and the internal
stamp-wrap unit test pass.

The standalone allocation probe prepares topology, query, dense activity, and
workspace before measuring 128 sample+analyze calls. All 15 topology/law records
(chain, square, cubic, sparse ER, power law x site/bond/mixed) report **0.000
allocations/sample** and **0.000 allocated bytes/sample**, compared with the F0
reference's 7 allocations and 86,016 bytes per sample. The feature-gated,
doc-hidden audit snapshot accounts for every owned workspace vector. At
`V=4,096`, one query, and labels prepared, allocator observation and summed
capacity both equal 196,609 bytes across six allocations:

| Workspace vector | Capacity bytes |
|---|---:|
| parent | 32,768 |
| component size | 32,768 |
| canonical identity | 32,768 |
| query root stamps | 32,768 |
| query outcomes | 1 |
| optional labels | 65,536 |
| **Total** | **196,609** |

This is `48.000 bytes/V` and, for this one-query probe, `196,609 bytes/query`.
The latter is total workspace normalized by query count, not marginal query
storage: each additional query outcome currently costs one byte while root
stamps remain vertex-sized. A no-label workload does not reserve label storage.
Struct headers and allocator metadata are excluded.

The Linux-only RSS target runs one case per fresh process and reads
`/proc/self/status`; it adds no dependency or unsafe code. From before topology
construction through 8 warmups plus 128 mixed sample+analyze calls, square
(`V=4,096`, `E=8,064`) measured `VmHWM 2,096 -> 2,796 KiB` (+700 KiB), while
sparse ER (`V=4,096`, `E=16,384`) measured `2,100 -> 3,572 KiB` (+1,472 KiB).
These are whole-process lifetime high-water deltas including topology,
configuration, RNG, runtime, and allocator effects; they are not workspace
payload or portable RSS guarantees.

`percolation_components` contains 63 short-run Criterion IDs: production
sample+analyze, production analyze-only, and equivalent independent flood-fill
sample+analyze for five topologies x three laws, plus three historical PR #4
square measurements and, since F4, fifteen `static-mc` sweep+measure records. Before timing each topology/law pair, 128 identical samples
assert equality of all six summary fields and query outcomes. The production and
equivalent-reference timed paths both use `TIMED_SAMPLE_SEED`, advance exactly
`TIMED_INITIAL_SAMPLES = 1` sample before Criterion, and advance exactly one
sample per logical iteration. A separate 129-state check (the initial sample plus
128 following samples) compares every vertex and physical-edge activity flag,
proving that corresponding iteration indices produce identical configuration
sequences. Criterion may request different total iteration counts for separate
benchmark functions, but each path starts from and advances through the same
indexed sequence. Labels are disabled on both timed paths. Criterion throughput
uses physical `E`, so it reports edges/s directly; samples/s below is
`1 / point_time`. Displayed edges/s is recomputed as `E * displayed samples/s`
and rounded to `0.01 M`, with an automated arithmetic check. Parameters are
`V=4,096`, `p_site=0.5927`, and `p_bond=0.5`.

| Topology | Law | Production samples/s | Production edges/s | Equivalent reference samples/s | Reference edges/s | Speedup |
|---|---|---:|---:|---:|---:|---:|
| chain | site | 13,493 | 55.25 M | 11,480 | 47.01 M | 1.175x |
| chain | bond | 14,799 | 60.60 M | 9,475 | 38.80 M | 1.562x |
| chain | mixed | 11,258 | 46.10 M | 8,319 | 34.07 M | 1.353x |
| square | site | 9,738 | 78.53 M | 4,250 | 34.27 M | 2.291x |
| square | bond | 8,683 | 70.02 M | 3,324 | 26.80 M | 2.612x |
| square | mixed | 7,577 | 61.10 M | 2,737 | 22.07 M | 2.769x |
| cubic | site | 7,608 | 87.64 M | 3,844 | 44.28 M | 1.979x |
| cubic | bond | 6,540 | 75.34 M | 2,963 | 34.13 M | 2.207x |
| cubic | mixed | 5,812 | 66.95 M | 2,260 | 26.04 M | 2.572x |
| sparse ER | site | 3,892 | 63.77 M | 2,935 | 48.09 M | 1.326x |
| sparse ER | bond | 4,220 | 69.14 M | 2,090 | 34.24 M | 2.019x |
| sparse ER | mixed | 4,026 | 65.96 M | 2,102 | 34.44 M | 1.915x |
| power law | site | 6,597 | 108.02 M | 3,122 | 51.12 M | 2.113x |
| power law | bond | 5,521 | 90.40 M | 2,098 | 34.35 M | 2.632x |
| power law | mixed | 4,806 | 78.69 M | 2,291 | 37.51 M | 2.098x |

Production analyze-only throughput for the same fixed configurations:

| Topology | Site samples/s / edges/s | Bond samples/s / edges/s | Mixed samples/s / edges/s |
|---|---:|---:|---:|
| chain | 17,440 / 71.42 M | 18,884 / 77.33 M | 15,536 / 63.62 M |
| square | 10,698 / 86.27 M | 10,187 / 82.15 M | 9,660 / 77.90 M |
| cubic | 8,237 / 94.89 M | 7,574 / 87.25 M | 7,255 / 83.58 M |
| sparse ER | 3,882 / 63.60 M | 4,780 / 78.32 M | 4,792 / 78.51 M |
| power law | 6,868 / 112.46 M | 6,625 / 108.48 M | 5,930 / 97.10 M |

The retained PR #4 square timing reproduces the historical per-call allocation
profile on the final analyzer; it is only a **historical non-equivalent migration
baseline**. It is not used to
claim speedup. The equivalent flood-fill table above is the F3 performance
comparison. Short-run point estimates are not CI thresholds. Run the evidence
with:

```bash
cargo test -p cmc-rs --test suite percolation_components --no-default-features
cargo bench -p cmc-rs --bench percolation_components --no-default-features -- --test
cargo bench -p cmc-rs --bench percolation_components --no-default-features
cargo bench -p cmc-rs --bench percolation_allocations \
  --features allocation-probe -- --test
SCUTTLE_RSS_CASE=square SCUTTLE_RSS_LAW=mixed \
  cargo bench -p cmc-rs --bench percolation_components_rss --no-default-features
SCUTTLE_RSS_CASE=sparse-er SCUTTLE_RSS_LAW=mixed \
  cargo bench -p cmc-rs --bench percolation_components_rss --no-default-features
```

F3 passed independent review and its engineering gates but remains
**Experimental**. F4 (below) adds the observable normalization,
finite-cluster susceptibility, and Carlo.rs adapter that turn this analyzer
into production-path observables; with the F4 review landed, the static family
is the platform's first Validated row.

#### `Vec<bool>` storage evidence

On this compiler, the reviewer concern is confirmed: these `Vec<bool>` values
are allocator-observed as byte-backed, despite `capacity()` being expressed in
logical `bool` elements. Each probe performed exactly one allocation:

| Logical len/capacity | Allocator-observed bytes |
|---:|---:|
| 1 | 1 |
| 63 | 63 |
| 64 | 64 |
| 65 | 65 |
| 4096 | 4096 |
| 8064 | 8064 |
| 11,520 | 11,520 |
| 16,384 | 16,384 |

This is empirical behavior of `rustc 1.98.0`, not a promised representation of
`Vec<bool>` across toolchains. The probe remains executable so future baselines
can detect a representation change instead of assuming bit packing or byte
storage from the type name.

For every topology and mode, the post-warmup sample+analyze measurement returned
**7.000 allocations/sample**, **86,016 allocated bytes/sample**, and a peak of
**7 live allocations / 86,016 live bytes** in the serial measurement region.
These are analyzer reference costs, not targets. Both bench targets' `--test`
commands completed successfully; Criterion exercised all 14 F2 law throughput
IDs listed above. Exact raw Criterion intervals remain in `target/criterion/` on
the measurement host.

## F4 scientific observables and production adapter

F4 passed independent review and all engineering gates. It composes
the F1-F3 foundation into the final scheduler surface and deletes the
provisional PR #4 public API in one migration (no compatibility wrapper, per
the roadmap rule for never-merged APIs). The provisional API PR #4 introduced
is superseded in place on this branch; PR #4 carries the final architecture as
its delivery vehicle. The static family is **Validated** and becomes
production-ready when the branch merges into `dev`.

### Observable definitions and Results semantics

`ObservablePlan` explicitly requests scalar observables and boundary queries.
Duplicate requests collapse (one sample per name per realization), and labels
are intentionally not part of the plan: every scalar follows from
`ComponentSummary`, whose canonical lowest-ID largest component realizes the
one-component tie policy. Each boundary query records `BoundaryCrossing0`,
`BoundaryCrossing1`, ... in plan order as a 0/1 indicator.

| Results name | Definition | Notes |
|---|---|---|
| `ActiveVertexCount` | occupied sites in the configuration | mixed mode counts sites only |
| `ActiveEdgeCount` | physical edges that are open under the law and have both endpoints occupied | site mode: edges between two occupied sites; bond mode: open bonds; mixed: both conditions |
| `ComponentCount` | tracked components among occupied sites | bond mode counts every isolate |
| `LargestSize` | `S_max`, `0` when nothing is occupied | |
| `GiantFraction` | `S_max / V`, denominator fixed to total topology vertices `V` | undefined and absent at `V = 0` |
| `RawSecondMoment` | `M2 = sum_i s_i^2`, no normalization, largest included | u128 accumulation |
| `FiniteClusterSusceptibility` | per sample, exclude one largest component (ties: canonical lowest identity), then `chi = (M2 - S_max^2) / (V_active - S_max)` | exclusion before averaging |
| `FiniteClusterSusceptibilityDefined` | 0/1 indicator recorded for every sample | 0 exactly when the chi denominator is zero |

Undefined samples never write the chi scalar: the `FiniteClusterSusceptibility`
row in `results.json` has an `n_bins` counting only defined samples, while the
indicator row counts every sample. This makes the defined fraction explicit
instead of silently emitting zero or NaN. The per-sample exclusion cannot be
reconstructed from separately aggregated `LargestSize`/`RawSecondMoment` means
(unit-tested counterexample against the ratio of means).

Moment observables convert u128 squares up to `V^2` into f64. Plans whose
topology could exceed the exact f64 integer range (`V^2 > 2^53`, i.e.
`V > 94,906,265`) are rejected at construction with
`ObservablePlanError::InexactMomentConversion`, for `RawSecondMoment` and for
the susceptibility numerator, so no silent precision loss is possible.

### Adapter and parameter schema

`StaticPercolationMC` composes an owned validated `CsrLattice`, the uniform
`StaticLaw` enum (site/bond/mixed, used only at this boundary),
`StaticConfiguration`, a prepared `ComponentWorkspace`, and the
`ObservablePlan`. `sweep()` only resamples occupancy; `measure()` analyzes
with the preallocated workspace and records exactly the planned scalars.
Configure `thermalization_sweeps = 0` (i.i.d. samples; the adapter cannot see
scheduler configuration to enforce it). The `FromParams` schema is:

- `mode`: `site` (default) / `bond` / `site-bond`; pure modes take `p`, mixed
  takes `p_site` + `p_bond`, and the wrong key for the mode is a typed error;
- lattice parameters reuse the standard builders (`chain`/`square`/`cubic`/
  `hypercubic`/...), `pbc` defaults to `false`;
- `observables`: comma-separated selection defaulting to the full plan;
- `spanning_from`/`spanning_to`: explicit comma-separated site lists, always
  both together; defaults are square left/right columns and chain endpoints,
  and any other parameter-built topology is rejected loudly (construct the
  adapter directly with a `BoundaryQuery` for arbitrary sets).

`validate_params` accepts exactly what `from_params` accepts. Borrowed graphs
and custom laws bypass the adapter through `analyze`/`analyze_with_labels` or
compose a runtime with `Run::from_parts()`.

### Validation

Every applicable PR #4 gate was migrated to the final API in
`tests/physics/percolation.rs` (14 tests, 2 `#[ignore]` long) and
`tests/physics/percolation_zscore.rs` (6 z-score tests): 2x2 exhaustive
enumeration vs hand-derived moments, site spanning polynomial, mixed closed
form `2 p_s^2 p_b - p_s^4 p_b^2`, mixed-to-pure reduction identities, chain
closed forms, union-find vs independent flood-fill across seven topology
families and the seeded random graph, crossing monotonicity (pure and mixed),
scheduler 2x2 exact moments and chain closed forms, fixed-seed bitwise
reproduction, `p = 0`/`p = 1` occupation extremes, the ignored 32x32 bond `p_c = 1/2` crossing and 16^3 cubic
critical-bracket long tests, and `BoundaryCrossing0`/`LargestSize` z-score
gates for site/bond/mixed (16 default seeds, `SCUTTLE_ZSCORE_SEEDS`-scalable,
verified at 64 seeds). New F4 unit coverage (10 tests in `observable.rs` +
`carlo.rs`): tie exclusion, interleaved duplicate collapse, no-active and
`V = 0`, single-giant-component undefined chi with indicator counts under the
scheduler, per-sample exclusion vs ratio-of-means counterexample, exact-f64
conversion limit, JSON schema, and the `FromParams` rejection matrix.

### Steady-state allocation gate

The `percolation_allocations` probe (feature `allocation-probe`) adds
`report_adapter_case`: it warms 8 full `sweep()+measure()` cycles, then
measures 128 more through `StaticPercolationMC` and a Carlo.rs `Context` with
a binsize strictly larger than warmup plus measurement, so no accumulator bin
completes inside the measured window. All **15 topology/law records** (chain,
square, cubic, sparse ER, power law x site/bond/mixed; 9 observables each)
report **0.000 allocations/sample**, **0.000 allocated bytes/sample**, and
zero peak live allocations/bytes:

```text
PERCOLATION_ADAPTER_ALLOC ... allocations_per_sample=0.000
allocated_bytes_per_sample=0.000 peak_live_allocations=0 peak_live_bytes=0
```

This gate required one behavior-preserving Carlo.rs change: with no
measurement namespace, `Context::measure*` previously allocated a fresh
`String` for every observable name on every call; `qualified_measurement_name`
now returns `Cow::Borrowed` in that case (`Cow::Owned(format!(...))` when a
parallel-tempering namespace is active, exactly as before). Recording nine
scalars per sample is otherwise allocation-free.

### Criterion data

`percolation_components` adds 15 `static-mc-sweep+measure` IDs (full adapter:
law sampling through the adapter, workspace analysis, and nine scalar
recordings) next to the F3 `sample+analyze` IDs. Recorded 2026-10-02 on the
same Xeon Gold 6148 / rustc 1.98.0 host, short-run configuration (10 samples,
100 ms warm-up, 250 ms measurement), **while the host carried an unrelated
sustained load of ~45/72 cores**, so absolute point estimates are inflated
relative to the 2026-09-03 F3 table by roughly 8-15% across all IDs —
including `analyze-only`, a path F4 does not touch and did not modify (the
F1-F3 core sources are byte-identical since the F3 record; the diff touches
only exports, the two new F4 files, benches, and tests). The load-independent
same-run comparison is therefore the meaningful F4 evidence: the adapter
composition costs the same as the bare core within contention noise (median
about -2%, scatter roughly +/-7%, one -23% power-law/site outlier where the
adapter measured faster).

| Topology | Law | Core sample+analyze (us) | Adapter sweep+measure (us) | Adapter samples/s | Adapter edges/s | Adapter vs core |
|---|---|---:|---:|---:|---:|---:|
| chain | site | 86 | 87 | 11,464 | 46.94 M | +1.9% |
| chain | bond | 85 | 79 | 12,695 | 51.99 M | -7.3% |
| chain | mixed | 100 | 101 | 9,940 | 40.71 M | +1.0% |
| square | site | 114 | 118 | 8,470 | 68.30 M | +3.5% |
| square | bond | 133 | 132 | 7,588 | 61.19 M | -1.2% |
| square | mixed | 155 | 145 | 6,891 | 55.57 M | -6.4% |
| cubic | site | 149 | 146 | 6,862 | 79.06 M | -2.4% |
| cubic | bond | 171 | 165 | 6,056 | 69.76 M | -3.2% |
| cubic | mixed | 191 | 188 | 5,308 | 61.15 M | -1.5% |
| sparse ER | site | 278 | 272 | 3,681 | 60.31 M | -2.4% |
| sparse ER | bond | 256 | 249 | 4,016 | 65.79 M | -2.8% |
| sparse ER | mixed | 274 | 265 | 3,774 | 61.83 M | -3.4% |
| power law | site | 218 | 168 | 5,966 | 97.68 M | -23.2% |
| power law | bond | 179 | 190 | 5,250 | 85.97 M | +6.5% |
| power law | mixed | 239 | 220 | 4,540 | 74.34 M | -8.0% |

For the F3 no-regression gate, the F3 core paths (`sample+analyze`,
`analyze-only`, `equivalent-reference`, laws) are unchanged code and their
same-run parity/sequence assertions still pass; the cross-day point-estimate
deltas under the current load are environmental, as evidenced by the equal
inflation of the untouched `analyze-only` IDs. Re-run on a quiet host to
reproduce absolute numbers:

```bash
cargo bench -p cmc-rs --bench percolation_components --no-default-features -- --test
cargo bench -p cmc-rs --bench percolation_components --no-default-features
cargo bench -p cmc-rs --bench percolation_allocations \
  --features allocation-probe -- --test
```

### Remaining F4 limitations

At F4 the adapter was uniform-only; N1 has since added the heterogeneous
`StaticLaw` variants on the same adapter (see the N1 section below), so the
uniform-only restriction no longer applies. The remaining limitations are
unchanged: the parameter path owns its lattice (borrowed graphs use
the core APIs); one boundary-query normalization per adapter plan; scalar
results only (vector observables need C0); wrapping/winding exists since W2
as Experimental `Wrapping`/`WindingRank` plan observables that require an
embedded periodic adapter (generator vectors remain analyzer-level); scalar
Newman-Ziff site/bond scans exist since Z1 (no scheduler adapter or vector
output until C0); no directed/process/continuum/FK/attack families; chi
tie conventions other than single canonical-largest exclusion are separate
named policies that do not exist.

## N1 heterogeneous Bernoulli laws

N1 promotes the F2 heterogeneous probability substrate to the production
adapter boundary. The analyzer, workspace, and observable plan are untouched
(the roadmap requires this); only the owned adapter surface grows.

### API

`StaticLaw` gains three owned heterogeneous variants:
`SiteHeterogeneous { probabilities: Vec<Probability> }`,
`BondHeterogeneous { probabilities: Vec<Probability> }`, and
`MixedHeterogeneous { site_probabilities, bond_probabilities }`. The
constructors `StaticLaw::site_heterogeneous` / `bond_heterogeneous` /
`mixed_heterogeneous` take `impl Into<Vec<Probability>>` (owned vectors or
borrowed slices) and cannot fail: every value is already a validated
`Probability`. Field lengths are checked against the adapter topology in
`StaticPercolationMC::new` — a mismatch returns the new typed
`StaticPercolationError::LawDomain(SamplingError::ProbabilityFieldLength {..})`,
reusing the F2 vocabulary — so a constructed adapter's `sweep()` can never hit
a domain error. `sweep()` borrows the fields (`ProbabilityField::Borrowed`),
so per-sweep sampling remains allocation-free. The `FromParams` schema stays
uniform-only (roadmap §9 N1 adds no probability-list input format);
heterogeneous adapters are scheduled through `Run::from_parts()` with a seeded
`Context`, exactly like other non-`FromParams` models. `StaticLaw::law()` now
returns a reference (the enum is no longer `Copy` because it owns vectors).

### Validation gates (roadmap §9 N1)

All gates live in `tests/physics/percolation_heterogeneous.rs` (7 tests) and
`tests/physics/percolation_heterogeneous_zscore.rs` (7 z-score tests), plus
3 adapter unit tests in `src/percolation/carlo.rs`:

- **Exact non-identically-distributed enumeration.** The 2x2 open square is
  enumerated with per-element weights `prod_i p_i^{x_i} (1-p_i)^{1-x_i}`
  (four distinct values per domain); scheduler statistics for
  `BoundaryCrossing0`, `LargestSize`, and `RawSecondMoment` match within
  |z| < 4 for site, bond, and mixed, over 4 seeds each at 200k sweeps. The
  weight products are checked to sum to one. Hand-derived heterogeneous chain
  closed forms `P = prod p_i` (site), `prod q_j` (bond), and the product of
  both (mixed) match through the full scheduler stack at L = 6 within |z| < 4.
  The 3x3 mixed enumeration reference of the z-score tests was independently
  reproduced during independent N1 review by a from-scratch Python
  flood-fill enumeration (different code and summation order) agreeing to
  ~1e-13-1e-14 relative difference.
- **Uniform field reduction, bitwise.** A heterogeneous field of identical
  `p` values consumes RNG in the same order as the uniform law (and consumes
  none at exact endpoints), so identical seeds reproduce identical `Results`
  means bit-for-bit across six observables — for site, bond, and mixed laws,
  at `p in {0.0, 0.3, 0.5, 1.0}`.
- **Mixed reduction identities, bitwise.** Heterogeneous mixed with an all-one
  edge field equals the heterogeneous site law, and with an all-one site field
  equals the heterogeneous bond law, bitwise at identical seeds.
- **Endpoint exactness.** All-0/1 fields produce deterministic configurations;
  hand-checked chain cases pin `ActiveVertexCount`, `ActiveEdgeCount`,
  `ComponentCount`, `LargestSize`, `RawSecondMoment`, and crossing exactly.
- **Length and invalid-probability rejection.** Wrong-length site/bond/mixed
  fields are typed `LawDomain` rejections at the adapter constructor (unit
  tests assert the exact payload); no unvalidated f64 can reach the
  heterogeneous constructors (`Probability::new` is the only entry).
- **Scheduler z-scores.** 16 default seeds (SCUTTLE_ZSCORE_SEEDS-scalable,
  verified at 64): `BoundaryCrossing0` and `LargestSize` against full
  heterogeneous enumeration on the 3x3 square for site (2^9), bond (2^12),
  and mixed (2^21, shared across its two tests), with per-seed |z| < 4,
  |z̄| < 1.5, and no one-sided bias; plus a Poisson-binomial activity gate
  (site `ActiveVertexCount` mean `sum p_i`, bond `ActiveEdgeCount` mean
  `sum q_j`, mixed `ActiveEdgeCount` mean `sum_j q_j p_u(j) p_v(j)`, |z| < 4
  per seed). The heterogeneous scheduler runs use the labeled seed base
  `0x4845_5445_0000` ("HETE").

### Uniform versus heterogeneous sampling benchmark

`percolation_laws` adds two Criterion IDs (heterogeneous site and bond next to
the existing heterogeneous mixed), for **16 IDs** total. Recorded 2026-10-08
on the same Xeon Gold 6148 / rustc 1.98.0 host family as the 2026-09-03 F2
record, short-run configuration (10 samples, 100 ms warm-up, 250 ms
measurement), open 64x64 square (`V=4,096`, `E=8,064`); heterogeneous fields
alternate mid-range values (0.25/0.75 per vertex, 0.2/0.6 per edge — no
endpoint skipping). Point estimates:

| Family | Uniform `p = 0.5` | Heterogeneous | Heterogeneous time delta |
|---|---:|---:|---:|
| Site | 7.009 us / 584.4 M vertices/s | 10.324 us / 396.7 M vertices/s | +47.3% |
| Bond | 13.543 us / 595.4 M edges/s | 20.453 us / 394.3 M edges/s | +51.0% |
| Mixed | 20.571 us / 591.1 M entities/s | 30.952 us / 392.9 M entities/s | +50.5% |

The ~50% cost is the per-element probability load and compare in the dense
heterogeneous loop versus the hoisted uniform scalar; the F2 record measured
the same effect (378.30 M entities/s for heterogeneous mixed on 2026-09-03).
These are point estimates, not CI thresholds. Run with:

```bash
cargo bench -p cmc-rs --bench percolation_laws --no-default-features -- --test
cargo bench -p cmc-rs --bench percolation_laws --no-default-features
```

### Heterogeneous adapter allocation gate

`percolation_allocations` adds `report_heterogeneous_adapter_case`: after 8
warmup `sweep()+measure()` cycles through `StaticPercolationMC` with
heterogeneous fields, 128 measured cycles report **0.000 allocations/sample,
0.000 allocated bytes/sample, and zero peak live allocations/bytes** for all
**15 topology/law records** (chain, square, cubic, sparse ER, power law x
heterogeneous site/bond/mixed; 9 observables each), exactly matching the
uniform adapter gate:

```text
PERCOLATION_ADAPTER_ALLOC ... law=heterogeneous-site|bond|mixed
allocations_per_sample=0.000 allocated_bytes_per_sample=0.000
peak_live_allocations=0 peak_live_bytes=0
```

### Remaining N1 limitations

The analyzer is untouched (heterogeneity enters only through activity
sampling, as the roadmap requires). The `FromParams` schema is uniform-only;
heterogeneous adapters are constructed directly and scheduled with
`Run::from_parts()`. Per-element probabilities must be validated upstream
(`Probability::new`); the adapter validates only field lengths. Network
facades over heterogeneous site fields (random failure, original-degree
observables) are provided by N2, and targeted attacks by N3; combinations
with long-range/correlated/augmented networks remain future work.

## W1 periodic embedding capability

Status: **Experimental** — independently reviewed; the embedding
capability is provisional-until-merge per the never-merged-API rule
(roadmap §16 W1). Engineering gates and tests below are green. W1 is the
data layer; W2 (the winding analyzer that consumes it) is Validated.

### API and semantics

- `topology::PeriodicEmbedding` is a read-only capability next to
  `GraphView`/`UndirectedGraphView`: `periodic_dimensions()` returns the
  embedding dimension `D >= 1`, and
  `edge_displacement(edge, from) -> Option<DirectedDisplacement>` returns the
  integer fundamental-cell displacement of one **directed incidence** — the
  pair (physical edge, endpoint) — or `None` when `from` is not an endpoint.
  The key is never the endpoint pair or a bond label, so parallel edges
  between the same pair keep independent displacements.
- `DirectedDisplacement<'_>` is a lazy signed view over the stored
  source-to-target vector; querying the same physical edge from the other
  endpoint negates it exactly, so antisymmetry
  (`d(v->u) = -d(u->v)`) is structural, not a validated copy. Accessors:
  `iter()`, `axis(a)`, `dimensions()`, `is_zero()`.
- `LatticeEmbedding` owns one vector per physical edge (flat, edge-major,
  `D` `i32` components, `16 + 4D` bytes per edge) plus the endpoint pair it
  was validated against. `try_from_edge_displacements` rejects, before
  storing anything: zero periodic dimensions, wrong table length, self-loop
  edges (the two directed incidences of a self-loop would need opposite
  displacements but share one `(edge, endpoint)` key), `i32::MIN` components
  (unnegatable), and unresolvable capacity. It does **not** re-verify cocycle
  consistency: that is a producer property, pinned for the built-in builders
  by the plaquette and winding tests below and exploited by the W2 analyzer
  rather than re-validated per construction.
- Open boundaries have no embedding at all: every `*_with_embedding` builder
  rejects `pbc = false` with `EmbeddingError::OpenBoundaries`. A non-periodic
  graph is represented by not having the capability, so a wrapping query on a
  non-embedded graph is an explicit consumer-side error, never a vacuous
  "no winding" answer.
- No-inference principle (the roadmap W1 constraint): displacements are
  recorded by the builder at the only stage where wrap information is
  unambiguous — while bonds are emitted. `BondType` and vertex IDs are never
  consulted for geometry: on a 2-wide torus the two parallel bonds between a
  pair of sites carry different displacements (`0` and `-1`), so endpoint
  geometry alone cannot identify them.
- `BorrowedUndirectedCsr` intentionally does not implement the capability: a
  CSR topology carries no displacement payload and post-hoc reconstruction
  would be the forbidden inference. An embedding is separate data supplied by
  whoever knows the geometry, beside the borrowed topology.

### Builder coverage

- `build_chain_with_embedding`, `build_square_with_embedding`,
  `build_hypercubic_with_embedding` (all `pbc = true`), via a shared
  `hypercubic_bonds` bond loop that tracks wraps; and
  `build_triangular_with_embedding` via `triangular_bonds`, whose diagonal
  bonds lift `(x, y)` to `(x + 1, y + 1)`. Each variant returns exactly the
  same lattice as its plain builder (pinned by test) plus the embedding.
- Honeycomb and kagome embeddings are delivered by the roadmap v2 T3 stage
  (see the T3 section below): both builders were restructured from
  symmetric-adjacency assembly to direct bond emission with recorded
  displacements, removing the W1 deferral reason — `from_adjacency` dropped
  the per-bond construction direction while deduplicating. The small-cell
  obstruction is now understood precisely (T3 review, algebraically
  verified): honeycomb `Lx = 2` genuinely DOUBLES horizontal bonds (two
  distinct physical bonds per pair with displacements `(0,0)` and
  `(-1,0)`), while kagome never has parallel bonds — its `Lx = 2`
  obstruction was a covering-image ATTRIBUTION ambiguity (the wrap pair
  `{A(0,y), B(1,y)}` and the non-wrap pair `{A(1,y), B(0,y)}` are distinct
  endpoint pairs whose displacements cannot be told apart from endpoints
  alone), not bond multiplicity. Post-hoc displacement inference from
  either would be guessing; emission-time recording removes the need.

### Validation gates

All in `CMC.rs` unit tests (`lattice::graph::tests`, `topology::embedding::tests`):

- **Hand-derived 2x2 torus:** every one of the 8 torus bonds is asserted in
  both directions against a hand-written displacement table; each of the four
  vertex pairs carries one in-cell and one wrapping bond.
- **Independent geometry re-derivation (property test):** for 8 dimension
  sets (1D L=2,3,5; 2x2, 3x3, 4x3; 3D 2x3x4 and 3x1x2), every bond's
  advancing axis and wrap flag are recomputed from the builder's documented
  mixed-radix coordinates and compared componentwise, including exact
  negation from the other endpoint and `is_zero`.
- **Parallel edges:** the L=2 ring's two parallel bonds carry `0` and `-1`;
  the 2x2 square pair `{0,1}` carries `(0,0)` and `(-1,0)`; the 2x2
  triangular double diagonals between sites 1 and 2 carry `(1,0)` and
  `(0,-1)`.
- **OBC semantics:** `pbc = false` is rejected with
  `EmbeddingError::OpenBoundaries` by the chain, square, and hypercubic
  embedding builders.
- **Cocycle consistency (Z^D embedding):** on the 3x3 square torus, all 9
  elementary plaquettes accumulate `(0,0)`; on the 3x3 triangular torus, all
  18 elementary triangles (both orientations) accumulate `(0,0)`.
- **Winding of fundamental cycles:** axis-ring walks on L=2/3/4 squares, a
  5-ring, and the 2x3x2 cubic torus accumulate exactly `+1` on the walked
  axis and `0` elsewhere, taking the bond emitted from each site (the L=2
  parallel-bond direction).
- **Table construction rejections:** zero dimensions, wrong length,
  self-loops, `i32::MIN`, non-endpoint queries, and generic-code acceptance
  of the `PeriodicEmbedding` trait; plus bit-identical lattice reproduction
  by every embedding builder.

### W2 consumption

The W2 consumer is exactly the pass this section sketched: one pass over
physical edges calling `embedding.edge_displacement(edge, from)` next to
`UndirectedGraphView::edge_endpoints`, running union-find with integer
potentials. The `(edge, endpoint)` lookup key means that pass never needs
`BondType` or endpoint-pair geometry. See the W2 section below for the
delivered analyzer, its gates, and its benchmark data.

### Limitations

No benchmarks were added: lookups are O(1) table reads with no allocation
(the storage note above is the capacity claim), and there is no analyzer to
bench. The cocycle property is a builder contract verified by tests, not
re-checked by `try_from_edge_displacements`. An embedding detached from its
lattice has no provenance guard — the same caller contract as the F1 view
IDs (both must stay unchanged together).

## W2 winding analyzer

Status: **Validated** — independent review passed 2026-10-22 (roadmap §16
W2). W2 is the analysis layer over
the W1 displacement data: one physical-edge scan computes connectivity and
winding together.

### API and semantics

- `analyze_winding(graph, embedding, vertex_activity, edge_activity, workspace)`
  (`percolation/analyzer/winding.rs`) is mode-free exactly like F3: it
  consumes an `UndirectedGraphView` plus `VertexActivity`/`EdgeActivity`
  (active edge = edge active with both endpoints active) and one
  `PeriodicEmbedding` table, and looks displacements up by the directed
  incidence `(edge, endpoint)`. It is an independent implementation on its
  own `WindingWorkspace` (`percolation/workspace/winding.rs`) rather than a
  layer over the F3 analyzer: potentials need their own path-compressed
  union-find, so a fused single scan is the roadmap's shared-traversal
  orientation; F3's semantics and workspace are untouched.
- The scan is a union-find with integer `Z^D` potentials: every vertex
  stores its displacement toward its union-find parent, unions across
  active edges update potentials with **checked `i32` arithmetic**
  (`i32` ownership is W2's: W1 stores raw `i32`, W2 guarantees no silent
  wrap — an accumulation that would leave `i32` returns the typed
  `WindingError::DisplacementOverflow { edge }` attributed to the physical
  edge being scanned; the check is conservative, never wrapping), and every
  same-root closure contributes its cycle-sum mismatch as a winding
  generator observation.
- Generators are folded incrementally into a per-root first-nonzero-pivot
  **echelon basis of at most `D` rows** (extended-gcd row combinations keep
  the generated subgroup invariant; leading components positive). Results
  therefore report per component the **winding rank and a generating set of
  the winding subgroup** — never only a Boolean: `WindingResult::generators`
  yields exactly `rank` vectors of `D` signed cell components (an echelon
  generating set, not a canonical Hermite normal form).
- `WindingSummary` carries the F3-parity connectivity scalars plus
  `wrapping_component_count` and the canonical largest component's winding
  rank; `ComponentWinding` lists every component in ascending canonical
  identity order (identity = lowest active vertex ID, the F3 policy; the
  canonical-largest tie policy is identical). The empty configuration
  yields the empty summary.
- Typed rejections before any scanning: vertex/edge activity length
  mismatches, zero embedding dimensions, and a missing displacement for an
  active-edge endpoint (`WindingError::MissingDisplacement` — a table that
  does not describe this graph). OBC/non-embedded graphs have no embedding
  capability at all, so wrapping requests there are explicit consumer-side
  errors (below), never a vacuous "no winding".
- Adapter integration: `StaticPercolationMC::with_embedding(lattice,
  embedding, law, plan)` pairs the table with the plan. Construction
  validates the edge count (`EmbeddingMismatch`), endpoint coverage
  (`EmbeddingEndpoints`), and a static magnitude gate `max |component| *
  (2V + 2) <= i32::MAX` (`EmbeddingMagnitude`), which bounds every `i32`
  potential sum the infallible per-sweep `measure()` can produce. A plan
  requesting `Wrapping`/`WindingRank` without an embedding is
  `WindingWithoutEmbedding`; the parameter schema accepts
  `observables = "wrapping,winding-rank"` with `pbc = true` on chain,
  square, cubic, hypercubic, and triangular topologies (built with their
  embeddings) and rejects open/honeycomb/kagome winding requests loudly.
  The winding scan runs only when the plan requests a winding observable.
- Scheduler-level scalars (plan-explicit): `Wrapping` is `1` when any
  tracked component has winding rank > 0 (`0` otherwise); `WindingRank` is
  the canonical largest component's rank in `0..=D` (an ensemble-averaged
  integer, not a probability). Generator vectors stay analyzer-level until
  Carlo.rs grows vector estimates (C0). The default observable set and
  `ObservablePlan::all` do not include either name.

### Validation gates

All in `tests/physics/percolation_winding.rs` (16 tests) plus 13 unit tests
(`percolation/analyzer/winding.rs`, `percolation/observable.rs`,
`percolation/carlo.rs`).

- **Independent reference (roadmap gate):** a finite-cover lift-BFS — the
  torus replicated into a `5^D` cover (radius 2), BFS from each component's
  seed pinned in the center cell, winding read off the cells where the seed
  reappears, rank by plain f64 Gaussian elimination (real rank = subgroup
  rank). This shares no algorithmic path with the production weighted
  union-find. Radius 2 provably suffices on the tested tori (a simple cycle
  winding k times around an axis of an L-torus needs >= kL edges). Every
  hand-constructed case below is parity-checked against it, and the full
  enumerations check it configuration-by-configuration. The Criterion
  target re-runs the parity on sampled tori, skipping cover-clipped
  samples.
- **Hand-constructed cycles:** an interior 3x3 plaquette is contractible
  (rank 0); a row ring wraps with generator exactly `[1, 0]`; row ring +
  column ring sharing a site is one component of rank 2 with generators
  `[1, 0]`, `[0, 1]`; a hand-built 4-cycle whose lift accumulates `(1, 1)`
  has rank 1 with generator `[1, 1]` — distinguished from a single-axis
  wrap by the generator, not the rank; the L=2 ring's two parallel bonds
  (displacements `0` and `-1`) wind only through the cycle (the wrap-only
  bond set is a connected open path, rank 0 — the `(edge, endpoint)` key
  path for parallel edges); the full 2x3x4 cubic torus reports rank 3 with
  the three unit generators.
- **Invariances (roadmap gate):** over 8 random occupation samples on the
  3x3 torus, gauge shifts `d'(s->t) = d(s->t) + phi(t) - phi(s)`
  (representative/origin choice; cycle sums telescope) and full edge
  direction reversal (negated table) leave summaries equal and generated
  lattices identical (compared through a test-side Hermite normal form);
  cyclic relabeling of a winding 4-ring preserves rank and the unit
  generator for every shift.
- **OBC/no-embedding explicit errors (roadmap gate):** adapter-level
  `WindingWithoutEmbedding`, parameter-path rejection mentioning the
  embedded-periodic-topology requirement (default-open square and
  pbc-only honeycomb), and `EmbeddingError::OpenBoundaries` from the W1
  builders.
- **Physical endpoints:** bond p=1 on the 3x3 torus wraps with rank 2 every
  sample, p=0 never wraps; site p=1 on the 2x2 torus wraps; fixed seeds
  reproduce `Wrapping`/`WindingRank`/`BoundaryCrossing0` bitwise (means,
  bin counts) and a shifted seed moves the stream.
- **Full enumeration vs reference:** every bond configuration of the 2x2
  torus (2^8 = 256) and the 3x3 torus (2^18 = 262144) is
  parity-checked component-by-component (identity, size, rank) against the
  lift-BFS reference. Exact results at p = 1/2: 2x2 P(wrapping) =
  187/256 = 0.73046875, mean largest winding rank 252/256 = 0.984375; 3x3
  P(wrapping) = 186684/262144 = 0.7121429443359375, mean largest rank
  259658/262144 = 0.9905166625976563. The 2x2 hand-check pins the
  two-disjoint-wrapping-2-cycle configuration (2 components, both rank 1)
  and the two-open-paths configuration (rank 0).
- **Statistical:** 16-seed scheduler z-scores (4000 sweeps/seed, seed base
  `0x5A5A_0000`) of `Wrapping` and `WindingRank` against the exact 2x2
  enumeration — per-seed |z| < 4, mean |z| < 1.5, no one-sided bias;
  seed count scales through `SCUTTLE_ZSCORE_SEEDS` (verified green at 64).
- **Cross-consistency with F3:** 32 random 4x4 samples agree on active
  vertex/edge counts, component count, and canonical largest identity/size.
- **Typed rejections and overflow:** activity-length mismatches, zero
  dimensions, missing-displacement embeddings; a 3-cycle of 2*10^9
  displacements (cycle sum exactly zero) is reported contractible while
  3*10^9 (sum 3e9 > i32) is a typed `DisplacementOverflow`, and an 8-cycle
  of 6*10^8 displacements (sums 4.8e9) errors at the analyzer and is
  rejected by the adapter's static magnitude gate at construction. Unit
  tests additionally pin extended-gcd reference values, workspace reuse
  across repeated analyses, and dimension switches (D=1 -> D=2) on one
  workspace.

### Performance

Criterion target `percolation_winding` (sample_size 10, 100 ms warmup,
250 ms measurement; throughput = physical edges/s; recorded on the Xeon
Gold 6148 host, p_bond = 0.5, 2026-10-09):

| Case | sample+f3 | sample+winding | sample+f3+winding | adapter plain | adapter winding |
|------|-----------|----------------|-------------------|---------------|-----------------|
| ring V=4096, E=4096, D=1 | 69.95 us (58.6 Melem/s) | 137.2 us (29.9 Melem/s) | 197.4 us | 74.1 us | 205.7 us |
| square 64x64 V=4096, E=8192, D=2 | 119.1 us (68.8 Melem/s) | 270.9 us (30.2 Melem/s) | 381.9 us | 118.7 us | 387.4 us |
| cubic 16x16x8 V=2048, E=6144, D=3 | 83.4 us (73.7 Melem/s) | 220.0 us (27.9 Melem/s) | 297.3 us | 84.4 us | 296.7 us |

The winding scan costs ~2.0-2.6x the F3 scan per edge (ring 1.96x through
cubic 2.64x; potential bookkeeping, potential-aware path compression, basis
folding), so an adapter requesting both analyses runs ~2.8-3.5x the plain
adapter per sample (ring 2.78x through cubic 3.52x); the winding-only path is strictly cheaper than F3+winding because
it fuses connectivity and winding into one scan. Reproduce with:

```bash
cargo bench -p cmc-rs --bench percolation_winding -- --test   # smoke + parity
cargo bench -p cmc-rs --bench percolation_winding
```

### Allocation probe

Six new `percolation_allocations` records (3 workspace + 3 adapter, 128
post-warmup samples each; `--features allocation-probe -- --test`): the
winding hot path (law sampling + `analyze_winding`) records **0.000
allocations and 0.000 bytes per sample** in steady state after one
`WindingWorkspace::prepare`, and the winding adapter (sweep + measure, 9
scalars + no boundary query) likewise 0.000 / 0.000 with zero peak live
allocations. Retained workspace capacity: ring (D=1) 327,696 B = 80.0
B/vertex, square (D=2) 393,248 B = 96.0 B/vertex, cubic (D=3) 245,808 B =
120.0 B/vertex — exactly `4*8 + 32 + 2*4 + 4D + 4D^2` bytes per vertex
(union-find arrays + worst-case component table + counters + potential +
basis blocks); 10 allocations at `prepare`, none afterwards.

### Limitations

Generators are an echelon generating set, not a canonical Hermite normal
form (two runs generating the same lattice may report different vectors;
the tests compare lattices, not rows). The overflow check is conservative:
it can reject an adversarial table whose intermediate sums leave `i32`
even when the mathematical winding subgroup is small, and never wraps
silently. The adapter magnitude gate bounds potential sums only; a crafted
table passing the gate could in principle overflow a checked basis-echelon
combination inside `measure()` and panic loudly (builder embeddings carry
unit displacements and cannot reach this) — use the core `analyze_winding`
for adversarial tables. An embedding table with fewer edges than the graph
panics at the infallible per-edge lookup, the same F1 caller contract
(pair via `LatticeEmbedding::edge_count`); the adapter checks it for you.
`Wrapping` is defined as "any tracked component has rank > 0" and
`WindingRank` tracks the canonical largest component only — other
conventions (seed-site wrapping, per-axis wrapping probabilities) are
separate named observables that do not exist. Vector/generator output is
analyzer-level until C0; honeycomb/kagome topologies gained embeddings in
T3 (below), so `analyze_winding` covers them at the analyzer level, but the
adapter parameter schema still accepts only the chain/square/cubic/
hypercubic/triangular builders; heterogeneous laws compose with winding
through the same adapter unchanged.

## Z1 Newman-Ziff parameter scans

Status: **Validated** — independent review passed 2026-10-22 (roadmap §10
Z1). Z1 is a core scalar API: it does not add
a Carlo.rs adapter or vector `Results`; those wait for C0.

### Semantics

A scan draws one uniformly random permutation of the tracked elements and
activates it incrementally:

- **Site:** the permutation is over vertices. `points[k]` describes the first
  `k` active vertices; a physical edge is active iff both endpoints are active.
  Each physical edge is counted once (self-loops once, parallel edges
  separately), matching the F3 active-edge definition.
- **Bond:** the permutation is over physical edges. Every vertex is tracked
  from the start, so isolates count as components; a self-loop activation is
  one active edge but does not change connectivity.

`ScanPoint` records `occupation`, `active_edge_count`, `component_count`,
`largest_component_size`, and `raw_second_moment = sum(s_i^2)`. `ScanCurve`
exposes the dense activation order and exactly `trials + 1` points
(`trials = V` site, `E` bond). One curve is one realization; independent
curves are i.i.d. No mixed, heterogeneous, deletion, crossing, or wrapping
process is part of Z1.

### API

- `NewmanZiffWorkspace::prepare(V, E)` reserves the union-find, activity,
  order, and point buffers; repeated `scan_site`/`scan_bond` calls for the same
  topology allocate nothing.
- `scan_site(graph, workspace, rng) -> ScanCurve<'_>` and
  `scan_bond(graph, workspace, rng) -> ScanCurve<'_>` work on any
  `UndirectedGraphView`; IDs/order follow the F1 unchanged-view provenance
  contract.
- `ScanObservable` selects `ComponentCount`, `ActiveEdgeCount`,
  `LargestComponentSize`, `GiantFraction`, or `RawSecondMoment`.
- `reweight(curve, Probability, observable, &mut ReweightWorkspace) -> f64`
  applies the exact binomial weights. `Probability` is the validated
  `[0, 1]` boundary; no raw `f64` is accepted.
- `ScanEstimator` accumulates independent per-scan reweighted values and
  returns `mean`, unbiased sample `variance`, `standard_error`, and `samples`.
- `ScanError` is capacity-only. `ReweightError` covers capacity, curve-length
  mismatch, undefined giant fraction at `V = 0`, and the exact-f64 moment
  limit.

### Reweighting

For `n` activations, a uniformly random permutation makes the size-`k` prefix
uniform over the `C(n,k)` subsets, so the Bernoulli expectation at `p` is
`sum_k C(n,k) p^k (1-p)^(n-k) f(point_k)`. `reweight` computes the weights by
starting at the binomial mode and recursing outward; all relative weights are
at most one, the common normalization cancels in the weighted mean, and no
per-`k` `exp`/`log` is needed. `p = 0` and `p = 1` return `point_0` and
`point_n` exactly.

`GiantFraction` is `largest_size / V` with the F4 total-topology denominator;
it is undefined at `V = 0`. `RawSecondMoment` is exact in `u128` on a point;
reweighting converts to `f64` and rejects `V > 94,906,265` (`V^2 > 2^53`) with
`ReweightError::InexactMomentConversion`, matching F4. The per-scan
`reweight` value is deterministic given the curve; the estimator's variance
and standard error describe independent-scan Monte Carlo uncertainty, not a
binomial per-`k` variance.

### Validation gates

`tests/physics/percolation_scan.rs` (6 tests) plus 11 unit tests
(`newman_ziff.rs` 5, `reweight.rs` 6):

- **Every occupation-count step vs naive recomputation:** for chain 8,
  4x4 square, 2x2x2 cubic, 3x3 triangular, 2x2 honeycomb, 2x2 kagome, a seeded
  random graph, and a multigraph with self-loop/parallel edges, four seeds,
  every prefix is compared field-by-field (`active_edge_count`,
  `component_count`, `largest`, `M2`, and the mode-specific active count) to a
  fresh F3 `analyze` over the same occupancy. F3 is itself independently
  flood-fill validated.
- **Exact enumeration and z-scores:** full 3x3 enumeration (site 2^9, bond
  2^12) at `p = 0.3, 0.5, 0.7` for all five observables. `R = 16 *
  zscore_seed_count(16)` independent scans (default 256,
  `SCUTTLE_ZSCORE_SEEDS`-scalable) feed `ScanEstimator`; each gate requires
  `|z| < 4.5` against the exact mean, with an exact equality path for the
  deterministic bond `ActiveEdgeCount` and a sign-balance check. Verified at
  `SCUTTLE_ZSCORE_SEEDS=64`.
- **Direct Bernoulli:** independent Bernoulli occupancy + F3 analysis at
  `p = 0.55` (8192 samples) agrees with the reweighted estimator within a
  combined 4-sigma interval for `ComponentCount` and `LargestComponentSize`.
- **Endpoints:** `p = 0/1` equal the first/last point exactly for every
  observable; an empty borrowed graph returns one empty point and rejects
  `GiantFraction`; fixed seeds reproduce order and points bitwise; a shifted
  seed changes the order; site/bond reuse of one workspace matches fresh
  workspaces.

### Performance

Criterion target `percolation_newman_ziff`, recorded 2026-10-21 on the Xeon
Gold 6148 / rustc 1.98.0 host, short-run configuration (10 samples, 100 ms
warm-up, 250 ms measurement). The comparison uses 17 probabilities and 64
independent replicas for both paths; "direct" is 64 independent Bernoulli
draws + F3 `analyze` at each `p`, and "combined" is 64 independent scans +
17 reweights per scan. Times are point estimates, not CI thresholds.

| Case | Scan | 17p reweight | Combined 64r | Direct 64r | Speedup |
|---|---:|---:|---:|---:|---:|
| square site (V=4096, E=8064) | 217.6 us | 278.8 us | 32.55 ms | 79.97 ms | 2.46x |
| square bond | 207.8 us | 520.3 us | 47.88 ms | 101.03 ms | 2.11x |
| cubic site (V=4096, E=11520) | 256.4 us | 279.8 us | 35.16 ms | 98.78 ms | 2.81x |
| cubic bond | 256.2 us | 724.0 us | 64.97 ms | 129.95 ms | 2.00x |
| sparse-ER site (V=4096, E=16384) | 321.7 us | 292.3 us | 38.69 ms | 179.79 ms | 4.65x |
| sparse-ER bond | 339.9 us | 1416.6 us | 112.12 ms | 205.59 ms | 1.83x |

The advantage grows with edge density for site scans and is strongest when the
same curve is reused across many probabilities; bond scans pay the larger
`n = E` reweight. The scan itself is a single incremental union-find pass; the
reweight is a weighted sum over the curve.

The `percolation_allocations` probe adds 10 records (5 topologies x site/bond,
128 post-warmup scan+reweight iterations). All report **0.000 allocations,
0.000 allocated bytes per sample, and zero peak live allocations/bytes** after
one `prepare` of each workspace. Retained capacity (scan workspace + reweight
weights) is:

| Topology | Mode | Total bytes | B/V | B/(V+E) |
|---|---:|---:|---:|---:|
| chain (V=4096, E=4095) | site | 335,927 | 82.0 | 41.0 |
| chain | bond | 335,919 | 82.0 | 41.0 |
| square (V=4096, E=8064) | site | 562,104 | 137.2 | 46.2 |
| square | bond | 593,848 | 145.0 | 48.8 |
| cubic (V=4096, E=11520) | site | 759,096 | 185.3 | 48.6 |
| cubic | bond | 818,488 | 199.8 | 52.4 |
| sparse-ER (V=4096, E=16384) | site | 1,036,344 | 253.0 | 50.6 |
| sparse-ER | bond | 1,134,648 | 277.0 | 55.4 |
| power-law (V=4096, E=16374) | site | 1,035,774 | 252.9 | 50.6 |
| power-law | bond | 1,133,998 | 276.9 | 55.4 |

The table above excerpts 8 of the 10 probe records (the power-law rows
are omitted for space; they are part of the probe output). The point
buffer (48 bytes per point) dominates; capacity excludes struct headers,
allocator metadata, and RSS effects. Reproduce with:

```bash
cargo test -p cmc-rs --test suite percolation_scan
SCUTTLE_ZSCORE_SEEDS=64 cargo test --release -p cmc-rs --test suite percolation_scan
cargo bench -p cmc-rs --bench percolation_newman_ziff -- --test
cargo bench -p cmc-rs --bench percolation_newman_ziff
cargo bench -p cmc-rs --bench percolation_allocations --features allocation-probe -- --test
```

### Limitations

Z1 is a core uniform-site/bond scan only. Its scan layer has no heterogeneous
probabilities, mixed `(p_site, p_bond)` surface, deletion process, boundary
crossing, wrapping, or targeted/network facade — heterogeneous probabilities
enter as a reweight-layer extension (T2 below) over the unchanged uniform
scan, never as a changed activation process. There is no Carlo.rs
`MonteCarlo`/`FromParams` adapter and no vector `Results`/`ArrayEstimate`
output until C0; callers consume the borrowed `ScanCurve` directly. The
`ScanEstimator` variance is scan-to-scan only. The benchmark's direct path is
raw Bernoulli draws plus F3 analysis, not the full law/scheduler stack, so the
speedup is indicative rather than a scheduler-level claim.

### T2: heterogeneous (Poisson-binomial) reweighting

Status: **Experimental** (roadmap v2 §2 T2; no independent review yet).

`reweight_heterogeneous(curve, field, observable, workspace)` reweights one
Z1 scan curve — `scan_site`/`scan_bond` are untouched — to a heterogeneous
Bernoulli field with one validated `Probability` per universe element
(vertices for site curves, physical edges for bond curves). Because a
uniformly random permutation's size-`n` prefix is a uniform `n`-subset,

```text
E_het[f] = sum_n C(N,n) * E_unif,n[ w * f ],
w(A)     = prod_{i in A} p_i * prod_{j not in A} (1 - p_j),
```

so the per-scan raw estimator `sum_n C(N,n) w(A_n) f_n` is unbiased by
linearity, and with `f = 1` it becomes `sum_n C(N,n) w(A_n)` whose
expectation is exactly the total Poisson-binomial mass, 1 — returned as
`HetReweightResult::unit_mass`, a self-diagnostic that must average to 1 over
scans. The returned value is deliberately the **raw** sum, not a weighted
average: the normalization fluctuates per scan under a heterogeneous field,
and dividing by it would introduce a ratio-estimator bias (the uniform
`reweight` normalizes only because its weights depend on `n` alone; uniform
`ProbabilityField`s therefore delegate to it outright with `unit_mass = 1`).

All arithmetic is log-space: `l_0 = sum_j ln(1 - p_j)`, entering element `e`
adds `ln(p_e/(1 - p_e))` along `activation_order()`, and `log C(N,n)` grows
by the `O(N)` recurrence from `logC(N,0) = 0` (an `ln(k)` table in
`HetReweightWorkspace`, no `lgamma` or new dependency). A running-max
log-sum-exp keeps every internal term in `(0, 1]` — positive terms, so no
catastrophic cancellation, and `N = 16384` mid-range fields stay finite. The
final `exp(l_max)` restores the linear scale; outside the envelope below the
reachable failure is underflow — a near-degenerate field (e.g. half 0.01 /
half 0.99 at `N = 10^4`) pushes the best single-scan log-mass below the
`f64` range and `value`/`unit_mass` return exactly `0.0` — while the mirror
side, saturation past `ln(f64::MAX) ~ 709.8`, needs a positive excursion
that large and has not been observed even at `N = 2*10^5` (review-verified).
Fields with entries `0.0`/`1.0` are typed rejections
carrying the index — boundary fields are the N1 direct-sampling path, and no
universe-splitting is attempted.

**Variance envelope (the central caveat).** The per-scan values are
heavy-tailed once the field spread times the universe size grows: the
log-weight of a size-`N/2` prefix fluctuates across scans with standard
deviation about `sqrt((N/4) * Var(logit p))`, so scan values behave roughly
lognormally with that parameter. Measured coefficient of variation of the
per-scan giant-fraction estimate on a column ramp field around `p_c` at
`L = 64` (`N = 4096`): 0.29 at zero width (the uniform limit), 1.25 at
half-width 0.012, 3.3 at 0.02, 45 at 0.09 — a bimodal half-0.2/half-0.9
field already gives cv 2.7 at `N = 8` and is unusable at `N = 4096`
(`Var(logit p)` is 3.21 for that bimodal field vs 0.047 for the
half-width-0.09 ramp — about 68x — and the envelope parameter
`sqrt((N/4) Var(logit p))` at `N = 4096` is 57 vs 7). Rule of thumb: the
reweight path is competitive while
`sqrt((N/4) Var(logit p))` stays of order one or below; beyond that use the
N1 direct sampler, which keeps light-tailed estimates at any spread. The
statistical gates below are parameterized inside that envelope by design.

### T2 validation gates

`tests/physics/percolation_scan_heterogeneous.rs` (5 tests) plus 6 unit
tests in `reweight.rs` (hand-worked two-order weights, uniform delegation,
constant-field reduction, typed rejections, workspace reuse/growth, empty
curve):

- **Exact enumeration (site N=8 seeded random graph, bond E=7 chain, R = 200
  * `zscore_seed_count(16)` = 3200 default scans):** full 2^8/2^7 enumeration
  under per-element weights with a test-local BFS reference (never production
  union-find). Fields: deterministic U(0.3,0.8) spread and bimodal
  half-0.2/half-0.9. All eight z-scores (two fields x two modes x
  GiantFraction/RawSecondMoment) within |z| <= 4; recorded max |z| = 1.35
  (bimodal site RawSecondMoment). Unit-mass z at the same settings:
  -1.03..+1.32.
- **Unit mass and finiteness:** the L=128 ramp field below at R = 400 gives
  unit-mass mean 1.000539 +/- 0.015284 (z = +0.04 vs exactly 1); every
  returned `value` and `unit_mass` is finite at `N = 10^4` (100x100 square)
  for spread, bimodal, and constant fields over 4 seeds each.
- **N1 cross-validation at scale (L = 64 and L = 128 open squares, site
  scans):** clipped linear ramp fields in the column coordinate centered at
  `p_c = 0.5927` with half-widths 0.004 (L=64) and 0.002 (L=128) — chosen
  inside the variance envelope above. Reweighted (R = 400, seed bases
  0x4845_5a32_0000 / 0x4845_5a33_0000) vs direct `SiteBernoulli` at the same
  field (4096 samples, seed 0x4449_5200+side): L=64 z = -1.21 (giant
  fraction, 0.27751 +/- 0.00608 vs 0.28509 +/- 0.00162) and -0.80 (raw
  second moment); L=128 z = -0.65 and -0.68. All |z| <= 4.
- **Uniform reduction:** an owned constant-valued field agrees with the v1
  binomial `reweight` within 1e-9 relative across p in {0.1, 0.25, 0.5,
  0.5927, 0.75, 0.9} on 4x4 square site and bond curves (the `Uniform`
  variant delegates, so it matches bitwise).
- **Typed rejections:** field length != trials, `p_i = 0` and `p_i = 1` with
  the offending index in the error; NaN is unconstructible through
  `Probability`.
- **Allocation probe:** `HET_SCAN_ALLOC` records (5 topologies x site/bond,
  128 post-warmup scan+reweight iterations) all assert **0 allocations and 0
  bytes per realization**; retained capacity is two `f64` words per element
  (log weights + ln table), e.g. 65,544 bytes site / 129,032 bytes bond on
  the 4096-vertex square.

### T2 performance

Criterion group `nz_het` (same short-run configuration as the Z1 table,
2026-10-22): one scan then reweighting at K = 8 graded fields vs 8 direct
heterogeneous sampling runs (per-element draws + F3 analysis) at those
fields, 64 replicas per iteration, `LargestComponentSize`:

| Case | Scan + 8 het reweights | Direct, 8 fields | Direct / reweight |
|---|---:|---:|---:|
| square site (V=4096, E=8064) | 88.7 ms | 56.3 ms | 0.63x |
| square bond (E=8064) | 152.4 ms | 62.7 ms | 0.41x |
| cubic site (V=4096, E=11520) | 87.1 ms | 69.4 ms | 0.80x |
| cubic bond (E=11520) | 215.2 ms | 81.4 ms | 0.38x |
| sparse-ER site (V=4096, E=16384) | 93.9 ms | 131.6 ms | **1.40x** |
| sparse-ER bond (E=16384) | 318.7 ms | 137.7 ms | 0.43x |

Unlike the uniform Z1 reweight (17 probabilities, 1.8-4.7x speedups), the
heterogeneous path pays two `ln` evaluations per universe element per call,
so at K = 8 fields it only wins where the direct path is expensive relative
to the scan universe — site scans on the sparsest topology. The advantage
grows with the number of distinct fields queried per scan (the ln table is
per workspace, and the per-call cost is field-dependent), and the variance
envelope above — not timing — is what limits where the method applies.
Times are point estimates, not CI thresholds.

### T2 limitations

One field per call, one observable per call (the curve is reusable). The
estimator's variance envelope is the binding constraint (see the variance
caveat); the unit-mass diagnostic detects envelope violations empirically
(scan averages drifting from 1 or fluctuating wildly), but no automatic
fallback to the N1 path exists — the caller chooses. No Carlo.rs adapter.
Independent review passed 2026-10-22
(`docs/reviews/2026-10-22-t2-het-reweight-review.md`): the identity
brute-forced over all permutations in exact rational arithmetic, the
ratio-bias finding independently confirmed (with v1 shown structurally
immune), a from-scratch sampler/union-find cross-check at L=64
(z = -0.18/-0.26), and the saturation story completed (the reachable
out-of-envelope failure is underflow to zero, not the +709.8 saturation).

## N2 random network failure and degree observables

Random failure removes vertices independently with keep probability `q`
(`RandomFailureMC`, `percolation/network/`). It is a domain facade over the
F2 uniform site law — `sweep()` delegates to `SiteBernoulli`, so reduction to
`StaticLaw::site(q)` is by construction and bitwise under one RNG stream —
never a duplicate sampler. Scalar component observables reuse the F4
[`ObservablePlan`] verbatim. Degree conditioning is always on **original**
degree cached once in `NetworkProfile::from_graph` (self-loops contribute two
incidences); residual-degree conditioning is a different process and is not
offered here.

### Degree arrays (C0 element-preserving output)

Recorded per sample with shape `[max_degree + 1]`, indexed by original degree
`k = 0 ..= max_degree`:

| Name | Definition |
|---|---|
| `DegreeVertexCount[k]` | `n_k`: vertices of original degree `k` (constant) |
| `DegreeOccupied[k]` | 1 when `n_k > 0`, else 0 (constant) |
| `ActiveFractionByDegree[k]` | active degree-`k` vertices over `n_k` |
| `GiantMembershipByDegree[k]` | degree-`k` vertices in the canonical largest component over `n_k` (F4 tie policy) |
| `DegreeFirstMoment[k]` | `sum_{active v, deg k} s(v) / n_k` |
| `DegreeSecondMoment[k]` | `sum_{active v, deg k} s(v)^2 / n_k` |

`s(v)` is the size of `v`'s component. The conditional mean component size of
a degree-`k` vertex — `S(k)` in Cohen et al. — is estimated by the ratio of
ensemble means `DegreeFirstMoment / ActiveFractionByDegree`, never by a
per-sample ratio (undefined when no degree-`k` vertex is active).

**Empty-bin semantics**: degrees absent from the graph (`n_k = 0`) keep the
full shape and record exactly `0.0` (never NaN, never skipped) in the four
ratio arrays; `DegreeOccupied[k] == 0` and `DegreeVertexCount[k] == 0`
identify them. Denominators are the constant `n_k`, so present bins are always
defined — there is no per-sample undefined case.

**Exactness gate**: degree second moments sum `s(v)^2` over vertices, bounded
by `sum_c s_c^3 <= V^3`; the constructor rejects `V^3 > 2^53`
(`V > 208_063`) exactly like the F4 moment gate rejects `V^2 > 2^53`.
Winding observables are rejected (`WindingUnsupported`) — the facade carries
no embedding.

### Validation gates (roadmap §9 N2)

All gates in `tests/physics/percolation_network.rs`:

- **Site-law reduction**: 32 sweeps of `RandomFailureMC` produce
  bitwise-identical vertex activity to `StaticPercolationMC::new(..,
  StaticLaw::site(q), ..)` and to a raw `SiteBernoulli` law on
  identically-seeded streams (star, `q = 0.61`).
- **Hand-graph exact expectations**: path5, cycle6, star6, complete4 at
  `keep in {0.3, 0.5, 0.7}` against exhaustive `2^V` subset enumeration whose
  per-subset quantities come from an **independent adjacency-BFS reference**
  (canonical-lowest largest-component identity, F4 susceptibility exclusion).
  Susceptibility is compared as a conditional mean over defined samples
  (`FiniteClusterSusceptibilityDefined` convention); a case where every
  sample is undefined skips that gate loudly. 16 seeds x 600 sweeps
  (verified at 64 seeds via `SCUTTLE_ZSCORE_SEEDS`), |z| < 4 per observable,
  |mean z| < 2 per (graph, keep): first run max |z| = 2.02, all mean |z| < 0.6.
- **Relabeling invariance**: star10 under the reversed labeling (statistical,
  2000 sweeps): constant arrays (`DegreeVertexCount`, `DegreeOccupied`)
  identical bitwise; every ratio array element and `GiantFraction` within
  4 combined sigma.
- **Empty bins**: star5 (degrees {4,1,1,1,1}) — full `[0..=4]` shape, bins
  {0,2,3} exactly zero with zero stderr, occupancy markers correct.
- **Literature anchors** (ER G(n,m) configuration graph, N = 100 000,
  z = 2.5, empirical degree distribution):
  - **Cohen et al., PRL 85, 4626 (2000)** (`docs/ref/Cohen_2000_PRL_internet_resilience.pdf`,
    arXiv:cond-mat/0007048): measured `kappa = <k^2>/<k> = 3.5025` gives
    `q_c = 1/(kappa-1) = 0.3996` against the ideal Poisson `0.4` (agreement
    0.0004). Bracket: at `q_c - 0.03` the giant is extinguished
    (`S = 0.0030`); at `q_c + 0.10` a giant exists (`S = 0.1851 ± 0.0001` vs theory 0.1855). Both points sit outside the finite-size rounding window
    `~N^{-1/3} ~ 0.02`.
  - **Callaway, Newman, Strogatz, Watts, PRL 85, 5468 (2000)**
    (`docs/ref/Callaway_2000_PRL_network_robustness.pdf`,
    arXiv:cond-mat/0007300): giant fraction `S = q(1 - G0(u))` with
    `u = 1 - q + q G1(u)` (their Eq. 5) evaluated on the graph's own
    empirical degree distribution. Measured vs theory:
    `q=0.5`: 0.1857(1) vs 0.1861; `q=0.65`: 0.4261(1) vs 0.4259;
    `q=0.8`: 0.6377(1) vs 0.6373 — every deviation `< 4 sigma + 0.006`
    (actual deviations <= 0.0005).
  - The platform's `FiniteClusterSusceptibility` (largest-excluded,
    active-denominator) is deliberately **not** compared to Callaway et al.'s
    `H0'(1)` mean cluster size, which uses a different convention
    (unoccupied vertices contribute zero, no giant exclusion); conflating the
    two would be a definitional error, not a validation.

Reproduce with:

```bash
cargo test -p cmc-rs --test suite percolation_network
SCUTTLE_ZSCORE_SEEDS=64 cargo test --release -p cmc-rs --test suite percolation_network
cargo bench -p cmc-rs --bench percolation_network
cargo bench -p cmc-rs --bench percolation_allocations --features allocation-probe
```

At 64 seeds every (graph, keep) gate stays inside |z| < 2 (max 1.98) with
|mean z| <= 0.81.

### Steady-state allocation gate

`PERCOLATION_RANDOM_FAILURE_ALLOC` (allocation-probe bench): after warmup,
`sweep + measure` (labels pass plus six degree arrays) performs **0**
allocations per sample on all five benchmark topologies; scratch
(`root_size`, degree counters, shared value buffer) is fully preallocated and
refilled in place.

### Criterion data

`cargo bench -p cmc-rs --bench percolation_network` — 64 sweeps per iteration
at `q = 0.5927` (V = 4096 for every topology; the facade also records the full
`ObservablePlan::all` scalar set, the baseline records two scalars, so the
ratio prices labels + degree arrays + scalars together):

| Topology | `random_failure_degree` | `scalar_baseline` | Ratio |
|---|---|---|---|
| chain (E=4095) | 12.4 ms | 5.3 ms | 2.33x |
| square (E=8064) | 12.7 ms | 7.4 ms | 1.71x |
| cubic (E=11520) | 14.0 ms | 9.8 ms | 1.43x |
| sparse-er (E=16384) | 22.9 ms | 17.7 ms | 1.29x |
| power-law (E=16374) | 19.8 ms | 11.5 ms | 1.72x |

The label materialization pass dominates the added cost; degree counting and
the six array records are two more vertex-order passes over preallocated
buffers.

### Limitations

N2 is a facade, not a process: it draws independent failure realizations and
has no removal trajectory (that is N3, which now unblocks on C0 and records
one vector realization per attack). There is no `FromParams` parameter schema
and no boundary-query defaults; construct `RandomFailureMC::new` directly or
through `Run::from_parts`. `DegreeSecondMoment` bounds graphs to
`V <= 208_063` for exact f64 conversion. Adaptive/residual-degree
conditioning is explicitly out of scope.

## N3 targeted degree attack and robustness curves

`DegreeAttackMC` (`percolation/network/attack.rs`) removes vertices in
descending **original** degree (the `NetworkProfile` cache). Ties resolve by
an explicit policy: `DegreeDescendingStableId` (ascending dense ID,
deterministic and seed-independent) or `DegreeDescendingSeededShuffle`
(Fisher-Yates over each equal-degree run, one realization per seed, bitwise
reproducible for a fixed seed). Adaptive residual-degree attack is a
different process and lives in its own T1 section below, on the same
adapter.

### One trajectory is one vector realization

One `sweep()` samples the removal order (tie randomness only); one
`measure()` walks the trajectory with a fresh F3 component analysis at every
checkpoint and records, per trajectory (never per checkpoint as independent
scalars):

| Name | Definition |
|---|---|
| `RobustnessCurve[i]` | `S_max / V` at checkpoint `i` (F4 `GiantFraction` convention, fixed denominator `V`) |
| `RobustnessSusceptibility[i]` | F4 finite-cluster susceptibility at checkpoint `i`; recorded only when defined at every checkpoint of the trajectory |
| `RobustnessSusceptibilityDefined` | scalar 0/1 indicator, always recorded |
| `RobustnessRemovedFraction[i]` | removed fraction `r_i / V`, constant coordinates |

Grids: `Full` (V+1 points, `O(V(V+E))` per trajectory), `Every { step }`
(always ending at `V`), `At { removals }` (user-exact, so a grid can stay
inside the susceptibility-defined region). `V^2 <= 2^53` is enforced
(susceptibility squares component sizes). Provenance: callers stamp
`results.set_attribute("attack_policy", mc.policy().as_provenance())` and
the grid; Carlo.rs stores and round-trips the values without interpreting
them.

### Validation gates (roadmap §9 N3)

All in `tests/physics/percolation_attack.rs`:

- **Per-step naive parity**: full-resolution trajectories on star5, path6,
  and a 12-vertex ER graph, both policies: every checkpoint's giant fraction
  and susceptibility equal an independent adjacency-BFS reference exactly
  (no statistics), and the survivor count is strictly monotone (recomputed
  independently from the survivor mask).
- **Tie reproduction**: the stable-ID order is seed-independent (and equals
  the hand-derived path8 order); seeded shuffle reproduces bitwise for a
  fixed seed, reshuffles tie groups for new seeds, and preserves the
  degree-descending group structure.
- **Coordinates and provenance**: `RobustnessRemovedFraction` is exactly
  `i * step / V`; full-grid trajectories hit undefined-susceptibility
  checkpoints (intact connected graph, empty end) so only the indicator is
  recorded; policy attributes survive `Results` JSON round-trip.
- **Literature anchor (exact)** — Callaway et al., PRL 85, 5468 (2000),
  Eq. (14): checkpoints at complete degree-group boundaries (removing
  exactly all vertices with degree > K) on a 60k-vertex ER configuration
  graph (z = 2.5) and a 60k power-law configuration graph
  (alpha = 2.5, m = 2, cutoff 400); theory `S = F0(1) - F0(u)`,
  `u = 1 - F1(1) + F1(u)` with `q_k = 1{k <= K}` evaluated on each graph's
  own degree distribution. Measured vs theory (deterministic stable-ID
  trajectories; deviations well inside the 0.006 finite-N slack):
  ER z = 2.5 — K=7: 0.8856/0.8851, K=6: 0.8641/0.8627, K=5: 0.7949/0.7954,
  K=4: 0.5671/0.5690 (max deviation 0.0019); power-law alpha = 2.5,
  m = 2, cutoff 400 — K=80: 0.9899/0.9898, K=40: 0.9694/0.9694,
  K=20: 0.9167/0.9175, K=10: 0.7678/0.7707 (max deviation 0.0029).
- **Literature anchor (bracket)** — Cohen et al., PRL 86, 3682 (2001):
  their Eqs. (2)/(7)/(9) attacked-distribution estimate with the
  `kappa = 2` criterion gives the breakdown fraction `f_c`. The discrete-
  empirical estimate on our finite N = 6e4 power-law graph lands at
  `f_c = 0.0211` versus `f_c = 0.2428` for ER z = 2.5 (the paper's own
  continuum equations give ~0.056 at alpha = 2.5 — the finite realized
  degree cutoff lowers our estimate; the review's independent solution of
  their Eqs. (8)+(11) established this, so the 0.02-0.03 range belongs to
  alpha ~ 3, not 2.5). Gates: `f_c(power-law) < 0.5 f_c(ER)`
  (0.087 ratio) and the measured giant survives at `0.3 f_c`
  (`S = 0.346 > 0.25`) but has collapsed at `1.5 f_c`
  (`S = 0.0016 < 0.05`). A statistical bracket, never an exact threshold
  claim (the estimate carries its own continuum approximations).

### Steady-state allocation gate

`PERCOLATION_DEGREE_ATTACK_ALLOC` (allocation-probe bench): after warmup,
one complete trajectory (order construction + 17 checkpoint analyses) is
allocation-free on all five benchmark topologies (0.000/trajectory).

### Criterion data

`cargo bench -p cmc-rs --bench percolation_attack` — 8 trajectories x ~17
checkpoints per iteration (V = 4096):

| Topology | stable-id | seeded-shuffle |
|---|---|---|
| chain (E=4095) | 5.30 ms | 9.06 ms |
| square (E=8064) | 7.56 ms | 11.23 ms |
| cubic (E=11520) | 9.37 ms | 13.24 ms |
| sparse-er (E=16384) | 23.8 ms | 25.4 ms |
| power-law (E=16374) | 13.9 ms | 15.0 ms |

The shuffle premium is the tie-group Fisher-Yates (largest on the chain,
whose 4094 interior vertices form one tie group).

### Limitations

Original-degree order only within THIS section: adaptive residual-degree
attack (recomputed during the attack) is the T1 section below. The
trajectory analysis is checkpoint-parallel-naive: every checkpoint is a
fresh exact analysis (no incremental removal structures), so `Full` grids
cost `O(V(V+E))` per trajectory; production curves should use `Every`/`At`.
No boundary queries, no degree-conditioned arrays (that is N2's
`RandomFailureMC`), no `FromParams` schema — construct directly or via
`Run::from_parts`.


## T1 adaptive residual-degree attack

Roadmap v2 §2 T1 (the deferred N3 tail). Two `AttackPolicy` variants extend
the SAME `DegreeAttackMC` adapter — grids, C0 vector output
(`RobustnessCurve` / `RobustnessSusceptibility` with the defined-indicator
convention / `RobustnessRemovedFraction` coordinates), and provenance
(`as_provenance()`: `residual-degree-descending/stable-id`,
`residual-degree-descending/seeded-key`) are exactly the N3 conventions:

- `ResidualDegreeDescendingStableId`: at every removal the victim is the
  maximum **residual** degree — degree within the still-active subgraph, a
  self-loop counting two incidences at its own vertex while it is active —
  with ties by ascending dense vertex ID. Deterministic, no RNG draws.
- `ResidualDegreeDescendingSeededKey`: same ordering, ties by ONE uniform
  `u64` tie-break key per vertex drawn at trajectory start (draw contract:
  ascending dense vertex order, one `u64` each, before the first removal).

**Tie-policy semantics.** For the original-degree policies the seeded rule
is a uniform shuffle of each fixed equal-degree run (computable offline).
The adaptive order emerges ONLINE — the equal-residual-degree groups change
after every removal — so a per-run shuffle is not well-defined; the fixed
per-vertex random tie-key is the standard online tie-break. This is the
roadmap's "tie policy explicit" requirement.

**Cost model.** Order generation uses a lazy max-heap keyed
`(residual degree desc, tie key asc)` with stale-entry skipping: a residual
degree only decreases, so every neighbor decrement re-pushes a fresh entry
and pops skip entries whose degree no longer matches the tracked residual
(or whose vertex is already removed). Every active vertex always has one
fresh entry in the heap, so the popped maximum dominates every active
vertex's current entry and the first fresh pop is exactly the argmax under
the key. Pushes are bounded by `V + E` (each non-loop physical edge
decrements exactly one endpoint once), giving `O((V + E) log V)` per
trajectory vs `O(V (V + E))` for the naive full rescan (every step
recompute every active residual degree from adjacency). The naive variant
is a test/bench reference only, never production code.

### Validation gates (roadmap v2 §2 T1)

All in `tests/physics/percolation_attack.rs`:

- **Naive parity (exact)**: on chain10, cycle12, star6, a dense 10-vertex
  ER graph, a hand-built multigraph with parallel edges and a self-loop,
  and two small configuration-model multigraphs — both tie policies, three
  seeds for the seeded policy — the heap-based removal order equals the
  naive full-rescan reference (same tie keys; the seeded keys reproduced
  from the draw contract) at EVERY step: 28 cases of exact full-order
  equality. Each trajectory additionally matches an independent per-step
  adjacency BFS on giant fraction and susceptibility at every checkpoint.
- **Monotonicity/invariants**: per step, recomputed independently from the
  adjacency: active count decreases by exactly one; the removed vertex's
  residual degree equals the scanned maximum; active edges equal exactly
  the edges with both endpoints active (removed edges exactly the incident
  ones).
- **Determinism**: stable order is seed-independent and equals the naive
  rescan under identity keys; seeded order reproduces bitwise per seed,
  differs across seeds only through the drawn keys (each order equals the
  naive rescan under its own replayed keys), and same seed reproduces
  identical curves.
- **N3 compatibility**: the original-degree behavior is untouched — the
  five N3 tests above are the regression and pass unchanged.

### Literature bracket (Holme et al. 2002, Fig. 7(e))

Anchor: Holme, Kim, Yoon, Han, PRE 65, 056109 (2002) (archived PDF +
markdown in `docs/ref/`). Their Fig. 7(e) configuration: BA m0 = 5, m = 3
(here N = 10^4 per the roadmap, E = 29 995, ⟨k⟩ ≈ 6; the paper itself uses
N = 1500). Statistical brackets only — never exact claims. R = 24 seeded
graph/trajectory pairs, one trajectory per process: RD = adaptive
recalculated degree (seeded-key ties), ID = initial degree (the N3 order,
seeded-shuffle ties), random failure = uniformly shuffled removal order.

| f | S_RD | S_ID | S_fail | RD/ID CIs disjoint |
|---|---|---|---|---|
| 0.005 | 0.9945 ± 0.0001 | 0.9945 ± 0.0001 | 0.9950 | no |
| 0.01 | 0.9884 ± 0.0002 | 0.9884 ± 0.0001 | 0.9900 | no |
| 0.02 | 0.9752 ± 0.0003 | 0.9752 ± 0.0003 | 0.9800 | no |
| 0.03 | 0.9613 ± 0.0004 | 0.9612 ± 0.0005 | 0.9700 | no |
| 0.05 | 0.9304 ± 0.0007 | 0.9305 ± 0.0007 | 0.9500 | no |
| 0.075 | 0.8868 ± 0.0010 | 0.8871 ± 0.0010 | 0.9249 | no |
| 0.10 | 0.8364 ± 0.0014 | 0.8380 ± 0.0014 | 0.8997 | no |
| 0.125 | 0.7793 ± 0.0015 | 0.7818 ± 0.0017 | 0.8743 | no |
| 0.15 | 0.7092 ± 0.0023 | 0.7174 ± 0.0022 | 0.8487 | yes |
| 0.20 | 0.4814 ± 0.0053 | 0.5371 ± 0.0045 | 0.7971 | yes |
| 0.25 | 0.0116 ± 0.0018 | 0.2189 ± 0.0114 | 0.7446 | yes |
| 0.30 | 0.0014 ± 0.0001 | 0.0072 ± 0.0012 | 0.6905 | yes |
| 0.35 | 0.0008 ± 0.0001 | 0.0025 ± 0.0003 | 0.6354 | yes |
| 0.40 | 0.0002 ± 0.0000 | 0.0010 ± 0.0001 | 0.5787 | yes |

Readings (bracket language, matching Holme's claims, no exact agreement):

- Early stage (f <= 0.125): RD ~ ID — 95% CIs overlap everywhere
  ("the differences among the removal procedures are not significant in
  the early attack stage", their Fig. 7(e) discussion).
- Late stage (f >= 0.15): S_RD < S_ID with non-overlapping 95% CIs
  (empirical standard-error intervals, `mean ± 1.96 s/sqrt(R)`), the
  separation growing to ~0.21 at f = 0.25 — recalculated degrees are the
  more harmful strategy.
- Both attack curves reach S ~ 0 (upper CI < 0.05) by f = 0.30 while
  random failure still sits at S = 0.69 — both attacks collapse the giant
  at a far smaller removal fraction than random failure.
- CI convention: the intervals gated on are empirical-SE intervals, which
  respect the actual across-graph variance of S(f). A literal
  binomial-Wilson interval with n = R treats each graph's fraction as a
  Bernoulli draw and overstates the variance by ~two orders of magnitude
  (its half-width is ~0.14 at these means); it is computed and reported in
  the test for transparency, never gated on.
- Secondary context (N3 Cohen kappa machinery) applies unchanged to
  adaptive orders; their WS-model crossover (RD least harmful near
  N_rm/N ~ 0.25) is out of scope — no WS generator on the platform.

### Steady-state allocation gate

`RESIDUAL_ATTACK_ALLOC` (allocation-probe bench, asserted like the
kcore/bootstrap probes): after warmup, one complete adaptive trajectory
(seeded-key draw + heap order generation + 17-checkpoint walk) performs 0
allocations and 0 allocated bytes on all five benchmark topologies
(chain/square/cubic/sparse-er/power-law, V = 4096).

### Criterion data

`cargo bench -p cmc-rs --bench percolation_attack -- percolation_attack_order`
(order generation only, BA N = 10^4, m0 = 5, m = 3, 8 trajectories per
iteration, stable-identity naive reference):

| Implementation | Time/iteration | Time/trajectory |
|---|---|---|
| lazy heap (production) | 46.7 ms | 5.84 ms |
| naive full rescan (reference) | 4.29 s | 536.7 ms |

Measured ratio ~92x at this size — the documented per-step cost model
(`O((V+E) log V)` vs `O(V(V+E))`) realized.

### Limitations

Independent review passed 2026-10-22
(`docs/reviews/2026-10-22-t1-attack-review.md`): 1492-case adversarial
fuzz (degenerates, self-loops, parallel edges, isolated vertices) against a
from-scratch reference with zero mismatches, the heap soundness argument
re-derived, the tie draw contract replayed bit-exactly, N3 orders shown
bit-identical to 53b2019, and the heap re-timed at 128x. Review follow-ups
applied: the BA builder's complete m0-core is now documented as a standard
implementation deviation (Holme Sec. IV E starts from m0 vertices and NO
edges), the early-overlap gate now asserts every grid point f <= 0.125
(not just the first), and the random-failure contrast is separated from
the Holme bracket claims. Order generation is interleaved
into `sweep()` (the complete order is materialized before the walk; the
trajectory walk itself is the unchanged N3 checkpoint analysis). Same
`O(V (V+E))` full-grid caveat as N3. Only degree-based adaptive strategies
(recalculated betweenness — their RB — is out of scope; it needs all-pairs
shortest paths).


## Network edge-list input

`csr_from_edge_list(spec, text)` (`percolation/network/input.rs`): one edge
per line as two whitespace-separated `usize` endpoints; `#` comments and
blank lines are skipped; physical edges are emitted in file order (CSR edge
IDs follow input order). The contract is explicit in `EdgeListSpec` —
declared vertex count (at least 1 — zero is a typed
`InvalidVertexCount`), `DuplicateEdgePolicy` (`Reject` errors on the first
repeated undirected pair, self-loops included; `AllowParallel` keeps
repetitions as distinct physical edges), `SelfLoopPolicy`
(`Reject`/`Allow`; a kept self-loop contributes two incidences, the CSR
convention). Every rejection is typed (`NetworkInputError`) with the
offending 1-based line number: `Parse` (token count, non-integer),
`Endpoint` (range violation with the value), `DuplicateEdge`, `SelfLoop`.
Nothing is silently dropped or merged.

Zero-copy ingestion of complete external CSR arrays remains the
independently reviewed F1 `BorrowedUndirectedCsr::new(offsets, neighbors,
edge_ids, edge_endpoints)`.

### Validation

`tests/physics/percolation_network_input.rs`: parsed graphs are structurally
identical to builder-made lattices (V/E, degree profile, full-activity
component summary); every rejection variant carries the exact line number;
AllowParallel keeps 4 distinct edges where Reject errors; mixed policies are
expressible; parsed graphs flow through the N2 facade and through the F1
borrowed view with identical analysis.

### Ingestion benchmark (format decision evidence)

`cargo bench -p cmc-rs --bench percolation_network_input` — parse + CSR
build of a G(n, m) edge list, text generated outside the timed region:

| Case | V | E | Time | Throughput |
|---|---|---|---|---|
| er16k | 4 096 | 16 384 | 3.0 ms | ~5.4 M edges/s |
| er1m | 250 000 | 1 000 000 | 814 ms | ~1.23 M edges/s (13.1 MB text, ~16 MB/s) |

A million edges ingest in well under a second, so the roadmap's bar for a
custom binary or memory-mapped format is not met; text remains the only
network input format.


## D1 directed topology capability

Public directed capabilities (roadmap §11 D1), promoted from the F1
crate-private sketch:

- `DirectedGraphView: GraphView` — `arc_count`, `arc_id`, `arc_ids`, and
  `arc_endpoints(arc) = [source, target]`; dense `ArcId`s follow the F1
  ID contract.
- `DirectedOutGraphView` / `DirectedInGraphView` — separate adjacency
  capabilities; no symmetric fallback exists and `CsrLattice` never
  masquerades as directed.

Owned `DirectedCsr` (`topology/directed_csr.rs`): one dense physical arc
table plus two CSR incidence structures (outgoing by source, transposed by
target), built by `from_arcs` (input order = arc ID order). Public fields
with a documented invariant list and `validate()`, matching the
`CsrLattice` contract. Parallel arcs keep distinct IDs; a self-loop is one
out- and one in-incidence at its vertex.

Borrowed zero-copy views (`topology/borrowed_directed_csr.rs`):
`BorrowedDirectedOutCsr` (out CSR + arc endpoint table — the minimum
carrying stable arc identity) and `BorrowedDirectedCsr` (+ transpose
arrays). One-time validation on both sides: offset shape, incidence
lengths, neighbor/arc ranges, per-arc sense (`OutArcMismatch` /
`InArcMismatch`), and exact multiplicity — `ArcOutMultiplicity` /
`ArcInMultiplicity` catch both duplicated and missing arcs.

### Validation

`tests/physics/topology_directed.rs` (5 tests): owned exposes both
directions with correct arc identity on a ring + parallel-arc + self-loop
graph; empty graphs and isolated vertices are valid; owned/borrowed parity
for every vertex in both directions; every malformed representation is
rejected with the exact typed variant (12 rejection cases across offset
shape, lengths, ranges, sense, and multiplicity on both sides);
`validate` catches mutated public fields (sense check fires before
multiplicity for inconsistent corruptions; sense-consistent duplicates
trip multiplicity).

### Trait-overhead benchmark (F1-style local gate)

`cargo bench -p cmc-rs --bench topology_view` — directed additions on a
V=16384, A=48832 mixed-orientation grid, same shapes as the F1 gates:

| Scan | direct | generic owned | generic borrowed |
|---|---|---|---|
| out incidences | 79.8 µs | 68.9 µs (−14%) | 75.1 µs |
| in incidences | 81.9 µs | 68.5 µs (−16%) | ~75 µs |
| arc table | 18.6 µs | 21.8 µs (+17%) | 22.0 µs |

The adjacency scans — the hot path D2 analyses iterate — are within the F1
band (generic even faster, as in F1's incidence scan). The arc-table scan's
+17% is a storage-layout effect, not abstraction cost: the direct fold over
the compact `Vec<Arc>` auto-vectorizes, while the generic per-element cost
(0.45 ns) is *below* the recorded F1 generic edge scan (0.72 ns per
element, 23.472 us / 32512 edges).
Recorded as data; the D2 analyzers iterate incidences.

### Limitations

No directed analysis yet (D2: reachability, weak components, SCC). No
directed lattice builders — construct via `DirectedCsr::from_arcs` or the
transpose-validated borrowed views. No arc weights or arc activity types;
D2/D3 define those against their own occupancy contracts.


## D2 static directed analysis

Three analyses over the D1 directed capabilities
(`percolation/directed/`), each with its own workspace and result type:

| Question | API | Algorithm |
|---|---|---|
| Weak components (direction ignored) | `analyze_weak_components` | union-find over the physical arc table; canonical lowest-member labels |
| Strongly connected components | `analyze_strong_components` | iterative Tarjan (explicit DFS stack — deep graphs cannot overflow the call stack); adjacency materialized into the workspace once per analysis, then plain-slice DFS |
| source -> target reachability | `reachable` | early-exit BFS with stamped visitation (repeated queries never clear the visited array); `reachable(v, v)` is true by definition |

Summaries follow the platform conventions: canonical component identity is
the lowest member vertex ID; equal-size largest components resolve to the
lowest identity; `raw_second_moment = sum(s_i^2)`. `StrongComponentSummary
::acyclic` is false exactly when a directed cycle exists — a size-1 SCC
with a self-loop is tracked separately from multi-vertex SCCs.

### Validation anchors

`tests/physics/directed_analysis.rs` (3 tests, all green):

- **Exhaustive definitional anchor**: every digraph on 3 vertices
  (2^9 = 512) and on 4 vertices (2^16 = 65,536), self-loops included —
  SCC labels, weak labels,
  and every (source, target) reachability pair equal ground truth computed
  from the mathematical definitions by transitive closure
  (`u ~ v iff u reaches v and v reaches u`), independent of any SCC
  algorithm. This is a stronger anchor than any published benchmark: it is
  exhaustive over the input space.
- **Tarjan vs Kosaraju**: an independently written test-side Kosaraju
  (own DFS, own reverse sweep) produces identical partitions on the
  exhaustive corpus and on 40 random digraphs (V in 5..45, up to 200 arcs
  with self-loops and parallel arcs, xorshift RNG independent of the
  production RNG).
- **Weak vs flood-fill**: test-side undirected BFS matches on both corpora.
- **Transpose invariance**: reversing every arc leaves the SCC partition
  unchanged (all 40 random rounds).
- **Closed forms**: directed chain (1 weak component, 10 singleton SCCs,
  acyclic, reachability = order relation), out/in-stars (acyclic), ring
  (one size-6 SCC, cyclic), complete DAG (acyclic, full downstream
  reachability), single self-loop (singleton SCC but `acyclic = false`).

Literature: Tarjan's SCC algorithm (R. E. Tarjan, SIAM J. Comput. 1(2),
146-160, 1972) predates arXiv, so no PDF can be archived; the D2 gates
instead anchor to the *definition* exhaustively (all digraphs up to 4
vertices) plus algorithm-cross-validation (Tarjan vs Kosaraju), which is
strictly stronger than reproducing published numbers. Kosaraju's
formalization (Sharir 1981) is likewise pre-arXiv for this purpose.

### Steady-state allocation gate

`DIRECTED_ANALYSIS_ALLOC`: one weak + one strong + one far-corner
reachability analysis per round, after warmup, on all five benchmark
topologies: **0 allocations per round** (stamped visitation and prepared
workspaces).

### Criterion data

`cargo bench -p cmc-rs --bench directed_analysis` (V=16384):

| Case | A | weak | strong | reachability (far corner) |
|---|---|---|---|---|
| directed grid | 40 576 | 279 µs | 498 µs | 117 µs |
| random digraph | 49 152 | 719 µs | 1.20 ms | 422 µs |

### Limitations

No occupancy/activity contract yet — these are pure-topology analyses over
the full arc set. Directed *percolation* (occupancy, observables,
space-time layers) is D3, which builds on this stage and on C0 vector
output. Reachability is single-pair; multi-source or all-pairs variants
wait for a demonstrated need.


## D3 layered space-time directed percolation

`LayeredDirectedPercolationMC` (`percolation/directed/process.rs`):
directed **bond** percolation as an absorbing-state process. Layer `t+1`
site `j` is active iff some spatial neighbor `i` is active in layer `t`
and the directed bond `i -> j` is open at time `t` (independent
Bernoulli(p) per time step and directed incidence). The spatial lattice is
any owned `CsrLattice` (chain/square/...; periodic or open), so the same
adapter runs 1+1 and 2+1 dimensional processes.

State is a pair of frontier bitsets (`u64` words, double-buffered):
memory is `ceil(L/64)` words (= L/8 bytes) plus `O(T)` history counters,
never `O(L * T)`.
Initial conditions: `DpSeedPolicy::Single` (survival experiments) or
`::Full` (density-decay experiments). One `sweep()` draws one complete
realization; one `measure()` records it as ONE C0 vector sample:
`DpLayerActive[t]` (active count), `DpSurvival[t]` (0/1 indicator), and
the scalar `DpExtinctionTime` — the layer count is the array shape, so a
time series is never misreported as T independent scalar samples.

**RNG draw contract** (bitwise reproducibility; the parity gate replays
it): per layer, active sites ascending, CSR incidence order, one `f64`
draw per incidence, open iff `< p`; bonds out of inactive sites are never
drawn (they cannot change the process law).

### Validation anchors

`tests/physics/directed_process.rs` (5 tests):

- **Frontier DP = explicit DAG path enumeration** (the roadmap gate): an
  independent oracle replays the draw contract with its own bookkeeping;
  the production layer counts match the oracle sets exactly, and the fully
  materialized layered DAG (open bonds as arcs) checked with the D2
  `reachable` engine agrees node-by-node at every layer.
- **One realization = one vector sample**: after 16 sweeps both arrays
  carry exactly 16 samples per element (binsize 1).
- **Critical power law** — Hinrichsen, Adv. Phys. 49, 815 (2000)
  (arXiv:cond-mat/0001070, archived in `docs/ref/`, Tables 1-2):
  L=512 periodic ring, p = p_c(bond, 1+1D) = 0.6447001, single center
  seed, T = 1200: P(64)/P(1024) gives delta_hat within 0.05 of
  delta = 0.159464 at 2000 sweeps (env-tunable `SCUTTLE_DP_SWEEPS`).
- **Density decay**: full occupation at p_c, rho(64)/rho(1024) recovers
  the same delta (tolerance 0.06) — a second, independent anchor with the
  same exponent.
- **Phase brackets**: p = 0.60 extinguishes (P(1024) < 0.01); p = 0.68 on
  the same ring stays active (P(1024) > 0.5; the eventual plateau is the
  percolation strength ~ (p-p_c)^beta ~ 0.4 and t = 1024 spans only a few
  correlation times xi_parallel ~ 200, so the measurement sits in the
  transient approach ~0.67 — a statistical bracket, never an exact
  threshold claim).

### Steady-state allocation gate

`DIRECTED_PROCESS_ALLOC`: one complete 256-layer realization per round on
all five benchmark topologies: **0 allocations** (frontier pair + history
buffers preallocated).

### Criterion data

`cargo bench -p cmc-rs --bench directed_process` — 8 realizations per
iteration, single seed at p_c: chain L=512, T=1200: 1.53 ms per 8 realizations (~401 M
site-layers/s in Criterion's element convention; ~3.2 G physical
site-layer updates/s); chain L=1024, T=2000: 3.29 ms per 8 (~637 M/s).

### Limitations

Bond DP only (site DP and the Domany-Kinzel stochastic cellular automaton
are different processes); no activity-conditioned observables (k-moments,
pair correlations) or winding variants yet; spatial lattice must be an
owned `CsrLattice` (borrowed views via a future adapter if a need
appears); the extinction-time scalar is a per-realization integer
recorded as f64.

## E3 k-core percolation

`kcore_decompose` + `KCoreWorkspace` (`percolation/process/kcore.rs`) and
the `KCoreMC` adapter: the k-core is the unique maximal induced subgraph
in which every vertex `v` keeps at least `t_v` neighbors (`CoreThresholds`
uniform or per-vertex `u32`). Degrees use the platform multigraph
convention (a self-loop contributes 2, parallel edges count separately).

### API

- `kcore_decompose(graph, thresholds, initial_active, workspace)`: the
  reusable, RNG-free core analysis, generic over `UndirectedGraphView`,
  driven through the E1 primitives — `ProcessQueue` (two-phase FIFO),
  `NeighborCounters` (residual degrees), `StampFlags` (deleted/visited
  marks). Deterministic (ascending initial scan, CSR incidence order).
  Typed `KCoreError` rejections (threshold/activity length, `u32` counter
  range, capacity); zero allocations after `prepare`.
- `KCoreSummary`: core vertex count, core edge count (physical edges with
  both endpoints in the core), largest core component (one BFS pass),
  peeling rounds, deleted total.
- Post-run workspace access: `deletion_order()` (pop order),
  `round_survivors()` (length rounds + 1; index 0 = initial activity),
  `in_core(v)`, `deletion_round(v)` (synchronous wave), `residual_degree(v)`.
- `KCoreMC::new(lattice, keep, thresholds, recorded_rounds)`: one sweep
  draws the initial activity with the F2 uniform site law (bitwise
  identical to `StaticPercolationMC` site sampling — verified gate; exact
  keep = 0/1 shortcuts) and peels deterministically. One measure records
  ONE C0 vector sample `KCoreSurvivorCount[r]` (fixed shape
  `recorded_rounds + 1`; early termination pads with the final core size,
  overrun truncates the recording only — the D3 fixed-layer convention)
  plus scalars `KCoreFinalSize`, `KCoreRounds`, `KCoreDeletedTotal`.
  Researchers sweeping `k` or the keep probability run one adapter per
  value (each point is its own ensemble); documented, deliberately no
  built-in sweep.

### Semantics: synchronous rounds under FIFO peeling

Peeling pops each queue segment (vertices enqueued during the previous
round) as one round; a vertex crossing below its threshold during round
`r` deletions is deleted in round `r + 1` — exactly the simultaneous-wave
process, at one-incidence-touch cost. The final core is order-independent
mathematically (the min-residual-degree subgraphs are closed under union,
so the maximal element is unique); the shuffled-relabeling gate verifies
this empirically rather than the code enforcing it.

### Validation gates (roadmap §12 E3)

`tests/physics/kcore_percolation.rs` (11 tests):

- **Exhaustive initial subgraphs vs repeated full scans**: every `2^V`
  initial activity subset on path5/cycle6/star6/complete4/complete5,
  uniform `k in {0,1,2,3,6}` plus two heterogeneous threshold patterns —
  1456 (graph, threshold, subset) runs — match an independently written
  repeated-full-scan reference on core membership, every summary field,
  and the round count.
- **Endpoints**: `k = 0` removes nothing on 288 enumerated subsets;
  `k > max_degree` removes everything in exactly one synchronous round
  with ascending deletion order and survivor curve `[V, 0]`.
- **Queue-order independence**: 12 random relabelings (with shuffled edge
  insertion order) per (graph, k, subset) across path/cycle/star/complete
  — 90 distinct deletion orders observed, identical core sets and
  summaries throughout.
- **Hand-worked heterogeneous thresholds** with full wave accounting
  (path/star cascades traced by hand, survivors and `deletion_round`
  asserted).
- **Classic deterministic facts**: the 1-core of any connected graph is
  the whole graph with zero rounds; `k = 1` removes only isolates; the
  k-core of a tree is empty for every `k >= 2`.
- **Adapter contract**: keep = 1 is bitwise the direct decomposition;
  initial activity is bitwise the F2 site law on the same stream; one
  realization records exactly one vector sample (fixed shape, one bin per
  element per realization, non-increasing means); the recorded-round
  budget truncates/pads without changing the exact scalars; typed
  constructor rejections.
- **Literature anchor** — Janson & Luczak, Random Structures Algorithms
  30, 50-62 (2007) (arXiv:math/0508453, archived in `docs/ref/`),
  Theorem 1.1 (after Pittel, Spencer & Wormald): an own solver computes
  `lambda_k = min_{mu>0} mu/psi_{k-1}(mu)` (`psi_j(mu) = P(Poisson(mu) >=
  j)`) and the above-threshold fraction `psi_k(mu_k(lambda))`. Measured
  (G(n = 5e4, m), k = 3, 12 realizations per point, env-tunable via
  `SCUTTLE_KCORE_REALIZATIONS`):
  - solver `lambda_3 = 3.350919` at `mu* = 1.7933` (identity check, not
    an MC claim);
  - below, c = 3.2009: 3-core fraction 0.00000 ± 0.00000, largest core
    0.00000 — extinguished;
  - above, c = 3.5009: 3-core fraction 0.45309 ± 0.00134 vs theory
    `psi_3(mu_3(c)) = 0.45167` (agreement within ~1 sigma + finite-N
    slack), largest-core fraction 0.45309 (the giant carries the core).
    A statistical bracket with real error bars, never an exact finite-N
    threshold claim.

### Steady-state allocation gate

`KCORE_PROCESS_ALLOC` (allocation probe, k = 3, keep 0.5927, 256 recorded
rounds, 128 measured realizations after warmup): **0 allocations and 0
allocated bytes per realization** on all five benchmark topologies
(chain/square/cubic/sparse-er/power-law at V = 4096).

### Criterion data

`cargo bench -p cmc-rs --bench kcore_percolation` (8 realizations per
iteration, V = 4096, keep 0.5927, quick mode):

| Case | k = 1 | k = 3 |
|---|---|---|
| chain (E 4095) | 1.04 ms | 1.04 ms |
| square (E 8064) | 1.73 ms | 1.67 ms |
| cubic (E 11520) | 2.34 ms | 2.55 ms |
| sparse-er (E 16384) | 3.45 ms | 3.63 ms |
| power-law (E 16374) | 2.94 ms | 3.24 ms |

Per realization that is ~0.13-0.45 ms for V = 4096 at ~9-32 M peeled
vertices/s in the Criterion elements convention (throughput unit: one
vertex per realization).

### Limitations

Vertex process only (no edge activity/thresholds — a bond k-core variant
would be a separate process); the recorded-round budget is fixed at
construction (per-realization round counts vary, so the C0 vector needs
the documented truncate/pad convention); no core-number decomposition
(every vertex's maximal k, not just one k per run) — a future analytic
layer if a workload justifies it; the adapter owns a `CsrLattice`
(borrowed graphs use `kcore_decompose` directly); no built-in k/keep
sweep by design (one adapter per ensemble point).

## E2 bootstrap percolation

`BootstrapWorkspace` + `BootstrapPercolationMC`
(`percolation/process/bootstrap.rs`): monotone activation-threshold
closure. An inactive vertex activates once its **active-neighbor
incidences** reach a per-vertex threshold; a single O(V + E) closure
replaces repeated full scans. Buffers are the E1 primitives only — there
is deliberately no shared physical transition type (roadmap §12 E1/E2).

### API

- `BootstrapThreshold::Homogeneous { value }` or
  `::Heterogeneous { values: Vec<u32> }` (one per dense vertex).
- `BootstrapWorkspace::synchronous_rounds(graph, seeds, thresholds,
  trajectory)` / `::queue_closure(...)`: RNG-free core analysis, generic
  over `G: UndirectedGraphView` with seeds supplied as **any
  `VertexActivity`** (researchers bring arbitrary graphs and
  caller-sampled activity). `trajectory` is a caller-owned buffer holding
  the per-step active counts; `BootstrapOutcome` carries `seed_count`,
  `final_active`, `activated_total = final - seeds`, and the clock time
  (`rounds` or `events`).
- `BootstrapPercolationMC::new(graph, thresholds, seeds, dynamics, steps)`:
  the Carlo.rs adapter, generic over the owned graph type. Seeds are
  `BootstrapSeeds::Density { p }` (each vertex independently a seed; draw
  contract: one `f64` per vertex, ascending dense index, seed iff
  `draw < p`) or `::Fixed { vertices }` (explicit, deterministic
  realizations). One `sweep()` is one complete realization; one
  `measure()` records ONE C0 element-preserving array sample plus closure
  scalars:
  - `BootstrapActiveCount[i]`, `i` in `0..steps`: active count at recorded
    step `i` (step 0 = seeds). Entries past the closure repeat the final
    count; a closure longer than `steps` truncates the array while the
    scalars still report the full closure.
  - `BootstrapFinalActive`, `BootstrapActivatedTotal` (rule activations
    beyond the seeds), and `BootstrapRounds` (synchronous) or
    `BootstrapEvents` (queue) as the clock scalar.
- `is_active(index)` exposes the last run's closure set on the workspace.

### Conventions

- **Two explicit clocks** (the roadmap's "specify ... whether recorded
  time means synchronous rounds or queue closure order"):
  `BootstrapDynamics::SynchronousRounds` is the literature's parallel
  update (`A_0` = seeds; round `t + 1` activates every inactive vertex
  with `|N(v) ∩ A_t| >= r_v`; trajectory `[|A_0|, |A_1|, ...]`, `rounds` =
  update rounds until no change). `BootstrapDynamics::QueueClosure` is an
  asynchronous FIFO monotone closure (seeds ascending, then threshold-0
  vertices ascending, then first crossings; `events` = rule activations
  after the seeds). The **final closure is order-independent** in the
  queue clock (monotone-process fixed point); only the event sequence
  depends on the order.
- Threshold-0 vertices activate at the first update step in both clocks
  and count as rule activations, not seeds. Hence `r = 0` activates every
  vertex and `r > max_degree` preserves exactly the seeds (both gates
  asserted verbatim below).
- Parallel edges count once per incidence (the F3 physical-edge
  convention); a self-loop never helps its own vertex activate.
- Typed `BootstrapError` rejections: empty graph, zero recorded steps,
  heterogeneous-threshold length, fixed-seed domain, seed-activity domain,
  seed density, capacity.

### Validation gates (roadmap §12 E2)

`tests/physics/bootstrap_percolation.rs` (13 tests, ~15 s):

- **Exhaustive 2^V initial states vs repeated full scans**: path P5,
  cycle C5, star (4 leaves), K4, K5; every seed mask x every homogeneous
  r in `0..=max_degree + 1` (720 mask x r combinations) x both clocks —
  final closure SETS equal an independent full-scan reference, plus
  scalar invariants (`seed_count`, `activated_total = final - seeds`,
  `events` in queue mode). Per-round trajectories additionally equal an
  independent simultaneous-sweep reference (528 mask x r combinations).
  Heterogeneous thresholds: 15 seeded threshold vectors (values
  `0..=degree+1`) x all seed masks x both clocks (672 combinations).
- **`r = 0` activates all vertices**; rounds = 0 only when the seeds
  already cover the graph.
- **`r > max_degree` preserves only initial seeds** (final set == seed
  mask, 0 rounds, 0 rule activations).
- **Order independence**: on seeded random graphs (simple + parallel-edge
  variants) with random seeds and random heterogeneous thresholds, the
  FIFO production queue closure equals a test-RNG one-random-qualifying-
  vertex-at-a-time closure and a test-RNG shuffled-sweep closure, and
  equals the synchronous-rounds closure (24 draws x 7 graphs).
- **Hand-worked heterogeneous example**: path 0-1-2-3-4, thresholds
  `[2,1,1,2,3]`, seeds `{0,4}` — trajectory `[2,3,4,5]`, 3 rounds / 3
  events, full closure.
- **Adapter contracts**: fixed seeds reproduce the core bitwise
  (outcome + trajectory); the density draw contract replays bitwise;
  `measure` records the fixed-shape array + scalars (deterministic seeds
  give exact means and zero stderr); typed domain rejections.
- **Literature anchor — Holroyd, Probab. Theory Relat. Fields 125,
  195-224 (2003)** (arXiv:math/0206132, archived in `docs/ref/`):
  Theorem 1 proves `I(L,p) -> 1` iff `p ln L > pi^2/18` for 2D `r = 2`
  bootstrap, i.e. the finite-size seed threshold obeys
  `p_c(L) ~ pi^2/(18 ln L)`. All gates are statistical brackets with real
  error bars (400 realizations per point by default, env-tunable via
  `SCUTTLE_BOOTSTRAP_SEEDS`), never an exact finite-L threshold claim,
  because Holroyd's convergence is notoriously slow (finite-L thresholds
  sit far below the asymptotic constant — the very effect that made the
  pre-Holroyd numerics of Adler-Stauffer-Aharony misleading):
  - *His own observable* `I(L,p)` = P(full L x L torus activation) at
    formula-derived densities: at `p+ = pi^2/(18 ln L)` it is
    `1.0000 +/- 0.0013` for L = 32, 64, 128 (the `p ln L > lambda` side
    is already certain at these sizes); at `p+/4` (`p ln L = lambda/4`)
    it is `0.0000 +/- 0.0013` (L=32), `0.0025 +/- 0.0025` (L=64),
    `0.0000 +/- 0.0013` (L=128).
  - *Half-majority density* `p_hat(L)` (bisected): `0.0692` (L=32),
    `0.0577` (L=64), `0.0480` (L=128) — monotonically decreasing (the
    1/ln L decay), with `p_hat * ln L = 0.240 / 0.240 / 0.233` strictly
    inside `(lambda/4, lambda) = (0.137, 0.548)`: a fixed ~0.43 lambda
    fraction at accessible sizes, approaching the asymptotic constant
    from below, consistent with the theorem and never claimed to reach
    it.
  - *Sharpness*: the metastable window is narrower than a factor 1.5 in
    p at every size: `P(p_hat/1.5) = 0.0175 / 0.0050 / 0.0000` and
    `P(1.5 p_hat) = 0.9925 / 1.0000 / 1.0000` for L = 32/64/128.

### Steady-state allocation gate

`BOOTSTRAP_PROCESS_ALLOC` (r = 2, density 0.5927, 32 steps, 128 measured
realizations): **0.000 allocations and 0 bytes per realization** on all
five benchmark topologies (chain, square, cubic, sparse-er, power-law);
the probe asserts the zero.

### Criterion data

`cargo bench -p cmc-rs --bench bootstrap_percolation -- --quick`
(8 realizations per iteration; throughput counts `V + E` per realization,
so times below are per 8 realizations):

| Topology | V | E | SynchronousRounds | QueueClosure |
|---|---:|---:|---:|---:|
| chain (open) | 4096 | 4095 | 293.17 us (~224 M elem/s) | 599.15 us (~109 M elem/s) |

(The chain row above was recorded at seed density 0.5927 during
calibration; the committed bench uses 0.1, where both clocks sit within
~6% — review follow-up 2026-10-23.)
| square (open) | 4096 | 8064 | 1.4461 ms (~67 M elem/s) | 1.4008 ms (~69 M elem/s) |
| cubic (open) | 4096 | 11520 | 2.3039 ms (~54 M elem/s) | 2.3558 ms (~53 M elem/s) |
| sparse-er | 4096 | 16384 | 2.7526 ms (~60 M elem/s) | 2.9785 ms (~55 M elem/s) |
| power-law | 4096 | 16374 | 2.4326 ms (~67 M elem/s) | 2.6060 ms (~63 M elem/s) |
| torus 128x128 (p = 0.048) | 16384 | 32768 | 2.7712 ms (~142 M elem/s) | 2.7179 ms (~145 M elem/s) |

At the committed sparse density (0.1) the two clocks sit within a few
percent on every topology; the calibrated dense chain run (0.5927) showed
a twofold queue deficit — at full activation the FIFO event order walks
the path non-contiguously while round frontiers stay contiguous — which
is a workload note, not a clock deficiency.

### Limitations

Monotone activation only (no de-activation; k-core peeling is E3); the
array horizon `steps` truncates long trajectories (scalars stay exact);
no per-vertex activation time observable, spanning/winding diagnostics,
or activity-conditioned observables yet; `BootstrapSeeds::Density` is the
only built-in stochastic seed law (arbitrary laws go through the core
with a caller `VertexActivity`); the adapter owns its graph (borrowed
graphs use the core functions directly); time scalars are per-realization
integers recorded as `f64`.

## R1 external bond-configuration analysis

Status: **Integrator-verified** — delivered in-branch (integrator checks 2026-10-23;
review (roadmap §15 R1).

### Decoupling contract

R1 is the entry to the FK/random-cluster track: analyze a
**caller-supplied** bond configuration through the shared component
analyzer, with no coupling to the sampler that produced the bonds. A bond
configuration is one `bool` per physical edge of any
`UndirectedGraphView` (`true` = bond present in the subset). Independent
bond Bernoulli, a future random-cluster q = 1 reduction (R2), or the
actual FK bond draw of an Edwards-Sokal/Swendsen-Wang adapter (R3) are
all just masks at this boundary; nothing on the analysis side imports,
constrains, or identifies a sampler. The worm high-temperature even
subgraph is not an FK configuration and is not relabeled as one.

Semantics follow F3/W2 exactly: every vertex is active (bond percolation
- the configuration is an edge subset), an edge is active iff its mask
entry is open, an open self-loop is a counted union no-op, and parallel
edges are distinct physical edges with independent mask entries.

### API

- `ExternalBondConfiguration::new(graph, &mask)` validates the mask
  length against the graph's physical edge count (typed
  `BondConfigurationError::MaskLength`; an *empty* mask on a non-empty
  graph is this rejection, while an all-`false` mask of correct length is
  the legal "nothing open" configuration leaving every vertex an
  isolate). The mask is borrowed, not copied. The type implements
  `EdgeActivity` (mask) and `VertexActivity` (all-active), so it feeds
  `analyze`, `analyze_with_labels`, and `analyze_winding` directly as
  both activity sources; the analyzer re-validates the pinned dimensions
  against the graph it is given.
- `analyze_bonds(graph, mask, queries, workspace)` is the one-call
  convenience from raw mask to `AnalysisResult`; a wrong length surfaces
  as the analyzer's own `AnalysisError::EdgeActivityLength`.
- `BondSample` is an owned record of *which* edges are open — a sorted,
  deduplicated `Vec<EdgeId>` (chosen over a bitset: the producing
  sampler — a future R3 SW adapter — naturally emits open edges as it
  draws them; membership is `O(log k)`; sparse draws cost `O(k)`, not
  `O(E)`). `from_mask`/`to_mask`/`write_mask` (reusable buffer) convert
  between the dense-mask and open-edge views; a foreign dense index is
  the typed `EdgeIndex` rejection. It stores no probability, RNG state,
  or sampler coupling. The R2 q=1 reduction and the R3 Edwards-Sokal FK
  exposure are delivered and consume this contract (their sections below).

### Validation gates

`tests/physics/external_bonds.rs` (plus module unit tests):

- **Reference parity**: summary (component count, largest, second
  moment), canonical labels, and query outcomes match an independent
  adjacency-BFS flood fill on the open-edge subgraph across
  chain/square/random-multigraph topologies, random masks at p in
  {0.1, 0.5, 0.9}, three seeds, sixteen draws each, through both entry
  points (`ExternalBondConfiguration` composition and `analyze_bonds`).
- **Typed rejections**: wrong mask length from the constructor,
  `BondSample::from_mask`, and `analyze_bonds`; empty mask on a
  non-empty graph rejected; all-`false` full-length mask legal with all
  isolates; `write_mask` buffer-length rejection. Exhaustive 2^E parity
  on a parallel-edge + self-loop multigraph (self-loop on an open bond =
  union no-op, counted as active).
- **Winding bridge**: on the 2x2 torus, `analyze_winding` with an
  `ExternalBondConfiguration` equals a `StaticConfiguration` holding the
  same bonds (read back entry by entry from the law's draw) — summary,
  per-component (identity, size, rank), and generators — over 64 draws,
  with both wrapping and contractible cases observed.
- **Sampler decoupling**: a uniform Bernoulli producer and a geometric
  checkerboard producer (provably different masks) analyze identically
  through the same API, each at reference parity, and `BondSample`
  round-trips either producer's mask back to the same analysis.

### Benchmark

`cargo bench -p cmc-rs --bench external_bonds -- --quick`
(quick mode, V = 4096, p = 0.5, one query, reused workspaces; the two
paths analyze bit-identical masks — asserted in the bench — drawn through
the law's documented one-f64-per-edge contract):

| Topology | E | analyze_bonds (external mask) | analyze (StaticConfiguration) |
|---|---:|---:|---:|
| chain | 4095 | 59.1 us | 66.3 us |
| square | 8064 | 101.3 us | 106.9 us |
| cubic | 11520 | 140.1 us | 142.1 us |
| sparse-er | 16384 | 215.0 us | 223.4 us |
| power-law | 16374 | 157.0 us | 165.3 us |

The two paths run the identical F3 analyzer over identical activity;
across repeated quick-mode runs the difference stays within roughly ±10%
and flips direction (chain +12% for static in one run, -11% in another),
i.e. quick-mode noise and code-layout effects between the two
`EdgeActivity` monomorphizations rather than algorithmic work — the
external wrapper adds one indexed load per edge at most.

### Limitations

Analysis only: no sampler exists here, by design (R2 below validates the
random-cluster q = 1 limit; R3 below exposes an actual SW/Edwards-Sokal bond
draw). The mask is per *physical* edge, so site
activity is fixed all-active — mixing heterogeneous vertex activity with
an external bond mask goes through the F3 capabilities directly with a
separate vertex source. `BondSample` is a record, not a distribution:
it carries no probability or RNG state and cannot generate draws.

## R2 random-cluster q = 1

Status: **Integrator-verified** — delivered in-branch (integrator checks 2026-10-23;
review (roadmap §15 R2).

### The exact statement validated

**The q = 1 FK/random-cluster measure on a graph `G = (V, E)` equals the
product Bernoulli(p) measure on its bonds.** The random-cluster measure
with edge parameter `p` and cluster weight `q` is

`phi_{p,q}(A) = p^{|A|} (1-p)^{E-|A|} q^{k(A)} / Z(p, q)` over bond
subsets `A ⊆ E`, with `k(A)` the component count of `(V, A)` under the
bond convention (every vertex active). At `q = 1` the cluster factor is
exactly `1^{k(A)} = 1` — `k(A)` drops out of the weight entirely, on any
graph including parallel edges and self-loops — so
`Z(p, 1) = Σ_A p^{|A|} (1-p)^{E-|A|} = (p + (1-p))^E = 1` and
`phi_{p,1}(A) = Π_e [p if e ∈ A, (1-p) otherwise]`: the independent bond
Bernoulli product measure, configuration by configuration. This is exact
mathematics, not an approximation to be anchored numerically — the
founding q = 1 observation of Fortuin and Kasteleyn (Physica 57, 536
(1972); pre-arXiv, so no PDF is archived) — and the gates below prove it
by computation on this codebase.

### API

- `ClusterQ` is the explicit cluster-weight parameter, pinned to `q = 1`:
  `ClusterQ::q1()` constructs it, `ClusterQ::new(value)` accepts exactly
  `1.0` and returns the typed `RandomClusterError::UnsupportedQ` for
  anything else (the message cites the roadmap requirement that generic
  real-`q` samplers — Sweeny, Chayes-Machta — need a separate plan after
  R3). The parameter is carried explicitly so the domain widens in one
  place.
- `RandomClusterMC` is the owned scheduler adapter. It duplicates **no
  sampling law**: `sweep()` delegates to the F2 `BondBernoulli` (uniform
  `p` via `new`, or a validated per-edge field via `new_heterogeneous` —
  at q = 1 the heterogeneous FK measure is `Π_e Bernoulli(p_e)`) and
  mirrors the draw into a reused dense mask; `measure()` analyzes that
  mask through the R1 external-bond path (`analyze_bonds`) and records
  the F3 `ObservablePlan` scalars. The analysis plumbing is therefore a
  different path than `StaticPercolationMC`'s `StaticConfiguration`
  analysis — the two adapters agree bitwise at matched seeds (same
  `BondBernoulli` stream) and statistically at independent seeds.
- After `sweep()`, `mask()` (dense view), `configuration()` (the law's
  draw target), and `realization()` (R1's owned `BondSample`, the R3
  contract surface — an SW/Edwards-Sokal adapter exposes its actual FK
  draw as the same record) describe one and the same realization.
  `realization()` materializes on demand so the sweep/measure path
  allocates nothing; constructor rejections mirror the static adapter
  (topology, winding plans, moment exactness, query domain, law domain,
  invalid `p`).

### Validation gates

`tests/physics/random_cluster_q1.rs` (plus module unit tests):

- **Configuration-by-configuration (roadmap gate)**: exhaustive
  enumeration of ALL `2^E` bond configurations on a 5-edge chain-ish
  graph (square + pendant tail, so `k(A)` is not a function of `|A|`
  alone) and a 6-edge multigraph (self-loop, parallel pair, 4-cycle), at
  uniform p ∈ {0, 0.25, 0.5, 0.75, 1} and one heterogeneous per-edge
  field each. Per configuration: the FK weight `p^{|A|} (1-p)^{E-|A|} q^{k(A)}`
  with the q = 1 cluster factor computed (`1.0f64.powi(k(A))`) and
  multiplied in equals the Bernoulli product weight computed edge by edge
  (≤ 1e-12 from power-versus-product floating association); the
  enumerated FK probabilities sum to `Z(p,1) = 1`; `k(A)` and `S_max(A)`
  come from an independent flood fill cross-checked against the
  production analyzer on every configuration; and the observable means
  ⟨k⟩, ⟨S_max⟩ under the two full distributions coincide.
- **Bitwise sampler identity (delegation contract)**: same seed → the
  adapter's mask and its `BondSample::to_mask()` equal the direct
  `BondBernoulli` draw read back from a `StaticConfiguration`, for
  uniform p ∈ {0, 0.25, 0.5, 0.73, 1} and a heterogeneous field, several
  seeds.
- **Statistical (roadmap gate)**: 8x8 open square at matched p = 0.5 (the
  bond threshold), `SCUTTLE_ZSCORE_SEEDS` seeds (default 16), 10 000
  i.i.d. realizations per seed and path with independent RNG streams:
  two-sample z-scores of GiantFraction and finite-cluster susceptibility,
  `RandomClusterMC` (Bernoulli delegate → external-mask analysis) versus
  the `StaticPercolationMC` bond path (`StaticConfiguration` analysis) —
  per-seed |z| < 4 and |mean z| < 2.
- **Typed rejections**: `q ≠ 1` (including `NaN` and `1 + eps`) rejected
  with the roadmap-citing error; module unit tests cover invalid `p`,
  winding plans, query-domain and heterogeneous-field-length rejections,
  endpoint-exact measures (p = 0, 1), and the edgeless graph.

### Benchmark

`cargo bench -p cmc-rs --bench random_cluster_q1 -- --quick`
(realization = BondBernoulli draw + analysis; the q = 1 path adds the
dense-mask mirror that feeds the R1 external-bond analysis; plan
recording excluded from both paths equally):

| Topology | V | E | q1 realization (external path) | static bond realization |
|---|---:|---:|---:|---:|
| chain | 4096 | 4095 | 74.9 us | 78.3 us |
| square | 4096 | 8064 | 126.8 us | 131.5 us |
| cubic | 4096 | 11520 | 159.5 us | 168.2 us |
| sparse-er | 4096 | 16384 | 257.0 us | 260.9 us |
| power-law | 4096 | 16374 | 187.2 us | 190.7 us |

The two realization paths execute the same analyzer over bitwise-identical
draws, and the external path lands at or slightly below the static path
across repeated quick-mode runs — the same ±10% quick-mode noise band R1
recorded for `analyze_bonds` versus `analyze` (code-layout effects between
the two `EdgeActivity` monomorphizations), with the mask mirror inside
that band rather than a systematic premium. A separate quick-mode timing
of `BondSample` materialization (the on-demand `realization()` record,
`O(E)` plus one allocation) runs 19.0 us (chain), 49.3 us (square),
66.6 us (cubic), 105.7 us (sparse-er), 104.8 us (power-law) — roughly a
quarter to half of a realization, paid by R3 only when the draw record is
actually consumed.

### Limitations

q = 1 only, by design: `q ≠ 1` is a typed rejection, and generic real-`q`
samplers (Sweeny, Chayes-Machta) need a separate plan after R3. The
adapter owns its lattice and carries no periodic embedding — winding an
FK configuration goes through the core `analyze_winding` with an
`ExternalBondConfiguration` directly. No `FromParams` parameter schema
(construct the adapter directly). The sweep/measure path is
allocation-free; only `realization()` allocates, on demand.

## R3 Edwards-Sokal FK exposure

Status: **Validated** — independent review PASS 2026-10-23 (reviewer addendum in docs/reviews/2026-10-23-final-round-reviews.md).

### The actual draw, exposed and analyzed

`EdwardsSokalMC` (`percolation/edwards_sokal.rs`) is the sampler side of
the R track: an owned-`CsrLattice` MonteCarlo adapter over ferromagnetic
couplings (`J_e = Bond.weight`, validated strictly positive and finite at
construction — the typed `Coupling` rejection — with the inverse
temperature `beta >= 0` validated finite). One `sweep()` is one
Swendsen-Wang / Edwards-Sokal step:

1. **Conditional FK bond draw** given the current spins: every physical
   edge in ascending dense `EdgeId` order whose endpoints carry equal
   spins opens independently with `p_e = 1 - exp(-2 beta J_e)`
   (precomputed per edge); opposite-spin edges are closed with probability
   one. The draw is recorded AS IT HAPPENS into the R1 `BondSample` (open
   edges are emitted in ascending dense order, exactly the sample's sorted
   storage) — `last_bond_sample()` returns the actual draw, never a
   reconstruction.
2. **Cluster spin assignment**: an internal union-find over the open bonds
   builds the FK clusters, and each cluster receives one independent
   uniform spin.

The following `measure()` runs the shared F3 `analyze` on that same
`BondSample` (as `EdgeActivity`, paired with the all-active vertex source
`AllActiveVertices`) and records the F4 `ObservablePlan` scalars — the
composition mirrors `StaticPercolationMC`, so the exposed and analyzed
configurations are one object. Boundary queries become
`BoundaryCrossing*` means; on a proper chain these are `P_FK(i ~ j)`.
A frozen-spins constructor pins a caller-supplied `±1` configuration
forever, reducing every sweep to the pure conditional draw — the
deterministic ensemble the exact gates exhaust. Unlike the i.i.d. static
family this is a Markov chain: thermalize before measuring.

### RNG contract

Documented and bitwise-replayed by tests: the bond pass consumes one `f64`
per same-spin edge with `0 < p_e < 1` (exactly the F2 endpoint convention
at `p_e = 0` — `beta = 0` — and at `p_e = 1`, reachable only through
`exp` underflow) and nothing for opposite-spin edges; the spin pass
consumes one `f64` per FK cluster, at the cluster's lowest-index site
under an ascending site scan (a cluster-set invariant, so the count and
order are determined by the bond draw alone).

### Validation gates (roadmap §15 R3)

All in `tests/physics/edwards_sokal.rs` (4 tests, draw counts tunable via
`SCUTTLE_EDWARDS_DRAWS`, default 60 000) plus 6 module unit tests
(typed rejections, bitwise RNG-contract replays frozen and unfrozen,
endpoint determinism, observable recording, accessor state):

- **Exact bond-subset enumeration** (frozen spins, a 4-edge heterogeneous
  graph, two 4-spin configurations): all `2^E` subsets of same-spin edges
  enumerated; the sampler's empirical frequencies z-gated against the
  product Bernoulli law (`|z| < 5` per subset; opposite-spin edges closed
  in every draw), and for EVERY subset the adapter's measure-time path
  (`BondSample` as edge activity) and the R1 `analyze_bonds` raw-mask path
  produce identical summaries and query outcomes.
- **Partition function**: on the heterogeneous graph, the 4-ring, and a
  3-vertex parallel-edge + self-loop multigraph (self-loops same-spin by
  definition, open self-loop a union no-op), at `beta in {0.3, 0.9}`:
  the Edwards-Sokal double sum `sum_s sum_{A compatible} p^{|A|}
  (1-p)^{E-|A|}` (exhaustive over spins AND bond subsets), the q = 2
  random-cluster single sum `sum_A p^{|A|}(1-p)^{E-|A|} 2^{k(A)}` with a
  test-local union-find `k(A)`, and the independent brute-force Ising
  `sum_s exp(-beta H)` agree: double sum = single sum, and
  `Z_Ising = exp(beta sum_e J_e) * double sum`, each to `1e-12` relative.
- **`P_FK(i ~ j) = <s_i s_j>`**: a thermalized proper ES chain on the
  4-ring at `beta = 0.9`; the boundary-crossing mean of the exposed draw
  and the spin correlator read from the same chain's spins are each
  z-gated against the exact brute-force correlator (`|z| < 5`) and to each
  other (per-sample conditional-mean relation
  `E[s_i s_j | A] = [i ~_A j]`, so the empirical gap has exact stderr
  `sqrt((1-p)/n)`, gated at `|z| < 6`).
- **Worm-representation separation and consistency**: the worm's
  high-temperature representation (even subgraphs with `tanh(beta J)`
  weights; two-defect odd subgraphs for correlators) is checked to agree
  with brute-force spin enumeration on the partition function and both
  correlators — the level at which the two representations MUST agree.
  Configuration-level agreement is impossible by construction and never
  claimed: the worm's even subgraph has `tanh(beta J)` weights on a
  parity-constrained support and is not an FK configuration; nothing in
  the R track analyzes it as one (the assertion lives in the module and
  test docs).

Edwards-Sokal (Phys. Rev. D 38, 2009 (1988)) and Swendsen-Wang (Phys.
Rev. Lett. 58, 86 (1987)) pre-date arXiv — no PDFs are archivable (the
Tarjan-1972 precedent in D2); the exact double-sum identities above are
the anchors.

### Benchmark (`benches/edwards_sokal.rs`)

Quick-mode medians on the five shared topologies (V = 4096, `beta = 0.5`,
`p_e = 1 - exp(-1) ~ 0.6327`, full `ObservablePlan::all` plan with one
boundary query; one `es_sweep` = FK draw + cluster spins + `BondSample`
exposure, `es_sweep+measure` adds the F3 analysis and scalar recording):

| Topology | E | es_sweep | es_sweep+measure |
|---|---:|---:|---:|
| chain | 4095 | 79.4 us | 249.7 us |
| square | 8064 | 156.1 us | 460.8 us |
| cubic | 11520 | 192.7 us | 595.6 us |
| sparse ER | 16384 | 321.0 us | 1004.6 us |
| preferential attachment | 16374 | 232.6 us | 811.0 us |

The draw+exposure path sustains roughly 50-70 M edges/s; the full cycle
16-21 M edges/s (the F3 analysis dominates, as in the static family).
One `Vec<EdgeId>` collection per sweep rebuilds the `BondSample` through
R1's public constructor — the single documented steady-state allocation
of the cycle (every other buffer — masks are implicit in the sample,
union-find, cluster spins, workspace — is reused).

### Limitations

Full Swendsen-Wang dynamics only (no Wolff single-cluster variant; the
spec makes the conditional bond draw the requirement, and the alternating
cluster resampling is what makes the chain's bond marginal the FK
measure — an i.i.d. per-site spin redraw would NOT produce `2^{k(A)}`).
Winding observables on the FK draw are a typed `WindingNotSupported`
rejection (FK winding analysis awaits an embedding-aware follow-up; the
core `analyze_winding` + `ExternalBondConfiguration` path remains
available to callers). No `FromParams` schema (construct directly or via
`Run::from_parts`). The `q = 2` cluster factor enters through the spin
resampling, not an explicit weight; generic real-`q` sampling is provided
by Q1 (Sweeny) and integer-`q` Chayes-Machta by Q2.

## P1 particle candidate-pair adapter

Experimental — delivered in-branch on `feat/cmc-percolation`, awaiting
independent review. This is the geometry/candidate-pair infrastructure for
the P track (roadmap §14); it deliberately implements **no** overlap or
percolation observables — those are P2.

### Borrow, not duplicate

`percolation/continuum.rs` adds no second geometry implementation. It
borrows the particle module's `ParticleConfiguration` (transient copy from
flat `&[f64]` coordinates when callers have none), `OrthorhombicCell`
minimum-image displacements (through the `SimulationCell` trait methods on
the concrete type), and the `CellList` packed stencil search (`new`,
`particle_cell`, `fill_candidates`, `validate`). The geometry is accepted
concretely as `OrthorhombicCell<D>` rather than generically over
`SimulationCell<D>` because the reused cell-list bucketing is defined only
over orthorhombic side lengths — the wider trait would admit geometries the
stencil search cannot support. Other cell shapes require generalizing the
particle module's cell list first (out of P1/P2 scope).

### Ephemeral `PairGraph` versus a persistent CSR

The roadmap forbids materializing a persistent CSR for every particle
snapshot: an equilibrium simulation would rebuild it every sweep. The
delivered design is the transient counterpart — each analysis builds an
owned `PairGraph` (vertex count plus a flat `Vec` of pairs), feeds it to
the F3 analyzer as an `UndirectedGraphView`, and drops it. No adjacency
index is built: `incidences` are derived by an on-demand scan (the F3
analyzer never calls it — it scans the physical edge table), so the
snapshot costs `O(pairs)` memory for the duration of one analysis only.

### Pair order contract

Pairs `(i, j)`, `i < j`, are emitted cell-major: packed cells in ascending
linear index (row-major over the `CellList` grid), then within each cell
first endpoints in ascending particle index, each with accepted partners in
ascending index. The sequence is sorted by `(cell(i), i, j)` — not globally
by `i`. An exact lattice-vector translation (bit-identical wrapped
coordinates) reproduces the identical sequence; a generic translation
preserves the pair set but may reorder across cell boundaries. Every pair
is emitted exactly once, from the first endpoint's own cell; the `3^D`
stencil provably contains every partner within the cutoff because the cell
width is at least the cutoff. The boundary is inclusive: minimum-image
squared distance `<=` squared cutoff. Zero particles are a typed rejection
(`ContinuumError::EmptyCoordinates`) rather than an empty analysis.

### API

- `CandidatePairGraph::new(&[f64], &OrthorhombicCell<D>, cutoff)` —
  borrowed flat-coordinate adapter; `pairs()` performs one transient
  analysis (coordinate copy + cell-list build + stencil walk).
- `candidate_pairs_from_configuration(&ParticleConfiguration<D>, cutoff)` —
  zero-copy path for particle-module callers.
- `candidate_pairs_with_cell_list(&ParticleConfiguration<D>, &CellList<D>)`
  — reuses a caller-prepared cell list (validated with the particle
  module's own `CellList::validate`; stale lists are typed errors).
- Typed `ContinuumError`: zero dimension, non-finite/non-positive cutoff,
  empty coordinates, non-multiple flat length, non-finite coordinate,
  cutoff beyond the minimum-image maximum (half the shortest side),
  inconsistent cell list.

### Validation gates (roadmap §14 P2 candidate-pair prerequisite)

All in `tests/physics/continuum_pairs.rs` plus module unit tests:

- **All-pairs parity**: exact pair-set equality against an independently
  coded exhaustive `O(N^2)` scan with minimum image on random 2D and 3D
  periodic configurations, three seeds each (the search is independent;
  the metric formula is shared — the gate targets the stencil, not
  floating-point noise).
- **Inclusive boundary**: dyadic configurations with pairs exactly at the
  cutoff, directly and through periodic images, are kept; a cutoff a few
  ulp shorter drops them.
- **Periodic-image pairs**: pairs close only through wrapped images
  (including through both coordinates' images) are found.
- **Degenerate sizes**: one particle yields a valid one-vertex empty-edge
  graph; zero particles is the typed rejection.
- **Determinism**: identical inputs reproduce the identical sequence, and
  the flat and configuration entry points agree.
- **Translation invariance**: exact lattice shift reproduces the identical
  sequence (dyadic quantization makes the wrap bit-exact); a generic shift
  preserves the pair set (guarded to sit a safe margin off the cutoff).
- **F3 bridge**: `analyze`/`analyze_with_labels` on the `PairGraph` with a
  random site mask match an independent adjacency-BFS reference on every
  `ComponentSummary` field and on every canonical label (lowest active
  index per component).
- **Prepared-list path**: matches the fresh build; a list built for
  different coordinates is rejected.
- **Throughput record** (debug-build test record): N = 10 000 in 2D at
  unit density, cutoff 1.0 — 15 885 pairs vs Poisson expectation 15 708
  (ratio 1.011), 2.2-2.6e5 pairs/s. Release figures below.

### Criterion data (`benches/continuum_pairs.rs`)

Unit number density, cutoff set from target mean coordination
`z = rho * V_D * r_c^D`; `enumerate` = one full transient
`CandidatePairGraph::pairs()` per iteration (copy + cell-list build +
stencil walk — the honest no-persistent-CSR per-analysis cost), `analyze` =
enumeration plus the F3 component pass, all sites active. Full-run
Criterion medians (measured pairs per particle track `z/2` as expected):

| System | pairs | enumerate | analyze |
|---|---|---|---|
| 2D N=4096 z=2 | 4 112 | 2.66 ms | 2.85 ms |
| 2D N=4096 z=6 | 12 270 | 2.13 ms | 2.25 ms |
| 2D N=4096 z=12 | 24 510 | 2.32 ms | 2.53 ms |
| 2D N=16384 z=2 | 16 473 | 11.88 ms | 12.45 ms |
| 2D N=16384 z=6 | 49 507 | 9.38 ms | 9.67 ms |
| 2D N=16384 z=12 | 98 687 | 9.90 ms | 10.83 ms |
| 3D N=4096 z=2 | 4 130 | 8.22 ms | 8.30 ms |
| 3D N=4096 z=6 | 12 166 | 5.19 ms | 5.35 ms |
| 3D N=4096 z=12 | 24 589 | 5.73 ms | 5.91 ms |
| 3D N=16384 z=2 | 16 316 | 36.91 ms | 38.82 ms |
| 3D N=16384 z=6 | 49 256 | 21.83 ms | 22.18 ms |
| 3D N=16384 z=12 | 98 182 | 24.42 ms | 24.77 ms |

Peak throughput ~1.0e7 pairs/s (2D, N=16384, z=12). Note the honest
low-coordination cost inversion: at small cutoffs the cell grid
`floor(L/r_c)^D` over-subdivides (sub-occupancy cells), so per-cell
fill/sort overhead dominates over distance checks — z=2 is slower than
z=12 at equal N despite 6x fewer pairs. This is the reused `CellList`
layout; callers analyzing many snapshots of one equilibrium chain at fixed
cutoff should keep a `CellList` alive and use
`candidate_pairs_with_cell_list`, which skips the per-analysis build.

### Limitations

Orthorhombic cells only (see the geometry judgment above); no overlap,
radius, or percolation observable lives here (P2); the flat-coordinate
path pays one transient `O(N)` copy per analysis (the configuration and
prepared-list paths do not); winding/wrapping over candidate pairs awaits
the continuum geometry and configuration definitions required by the W
track contract.

## P2 Boolean percolation

Experimental — delivered in-branch on `feat/cmc-percolation`, awaiting
independent review. Initial scope (roadmap §14, verbatim): orthorhombic
cells; disks/spheres; monodisperse and finite-polydisperse radii; binomial
and Poisson point processes.

### Composition on P1

`percolation/boolean.rs` adds no second geometry implementation. Overlap
enumeration composes the P1 candidate-pair adapter at the scalar cutoff
`2 r_max` — the smallest cutoff that provably contains every overlapping
pair, since `r_i + r_j <= 2 r_max` — and then keeps exactly the pairs with
minimum-image `d^2 <= (r_i + r_j)^2` (inclusive, squared, no epsilon,
matching P1's boundary discipline). Monodisperse radii never pay the
filter: `r_i + r_j = 2r` IS the enumeration cutoff, and `(r+r)^2` is
bit-identical to P1's `(2r)^2`, so the P1 pair list passes through
untouched (unit-tested against the filtered all-equal-radii path).
Per-particle radii run one in-place `retain` recomputing each candidate
pair's squared distance with the particle module's own metric — one extra
metric call per candidate pair, no duplicated geometry. The P1
minimum-image ceiling carries over unchanged: `2 r_max <= L_min / 2`.

### API

- `Radii::monodisperse(r)` / `Radii::per_particle(...)` — typed validation
  (finite, strictly positive, non-empty; length checked against the
  particle count at every enumeration entry point, and re-validated so a
  direct enum literal cannot smuggle an invalid radius into geometry).
- `PointProcess::binomial(count)` / `PointProcess::poisson(mean)` — the
  SAMPLING side. Positions are uniform in the cell under a documented
  draw contract (particle-major, axis-ascending, one
  `random_range(0.0..L_axis)` per coordinate). The Poisson count draw is
  Knuth's multiplication method in chunks of at most 512 (Poisson
  additivity; keeps `exp(-part)` far from underflow for every accepted
  mean, cap 1e6); the binomial count consumes no randomness.
- `overlap_pairs(&[f64], &OrthorhombicCell<D>, &Radii)` and
  `overlap_pairs_from_configuration(&ParticleConfiguration<D>, &Radii)`
  — the MEASUREMENT side for existing coordinates (the equilibrium-chain
  case): the caller's accepted coordinates are analyzed, never resampled.
  Output is an ephemeral `OverlapGraph` (vertex count + flat pair `Vec`,
  `UndirectedGraphView` with derived incidences) following the documented
  P1 transient pattern.
- `BooleanPercolationMC::sampling(...)` — Carlo.rs adapter drawing one
  point-process realization per sweep and recording the F4
  `ObservablePlan` scalars (GiantFraction, susceptibility, ...) per
  measure. `thermalization_sweeps = 0` documented (i.i.d. samples). A
  Poisson realization may draw zero particles — a legitimate vacuum
  sample recorded as the all-zero summary (an empty-buffer REQUEST
  through `overlap_pairs` remains the P1 typed rejection: caller bug).
- `BooleanPercolationMC::measuring(cell, coordinates, radii, plan)` —
  frozen adapter owning an immutable coordinate copy: `sweep()` is a
  documented no-op and every `measure()` analyzes the identical
  configuration. Live equilibrium chains should call
  `overlap_pairs_from_configuration` per accepted snapshot instead.
- Typed `BooleanError`: P1 rejections wrapped, invalid radii, empty
  radii, radii-length mismatch, zero binomial count, invalid/oversized
  Poisson mean, Poisson-with-per-particle-radii (the radii pin the
  count), boundary queries with a random count, query-domain mismatch,
  winding observables (no continuum embedding — W-track contract),
  capacity.

Deliberate scope cuts, documented in the code: per-particle radii cannot
pair with a Poisson process (they pin the particle count; species-per-
particle DRAWN radii are a point-process extension beyond the roadmap's
initial scope); boundary queries are rejected for Poisson adapters (their
vertex domain varies per realization); no degree-array observable
(unbounded support would force an arbitrary bin cap — the moment scalars
carry the connectivity information); winding observables are rejected.

### Validation gates (`tests/physics/boolean_percolation.rs`)

- **All-pairs overlap parity**: exact pair-set equality against an
  independently coded exhaustive `O(N^2)` overlap scan with minimum image
  on random 2D and 3D periodic polydisperse configurations (species
  0.4/0.65/0.9 and 0.3/0.5/0.55), three seeds each.
- **Translation invariance**: a generic shift (then wrap) preserves the
  pair set and every label-free summary field (fixture guarded to sit a
  safe margin off all contact distances); an exact lattice shift of a
  dyadic-quantized configuration reproduces the identical canonical pair
  sequence.
- **Particle-permutation invariance**: a Fisher-Yates relabeling of
  positions+radii maps the pair set exactly (relabeled identities,
  identical connectivity) and leaves the summary's label-free fields
  unchanged.
- **Periodic-representative invariance**: shifting individual particles
  by random integer box vectors changes neither pair set nor summary.
- **Zero-radius limit (exact)**: at `r = 5e-3` in a 6x6 box (the smallest
  scale whose `(L/2r)^D` cell grid allocates comfortably), distinct
  particles (distances >= 1) yield zero pairs and exactly coincident
  particles yield exactly one pair — decided exactly in floating point.
- **Fully-overlapping limit**: at the largest permitted contact cutoff
  (exactly `L/2`) a 1D ring is the complete graph `K_6` — all 15 pairs,
  one component, second moment 36 — verified against the brute-force
  reference; in 2D (where the box diagonal exceeds `L/2`, so a complete
  graph is unreachable under minimum image) an 8x8 grid at `r = L/4`
  collapses into a single component (saturated connectivity).
- **Poisson count statistics**: the chunked-Knuth count draw (direct, n
  = 6000 draws at lambda = 37.7 and n = 800 at lambda = 1500 exercising
  the chunked path) satisfies `|mean - lambda| <= 4 sqrt(lambda/n)` and
  `|s^2 - lambda| <= 4 sqrt((lambda + 2 lambda^2)/n)`;
  the same gates hold through the adapter's `ActiveVertexCount`
  observable (3000 realizations). Measured: lambda = 37.7 -> mean
  37.818 +- 0.317, s^2 = 37.661 +- 2.771; lambda = 1500 (chunked) ->
  mean 1499.605 +- 5.477, s^2 = 1549.09 +- 300.05; adapter lambda = 40
  -> mean 39.946, s^2 = 38.41.
- **Equilibrium-measurement path**: a frozen 200-particle polydisperse
  configuration measured across 12 sweeps records zero-variance
  estimates for every always-defined observable (integer observables
  exactly zero; the rational GiantFraction at most one bin-mean-rounding
  ulp) with the coordinate buffer bit-identical throughout, and the
  zero-copy configuration entry point enumerates the identical overlap
  set as the flat path.
- **Determinism**: identical seeds reproduce bitwise-identical scheduler
  estimates (3D, 60 sweeps); different seeds differ.
- **Typed rejections** at integration level: empty coordinates (P1
  wrap), radii-length mismatch, cutoff beyond the minimum-image maximum;
  adapter-level rejections in the module unit tests.
- **Monodisperse 2D disk threshold bracket**: N = 50 000 binomial disks
  on the 400x400 torus, percolation probability = fraction of
  realizations with a component whose UNWRAPPED x-extent exceeds 80% of
  the side (spanning judged on minimum-image-unfolded coordinates — raw
  extents are meaningless on a torus, where a tiny seam-straddling
  cluster has raw extent ~L; the unfolding BFS also catches every
  x-wrapping cluster). Against the Mertens & Moore Table I anchor
  `eta_c = 1.12808737(6)` (PRE 86, 061109 (2012), archived in
  `docs/ref/`, BIBLIOGRAPHY.md row): the Wilson-95% interval at
  `eta_c/1.15` lies strictly below 1/2 and at `1.15 eta_c` strictly
  above 1/2, 12 realizations per side — the measured crossing of 1/2
  sits inside the bracket. Measured: eta = 0.98095 -> 0/12 crossings
  (Wilson-95% [0, 0.250]), giant fraction 0.030; eta = 1.29730 -> 12/12
  crossings ([0.750, 1]), giant fraction 0.950. Statistical bracket with
  real error bars; never labeled exact.

### Steady-state allocation (probe `report_boolean_case`)

The adapter allocates nothing per realization beyond the documented P1
transient: positions land in one retained buffer reserved for the count
ceiling, the F3 workspace and query-outcome scratch are reused, and the
observable plan records scalars under the standard first-measure/binsize
discipline. The transient (coordinate copy, cell-list buckets, pair
list) allocates EVERY realization by design under the no-persistent-CSR
rule, so warmup cannot exclude it; the probe instead replays the same
construction in isolation over identical coordinate snapshots (copied
into a preallocated arena inside the measured window) and asserts
allocation-count equality — every steady-state allocation is
attributable to the documented transient and the adapter adds exactly
zero. Both the sampling path (monodisperse and polydisperse, N = 2048,
128 measured realizations) and the frozen measuring path are asserted;
binomial counts keep the realized N fixed (a Poisson adapter differs
only by rare tail excursions growing the retained buffer once).
Measured (release, N = 2048, 128 samples, unit-density 2D box):
monodisperse (eta ~ 0.69) 3 691.7 transient allocations per realization
(523.9 KB), polydisperse 0.25/0.5/0.75 species (eta ~ 0.92; the larger
`2 r_max = 1.5` cutoff coarsens the cell grid) 1 808.6 per realization
(473.8 KB), frozen measuring mode 3 699.0 per measurement — and
`adapter_added_allocations = 0` asserted in all three records. The
transient is dominated by cell-list buckets, exactly the documented P1
pattern.

### Criterion data (`benches/boolean_percolation.rs`)

Unit number density, filling factor `eta` as the overlap knob;
`realize` = full sampling realization through the production adapter
(draw + overlap + F3 + recording), `measure` = frozen-configuration
analysis (the equilibrium path, comparable to P1's `analyze` figure).
`eta = 1.1281` is the disk-threshold operating point, not an exact claim.
Representative quick-run Criterion medians (sample-size 10, warm-up 0.5 s,
measurement 2 s; the bench records all 24 systems — measured pairs per
particle track `2^D eta / 2` as expected):

| System | pairs | measure | realize |
|---|---|---|---|
| 2D N=4096 eta=0.700 mono | 5 728 | 2.33 ms | 2.44 ms |
| 2D N=4096 eta=1.128 mono | 9 279 | 2.22 ms | 2.24 ms |
| 2D N=4096 eta=1.800 mono | 14 858 | 2.21 ms | 2.14 ms |
| 2D N=4096 eta=1.128 poly | 8 667 | 2.60 ms | 2.55 ms |
| 2D N=4096 eta=1.800 poly | 13 818 | 3.11 ms | 3.19 ms |
| 2D N=16384 eta=1.128 mono | 36 944 | 9.67 ms | 9.55 ms |
| 2D N=16384 eta=1.128 poly | 34 256 | 10.99 ms | 10.77 ms |
| 3D N=4096 eta=1.128 mono | 18 408 | 5.14 ms | 5.22 ms |
| 3D N=4096 eta=1.128 poly | 15 499 | 8.79 ms | 8.62 ms |
| 3D N=16384 eta=1.128 mono | 74 431 | 24.12 ms | 22.90 ms |
| 3D N=16384 eta=1.800 poly | 98 699 | 55.50 ms | 53.64 ms |

Peak ~4.1e6 overlapping pairs/s (3D N=16384 eta=1.128 mono). The
polydisperse per-candidate contact filter costs ~30-70% over the
monodisperse pass-through at equal eta — one metric call per candidate
pair on top of P1's enumeration-dominated cost, largest where the
candidate list has the most non-survivors (3D, mixed radii).

### Limitations

Orthorhombic cells only (the reused cell list); contact cutoff
`2 r_max <= L_min/2` (minimum image — the fully-overlapping 2D limit is
therefore connectivity saturation, not a complete graph); radii are
caller-assigned per particle (a finite-species DRAWN-radius point process
is an extension); no continuum winding/wrapping at the P2 adapter level
(the T4 displacement source below supplies the analysis layer);
known thresholds are statistical brackets only.

## T4 continuum winding embedding

Status: **Experimental** (roadmap v2 §2 T4; engineering gates green,
awaiting independent review). The v1 W-track note — continuum wrapping
"may reuse this analysis contract only after their own geometry and
physical configurations are defined" (v1 roadmap §16, quoted verbatim; the
v1 document stays untouched) — is discharged by this section: the geometry
is now defined and the W2 contract is reused unchanged.

### The definition (the deliverable)

`percolation/continuum_winding.rs` adds no analyzer. It defines how a
continuum pair edge acquires its integer `Z^D` cell vector and stores the
result in a [`PeriodicEmbedding`] table:

- **Derivation:** for a pair `(i, j)` at caller coordinates `x_i`, `x_j`,
  per axis `c` the raw in-buffer delta `r_c = x_j[c] - x_i[c]` and the
  particle module's minimum-image displacement
  `d_c = cell.displacement(x_i, x_j)[c]` differ by an integer number of box
  lengths. The stored component is `n_c = round((d_c - r_c) / L_c)`, the
  boundary crossings of the minimum-image connection (for in-box
  representatives; a general representative choice shifts it by a
  telescoping gauge amount, exactly W2's gauge invariance — re-wrapping any
  particle changes incident vectors but no cycle sum).
- **Tolerance, never silent rounding:** the dimensionless ratio must lie
  within `1e-9` of an integer. Consistent coordinates land within a few
  ulps; the tolerance only fires on inconsistent representatives (~1e8 box
  lengths apart, verified reachable through public construction), on
  non-finite ratios, or on magnitudes above `2e9` (the exact-negation bound
  for `i32`) — all one typed `ContinuumWindingError::NonIntegerCellOffset`
  naming the edge and axis. No saturation, no wrap.
- **Negation:** vectors are stored once in the pair's `i -> j` direction;
  the reverse query returns the exact negation through the lazy
  `DirectedDisplacement` view (structural antisymmetry, the W1 contract).
- **Precondition `r_max < L_c/2` on every axis** (strict; P1's cell list
  accepts the `L/2` tie, T4 deliberately does not), rejected at
  construction with the offending axis named. Why it is load-bearing: below
  half a box side, a pair within cutoff has a unique closest periodic
  image, so the minimum-image displacement is THE displacement of the one
  actual physical connection and `n` records exactly how many boundaries it
  crosses — the winding subgroup of these vectors is the winding of the
  cluster as a geometric object. At `r_max >= L_c/2` two images of the same
  pair can both lie within cutoff: the pair graph records one edge whose
  two geometric crossings are distinct and the minimum-image tie silently
  selects one, so integer winding stays a well-defined cochain but stops
  describing a definite physical connection (sub-ulp coordinate noise at
  the tie can flip verdicts).

### API

- `ContinuumEmbedding::from_pairs(pairs, coordinates, cell, cutoff)`
  accepts any P1 `PairGraph::pairs()` / P2 `OverlapGraph::pairs()` table
  plus the SAME flat coordinate buffer the pairs were enumerated from and
  the enumeration cutoff (the precondition needs it as data). Typed
  rejections before any state is kept: zero dimension, invalid cutoff,
  buffer length not a multiple of `D`, empty buffer, non-finite
  coordinates, the half-box precondition (axis named), out-of-range and
  self-loop endpoints, non-integer offsets (edge and axis named),
  capacity.
- `ContinuumEmbedding` implements `PeriodicEmbedding` (`edge_count` for
  pre-pairing), so the unchanged `analyze_winding(graph, embedding, ...)`
  runs on continuum configurations. Construction is a per-realization
  transient like the `PairGraph` itself (`16 + 4D` bytes per edge);
  steady-state analysis on a fixed graph and embedding allocates nothing.

### Validation gates

Unit tests in `percolation/continuum_winding.rs` (hand 1D/2D tables,
negation, precondition axis naming, helper tolerance cases, structural
rejections, gauge re-wrapping) plus `tests/physics/continuum_winding.rs`:

- **Exact negation + independent recomputation:** every edge of random 2D
  and 3D configurations (sizes 40/160/500, 3 seeds, cutoffs 0.3 and 0.95 of
  `L/2`) satisfies `d(i->j) = -d(j->i)` per axis AND equals an
  independently recomputed `round((d - r)/L)` with a division-form
  minimum image: 515,079 edges, 1,222,184 axis entries, 0 disagreements;
  non-endpoint queries return `None`.
- **Hand-derived tables:** the five-particle x-ring (four in-box edges at
  exactly the cutoff, one wrap edge at exactly the cutoff with stored
  vector `-1`) analyzes to one component, rank 1, generator `[1, 0]`; a
  two-particle wrap pair and a three-particle open chain containing a wrap
  edge are connected but rank 0 (a tree has no cycle; crossing is not
  wrapping).
- **P2 independent reference (the physics gate):** 80 Poisson Boolean-model
  configurations — 20 seeds x {spanning, sub-spanning} x {2D (L=40, mean
  N=2500, eta 1.5 / 0.7 against eta_c ~ 1.128), 3D (L=15, mean N=4000,
  eta 0.6 / 0.15 against eta_c ~ 0.34)} — 259,068 particles, 50,484
  components. Reference: a test-local BFS unfolding (float minimum-image
  accumulation along BFS tree edges, a division-form image, no production
  cell call) with the winding subgroup generated by the cycle sums of the
  non-tree edges; extent-based spanning (max-min unwrapped coordinate
  `>= L - tol`) recorded alongside. Every component: production rank ==
  reference rank and wraps == (reference rank > 0), zero disagreements.
  Winding and extent-spanning are logically independent in BOTH directions
  (review follow-up 2026-10-23: the boundary-exact hand ring winds with
  extent 8 < 10, and the review found 3/20 fresh small-box configurations
  with winding-but-not-spanning components — the previous claim that
  winding implies extent `>= L` held only by luck of these far-from-
  threshold densities): both events are census-counted, never asserted;
  at these densities spanning-without-wrapping and
  wrapping-without-spanning both read 0. Spanning regimes wrap in 20/20
  configurations each; sub-spanning regimes wrap 0/20.
- **Precondition boundary:** a cutoff of exactly `L_x/2` on a `[4, 10]`
  box is rejected naming axis 0 although P1's cell list accepts the tie;
  a cutoff just under `L_x/2` runs end to end (a three-particle boundary
  ring winds x exactly once, generator `[1, 0]`).
- **Tolerance rejections:** hand-built pair lists with representatives
  1e8 box lengths apart (residual ~1.5e-8) and with an integer-valued but
  unrepresentable offset (1e300 boxes) both fail with the typed error.

### Allocation probe (`report_continuum_winding_case`)

Fixed 2D Poisson configuration (N = 2048, E = 1588): construction
(pair enumeration + embedding) measures 5,735 allocations / 719,104 bytes
per realization — the documented P1-transient accounting, never asserted
zero; 128 steady-state samples of bond redraw + `analyze_winding` assert
0 allocations (workspace 196,640 bytes, 96.0 bytes/vertex).

### Limitations

Orthorhombic cells only (inherits P1); the strict `r_max < L/2` cutoff is
deliberately stronger than P1's `<= L_min/2` (documented above); no
scheduler adapter or plan observable wraps this yet (the analyzer is
direct-call, like F3's exact-analysis path); `analyze_winding`'s adapter
magnitude gate is not applied to continuum embeddings because each vector
is bounded by construction to `|n| <= 2e9` with small physical values.

Independent review passed 2026-10-23
(`docs/reviews/2026-10-23-t4-continuum-winding-review.md`): a third
arithmetic convention (own `rem_euclid` minimum image) over 1,156 fresh
randomized configurations / 1,831,913 edges with zero disagreements and
proven tie-unreachability at legal cutoffs; negation exact to `+-2e9`
with clean rejection beyond; gauge invariance over 90 re-wrap trials /
241,230 checks; 47 fresh-seed configurations / 10,586 components with
zero verdict disagreements against an independent BFS unfolding; the
`L/2` double-image ambiguity demonstrated (`0 -> -1 -> +1` under `1e-15`
noise). The review's one Major — the gate asserting `winding implies
extent >= L`, false in general (the hand ring itself: rank 1, extent
8.5; 3/20 fresh small-box configs) — was corrected as a follow-up: both
independence directions are now census-counted, the hand ring is pinned
as the compact-wrapping counterexample in a regression test, the census
total was corrected to 259,068, and the per-edge footprint to
`16 + 4D` bytes. Production verdicts were correct everywhere; the
finding was an overstrong gate claim, not an engine defect.

## T3 honeycomb and kagome periodic embeddings

Status: **Validated** (roadmap v2 §2 T3, the W1 builder-coverage tail;
independent review passed 2026-10-23 — see the end of this section and
`docs/reviews/2026-10-23-t3-embeddings-review.md`).

### The construction decision

Both plain builders historically assembled symmetric adjacency lists and
let `CsrLattice::from_adjacency` deduplicate directed entries into physical
bonds, dropping the per-bond construction direction — the W1 reason these
families had no embeddings. T3 restructures both builders to a direct
bond-emission pass shared by the plain and embedding variants
(`honeycomb_bonds`/`kagome_bonds` with a `track_displacements` flag, exactly
the `hypercubic_bonds`/`triangular_bonds` pattern): the plain builders now
consume the same emission, so plain and with-embedding construction cannot
drift by design. Each physical bond is emitted exactly once, from its
smaller-index endpoint, visiting pairs in ascending `(source, target)`
order — the order `from_adjacency`'s BTreeMap produced historically, so the
refactored plain builders stayed bit-identical: the in-repo test pins
plain-vs-embedding equality, and the historical pre-refactor parity was
verified out-of-repo by the independent review (which ran the verbatim
pre-T3 builders from `2a9165f` against production — bit-identical at all
25 swept sizes). That same order doubles as the displacement-table row order.

### Documented coordinate schemes (public contracts)

- Honeycomb (brick-wall): site `(x, y)` carries index `y*lx + x`; neighbours
  of every site in family order: right `(x+1, y)`, left `(x-1, y)` — both
  `HoneyX` — and the matched inter-row bond `(x+1, y+1)` for even `x`,
  `(x-1, y-1)` for odd `x` (`HoneyY`), coordinates mod `lx`/`ly`.
- Kagome: unit cell `(x, y)` hosts sublattices A/B/C at index
  `s + 3*(x + y*lx)`; neighbours `A: B C B(x-1,y) C(x,y-1)`,
  `B: A C A(x+1,y) C(x-1,y)`, `C: A B A(x,y+1) B(x+1,y)`; every bond
  `Kagome`.
- Displacement convention (identical to the hypercubic/triangular
  builders): a `+1` step along an axis from coordinate `c` contributes `+1`
  on that axis iff `c+1` equals the axis length (a fundamental-cell
  boundary crossing), a `-1` step contributes `-1` iff `c == 0`, else `0`;
  the vector is stored source-to-target and the reverse query negates
  exactly through the W1 `DirectedDisplacement`.

### Small-cell degeneracy (the point of the stage)

- Honeycomb `Lx = 2`: right and left of one site reach the same neighbour,
  so every row carries one doubled horizontal pair — two indistinguishable
  physical bonds carrying distinct displacements `(0,0)` and `(-1,0)`;
  within a pair the non-wrapping bond is emitted first (a documented
  emission-order contract fixing the table rows). `2x2` has 2 doubled pairs
  (6 edges on 4 distinct pairs); `2x3` has 3; `Lx >= 4` doubles nothing.
- Kagome: no parallel bonds exist at any `Lx, Ly >= 2` (each site's four
  neighbours are distinct). The W1 "collapsed covering bonds" are a
  reconstructability collapse, not a multiplicity one: at `Lx = 2` the
  left/right cross families join the same sublattice pairs with swapped
  cells — e.g. `{A(0,y), B(1,y)}` wraps (`(-1,0)`) while `{A(1,y), B(0,y)}`
  does not (`(0,0)`) — so endpoint data alone cannot attribute which pair
  wraps; the displacement is recorded at emission instead. `2x2` has 6
  wrapping bonds (2 A–B and 2 B–C left steps wrapping x, 2 A down steps
  wrapping y) out of 24; the wrapping-bond count is the closed form
  `2*ly + lx` at every size.

### Validation gates

Unit tests in `lattice::graph::tests`; winding integration in
`tests/physics/percolation_winding.rs`.

- **Plain-builder parity (gate 1):** with-embedding builders return
  bit-identical lattices to the plain builders — `edges` (endpoints, kinds,
  weights), `offsets`, `neighbors`, `edge_ids`, `n_sites`, `n_bonds` — over
  10 honeycomb sizes (2x2..8x4, even `Lx` only) and 10 kagome sizes
  (2x2..6x4, including 5x3 and 3x2).
- **Hand-derived 2x2 tables (gate 2):** every bond asserted in both
  directions — honeycomb 2x2: all 6 edges including both doubled pairs;
  kagome 2x2: all 24 edges including the 6 wrapping bonds. Spot-checked 3x2
  and 4x3 tables for kagome (all 7 wrapping bonds of 3x2, both directions)
  and the analogous even-`Lx` honeycomb 4x2 (full 12-bond table) and 4x3.
- **Independent geometry re-derivation (gate 3, property test):** for 8
  honeycomb and 9 kagome sizes, every bond's displacement is re-derived by
  an independently written checker that enumerates the documented
  neighbour-family steps of every site (per-directed-incidence bookkeeping,
  mixed-radix endpoint decomposition, wrapping cell arithmetic — no call
  into the emission logic): per unordered pair, the two directions must be
  exact negations, and the min-normalized displacement multiplicity must
  halve exactly to the emission's table. Every edge additionally checks
  exact reverse negation and the `is_zero` convention.
- **Table construction (gate 4):** both builders route their tables through
  `LatticeEmbedding::try_from_edge_displacements`, whose validation
  (length, self-loops, `i32::MIN`, capacity) runs on every construction in
  every test above.
- **Winding (gate 5):** fully occupied tori at sizes >= 2 (honeycomb
  2x2/2x3/4x2/4x3/6x4, kagome 2x2/2x3/3x2/4x3/5x4) report one component of
  winding rank 2 with determinant-±1 generators (the full-torus winding
  subgroup is all of `Z^2`), each case parity-checked component-by-component
  against the W2 lift-BFS reference; deactivating every nonzero-displacement
  bond (a cut of all wrapping bonds, including the wrapping half of each
  honeycomb doubled pair) leaves no wrapping component at any size, also
  reference-checked.
- **Bond content (gate 6):** closed forms asserted across sizes — honeycomb
  `3*lx*ly/2` edges with `ly` doubled pairs at `Lx = 2` (none at `Lx >= 4`)
  and `ly + lx/2` wrapping bonds; kagome `6*lx*ly` edges, zero parallel
  pairs, `2*ly + lx` wrapping bonds.

### Limitations

The adapter parameter schema and `StaticPercolationMC::with_embedding` are
unchanged: winding on these families is analyzer-level
(`analyze_winding`), and plan-level `Wrapping`/`WindingRank` still reject
the honeycomb/kagome builders loudly (the W2 behavior). Embeddings exist
only for the periodic builders (`Lx, Ly >= 2`, honeycomb `Lx` even), the
plain builders' domains. No new benchmarks: the emission replaces the
adjacency assembly with the same asymptotic work, and the W1 storage note
(`16 + 4D` bytes per edge, O(1) lookups) covers the tables.

Independent review passed 2026-10-23
(`docs/reviews/2026-10-23-t3-embeddings-review.md`, verdict PASS, zero
blockers/majors): the hand tables were re-derived from the rustdoc alone
BEFORE reading the tests (all rows matched); an independently written
generator (lifted coordinates, own division-form wrap arithmetic) matched
production exactly per bond, in emission order, at all 25 swept sizes;
the pre-T3 builders were run verbatim via a path dependency and produced
bit-identical lattices; the no-parallel-bonds claim was proven
algebraically (a repeated neighbour requires `lx|1` or `ly|1`) and the
covering-image reinterpretation of the W1 kagome note endorsed; winding
was re-verified with an own lift-BFS on ten fresh sizes plus adversarial
cut configurations (single-axis rank 1, wrapping-cut rank 0). Doc-only
follow-ups applied: the W1 bullet now states the honeycomb-doubles /
kagome-attribution-ambiguity distinction correctly, the parity wording
distinguishes the in-repo pin from the externally verified historical
parity, the subsumption of the roadmap's W1-style plaquette/cocycle gates
is stated, and one garbled kagome test comment was rewritten against the
asserted table.

## G1 incremental component process

Status: **Validated** (independent review passed 2026-10-23).

### Contract: infrastructure, not a law

Achlioptas (G2) and invasion (G3) share infrastructure, not a fabricated
common law trait — and the infrastructure they share is this:
**incremental union operations with online observables over any
`UndirectedGraphView`**. `IncrementalUnionWorkspace`
(`percolation/process/incremental.rs`) knows nothing about candidate edges,
frontiers, thresholds, or rules. The Z1 step (one uniformly random vertex
activation), the future G2 step (one accepted candidate edge), and the
future G3 step (one frontier growth) are all just sequences of
`activate_vertex`/`activate_edge` calls on this substrate.

Every accepted operation returns an online `IncrementalSummary` — component
count (active isolates included), largest size with its canonical identity
(the lowest dense vertex index in the component; equal-size ties resolved
by the lowest identity — the F3 convention), raw second moment
`sum(s_i^2)`, and active vertex/edge counts — and appends
`(step, summary)` to the workspace's `ProcessTrajectory`. **The step
coordinate is caller-defined and never interpreted**: an activation index,
an accepted-edge count, an invasion time, or anything else the driving law
means by "step" (roadmap §13: "record process coordinates without forcing
a specific law").

Semantics follow F3 exactly: an edge is active iff it was activated AND
both endpoints are active (the transition happens once, at the later of
the two events; a waiting edge goes live at the second endpoint's
activation); a self-loop counts once; parallel edges count separately.
Activation is order-independent — the summary after any sequence depends
only on the sets of activated vertices and edges.

### API

- `prepare(vertices, max_edges)`: reserve every buffer and reset all state;
  re-preparing with a different size is legal and leaves no stale state.
  The trajectory budget is exactly `vertices + max_edges` (each vertex and
  each physical edge is accepted at most once), so accepted operations can
  never overflow it.
- `activate_vertex(graph, vertex, step)` / `activate_edge(graph, edge,
  step)`: mutate, return the post-operation summary, record `(step,
  summary)`. Typed refusals (`VertexAlreadyActive`, `EdgeAlreadyActive`,
  index-out-of-domain, prepared-size vs graph mismatch) leave the state
  untouched and record nothing.
- `reset()`: deactivate everything and clear the trajectory for the next
  realization (O(V + E) byte clears; capacity retained). Per-element
  deactivation is out of scope — union-find cannot split a component.
- `summary()`: the live snapshot without recording (e.g. the state at step
  zero). `trajectory().steps()` / `.summaries()` expose the record as
  slices.
- `component_size(vertex)`: the scoring primitive for future rule-driven
  processes (G2 candidate selection) — read-only root walk, `None` for
  inactive vertices.

### Validation gates

`tests/physics/incremental_process.rs` (5 tests, plus 7 module unit tests):

- **Random-activation prefix parity** (the Z1 convention on the new
  substrate): random interleavings of ALL vertex and ALL edge activations
  on chain/square/random/multigraph topologies, four seeds each; every
  prefix summary equals a fresh F3 `analyze` pass over the same active
  sets (all six fields, including the canonical identity), and the
  trajectory mirrors the returned summaries with the caller-defined step
  values verbatim.
- **Online observable correctness, exhaustive orderings**: EVERY
  permutation of all 7 activation operations on the 3-vertex
  self-loop+parallel multigraph (5040 orders) and the 4-vertex star
  (5040), and all 6 on a path-with-self-loop (720) — each prefix against
  an independent from-scratch flood-fill reference (components, largest,
  identity, second moment, active edges). Subsets and order independence
  in one gate.
- **Canonical identity**: at every prefix, the G1 identity equals the F3
  `analyze_with_labels` selected largest component, the F3 label of that
  vertex is itself, the whole component carries the label, and
  `component_size` agrees.
- **Z1 compatibility demonstration (test-only)**: activating every edge
  first, then vertices in the `scan_site` activation order, reproduces the
  `ScanPoint` prefix summaries (components, largest, second moment, active
  edges) exactly on chain/square/random/multigraph; Z1 keeps its own fused
  scan in production.
- **Reuse**: reset cycles and grow/shrink/grow `prepare` cycles leave no
  stale state; eight random orders on the multigraph reach bit-identical
  terminal summaries; the trajectory budget fits `V + E` exactly.

### Allocation probe

`report_incremental_case` in `benches/percolation_allocations.rs`
(`cargo bench -p cmc-rs --bench percolation_allocations --features
allocation-probe -- INCREMENTAL`): one realization = reset + all E edge
activations + all V vertex activations in a fresh random order (V + E
accepted operations, each appending one trajectory entry). **Zero
allocations per realization are asserted** on all five benchmark
topologies (chain/square/cubic/sparse-er/power-law at V = 4096; measured
0.000 allocations and bytes, 128 realizations each).

### Benchmark

`cargo bench -p cmc-rs --bench incremental_process` (one iteration = 4
complete realizations on the shared topologies, V = 4096; throughput unit
= accepted operations, V + E per realization; workspace reused):

| Topology | V | E | Ops/realization | Time/iteration | Throughput |
|---|---|---|---|---|---|
| chain | 4096 | 4095 | 8191 | 1.18 ms | 27.83 M ops/s |
| square | 4096 | 8064 | 12160 | 1.91 ms | 25.44 M ops/s |
| cubic | 4096 | 11520 | 15616 | 2.41 ms | 25.97 M ops/s |
| sparse ER | 4096 | 16384 | 20480 | 3.16 ms | 25.94 M ops/s |
| preferential attachment | 4096 | 16374 | 20470 | 3.04 ms | 26.94 M ops/s |

### Limitations

Deliberately no process law lives here: no Achlioptas candidate rule, no
invasion frontier, no probabilities, no scheduler adapter — G2/G3 build
those on this substrate. Deletion/deactivation of individual elements is
not supported (union-find cannot split); the trajectory is a dense
`(step, summary)` record per accepted operation, not a per-vertex trace.
The largest-selection bookkeeping is exact but recomputed incrementally;
callers needing full component-size distributions at every step must
derive them from the summary stream or wait for a process that records
them.

## G2 Achlioptas processes

Status: **Validated** — independent review PASS WITH NOTES 2026-10-23 (addendum 3 in docs/reviews/2026-10-23-final-round-reviews.md).

### Contract: the process and its explicit parameters

An Achlioptas process turns random edge addition into a competition: every
step proposes `m` candidate edges, scores each against the current
component structure, and adds exactly ONE — the minimum-score winner
(giant suppression). `AchlioptasProcessMC`
(`percolation/process/achlioptas.rs`) drives the G1
`IncrementalUnionWorkspace` unchanged: every vertex is active from step
zero (an edge-addition process on a fixed vertex set; the initial state is
V isolates, largest component size 1), candidates are scored through
`component_size` (the substrate's designated scoring primitive), and each
accepted edge is one `activate_edge` call. The process runs until EVERY
permitted edge is added (terminal state = complete occupation of the
universe; roadmap gate).

Roadmap §13 demands that the physical parameters be explicit because "a
generic `score(size_a, size_b)` callback is insufficient to define the
process" — the process is defined by named, closed enums:

- **`AchlioptasRule`**: `Product` (`s_u * s_v`, the classic product rule)
  or `Sum` (`s_u + s_v`), evaluated as `u128` over the two endpoint
  component sizes. Selection is minimum-score in both.
- **`CandidateSampling`**: `WithoutReplacement` draws `m` distinct edges
  from the not-yet-added pool (partial Fisher-Yates; every step adds a
  fresh edge, exactly `|U|` steps total) or `WithReplacement` draws `m`
  iid uniform edges from the static universe — the two conventions the
  Riordan-Warnke paper itself describes ("from all possible edges [or]
  from those edges not present").
- **`EdgeUniverse`**: all physical edges, or a validated subset of dense
  edge indices (sorted, deduplicated; out-of-domain and duplicate entries
  are typed errors).
- **`TiePolicy`**: a tie is two candidate edges with EQUAL score in one
  step; `LowestEdgeIndex` (the default convention) selects the tied
  candidate with the lowest dense edge index, `EarliestDraw` the one
  drawn earliest. Repeated identical candidates (possible only under
  with-replacement) are one candidate.

**Null steps are decided, documented, and tested**: under
`WithReplacement` an already-added candidate competes normally (its score
`s*s` resp. `2*s` is legal), and when it wins the step consumes its draws,
advances the step counter, and adds NOTHING — the candidate pool is never
filtered or redrawn (filtering would silently change the law). The process
then terminates at coupon-collector time, a.s.

**RNG draw contract** (bitwise reproducible; verbatim in the module docs):
initialization consumes no randomness (vertices activate in ascending
dense order); each step consumes `k = m` draws (`WithReplacement`, each
`random_range(0..U)` through the static universe) resp.
`k = min(m, L)` draws (`WithoutReplacement`: the j-th draw is
`random_range(0..L-j)`, partial Fisher-Yates into the pool tail); after
selection the winner leaves the pool via a documented compaction (its tail
slot `L-1-r` swaps with the last slot `L-1`). `step_with_candidates`
drives one step from a caller-supplied candidate list with NO RNG —
deterministic replay of a recorded trace.

**Output (C0)**: one realization is ONE vector sample over the
occupation-count axis — `AchlioptasLargestSize[i]` (largest component size
at checkpoint `i`; 1 at zero occupation: the V active isolates) plus
`AchlioptasOccupiedEdges[i]` (fixed checkpoint coordinates, so merged rows
are self-describing) — plus scalars: `AchlioptasGiantReached` (0/1) and
`AchlioptasStepsToGiant` (the EXPLICIT percolation criterion: the first
accepted edge whose summary satisfies `2 * largest > V`, recorded as the
occupied-edge count), `AchlioptasSteps` (proposal steps incl. nulls),
`AchlioptasNullSteps`. Provenance strings (`as_provenance`) for run
metadata.

### Validation gates

`tests/physics/achlioptas_process.rs` (10 tests; plus 6 module unit tests):

- **Deterministic replay of a fixed candidate trace** (roadmap gate): a
  test-side driver reimplements the DOCUMENTED draw contract and drives the
  process through the public `step_with_candidates` API only — same seed
  reproduces the drawing sweep's accepted sequence BIT-TO-BIT (accepted
  edges, curve, giant step, terminal summary) on chain/square/random for
  m in {2, 4, 3} across both sampling modes; re-replaying the recorded
  trace on a fresh adapter (no RNG anywhere) reproduces it again. Frozen
  seeds give identical sweeps.
- **Every step vs naive BFS** (roadmap gate): after EVERY accepted edge of
  every realization, component count / largest size / raw second moment /
  active edge count recomputed from scratch by an independent flood fill
  and compared against the G1 trajectory summaries — chain/square/random/
  multigraph x {product m=2, sum m=3, product-with-replacement m=2} x 3
  seeds; the trajectory step coordinates (`0` x V, then 1..=|U|) and the
  `steps_to_giant` scalar are re-derived from the reference; the
  `AchlioptasLargestSize` curve mirrors the trajectory at every Full-grid
  checkpoint.
- **m = 1 reduces to ordinary random edge addition** (roadmap gate):
  bitwise — the accepted order equals a test-side uniform edge permutation
  implementing the documented m=1 contract under the same seed (the
  candidate IS the added edge; scoring cannot reroute it); statistically —
  mean largest sizes at four fixed occupation counts agree with the Z1
  `scan_bond` activation order within combined error bars (|z| <= 1.1, 40
  seeds, both sampling modes — with-replacement reduces in distribution
  because the first-appearance order of iid draws IS a uniform
  permutation); m=1 with-replacement total steps match the
  coupon-collector mean `U H_U` (measured 87.7 +/- 3.4 vs 90.6, z=-0.87,
  60 realizations) with exact bookkeeping `steps - null_steps = U`.
- **Terminal state = every permitted edge** (roadmap gate): exact —
  accepted set == universe set, all permitted edges live, pool empty, for
  AllEdges AND Subset universes under both sampling modes.
- **Product vs sum by hand**: components of sizes 5, 2, 3, 1 built by
  hand; candidates joining (5,1) and (2,3) score product 5 vs 6 but sum
  6 vs 5 — the two rules accept DIFFERENT edges from the same candidate
  set (hand-computed, unit + integration).
- **Replacement semantics**: occupied candidates can win (null step),
  leave the summary untouched, and the process still terminates complete.
- **C0 recording**: `Every`-grid vector shape, one vector sample per
  realization, constant coordinates, all four scalars, provenance strings.

### Literature anchor: Riordan & Warnke continuity

O. Riordan, L. Warnke, "Achlioptas Process Phase Transitions Are
Continuous", Science 335, 1189-1191 (2012) (arXiv:1102.5306, archived in
`docs/ref/`; the archived v3 is the extended Ann. Appl. Probab. version):
Theorems 1-2 prove that for any merging rule — in particular any
Achlioptas process, including the product rule — the rescaled giant
component evolves CONTINUOUSLY (no discontinuous jump, `L1(m) >= N_>=K(m) -
eps*n` whp at every step).

Anchor gate (product rule, m = 2, without replacement, fixed sparse ER
graphs of mean degree 4; every number a mean with real error bars over
independent seeded realizations; statistical claims, never exact):

- **Finite-cluster susceptibility** `chi_fin = (sum s_i^2 - S_max^2) /
  (N - S_max)` in the transition window 0.4-1.6 edges/vertex PEAKS at
  **363.4 +/- 16.3** (N = 2^12, 48 realizations) vs
  **1898.1 +/- 116.0** (N = 2^15, 24 realizations): growth x5.22 with a
  13.1-sigma separation — the signature of a CONTINUOUS transition (a
  discontinuous one would keep the susceptibility O(1)).
- **No O(1) jump**: the max single-window increase of the mean
  `S_max/N` curve SHRINKS as the window is refined (a discontinuous
  transition would saturate at a fixed jump): N = 2^12 — 0.224 over
  32-edge windows vs 0.626 over 128-edge (shrink ratio 0.36); N = 2^15 —
  0.401 over 256-edge vs 0.722 over 1024-edge (shrink ratio 0.56). All
  absolute jumps < 1. The 2^15 curve is genuinely steep — "explosive" at
  finite N, exactly why the Riordan-Warnke proof was needed; the window
  refinement and the susceptibility growth together are the
  accessible-size continuity evidence.

### Allocation probe

`report_achlioptas_case` in `benches/percolation_allocations.rs`
(`cargo bench -p cmc-rs --bench percolation_allocations --features
allocation-probe`): one realization = V vertex re-activations + E
candidate steps (2 draws, 2 `component_size` scorings, 1 union each) to
complete occupation. **Zero allocations per realization are asserted** on
all five benchmark topologies (chain/square/cubic/sparse-er/power-law at
V = 4096, 128 realizations each; measured 0.000 allocations and bytes).

### Benchmark

`cargo bench -p cmc-rs --bench achlioptas_process` (one iteration = 2
complete realizations, product rule m = 2 without replacement, `Every`
grid of 64 windows; throughput unit = process steps):

| Topology | V | E | Steps/realization | Throughput |
|---|---|---|---|---|
| chain | 4096 | 4095 | 4095 | ~9.3 M steps/s |
| square | 4096 | 8064 | 8064 | ~8.9-9.5 M steps/s |
| cubic | 4096 | 11520 | 11520 | ~9.9-10.4 M steps/s |
| sparse ER | 4096 | 16384 | 16384 | ~11.4-11.8 M steps/s |
| preferential attachment | 4096 | 16374 | 16374 | ~10-11 M steps/s |

(Measured on a busy shared machine; ~0.9-2.9 ms per iteration. Each step
= 2 candidate draws + 2 root-walk scorings + 1 union-find union + curve
bookkeeping.)

### Limitations

The percolation scalar is the explicit half-of-V criterion only; a
boundary-query arm (a component touching caller-defined vertex sets)
needs cross-component identity queries the G1 substrate does not expose
(`component_size` gives sizes; two distinct components may share a size),
and fabricating a second union-find to answer them would duplicate G1 —
deferred until the substrate grows a `same_component` query. The
without-replacement pool is an unordered set (draws are position-uniform),
so pool ARRANGEMENT after `step_with_candidates` may differ from a
drawing sweep at the same point (same content, same law, bitwise-identical
accepted sequences; only a sweep continuing to draw after a partial replay
would notice). Self-loops and parallel edges follow the platform
convention (a self-loop is a permitted edge that never merges anything).
The adapter owns its lattice (clone for repeated use); driving borrowed
views means calling the G1 substrate plus this module's step logic
manually.

## G3 invasion percolation

Status: **Validated** — independent review PASS 2026-10-23 (reviewer addendum in docs/reviews/2026-10-23-final-round-reviews.md).

One realization freezes a threshold field — one `f64` per vertex (site
mode) or per physical edge (bond mode) — and grows an invasion cluster
from explicit seeds through the **minimum-threshold frontier**: every step
accepts the frontier's smallest entry, its un-invaded endpoint joins the
cluster, and the entry's edge becomes the invasion-tree edge. The founding
reference is Wilkinson & Willemsen, J. Phys. A 16, 3365 (1983) — it
predates arXiv, so no PDF is archivable (the same Tarjan precedent as D2);
the delivered rule set is instead anchored to exact in-repo verification
and to the archived Sweeney & Middleton, Phys. Rev. E 88, 032129 (2013)
(`docs/ref/SweeneyMiddleton_2013_PRE_mst_percolation_threshold.pdf`),
which states the invasion-percolation/MST correspondence.

### API

- `invasion_run(graph, mode, thresholds, seeds, stop, trapping, workspace)`:
  the RNG-free core. Any finite thresholds are legal (the law does not
  require uniformity); NaN/infinity are typed rejections (`total_cmp`
  needs a total order). The workspace is re-prepared idempotently, so a
  correctly sized workspace never allocates in steady state. Returns an
  `InvasionOutcome` (seed count, accepted total, invaded vertices,
  trapped sites, max accepted threshold — the invasion barrier — stop
  reason and stop step); the workspace exposes `invaded_order()`,
  `parent()`, `accepted_edges()`, `accepted_thresholds()`, and the driven
  G1 `substrate()`.
- `InvasionMode::Site` / `InvasionMode::Bond`: which element carries the
  frozen field. Site mode queues each site once — when the first invader
  reaches it — carrying that invader's first incidence as the tree edge;
  bond mode queues boundary edges (exactly one invaded endpoint) and
  discards stale pops whose far endpoint was invaded meanwhile; parallel
  edges are independent frontier entries; self-loops never enter.
- `InvasionStop`: `GrowAll` (until frontier exhaustion — with non-trapping
  growth this invades the seeds' whole components), `FirstSpanning { from,
  to }` (checked after planting and after every accepted step, so two
  seeds already satisfying it stop at step zero), `Steps(n)` (exactly n
  accepted steps, zero legal). A frontier exhausted first stops with
  `FrontierExhausted` and the full accepted count as the stop step.
- `InvasionTrapping`: `NonTrapping` (default) or `Trapping { exit }` — the
  Wilkinson-Willemsen rule. A candidate site is invadable iff it is
  connected to an un-invaded exit vertex through un-invaded sites (an
  invaded outlet no longer provides escape; exit vertices are invadable
  through the normal minimum rule). The defending phase only shrinks, so a
  trapped verdict is final: sites retire permanently, each verdict costing
  one early-exit DFS of the defending region (O(V·(V+E)) worst case — the
  honest cost model, stated in the module docs rather than hidden).
- `InvasionPercolationMC`: the scheduler adapter over an owned
  `UndirectedGraphView`. One `sweep()` draws the field (bitwise contract:
  one `f64` per element in ascending dense order, uniform on [0,1)) and
  grows one realization; `measure()` records ONE C0 vector sample plus
  scalars. Configure `thermalization_sweeps = 0`; realizations are i.i.d.

### Explicit tie policy

Equal thresholds resolve by the **lowest dense ID** — the vertex index in
site mode, the physical-edge index in bond mode. The heap key is
`f64::total_cmp` on the threshold, then dense ID ascending; keys are
unique per entry (one heap entry per site / per edge), so the pop order is
deterministic even though `BinaryHeap` itself is unstable. `total_cmp`
orders `-0.0 < +0.0` as distinct keys (no tie) — documented and pinned by
a gate that would fail under an ID tie-break.

### Validation gates (roadmap §13, verbatim)

`tests/physics/invasion_percolation.rs` (16 tests):

- **Fixed-weight hand-worked graphs**: a 7-vertex, 9-edge multigraph
  (parallel pair 1-3, self-loop at 3) with hand-assigned thresholds. The
  complete sequence, the live-frontier evolution (a hand table of every
  frontier after planting and after each step — including the three-entry
  race for vertex 3 that the parallel edge wins, the never-queued
  self-loop, and the stale pops), and the parent tree are asserted against
  production in BOTH modes; the hand tables are simultaneously asserted
  against the scan reference, closing the loop hand → reference →
  production.
- **Binary heap versus scanning the complete frontier each step**: an
  independent test-side reference rebuilds the entire frontier from
  scratch before every acceptance (all edges in bond mode, all un-invaded
  frontier sites in site mode), sorts by the documented key, and walks it
  (retirement included). 576+ grid comparisons — chain / square /
  self-loop+parallel random multigraph / two-component graph × both modes
  × continuous and five-level-quantized (tie-saturated) threshold draws ×
  single and double seeds × GrowAll / Steps(3) / FirstSpanning — plus a
  trapping grid (both modes, exit at the far vertex). Every comparison
  checks order, accepted edges, accepted thresholds, parent array, stop
  reason/step, and every outcome scalar.
- **Tree uniqueness**: the invasion graph is a forest — exactly
  `invaded - seeds` accepted edges, no cycle (independent union-find
  replay: every accepted edge must merge two distinct classes), one
  parent per invaded non-seed vertex, one tree per seed (trees never
  merge, verified by BFS partition on two-seed same-component runs), and
  the forest spans exactly the invaded vertices.
- **Non-trapping bond invasion agrees with the applicable Prim/MST prefix
  relation when thresholds are unique** — stated precisely and verified:
  single seed ⇒ the accepted-edge sequence is exactly Prim's growth order
  from the seed (independent test-side Prim, all-edges cut scan) and the
  terminal tree is the unique MST of the seed's connected component
  (independent test-side Kruskal, edge-set equality); k seeds ⇒ Prim's
  order from a virtual super-root joined to every seed by `-inf` edges,
  terminal forest = the min-weight spanning forest with k trees, one seed
  each; `Steps(n)` ⇒ the exact n-edge prefix. Chain, square, and
  self-loop+parallel multigraph × four unique-threshold draws each
  (dyadic values, exact comparisons), all agreeing.
- **Explicit-policy gates**: tie hand examples (star graph, full and
  partial ties in both modes, the `-0.0 < +0.0` pin), stop policies
  (Steps 0/3/6/over-budget; FirstSpanning at the BFS-verified spanning
  step for every target — reachable at the stop step, unreachable one
  step earlier — plus the two-seed step-zero stop and the exhausted
  unreachable case), Wilkinson-Willemsen trapping hand traces (outlet
  invaded then branch sealed; near side trapped while the chain drains to
  a live outlet; bond mode where retiring a site consumes BOTH parallel
  frontier entries targeting it), the G1 substrate drive (trajectory
  steps are the invasion times, `[0, 1, 1, 2, 2, ...]`; online vertex and
  edge counts after every complete growth; two seeds in one component
  remain two substrate components — the substrate tracks the forest, not
  graph connectivity), the bitwise adapter draw contract (dense-order
  replay), C0 history compression (truncate / final-threshold padding /
  all-zero no-acceptance), typed rejections (length mismatch, NaN,
  infinity, empty/duplicate/out-of-domain seeds, spanning and exit
  domains, empty graph, zero-length recorded array), and workspace reuse
  (determinism, trapping-marks clearing, cross-topology re-prepare).

### Allocation probe

`report_invasion_case` in `benches/percolation_allocations.rs`: one
complete realization (frozen-field draw, heap-frontier growth to
exhaustion, substrate drive, C0 history compression, measurement) —
site, bond, and the trapping variant — **asserts zero allocations per
realization** on all five shared topologies at V = 4096 (128 measured
realizations each; measured 0.000 allocations and bytes, matching the
kcore/bootstrap/incremental/achlioptas probe policy). One subtlety fixed
during delivery: `InvasionWorkspace::prepare` now drains the heap BEFORE
`try_reserve_exact` (which budgets `len + additional`), so reuse demands
exactly `max(V, E)` capacity instead of transiently twice that.

### Criterion data

`cargo bench -p cmc-rs --bench invasion_percolation` (one iteration = 4
complete GrowAll realizations from seed 0 on the shared topologies,
V = 4096; every realization accepts the full seed component = 4095
growth steps; throughput unit = accepted growth steps; quick mode):

| Topology | V | E | Mode | Time/iteration | Throughput |
|---|---|---|---|---|---|
| chain | 4096 | 4095 | site | 1.010 ms | 16.2 M steps/s |
| chain | 4096 | 4095 | bond | 1.021 ms | 16.0 M steps/s |
| square | 4096 | 8064 | site | 3.279 ms | 5.00 M steps/s |
| square | 4096 | 8064 | bond | 5.311 ms | 3.08 M steps/s |
| cubic | 4096 | 11520 | site | 3.616 ms | 4.53 M steps/s |
| cubic | 4096 | 11520 | bond | 8.073 ms | 2.03 M steps/s |
| sparse-er | 4096 | 16384 | site | 4.817 ms | 3.40 M steps/s |
| sparse-er | 4096 | 16384 | bond | 12.484 ms | 1.31 M steps/s |
| power-law | 4096 | 16374 | site | 4.196 ms | 3.90 M steps/s |
| power-law | 4096 | 16374 | bond | 10.567 ms | 1.55 M steps/s |

Bond mode is 1.6-2.6x slower than site mode at equal E: every accepted
growth pushes all incidences of the new vertex as independent heap
entries (bounded by E total), versus site mode's one-entry-per-site
guarantee (bounded by V). Heap-versus-full-frontier-scan is a
correctness gate, not a benchmark axis (the scan is O(frontier) per step
by construction and is not production code).

### Limitations

`InvasionMaxAcceptedThreshold` is 0.0 for a zero-acceptance realization
(the adapter draws [0,1), so negative maxima arise only through the core
with caller-supplied negative thresholds — handled exactly via `reduce`).
The trapping escape check is an early-exit DFS with monotone retirement
(O(V·(V+E)) worst case; a split-maintaining union-find has the same
asymptotic worst case, so the simpler mechanism is the documented
choice). Site-mode tree edges are fixed at first push by the documented
carrying rule — a modeling choice, not a physical asymmetry. The adapter
owns its lattice (clone it per process); driving borrowed views means
calling `invasion_run` directly.

## F2 first-passage observables

Roadmap v2 §3 F2 (Experimental). `percolation/first_passage_observables.rs`
builds the FPP observables on the v2 F1 weighted shortest-path pass; the
anchor is Auffinger-Damron-Hanson, *50 years of first-passage percolation*,
AMS ULS 68 (2017) (archived, `docs/ref/`).

### Observable definitions

Two Z^2 source/target layouts, both on the open box (dense vertex
`x + y*width`, the `build_square` row-major convention) plus ONE
super-source vertex joined to every source-set vertex by structural
zero-weight edges — one F1 pass from the super-source realizes
`min over the source set` exactly (zero weights are legal Dijkstra
inputs):

- **Corridor** `PassageGeometry::corridor(length r, width w)`: source set
  the left column x = 0, target set the right column x = r. The profile
  has one entry per row: `T_y = min over paths from any left-column
  vertex to (r, y)`.
- **Segment** `PassageGeometry::segment(steps, margin, half_height)`: one
  origin vertex, explicit target points `s in steps` along the +x row of
  the origin, inside a box with symmetric longitudinal margin and
  transversal half-height. One profile entry per step — the
  point-to-point `T(0, x)` of the literature, several distances at once
  in one box.

One realization draws one weight per box edge (law below) and yields ONE
C0 vector sample: the passage-time profile along the target boundary. The
scalar passage time `T_r = min_y T_y` is a derived accessor, never the
primary sample. Path degeneracy never reaches the observable — profile
entries are unique minima over path weights — and tied minima across
targets keep all entries (no row is singled out).

**Weight laws and the bitwise draw contract.** `EdgeWeightLaw`:
Exponential(rate) `w = -ln(1-u)/rate` (rate 1 is the classical model) and
Uniform `w = u`, both transforms of one uniform draw `u in [0, 1)`. The
contract: one `rng.random::<f64>()` per box edge in dense `build_square`
order, then structural zeros consuming no draws. The exponential is
written `0.0 - ln(1-u)/rate` so `u = 0` yields `+0.0` (never `-0.0`,
which would pass F1 validation but bitwise-compare unequally against a
naive `-x` replay). `PassageSampler::weights()` exposes the last slice
(replay handle); `realize_with_weights` runs production on a caller-owned
field (the exact-enumeration path — weights are data, exactly as in F1).

**Estimator.** `PassageEstimator`/`PassageEstimate` (the `ScanEstimator`
pattern plus vector accumulation): scalar mean, unbiased variance, and
standard error of `T_r` across realizations, and per-target profile mean
and variance. The fluctuation observable is the scalar sample variance.

### Why subadditivity runs on segments, not corridors

The anchor inequality `E[T_{2r}] <= 2 E[T_r]` is exact in expectation for
the point-to-point passage time on the infinite lattice
(geodesic concatenation + translation invariance). On a reflection-symmetric
segment box the same argument is exact on a FINITE graph: the reflection
about the box's center line swaps origin and 2r-target, fixes the
r-target, and preserves the i.i.d. law, so `E T(0, re) = E T(re, 2re)` and
concatenation gives the inequality with no boundary correction. The
CORRIDOR observable instead carries entry/exit surface corrections (the
free min over w rows on each side), so its `E[T_r]/r` converges to the
cylinder constant from BELOW — the two families bracket complementary
physics and both are exposed. The roadmap's alternative sphere-hitting
formulation was not chosen: the sphere minimum does not satisfy the exact
finite-box inequality, and a high-power gate could then fail honestly.

### Validation gates (roadmap v2 §3 F2)

All in `tests/physics/first_passage_observables.rs`; every statistical
statement is a bracket, never an exact claim against 1/3.

**Exact enumeration.** Two-level dyadic weights {0.5, 1.0} (every path
sum exact in f64), ALL 2^E box-edge assignments, production profile vs a
test-local Bellman-Ford over test-built topologies (edge order matched to
production's dense order so one mask indexes the same edge in both):

- Corridors 2x4 (E=17), 3x3 (E=17), 3x2 (E=10): 263,168 assignments,
  919,552 profile entries, bit-exact, 0 mismatches. 128,272 assignments
  have a tied profile minimum (non-unique optimal paths) — e.g. the
  all-0.5 field makes every row attain `0.5 r` — and the profile content
  is asserted in every one of them (the tie policy: values are unique
  minima; no row wins).
- Segments ([2] m0 h1, [2,3] m1 h0, [1,2] m1 h1): 4,198,432 assignments,
  0 mismatches.
- Production-vs-test topology pinned by sorted endpoint-pair multiset in
  every geometry.

**Draw-contract replay.** 474 weight entries reproduced bitwise from
re-seeded streams (Exponential rates 1 and 2.5, Uniform; corridor and
segment), structural zeros exact.

**Determinism and error scaling.** Same seed -> bitwise-identical profiles
on two samplers (2 laws x 50 corridor runs, 2 x 20 segment runs); a
re-seeded stream reproduces chain one bitwise on a history-laden sampler
(no workspace-state leak). SEM ratio N=512 vs N=128: 0.531, inside the
1/sqrt(4) = 0.5 band (0.35, 0.65).

**Fluctuation bracket (chi).** Aspect-1 corridors (w = r, so the
transversal wandering `r^{2/3}` never meets the walls), r in
{8, 16, 32, 64, 128}, R = 600 realizations per size per law; OLS fit of
ln Var vs ln r with a t_{0.975}(df 3) CI on the slope/2:

| law | fitted chi | 95% CI | R^2 |
|---|---|---|---|
| Exponential(1) | 0.2476 | (0.2028, 0.2923) | 0.9904 |
| Uniform[0,1) | 0.2658 | (0.2198, 0.3118) | 0.9912 |

Both CIs sit strictly inside the rigorous open interval (0, 1/2)
(Kesten's bound; the survey's prediction chi = 1/3 in d = 2 is NOT
claimed). The effective exponent sits below 1/3 at these sizes — the
documented small-r crossover, not a disagreement; no larger-r rebracket
was needed (the CI does not clip an endpoint).

**Subadditivity bracket.** `E[T_{2r}] <= 2 E[T_r]` on reflection-symmetric
segment boxes (margin = r, half-height = r), Exponential(1), per-target
means from one box per pair:

| r | N | mu_r = E[T_r]/r | mu_2r | diff | z |
|---|---|---|---|---|---|
| 8 | 1200 | 0.555571 | 0.510286 | -0.045286 ± 0.004350 | -10.41 |
| 16 | 1200 | 0.510530 | 0.474303 | -0.036227 ± 0.002670 | -13.57 |
| 32 | 1200 | 0.473845 | 0.450111 | -0.023733 ± 0.001572 | -15.10 |
| 64 | 500 | 0.449508 | 0.433615 | -0.015894 ± 0.001535 | -10.35 |

No apparent violation anywhere (all z strongly negative; the consistency
threshold z < 3 never approached).

**Time-constant sequence.** One fixed segment box (steps 8..128, margin
64, half-height 64), N = 800: mu(r) = E[T_r]/r = 0.562055 ± 0.004709 ->
0.510794 ± 0.002721 -> 0.474296 ± 0.001633 -> 0.448849 ± 0.000982 ->
0.434147 ± 0.000602 — strictly decreasing with shrinking steps
(convergence from above, as subadditivity demands) toward the box-family
time constant (the literature's mu for Exp(1) weights on Z^2 is ~0.40;
the finite-box and finite-r corrections keep the sequence above it,
consistent, not claimed).

**Allocation probe** (`benches/percolation_allocations.rs`,
`F2_OBSERVABLES_ALLOC`, asserted): corridor 64x16 (V=1041, E=2015) and
128x32 (V=4129, E=8127) — 0 allocations and 0 bytes per realization
(law draw + F1 pass + profile copy + estimator add) over 128 measured
realizations; construction (geometry + sampler + estimator) measured
separately at 1077/4167 allocations (688,248 / 2,769,400 bytes) on the development machine — machine-dependent allocator counts; the independent review measured 1075/4165 (496,344 / 1,992,280 bytes), same order, informational only.

**Throughput** (eprintln-timed test, debug profile): aspect-1 corridors
64x64 (V=4161, E=8255) 4,283 us/realization = 233 realizations/s;
128x128 (V=16513, E=32895) 18,258 us/realization = 55 realizations/s.

### Limitations

Corridor fluctuation brackets use aspect-1 (w = r) geometry: at fixed
narrow width the geodesic is width-constrained and Var(T_r) crosses over
toward r^1 behavior, so width must scale with r (documented; the ladder
does). The segment boxes truncate the infinite lattice; the margins
(r longitudinal, r transversal) keep the truncation exponentially small
in `(margin / r^{2/3})^2` but never zero — the subadditivity inequality
itself is exact regardless (reflection), only the proximity to the
infinite-lattice mu carries the correction. Exponential draws use the
`-ln(1-u)` inverse transform (one transcendental per edge); a Ziggurat
would be faster but is not needed at these throughputs. No MonteCarlo
adapter and no scheduler integration — F3 (shape theorem, geodesic
aging) builds on this layer next.

## F3 shape theorem brackets and geodesics

Status: **Validated** — independent review PASSED 2026-10-24 (no Blocker;
the Major overflow-to-infinity tight-edge miscount is fixed with a
regression test, and the aging-convergence wording is corrected; the
remaining minors are recorded in `docs/reviews/2026-10-24-f3-review.md`)
 — all stage-F3 gates pass 2026-10-24; independent
review pending. The probe `tests/physics/first_passage_shape.rs` holds six
tests (uniform shape/counts, two-level counts vs brute force, directional
time constants, exploratory aging, draw-contract determinism, typed
rejections).

### Anchor statements (Auffinger-Damron-Hanson 2017)

Archived at
`docs/ref/md/Auffinger_2017_50yrs_first_passage_percolation.md`. The
following lines are quoted verbatim in the module rustdoc and drive the
gates:

- **Time constant** (Theorem 2.1, ~lines 392-405): under
  `E min[t_1, ..., t_{2d}] < infinity` the time constant `mu(e_1)` exists
  with `T(0, n e_1)/n -> mu(e_1)` a.s. and in `L1`. The subadditive proof
  gives `T(0, n e_1) <= T(0, m e_1) + T(m e_1, n e_1)` (~lines 442-464)
  and `lim E[T(0, n e_1)]/n = inf_n E[T(0, n e_1)]/n` (~lines 424-433):
  convergence from **above**.
- **Bounds** (~lines 518-526): `0 <= mu(e_1) <= E tau_e`, with Theorem
  2.4 (Hammersley-Welsh) making the upper inequality strict for any
  non-trivial passage-time law.
- **Shape theorem** (Cox-Durrett, Theorem 2.16, ~lines 745-755): a
  deterministic, convex, compact, axis-symmetric `B_nu` with non-empty
  interior satisfies `(1 - eps) B_nu subset B(t)/t subset (1 + eps) B_nu`
  eventually a.s.
- **Geodesics** (Definition 4.1, ~lines 3904-3914) and the uniquely
  geodesic spanning tree (~lines 4104-4117).
- **Overlap / persistence context**: Question 16 (~lines 2280-2283) asks
  for a positive power lower bound on the expected overlap of two
  independently sampled geodesics; Question 30 (~lines 5044-5046) asks
  whether the number of edges that belong to infinitely many geodesics
  `Gamma(v)` diverges on circles `|x| = r`. Neither gives a crisp
  finite-horizon aging law.

### Observable and reuse contract

`ShapeGeometry::new(directions, distances, margin)` builds one open
centered `side x side` Z^2 box (`side = 2*(max n + margin) + 1`) with the
origin at the center and one target per `(Direction, distance)` pair in
directions-major order; the four directions are `(1,0)`, `(0,1)`, `(1,1)`
and `(1,-1)`. All targets live in one box, so ONE F1 `first_passage_run`
per realization yields the full C0 passage-time profile. `ShapeSampler`
draws one uniform per box edge (F2 bitwise contract) under `EdgeWeightLaw`
and `realize_tracked` records the F1 predecessor tree; `ShapeEstimator`
keeps per-target mean, unbiased variance and SEM. `SHAPE_SIDE_CAP = 1024`
is a guard-rail, not a feasibility claim.

`geodesic_count_run` runs one F1 pass, sorts vertices by
`(distance, id)`, and propagates `ways[source] = 1` along every tight edge
`dist[u] + w == dist[v]` exactly. Counts are `u128` with explicit
saturation, parallel edges are distinct geodesics, self-loops cannot be
tight, unreachable targets count 0 and `source == target` counts 1. Zero
weights are rejected: a zero-cost tight cycle makes the geodesic count
infinite (F1 allows zero weights; this counting layer deliberately does
not). `geodesic_path` reconstructs the deterministic predecessor path
into caller buffers. `geodesic_aging` computes the shared leading-prefix
fraction of the geodesics to `n u` and `2 n u` in the same environment.

### Validation gates

**Exact shape and counts.** Uniform weights 1.0: all 36 vertex distances
of a 6x6 open box equal the l1 norm, and all 12 directional targets equal
`n*(|dx|+|dy|)`. Uniform corner-to-corner geodesic counts on 6x6 / 10x10
/ 12x12 are exactly C(10,5) = 252, C(18,9) = 48,620, C(22,11) = 705,432
(interior targets 126 / 24,310 / 352,716). Two-level dyadic {0.5, 1.0}
fields on 3x3 and 4x3 boxes, 8 seeds each (16 graphs): the tight-edge DAG
count equals a test-local brute-force enumeration of all simple paths
(pruned only by the target distance, sharing no code with the DAG
method), 7 graphs have non-unique geodesics, 34 geodesics counted in
total; all-0.5 fields give exactly the monotone-path counts C(4,2) = 6
and C(5,3) = 10. Every reconstructed path is simple and its weight sum is
bitwise equal to the F1 distance.

**Directional time constants.** Exponential rate 1, 97x97 centered box,
margin 16, n = 8/16/32, R = 150 (SEM in parentheses):

| direction | mu(8) | mu(16) | mu(32) |
|---|---|---|---|
| east (1,0) | 0.5653 (108) | 0.4978 (71) | 0.4750 (39) |
| north (0,1) | 0.5513 (107) | 0.5077 (65) | 0.4723 (39) |
| northeast (1,1) | 0.7637 (107) | 0.6995 (61) | 0.6551 (44) |
| southeast (1,-1) | 0.7715 (111) | 0.6997 (62) | 0.6539 (38) |

Every `mu` CI is positive and strictly below `E tau_e = 1`, every
consecutive change is negative within 3 combined SEM (subadditivity's
convergence from above), the 8->16 decrement exceeds the 16->32 decrement
within noise (finite-ladder flattening) and the latter stays at or below
0.06. Axis symmetry `|mu(east) - mu(north)|` is within 3 noise at n = 32.
Paired per-realization norm inequalities at n = 32:
`mu(east) + mu(north) - mu(northeast) = 0.29221 ± 0.00586` and
`mu(east) + mu(north) - mu(southeast) = 0.29345 ± 0.00599`, both
CI-positive. NO closed-form shape is claimed: the survey's shape theorem
is quoted as the reason to expect a norm, not as an object measured here.

**Exploratory aging.** Exponential rate 1, east, R = 60 per radius:
persistence 0.6551 ± 0.0427 (r = 8), 0.7294 ± 0.0405 (r = 16),
0.6368 ± 0.0445 (r = 24), 360 exact geodesic path checks. The drift
between the two largest radii is 0.0925 (noise 0.0602) inside the
deliberately generous 0.10 + 3*noise exploratory band. The survey has no
crisp finite-horizon aging law, so the statistic is labelled EXPLORATORY,
is not the Question 16 independent-overlap object, and no asymptotic
claim (monotone, convergent or otherwise) is made.

**Typed rejections.** Geometry (5), sampler (7), estimator (3), geodesic
counting (6) and path reconstruction (3) cases, with the previous
profile / counts / distances asserted untouched after every counting
rejection.

### Limitations

Finite boxes with free boundaries and a fixed margin measure
`T_box(nu)/n >= mu_infinite(u)` only as a statistical bracket; the
flattening band is a three-point ladder statement, not a limit. Geodesic
counts use exact f64 equality on the tight-edge relation: dyadic fields
make that exact, and arbitrary continuous weights make non-unique
geodesics vanish, so the count is then 1 per target (the interesting
counting regime is multi-level weights). The DAG sort and way buffer are
reused but no allocation probe was added at this stage. The aging
statistic is one tie-broken path per target in one environment; a different
tie policy or an all-geodesic overlap would define a different statistic.
No MonteCarlo adapter exists at this layer.

## L1 long-range power-law edge law

Status: **Validated** — all stage-L1 gates pass 2026-10-23; independent
review passed 2026-10-24 (the review's `coordinates()` d >= 2 Major and
the O(log V) doc Minor are fixed as follow-ups; see
`docs/reviews/2026-10-24-l1-review.md`).

### Contract: Biskup's edge law on a finite box

The anchor (archived
`docs/ref/md/Biskup_2003_long_range_chemical_distance.md`) defines
long-range percolation through independent bond variables with
`p_xy = |x-y|^{-s+o(1)}`, `|.|` any norm equivalent to ℓ2 (~line 134:
"for definiteness, we can take |.| to be the usual ℓ2-norm; however, any
other equivalent norm will do"), and works with the explicit
representative `p_xy = 1 - exp{-β(1+|x-y|)^-s}` (Eq. 1.7). Stage L1
realizes the explicit power-law member of that family,
`p(r) = min(1, beta * r^{-s})` (first-order equivalent to Eq. 1.7 at
small `beta r^{-s}`), on a finite box of Z^d; the Biskup regime structure
(s versus d and 2d) is the L2/L3 target, not claimed here.

Every semantic choice is explicit:

- **Norm: ℓ1** (`LongRangeBox::distance`, a `usize`). Integer-exact (the
  occupation probability depends on floating rounding only through one
  `powf`), equivalent to Biskup's ℓ2 so the regime statements transfer,
  and identical to ℓ2 in d = 1, the stage-L1 target dimension.
- **Box and boundary**: `LongRangeBox::new(&sides)` builds a product of
  per-dimension side lengths (each >= 1; vertex/pair counts must fit a
  usize — typed rejections). The pair universe is every unordered pair of
  DISTINCT box vertices — the finite-box restriction of the
  infinite-measure law: no edge leaves the box, each pair keeps the same
  `p(r)`, and boundary effects are confined to the box surface (a vertex
  at depth k misses exactly the partners the infinite lattice would offer
  beyond it).
- **Consideration order**: lexicographic by dense (x, y), x < y —
  ascending pair index (`pair_index`/`pair_endpoints` are inverses; dense
  order is coordinate-lexicographic order).
- **Draw contract** (bitwise reproducible on one platform): one uniform
  `f64` per pair in ascending pair order — saturated and near-zero pairs
  included, nothing skipped. The pair at index k is occupied iff
  `uniforms[k] < pair_probability(r, beta, s) = min(1, beta * r^(-s))`
  (the pinned expression; same-platform replay mirrors its bits).
- **Saturation**: `p = 1` exactly when `beta * r^{-s} >= 1` — in
  particular every r = 1 pair is surely occupied at `beta >= 1` (`u < 1`
  holds for every uniform in `[0,1)`), and a huge beta saturates every
  pair within radius `beta^{1/s}` (the complete graph in the limit).
  With `beta` finite, `s >= 0`, `r >= 1` the factor `r^{-s}` lies in
  `(0, 1]`, so the product never overflows past beta and no (s, d)
  combination is rejected. `beta = 0` (empty graph) and `s = 0`
  (distance-independent `min(1, beta)`) are legal but degenerate.
- **Process coordinate**: one realization activates all `n` box vertices
  first (dense ascending, step coordinate 0 — the G2 convention; the
  fixed vertex set is the stage) and then every occupied pair in
  consideration order with step coordinate the 1-based ADDED-EDGE COUNT.
  The realized count is random; the trajectory steps are exactly
  `[0; n] ++ [1..=m]`, one `IncrementalSummary` per accepted operation.

### API

- `long_range_run(box, beta, s, uniforms, workspace)`: the RNG-free core.
  `uniforms` is one `f64` per pair in consideration order, each finite
  and in `[0, 1]` (typed rejections otherwise; rejections leave the
  workspace untouched). The workspace is re-prepared idempotently.
  Returns a `LongRangeOutcome` (added edges, final component count,
  largest size, raw second moment); the workspace exposes
  `occupied_pairs()` (consideration order), `component_size(vertex)`, and
  the driven G1 `substrate()`.
- `LongRangeBox`: `new(&sides)`, `sides()`, `dimension()`, `vertices()`,
  `pairs()`, `vertex_index(coords)`, `coordinates(index)`, `distance`,
  `pair_index` / `pair_endpoints`, `pair_probability(distance, beta, s)`
  (the contract expression).
- `LongRangeWorkspace`: `prepare(vertices, pairs)` / `reset()` /
  accessors; `LongRangePercolationMC::new(box, beta, s, steps)` with
  `realize(rng)` and the Carlo.rs adapter (`sweep` + `measure`).

### C0 recording

One realization is ONE element-preserving vector sample plus scalars (the
G3 padding convention for random-length growth): `LongRangeLargestSize`
has fixed length `steps`, entry i = largest component size after the
(i+1)-th added edge, entries past the realized length repeat the final
value (1 when nothing was added); the scalars `LongRangeAddedEdges`,
`LongRangeComponentCount`, `LongRangeSecondMoment` carry the true values.
Configure `thermalization_sweeps = 0`; realizations are i.i.d.

### Validation gates

`tests/physics/long_range.rs` (8 tests):

- **Threshold replay (law + pair order)**: for each of the 28 pairs of a
  1x8 box, three (beta, s) laws: a uniform bit-adjacent BELOW the
  recomputed `min(1, beta r^{-s})` occupies exactly that pair (all other
  uniforms 1.0, which never occupy — the strict `<` at the saturation
  boundary), one bit-adjacent ABOVE does not; the pair-order pin is that
  the tailored uniform sits at index k. `p` also cross-checked against
  the direct `beta / r^s` reading (relative < 1e-12).
- **Law frequencies + independence**: per-pair occupancy over 600
  realizations on 1x8 at (0.8, 0.5), (0.3, 1.0), (2.5, 1.7) — the last
  with 7 saturated pairs (asserted ALWAYS occupied) — and on 4x4 at
  (0.5, 1.2) exercising the 2D ℓ1 distance. Max per-pair |z| = 2.50 /
  3.54 / 2.59 / 2.25 (all < 5); pooled Pearson chi-square 33.3 / 36.0 /
  23.2 / 110.2 against df = 28 / 28 / 21 / 120, Wilson-Hilferty critical
  (alpha = 0.001) 57.0 / 57.0 / 46.9 / 173.7. Joint-occupancy
  independence probes (3 pair-of-pairs per box) z in [-2.26, 1.96]
  against the exact product law.
- **Degree moments vs the LRP measure**: closed-form
  `E[deg(x)] = sum_{y != x} min(1, beta |x-y|^{-s})` (summed explicitly
  and cross-checked against the d=1 hand shell formula) vs the empirical
  mean center degree on 1x16, R = 4000, s < 1 / s = 1 / s > 2 laws:
  E = 6.7115 / 1.5932 / 0.6415, z = -0.41 / -0.59 / 0.66 (all < 4).
- **Bitwise draw contract**: same seed twice reproduces history,
  occupied pairs, trajectory steps and summaries bitwise; the adapter
  replay equals a manual one-f64-per-pair draw through the RNG-free core
  (the documented order pinned); different seeds differ; the C0 padding
  convention verified entry-by-entry against the substrate trajectory.
- **Trajectory parity**: 6 (box, law) configs x 5 seeds — 773 added-edge
  steps checked against a test-local independent union-find (all six
  summary fields, per step), the vertex stage checked entry-by-entry,
  steps == `[0; n] ++ [1..=m]`, and 30 final structures checked against a
  test-local flood fill (component count, largest size and canonical
  identity, second moment, `component_size` for every vertex, occupied
  set equality).
- **Saturation**: r = 1 pairs ALWAYS occupied at beta >= 1 (three laws x
  40 seeds); the 38 pairs with r <= 4 at (beta=4, s=1) always occupied
  (40 seeds); beta = 1e6, s = 2 gives the complete 1x8 graph every time
  (28 edges, 1 component, second moment 64); beta = 0 adds nothing; a
  1x1 box has an empty pair universe.
- **Typed rejections**: 13 refusal cases (empty box, zero side, negative
  / NaN / ±inf beta, negative / NaN / inf s, wrong-length uniform field,
  uniform out of [0,1] three ways, zero steps) — each leaving the last
  valid realization untouched.
- **Throughput report** (eprintln-timed test, debug profile): below.

### Allocation probe

`report_long_range_case` in `benches/percolation_allocations.rs`
(`LONG_RANGE_PROCESS_ALLOC`, asserted): 1x256 (32,640 pairs), 1x512
(130,816 pairs), and 16x16 (32,640 pairs) boxes — 0.000 allocations and
0.000 bytes per realization (uniform draw + per-pair decision + G1
substrate drive + C0 compression + measurement) over 128 measured
realizations each.

### Throughput

(eprintln-timed test, debug profile, beta = 0.8, s = 1.0): 1x256 —
8,229 us/realization = 122 boxes/s (3.97M pair decisions/s); 1x512 —
32,154 us/realization = 31 boxes/s (4.07M pair decisions/s).

### Limitations

The trajectory vertex-stage activations cost an O(n^2) incidence scan on
the pair view (the substrate's `activate_vertex` semantics); at the
current box sizes this is a small constant next to the pair loop, but a
much larger box would want a specialized substrate entry point. The law
is the pure power-law member of Biskup's family — Eq. 1.7's
`1 - exp{-beta(1+r)^-s}` is NOT offered (first-order equivalent at small
`beta r^{-s}`; the difference is a documentation choice, not an
approximation inside the code). L2 chemical-distance regimes and L3
diameter transitions are delivered (see their sections); no percolation
(infinite-cluster) claim is made on finite boxes — the observable set is component
structure only.

## L2 chemical-distance regimes

Status: **Validated** — all stage-L2 gates pass 2026-10-24; independent
review passed 2026-10-24 (no Blocker; the review's two Majors — a
seed-fragile s=2.5 gate and a pre-run `occupied_edge_count` bug — are
fixed as follow-ups, along with the estimator-doc, public-overflow,
capacity-state, and wording minors; see
`docs/reviews/2026-10-24-l2-review.md`).

### Contract: graph distance on one L1 realization

The chemical distance `D(x, y)` is the unweighted graph distance on the
occupied edge set of one L1 realization (this layer has no edge weights;
weighted shortest paths are the separate F1/F2 track). One realization
produces:

- **Distance profile** (`ChemicalDistanceWorkspace::profile()`): one
  `f64` per box vertex, entry `y` = `min_{x in sources} D(x, y)` for a
  caller-chosen non-empty source set. One realization is ONE
  element-preserving vector sample of fixed shape `vertices` — never a
  scalar distance alone. Unreachable entries carry the F1 sentinel
  `f64::INFINITY`; reachability is the derived `is_finite` test.
- **Exact box diameter** (optional, `DiameterTracking::AllVertexBfs`):
  the maximum finite distance over all vertex pairs including the
  diagonal — the maximum eccentricity. On a connected box graph this is
  the graph diameter; on a disconnected realization it is the largest
  within-component distance (isolates contribute 0). Always finite,
  always defined.

### Occupancy: the L1 draw contract replayed in-layer

The occupied set is produced by replaying the L1 contract, never by
re-running the L1 substrate: one uniform per pair in ascending
consideration order (saturated and near-zero pairs included) with the
pair at index k occupied iff
`uniforms[k] < LongRangeBox::pair_probability(r, beta, s)`. The pinned
expression is evaluated once per ℓ1 distance into a table — the same
deterministic function and thus the same bits as a per-pair evaluation.
For identical uniform fields the realized graph is bitwise the one
`long_range_run` exposes through `occupied_pairs()`; the sampler draws the
same stream as `LongRangePercolationMC` at the same seed. Both parities
are validation gates.

The CSR adjacency is built in two passes over the pair universe: degrees
first, then edge placement into the prepared offsets. The
`ChemicalDistanceSampler` reserves the complete-graph worst case
(`2 * pairs` slots, the L1 budget convention) at construction, so
realizations perform zero steady-state allocations.

### Regime semantics (regime-correct, anchor-faithful)

The anchor is Biskup, *On the scaling of the chemical distance in
long-range percolation*, Ann. Probab. 32, 2938 (2004), archived at
`docs/ref/md/Biskup_2003_long-range_chemical_distance.md`. In `d = 1`
the boundaries fall at `s = 1` and `s = 2`; the statements differ by
observable, and this layer keeps them strictly apart:

| Regime | Proved statement | Gate here |
|---|---|---|
| `s in (1, 2)` | Typical distance `D(x,y) = (log|x-y|)^(Delta + o(1))`, `Delta^-1 = log2(2/s)` | Typical-distance profile brackets: growth, slower than any power, fitted exponent below and ordered by `Delta` |
| `s = 1` | Diameter `D_N = Theta(log N / log log N)`; typical distance "at most logarithmic", lower bound explicitly missing | Diameter upper envelope only (`<= ln N`), no typical-distance gate |
| `s < 1` | Sup distance `sup D(x,y) = d/(d-s)` a.s. (constant); box diameter converges to it | Box-diameter constant envelope, no typical-distance gate |
| `s > 2` | Box diameter `D_N` expected linear in N (proved in d = 1); typical distance linear with a positive time constant | Paired log-log slope plus `D/r = A + c/r` extrapolation on the typical distance |
| `s = 2` | Explicitly "not very well understood" | **No gate, no run** |

The gates are finite-box statistical brackets of infinite-lattice
asymptotics, not proofs; every reported exponent is quoted against the
paper value with its finite-size direction. Biskup's statements are for
`x, y` in the infinite cluster, so the estimator conditions each profile
entry on reachability and reports the reachable fraction — at `beta >= 1`
in `d = 1` every `r = 1` pair is surely occupied, the box graph is
connected, and the reachable fraction is identically 1 (asserted per
realization in the statistical gates).

### API

- `chemical_distance_run(box, beta, s, uniforms, sources, tracking, workspace)`:
  the RNG-free core. Validates `(beta, s)`, sources, the diameter cap,
  and the uniform field before touching the workspace; rejections leave
  the previous profile/adjacency intact. Returns
  `ChemicalDistanceOutcome { occupied_edges, source_component_size,
  unreachable_count, diameter: Option<usize> }`.
- `ChemicalDistanceWorkspace`: `prepare(vertices, pairs, distance_range)`
  (idempotent, reserve-only growth), `reset()`, and accessors
  `profile`, `adjacency_offsets`, `adjacency_neighbors`,
  `occupied_edge_count`.
- `ChemicalDistanceSampler::new(box, beta, s, sources)` with
  `realize(rng)` (profile) and `realize_with_diameter(rng)` (profile +
  exact diameter); accessors for geometry, parameters, sources, profile,
  uniforms, workspace, and the last outcome.
- `ChemicalDistanceEstimator::new(vertices)` with `add_profile`,
  `add_diameter`, and `estimate() -> Option<ChemicalDistanceEstimate>`
  (per-entry mean/SEM/reachable fraction, diameter mean/SEM; NaN/0.0
  conventions documented on the type). `estimate()` allocates its
  vectors at report time only.

There is no Carlo.rs adapter at this stage (as with F1/F2): the estimator
exposes the per-vertex C0-shaped vectors directly, and L3 composes the
diameter track. `DIAMETER_VERTEX_CAP = 4096` bounds the
`O(vertices * (vertices + m))` all-vertex diameter pass; the single-source
profile path has no cap beyond memory.

### Validation gates

`tests/physics/long_range_distances.rs` (9 tests):

- **Pre-run workspace state**: after `prepare` (and in a freshly built
  sampler) the adjacency is empty — `occupied_edge_count() == 0`,
  offsets all zero, no neighbor slots — while the profile is the
  documented all-INFINITY buffer. (Review M2 follow-up: `prepare` now
  reserves adjacency CAPACITY without resizing, so the accessor cannot
  report the pair count before a run.)
- **BFS and diameter correctness**: 3 configs — 1x16 at `beta = 0.3`
  (disconnections and INFINITY live), 1x24 at `beta = 1`, 5x5 at
  `(beta, s) = (2.2, 0.6)` — x 5 seeds x 1-3 source sets. The occupied
  set is checked three ways (production CSR, test-local law replay,
  L1 core `occupied_pairs`); the profile matches a test-local multi-source
  BFS entry-by-entry bitwise (35 profiles), the diameter matches a
  test-local all-pairs reference (35 diameters), and
  `source_component_size`/`unreachable_count` match the reference counts.
- **Estimator statistics**: 40 `beta = 0.3` realizations from one source;
  mean, unbiased SEM, and reachable fraction match manual arithmetic at
  every vertex; the far corner is unreachable in some realizations
  (minimum reachable fraction 0.40 — the conditioning is live), the
  source is always reachable, never-reached entries carry NaN mean and
  0.0 error (directly exercised with `beta = 0`), and the sources' NaN
  conventions are exercised.
- **Bitwise contract**: same seed twice gives bit-identical profile,
  outcome, and adjacency; the L1 adapter at the same seed realizes the
  same occupied graph; a different seed changes the profile.
- **(1, 2) typical distance** (`beta = 1`, `n = 2048`, 150 realizations):
  for `s = 1.3 / 1.5 / 1.8` the pooled two-direction mean grows strictly
  at every r in `[8, 512]` against combined SEM; the `D`-vs-`r` log-log
  slope is 0.172 / 0.219 / 0.329 with CI upper bounds below 0.5 (a
  finite-window empirical bound, not a proof of sub-polynomial growth);
  the fitted `ln D`-vs-`ln ln r` exponents are 0.658 / 0.838 / 1.258,
  each strictly below `Delta = 1.609 / 2.409 / 6.579` with CI-disjoint
  increasing ordering in `s` — the documented finite-window approach
  from below at `ln r <= 7`, never an equality claim.
- **s = 2.5 linear regime** (`n = 8192`, 24 realizations, tail
  `r in [256, 2048]`): the finite-window log-log slope sits at or below
  1, so the gate is PAIRED — each realization contributes its own slope
  and its own `D/r = A + c/r` fit, and the aggregate 95% CIs carry the
  claim. D = 105.6, 202.1, 404.7, 747.2; slope 0.949 +/- 0.043, CI
  (0.907, 0.992); `A = 0.3713 +/- 0.0199`, CI (0.351, 0.391), strictly
  inside `(0, 1)`; `c = 11.18 +/- 8.55`, CI (2.63, 19.73), strictly
  positive; the measured ratio is non-increasing (0.413 -> 0.395 ->
  0.395 -> 0.365) and strictly decreases over the window. The paired
  design was adopted after the independent review showed the previous
  single-aggregate R=12 gate failed on 3/10 fresh seed bases.
- **s = 0.5 diameter envelope** (`beta = 1`, n = 256/512/1024/2048,
  8 realizations each): a BOUNDEDNESS gate, not a measurement of the
  constant 2. Every realization gives diameter 3 at every size (the
  constant sequence is reported explicitly instead of a degenerate
  zero-residual fit); means stay in [2, 4] with no positive power growth.
  The box graph is a subgraph of Z^d without outside shortcut vertices,
  so finite-box diameters may sit above the infinite-lattice a.s. sup
  `d/(d-s) = 2`; Biskup's `D_N` statement is quoted qualitatively. No
  typical-distance claim.
- **s = 1 diameter upper envelope** (n = 512/1024/2048/4096,
  3 realizations each): diameters 5.0 / 5.33 / 5.33 / 6.0 stay below
  `ln N = 6.24 / 6.93 / 7.62 / 8.32`; the diameter-vs-N slope is
  0.079 +/- 0.096, CI (-0.017, 0.175) — consistent with the
  `Theta(log N / log log N)` envelope and far below linear, with NO
  lower-bound claim. `s = 2` is not run anywhere (documented above).
- **Typed rejections**: 18 refusal cases — non-finite/negative beta
  (value bits preserved), negative/NaN s, empty sources, out-of-range
  source (position and value carried), wrong-length uniforms, out-of-range
  uniform (index carried), the diameter cap on a 4097-vertex box, sampler
  source domains, zero-vertex estimator, profile length mismatch, and the
  three public dimension-overflow paths (`vertices + 1`,
  `distance_range + 1`, `2 * pairs`; all typed before any buffer
  mutation), plus a `usize::MAX` diameter accumulated without overflow
  (the square is taken in f64) — all with the state-untouched assertion
  for the workspace and the empty-estimator `None`.
- **Throughput report** (eprintln-timed test, debug profile): below.

### Allocation probe

`report_chemical_distance_case` in
`benches/percolation_allocations.rs` (`CHEMICAL_DISTANCE_ALLOC`,
asserted): 1x256 (32,640 pairs), 1x512 (130,816 pairs), and 16x16
(32,640 pairs) profile realizations — 0.000 allocations and 0.000 bytes
per realization over 128 measured realizations each (uniform replay +
CSR build + multi-source BFS + estimator add), plus the 1x256
all-vertex-diameter mode (eccentricity scratch reuse) at 0.000
allocations and 0.000 bytes over 128 measured realizations (both
counters asserted on every path).

### Throughput

(eprintln-timed test, debug profile, `beta = 1`, `s = 1.5`): 1x256 —
5,197 us/realization = 192 realizations/s (6.28M pair decisions/s);
1x1024 — 82,414 us/realization = 12.1 realizations/s (6.36M pair
decisions/s); the 1x256 all-vertex diameter mode takes 12.5-12.8 ms per
realization (78-80/s).

### Limitations

The draw contract forces one full pair scan per realization (`O(pairs)`
occupancy decisions regardless of the realized edge count), so the layer
is quadratic in the box side lengths in `d = 1`. The exact diameter needs
one BFS per vertex; `DIAMETER_VERTEX_CAP = 4096` is a hard guard-rail, not
a feasibility guarantee (`O(V^3)` in the complete-graph worst case), and
the profile path is not capped. The sampler reserves the complete-graph
adjacency as capacity only (`2 * pairs` slots): the n = 8192 statistical
gate reserves ~805 MB of address space (the realized sparse adjacency
touches only a small fraction), which is worth knowing for CI. Distances
are unweighted integers carried in `f64` (exact up to `2^53`, unreachable
at the documented cap) — weighted first-passage distances are F1/F2, not
this layer. All regime gates are statistical brackets on finite boxes:
the asymptotic exponents are quoted from Biskup, never claimed as measured
equalities, and `s = 2` is deliberately gated not at all. At `beta < 1`
profile means condition on reachability; the reachable fraction is
reported alongside so conditioning is never implicit. Validation
rejections leave the workspace untouched; a capacity failure during
re-preparation may leave earlier buffers re-prepared but never returns a
result. There is no MonteCarlo adapter or scheduler integration yet
(L3's layer).

## L3 diameter transitions

Status: **Validated** — all stage-L3 gates pass 2026-10-24;
independent review PASSED 2026-10-24 (no Blocker/Major; the four minors
are fixed as follow-ups; see `docs/reviews/2026-10-24-l3-review.md`).

### Contract: the box diameter across Biskup's regimes

L3 is an analysis layer over L2, not a second distance engine: every
realization is produced by
`ChemicalDistanceSampler::realize_with_diameter`, so the L1 draw
contract, the O(N^2) pair replay, the CSR build, the BFS passes, and
`DIAMETER_VERTEX_CAP = 4096` are exactly L2's. The recorded observable
is the exact finite-box graph diameter `D_N` = maximum finite chemical
distance over all vertex pairs (isolates contribute 0; on a
disconnected realization it is the largest within-component distance).

`DiameterRegime::classify(s)` fixes the d = 1 boundaries and makes the
stage's claim scope explicit:

| `s` | Regime | Diameter statement (anchor) | L3 gate |
|---|---|---|---|
| `s < 1` | `Bounded` | Eq. (1.8) sup = `d/(d-s)` a.s.; box graph may sit above it | bounded, no power growth |
| `s = 1` | `Logarithmic` | Eq. (1.9) `D_N = Theta(log N / log log N)` | upper envelope relevant to `ln N` plus bounded ratio |
| `1 < s < 2` | `Polylogarithmic` | diameter bounds are explicitly NOT sharp; only the typical distance is pinned | sub-polynomial growth only |
| `s = 2` | `Critical` | "not very well understood" | **no gate** (`supports_diameter_gate() == false`) |
| `s > 2` | `Linear` | Biskup expects `D_N` linear in N (proved in d=1) | `D_N/N` approaches a positive constant from above |

These are finite-box brackets of infinite-lattice asymptotics. The
`(1,2)` row is the one place where it is essential not to transfer the
sharp typical-distance law to the diameter: Biskup states that the
available diameter bounds are not sharp enough for a definitive
asymptotic statement, and the module/test docs repeat that limitation.

### API and draw contract

- `DiameterRegime::{classify, supports_diameter_gate, label}`.
- `measure_diameter_ladder(sizes, exponents, beta, realizations, seed)`
  measures `(s, N)` rows (`s` outer, `N` inner), each with
  `realizations >= 2` independent draws from one source at `N / 2`.
  Row `r` uses `seed + r + i` for realization `i`; the per-realization
  integer diameters and their sample mean/SEM are recorded on
  `DiameterRow` (`diameters`, `mean`, `standard_error`, `ratio`,
  `ratio_standard_error`). `DiameterLadder` keeps the rows in order and
  offers `row(s, N)` / `rows_for(s)`.
- Typed `DiameterLadderError`: empty size/exponent ladder, zero size,
  size above the L2 cap, fewer than two realizations, invalid exponent,
  and transparent geometry/L2 error wrappers. All validation happens
  before any sampler is built.

### Validation gates

`tests/physics/long_range_diameter.rs` (7 tests):

- **Regime classification**: all five boundaries, `None` for
  negative/NaN `s`, and `supports_diameter_gate() == false` exactly at
  `s = 2`.
- **Replay parity**: 10 ladder diameters on two small boxes match
  direct L2 `realize_with_diameter` draws bitwise at the documented
  row seeds.
- **Linear (`s = 2.5`)**: N = 128/256/512/1024/2048, 24 realizations.
  Means 57.0/102.5/219.7/396.8/750.9; ratios 0.445/0.400/0.429/0.387/
  0.367 — non-increasing within two combined SEM, strictly lower at the
  ends; log-log slope 0.939 +/- 0.021, CI (0.871, 1.007); the `1/N`
  extrapolation `D_N/N = A + c/N` gives A = 0.3805 +/- 0.0153, CI
  (0.332, 0.429), strictly inside (0, 1), with c = 8.35 (the finite-size
  approach from above).
- **(1, 2) (`s = 1.5`)**: N = 128..1024, 12 realizations; diameters
  9.92/10.92/12.08/13.83, strictly increasing; log-log slope
  0.159 +/- 0.009 with upper CI 0.199, far below 0.5 — a finite-window
  sub-polynomial bracket, never a sharp exponent.
- **Logarithmic boundary (`s = 1`)**: N = 512/1024/2048, 4 realizations;
  diameters 5.0/5.0/5.5 against `ln N` = 6.24/6.93/7.63; every mean is
  below `ln N`, the ratio to `ln N / ln ln N` is 1.47/1.40/1.47 (bounded),
  and the sequence is non-decreasing within one integer step.
- **Bounded (`s = 0.5`)**: N = 256/512/1024, 4 realizations; the diameter
  is exactly 3 at every size, inside [2, 4], with no growth.
- **Typed rejections**: the seven refusal classes above.

### Limitations

L3 is one-dimensional (`LongRangeBox::new(&[N])`): Biskup's regime
statements and the recorded gates are the d = 1 ones; a d >= 2 ladder
would need its own anchors and is not gated here. Every realization pays
L2's all-vertex diameter cost (`O(N (N + m))`, `O(N^3)` in the
complete-graph worst case) and is capped at 4096 vertices; the tests use
N <= 2048 for that reason. The `s = 2` critical case is deliberately
ungated, and in `(1, 2)` no sharp diameter exponent is claimed. `D_N` is
integer-valued, so individual N steps can tie; the s = 1 gate therefore
allows a one-integer-step tolerance rather than strict per-step growth.
The layer has no Carlo.rs adapter (as with L2; observables are exposed
through the ladder records).

## C correlated percolation (C1–C3)

Status: **C1 and C2 Validated; C3 Experimental.** Independent review
(docs/reviews/2026-10-24-cq-review.md) passed C1; the C2 bootstrap-CI,
anchor, and estimator-narrative Majors are fixed, and the 256-replica
ladder gives `p_c(inf)=0.5457`, CI (0.5355, 0.5558), ~0.65 sigma from
Pruessner's 0.54299. C3's analytic map is exact and its `nu` estimator is
now calibrated on the exact iid control, but the correlated branch flip
is not statistically certified at L <= 1024 (a precise negative result).

### Anchors and conventions

Schrenk et al., PRE 88, 052102 (2013): `S(f) ~ |f|^{-beta_c}`,
`beta_c = 2(H+1)`, `c(r) ~ |r|^{2H}`. Pruessner et al., PRE 97, 062141
(2018): `C(r) ~ |r|^{-a}`, target `C(r) = (1+r^2)^{-a/2}`, discrete FFM
with zero-cutoff, and the threshold table used by C2
(a=inf 0.592746, a=1.5 0.54299(1), ..., a=0.5 0.5239(3)). The module
maps `a = -2H`. Consequence: `a -> infinity` is white noise / iid;
`a -> 0` is a constant field, NOT iid. This corrects the roadmap's
earlier `a = 0` wording.

### C1 — correlated fields

`CorrelatedGaussianField` builds the design correlation
`C(r) = (1+r^2)^{-a/2}` on a periodic power-of-two lattice, forward-
transforms it, applies the zero-cutoff, normalizes to unit site variance,
draws the Hermitian modes with a fixed Box-Muller contract, and
inverse-transforms. `CorrelatedSiteLaw` occupies `f_x > Phi^{-1}(1-p)`,
so the marginal is exactly p; `a = infinity` delegates bitwise to
`SiteBernoulli`; `a = 0` is the constant-field limit. The adapter
`CorrelatedPercolationMC` reuses the F3 analyzer and `ObservablePlan`.
Unit gates cover the design correlation, unit variance, marginal,
super-binomial fraction variance, adapter parity, determinism, the two
limits, and typed rejections. The thresholded binary correlation is
gated separately: a=1.5 gives 0.0993 at r=1 and 0.0212 at r=4, a=2.5
gives 0.0049 at r=4, and white noise gives -0.0004/-0.0002 — positive,
decaying, exponent-ordered, and vanishing in the iid limit.

### C2 — shifted thresholds

`ThresholdLadder` evaluates a probability grid on periodic `L x L`
lattices by sorting the field once per replica and activating vertices in
descending order through one incremental union-find pass; it records the
order-parameter susceptibility, the F3 finite-cluster susceptibility and
the largest-cluster fraction, then fits `p_c(L) = p_c(inf) + A L^{-1/nu}`
with bootstrap intervals and a log-derivative `nu` fit. The fast ladder
is a smoke test; the quantitative gate is the ignored long test at
L=128/256/512 with 256 replicas. The peaks are 0.5139/0.5268/0.5344 and
the `L^{-1/nu}` extrapolation with nu=4/3 gives `p_c(inf)=0.5457`,
CI (0.5355, 0.5558), chi2/dof=0.00. The 0.0027 offset from Pruessner's
0.54299(1) is ~0.65 bootstrap sigma, so the reference value lies inside
the CI; a +/-0.005 tolerance is also retained. Note the observable
difference: Pruessner uses a wrapping threshold on L=64-8192 with 10^5
replicas, while this estimator is the order-parameter susceptibility
peak `N Var(S_max/N)`. The
largest-fraction crossing estimator is separately calibrated on white
noise (0.5987 vs 0.592746). C2 remains Experimental pending independent
review.

### C3 — Weinrib–Halperin map

The analytic map implements `a nu_iid = 2`, `a_H = 2/nu_iid`, and the
two branches (`nu_iid` above, `2/a` below). For the 2D square,
`nu_iid = 4/3`, `a_H = 1.5`, and the predicted sequence is
4/3, 4/3, 8/5, 2, 8/3, 4 at a = 2, 1.5, 1.25, 1, 0.75, 0.5 — exactly
Pruessner's `nu_a` column.

The earlier log-derivative `nu` estimator is diagnosed and retired. Its
maximum of the central difference of `ln chi_op` is selected on a noisy
flank with a one-grid-spacing difference window, so it is not scale
covariant: on the exact iid control it returns 1.215 (0.739, 1.690) with
`chi2/dof = 11.7` on a pooled ladder, and 3.04 on the coarse default
window. The old `0.934/1.020/0.804` values are that failure mode, not a
physical branch shift.

The replacement is `ThresholdLadder::nu_from_pc_shift` on a
`ThresholdLadder::pool` of independent ladders: refine each order-
parameter susceptibility peak with a local log-quadratic fit over
`PEAK_FIT_RADIUS = 0.004`, then fit
`p_c(inf) - p_c(L) = A L^{-1/nu}` with `p_c(inf)` fixed. Eight pooled
ladders of 512 replicas (4096 per size) at L=128/256/512/1024 give:

| a | fixed p_c(inf) | refined peaks L=128..1024 | nu (95% CI) | predicted |
|---|---|---|---|---|
| +inf (control) | 0.592746 exact | 0.582975, 0.586747, 0.589137, 0.590577 | 1.379 (1.325, 1.433), chi2/dof 0.11 | 4/3 |
| 1.5 | 0.54299 | 0.51481, 0.52226, 0.53057, 0.53499 | 1.592 (1.413, 1.771), chi2/dof 1.42 | 4/3 |
| 1.25 | 0.53895 | 0.50223, 0.51571, 0.52093, 0.52754 | 1.819 (1.542, 2.095), chi2/dof 1.47 | 8/5 |
| 1.0 | 0.53452 | 0.50230, 0.50922, 0.51448, 0.52117 | 2.529 (1.741, 3.318), chi2/dof 0.58 | 2 |

The iid control now brackets 4/3 (bias +3.4%), which the old estimator
could not do. The correlated values are ordered correctly but all sit
high, increasingly so as `a` decreases. The fit is prior-dominated:
`d nu / d p_c(inf) ~ 0.12` per 0.001, so an effective `p_c(inf)` only
0.002-0.003 below the table would bring all three rows onto the
predictions; that offset is comparable to the finite-L zero-cutoff /
reference-scale spread C2 already records, but it is not independently
measured. A free `(p_c(inf), A, nu)` fit is ill-conditioned here (on the
iid control it returns 0.5932/1.43 instead of 0.592746/1.333). The gap
`nu(1.25) - nu(1.5) = 0.226 +/- 0.208` is 1.1 sigma, and the a=1.0
interval (1.741, 3.318) contains 2.

**Status: Experimental, flip NOT certified.** The estimator is
calibrated and usable, but L <= 1024 with a fixed table `p_c(inf)` cannot
distinguish the branches at 95%. A one-off cross-check extending to
L=2048 (L=256/512/1024/2048, four pooled ladders of 128 replicas) gives
`nu(1.5)=1.620` (CI 1.356-1.884, chi2/dof 1.76) and an iid control that
ranges 1.334-1.512 with the fit radius, so the larger box does not
certify the flip either. Certification needs L >= 2048-8192
(Pruessner used 8192 with 10^5 replicas) to suppress the zero-cutoff and
correction terms, or an independent 0.001-level `p_c(inf)` for this
finite-L FFM (a Binder-cumulant or very-large-L crossing), or both. The
nightly control gate
`c3_nu_pc_shift_calibrates_on_the_white_noise_control` asserts only the
iid calibration; no branch-value gate is claimed.

### Limitations

C1 gates the Gaussian field correlation, not yet the thresholded binary
correlation. C2/C3's literature comparison is currently a recorded
limitation, not a validated bracket. The FFT is restricted to
power-of-two axes; `a = 0` has no percolation transition. No Carlo.rs
adapter is provided for C2/C3 (analysis layers); C1 has one.

## Q generic-q random cluster (Q1–Q3)

Status: **Q1 and Q2 Validated; Q3 Experimental.** Q1's exact-FK
marginals (q=0.5-3), q=1 reduction, fuzz, and rate-level detailed-balance
audit pass; the q<1 shortcut Blocker found by the independent review is
fixed. Q2 passes its invariance/R3 comparisons and an exact joint
stationarity gate (706/1458 states at q=2/3, max deviation 1.3e-17);
Q3's phase scan is explicitly exploratory.

### Anchors

Geyer, J. Stat. Mech. (2015) P08010: Sweeny heat-bath Eq. (4),
Metropolis Eq. (5), the `1/2 <= pHB/pMH < 1` bound, and the connectivity
cost discussion. Fortuin-Kasteleyn (1972) defines the measure; the exact
2D self-dual curve `p_c(q) = sqrt(q)/(1+sqrt(q))` is the Q3 anchor (Wu,
Rev. Mod. Phys. 54, 235 (1982); see the C/Q reference notes).

### Q1 — Sweeny sampler

`SweenyMC` proposes every edge once per sweep in ascending order, one
uniform per proposal, with heat-bath or Metropolis acceptance. The forest
supports insertion by union-by-size, deletion of non-tree edges in O(1),
and tree-edge removal by an interleaved two-sided BFS whose cost is
proportional to the smaller side (worst case O(V)); it is not the
asymptotically optimal dynamic-connectivity structure and says so.
Gates: exact `2^E` enumeration on a 3x3 square for q=1,2,3 (own
union-find, full Z) matches the empirical |A| and k(A) marginals with
|z| <= 1.36; q=1 one all-closed heat-bath sweep equals the F2
`BondBernoulli` draw bitwise; update fuzz against full recomputation;
the acceptance-ratio bound; typed rejections; allocation note; and a
rate-level detailed-balance audit over every mask/edge on the 3x3 for
q=1,2,3, matching the exact FK two-state conditional to 1e-12.

### Q2 — Chayes–Machta for integer q

`ChayesMachtaMC` activates one colour class per step, resamples its
clusters, and updates awake bonds with independent `p` draws; the colour
index cycles deterministically. At p=0.45 on a periodic 3x3, q=2/3 agree
with Q1 within 4 sigma (differences 0.075/0.0007 vs noise 0.173/0.164),
the exact-FK |A| histograms have max |z| 2.96/3.74, and the colour
marginals have max |z| 1.48/2.93. At q=2 the CM chain also agrees with
R3's Edwards–Sokal chain (difference 0.041 vs noise 0.191). Exact joint
stationarity is gated: on a periodic 2x2 the full (bond-mask, colouring)
space has 706 states at q=2 and 1458 at q=3, and one CM step preserves
the Edwards-Sokal measure to `max |mu P - mu| <= 1.3e-17` for every active
colour.

### Q3 — q phase diagram (exploratory)

`critical_probability(q)` returns the exact self-dual curve (validated at
q=1,2,4,8). `phase_scan` runs Q1 at that probability on a periodic
square and reports bond-density mean/variance, skewness/kurtosis and the
bimodality coefficient. At L=8 the exploratory values are 0.362 (q=2)
and 0.505 (q=8). The ignored nightly L=16 scan gives 0.363/0.400/0.519/
0.654 for q=2/4/6/8: below the 5/9 bimodality threshold for q<=4, above
it for q=8, with q=6 unresolved at L=16 and monotone ordering in q. No
latent-heat or finite-size-scaling claim is made.
A production phase diagram needs larger L, a two-mode estimator and a
first-order/continuous classification gate — recorded as future work.

### Limitations

The Sweeny forest is not asymptotically optimal; Q2/Q3 use the same
connectivity substrate. Q3's bimodality coefficient is a descriptive
diagnostic, not a latent heat. Q2's detailed balance is argued
algorithmically but not yet exhaustively audited; Q1 has the exact-FK
and fuzz gates. No C0 array adapter is added for the QR tracks yet.