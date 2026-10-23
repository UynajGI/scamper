# T2 (Heterogeneous Newman-Ziff Reweighting) — Independent Reviewer Report

> Independent review of commit `7a9c5ef` ("feat(CMC): add heterogeneous
> Newman-Ziff reweighting (T2)", branch `feat/cmc-percolation`), 2026-10-22.
> Reviewer did not write the code; all numbers below are the reviewer's own
> computations unless flagged as a re-run of a committed gate. No code was
> modified. Host: the same Xeon Gold 6148-class Linux box, `rustc 1.98.0
> (88d9e12ae 2026-08-18)`, `cargo 1.98.0`.

**VERDICT: PASS** (three minor/notes for follow-up; none block)

Scope reviewed: `CMC.rs/src/percolation/scan/reweight.rs`
(`reweight_heterogeneous`, `HetReweightResult`, `HetReweightWorkspace`, new
`ReweightError` variants, factored `validate_curve`), `scan/mod.rs` exports +
crate-root re-export in `src/lib.rs`, `tests/physics/percolation_scan_heterogeneous.rs`,
the `HET_SCAN_ALLOC` probe in `benches/percolation_allocations.rs`, the `nz_het`
group of `benches/percolation_newman_ziff.rs`, and the PERCOLATION.md T2 row +
T2 sections under Z1. `git show 7a9c5ef` confirms `newman_ziff.rs` is untouched
by the commit, so the Z1 prefix-parity gate carries over to the unchanged scan
layer and the roadmap's "prefix parity as v1" gate is satisfied by
construction.

## Required verifications

| # | Claim under review | Result | Reviewer's independent evidence |
|---|---|---|---|
| 1 | Estimator identity `E_het[f] = Σ_n C(N,n) E_unif,n[w·f]`, per-scan raw sum unbiased, unit mass exactly 1 in expectation | **REPRODUCED (exact)** | From-scratch re-derivation (below) + exact rational brute force over ALL N! permutations at N=4 and N=5, random rational fields, f = \|A\| and f = largest component on a random graph: `E_perm[V] == 2^N-enumeration` as exact `Fraction`s, and `E_perm[S] == 1` exactly, all 4 cases. Stronger than machine precision. |
| 2 | Raw sum returned, not per-scan ratio v/s, because the ratio is biased | **REPRODUCED** | Own MC, own graph (G(9,0.4), 12 edges) and own spread field, own permutation sampler + own union-find, R = 200 000: RAW estimator z = **+0.12** vs exact 2^9 enumeration; RATIO estimator v/s z = **−24.01** (stage reported −17 at N=9/R=200k — same phenomenon, magnitude depends on the field). Production returns the raw sum. Uniform-field control: with exact `Fraction`s, `s = Σ_n C(N,n) p^n (1−p)^{N−n} = 1` **identically on 500/500 random permutations** — the v1 `reweight` denominator depends only on (N, p), never on the trajectory, so v1 has no ratio-bias exposure. |
| 3 | Log-space numerics: prefix log-likelihood, ln(k)-table binomial recurrence, running-max log-sum-exp; stability envelope | **REPRODUCED with one doc gap** | Line-by-line check of reweight.rs:366–398 against the derivation: all three recurrences correct (see below). Own Python IEEE-double mirror at N=10^4: constant field p=0.3/0.5927 vs two independent references (v1 mode-anchored recurrence; lgamma log-space) — rel. 8.7e-10 / 1.65e-9, and the production crate reproduces the same magnitudes (8.673e-10 / 1.651e-9). Spread U(0.3,0.8): per-scan unit mass ~1e-203 (log-max ≈ −467), finite, 6217/10001 terms underflow harmlessly relative to the max. Near-degenerate half-0.01/half-0.99 (envelope σ = 229.8): **production returns value = 0, unit_mass = 0** (finite) on ordinary scans — see Minor finding M1 on the doc's one-sided saturation story. |
| 4 | Variance envelope: per-scan log-values fluctuate with sd ≈ `sqrt((N/4) Var(logit p))`, roughly lognormal | **REPRODUCED** | Re-derived: for a uniform size-n subset, `Var(Σ_{i∈A} logit p_i) = n(N−n)/(N−1)·Var_pop(logit p)` → `N²/(4(N−1)) ≈ N/4` at n = N/2. Own MC on the unit-mass estimator (f = 1, pure weight fluctuation) over N ∈ {64,256,1024,4096} × 5 fields: measured sd(ln S) = 0.0347/0.1324/2.311/… vs formula 0.0383/0.1435/2.516/… — the formula tracks every sqrt(N) and sqrt(spread) step and is 4–10% **conservative** (overestimates the tail). Uniform field: sd(ln S) = 0 exactly. Measured cv(S) matches the lognormal prediction `sqrt(exp(σ²)−1)` inside the envelope (e.g. ramp h=0.03 at N=4096: measured 1.63 vs 1.65). The envelope is drawn honestly. |
| 5 | Gates re-run (release + debug) | **REPRODUCED** | Note: the tests live under the `suite` target (`tests/suite.rs` wires `#[path]` modules), so the invocation is `cargo test -p cmc-rs --test suite percolation_scan_heterogeneous …` (as PERCOLATION.md documents). Release: **5/5 passed**, exact-enumeration z's −1.135…+1.354, **max \|z\| = 1.354** (doc: 1.35, bimodal site RawSecondMoment), unit-mass z +0.035 at L=128 (doc: +0.04); at-scale L=64 z = −1.205/−0.802, L=128 z = −0.646/−0.676 (doc: −1.21/−0.80/−0.65/−0.68). Debug: **5/5 passed** (27.2 s). Every documented z-score regenerated to two decimals. |
| 6 | Uniform reduction + typed rejections | **REPRODUCED** | Own scan (6×6 square, own seed): owned constant field vs v1 `reweight` across p ∈ {0.1, 0.25, 0.5, 0.5927, 0.75, 0.9}: max relative difference **2.3e-14** (≪ 1e-9 claimed), unit_mass = 1 ± 1e-15; `ProbabilityField::Uniform` delegation **bitwise** identical (`to_bits` equal) with unit_mass exactly 1.0. Rejections on my own curve: `ProbabilityFieldLength{expected:36, actual:34}`; `BoundaryProbability{index:7, value:0.0}`; `{index:3, value:1.0}`; `Probability::new(NaN/±inf)` errors; reviewer extra probe: a **−0.0** entry is caught (IEEE −0.0 == 0.0) at index 0 — boundary rejection is airtight. |
| 7 | Independent cross-check at scale (load-bearing physics gate) | **REPRODUCED** | Reviewer's own direct sampler — own splitmix64 RNG, own 53-bit unit draws, own weighted-union-find, own moment bookkeeping, sharing no code with `SiteBernoulli`/`analyze` — on the same open 64×64 square with the same graded ramp (half-width 0.004 centered at 0.5927), 16 384 direct samples vs 512 production scans+reweights with reviewer-derived seeds: giant fraction **z = −0.181** (0.285005 ± 0.005304 vs 0.285977 ± 0.000802), raw second moment **z = −0.259**, unit mass z = −0.509 vs 1. Bonus bond-mode gate (chain-10, ramp field over the edge universe, R = 200 000 reweighted vs 65 536 own direct): **z = −0.247**, unit-mass 1.000687 ± 0.001504. |
| 8 | Allocation probe + nz_het performance claims | **REPRODUCED** | `cargo bench -p cmc-rs --bench percolation_allocations --features allocation-probe -- --test`: **10/10 `HET_SCAN_ALLOC` records** (5 topologies × site/bond) with `allocations_per_sample=0.000`, `allocated_bytes_per_sample=0.000`, zero peak live; retained capacity matches the doc exactly (square site 65 544 B, bond 129 032 B — two f64 words per element). `nz_het` Criterion re-run (same short config): direct/reweight = square site **0.62x**, square bond **0.41x**, cubic site **0.77x**, cubic bond **0.40x**, sparse-ER site **1.28x**, sparse-ER bond **0.43x** vs doc 0.63/0.41/0.80/0.38/**1.40**/0.43 — every dense sub-1 and the single sparse-site above-1 reproduce in direction; the 1.40x point estimate re-measured 1.28x (10%-level short-run noise; the doc flags point-estimate status). |
| 9 | Docs honesty (PERCOLATION.md T2 row + T2 sections, rustdoc vs code) | **REPRODUCED with nits** | Status honestly **Experimental** with "No independent review yet". The T2-row claims (max z 1.35 @ R=3200, L=64/128 max z 1.21, 1e-9 reduction, typed rejections, 0-alloc on five topologies) all regenerated. The envelope caveat is prominent and matches my V4 numbers; sub-1 speedups are disclosed in the row itself. Rustdoc derivation, recurrences, "one ln per entry", and the raw-vs-ratio story all match the code (verified line by line). Nits N1–N3 below. |

## Independent computations (detail)

**Identity, from scratch.** Split the heterogeneous product measure by subset
size: `E_het[f] = Σ_A w(A) f(A) = Σ_n C(N,n) · (1/C(N,n)) Σ_{|A|=n} w(A) f(A)`.
A uniformly random permutation's size-n prefix is exactly a uniform n-subset,
so `w(A_n) f(A_n)` estimates `E_unif,n[w f]` for every n at once; the
per-trajectory sum `Σ_n C(N,n) w(A_n) f_n` is unbiased by linearity regardless
of the dependence between prefixes of one permutation. Verified **exactly**
(Fraction arithmetic, probabilities with denominator 97): at N = 4 and 5, the
average of `V(π)` over **all** N! permutations equals the full-2^N enumeration
value as exact rationals (e.g. N=5, largest-component f: both
`21413516033/8587340257`), and `E_perm[S] = 1` exactly, for both a trivial and
a graph functional. The committed hand-worked test was also re-derived on
paper (order [1,0]: mass 1.5, value 1.875; order [0,1]: 0.5/0.875; average
1.375 = exact heterogeneous mean).

**Bias of the per-scan ratio.** Own field/graph/sampler at N=9, R=200k (own
incremental union-find along each permutation): raw mean 0.394028 ± 0.00078 vs
exact 0.393934 (z = +0.12); self-normalized ratio mean 0.389308 ± 0.000193
(z = −24.01). The ratio's sem looks deceptively small — the bias is real and
tens of sigma. The production path never divides by the per-scan mass.
Uniform case: `S ≡ 1` identically (binomial theorem; 500/500 exact-rational
random permutations at N=6), so v1 `reweight`'s division by
`Σ_k B(k)` — a deterministic function of (N, p) only — cannot introduce ratio
bias; the two paths are structurally different in exactly the way the rustdoc
says.

**Numerics.** The three recurrences in reweight.rs:366–398:
`l_0 = Σ_j ln(1−p_j)` over the full universe (line 367–370); entering element
`e` adds `ln(p_e/(1−p_e))` (one transcendental, line 378);
`logC(N,n) += ln(N−n+1) − ln(n)` via the `ln(k)` table (`log_integers[trials−n]
− log_integers[n−1]`, line 379) — each is the correct log of its exact
counterpart; `log_weights[n] = l_n + logC(N,n)`; running max, `u_n ∈ (0,1]`,
final `exp(l_max)` rescale (383–398). The odds ratio cannot produce NaN/inf
for validated `p ∈ (0,1)`: `1−p` is exact for p ≥ 0.5 (Sterbenz) and the worst
underflow of the odds (p = min subnormal) still logs to ≈ −744. My IEEE-double
mirror matches the production crate's outputs at N=10^4 to the reported digits
(constant-field deviation vs v1: 8.730e-10 mirror vs 8.673e-10 production at
p=0.3 — the tiny difference is permutation-order noise, same magnitude).
Error budget: at N=10^4 the log-space partial sums reach magnitude ~7000, so
~2·10^4 accumulated roundings give ~1e-9 relative — exactly what is observed;
the 1e-12 unit-test tolerance is an N=10 statement, not an N=10^4 one (Note
N4). Saturation: across bimodal 0.01/0.99 fields from N=10^4 to 2·10^5 the
log-mass excursion stays hugely **negative** (−1.6·10^4 … −3.2·10^5): the
documented +709.8 overflow side is essentially unreachable, while the −745
underflow side occurs almost surely outside the envelope (M1).

**Envelope.** Own MC grid (20k/8k scans per point) on `S(π)` (f = 1):
sd(ln S) tracks `sqrt((N/4) Var_pop(logit p))` at 0.90–0.96× across every
(N, spread) combination — e.g. U(0.3,0.8): 2.31/4.78/9.73/19.63 at
N = 64/256/1024/4096 vs formula 2.52/5.01/10.02/20.03; ramps half=0.03:
0.13/0.27/0.56/1.13 vs 0.14/0.29/0.57/1.15. Uniform → exactly 0. cv matches
the lognormal prediction inside the envelope and is *under*-estimated by any
finite sample far outside it — i.e. the caveat is, if anything, optimistic in
the user's favor only inside σ ≲ 1, and the doc's rule of thumb ("competitive
while σ is of order one or below") is well-calibrated. The doc's specific
anchor numbers (bimodal σ = 57.3 at N=4096) regenerate exactly.

**Independent at-scale gate.** Fully independent implementation (own RNG,
own union-find, own sampler; only the graph object and the code under review
are shared): L=64 ramp — z = −0.181 (giant), −0.259 (moment); chain-10 bond
ramp at R = 200k — z = −0.247. Combined with the committed gates (which use
`SiteBernoulli`+F3 as reference and a test-local BFS for enumeration), the
physics is confirmed by two disjoint reference stacks.

## Findings

**Blockers: none. Majors: none.**

- **M1 (minor, docs)** — reweight.rs:296–300 rustdoc "Numerics" (echoed at
  PERCOLATION.md:1396–1398) tells only one side of the out-of-envelope story:
  `exp(l_max)` "saturates to infinity only when a single scan's log-mass
  excursion exceeds ln(f64::MAX) ~ 709.8". In practice that side is
  unreachable (reviewer's scan over bimodal fields up to N=2·10^5 never saw a
  positive excursion), while the **reachable** failure is the opposite:
  `exp(l_max)` underflows to 0.0 and the call returns `value = 0.0,
  unit_mass = 0.0` — verified through the production API on a half-0.01/
  half-0.99 field at N=10^4 (σ = 229.8, far outside the envelope). Not silent
  if the caller heeds `unit_mass` (its scan average collapses to 0 ≠ 1), and
  strictly outside the documented usable envelope, but a user scanning the
  "Numerics" section would expect `inf`, not a quiet zero. Suggest one added
  sentence: beyond ~−745 the per-scan value and unit mass underflow to
  exactly 0 — treat unit-mass averages of 0 as an envelope violation.

- **N1 (note, docs)** — PERCOLATION.md:1419 says "plus **7** unit tests in
  reweight.rs" and lists six; the commit adds exactly 6 new `#[test]`
  functions (verified by diff). Trivial count fix.

- **N2 (note, docs)** — PERCOLATION.md:1410–1412: "the bimodal logit variance
  is ~70x larger still" (than the wide-ramp cv-45 field). Not reproducible
  under any convention the reviewer tried: a uniform ±0.09 ramp at p_c gives
  Var(logit) = 0.0117 vs bimodal 3.2104 → **275x** (two-point convention 22x,
  cv-back-solve 432x). Direction and conclusion (bimodal is catastrophically
  outside the envelope; σ = 57.3 at N=4096) are correct and the error is
  conservative, but the multiplier as printed does not regenerate.

- **N3 (note, docs)** — the task brief's `cargo test -p cmc-rs --test
  percolation_scan_heterogeneous` is not a real target; the physics tests are
  wired through `tests/suite.rs` `#[path]` modules, so the working invocation
  is `cargo test -p cmc-rs --test suite percolation_scan_heterogeneous`
  (which is what PERCOLATION.md:1346 documents). No action needed beyond
  keeping the documented form.

- **N4 (note, precision envelope)** — the constant-field↔v1 agreement
  tolerance demonstrated by the unit test is 1e-12 at N=10, but at N=10^4 the
  same comparison is ~1.7e-9 (both the reviewer's mirror and the production
  crate). This is exactly the documented "accumulated-rounding scale of the
  log recurrences" (log-space partial sums reach ~7000 before the rescale),
  so the rustdoc is honest in kind — but no test pins the N-scaling, and a
  future contributor could mistake 1e-12 for the general guarantee. Consider
  one constant-field-vs-v1 assertion at N=10^4 with a ~1e-8 tolerance to pin
  it.

- **N5 (note, observed and dismissed)** — spread U(0.3,0.8) at N=10^4 gives
  per-scan unit masses around 1e-203 with log-max ≈ −467 (envelope σ ≈ 37):
  the estimator's mean is carried by roughly one-in-10^200 tail scans. This
  is the documented heavy tail behaving exactly as advertised (σ ≫ 1 → use
  N1), and the committed finiteness test explicitly labels itself "a
  numerics claim, not an accuracy claim". No action.

Additional things probed and dismissed: `prepare` shrink/regrow semantics
(`prepared_integers` refills only on trial-count change; content depends only
on the count — no stale reads; every `log_weights` slot is rewritten each
call); boundary scan order (first offending index in dense order, as tested);
odds-ratio overflow (impossible for validated p; `1−p` exact for p ≥ 0.5 by
Sterbenz); `−0.0` entries (caught: IEEE equality); NaN (unconstructible
through `Probability`); v1-path behavior after the `validate_curve` factoring
(diffed against the pre-commit body — identical checks in identical order,
endpoints p=0/1 unchanged); `reserve_len` arithmetic (over-reserves
`len − len()` when shrinking capacity is short — never under-reserves);
`HetReweightWorkspace` cloning; bench-probe independence (borrowed
heterogeneous field with alternating values, so the probe exercises the het
path, not the uniform delegation).

## Recommendation

**PASS.** The estimator is exactly unbiased (proven in exact rational
arithmetic, not just statistically), the deliberate raw-sum API design is
vindicated by a −24σ ratio bias in the alternative, the log-space
implementation is faithful to its derivation and stable across the documented
envelope, the envelope formula itself is honest and slightly conservative,
the 0-allocation claim holds, the performance story (including the sub-1
speedups) is disclosed accurately, and — the load-bearing check — a fully
independent sampler/union-find of the reviewer's own agrees with the
production reweighter at L=64 within z = 0.26 and at R = 200k on a bond-mode
chain within z = 0.25. The three minor/notes are documentation polish; M1 is
the only one worth applying before the stage is promoted past Experimental.
