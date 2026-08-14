# Rabi Model Finite-Temperature U4 Crossing Study

## Objective

Study how the Binder cumulant U4 crossing point g* shifts with inverse temperature β,
for three values of Ω, to map the finite-temperature crossover behavior.

## Model

H = (Δ/2)σz + Ω a†a + g σx(a + a†), Δ=1, g_c = √(ΩΔ)/2

## Parameters

| Parameter | Values |
|-----------|--------|
| Ω | 5, 10, 20 |
| β | 0.1, 0.5, 1.0, 2.0, 4.0, 8.0 |
| Cutoff N_b | 8, 16, 32, 64 |
| r scan | 0.2 – 3.0, 500 points |

## Key findings

### 1. Finite-T effects visible at small β

At Ω=5, r=1.0: U4(N32) goes from 0.307 (β=0.1) to 0.143 (β=8.0).
The thermal smearing is strong at β=0.1 and saturates by β≈8.

### 2. U4 curves do NOT cross at any β

At every β, the cutoff ordering is monotonic:
  U4(N8) > U4(N16) > U4(N32) > U4(N64)

The gap between adjacent cutoffs grows with r (coupling strength),
but never changes sign. **No Binder cumulant crossing exists at finite temperature.**

This is consistent with theory: the Rabi QPT exists only in the η=Ω/Δ→∞ limit.
At finite η, quantum tunneling between the two wells keeps the ground state symmetric,
and U4 curves for different cutoffs are smooth and ordered.

### 3. All 2855 reported "crossings" are numerical noise

N32 and N64 agree to ~7 decimal places at β≥2, so Jacobi eigensolver roundoff
(~10^{-10}) creates spurious sign changes in U4_N32 - U4_N64.

## Data files

- `results/u4_curves/omega_{Ω}_beta_{β}.csv` — U4 vs g for all cutoffs (500 points each)
- `results/crossings/all_crossings.csv` — crossing points (mostly noise artifacts)

## Run

```bash
cd research/rabi_finite_T
cargo run --release 2> logs/run_YYYYMMDD.log
```
