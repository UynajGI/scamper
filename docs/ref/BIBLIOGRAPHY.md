# Percolation Platform Reference Bibliography

Every physics or algorithm correctness claim in the percolation platform is
anchored to a published paper. The paper is archived here as the arXiv PDF,
converted to Markdown under `md/` (incrementally, via `convert_to_md.py`), and
the comparison evidence lives in the stage section of `CMC.rs/PERCOLATION.md`.
A claim counts as correct only when our measurement agrees with the paper's
result within its stated statistical convention.

| Stage | Paper | Local copy | What is compared | Agreement |
|---|---|---|---|---|
| Z1 Newman-Ziff scans | M. E. J. Newman, R. M. Ziff, *Efficient Monte Carlo Algorithm and High-Precision Results for Percolation*, Phys. Rev. Lett. **85**, 4104 (2000) | [pdf](NewmanZiff_2000_PRL_percolation.pdf) · [md](md/NewmanZiff_2000_PRL_percolation.md) · arXiv:cond-mat/0005264 | Incremental-activation scan + binomial reweighting reproduces direct Bernoulli estimates; square site threshold `p_c = 0.592746...` anchors benchmark operating points (`p_site = 0.5927`) | Scan-vs-direct CI overlap at 17 comparison `p` values; 1.83-4.65x scan speedup (PERCOLATION.md, Z1) |
| N2 random failure (giant fraction) | D. S. Callaway, M. E. J. Newman, S. H. Strogatz, D. J. Watts, *Network Robustness and Fragility: Percolation on Random Graphs*, Phys. Rev. Lett. **85**, 5468 (2000) | [pdf](Callaway_2000_PRL_network_robustness.pdf) · [md](md/Callaway_2000_PRL_network_robustness.md) · arXiv:cond-mat/0007300 | Eq. (5) generating-function giant fraction `S = q(1-G0(u))`, `u = 1-q+q G1(u)`, evaluated on the graph's own degree distribution | ER N=1e5, z=2.5: `\|S_mc - S_theory\| <= 0.0005` at q = 0.5/0.65/0.8 (PERCOLATION.md, N2) |
| N2 random failure (threshold) | R. Cohen, K. Erez, D. ben-Avraham, S. Havlin, *Resilience of the Internet to Random Breakdowns*, Phys. Rev. Lett. **85**, 4626 (2000) | [pdf](Cohen_2000_PRL_internet_resilience.pdf) · [md](md/Cohen_2000_PRL_internet_resilience.md) · arXiv:cond-mat/0007048 | Critical keep fraction `q_c = 1/(kappa-1)`, `kappa = <k^2>/<k>`; used as a statistical bracket (never an exact finite-N threshold) | Empirical kappa = 3.5025 -> q_c = 0.3996 vs Poisson ideal 0.4; giant extinguished below, present + theory-consistent above (PERCOLATION.md, N2) |

| N3 targeted attack (breakdown threshold) | R. Cohen, K. Erez, D. ben-Avraham, S. Havlin, *Breakdown of the Internet under Intentional Attack*, Phys. Rev. Lett. **86**, 3682 (2001) | [pdf](Cohen_2001_PRL_intentional_attack.pdf) · [md](md/Cohen_2001_PRL_intentional_attack.md) · arXiv:cond-mat/0010251 | Their Eqs. (2), (7), (9) attacked-distribution estimate with the `kappa = 2` criterion give the breakdown fraction `f_c`; used as a statistical bracket | `f_c(power-law alpha=2.5) << f_c(ER)`, and the measured giant survives at `0.5 f_c` / collapses at `1.5 f_c` (PERCOLATION.md, N3) |
| N3 targeted attack (cutoff curve values) | D. S. Callaway, M. E. J. Newman, S. H. Strogatz, D. J. Watts, PRL **85**, 5468 (2000) | [pdf](Callaway_2000_PRL_network_robustness.pdf) · [md](md/Callaway_2000_PRL_network_robustness.md) · arXiv:cond-mat/0007300 | Their Eq. (14) `S = F0(1) - F0(u)` with degree-cutoff occupation `q_k = 1{k <= K}` at complete degree-group boundaries (an exact configuration-model prediction) | Measured vs theory at 4+ cutoffs on ER and power-law configuration graphs (PERCOLATION.md, N3) |

| D3 layered directed percolation (threshold + decay exponent) | H. Hinrichsen, *Nonequilibrium Critical Phenomena and Phase Transitions into Absorbing States*, Adv. Phys. **49**, 815 (2000) | [pdf](Hinrichsen_2000_AdvPhys_absorbing_states.pdf) · [md](md/Hinrichsen_2000_AdvPhys_absorbing_states.md) · arXiv:cond-mat/0001070 | Table 1: bond DP (1+1)D `p_c = 0.6447001(1)`; Table 2: `delta = 0.159464(6)` | Survival P(64)/P(1024) and full-seed density rho(64)/rho(1024) both recover delta within 0.05-0.06 at 2000 realizations; sub/super-critical brackets (PERCOLATION.md, D3) |

| E2 bootstrap percolation (metastability threshold) | A. E. Holroyd, *Sharp Metastability Threshold for Two-Dimensional Bootstrap Percolation*, Probab. Theory Relat. Fields **125**, 195-224 (2003) | [pdf](Holroyd_2003_PTRF_bootstrap_sharp_threshold.pdf) · [md](md/Holroyd_2003_PTRF_bootstrap_sharp_threshold.md) · arXiv:math/0206132 | Theorem 1: `I(L,p) -> 1` when `p ln L > lambda` and `-> 0` when `p ln L < lambda`, with `lambda = pi^2/18` — i.e. finite-size seed threshold `p_c(L) ~ pi^2/(18 ln L)`; used as a statistical bracket (never an exact finite-L threshold) | On the L x L torus (r = 2) the majority-activation probability at `p* = pi^2/(18 ln L)` stays strictly inside (0.05, 0.95) for L = 32/64/128 and is bracketed by `p*/1.5` (low) and `1.5 p*` (high) with real error bars (PERCOLATION.md, E2) |

| E3 k-core percolation (emergence threshold) | S. Janson, M. J. Luczak, *A Simple Solution to the k-Core Problem*, Random Structures Algorithms **30**, 50-62 (2007) | [pdf](Janson_2007_kcore_simple_solution.pdf) · [md](md/Janson_2007_kcore_simple_solution.md) · arXiv:math/0508453 | Theorem 1.1 (after Pittel, Spencer & Wormald): for G(n, m) with mean degree `lambda = 2m/n -> c`, the k-core (k >= 3) is empty whp below `lambda_k = min_{mu>0} mu/psi_{k-1}(mu)` and giant above it with `v(Core_k)/n -> psi_k(mu_k(lambda))`, where `psi_j(mu) = P(Poisson(mu) >= j)` and `mu_k(lambda)` is the largest solution of `mu = lambda psi_{k-1}(mu)`; own solver in the test; used as a statistical bracket with real error bars | Solver reproduces `lambda_3 = 3.350919`; G(n=5e4): at `c = lambda_3 - 0.15` the 3-core fraction is consistent with zero, at `c = lambda_3 + 0.15` it agrees with `psi_3(mu_3(c)) ~ 0.433` within 4 sigma + finite-N slack (PERCOLATION.md, E3) |

## Conversion

`convert_to_md.py` converts new PDFs incrementally with PyMuPDF (any local
pdf→md tool is acceptable; this replaced the MinerU-only flow until a token
is configured):

```bash
python3 -m venv /tmp/pdfenv && /tmp/pdfenv/bin/pip install pymupdf
/tmp/pdfenv/bin/python docs/ref/convert_to_md.py
```

The MinerU-based `literature-extract` skill remains the higher-fidelity path
for tables/formulas once `mineru-open-api` and a token are available.
