# Explicit boson-occupation solver

Location: `QMC.rs/src/impurity/spin_boson/occupation/`

This backend targets one or a few explicitly retained cavity modes with finite occupation cutoffs. It supports the Rabi and Jaynes-Cummings Hamiltonians and keeps the photon occupations in the sampled configuration.

## Algorithm

The solver diagonalizes the finite-cutoff Hamiltonian once and constructs the exact link propagator

`T = exp[-(beta/M)(H-E0)]`.

A complete periodic occupation worldline is then drawn with a closed-path bridge heat bath. Since `T^M = exp[-beta(H-E0)]`, changing `M` does not introduce Trotter error. The controlled approximation is the bosonic occupation cutoff. The method is intended for single/few-mode benchmarks; its dense setup scales cubically in the total finite Hilbert-space dimension.

The complete-path bridge update is important for ergodicity: local single-slice updates alone cannot cross conserved Jaynes-Cummings excitation sectors or Rabi parity sectors.

## Measurements

- total and mode-resolved photon occupations;
- `<n_m^2>`, `<n_m(n_m-1)>`, and `g_m^(2)(0)`;
- mode-mode equal-time number correlations;
- `sigma_z` and an off-diagonal `sigma_x` transfer estimator;
- energy transfer estimator;
- Rabi parity `sigma_z (-1)^(sum n)`;
- connected `sigma_z`--mode-number covariance;
- finite-cutoff reduced-spin purity as an entanglement/mixing proxy.

## Deliberate boundary

This is not a replacement for the retarded wormhole or continuous-time cluster solvers. It is a small-Hilbert-space cavity benchmark backend. Occupation-number SSE and directed-loop implementations can later live beside it under the same `occupation/` namespace.
