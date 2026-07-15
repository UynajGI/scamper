# Priority 2 cluster-solver review

Reviewed scope: `QMC.rs/src/impurity/spin_boson/cluster` only.

## Physics checks

- Hamiltonian convention is consistently
  `H = -(Delta/2) sigma_x + (epsilon/2) sigma_z + sigma_z sum g_k(b_k+b_k^dagger) + sum omega_k n_k`.
- Poisson auxiliary-cut rate `Delta/2` is consistent with the continuous-time expansion of the transverse term.
- Positive bias correctly favors `sigma_z=-1`: the cluster heat bath uses
  `p(+)=1/(1+exp(epsilon L))`, and the conditional mean is `-tanh(epsilon L/2)`.
- Retarded bonds are only proposed between equal-spin segments with
  `p=1-exp(-2 integral_I integral_J K_beta)`, consistent with the ferromagnetic retarded Ising representation used by this implementation.
- Single-mode normalization `lambda=g^2/omega` is consistent with the code's direct `sigma_z` coupling convention.
- The expansion-order estimator `E_x=-<N_kink>/beta` and
  `<sigma_x>=2<N_kink>/(beta Delta)` are mutually consistent.
- Existing free-spin and truncated single-mode ED tests exercise the most important normalization and sign conventions.

## No correctness change made

The initially suspected finite-bias sign mismatch is not present. The sampling probability and improved estimator have the same convention, and the existing positive-bias regression test covers it.

## Remaining engineering limitation

Retarded bond construction currently scans every segment pair, so one sweep costs `O(N_segment^2)`. This is correct but may dominate at large beta or large tunnelling. A future performance stage can add stochastic bond generation / cumulative-kernel sampling without changing the public solver contract.
