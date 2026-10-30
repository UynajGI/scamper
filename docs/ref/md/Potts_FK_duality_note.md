# 2D Fortuin-Kasteleyn / Potts self-dual critical curve (Q3 anchor note)

Anchor for Q3. The `q`-state ferromagnetic Potts model on the square
lattice maps to the Fortuin-Kasteleyn random-cluster measure

```text
phi_{p,q}(A) = p^{|A|} (1-p)^{E-|A|} q^{k(A)} / Z(p,q)
```

(Fortuin & Kasteleyn, Physica 57, 536 (1972)). The square lattice is
self-dual under `p <-> 1-p` with the cluster weight `q`, and the
self-dual (critical) bond probability is

```text
p_c(q) = sqrt(q) / (1 + sqrt(q)).
```

This is the critical curve of the `q`-state Potts model: `p_c(1) = 1/2`
(percolation), `p_c(2) = sqrt(2)/(1+sqrt(2)) = 0.585786...` (Ising),
`p_c(4) = 2/3` (the marginal case), and `p_c(q) > 2/3` for `q > 4`.
The transition is continuous for `1 <= q <= 4` and first order for
`q > 4` (Baxter; Wu).

Primary references (pre-arXiv, no PDF archived):

- P. W. Kasteleyn and C. M. Fortuin, J. Phys. Soc. Japan 26 (Suppl.),
  11 (1969); C. M. Fortuin and P. W. Kasteleyn, Physica 57, 536 (1972).
- F. Y. Wu, "The Potts model", Rev. Mod. Phys. 54, 235 (1982),
  doi:10.1103/RevModPhys.54.235.
- R. J. Baxter, *Exactly Solved Models in Statistical Mechanics*,
  Academic Press (1982), Ch. 12.

The Q3 module uses this curve as an anchor input and never fits a
critical point. The finite-size phase-gap scan there is explicitly
exploratory.
