# VMC-Family Paper Corpus (arXiv) — Manifest

Generated 2026-08-20. All PDFs fetched from `https://arxiv.org/pdf/<id>`; markdown produced with
`mineru-open-api v0.5.9 flash-extract <pdf> -o md/<slug>/ --language en` (token-free flash mode,
serialized with 10–15 s spacing; no HTTP 429 rate-limit hits encountered).

Layout: `pdf/<slug>.pdf` = source; `md/<slug>/<slug>.md` = converted markdown.
Page counts are estimates parsed from PDF internals (`/Type /Page` objects).

| file | arXiv id | title | authors | year | topic tags | pages | status |
|---|---|---|---|---|---|---|---|
| carleo-troyer-2017-nqs | 1606.02318 | Solving the Quantum Many-Body Problem with Artificial Neural Networks (Science 355, 602) | G. Carleo, M. Troyer | 2017 | nqs, vmc | ~10 | md-ok |
| sorella-1998-sr | cond-mat/9803107 | Green Function Monte Carlo with Stochastic Reconfiguration (PRL 80, 4558) | S. Sorella | 1998 | vmc, sr | ~13 | md-ok |
| toulouse-assaraf-umrigar-2015-vmc-dmc-intro | 1508.02989 | Introduction to the variational and diffusion Monte Carlo methods (Adv. Quantum Chem. 73) | J. Toulouse, R. Assaraf, C. J. Umrigar | 2015 | dmc, vmc, review | ~26 | pdf-only |
| medvidovic-carleo-2024-nqs-notes | 2402.11014 | Neural-network quantum states for many-body physics (EPJ Plus 139, 638) | M. Medvidović, G. Carleo | 2024 | nqs, review, vmc | ~26 | pdf-only |
| pfau-2020-ferminet | 1909.02487 | Ab initio solution of the many-electron Schrödinger equation with deep neural networks (Phys. Rev. Research 2, 033429) | D. Pfau, J. S. Spencer, A. G. D. G. Matthews, W. M. C. Foulkes, et al. | 2020 | nqs, vmc | ~21 | pdf-only |
| baroni-moroni-1998-reptation | cond-mat/9808213 | Reptation Quantum Monte Carlo (PRL 82, 2130) | S. Baroni, S. Moroni | 1998 | reptation | ~29 | pdf-only |
| carleo-2010-reptation-lattice | 1003.3696 | Reptation quantum Monte Carlo for lattice Hamiltonians with a directed-update scheme | G. Carleo, F. Becca, S. Moroni, S. Baroni | 2010 | reptation | ~11 | md-ok |
| carleo-2016-tvmc | 1612.06392 | Unitary dynamics of strongly-interacting Bose gases with time-dependent variational Monte Carlo in continuous space (PRX 7, 031026) | G. Carleo, L. Cevolani, L. Sanchez-Palencia, M. Holzmann | 2017 | tvmc | ~11 | md-ok |
| kwon-1998-backflow-eg | cond-mat/9803092 | Effects of Backflow Correlation in the Three-Dimensional Electron Gas: Quantum Monte Carlo Study (PRB 58, 6800) | Y. Kwon, D. M. Ceperley, R. M. Martin | 1998 | backflow, vmc, electron-gas | ~19 | md-ok |
| holzmann-2019-backflow | 1910.07167 | Orbital-dependent backflow wave functions for real-space quantum Monte Carlo | M. Holzmann, S. Moroni | 2019 | backflow, vmc | ~5 | md-ok |
| calcavecchia-2016-shadow | 1604.05804 | Metal-Insulator Transition of Solid Hydrogen by the Antisymmetric Shadow Wave Function | F. Calcavecchia, T. D. Kühne | 2016 | shadow, vmc | ~13 | md-ok |

## Coverage

The corpus spans the variational-QMC family end to end: stochastic-reconfiguration optimization
(Sorella 1998), real-space VMC/DMC methodology including optimization and fixed-node projector
technique (Toulouse-Assaraf-Umrigar 2015), reptation QMC (the original Baroni-Moroni preprint plus
the lattice directed-update variant), t-VMC real-time dynamics (Carleo-Cevolani-Sanchez-Palencia-
Holzmann), backflow-correlated Jastrow-Slater wave functions for the electron gas (Kwon-Ceperley-
Martin) and for real-space QMC (Holzmann-Moroni), shadow wave functions (Calcavecchia-Kühne), and
the neural-quantum-state line from the Carleo-Troyer foundation through FermiNet to the 2024
Medvidović-Carleo lecture notes. 7/11 items are converted to markdown; 4 long items (21–29 pages)
exceed the 20-page flash-extract cap and remain `pdf-only` — converting them requires a MinerU
token (`mineru-open-api extract`) or `--pages` chunked extraction.

## Classics NOT on arXiv (source separately)

- W. M. C. Foulkes, L. Mitas, R. J. Needs, G. Rajagopal, "Quantum Monte Carlo simulations of solids", Rev. Mod. Phys. 73, 33 (2001) — no arXiv deposit exists; the id suggested in the brief (cond-mat/0008185) belongs to an unrelated EPR paper. Free PDF mirror: aquila.infn.it (CPierleoni QMC bibliography).
- W. L. McMillan, Phys. Rev. 138, A442 (1965) — original Jastrow VMC for He (pre-arXiv).
- J. B. Anderson, J. Chem. Phys. 63, 1499 (1975) — original fixed-node DMC (pre-arXiv).
- D. M. Ceperley, B. J. Alder, Phys. Rev. Lett. 45, 566 (1980) — fixed-node DMC of electron gas (pre-arXiv).
- M. H. Kalos et al., GFMC series (Phys. Rev. 128, 179 (1962); Kalos-Whitlock book) — pre-arXiv.
- C. J. Umrigar, M. P. Nightingale, K. J. Runge, J. Chem. Phys. 99, 2865 (1993) — SR-precursor variance/energy minimization (not deposited).
- S. Baroni, S. Petrosyan, "A guide to reptation quantum Monte Carlo" (book chapter) — no arXiv version found; nearest arXiv stand-ins included here (cond-mat/9808213, 1003.3696).
- F. Becca, S. Sorella, "Quantum Monte Carlo Approaches for Correlated Systems" (Cambridge Univ. Press, 2017) — book, not on arXiv; the Sorella-Casula-Tavana "Optimized trial wave functions"-type review was also not locatable on arXiv.
