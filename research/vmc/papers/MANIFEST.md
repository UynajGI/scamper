# VMC-Family Paper Corpus (arXiv) — Manifest

Generated 2026-08-20. All PDFs fetched from `https://arxiv.org/pdf/<id>`; markdown produced with
`mineru-open-api v0.5.9 flash-extract <pdf> -o md/<slug>/ --language en` (token-free flash mode,
serialized with 10–15 s spacing; no HTTP 429 rate-limit hits encountered).
Updated 2026-08-21: three non-arXiv classics added from legitimate open mirrors (see table notes);
no paywalled or shadow-library sources were used.

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
| ceperley-alder-1980-electron-gas | — eScholarship (LBL preprint) | Ground State of the Electron Gas by a Stochastic Method (PRL 45, 566) | D. M. Ceperley, B. J. Alder | 1980 | qmc, electron-gas, dmc, green-oa | ~15 | md-ok |
| umrigar-1993-dmc-small-timestep | — URI DigitalCommons | A diffusion Monte Carlo algorithm with very small time-step errors (J. Chem. Phys. 99, 2865) | C. J. Umrigar, M. P. Nightingale, K. J. Runge | 1993 | dmc, time-step, population-control | ~29 | pdf-only |
| foulkes-2001-rmp-qmc-solids | — INFN L'Aquila mirror | Quantum Monte Carlo simulations of solids (Rev. Mod. Phys. 73, 33) | W. M. C. Foulkes, L. Mitas, R. J. Needs, G. Rajagopal | 2001 | review, vmc, dmc, fixed-node | ~51 | pdf-only |

Non-arXiv additions (2026-08-21), fetched via legitimate open-access routes only:

- `ceperley-alder-1980-electron-gas`: green-OA author preprint (LBL-11180, May 1980) from the UC
  eScholarship repository — https://escholarship.org/content/qt2d7023jm/qt2d7023jm_noSplash_41dd75deecca357d03aed90236b6dcd2.pdf
  (scanned 1980 typescript; OCR has era-typical glitches, physics-term density verified).
- `umrigar-1993-dmc-small-timestep`: institutional-repository deposit (Runge/URI, incl. required AIP
  citation attribution) — https://digitalcommons.uri.edu/cgi/viewcontent.cgi?article=1279&context=phys_facpubs
- `foulkes-2001-rmp-qmc-solids`: INFN L'Aquila open mirror (CPierleoni QMC bibliography) —
  https://www.aquila.infn.it/cpierleo/MCMC/QMC/BIBLIO/FoulkesMitasRPM2001.pdf
  The same mirror also holds Ceperley's 1995 helium-path-integrals RMP (`Ceperley95RMP67_279.pdf`, 15 MB,
  >20 pp flash cap) if ever wanted.
- Both `umrigar-1993` (29 pp) and `foulkes-2001` (51 pp) exceed the 20-page flash-extract cap and are
  `pdf-only`; converting them requires a MinerU token or chunked `--pages` extraction.

## Coverage

The corpus spans the variational-QMC family end to end: stochastic-reconfiguration optimization
(Sorella 1998), real-space VMC/DMC methodology including optimization and fixed-node projector
technique (Toulouse-Assaraf-Umrigar 2015), reptation QMC (the original Baroni-Moroni preprint plus
the lattice directed-update variant), t-VMC real-time dynamics (Carleo-Cevolani-Sanchez-Palencia-
Holzmann), backflow-correlated Jastrow-Slater wave functions for the electron gas (Kwon-Ceperley-
Martin) and for real-space QMC (Holzmann-Moroni), shadow wave functions (Calcavecchia-Kühne), and
the neural-quantum-state line from the Carleo-Troyer foundation through FermiNet to the 2024
Medvidović-Carleo lecture notes — plus, since 2026-08-21, the founding stochastic-method electron-gas
calculation (Ceperley-Alder 1980), the small-time-step DMC/population-control algorithm paper
(Umrigar-Nightingale-Runge 1993), and the standard QMC-for-solids review (Foulkes et al. 2001).
8/14 items are converted to markdown; 6 long items (21–51 pages) exceed the 20-page flash-extract
cap and remain `pdf-only` — converting them requires a MinerU token (`mineru-open-api extract`)
or `--pages` chunked extraction.

## Classics NOT on arXiv (source separately)

Sourced legitimately 2026-08-21 and moved into the table above: Foulkes-Mitas-Needs-Rajagopal RMP 2001
(INFN mirror), Ceperley-Alder PRL 1980 (eScholarship preprint), Umrigar-Nightingale-Runge JCP 1993
(URI DigitalCommons). Still missing — closed-access everywhere, no legitimate free copy found:

- W. L. McMillan, Phys. Rev. 138, A442 (1965) — original Jastrow VMC for He. Tried: Unpaywall /
  Semantic Scholar / OpenAIRE all report `closed`; APS paywall; author deceased 1984 so no author
  page; NASA ADS scanned-article endpoint session-gated for non-browser fetchers; INFN mirror holds
  only the Foulkes/Ceperley RMPs; DOI-string web search surfaced no institution-hosted copy.
- J. B. Anderson, J. Chem. Phys. 63, 1499 (1975) "A random-walk simulation of the Schrödinger
  equation: H3+" — original DMC. Tried: same checks as above; AIP paywall; Penn State personal page
  unreachable; course servers host only *about* it, not the scan; ADS gated.
- M. H. Kalos et al., GFMC series (Phys. Rev. 128, 1791 (1962); Kalos-Whitlock book) — pre-arXiv.
  Tried: same closed-access checks; OSTI record is bibliography-only; INSPIRE hosts metadata only;
  ADS gated. Nearest legitimately-free stand-ins already in corpus: reptation QMC papers and the
  Toulouse-Assaraf-Umrigar 2015 introduction.
- C. J. Umrigar, M. P. Nightingale, K. J. Runge, J. Chem. Phys. 99, 2865 (1993) — obtained
  (see table); a green-OA deposit exists and should be preferred over any aggregator mirror.
- S. Baroni, S. Petrosyan, "A guide to reptation quantum Monte Carlo" (book chapter) — no arXiv version found; nearest arXiv stand-ins included here (cond-mat/9808213, 1003.3696).
- F. Becca, S. Sorella, "Quantum Monte Carlo Approaches for Correlated Systems" (Cambridge Univ. Press, 2017) — book, not on arXiv; the Sorella-Casula-Tavana "Optimized trial wave functions"-type review was also not locatable on arXiv. Tried (2026-08-21): the historically free author draft is no longer posted on any live page (people.sissa.it/~becca, people.sissa.it/~sorella incl. lecture/research pages, turborvb.sissa.it, cm.sissa.it); the SISSA IRIS record is Cloudflare-gated; web.archive.org unreachable from this host, so even the archived author-shared draft could not be retrieved. Publisher version is paid; no pirated copies used.
