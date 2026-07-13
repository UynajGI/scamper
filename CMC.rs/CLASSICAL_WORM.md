# Persistent classical worm infrastructure

Stage 5 adds a reusable extended-configuration-space driver and one fully
validated model: the ferromagnetic zero-field Ising high-temperature graph
representation.

## State and model split

The reusable layer consists of:

```text
WormState<Configuration, Defect>
    Physical sector: no endpoints
    Worm sector: ordered tail/head endpoints

WormModel
    opening-defect measure
    available local step
    local log-weight ratio
    transactional patch/commit
    defect/cache validation
    optional endpoint-bin mapping

WormKernel<Model>
    open / close / move-head proposal selection
    all Hastings corrections
    log-domain acceptance
    sector and transition statistics
    optional endpoint-pair histogram
    periodic cache audit
```

`WormKernel` is persistent: a sweep is a configurable number of local extended-
space transitions, and the chain may remain in the worm sector across sweep or
checkpoint boundaries. This is intentional. It makes the sampled extended
distribution explicit and gives endpoint observables a well-defined stationary
measure.

## Open, close and local-step detailed balance

The worm sector has relative fugacity `eta = exp(log_worm_fugacity)`. Opening
selects one of `N_open` defects uniformly. When head and tail coincide, closure
is proposed with probability `p_close`; otherwise a local step is mandatory.

The open and close log acceptance ratios are:

```text
open:  log(eta) + log(p_close) + log(N_open)
close: -log(eta) - log(p_close) - log(N_open)
```

For a local move `head -> head'`, the model provides the configuration log-
weight ratio and the local proposal ratio. The generic kernel additionally
includes the close-versus-step branch correction:

```text
log A = log(W_new / W_old)
      + log(q_reverse / q_forward)
      + log(p_step(new) / p_step(old))
```

Every stochastic decision uses `carlo_rs::accept_log_probability`; no
exponentiated acceptance probability is formed.

## Ising high-temperature graph model

For physical edge coupling `J_e = J * edge.weight`, define:

```text
K_e = beta * J_e
 t_e = tanh(K_e)
```

A graph configuration stores one occupation bit `n_e` per physical edge. The
physical sector contains even-degree subgraphs with reduced weight:

```text
W(graph) = product_e t_e ^ n_e
```

A worm-sector graph has odd parity at exactly the distinct head and tail
vertices. A local head move toggles one incident physical edge. On an irregular
graph, the edge-incidence proposal includes the exact degree Hastings factor:

```text
log(q_reverse / q_forward) = log(degree(old_head) / degree(new_head))
```

The configuration caches:

- occupied-edge count;
- vertex parity bits;
- reduced log graph weight.

Trial evaluation is read-only. Accepted toggles patch all three caches once.
The `cache-audit` policy recomputes and checks every invariant without repairing
it first.

### Supported graph class

- arbitrary loop-free `CsrLattice`;
- weighted and parallel physical edges;
- isolated vertices are allowed and produce a valid open-sector bounce;
- only non-negative effective couplings `J * edge.weight` are supported.

Negative graph weights would create a sign problem and are rejected. Self-loops
are also rejected in this first representation because they do not create the
two-defect parity transition used by the local driver.

## Observables

`IsingGraphWormMC` is directly constructible by Carlo.rs `Scheduler`. It records:

- `WormSector` and `PhysicalSector` indicators;
- open, close and step acceptance fractions;
- last completed worm length;
- physical-sector graph occupation, edge density and canonical energy estimator;
- current head and tail in the worm sector;
- optional dense one-hot `WormEndpointPairs` samples.

The cumulative `EndpointPairHistogram` is available from the kernel. For equal
endpoint fugacity, the ratio

```text
count(tail, head) / count(tail, tail)
```

estimates the Ising two-point correlation for that ordered endpoint pair.

The canonical physical-sector energy estimator is:

```text
E = -sum_e J_e tanh(K_e)
    -sum_e n_e J_e [1 / tanh(K_e) - tanh(K_e)]
```

Physical observables are emitted only while the chain is in the physical
sector; Carlo.rs accumulates each observable with its own sample count.

## Parameters

`IsingGraphWormMC::from_params` accepts the standard lattice parameters plus:

| Parameter | Default | Meaning |
|---|---:|---|
| `beta` | `1.0` | inverse temperature |
| `J` | `1.0` | global Ising coupling |
| `worm_updates_per_sweep` | physical edge count | local transitions per sweep |
| `worm_close_probability` | `0.25` | closure proposal probability at coincident endpoints |
| `worm_fugacity` | `1 / n_sites` | positive relative worm-sector fugacity |
| `log_worm_fugacity` | unset | log fugacity alternative; mutually exclusive with `worm_fugacity` |
| `worm_track_endpoint_pairs` | `false` | enable dense endpoint-pair samples and cumulative histogram |
| `worm_cache_audit_interval` | `0` | explicit audit cadence; zero uses build-mode policy |

## Checkpoint contract

`IsingGraphWormMC::save_snapshot` writes `cmc-rs-ising-worm-v1`, including:

- exact topology and model parameters;
- transition configuration;
- graph occupations and sector endpoints;
- transition counters, current worm length and sweep count;
- optional endpoint-pair counts.

`load_snapshot` validates constructor-owned model and transition parameters before
restoring. Carlo.rs remains responsible for concrete RNG and measurement
checkpoint state. Restoring both produces the exact future trajectory.

## Scope boundary

Stage 5 does not claim a universal worm representation. Integer-current, dimer,
loop-gas and other defect spaces can implement `WormModel`, but they are not
included in this version. Quantum directed-loop/wormhole code remains a
separate QMC.rs concern.
