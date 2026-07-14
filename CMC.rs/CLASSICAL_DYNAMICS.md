# Classical dynamics and rejection-free algorithms

Stage 6 separates equilibrium sweep bookkeeping from physical/event time.  The
new kernels share CMC.rs states and cache audits, but they do not pretend that a
Monte Carlo sweep is a physical clock.

## Clock contract

`carlo_rs::Context` records four distinct clocks:

- scheduler sweeps;
- attempted elementary transitions;
- accepted/executed transitions;
- continuous event time.

`SimulationClock` exposes these values without converting one unit into another.
Serde/JSON context checkpoints preserve all four values and default the three
new clocks to zero when loading an older checkpoint.  The legacy HDF5 backend
still requires its separate API migration before these fields can be added to
that format.

## Kawasaki exchange

`KawasakiCore` selects a physical graph edge and proposes to exchange unlike
Ising spins.  The two-site move uses the existing transactional batch-energy
path and canonical log-domain Metropolis-Hastings acceptance.  Signed
magnetization is conserved exactly.  One sweep attempts either `N` exchanges or
an explicitly configured number.

This is stochastic conserved-order-parameter dynamics.  Its sweep count is an
algorithmic clock, not automatically a calibrated physical time.

## Direct Gillespie driver

`RejectionFreeModel` describes a finite event catalog:

```text
event_count -> event_rate -> prepare_event -> commit_event
```

`GillespieKernel` recomputes the complete rate catalog, samples one event with
probability `r_i / sum(r)`, and advances time by an exponential waiting time of
rate `sum(r)`.  `advance_by(dt)` stops exactly on the requested observation-time
boundary.  A zero-rate absorbing state advances the observation clock without
executing an event.

`KineticIsingModel` is the reference backend.  It supplies single-spin-flip
Glauber or continuous-time Metropolis rates and validates the detailed-balance
ratio against the canonical Ising energy change.

## BKL / n-fold way

`BklIsingKernel` caches all Ising flip rates in a Fenwick tree.  Event selection
and rate updates are `O(log N)` plus the degree of the flipped site.  On uniform
lattices this plays the same rejection-free role as traditional rate-class
n-fold-way tables; unlike fixed classes it also supports arbitrary weighted
CSR graphs and parallel edges.

The versioned snapshot contains:

- spins and cached physical energy;
- event time and event count;
- per-site rates;
- the exact Fenwick tree and values;
- model/rate-law identity and audit interval.

Storing the tree itself preserves the exact future trajectory after JSON
checkpoint restoration instead of merely reconstructing a numerically
close summation order.

`KineticIsingBklMC` samples at fixed event-time windows and reports energy,
magnetization, events per window, total rate and event time through Carlo.rs.

## Hard-sphere event-chain Monte Carlo

`HardSphereEventChain<D>` implements straight, rejection-free lifted chains for
identical hard spheres in an orthorhombic periodic cell.  A chain chooses an
active particle, Cartesian axis and direction.  Motion continues to the nearest
contact, transfers the lifting variable to the collided particle and proceeds
until the configured chain length is exhausted.

The first implementation deliberately uses an exact `O(N)` collision search.
It is a correctness backend and benchmark baseline; a cell-list accelerated
collision structure can be added later without changing the public lifting
contract.  Diameter, box, overlap and snapshot invariants are validated.

`HardSphereEventChainMC<D>` records chains as attempted/executed transitions and
reports packing fraction, collisions and cumulative lifted distance.  Lifted
distance is not conflated with kinetic event time.

## Cache audit

The existing `cache-audit` policy extends to Stage 6:

- Kawasaki checks cached lattice energy;
- BKL checks physical energy, every cached rate and the Fenwick sums;
- event-chain checks finite wrapped positions and absence of hard-sphere overlap.

Release builds keep automatic audits disabled unless the feature or an explicit
interval is enabled.

## Benchmarks and validation

`benches/dynamics_bench.rs` reports:

- Kawasaki attempted exchanges/s;
- direct Gillespie events/s;
- Fenwick BKL events/s;
- event-chain lifted distance/s;
- BKL checkpoint serialization;
- BKL integrated autocorrelation time, ESS and ESS/s.

The Stage 6 tests cover rate-weighted event selection, exponential waiting
times, exact fixed-time windows, Kawasaki conservation, BKL cache/checkpoint
identity, a four-site exact Ising energy comparison, analytic hard-sphere
collisions, periodic wrapping and scheduler clock separation.

## Deliberate boundaries

Stage 6 does not claim calibrated real-time dynamics for arbitrary Metropolis
sweeps.  It also does not yet include reaction networks with dynamically sized
event catalogs, event-chain cell lists, anisotropic particles or irreversible
non-equilibrium rates.  Those additions should implement explicit model/rate
contracts rather than overload canonical sampling APIs.
