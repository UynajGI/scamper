# Performance and statistical-efficiency benchmarks

The cross-cutting benchmark suite is registered as a Criterion target:

```bash
cargo bench -p cmc-rs --bench performance_bench
```

It covers:

- Ising Metropolis attempted updates per second;
- Wolff cluster sites per second;
- Swendsen-Wang physical edges per second;
- Lennard-Jones trial translations per second;
- cell-list neighbor-query cost;
- batch spin-move delta-energy cost;
- Wang-Landau JSON checkpoint serialization cost.

At benchmark startup, fixed-seed pilot chains print records such as:

```text
STAT_EFF ising_metropolis tau_int=... ess=... ess_per_second=...
```

The same report is produced for Metropolis, Wolff, Swendsen-Wang and Lennard-Jones translation. `tau_int` uses the convention that independent samples have integrated autocorrelation time 0.5, so:

```text
ESS = N / (2 * tau_int)
ESS/s = ESS / elapsed_seconds
```

Throughput and statistical efficiency answer different questions. A change should not be accepted solely because attempted updates per second increased; compare both Criterion timing and ESS/s on the same hardware, model size, temperature and seed policy.

To compile every benchmark without executing the timing loops:

```bash
cargo bench -p cmc-rs --bench performance_bench --no-run
```

## Stage 5 worm benchmark

The persistent classical worm path has a separate Criterion target:

```bash
cargo bench -p cmc-rs --bench worm_bench
```

It reports local extended-space transitions, full sweeps, endpoint-pair tracking overhead, JSON snapshot serialization and a fixed-seed physical-sector occupied-edge `tau_int` / ESS / ESS-per-second pilot. The endpoint benchmark is separate because dense correlation tracking is optional and should not be hidden in the core transition number.

## Stage 6 classical-dynamics benchmark

Dynamic and rejection-free kernels have a separate Criterion target:

```bash
cargo bench -p cmc-rs --bench dynamics_bench
```

It reports Kawasaki attempted exchanges, direct Gillespie events, Fenwick BKL
events, hard-sphere event-chain lifted distance and BKL JSON serialization.
A fixed-event-time BKL pilot also prints `tau_int`, ESS and ESS/s.  Event rate,
attempt throughput and statistical efficiency are intentionally kept separate:
a rejection-free method can execute fewer but more statistically useful state
changes, while event-chain distance is not a physical-time unit.
