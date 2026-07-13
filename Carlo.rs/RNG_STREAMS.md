# Deterministic RNG stream contract

`RngStreamKey` derives domain-separated seeds from logical simulation identity:

```text
base seed
+ task ID
+ run ID
+ chain ID
+ replica ID
+ logical thread ID
+ lifecycle phase
+ substream ID
```

Use builder methods and then construct the concrete generator:

```rust,ignore
use carlo_rs::{RngPhase, RngStreamKey};
use rand_xoshiro::Xoshiro256PlusPlus;

let rng: Xoshiro256PlusPlus = RngStreamKey::new(master_seed)
    .with_task(task_id)
    .with_chain(chain_id)
    .with_replica(replica_id)
    .with_phase(RngPhase::Measurement)
    .seeded();
```

Rules:

1. IDs are logical and stable. Do not use an opportunistic Rayon worker index when task scheduling can change.
2. A module must not call `thread_rng()` or invent arithmetic offsets such as `base_seed + task_id * 10000`.
3. Independent chains and replicas derive their streams from their IDs, not by consuming another stream with `next_u64()`.
4. Checkpoints store the concrete RNG state. Resumption continues that exact state rather than deriving a new seed.
5. `RngPhase` separates initialization, thermalization, measurement, exchange, checkpoint and backend-task domains when distinct streams are intentionally required.

Carlo.rs schedulers, Rayon tasks, MPI run identities, parallel tempering, MCMC multi-chain execution and MCMC replica exchange now use this contract. The workspace enables serde_json's `float_roundtrip` feature so JSON checkpoints preserve all serialized floating-point state and traces exactly as well as the integer RNG state.
