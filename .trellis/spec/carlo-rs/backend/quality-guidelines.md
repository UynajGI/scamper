# Carlo.rs Quality Guidelines

## Required patterns

### Implement MonteCarlo + FromParams together

Every simulation type must implement both:

```rust
impl MonteCarlo for MyModel {
    type Rng = Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        // Access RNG via ctx.rng (pub field)
        // Call ctx.measure() during sweep if measuring inline
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        ctx.measure("Observable", value);
        // ctx.measure_array("Correlation", &array);
        // ctx.measure_complex("OrderParam", re, im);
    }
}

impl FromParams for MyModel {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        // Parse parameters, initialize state with RNG
        // This is the RIGHT place for random initialization — RNG is provided here
    }
}
```

### Use ctx.rng directly — it's a pub field

```rust
fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
    let x: f64 = ctx.rng.random();          // rand 0.10: random() not gen()
    let j = ctx.rng.random_range(0..=i);    // rand 0.10: random_range() not gen_range()
}
```

Requires `use rand::RngExt;` in scope for convenience methods.

### FromParams receives RNG for a reason

Carlo.rs deliberately passes `rng: &mut Self::Rng` to `from_params`. This is the designated place for random initialization (random spins, random lattice). Do not defer initialization to a separate method unless you also update the Scheduler to call it.

### Keep measure() lightweight

`measure()` is called after every measurement sweep. Don't compute O(N^2) quantities here. Pre-compute during sweep if possible. `ctx.measure()` is cheap (just pushes to a bin accumulator).

### Type alias for composed types

For crates that provide pre-composed MonteCarlo impls, provide type aliases:

```rust
type IsingMetropolis = ClassicalMC<IsingModel, MetropolisCore>;
type IsingWolff = ClassicalMC<IsingModel, WolffCore>;
```

## Forbidden patterns

### Don't store RNG in your model

The Context owns the RNG. Use `ctx.rng`, don't create and store a separate RNG.

### Don't call ctx.advance_sweep()

The Scheduler calls `ctx.advance_sweep()` after each sweep. Your sweep() should only update configuration.

### Don't panic in from_params

Return `Err(CarloError::InvalidConfig { field, reason })` for bad parameters. Let the Scheduler handle errors.

### Don't use gen() / gen_range() with rand 0.10

Rand 0.10 renamed these to `random()` / `random_range()`. Requires `use rand::RngExt;`.

## Testing MonteCarlo impls

### Unit tests (no Scheduler)
- Test sweep logic directly: create Context with `seed_from_u64(42)`, call sweep() in loop, verify state converged

### Integration tests (with Scheduler)
- `Scheduler.run_one::<YourType>(&params)` — end-to-end
- Verify `results.get("Observable")` returns Estimate with reasonable mean and positive stderr
- For known exact solutions (e.g., Onsager), validate within tolerance + 3σ error

## CarloError conventions

```rust
// Bad parameter
CarloError::InvalidConfig { field: "beta".into(), reason: "must be positive".into() }

// Missing required measurement
CarloError::MeasurementNotFound { name: "Energy".into() }
```
