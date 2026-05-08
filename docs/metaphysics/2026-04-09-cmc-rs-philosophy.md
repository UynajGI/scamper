# CMC.rs Philosophy

> This document defines core principles and design criteria, serving as the judgment basis for all architectural decisions.
>
> Established through Socratic dialogue on 2026-04-09.

---

## I. Ontology: What is CMC.rs?

### Core Definition

CMC.rs is a **Classical Monte Carlo algorithm toolbox** — a component library between the Carlo.rs framework layer and concrete model implementations. It is the **peer sibling** of QMC.rs: symmetric architecture, independent algorithm domains.

### Position in the Stack

```
┌─────────────────────────────────────────────────────┐
│  User's Model Implementation                        │
│  (implements ModelMC trait, defines Hamiltonian)    │
├─────────────────────────────────────────────────────┤
│  CMC.rs — Algorithm Toolbox                         │
│  ├─ MetropolisCore<MC>  (single-spin flip)          │
│  ├─ WolffCore<MC>       (cluster growth)            │
│  └─ SWCore<MC>          (Swendsen-Wang)             │
├─────────────────────────────────────────────────────┤
│  Carlo.rs — Core Framework                          │
│  ├─ MonteCarlo trait (base)                        │
│  ├─ Backend trait (parallel abstraction)           │
│  ├─ Scheduler (thermalization → measurement)       │
│  └─ Measurements + Results                         │
└─────────────────────────────────────────────────────┘
```

### Relationship with QMC.rs

CMC.rs and QMC.rs are **equal siblings** — same architectural patterns, same granularity, same trait hierarchy approach. They share no code dependencies; each defines its own lattice topology types independently.

### What Users Implement

Users implement **Hamiltonian decomposition**: defining energy changes, coupling constants, and spin configurations via `ModelMC` trait. The engine automatically handles update logic, acceptance probabilities, and measurement accumulation.

**User responsibility**: Define the model (energy function, lattice, parameters).
**Framework responsibility**: Handle update algorithms, RNG management, measurement accumulation, checkpointing.

### Granularity: Mid-level

Same as QMC.rs — users implement the model interface, engines handle algorithm internals. Users do not touch algorithm internals unless they choose Level 3 (custom engine).

### Domain Coverage (Long-term Vision)

| Phase | Domain | Algorithms |
|-------|--------|------------|
| Phase 1 | Lattice Classical MC | Metropolis, Wolff, Swendsen-Wang (Ising, Potts, XY, Heisenberg) |
| Phase 2+ | Continuous Space MC | Molecular dynamics, path integral classical limit |
| Phase 3+ | MCMC / Bayesian Inference | Separate package (not CMC.rs) |

### Trait Hierarchy

```rust
// Carlo.rs (base)
trait MonteCarlo { fn sweep(&mut self, ctx: &mut Context); }

// CMC.rs (domain layer)
trait LatticeMC { fn lattice(&self) -> &Lattice; }

// CMC.rs (method layer)
trait ModelMC: LatticeMC {
    fn spin_dim(&self) -> usize;
    fn coupling(&self) -> f64;
    fn beta(&self) -> f64;
    fn local_energy_change(&self, site: usize, old: f64, new: f64) -> f64;
    fn total_energy(&self) -> f64;
    fn spins(&self) -> &[f64];
    fn spins_mut(&mut self) -> &mut [f64];
}
```

---

## II. Teleology: What is Success?

### Priority Hierarchy

**C > B > A > D** (same as QMC.rs)

| Priority | Criterion | Description |
|----------|-----------|-------------|
| **C (Foundation)** | Correctness | Results match exact solutions and Carlo.jl reference within statistical error (≤ 3σ). Without correctness, nothing else matters. |
| **B (Necessity)** | Performance | Must achieve competitive sweep rates (not slower than reference implementations). Without performance, users won't adopt. |
| **A (Goal)** | Efficiency | Users implement new models in days, not weeks. This is the ultimate user value. |
| **D (Optional)** | Education | Code readability as teaching material — nice-to-have, not required. |

### Success Metrics

1. **Correctness**: 2D Ising results match Onsager exact solution within ≤ 3σ. Models without exact solutions match Carlo.jl reference within ≤ 3σ.
2. **Performance**: Sweep rate ≥ 80% of optimized C++/Fortran implementations for standard models.
3. **Efficiency**: New model implementation time ≤ 3 days for users familiar with the trait system.

### Validation Strategy

- **Onsager exact solution** for 2D Ising (primary benchmark)
- **Carlo.jl reference** comparison for all models
- **Both required** — exact solution where available, Carlo.jl where not

---

## III. Methodology: How to Achieve It?

### Core Design Principle: Progressive Complexity (3 Layers)

Simplified from QMC.rs's 4 layers — more pragmatic, fewer abstraction levels.

```
Level 1: Use built-in IsingModel + MetropolisCore → run simulation immediately
Level 2: Implement custom ModelMC → use built-in engines with your physics
Level 3: Replace engine entirely → custom algorithm for special cases
```

Users can start at Level 1 and progressively customize. No user should need Level 3 for standard use cases.

### Cross-Package Symmetry

CMC.rs and QMC.rs share the same architectural patterns:
- Wrapper pattern: `*Core<MC>` implements `MonteCarlo`
- Trait hierarchy: base → domain → method
- Lattice topology: adjacency list representation (independently defined in each package)
- Model placement: inside the package, not as a separate crate

### Model Placement: Inside CMC.rs

Concrete models live inside CMC.rs:
- `cmc::models::IsingModel`
- `cmc::models::PottsModel` (future)
- `cmc::models::XYModel` (future)

**Rationale**: Same as QMC.rs — models demonstrate trait usage, provide test coverage, reduce user friction.

### Implementation Order

Phase 1: Lattice infrastructure (LatticeMC, Lattice, BondType, builders)
Phase 2: ModelMC trait + IsingModel
Phase 3: MetropolisCore + WolffCore
Phase 4: SWCore + validation against Onsager/Carlo.jl
Phase 5: XY/Heisenberg models

---

## IV. Boundaries: What is CMC.rs NOT?

### Explicit Exclusions

| Exclusion | Reason | Where it belongs |
|-----------|--------|------------------|
| **MCMC / Bayesian Inference** | Different problem class (statistical inference vs. statistical physics) | Future separate package (e.g., MCMC.rs) |
| **Visualization/Analysis** | Post-processing is domain-specific | Python scripts, Jupyter notebooks |
| **Parallel Scheduling** | Framework-level responsibility | Carlo.rs Backend trait |
| **Continuous Space MC** | Different abstraction (not lattice-based) | Future package extension |

### What it CAN Include

| Included | Reason |
|----------|--------|
| **Concrete models** | Demonstrate usage, provide test coverage |
| **Algorithm components** | Core value proposition |
| **Measurement helpers** | Classical MC-specific observables |

---

## V. Decision Criteria

When facing design choices, test with these questions:

### Ontology Check
1. Is this a framework-level feature (belongs to Carlo.rs) or CMC-specific?
2. What granularity level is appropriate? (User responsibility vs. engine responsibility)
3. Does this respect the sibling relationship with QMC.rs?

### Teleology Check
4. Does this design preserve correctness guarantees?
5. Does this design enable competitive performance?
6. Does this design reduce user implementation effort?

### Methodology Check
7. Can users start simple (Level 1) and progressively customize?
8. Is the trait hierarchy clean and minimal?

### Boundary Check
9. Is this feature within CMC.rs scope (not Carlo.rs, not MCMC.rs, not post-processing)?

---

## VI. Key Architectural Decisions Summary

| Decision Point | Choice | Rationale |
|----------------|--------|-----------|
| **Position** | Algorithm toolbox (component library) | Between framework and models, maximum flexibility |
| **Relationship with QMC** | Equal siblings, independent domains | Symmetric architecture, no code sharing |
| **Granularity** | Mid-level (model interface → engine handles updates) | Balance ease-of-use with customization |
| **Lattice types** | Independently defined (no cross-package dependency) | Enforces dependency direction: CMC→Carlo, QMC→Carlo |
| **Progressive complexity** | 3 layers (simplified from QMC's 4) | More pragmatic for classical MC |
| **Validation** | Onsager exact + Carlo.jl reference (both required) | Gold standard + practical cross-check |
| **Success priority** | Correctness → Performance → Efficiency | Foundation before value |
| **Model placement** | Inside CMC.rs domain modules | Reduce friction, provide examples |

---

*This philosophy document was established through Socratic dialogue. It serves as the design criterion for all CMC.rs architectural decisions.*
