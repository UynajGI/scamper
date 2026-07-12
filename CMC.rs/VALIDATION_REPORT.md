# Validation report

Date: 2026-07-12

## Completed checks

- Parsed all 94 Rust source files in the final Scuttle workspace with the Tree-sitter Rust grammar: no syntax errors or missing syntax nodes.
- Parsed every `Cargo.toml` with Python `tomllib`.
- Parsed all 13 CMC.rs Rust source/test files independently.
- Checked the new periodic honeycomb construction at `(Lx, Ly) = (2,3), (4,4), (6,5)`: symmetric adjacency, degree 3 at every site, and `3N` incidences.
- Scanned CMC.rs for removed cluster APIs (`fk_bond_probability`, `random_cluster_spin`, scalar/vector branch shortcuts) and the old physical-energy `sum / 2` convention: no residual uses.
- Added internal tests for physical-edge counting, weighted/parallel edges, arbitrary O(N), endpoint-dependent cluster probability, cached-energy consistency, independent SW cluster assignments, Carlo.rs scheduler integration, strict parameter parsing, snapshot topology mismatch, and packed-replica PT energy.
- The delivered ZIP is tested with `unzip -t` after creation.

## Toolchain limitation

This execution environment did not contain `rustc`, `cargo`, or `rustfmt`. An attempt to install the stable Rust toolchain failed because outbound DNS/network access to `static.rust-lang.org` is disabled. Therefore `cargo fmt`, `cargo check`, `cargo clippy`, and `cargo test` were **not** executed here.

The first local verification command should be:

```bash
cargo fmt --all --check
cargo check --workspace --all-targets
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```
