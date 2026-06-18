# Lefthook helpers for Scuttle

Thin scripts invoked by `lefthook.yml`. They encapsulate logic that's awkward
to express inline in YAML (file-list → crate mapping, Conventional Commits
regex). All are plain bash, executable, and testable standalone.

## Files

| Script                 | Called from      | Purpose                                                 |
| ---------------------- | ---------------- | ------------------------------------------------------- |
| `fmt-check.sh`         | pre-commit / fmt | `cargo fmt --check` with a helpful failure hint.        |
| `clippy.sh <files...>` | pre-commit / clippy | `cargo clippy -D warnings` scoped to affected crates.   |
| `affected-crates.sh`   | (by clippy.sh)  | Map file paths → deduped `-p <crate>` flags.            |
| `commit-msg-check.sh <msg-file>` | commit-msg | Conventional Commits validation.              |
| `pre-push-test.sh`     | pre-push / test  | `cargo test` scoped to affected crates of push range.   |

## Crate mapping

```
Carlo.rs/* → -p carlo-rs
QMC.rs/*   → -p qmc-rs
CMC.rs/*   → -p cmc-rs
```

Files outside these dirs (or no Rust files staged) fall back to `--workspace`.

## Standalone testing

```bash
./.lefthook/affected-crates.sh CMC.rs/src/x.rs QMC.rs/y.rs   # → "-p cmc-rs -p qmc-rs "
./.lefthook/commit-msg-check.sh /path/to/COMMIT_EDITMSG       # exit 0/1
./.lefthook/clippy.sh CMC.rs/src/lib.rs                       # runs clippy on cmc-rs
./.lefthook/fmt-check.sh                                      # runs cargo fmt --check
./.lefthook/pre-push-test.sh                                  # runs cargo test on push range
```

## See also

- [`../lefthook.yml`](../lefthook.yml) — the hook configuration.
- [`../CLAUDE.md`](../CLAUDE.md) — top-level project notes (hooks section).
