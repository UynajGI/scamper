# Scuttle Git Hooks

Project-local git hooks enforcing formatting, lints, tests, and commit-message
convention. Lives under `.githooks/` (set via `core.hooksPath`, so nothing is
installed into `.git/hooks/`).

## Setup

```bash
git config core.hooksPath .githooks
```

Verify:

```bash
git config --get core.hooksPath   # → .githooks
```

Or via just:

```bash
just hooks          # sets core.hooksPath
just hooks --status # shows current config
```

## Hooks

| Hook         | When          | What it does                                                              |
| ------------ | ------------- | ------------------------------------------------------------------------- |
| `pre-commit` | `git commit`  | `cargo fmt --check` + `cargo clippy -D warnings` on affected crates only. |
| `commit-msg` | `git commit`  | Enforces [Conventional Commits](https://www.conventionalcommits.org).     |
| `pre-push`   | `git push`    | `cargo test` on affected crates (pushed commit range).                    |

### Design rationale

- **pre-commit is fast.** Only fmt-check + clippy on the crates you actually
  touched — seconds, not minutes. `cargo fmt --check` is workspace-wide (cargo
  has no per-package fmt filter), but it's near-instant.
- **pre-push runs tests.** Full `cargo test` runs once on push, so local
  history rewriting and quick WIP commits stay snappy.
- **Affected-crate scoping.** Staged files under `Carlo.rs/`, `QMC.rs/`,
  `CMC.rs/` map to crates `carlo-rs`, `qmc-rs`, `cmc-rs` and get passed as
  `-p` flags. Files outside those dirs fall back to `--workspace`.
- **Shared helpers** live in `_common.sh` (sourced by the executable hooks).

## Skipping

Hooks honor the standard `SKIP` env var (mirrors `pre-commit(1)`):

```bash
SKIP=clippy  git commit        # skip clippy only
SKIP=test    git commit        # (commit-msg has no test; harmless)
SKIP=all     git commit        # skip everything
SKIP=test    git push          # skip pre-push tests
```

Emergency overrides:

```bash
HOOK_SKIP=1        git commit   # disable pre-commit + commit-msg
HOOK_SKIP=1        git push     # disable pre-push
PUSH_HOOK_SKIP=1   git push     # disable pre-push only
```

`commit-msg` honors `SKIP=msg`.

## Conventional Commits

Enforced format:

```
<type>(<scope>)!: <subject>

<body>
```

- **type** (required): `feat` `fix` `docs` `style` `refactor` `perf` `test`
  `build` `ci` `chore` `revert`
- **scope** (optional): `Carlo` `QMC` `CMC` `carlo-rs` `qmc-rs` `cmc-rs`
  `spec` `task` `docs` `deps` `release` (unknown scopes warn, don't block)
- **`!`** marks a breaking change (e.g. `refactor(qmc-rs)!: …`)
- **subject**: imperative mood, ≤ 72 chars, no trailing period
- **body** (optional): separated by one blank line, wrap at 100 chars
  (warnings only; trailers like `Co-Authored-By:` are exempt)

Examples:

```
feat(CMC): add Wolff cluster algorithm
fix(CMC): correct beta sign in metropolis acceptance
docs: update CLAUDE.md architecture table
refactor(qmc-rs)!: replace SSE with worldline QMC
chore: bump version to 0.1.0-dev1
```

Pass-through (never blocked): `Merge …`, `Revert …`, `squash!`/`fixup!`/`amend!`.

## Troubleshooting

- **`cargo fmt --check` fails** → run `cargo fmt`, then `git add -u`.
- **clippy fails** → fix the warning, or `SKIP=clippy git commit` to defer.
- **commit-msg rejects a valid-looking message** → check for trailing period,
  > 72-char subject, missing `:` after type/scope, or an unknown shape.
- **pre-push tests the whole workspace** → happens when pushed files are
  outside the three crate dirs, or on the first push of a new branch.

## Files

```
.githooks/
├── _common.sh     # shared helpers (sourced, not executed)
├── pre-commit     # fmt-check + clippy
├── commit-msg     # Conventional Commits
├── pre-push       # cargo test
└── README.md      # this file
```
