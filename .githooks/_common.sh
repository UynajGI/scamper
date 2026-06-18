#!/usr/bin/env bash
# Shared helpers for Scuttle git hooks. Sourced, not executed.
#
# Provides:
#   repo_root        — echo absolute path to the repo root
#   staged_files     — echo staged .rs files (added/copied/modified)
#   affected_crates  — echo deduped `-p <crate>` flags for staged Rust files
#   hook_skip NAME   — returns 0 (true) if NAME is in $SKIP (e.g. SKIP=clippy)
#
# Scuttle crate layout:
#   Carlo.rs/ -> crate `carlo-rs`
#   QMC.rs/   -> crate `qmc-rs`
#   CMC.rs/   -> crate `cmc-rs`
#
# Hook behavior is opt-out via the standard `SKIP` env var, mirroring
# pre-commit(1): `SKIP=clippy git commit` skips the clippy check.

# Resolve repo root defensively: hooks run with cwd = repo root by default,
# but GIT_DIR / worktree configs can change that.
repo_root() {
    git rev-parse --show-toplevel
}

# Staged .rs files (added/copied/modified — not deleted).
staged_files() {
    git diff --cached --name-only --diff-filter=ACM -- '*.rs'
}

# Deduped `-p <crate>` flags derived from staged paths.
# Empty output (no Rust changes in known crates) → caller should early-exit.
affected_crates() {
    local crates=""
    local f
    while IFS= read -r f; do
        [ -z "$f" ] && continue
        case "$f" in
            Carlo.rs/*) crates="${crates}carlo-rs " ;;
            QMC.rs/*)   crates="${crates}qmc-rs " ;;
            CMC.rs/*)   crates="${crates}cmc-rs " ;;
        esac
    done <<EOF
$(staged_files)
EOF
    # dedupe while preserving order
    local seen="" out="" c
    for c in $crates; do
        case " $seen " in
            *" $c "*) ;;
            *) seen="$seen $c"; out="${out}-p ${c} " ;;
        esac
    done
    echo "$out"
}

# True (0) if the given hook/check name appears in $SKIP (whitespace/comma sep).
hook_skip() {
    local name="$1"
    local skip="${SKIP-}"
    [ -z "$skip" ] && return 1
    case " $skip " in
        *[[:space:]]${name}[[:space:]]*) return 0 ;;
        *\",${name},\"*) return 0 ;;
        *) return 1 ;;
    esac
}

# Print a labelled section header so hook output is greppable.
hook_header() {
    printf '\n\033[1;36m▸ %s\033[0m\n' "$1"
}

hook_fail() {
    printf '\033[1;31m✗ %s\033[0m\n' "$1" >&2
    exit 1
}

hook_pass() {
    printf '\033[1;32m✓ %s\033[0m\n' "$1"
}
