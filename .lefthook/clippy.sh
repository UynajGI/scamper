#!/usr/bin/env bash
# Run cargo clippy scoped to the affected crates. File paths come in as args
# (lefthook expands {{staged_files}}). Empty mapping → --workspace fallback.
set -euo pipefail

PKG_FLAGS="$(./.lefthook/affected-crates.sh "$@")"
if [ -z "$PKG_FLAGS" ]; then
    PKG_FLAGS="--workspace"
    SCOPE="workspace"
else
    SCOPE="(crates: $(echo "$PKG_FLAGS" | sed 's/-p //g'))"
fi

printf '\033[1;36m▸ cargo clippy --all-targets  %s\033[0m\n' "$SCOPE"
cargo clippy --all-targets $PKG_FLAGS -- -D warnings
