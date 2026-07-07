#!/usr/bin/env bash
# cargo check scoped to the affected crates — a fast type-check gate that runs
# between fmt and clippy. Catches obvious typos / type errors in seconds,
# before the slower clippy pass.
#
# Invoked by lefthook from the pre-commit group; file paths come in as args
# (lefthook expands {{staged_files}}). Empty mapping → --workspace fallback.
set -euo pipefail

PKG_FLAGS="$(./.lefthook/affected-crates.sh "$@")"
if [ -z "$PKG_FLAGS" ]; then
    PKG_FLAGS="--workspace"
    SCOPE="workspace"
else
    SCOPE="(crates: $(echo "$PKG_FLAGS" | sed 's/-p //g'))"
fi

printf '\033[1;36m▸ cargo check  %s\033[0m\n' "$SCOPE"
cargo check $PKG_FLAGS
