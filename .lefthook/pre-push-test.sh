#!/usr/bin/env bash
# Pre-push test runner for Scuttle (invoked by lefthook).
#
# Runs `cargo test` scoped to the crates touched in the pushed range.
# Lefthook doesn't pass push metadata to pre-push run commands reliably, so
# this script derives the range from git itself.
#
# Skip:  LEFTHOOK=0 git push
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

# Determine the pushed range from the remote-tracking state of the current
# branch. Falls back to the last 10 commits if no upstream is set.
REMOTE_REF="$(git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}' 2>/dev/null || true)"

if [ -n "$REMOTE_REF" ]; then
    BASE="$(git merge-base "$REMOTE_REF" HEAD)"
    RANGE="${BASE}..HEAD"
else
    # No upstream yet (first push of branch): test what's been committed since
    # divergence from the default branch, or the last 10 commits as a last resort.
    BASE="$(git merge-base HEAD origin/main 2>/dev/null || git merge-base HEAD origin/master 2>/dev/null || echo "")"
    if [ -n "$BASE" ]; then
        RANGE="${BASE}..HEAD"
    else
        RANGE="HEAD~10..HEAD"
    fi
fi

PUSHED_RS="$(git diff --name-only --diff-filter=ACMR "$RANGE" -- '*.rs' || true)"

# Map to crates.
crates=""
while IFS= read -r f; do
    [ -z "$f" ] && continue
    case "$f" in
        Carlo.rs/*) crates="${crates}carlo-rs " ;;
        QMC.rs/*)   crates="${crates}qmc-rs " ;;
        CMC.rs/*)   crates="${crates}cmc-rs " ;;
    esac
done <<EOF
$PUSHED_RS
EOF

# Dedupe.
seen="" PKG_FLAGS=""
for c in $crates; do
    case " $seen " in
        *" $c "*) ;;
        *) seen="$seen $c"; PKG_FLAGS="${PKG_FLAGS}-p ${c} " ;;
    esac
done

if [ -z "$PKG_FLAGS" ]; then
    PKG_FLAGS="--workspace"
    SCOPE="workspace"
else
    SCOPE="(crates: $(echo "$PKG_FLAGS" | sed 's/-p //g'))"
fi

printf '\033[1;36m▸ cargo test  %s\033[0m\n' "$SCOPE"
cargo test $PKG_FLAGS
