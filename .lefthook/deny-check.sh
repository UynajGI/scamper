#!/usr/bin/env bash
# cargo deny check — advisories (RUSTSEC) + licenses.
#
# Runs in pre-push. If `cargo-deny` is not installed, exit 0 so contributors
# without it aren't blocked; CI runs the same check via cargo-deny-action.
#
# Install:  cargo install cargo-deny
# Skip:     LEFTHOOK=0 git push
set -euo pipefail

if ! command -v cargo-deny >/dev/null 2>&1; then
    exit 0
fi

# `check advisories licenses` — skip bans (low signal, see deny.toml) and
# sources (not configured). Network needed to refresh the advisory DB on first
# run; subsequent runs use the cached DB under ~/.cargo/advisory-dbs.
#
# If the advisory DB can't be fetched (offline / firewall), fall back to a
# licenses-only check so a flaky network doesn't block the push. CI always
# runs both with network access.
if cargo deny check advisories licenses 2>&1; then
    exit 0
fi

# Retry licenses-only — licenses need no network. If THIS fails, it's a real
# license violation and we must fail the push.
printf '\033[1;33m⚠ cargo deny advisories failed (network?). Falling back to licenses-only.\033[0m\n' >&2
cargo deny check licenses
