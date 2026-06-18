#!/usr/bin/env bash
# cargo fmt --check (workspace-wide — cargo has no per-package fmt filter).
# Invoked by lefthook from the pre-commit group; the glob filter in
# lefthook.yml ensures this only runs when *.rs files are staged.
set -euo pipefail

if ! cargo fmt --check; then
    echo "✗ cargo fmt --check failed. Run \`cargo fmt\` then re-stage." >&2
    exit 1
fi
