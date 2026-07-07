#!/usr/bin/env bash
# typos — pure spelling check across the repo (config: .typos.toml).
#
# Runs in pre-commit. If `typos` is not installed, exit 0 so contributors
# without it aren't blocked; CI installs it as the source of truth.
#
# Install:  cargo install typos-cli   (or: brew install typos-cli)
set -euo pipefail

if ! command -v typos >/dev/null 2>&1; then
    exit 0
fi

typos
