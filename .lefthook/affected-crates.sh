#!/usr/bin/env bash
# Print `-p <crate> -p <crate>` flags for the affected crates of the given file
# list (one path per arg). If no file maps to a known crate, prints nothing —
# caller falls back to --workspace.
#
# Usage:  affected-crates.sh <file> <file> ...
#   or:   affected-crates.sh "$@"     # when lefthook passes {staged_files}
#
# Crate layout:
#   Carlo.rs/ -> carlo-rs   QMC.rs/ -> qmc-rs   CMC.rs/ -> cmc-rs
#
# Used by lefthook.yml to scope `cargo clippy` / `cargo test` to touched crates.
set -eu

# Collect args (may be empty when nothing staged in a known crate dir).
crates=""
for f in "$@"; do
    case "$f" in
        Carlo.rs/*) crates="${crates}carlo-rs " ;;
        QMC.rs/*)   crates="${crates}qmc-rs " ;;
        CMC.rs/*)   crates="${crates}cmc-rs " ;;
    esac
done

# Dedupe while preserving order.
seen="" out=""
for c in $crates; do
    case " $seen " in
        *" $c "*) ;;
        *) seen="$seen $c"; out="${out}-p ${c} " ;;
    esac
done

echo "$out"
