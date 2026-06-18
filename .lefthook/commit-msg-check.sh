#!/usr/bin/env bash
# Validate a commit message against Scuttle's Conventional Commits convention.
# The message file path is $1 (lefthook passes {1}).
#
# Format:  <type>(<scope>)!: <subject>
#   subject ≤ 72 chars, imperative, no trailing period
#   body (optional) separated by one blank line, wrap ≤ 100 (warning only)
#
# Allowed types:  feat fix docs style refactor perf test build ci chore revert
# Allowed scopes: Carlo QMC CMC carlo-rs qmc-rs cmc-rs spec task docs deps release
#
# Pass-through (never blocked): Merge / Revert / squash! / fixup! / amend!
set -euo pipefail

MSG_FILE="${1:?usage: commit-msg-check.sh <msg-file>}"
SUBJECT="$(sed -n '1p' "$MSG_FILE")"

red()    { printf '\033[1;31m%s\033[0m\n' "$1" >&2; }
yellow() { printf '\033[1;33m%s\033[0m\n' "$1" >&2; }

# --- pass-through commits --------------------------------------------------
case "$SUBJECT" in
    "Merge "*|"Revert "*) exit 0 ;;
    "squash! "*|"fixup! "*|"amend! "*) exit 0 ;;
esac

# --- regex -----------------------------------------------------------------
TYPES='feat|fix|docs|style|refactor|perf|test|build|ci|chore|revert'
SCOPES='Carlo|QMC|CMC|carlo-rs|qmc-rs|cmc-rs|spec|task|docs|deps|release'
HEADER_RE="^(${TYPES})(\((${SCOPES})\))?(!)?: .+"

fail() {
    red "✗ Commit message rejected."
    cat >&2 <<'EOF'

  Expected Conventional Commits format:

      <type>(<scope>)!: <subject>

  Examples:
      feat(CMC): add Wolff cluster algorithm
      fix(CMC): correct beta sign in metropolis acceptance
      docs: update CLAUDE.md architecture table
      refactor(qmc-rs)!: replace SSE with worldline QMC
      chore: bump version to 0.1.0-dev1

  Allowed types:  feat fix docs style refactor perf test build ci chore revert
  Allowed scopes: Carlo QMC CMC carlo-rs qmc-rs cmc-rs spec task docs deps release
                   (scope optional)

  Subject rules: imperative mood, ≤ 72 chars, no trailing period.
  Skip:          LEFTHOOK=0 git commit
EOF
    exit 1
}

# 1. Header shape.
if ! [[ "$SUBJECT" =~ $HEADER_RE ]]; then
    fail
fi

SCOPE="${BASH_REMATCH[3]}"

# 2. Subject length.
if ((${#SUBJECT} > 72)); then
    red "✗ Subject line is ${#SUBJECT} chars (max 72)."
    cat >&2 <<'EOF'

  Wrap the subject at 72 characters. Move detail to the body.

EOF
    exit 1
fi

# 3. No trailing period.
if [[ "$SUBJECT" =~ \.$ ]]; then
    red "✗ Subject ends with a period — remove it."
    exit 1
fi

# 4. Unknown scope → warn only (forward-compatible with future scopes).
if [ -n "$SCOPE" ]; then
    case "$SCOPE" in
        Carlo|QMC|CMC|carlo-rs|qmc-rs|cmc-rs|spec|task|docs|deps|release) ;;
        *) yellow "⚠ Unknown scope '$SCOPE' (allowed but unusual). Continuing." ;;
    esac
fi

# 5. Body wrap warning (trailers exempt).
if [ "$(sed -n '2,$p' "$MSG_FILE" | wc -l)" -gt 0 ]; then
    long_lines="$(awk 'NR>1 && length($0)>100 && $0 !~ /^(Signed-off-by|Co-Authored-By|Reviewed-by|Refs|Fixes|Closes|Resolves):/ {printf "  line %d (%d chars)\n", NR, length($0)}' "$MSG_FILE")"
    if [ -n "$long_lines" ]; then
        yellow "⚠ Some body lines exceed 100 chars (trailers ignored):"
        yellow "$long_lines"
    fi
fi

exit 0
