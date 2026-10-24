# Install lefthook git hooks (pre-commit, commit-msg, pre-push)
hooks:
    lefthook install
    @echo "lefthook hooks installed (verify: lefthook check-install)"

# Show lefthook install status
hooks-status:
    @lefthook check-install 2>&1 || echo "(lefthook not installed — run: just hooks)"

# Uninstall lefthook git hooks
hooks-disable:
    lefthook uninstall
    @echo "lefthook hooks uninstalled"

# Quick feedback (format + lint + test). Lint level is set in
# [workspace.lints] (Cargo.toml) — clippy here picks it up automatically.
check:
    cargo fmt --check
    cargo clippy --all-targets --workspace
    cargo test --workspace

# Format code
fmt:
    cargo fmt

# Build release
build:
    cargo build --release

# Build with MPI support (requires MPI installed)
build-mpi:
    cargo build --release --features mpi

# Test all quick (non-ignored) workspace tests. Long #[ignore] physics tests
# are CI-only by default; run them locally with `just test-long`.
test: test-fast

# Quick tests, parallel across crates: build every test binary once, then run
# each crate's suite concurrently. Tests inside a binary already use libtest's
# thread pool, so this only adds the cross-crate dimension.
test-fast:
    #!/usr/bin/env bash
    set -uo pipefail

    if ! cargo test --workspace --no-run; then
        exit 1
    fi

    crates=(carlo-rs cmc-rs qmc-rs mcmc-rs scamper)
    pids=()
    for crate in "${crates[@]}"; do
        cargo test -p "$crate" --no-fail-fast &
        pids+=("$!")
    done

    status=0
    for pid in "${pids[@]}"; do
        wait "$pid" || status=1
    done
    exit "$status"

# Quick tests for the crates touched since the upstream merge base. Falls back
# to the full quick suite when no crate maps (lockfile-only or docs changes).
test-changed:
    #!/usr/bin/env bash
    set -euo pipefail

    base="$(git merge-base HEAD '@{upstream}' 2>/dev/null || git merge-base HEAD origin/main 2>/dev/null || true)"
    if [ -n "$base" ]; then
        changed="$(git diff --name-only --diff-filter=ACMR "$base" -- '*.rs' '*/Cargo.toml' || true)"
    else
        changed="$(git diff --name-only --diff-filter=ACMR HEAD -- '*.rs' '*/Cargo.toml' || true)"
    fi
    changed="$(printf '%s\n%s\n' "$changed" "$(git ls-files --others --exclude-standard -- '*.rs' '*/Cargo.toml')" | sed '/^$/d')"

    flags="$(./.lefthook/affected-crates.sh $changed)"
    if [ -z "$flags" ]; then
        echo "No crate-scoped changes; running the full quick suite."
        just test-fast
    else
        echo "Testing affected crates: ${flags//-p /}"
        cargo test $flags --no-fail-fast
    fi

# Long #[ignore] physics tests (release, all features). CI runs these in
# nightly.yml; locally they need libhdf5/libopenmpi (`just install-deps`).
test-long:
    cargo test --release --workspace --all-features --no-fail-fast -- --ignored --nocapture

# Quick tests with every optional feature enabled (hdf5, mpi, allocation
# probe). Requires the system libraries installed by `just install-deps`.
test-all:
    cargo test --workspace --all-features --no-fail-fast

# Test unit tests only
test-unit:
    cargo test --lib

# Test integration tests only
test-integration:
    cargo test --test '*'

# Test the MPI backend (requires MPI and mpirun; this machine uses MPICH).
# Each MPI test runs in its OWN mpirun invocation with an --exact filter:
# MPI cannot be initialized twice in one process, and a single cargo test
# process per rank would otherwise run every matched #[test] back-to-back.
# Usage: just mpi-test [np=2] — np 1/2/4 recommended. Rank-count notes:
# np 1 runs only the singleton-safe suites; the pt_exchange end-to-end test
# hardcodes a 2-chain config (its entry point owns the MPI init, so it
# cannot probe the world size first) and therefore runs at np 2 only.
mpi-test np="2":
    #!/usr/bin/env bash
    if command -v mpirun &> /dev/null; then
        if [ "{{ np }}" -eq 2 ]; then
            tests=(
                "mpi_backend_distributed::mpi_backend_distributed_suite"
                "mpi_test::mpi_tests::mpi_backend_smoke_suite"
                "mpi_distributed::distributed_tests::dynamic_scheduler_completes_and_merges_runs"
                "mpi_pt_exchange::pt_exchange_completes_and_returns_results"
                "mpi_pt_dynamics::pt_exchange_dynamics_suite"
            )
        elif [ "{{ np }}" -gt 2 ]; then
            tests=(
                "mpi_backend_distributed::mpi_backend_distributed_suite"
                "mpi_test::mpi_tests::mpi_backend_smoke_suite"
                "mpi_distributed::distributed_tests::dynamic_scheduler_completes_and_merges_runs"
                "mpi_pt_dynamics::pt_exchange_dynamics_suite"
            )
        else
            tests=(
                "mpi_backend_distributed::mpi_backend_distributed_suite"
                "mpi_pt_dynamics::pt_exchange_dynamics_suite"
            )
        fi
        for test in "${tests[@]}"; do
            echo "=== mpirun -np {{ np }} :: ${test} ==="
            mpirun -np {{ np }} cargo test -p carlo-rs --features mpi,hdf5 --test suite -- --ignored --nocapture --exact "${test}" || exit 1
        done
    else
        echo "MPI not installed. Install with:"
        echo "  Ubuntu/Debian: sudo apt-get install libopenmpi-dev openmpi-bin"
        echo "  macOS: brew install open-mpi"
        exit 1
    fi

# Multi-seed z-score monitoring — local equivalent of the nightly.yml
# `zscore-monitor` job (P2.8). Raises the seed count of every z-score test
# via SCUTTLE_ZSCORE_SEEDS. Usage: just nightly-zscore [seeds=64]
nightly-zscore seeds="64":
    SCUTTLE_ZSCORE_SEEDS={{seeds}} cargo test --release -p cmc-rs --all-features --test suite zscore -- --nocapture --include-ignored
    SCUTTLE_ZSCORE_SEEDS={{seeds}} cargo test --release -p qmc-rs --all-features --test suite zscore -- --nocapture --include-ignored
    SCUTTLE_ZSCORE_SEEDS={{seeds}} cargo test --release -p qmc-rs --all-features --test suite ergodicity_multi_seed -- --nocapture --include-ignored

# Generate docs
doc:
    cargo doc --workspace --no-deps --open

# Check docs without opening
doc-check:
    cargo doc --workspace --no-deps

# Clean build artifacts
clean:
    cargo clean

# Dependency audit (advisories + licenses). Requires cargo-deny:
#   cargo install cargo-deny
# CI runs the same via cargo-deny-action.
deny:
    cargo deny --all-features check advisories licenses

# Spelling check (config: .typos.toml). Requires typos:
#   cargo install typos-cli
typos:
    typos

# Cargo rewrites workspace paths to registry dependencies during packaging, so
# dependent dev versions cannot be packaged until carlo-rs has landed. The tag
# workflow publishes and verifies the remaining dependency layers in order.
publish-dry:
    cargo package -p carlo-rs --no-verify

# Run benchmarks
bench:
    cargo bench --manifest-path Carlo.rs/Cargo.toml

# Compare with Julia baseline (requires Julia)
bench-compare:
    cargo bench --manifest-path Carlo.rs/Cargo.toml
    @echo "To compare with Carlo.jl, run:"
    @echo "  cd Carlo.jl && julia --project -e 'using Pkg; Pkg.instantiate(); include(\"benchmark/bench.jl\")'"

# Quick compile check
check-fast:
    cargo check

# Run single-threaded test
test-st:
    cargo test -- --test-threads=1

# Install system dependencies (Ubuntu/Debian)
install-deps:
    #!/usr/bin/env bash
    if command -v apt-get &> /dev/null; then
        sudo apt-get update
        sudo apt-get install -y libhdf5-dev openmpi-bin libopenmpi-dev
    else
        echo "Please install dependencies manually:"
        echo "  libhdf5-dev: HDF5 library"
        echo "  openmpi-bin libopenmpi-dev: MPI library"
    fi