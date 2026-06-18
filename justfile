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

# Quick feedback (format + lint + test)
check:
    cargo fmt --check
    cargo clippy --all-targets -- -D warnings
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

# Test all
test:
    cargo test --workspace

# Test unit tests only
test-unit:
    cargo test --lib

# Test integration tests only
test-integration:
    cargo test --test '*'

# Test MPI backend (requires MPI and mpirun)
test-mpi:
    #!/usr/bin/env bash
    if command -v mpirun &> /dev/null; then
        mpirun -np 4 cargo test --features mpi --test mpi_test
    else
        echo "MPI not installed. Install with:"
        echo "  Ubuntu/Debian: sudo apt-get install libopenmpi-dev openmpi-bin"
        echo "  macOS: brew install open-mpi"
        exit 1
    fi

# Generate docs
doc:
    cargo doc --workspace --no-deps --open

# Check docs without opening
doc-check:
    cargo doc --workspace --no-deps

# Clean build artifacts
clean:
    cargo clean

# Security audit
audit:
    cargo audit 2>/dev/null || echo "cargo-audit not installed"

# Publish dry-run
publish-dry:
    cd Carlo.rs && cargo publish --dry-run

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