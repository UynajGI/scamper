//! Carlo.rs - Monte Carlo simulation framework CLI

use carlo_rs::cli;

fn main() {
    // Initialize logging
    tracing_subscriber::fmt::init();

    if let Err(e) = cli::run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
