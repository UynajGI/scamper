// CLI tests - using clap Parser directly
use clap::Parser;

// Recreate the CLI structure for testing
#[derive(clap::Parser)]
#[command(name = "carlo-rs")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(clap::Subcommand)]
enum Commands {
    Run {
        #[arg(short, long)]
        single: bool,
        #[arg(short, long)]
        restart: bool,
    },
    Status,
    Merge,
    Delete,
}

#[test]
fn test_cli_run_command() {
    let cli = Cli::try_parse_from(["carlo-rs", "run"]).unwrap();
    assert!(matches!(cli.command, Commands::Run { .. }));
}

#[test]
fn test_cli_status_command() {
    let cli = Cli::try_parse_from(["carlo-rs", "status"]).unwrap();
    assert!(matches!(cli.command, Commands::Status));
}

#[test]
fn test_cli_run_with_options() {
    let cli = Cli::try_parse_from(["carlo-rs", "run", "--single", "--restart"]).unwrap();
    match cli.command {
        Commands::Run { single, restart } => {
            assert!(single);
            assert!(restart);
        }
        _ => panic!("Expected Run command"),
    }
}
