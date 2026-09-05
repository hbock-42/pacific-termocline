//! CLI entry point for the engine.
//!
//! The binary is a set of subcommands over one run directory. `inspect` is
//! here (T-06.4); the scenario runner behind `run` lands with the rest of
//! Epic 06.
//!
//! Anything the user got wrong — a run that is not there, a header this build
//! cannot read — leaves through [`std::process::exit`] with a message on
//! stderr naming the run, per CODING_STANDARDS.md § *Correctness and failure*:
//! invalid input is a `Result` all the way up, and never a stack trace.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, Subcommand};

/// The engine: run a scenario, or look at a run one has already produced.
#[derive(Debug, Parser)]
#[command(name = "termocline", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Print a run's header — grid, physical parameters, scenario, frame
    /// count — without reading its frames.
    Inspect {
        /// The run directory to report on: the one holding `header.json` and
        /// `frames.bin`.
        #[arg(long = "run", value_name = "DIR")]
        run: PathBuf,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match cli.command {
        Command::Inspect { run } => match engine::inspect_run(&run) {
            Ok(summary) => {
                print!("{summary}");
                ExitCode::SUCCESS
            }
            Err(error) => {
                eprintln!("could not inspect {}: {error}", run.display());
                ExitCode::FAILURE
            }
        },
    }
}
