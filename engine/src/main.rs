//! CLI entry point for the engine.
//!
//! The binary is a set of subcommands over one run directory: `run` takes a
//! scenario file and produces a run (T-06.1), and `inspect` reports on one that
//! already exists (T-06.4).
//!
//! Anything the user got wrong — a scenario that is not there, a timestep the
//! scheme refuses, a run this build cannot read — comes back as an
//! [`ExitCode::FAILURE`] and a message on stderr naming what was asked for, per
//! CODING_STANDARDS.md § *Correctness and failure*: invalid input is a `Result`
//! all the way up, and never a stack trace.

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
    /// Run a scenario forward in time and write the result as a run
    /// directory.
    Run {
        /// The scenario to run: the TOML file described by
        /// `docs/scenario-config-reference.md`.
        #[arg(long = "config", value_name = "FILE")]
        config: PathBuf,
        /// Where to write the run: a directory, created if it is not there,
        /// that will hold `header.json` and `frames.bin`.
        #[arg(long = "out", value_name = "DIR")]
        out: PathBuf,
    },
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
        Command::Run { config, out } => match engine::run_scenario_file(&config, &out) {
            Ok(report) => {
                // One line, at the end: what the run wrote and where. Progress
                // during the run is T-06.2's, not this command's.
                println!(
                    "{}: {} steps, {} frames written to {}",
                    config.display(),
                    report.steps_taken(),
                    report.frames_written(),
                    out.display()
                );
                ExitCode::SUCCESS
            }
            Err(error) => {
                eprintln!("could not run {}: {error}", config.display());
                ExitCode::FAILURE
            }
        },
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
