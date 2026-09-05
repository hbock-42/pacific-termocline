//! CLI entry point for the engine.
//!
//! The binary is a set of subcommands over one run directory: `run` takes a
//! scenario file and produces a run (T-06.1), and `inspect` reports on one that
//! already exists (T-06.4).
//!
//! A run reports its progress to stderr as it goes (T-06.2), so a multi-minute
//! run is not a black box; `--quiet` suppresses that for a scripted or CI run,
//! and `--verbose` adds the run's per-frame detail. stdout carries the run's
//! one-line summary and nothing else, so a script may read it without
//! filtering.
//!
//! Anything the user got wrong — a scenario that is not there, a timestep the
//! scheme refuses, a run this build cannot read — comes back as an
//! [`ExitCode::FAILURE`] and a message on stderr naming what was asked for, per
//! CODING_STANDARDS.md § *Correctness and failure*: invalid input is a `Result`
//! all the way up, and never a stack trace.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand};

use engine::progress::{RunProgress, Verbosity};

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
        #[command(flatten)]
        reporting: Reporting,
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

/// How loudly a run reports on itself.
///
/// The two flags are contradictory instructions, so clap refuses them together
/// rather than the binary picking one silently.
#[derive(Debug, Args)]
struct Reporting {
    /// Suppress progress and log output, for a scripted or CI run.
    #[arg(long = "quiet", short = 'q')]
    quiet: bool,
    /// Log the run's per-frame detail as well as its progress.
    #[arg(long = "verbose", short = 'v', conflicts_with = "quiet")]
    verbose: bool,
}

impl Reporting {
    /// The verbosity these flags ask for.
    const fn verbosity(&self) -> Verbosity {
        if self.quiet {
            Verbosity::Quiet
        } else if self.verbose {
            Verbosity::Verbose
        } else {
            Verbosity::Normal
        }
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match cli.command {
        Command::Run {
            config,
            out,
            reporting,
        } => {
            // Progress goes to stderr, so the summary below stays the only
            // thing on stdout.
            let mut progress = RunProgress::to_stderr(reporting.verbosity());
            match engine::run_scenario_file_observed(&config, &out, &mut progress) {
                Ok(report) => {
                    // One line, at the end: what the run wrote and where.
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
            }
        }
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
