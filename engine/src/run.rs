//! Running a scenario forward in time and writing what it produced.
//!
//! This is the join between the five epics behind it: [`Scenario`] loads the
//! config (T-03.4), [`RunLoop`] builds the basin, the solver and the wind
//! forcing and takes the steps (Epics 01-04), and [`RunWriter`] saves the ones
//! the schedule asks for (T-05.2). Nothing physical or numerical is decided
//! here — every bound was checked by whichever constructor owns it, and this
//! module's whole job is to drive the loop into a directory.
//!
//! # The loop
//!
//! The loop itself is [`RunLoop`]'s, and lives there because the browser takes
//! the same steps without a filesystem to write them to (ADR-0012). What this
//! module adds is the two things a *file* needs: the writer the frames go to,
//! and the observer that reports progress to a terminal (T-06.2). A run of
//! `total_steps` steps therefore visits `total_steps + 1` states here exactly
//! as it does there, and the frame at index `k` is the state after `k · N`
//! steps at model time `k · N · dt`.
//!
//! # What is allocated, and when
//!
//! The state, the solver's stage buffers, the composite wind and the stress
//! field written into the frames are all built by [`RunLoop::of_scenario`]
//! before the first step and reused (CODING_STANDARDS.md § *Performance*).
//! What remains per *frame* is the frame itself, which [`RunWriter::append`]
//! builds — at the output cadence, not in the time-stepping loop.
//!
//! # Failure
//!
//! Every way this can fail is the user's scenario or the user's filesystem, so
//! all of them are [`RunError`] rather than a panic
//! (CODING_STANDARDS.md § *Correctness and failure*). The scenario is built,
//! and the solver constructed, before the run directory is opened, so a
//! scenario the engine refuses leaves no half-written run behind.

use std::fmt;
use std::path::Path;

use termocline_format::FormatError;

use crate::progress::{RunObserver, RunReport};
use crate::run_loop::{RunLoop, RunLoopError};
use crate::run_writer::{RunWriteError, RunWriter};
use crate::scenario::{Scenario, ScenarioError};
use crate::solver::SolverError;

/// Why a run could not be made.
///
/// Each variant names the stage that objected and carries that stage's own
/// error, which already names the offending value and the bound it violated;
/// nothing is re-worded here, because the message the constructor wrote is the
/// actionable one.
#[derive(Debug)]
pub enum RunError {
    /// The scenario could not be read, or is not one the engine will run.
    Scenario(ScenarioError),
    /// The scenario's grid is not one the output format can describe.
    Grid(FormatError),
    /// The scenario asked for a timestep the scheme cannot take.
    Solver(SolverError),
    /// The run directory, its header or one of its frames could not be
    /// written.
    Write(RunWriteError),
}

impl fmt::Display for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Scenario(source) => source.fmt(f),
            Self::Grid(source) => write!(f, "[basin]: {source}"),
            Self::Solver(source) => write!(f, "[run]: {source}"),
            Self::Write(source) => source.fmt(f),
        }
    }
}

impl std::error::Error for RunError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Scenario(source) => Some(source),
            Self::Grid(source) => Some(source),
            Self::Solver(source) => Some(source),
            Self::Write(source) => Some(source),
        }
    }
}

impl From<ScenarioError> for RunError {
    fn from(source: ScenarioError) -> Self {
        Self::Scenario(source)
    }
}

impl From<FormatError> for RunError {
    fn from(source: FormatError) -> Self {
        Self::Grid(source)
    }
}

impl From<SolverError> for RunError {
    fn from(source: SolverError) -> Self {
        Self::Solver(source)
    }
}

impl From<RunWriteError> for RunError {
    fn from(source: RunWriteError) -> Self {
        Self::Write(source)
    }
}

/// The two ways a run refuses to start, as the `run` command reports them.
///
/// [`RunLoop`] is the portable half of this module (ADR-0012) and objects
/// before a directory is opened; its two errors are two of this one's four,
/// under the same names.
impl From<RunLoopError> for RunError {
    fn from(source: RunLoopError) -> Self {
        match source {
            RunLoopError::Grid(source) => Self::Grid(source),
            RunLoopError::Solver(source) => Self::Solver(source),
        }
    }
}

/// Run the scenario in `config_path` and write it into `directory`.
///
/// The whole `run` command: load, run, write. The run is named after the
/// config file it came from — the file stem, or the path as written if it has
/// none — which is what the header carries as its scenario description and
/// what `inspect` prints back.
///
/// # Errors
/// A [`RunError`] naming the stage that objected: a scenario that could not be
/// read or built, a grid the output format cannot describe, a timestep the
/// scheme refuses, or a run directory that could not be written.
pub fn run_scenario_file(config_path: &Path, directory: &Path) -> Result<RunReport, RunError> {
    run_scenario_file_observed(config_path, directory, &mut ())
}

/// [`run_scenario_file`], reporting what it does to `observer`.
///
/// The variant the CLI calls: the observer is where progress and logging live
/// (T-06.2), so that this module decides nothing about a terminal.
///
/// # Errors
/// The errors of [`run_scenario_file`].
pub fn run_scenario_file_observed(
    config_path: &Path,
    directory: &Path,
    observer: &mut dyn RunObserver,
) -> Result<RunReport, RunError> {
    let scenario = Scenario::load(config_path)?;
    let description = config_path.file_stem().map_or_else(
        || config_path.display().to_string(),
        |stem| stem.to_string_lossy().into_owned(),
    );
    run_scenario_observed(&scenario, &description, directory, observer)
}

/// Run `scenario` to the end of its output schedule, writing it into
/// `directory` under the name `description`.
///
/// `directory` is created if it is not there, and its two files are truncated
/// if they are: a run writes a run directory, it does not merge into one.
///
/// # Errors
/// A [`RunError`] naming the stage that objected. The solver is built before
/// the directory is opened, so a scenario refused for its timestep leaves
/// nothing on disk.
pub fn run_scenario(
    scenario: &Scenario,
    description: &str,
    directory: &Path,
) -> Result<RunReport, RunError> {
    run_scenario_observed(scenario, description, directory, &mut ())
}

/// [`run_scenario`], reporting what it does to `observer`.
///
/// The run tells the observer four things: that it started, that a step was
/// taken, that a frame was written, and that it finished. What any of that
/// looks like — a progress bar, a log line, nothing at all — is the observer's
/// business, not this module's.
///
/// # Errors
/// The errors of [`run_scenario`].
pub fn run_scenario_observed(
    scenario: &Scenario,
    description: &str,
    directory: &Path,
    observer: &mut dyn RunObserver,
) -> Result<RunReport, RunError> {
    let mut run = RunLoop::of_scenario(scenario, description)?;
    let schedule = run.schedule();

    let mut writer = RunWriter::create(directory, run.header())?;
    observer.run_started(description, schedule);
    let mut frames_written = 0;
    loop {
        if let Some(saved) = run.take_frame() {
            writer.append(saved.t_s, saved.state, saved.wind_stress)?;
            observer.frame_written(frames_written, saved.t_s);
            frames_written += 1;
        }
        if !run.take_step() {
            break;
        }
        // The step just taken reached the model time the loop is now at —
        // which is what the observer reports, so the time on screen is the
        // time the state is at.
        observer.step_taken(run.steps_taken(), run.model_time_s());
    }
    writer.finish()?;

    let report = RunReport::new(schedule.total_steps(), frames_written);
    observer.run_finished(&report);
    Ok(report)
}
