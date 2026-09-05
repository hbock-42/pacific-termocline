//! Running a scenario forward in time and writing what it produced.
//!
//! This is the join between the five epics behind it: [`Scenario`] loads the
//! config (T-03.4) and builds the [`Basin`](crate::Basin) (T-04.1) and the
//! wind forcing (Epic 03), [`Solver`] takes the RK4 steps at the CFL-checked
//! timestep (Epics 01–02), and [`RunWriter`] saves the ones the
//! [`OutputSchedule`] asks for (T-05.2). Nothing physical or numerical is
//! decided here — every bound was checked by whichever constructor owns it,
//! and this module's whole job is the order the pieces go in.
//!
//! # The loop
//!
//! A run of `total_steps` steps visits `total_steps + 1` states: the initial
//! one and the one after each step. Each is offered to the schedule, which
//! keeps the multiples of the output cadence, so the frame at index `k` is the
//! state after `k · N` steps at model time `k · N · dt`. The step is taken
//! *after* the state is offered, which is what puts the initial condition in
//! the file as frame zero and leaves the final state as the last frame.
//!
//! Model time is `step · dt` rather than an accumulator advanced by `dt` each
//! iteration: the product is one rounding, the accumulator is one per step,
//! and the frame times a reader plots against have to be the times the forcing
//! was actually sampled at.
//!
//! # The SST coupling
//!
//! A scenario's `[sst]` section is read here and nowhere else: it decides
//! which solver is built and which state is allocated, and the time loop is
//! the same either way. `T'` is *not* written to the run's frames — the
//! interchange format of [ADR-0004] describes the three variables of the
//! linear core, and extending it is a change to the contract the visualizer
//! reads, not a side effect of adding a term (T-12.3's business).
//!
//! [ADR-0004]: ../../docs/planning/adr/0004-data-interchange-format.md
//!
//! # What is allocated, and when
//!
//! The state, the solver's stage buffers, the composite wind and the stress
//! field written into the frames are all built before the first step and
//! reused (CODING_STANDARDS.md § *Performance*). What remains per *frame* is
//! the frame itself, which [`RunWriter::append`] builds — at the output
//! cadence, not in the time-stepping loop.
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

use termocline_format::{FormatError, GridSpec, RunHeader};

use crate::coriolis::BetaPlane;
use crate::forcing::WindForcing;
use crate::progress::RunObserver;
use crate::run_writer::{RunWriteError, RunWriter};
use crate::scenario::{Scenario, ScenarioError};
use crate::solver::{Solver, SolverError};
use crate::state::OceanState;

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

/// What a finished run wrote: the numbers the CLI reports back and a test
/// asserts on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunReport {
    /// Steps the run took from its initial state.
    steps_taken: u64,
    /// Frames written, including the one holding the initial state.
    frames_written: u64,
}

impl RunReport {
    /// Steps the run took from its initial state.
    #[must_use]
    pub const fn steps_taken(self) -> u64 {
        self.steps_taken
    }

    /// Frames written, including the one holding the initial state.
    #[must_use]
    pub const fn frames_written(self) -> u64 {
        self.frames_written
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
    let basin = scenario.basin();
    let grid = basin.grid();
    let params = scenario.physical_params();
    let schedule = scenario.output_schedule();

    // The beta-plane is placed by the basin's southern boundary rather than
    // centred on the grid, so a scenario that does not straddle the equator
    // gets the rotation of the latitudes it actually covers. It is the same
    // plane the scenario loader checked the rotation bound against, because
    // both ask the basin for it.
    let plane = BetaPlane::of_basin(params, basin);
    // The Epic 12 switch, and the whole of it: a scenario with an `[sst]`
    // section integrates the mixed-layer anomaly `T'` alongside the core, one
    // without integrates the three-variable model of Epics 01-07. The two
    // branches differ in one term of the right-hand side and one field of the
    // state; everything below is the same loop (T-12.1).
    let mut solver = match scenario.sst_params() {
        Some(sst) => {
            Solver::coupled_to_sst(grid, basin.spacing(), params, plane, schedule.dt_s(), sst)?
        }
        None => Solver::new(grid, basin.spacing(), params, plane, schedule.dt_s())?,
    };

    let header = RunHeader::new(
        GridSpec::new(grid.nx(), grid.ny(), scenario.bounds().into())?,
        params.into(),
        description,
        schedule.timing(),
    );

    let mut state = if solver.couples_sst() {
        OceanState::at_rest_with_sst_anomaly(grid)
    } else {
        OceanState::at_rest(grid)
    };
    // The run's forcing: the scenario's wind and the one field it is sampled
    // into, held here rather than inside the solver so that it survives the
    // whole time loop. A steady wind is therefore sampled once for the run
    // (T-10.5, `docs/performance-notes.md`), and the field a frame records is
    // the very field that stage of the integration read — the same instant,
    // and now literally the same buffer.
    let mut forcing = WindForcing::new(basin, scenario.wind());

    let mut writer = RunWriter::create(directory, &header)?;
    observer.run_started(description, schedule);
    let mut frames_written = 0;
    for step in 0..=schedule.total_steps() {
        let t_s = schedule.model_time_at_step(step);
        if schedule.writes_at_step(step) {
            writer.append(t_s, &state, forcing.at(t_s))?;
            observer.frame_written(frames_written, t_s);
            frames_written += 1;
        }
        if step < schedule.total_steps() {
            solver.step_with_forcing(&mut state, t_s, &mut forcing);
            // The step just taken is `step + 1` of the run, and it reached the
            // model time of the *next* iteration — which is what the observer
            // reports, so the time on screen is the time the state is at.
            observer.step_taken(step + 1, schedule.model_time_at_step(step + 1));
        }
    }
    writer.finish()?;

    let report = RunReport {
        steps_taken: schedule.total_steps(),
        frames_written,
    };
    observer.run_finished(&report);
    Ok(report)
}
