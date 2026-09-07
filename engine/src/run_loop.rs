//! A run taken one step at a time, without a filesystem.
//!
//! [`crate::run`] is a run driven to completion into a directory. This is the
//! same run with the loop turned inside out: the caller holds the [`RunLoop`]
//! and decides when the next step is taken, which is what [ADR-0012] asks for
//! on the web. A browser tab cannot block its main thread for the minutes a
//! full run takes, so it steps in chunks and yields between them — and a
//! chunk is `take_step` called a few times, nothing more.
//!
//! # Why the loop lives here and not in the visualizer
//!
//! Everything before the first step is a decision: which solver a scenario's
//! `[sst]` section asks for, which state that solver needs, which shape the
//! winds take, and what the header therefore promises. ADR-0012's
//! reproducibility argument is that a browser run and a native run are the
//! same computation, and that holds only while there is one place those
//! decisions are made. So [`crate::run`] drives this loop too: the native
//! `run` command and the browser take the same steps in the same order, and
//! the only thing that differs is where the frames go.
//!
//! # The loop
//!
//! A run of `total_steps` steps visits `total_steps + 1` states: the initial
//! one, and the one after each step. Each is offered to the [`OutputSchedule`],
//! which keeps the multiples of the output cadence. A driver therefore reads
//!
//! ```text
//! loop {
//!     if let Some(frame) = run.take_frame() { save(frame) }
//!     if !run.take_step() { break }
//! }
//! ```
//!
//! and the frame at index `k` is the state after `k · N` steps at model time
//! `k · N · dt`, which is what puts the initial condition in the run as frame
//! zero and leaves the final state as the last frame.
//!
//! [ADR-0012]: ../../docs/planning/adr/0012-the-browser-runs-the-engine.md

use core::fmt;

use termocline_format::{FormatError, Frame, GridSpec, RunHeader};

use crate::basin::Basin;
use crate::coriolis::BetaPlane;
use crate::forcing::{CompositeWind, StageForcing, WindForcing, WindStressField};
use crate::run_writer::OutputSchedule;
use crate::scenario::Scenario;
use crate::solver::{Solver, SolverError};
use crate::state::OceanState;
use crate::wind_response::{CoupledWind, SstWindResponse};

/// Why a run could not be started.
///
/// Both variants describe a scenario the engine will not run — a grid the
/// output format cannot describe, or a timestep the scheme refuses — so they
/// are returned rather than panicked, and each carries the error of whichever
/// constructor objected (CODING_STANDARDS.md § *Correctness and failure*).
#[derive(Debug)]
pub enum RunLoopError {
    /// The scenario's grid is not one the output format can describe.
    Grid(FormatError),
    /// The scenario asked for a timestep the scheme cannot take.
    Solver(SolverError),
}

impl fmt::Display for RunLoopError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Grid(source) => write!(f, "[basin]: {source}"),
            Self::Solver(source) => write!(f, "[run]: {source}"),
        }
    }
}

impl std::error::Error for RunLoopError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Grid(source) => Some(source),
            Self::Solver(source) => Some(source),
        }
    }
}

impl From<FormatError> for RunLoopError {
    fn from(source: FormatError) -> Self {
        Self::Grid(source)
    }
}

impl From<SolverError> for RunLoopError {
    fn from(source: SolverError) -> Self {
        Self::Solver(source)
    }
}

/// One saved timestep, as the run holds it: the state, and the wind stress
/// that drove the run to it.
///
/// Borrowed rather than copied, and the two together rather than either alone:
/// a frame records the stress the step actually read, at the instant it read
/// it (`crate::run_writer::frame_of`), so handing out one without the other
/// would let a caller pair a state with a stress from another moment.
#[derive(Debug)]
pub struct SavedStep<'a> {
    /// Model time of this state, in seconds since the start of the run.
    pub t_s: f64,
    /// The prognostic state: `h`, `u`, `v`, and `T'` where the run couples it.
    pub state: &'a OceanState,
    /// The wind stress field the run is being forced by at `t_s`, in pascals.
    pub wind_stress: &'a WindStressField,
}

impl SavedStep<'_> {
    /// This step as the frame a run records for it, on `grid`.
    ///
    /// The three fields travel together because a frame needs all three at
    /// once, so this is where they are handed over together as well:
    /// [`crate::frame_of`] is the conversion, and this is it applied to the
    /// step the loop is on.
    ///
    /// # Errors
    /// [`FormatError`] if the run's state does not cover the basin `grid`
    /// describes — which would mean the loop and the header it wrote disagree.
    pub fn frame(&self, grid: &GridSpec) -> Result<Frame, FormatError> {
        crate::run_writer::frame_of(self.t_s, grid, self.state, self.wind_stress)
    }
}

/// A run in progress: everything it needs to take its next step, and how far
/// it has got.
///
/// Built from a [`Scenario`] before the first step, and then driven by
/// [`RunLoop::take_frame`] and [`RunLoop::take_step`]. Every buffer the run
/// uses is allocated here and reused (CODING_STANDARDS.md § *Performance*);
/// what remains per *frame* is the frame itself, which the driver builds at
/// the output cadence and not inside the loop.
pub struct RunLoop {
    /// The scheme, with its stage buffers.
    solver: Solver,
    /// The prognostic state, stepped in place.
    state: OceanState,
    /// The run's winds, and the one field they are sampled into.
    forcing: RunForcing,
    /// How long the run is and which of its steps are saved.
    schedule: OutputSchedule,
    /// What the run promises about itself, built before the first step because
    /// a reader takes the frame count from it (`crate::run_writer`).
    header: RunHeader,
    /// Steps taken so far, so the state is the one after this many.
    steps_taken: u64,
    /// Whether the frame of the current step has been handed out already.
    ///
    /// The schedule says which steps are saved; this says whether *this* visit
    /// to a saved step has been served, so a driver that calls
    /// [`RunLoop::take_frame`] twice between steps — which a chunked driver
    /// does at every chunk boundary — gets the frame once.
    frame_taken: bool,
    /// Frames handed out so far.
    frames_taken: u64,
}

impl RunLoop {
    /// Start `scenario`, under the name `description` its header will carry.
    ///
    /// Nothing is stepped: the returned loop sits at the initial state, whose
    /// frame is the run's first.
    ///
    /// # Errors
    /// [`RunLoopError::Grid`] if the scenario's grid is not one the output
    /// format can describe, and [`RunLoopError::Solver`] if the scenario asked
    /// for a timestep the scheme refuses.
    pub fn of_scenario(scenario: &Scenario, description: &str) -> Result<Self, RunLoopError> {
        let basin = scenario.basin();
        let grid = basin.grid();
        let params = scenario.physical_params();
        let schedule = scenario.output_schedule();

        // The beta-plane is placed by the basin's southern boundary rather
        // than centred on the grid, so a scenario that does not straddle the
        // equator gets the rotation of the latitudes it actually covers. It is
        // the same plane the scenario loader checked the rotation bound
        // against, because both ask the basin for it.
        let plane = BetaPlane::of_basin(params, basin);
        // The Epic 12 switch, and the whole of it: a scenario with an `[sst]`
        // section integrates the mixed-layer anomaly `T'` alongside the core,
        // one without integrates the three-variable model of Epics 01-07. The
        // two branches differ in one term of the right-hand side and one field
        // of the state; the loop below is the same either way (T-12.1).
        let solver = match scenario.sst_params() {
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
        // The header's variable list and the state's fields come from the same
        // switch, one line apart, so a run cannot promise a `T'` it does not
        // integrate or integrate one it does not promise.
        let header = if solver.couples_sst() {
            header.with_sst_anomaly()
        } else {
            header
        };

        let state = if solver.couples_sst() {
            OceanState::at_rest_with_sst_anomaly(grid)
        } else {
            OceanState::at_rest(grid)
        };
        // The run's forcing: the scenario's wind and the one field it is
        // sampled into, held here rather than inside the solver so that it
        // survives the whole run. A steady wind is therefore sampled once for
        // the run (T-10.5, `docs/performance-notes.md`), and the field a frame
        // records is the very field that stage of the integration read — the
        // same instant, and now literally the same buffer.
        let forcing = match scenario.wind_response_params() {
            Some(response) => RunForcing::Coupled(Box::new(CoupledWind::new(
                basin,
                scenario.wind(),
                SstWindResponse::new(basin, response),
            ))),
            None => RunForcing::Prescribed(WindForcing::new(basin, scenario.wind())),
        };

        Ok(Self {
            solver,
            state,
            forcing,
            schedule,
            header,
            steps_taken: 0,
            frame_taken: false,
            frames_taken: 0,
        })
    }

    /// What this run promises about itself: its grid, its parameters, and the
    /// number of frames it will produce.
    #[must_use]
    pub const fn header(&self) -> &RunHeader {
        &self.header
    }

    /// How long this run is and which of its steps are saved.
    #[must_use]
    pub const fn schedule(&self) -> OutputSchedule {
        self.schedule
    }

    /// Steps taken so far; the state is the one after this many.
    #[must_use]
    pub const fn steps_taken(&self) -> u64 {
        self.steps_taken
    }

    /// Frames handed out so far, of the [`RunHeader`]'s promised count.
    #[must_use]
    pub const fn frames_taken(&self) -> u64 {
        self.frames_taken
    }

    /// Model time the run has reached, in seconds.
    ///
    /// The product `step · dt` rather than an accumulator advanced once per
    /// step, so the time a frame records carries one rounding rather than one
    /// per step.
    #[must_use]
    pub fn model_time_s(&self) -> f64 {
        self.schedule.model_time_at_step(self.steps_taken)
    }

    /// Whether the run has taken every step it was asked for.
    #[must_use]
    pub const fn is_finished(&self) -> bool {
        self.steps_taken >= self.schedule.total_steps()
    }

    /// The frame of the step the run is on, if the schedule saves it and it
    /// has not been handed out yet.
    ///
    /// Called before each step and once more after the last one, which is what
    /// puts the initial state in the run as frame zero and the final state at
    /// the end of it.
    pub fn take_frame(&mut self) -> Option<SavedStep<'_>> {
        if self.frame_taken || !self.schedule.writes_at_step(self.steps_taken) {
            return None;
        }
        self.frame_taken = true;
        self.frames_taken += 1;
        let t_s = self.schedule.model_time_at_step(self.steps_taken);
        // Destructured so the state and the stress are borrowed from disjoint
        // fields: the stress is sampled through `&mut forcing` and the state
        // read beside it, and a frame needs both at once.
        let Self { state, forcing, .. } = self;
        Some(SavedStep {
            t_s,
            wind_stress: forcing.at(t_s, state),
            state,
        })
    }

    /// Take one step of the run, and say whether there was one left to take.
    ///
    /// `false` means the run has reached the end of its schedule and nothing
    /// was stepped; the frame of the final state is [`RunLoop::take_frame`]'s
    /// to give, and it is given before this returns `false`.
    pub fn take_step(&mut self) -> bool {
        if self.is_finished() {
            return false;
        }
        let t_s = self.schedule.model_time_at_step(self.steps_taken);
        self.solver
            .step_with_forcing(&mut self.state, t_s, &mut self.forcing);
        self.steps_taken += 1;
        self.frame_taken = false;
        true
    }
}

/// The forcing of one run: the scenario's prescribed winds, and — when the
/// `[sst]` section switched the Epic 12 coupling on — the atmospheric response
/// added to them.
///
/// A closed enum rather than a `Box<dyn StageForcing>` for the same reason
/// [`ScenarioWind`](crate::ScenarioWind) is one: a run has exactly these two
/// shapes, and which one it has is decided once, before the first step. The
/// coupled arm is boxed because it is much the larger of the two and only one
/// scenario in the format has it.
enum RunForcing {
    /// The three-variable model of Epics 01-07: whatever the `[[wind]]`
    /// entries prescribe, and nothing else.
    Prescribed(WindForcing<CompositeWind>),
    /// The coupled model of Epic 12: the same winds, plus the wind response to
    /// the SST anomaly of the stage being evaluated.
    Coupled(Box<CoupledWind<CompositeWind>>),
}

impl StageForcing for RunForcing {
    fn basin(&self) -> Basin {
        match self {
            Self::Prescribed(forcing) => forcing.basin(),
            Self::Coupled(forcing) => forcing.basin(),
        }
    }

    fn at(&mut self, t_s: f64, state: &OceanState) -> &WindStressField {
        match self {
            Self::Prescribed(forcing) => forcing.at(t_s),
            Self::Coupled(forcing) => forcing.at(t_s, state),
        }
    }
}
