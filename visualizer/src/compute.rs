//! Computing a run in the tab that shows it.
//!
//! Per [ADR-0012] the web build has nothing to download: it holds a
//! [`Scenario`], steps it with the engine, and renders the frames it produces.
//! This module is that loop and the two limits it runs against.
//!
//! # The main thread must not block
//!
//! A run is 17 520 steps. Taking them in one call would freeze the tab for as
//! long as it takes — seconds at best — and a frozen tab shows no progress,
//! answers no clicks and, past a few seconds, is offered to the user for
//! killing. So [`ComputedRun::advance_within`] takes steps until a wall-clock
//! deadline and returns, and the shell calls it once per displayed frame. What
//! the reader sees is the run developing, which is the *watch the simulation
//! run* experience ADR-0001 traded away and ADR-0012 buys back.
//!
//! The deadline is checked every [`STEPS_PER_SLICE`] steps rather than every
//! step, so a clock read is amortised over a slice; the slice is what bounds
//! the overshoot past the deadline. Both constants say what they were measured
//! against.
//!
//! # Memory is the limit, and it is enforced here
//!
//! ADR-0012's other consequence: with nothing downloaded, what bounds a run is
//! the memory its frames occupy. A frame of the 0.5° control basin is 1.29 MB
//! and the control run is 731 of them — 941 MB, which no tab should be asked
//! to hold. [`FrameBudget`] is that limit written down: it is checked against
//! the header *before* the first step, so a scenario too big for a browser is
//! refused with an explanation rather than discovered when the tab dies.
//!
//! # The scale, and why the run is drawn before it is finished
//!
//! A run-wide colour scale needs every frame (T-08.2: it is what lets a
//! collapsing tilt be seen to collapse). A computed run does not have every
//! frame until it is over, so [`LoadedRun`] widens its scales as frames arrive
//! and reports [`LoadedRun::is_complete`] until then. The alternative —
//! compute everything, then draw — is never misleading but withholds exactly
//! what ADR-0012 bought. Widening is what the coordinator's own render tool
//! did, and what the shell shows is honest as long as it says the scale is
//! provisional, which is what [`crate::app`] does.
//!
//! [ADR-0012]: ../../docs/planning/adr/0012-the-browser-runs-the-engine.md

use core::fmt;

use engine::{RunLoop, RunLoopError, Scenario, ScenarioError};
use termocline_format::{GridSpec, RunHeader, Variable};
use web_time::{Duration, Instant};

use crate::run::FrameAppendError;
use crate::{LoadedRun, ScenarioPreset};

/// Steps taken between two reads of the clock.
///
/// The deadline of [`ComputedRun::advance_within`] can only be honoured to
/// within one slice, so this is how far past it a chunk may run. One step of
/// the presets' 80 × 25 grid takes **41 µs** natively (measured:
/// 17 520 steps of the control preset ([`ScenarioPreset::ALL`]) in 0.72 s, release
/// build, Apple M1 Pro — the same 0.7 s the halted T-08.4 work reported), so
/// eight steps is 0.33 ms there and 1.3 ms even at four times that cost, which
/// is what a wasm build is assumed to be until one is profiled. Both are
/// comfortably inside the budget below; a slice of one would pay a
/// `performance.now()` call every 41 µs to buy accuracy nothing can see.
pub const STEPS_PER_SLICE: u64 = 8;

/// How long a displayed frame may spend stepping the run.
///
/// Half of a 60 Hz frame's 16.7 ms, leaving the other half for what the tab
/// exists to do — decode the newest frame, colour-map it, and draw. At the 41
/// µs a step costs natively that is around 195 steps a frame, so the 17 520 of
/// the control preset take about 90 displayed frames, or 1.5 s. A browser is
/// slower than that by some factor nothing here has measured, and what the
/// factor changes is how long the run takes, not whether the tab keeps
/// drawing: that is what taking the budget in *time* buys.
pub const STEP_BUDGET: Duration = Duration::from_micros(8_000);

/// Bytes an `f64` occupies in an encoded frame.
///
/// The frame encoding is `bincode`'s fixed-width standard configuration
/// (`termocline_format::frame_encoding`), so every value of every field is
/// eight bytes wherever it sits.
const BYTES_PER_VALUE: u64 = 8;

/// The most memory one run's frames may occupy in a tab.
///
/// A limit rather than a hope, and per *run* rather than per tab: the shell
/// shows at most two (T-09.5), so the worst case is twice this. Thirty-two
/// mebibytes leaves the control preset — 244 frames of an 80 × 25 basin,
/// 19.9 MB — at 59 % of its budget, and the wind burst, which buys a daily
/// frame so the wave can be followed, at 89 % of it (`crate::presets`). That
/// is room for a longer or a more finely sampled scenario without room for a
/// careless one: the 941 MB control run is refused twenty-eight times over.
///
/// It is deliberately far below what a tab can technically allocate. What a
/// browser will let a page keep before it kills the tab is neither documented
/// nor constant across devices, and a limit tuned to the most generous desktop
/// is a limit that fails on a phone.
const BROWSER_FRAME_BUDGET_BYTES: u64 = 32 * 1024 * 1024;

/// How much of a run's frames a tab will hold.
///
/// The check is on the run's *header*, which states the grid and the frame
/// count, so it is answerable before a single step is taken — which is the
/// point. A scenario is refused for its size or it is run to the end; nothing
/// is half-computed and then abandoned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FrameBudget {
    /// The most bytes of encoded frames one run may retain.
    max_bytes: u64,
}

impl FrameBudget {
    /// The budget a browser tab is held to.
    #[must_use]
    pub const fn browser() -> Self {
        Self {
            max_bytes: BROWSER_FRAME_BUDGET_BYTES,
        }
    }

    /// The most bytes of encoded frames one run may retain.
    #[must_use]
    pub const fn max_bytes(&self) -> u64 {
        self.max_bytes
    }

    /// The frames `header` promises, as bytes of encoded frame data.
    ///
    /// The sum over the run's variables of one `f64` per point of that
    /// variable's staggered position, times the frames the header promises.
    /// The encoding adds a few bytes a frame on top — a length prefix per
    /// field and the frame's own model time — which for the control preset
    /// is 24 bytes in 81 704, or 0.03 %: far inside the headroom the budget
    /// leaves, and not worth modelling `bincode`'s varints in the visualizer
    /// to recover.
    #[must_use]
    pub fn bytes_of(header: &RunHeader) -> u64 {
        bytes_of_fields(
            header.grid,
            header.variables.iter().map(|spec| spec.variable),
            header.output.frame_count,
        )
    }

    /// The same, for a run that has not started: `frame_count` frames over
    /// `grid`, carrying the five variables of the linear core.
    ///
    /// What lets a [`ScenarioPreset`] state its size before it is pressed
    /// (`crate::presets`): a header exists only once a run has been started,
    /// and a visitor deserves the number before that. The variables are
    /// [`Variable::LINEAR_CORE`] because a preset is an uncoupled scenario;
    /// `tests/scenario_presets.rs` holds the estimate against
    /// [`FrameBudget::bytes_of`] of the header each preset actually produces.
    #[must_use]
    pub fn bytes_of_linear_core(grid: GridSpec, frame_count: u64) -> u64 {
        bytes_of_fields(grid, Variable::LINEAR_CORE.into_iter(), frame_count)
    }

    /// Whether a run with this header fits, and by how much it does not.
    ///
    /// # Errors
    /// [`BudgetExceeded`] naming what the run would cost and what it is
    /// allowed, so the refusal a reader sees is a number they can act on
    /// rather than "too big".
    pub fn admits(&self, header: &RunHeader) -> Result<(), BudgetExceeded> {
        let needed_bytes = Self::bytes_of(header);
        if needed_bytes <= self.max_bytes {
            return Ok(());
        }
        Err(BudgetExceeded {
            needed_bytes,
            budget_bytes: self.max_bytes,
            frame_count: header.output.frame_count,
            nx: header.grid.nx(),
            ny: header.grid.ny(),
        })
    }
}

/// Bytes of encoded frame data in `frame_count` frames over `grid` carrying
/// `variables`.
///
/// The one place the format's arithmetic is done — a sum over the variables of
/// one `f64` per point of that variable's staggered position — so the size a
/// run is admitted on and the size a preset states before it starts cannot be
/// two different counts.
fn bytes_of_fields(
    grid: GridSpec,
    variables: impl Iterator<Item = Variable>,
    frame_count: u64,
) -> u64 {
    let per_frame: u64 = variables
        .map(|variable| grid.field_len(variable) as u64 * BYTES_PER_VALUE)
        .sum();
    per_frame * frame_count
}

/// A scenario whose frames would not fit in the tab that was asked to compute
/// it.
///
/// Carries the run's shape as well as its size: "941 MB" says the scenario is
/// too big, and "731 frames of 320 × 100" says which of the three dials —
/// resolution, length, output cadence — to turn.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BudgetExceeded {
    /// Bytes of frames the run would retain.
    pub needed_bytes: u64,
    /// Bytes it is allowed.
    pub budget_bytes: u64,
    /// Frames the run promises.
    pub frame_count: u64,
    /// Cells across the basin.
    pub nx: usize,
    /// Cells up the basin.
    pub ny: usize,
}

impl fmt::Display for BudgetExceeded {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "this scenario's {frames} frames of {nx} × {ny} cells would hold {needed} in memory, \
             and a run in a browser is allowed {budget}; use a coarser grid, a shorter run, or a \
             longer output interval",
            frames = self.frame_count,
            nx = self.nx,
            ny = self.ny,
            needed = InMegabytes(self.needed_bytes),
            budget = InMegabytes(self.budget_bytes),
        )
    }
}

impl std::error::Error for BudgetExceeded {}

/// A size in bytes, written the way the UI states one.
///
/// Megabytes rather than mebibytes: the sizes this reports are memory a reader
/// compares against a scenario file's frame count, not against a page size.
pub struct InMegabytes(
    /// The size, in bytes. Megabytes are how it is *written*, not how it is
    /// held: rounding it before formatting would put the rounding in two
    /// places.
    pub u64,
);

impl fmt::Display for InMegabytes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        #[allow(clippy::cast_precision_loss)]
        let mb = self.0 as f64 / 1e6;
        write!(f, "{mb:.1} MB")
    }
}

/// Why a run could not be computed in the browser.
///
/// Each variant is something the *scenario* asked for that this build will not
/// do, so all of them are returned rather than panicked and each carries the
/// message of whatever objected (CODING_STANDARDS.md § *Correctness and
/// failure*).
#[derive(Debug)]
pub enum ComputeError {
    /// The scenario text is not one the engine will run.
    Scenario(ScenarioError),
    /// The engine refused to start the run — a grid the format cannot
    /// describe, or a timestep outside the CFL bound.
    Engine(RunLoopError),
    /// The run's frames would not fit in the tab.
    Budget(BudgetExceeded),
    /// A computed frame did not fit the header the engine wrote for it, which
    /// would mean the engine and the format disagree.
    Frame(FrameAppendError),
}

impl fmt::Display for ComputeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Scenario(error) => error.fmt(f),
            Self::Engine(error) => error.fmt(f),
            Self::Budget(error) => error.fmt(f),
            Self::Frame(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for ComputeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Scenario(error) => Some(error),
            Self::Engine(error) => Some(error),
            Self::Budget(error) => Some(error),
            Self::Frame(error) => Some(error),
        }
    }
}

impl From<ScenarioError> for ComputeError {
    fn from(error: ScenarioError) -> Self {
        Self::Scenario(error)
    }
}

impl From<RunLoopError> for ComputeError {
    fn from(error: RunLoopError) -> Self {
        Self::Engine(error)
    }
}

impl From<BudgetExceeded> for ComputeError {
    fn from(error: BudgetExceeded) -> Self {
        Self::Budget(error)
    }
}

impl From<FrameAppendError> for ComputeError {
    fn from(error: FrameAppendError) -> Self {
        Self::Frame(error)
    }
}

/// A run being computed in the tab, and the frames of it produced so far.
///
/// The two halves of the browser's answer to ADR-0012: the engine's
/// [`RunLoop`], stepped a slice at a time, and the [`LoadedRun`] its frames go
/// into — which is the same `LoadedRun` a file would have produced, so every
/// view reads it without knowing where it came from.
pub struct ComputedRun {
    /// The engine's loop, mid-run.
    stepping: RunLoop,
    /// The frames produced so far, as a run the views can draw.
    run: LoadedRun,
}

impl ComputedRun {
    /// Start computing `scenario`, under the name `description`, holding it
    /// to `budget`.
    ///
    /// Nothing is stepped: the returned run holds no frames yet, and the first
    /// call to [`ComputedRun::advance_within`] produces the initial state as
    /// frame zero.
    ///
    /// # Errors
    /// [`ComputeError::Engine`] if the engine refuses to start the run, and
    /// [`ComputeError::Budget`] if its frames would not fit in the tab — the
    /// last checked before a single step is taken.
    pub fn of_scenario(
        scenario: &Scenario,
        description: &str,
        budget: FrameBudget,
    ) -> Result<Self, ComputeError> {
        let stepping = RunLoop::of_scenario(scenario, description)?;
        budget.admits(stepping.header())?;
        let run = LoadedRun::computing(description, stepping.header().clone());
        Ok(Self { stepping, run })
    }

    /// Start computing what `preset` describes, under the preset's own name.
    ///
    /// The whole of what pressing a scenario button does (T-13.3), and the
    /// reason nothing carries over between two of them: a preset produces a
    /// *new* [`ComputedRun`], holding a new [`RunLoop`] over a new
    /// [`LoadedRun`], so the run that was being computed is dropped rather
    /// than reset.
    ///
    /// # Errors
    /// [`ComputeError::Scenario`] if the preset does not describe a scenario
    /// the engine will run, plus the errors of
    /// [`ComputedRun::of_scenario`].
    pub fn of_preset(preset: ScenarioPreset, budget: FrameBudget) -> Result<Self, ComputeError> {
        Self::of_scenario(&preset.scenario()?, preset.name(), budget)
    }

    /// The run as far as it has been computed.
    ///
    /// A [`LoadedRun`] like any other: it holds the frames produced so far,
    /// its scales cover those frames, and it says so through
    /// [`LoadedRun::is_complete`] until the last one has arrived.
    #[must_use]
    pub const fn run(&self) -> &LoadedRun {
        &self.run
    }

    /// Frames produced, of the frames the run promises.
    #[must_use]
    pub const fn progress(&self) -> (u64, u64) {
        (
            self.stepping.frames_taken(),
            self.stepping.header().output.frame_count,
        )
    }

    /// Whether every step of the run has been taken.
    #[must_use]
    pub const fn is_finished(&self) -> bool {
        self.stepping.is_finished()
    }

    /// Take up to `max_steps` steps, saving whichever frames fall inside them.
    ///
    /// The deterministic half of the loop, and the one the tests drive: where
    /// a chunk ends changes when the caller is handed control, never what the
    /// solver did with the state (`engine/tests/run_loop.rs`).
    ///
    /// # Errors
    /// [`ComputeError::Frame`] if a computed frame does not fit the header the
    /// engine wrote for it.
    pub fn advance_steps(&mut self, max_steps: u64) -> Result<(), ComputeError> {
        for _ in 0..max_steps {
            self.save_due_frame()?;
            if !self.stepping.take_step() {
                return Ok(());
            }
        }
        // The frame of the step the chunk ended on, so a run that stops
        // between two calls is never one frame behind what it has computed.
        self.save_due_frame()
    }

    /// Step until `budget` of wall-clock time has been spent, checking the
    /// clock every [`STEPS_PER_SLICE`] steps.
    ///
    /// What the shell calls once per displayed frame. It returns as soon as
    /// the budget is spent or the run is over, whichever comes first, so a
    /// finished run costs a comparison and nothing else.
    ///
    /// # Errors
    /// The errors of [`ComputedRun::advance_steps`].
    pub fn advance_within(&mut self, budget: Duration) -> Result<(), ComputeError> {
        let deadline = Instant::now() + budget;
        while !self.is_finished() {
            self.advance_steps(STEPS_PER_SLICE)?;
            if Instant::now() >= deadline {
                break;
            }
        }
        Ok(())
    }

    /// Save the frame of the step the run is on, if the schedule saves it and
    /// it has not been saved already.
    fn save_due_frame(&mut self) -> Result<(), ComputeError> {
        let grid = self.run.header().grid;
        let Some(saved) = self.stepping.take_frame() else {
            return Ok(());
        };
        let frame = saved
            .frame(&grid)
            .map_err(|error| ComputeError::Frame(FrameAppendError::Mismatch(error)))?;
        self.run.append_frame(&frame)?;
        Ok(())
    }
}
