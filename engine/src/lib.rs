//! The simulation core: pure computation, no rendering and no UI.
//!
//! The engine takes a scenario — grid, physical parameters, wind forcing, run
//! length — integrates the 1.5-layer reduced-gravity shallow-water equations
//! forward in time, and writes the resulting ocean state through
//! [`termocline_format`].
//!
//! The physics lands in Epics 01–04; so far the crate carries the time
//! integrator, the prognostic [`OceanState`], the [`PhysicalParams`] the
//! equations are written in terms of, the [`Solver`] that puts them together
//! into one time step of the linear shallow-water core, the [`WindStress`]
//! forcing that drives it — steady, seasonal, or a burst stacked on either —
//! the [`Scenario`] that names all of those in one TOML file, and the
//! [`RunWriter`] that saves the result at a configurable output cadence; and,
//! on the CLI side, the scenario runner behind `run` and the `inspect` command
//! that reports a written run's header back to a terminal.
//!
//! A run is driven either way through [`RunLoop`], which holds everything a
//! scenario needs to take its next step and hands out the frames the schedule
//! saves. `run` is that loop driven to completion into a directory; a browser
//! is the same loop driven a chunk at a time so the tab stays live (ADR-0012).
//!
//! [`benchmark`], [`profiling`] and [`precision`] are the odd ones out: they
//! compute nothing the simulation needs. [`benchmark`] holds the workloads
//! `benches/` measures (`docs/benchmarks.md`), [`profiling`] the instrument
//! that says where a timestep's time goes, and [`precision`] the one that says
//! what storing the fields at a narrower width would do to the answer (both
//! `docs/performance-notes.md`). All three are modules of the library rather
//! than helpers local to a `benches/` or `examples/` target, so that a test
//! can assert on them — which for a measurement is the difference between a
//! number and a claim.
//!
//! # The `fs` feature
//!
//! [ADR-0012] puts this crate in the browser: on the web the visualizer links
//! the engine and computes runs itself, so the engine has to build for
//! `wasm32-unknown-unknown`, where there is no filesystem. What a filesystem
//! buys is behind the `fs` feature — [`Scenario::load`], [`RunWriter::create`],
//! [`run`], [`inspect`] and [`benchmark`] — and the browser does not enable it.
//! Nothing physical is: [`solver`], [`shallow_water`], [`coriolis`],
//! [`boundary`], [`forcing`], [`state`] and [`params`] are compiled either way,
//! so a browser run and a native run are the same computation, which is what
//! makes ADR-0012's reproducibility argument hold. The portable path is
//! [`Scenario::from_toml`] over TOML *text*, [`Solver::new`], and
//! [`Solver::step_with_forcing`]; `engine/tests/filesystem_free_api.rs` walks
//! it, and CI builds the crate for `wasm32` with the feature off so it cannot
//! rot.
//!
//! Two modules compile for the browser without being usable there:
//! [`progress`] and [`profiling`] read the clock, and `Instant::now` panics on
//! `wasm32-unknown-unknown`. They are native instruments — a progress bar and a
//! timing probe — rather than filesystem access, so the `fs` feature is the
//! wrong name to hide them behind; what keeps them out of a browser build is
//! that the browser has no reason to call them.
//!
//! The `termocline` binary's own dependencies are the separate `cli` feature,
//! because an argument parser is not a filesystem. Both are on by default, so
//! a native build and the CLI are exactly what they were.
//!
//! [ADR-0012]: ../../docs/planning/adr/0012-the-browser-runs-the-engine.md

// The acceptance criterion of T-02.1 is that every field states its unit; the
// lint is what keeps that true as the crate grows.
#![warn(missing_docs)]

pub mod basin;
#[cfg(feature = "fs")]
pub mod benchmark;
pub mod boundary;
pub mod coriolis;
pub mod forcing;
#[cfg(feature = "fs")]
pub mod inspect;
pub mod integrator;
pub mod params;
pub mod precision;
pub mod profiling;
pub mod progress;
#[cfg(feature = "fs")]
pub mod run;
pub mod run_loop;
pub mod run_writer;
pub mod scenario;
pub mod shallow_water;
pub mod solver;
pub mod sst;
pub mod state;
pub mod wind_response;

pub use coriolis::{BetaPlane, BetaPlaneError, CoriolisTerm};
pub use integrator::{Rk4, StateVector};

/// Re-exported so the `inspect` command and its tests name one rendering of a
/// run's header.
#[cfg(feature = "fs")]
pub use inspect::{inspect_run, render_header};

/// Re-exported so the `run` command, its tests and any other caller name one
/// scenario runner: [`run::run_scenario_file`] is the whole command — load,
/// run, write — and [`run::run_scenario`] the same run from a scenario already
/// in hand.
#[cfg(feature = "fs")]
pub use run::{
    run_scenario, run_scenario_file, run_scenario_file_observed, run_scenario_observed, RunError,
};

/// Re-exported so the `run` command and its tests name one progress reporter:
/// [`progress::RunObserver`] is what a run tells, and [`progress::RunProgress`]
/// the reporter that turns it into a progress line and structured logs.
pub use progress::{
    LogLevel, ProgressReport, ProgressStyle, RunObserver, RunProgress, RunReport, Verbosity,
};

/// Re-exported so the CLI and the visualizer name one run loop: [`RunLoop`]
/// is a run taken a step at a time, which is what the browser holds
/// (ADR-0012), and [`run_loop::SavedStep`] one of the timesteps it saves.
pub use run_loop::{RunLoop, RunLoopError, SavedStep};

pub use run_writer::{frame_of, OutputSchedule, OutputScheduleError, RunWriteError, RunWriter};

/// Re-exported so a scenario, the solver and its tests name one set of
/// physical parameters and one state type.
pub use params::{
    PhysicalParams, PhysicalParamsError, EQUATORIAL_BETA_PER_M_PER_S,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
};
pub use precision::{StorageWidth, FIELD_STORAGE};
pub use shallow_water::{shallow_water_rhs, ShallowWaterRhs};

pub use solver::{check_rotation_timestep, step, RotationLimitError, Solver, SolverError};
/// Re-exported so a scenario, the solver and its tests name one mixed-layer
/// coupling: the [`sst::SstParams`] a `[sst]` section builds and the
/// [`sst::SstTerm`] that adds `∂T'/∂t` to a tendency. The extension is
/// Epic 12's and switched on per scenario; `CONTEXT.md` keeps `T'` out of the
/// linear ocean core, and so does the engine.
pub use sst::{SstParams, SstParamsError, SstTerm, SurfaceLayer, DEFAULT_SURFACE_DRAG_PER_S};
pub use state::OceanState;

/// Re-exported so a scenario, the forcing and the rotation all name one
/// basin geometry.
pub use basin::{Basin, BasinBounds, BasinBoundsError, BasinError};

/// Re-exported so the solver and its tests name one boundary condition: the
/// closed basin's no-normal-flow walls.
pub use boundary::NoNormalFlow;

/// Re-exported so a scenario names one wind forcing: the [`forcing::WindStress`]
/// trait a scenario implements, the scenarios that implement it, the
/// [`forcing::CompositeWind`] that stacks them, the
/// [`forcing::TimeDependence`] each of them declares, and the
/// [`forcing::WindStressField`] the solver actually reads — held across a run
/// by the [`forcing::WindForcing`] a time loop steps with.
pub use forcing::{
    CompositeWind, SeasonalTradeWinds, StageForcing, SteadyTradeWinds, TimeDependence,
    WindBurstAnomaly, WindForcing, WindStress, WindStressError, WindStressField, TROPICAL_YEAR_S,
};

/// Re-exported so a coupled scenario names the wind half of the Bjerknes
/// feedback in one place: the [`wind_response::SstWindResponse`] the atmosphere
/// answers with, the [`wind_response::WindResponseParams`] a scenario tunes it
/// by, and the [`wind_response::CoupledWind`] that adds it to the prescribed
/// winds of a run.
pub use wind_response::{
    CoupledWind, SstWindResponse, WindResponseError, WindResponseParams,
    ATMOSPHERIC_GRAVITY_WAVE_SPEED_M_PER_S, DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M,
};

/// Re-exported so the CLI, the tests and the example files name one scenario
/// format: [`scenario::ScenarioConfig`] is the TOML record on disk, and
/// [`scenario::Scenario`] the validated result the engine runs.
pub use scenario::{Scenario, ScenarioConfig, ScenarioError, ScenarioWind};

/// Re-exported so binaries and the visualizer agree on one format version, on
/// where a run's two files live, and on the one reader that opens them
/// (ADR-0004: the format crate is the one place any of them is defined).
pub use termocline_format::{
    RunReadError, RunReader, FORMAT_VERSION, FRAME_FILE_NAME, HEADER_FILE_NAME,
    OLDEST_READABLE_FORMAT_VERSION,
};

/// Re-exported so the solver, its tests and the scenario loader all share one
/// definition of what a grid cell is and where each variable sits on it.
pub use termocline_grid::{Field2D, Grid, Staggering, H_STAGGERING, U_STAGGERING, V_STAGGERING};

/// Re-exported so a scenario, the solver and the CLI all reach the same CFL
/// bound: [`check_timestep`] is the runtime check a run passes before its
/// first step, and it refuses an unstable timestep rather than clamping it.
/// [`Spacing`] is the one cell-spacing type both the bound and the operators
/// are stated in.
pub use termocline_numerics::{
    check_timestep, max_stable_dt, CflError, Spacing, SpacingError, WaveSpeed, CFL_SAFETY_FACTOR,
};

/// Re-exported so the right-hand side differentiates with the same C-grid
/// operators the numerics crate defines, rather than open-coding stencils.
pub use termocline_numerics::CGridOperators;

#[cfg(test)]
mod tests {
    #[test]
    fn workspace_links_the_format_crate() {
        assert_eq!(crate::FORMAT_VERSION, termocline_format::FORMAT_VERSION);
        assert_eq!(
            crate::OLDEST_READABLE_FORMAT_VERSION,
            termocline_format::OLDEST_READABLE_FORMAT_VERSION
        );
    }

    #[test]
    fn workspace_links_the_numerics_crate() {
        // The engine takes the CFL bound from `termocline-numerics` rather
        // than keeping a second copy of the formula. What the check does with
        // an unsafe timestep is covered in `tests/cfl_timestep.rs`.
        assert_eq!(
            crate::CFL_SAFETY_FACTOR,
            termocline_numerics::CFL_SAFETY_FACTOR
        );
    }

    #[test]
    fn workspace_links_the_grid_crate() {
        // The engine takes its grid vocabulary from `termocline-grid` rather
        // than defining a second one, per CODING_STANDARDS.md § Scope guards.
        assert_eq!(crate::H_STAGGERING, termocline_grid::Staggering::CellCenter);
    }
}
