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
//! [`benchmark`] and [`profiling`] are the odd ones out: they compute nothing
//! the simulation needs. [`benchmark`] holds the workloads `benches/` measures
//! (`docs/benchmarks.md`), and [`profiling`] the instrument that says where a
//! timestep's time goes (`docs/performance-notes.md`). Both are modules of the
//! library rather than helpers local to a `benches/` or `examples/` target, so
//! that a test can assert on them — which for a measurement is the difference
//! between a number and a claim.

// The acceptance criterion of T-02.1 is that every field states its unit; the
// lint is what keeps that true as the crate grows.
#![warn(missing_docs)]

pub mod basin;
pub mod benchmark;
pub mod boundary;
pub mod coriolis;
pub mod forcing;
pub mod inspect;
pub mod integrator;
pub mod params;
pub mod profiling;
pub mod progress;
pub mod run;
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
pub use inspect::{inspect_run, render_header};

/// Re-exported so the `run` command, its tests and any other caller name one
/// scenario runner: [`run::run_scenario_file`] is the whole command — load,
/// run, write — and [`run::run_scenario`] the same run from a scenario already
/// in hand.
pub use run::{
    run_scenario, run_scenario_file, run_scenario_file_observed, run_scenario_observed, RunError,
    RunReport,
};

/// Re-exported so the `run` command and its tests name one progress reporter:
/// [`progress::RunObserver`] is what a run tells, and [`progress::RunProgress`]
/// the reporter that turns it into a progress line and structured logs.
pub use progress::{LogLevel, ProgressReport, ProgressStyle, RunObserver, RunProgress, Verbosity};

pub use run_writer::{OutputSchedule, OutputScheduleError, RunWriteError, RunWriter};

/// Re-exported so a scenario, the solver and its tests name one set of
/// physical parameters and one state type.
pub use params::{
    PhysicalParams, PhysicalParamsError, EQUATORIAL_BETA_PER_M_PER_S,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
};
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
