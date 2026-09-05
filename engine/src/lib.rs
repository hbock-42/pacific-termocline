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
//! and the [`RunWriter`] that saves the result at a configurable output
//! cadence; and, on the CLI side, the `inspect` command that reports a written
//! run's header back to a terminal.

// The acceptance criterion of T-02.1 is that every field states its unit; the
// lint is what keeps that true as the crate grows.
#![warn(missing_docs)]

pub mod basin;
pub mod coriolis;
pub mod forcing;
pub mod inspect;
pub mod integrator;
pub mod params;
pub mod run_writer;
pub mod shallow_water;
pub mod solver;
pub mod state;

pub use coriolis::{BetaPlane, BetaPlaneError, CoriolisTerm};
pub use integrator::{Rk4, StateVector};

/// Re-exported so the `inspect` command and its tests name one rendering of a
/// run's header.
pub use inspect::{inspect_run, render_header};

pub use run_writer::{OutputSchedule, OutputScheduleError, RunWriteError, RunWriter};

/// Re-exported so a scenario, the solver and its tests name one set of
/// physical parameters and one state type.
pub use params::{
    PhysicalParams, PhysicalParamsError, EQUATORIAL_BETA_PER_M_PER_S,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
};
pub use shallow_water::{shallow_water_rhs, ShallowWaterRhs};
pub use solver::{step, Solver, SolverError};
pub use state::OceanState;

/// Re-exported so a scenario, the forcing and the rotation all name one
/// basin geometry.
pub use basin::{Basin, BasinError};

/// Re-exported so a scenario names one wind forcing: the [`forcing::WindStress`]
/// trait a scenario implements, the scenarios that implement it, the
/// [`forcing::WindStressField`] the solver actually reads, and the
/// [`forcing::CompositeWind`] that stacks an anomaly such as
/// [`forcing::WindBurstAnomaly`] on a base scenario.
pub use forcing::{
    CompositeWind, SeasonalTradeWinds, SteadyTradeWinds, WindBurstAnomaly, WindStress,
    WindStressError, WindStressField, TROPICAL_YEAR_S,
};

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
