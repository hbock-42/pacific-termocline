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
//! into one time step of the linear shallow-water core, and the [`RunWriter`]
//! that saves the result at a configurable output cadence.

// The acceptance criterion of T-02.1 is that every field states its unit; the
// lint is what keeps that true as the crate grows.
#![warn(missing_docs)]

pub mod coriolis;
pub mod integrator;
pub mod params;
pub mod run_writer;
pub mod shallow_water;
pub mod solver;
pub mod state;
pub mod wind;

pub use coriolis::{BetaPlane, BetaPlaneError, CoriolisTerm};
pub use integrator::{Rk4, StateVector};

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
pub use wind::WindStress;

/// Re-exported so binaries and the visualizer agree on one format version, and
/// on where a run's two files live (ADR-0004: the format crate is the one
/// place either is defined).
pub use termocline_format::{FORMAT_VERSION, FRAME_FILE_NAME, HEADER_FILE_NAME};

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
