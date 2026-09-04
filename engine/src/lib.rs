//! The simulation core: pure computation, no rendering and no UI.
//!
//! The engine takes a scenario — grid, physical parameters, wind forcing, run
//! length — integrates the 1.5-layer reduced-gravity shallow-water equations
//! forward in time, and writes the resulting ocean state through
//! [`termocline_format`].
//!
//! The physics lands in Epics 01–04; so far the crate carries the time
//! integrator, the prognostic [`OceanState`] and the [`PhysicalParams`] the
//! equations are written in terms of.

// The acceptance criterion of T-02.1 is that every field states its unit; the
// lint is what keeps that true as the crate grows.
#![warn(missing_docs)]

pub mod integrator;
pub mod params;
pub mod state;

pub use integrator::{Rk4, StateVector};

/// Re-exported so a scenario, the solver and its tests name one set of
/// physical parameters and one state type.
pub use params::{
    PhysicalParams, PhysicalParamsError, EQUATORIAL_BETA_PER_M_PER_S,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
};
pub use state::OceanState;

/// Re-exported so binaries and the visualizer agree on one format version.
pub use termocline_format::FORMAT_VERSION;

/// Re-exported so the solver, its tests and the scenario loader all share one
/// definition of what a grid cell is and where each variable sits on it.
pub use termocline_grid::{Field2D, Grid, Staggering, H_STAGGERING, U_STAGGERING, V_STAGGERING};

/// Re-exported so a scenario, the solver and the CLI all reach the same CFL
/// bound: [`check_timestep`] is the runtime check a run passes before its
/// first step, and it refuses an unstable timestep rather than clamping it.
pub use termocline_numerics::{
    check_timestep, max_stable_dt, CflError, Spacing, WaveSpeed, CFL_SAFETY_FACTOR,
};

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
