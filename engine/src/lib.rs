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

#[cfg(test)]
mod tests {
    #[test]
    fn workspace_links_the_format_crate() {
        assert_eq!(crate::FORMAT_VERSION, termocline_format::FORMAT_VERSION);
    }

    #[test]
    fn workspace_links_the_grid_crate() {
        // The engine takes its grid vocabulary from `termocline-grid` rather
        // than defining a second one, per CODING_STANDARDS.md § Scope guards.
        assert_eq!(crate::H_STAGGERING, termocline_grid::Staggering::CellCenter);
    }
}
