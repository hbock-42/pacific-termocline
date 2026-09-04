//! The prognostic state of the 1.5-layer reduced-gravity ocean.
//!
//! [`OceanState`] holds the three prognostic variables of the shallow-water
//! core — the thermocline depth anomaly `h` and the current anomalies `u` and
//! `v` — each allocated at its own position on the Arakawa C-grid fixed in
//! [ADR-0003]. The staggering comes from `termocline-grid`'s named constants
//! rather than from index arithmetic here, so this module never has to know
//! that a face field carries one extra line of points.
//!
//! The same type serves as a *tendency*: a right-hand-side evaluation (T-02.3)
//! produces an `OceanState` whose components are `∂h/∂t`, `∂u/∂t` and `∂v/∂t`
//! in SI-per-second, and [`StateVector`] is what lets [`crate::Rk4`] combine
//! those with a state.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use crate::integrator::StateVector;
use crate::params::PhysicalParams;
use termocline_grid::{Field2D, Grid, H_STAGGERING, U_STAGGERING, V_STAGGERING};

/// Value every prognostic variable takes in the undisturbed ocean.
///
/// All three are anomalies (`CONTEXT.md`), so the resting ocean is exactly
/// zero rather than some reference profile: at rest the thermocline sits at
/// its mean depth `H` and nothing is moving.
const AT_REST: f64 = 0.0;

/// The prognostic state of the upper layer over one basin.
///
/// The three fields are anomalies in SI units: `h` in metres, `u` and `v` in
/// m/s. `h` is a departure from the mean thermocline depth `H`, **not** a
/// total depth — the total is `H + h`, which
/// [`OceanState::total_thermocline_depth_m`] is there to spell out.
#[derive(Debug, Clone, PartialEq)]
pub struct OceanState {
    /// Shape of the basin the three fields cover.
    grid: Grid,
    /// Thermocline depth anomaly `h`, in metres, at cell centers.
    h: Field2D<f64>,
    /// Zonal current anomaly `u`, in m/s, at east/west faces.
    u: Field2D<f64>,
    /// Meridional current anomaly `v`, in m/s, at north/south faces.
    v: Field2D<f64>,
}

impl OceanState {
    /// The ocean at rest over `grid`: every anomaly exactly zero.
    ///
    /// This is the initial condition of every unforced run and the state a
    /// damped, unforced basin relaxes back to.
    #[must_use]
    pub fn at_rest(grid: Grid) -> Self {
        Self {
            grid,
            h: grid.allocate(H_STAGGERING, AT_REST),
            u: grid.allocate(U_STAGGERING, AT_REST),
            v: grid.allocate(V_STAGGERING, AT_REST),
        }
    }

    /// Shape of the basin this state covers.
    #[must_use]
    pub const fn grid(&self) -> Grid {
        self.grid
    }

    /// Thermocline depth anomaly `h`, in metres, at cell centers.
    #[must_use]
    pub const fn h(&self) -> &Field2D<f64> {
        &self.h
    }

    /// Zonal current anomaly `u`, in m/s, at east/west faces.
    #[must_use]
    pub const fn u(&self) -> &Field2D<f64> {
        &self.u
    }

    /// Meridional current anomaly `v`, in m/s, at north/south faces.
    #[must_use]
    pub const fn v(&self) -> &Field2D<f64> {
        &self.v
    }

    /// Mutable access to the thermocline depth anomaly `h`, in metres.
    pub fn h_mut(&mut self) -> &mut Field2D<f64> {
        &mut self.h
    }

    /// Mutable access to the zonal current anomaly `u`, in m/s.
    pub fn u_mut(&mut self) -> &mut Field2D<f64> {
        &mut self.u
    }

    /// Mutable access to the meridional current anomaly `v`, in m/s.
    pub fn v_mut(&mut self) -> &mut Field2D<f64> {
        &mut self.v
    }

    /// Total thermocline depth `H + h` at the cell center `(i, j)`, in metres,
    /// or `None` if that cell lies outside the basin.
    ///
    /// The one place the anomaly is turned into an absolute depth. Nothing in
    /// the solver needs this — the equations are written in `h` — but output
    /// and validation do, and doing the addition here keeps the distinction
    /// `CONTEXT.md` draws from being re-derived at each call site.
    #[must_use]
    pub fn total_thermocline_depth_m(
        &self,
        params: &PhysicalParams,
        i: usize,
        j: usize,
    ) -> Option<f64> {
        self.h
            .get(i, j)
            .map(|anomaly_m| params.mean_depth_m() + anomaly_m)
    }

    /// Panic unless `other` covers the same basin as `self`.
    ///
    /// An `OceanState`'s shape is not in its type, so [`StateVector`] requires
    /// a panic rather than a truncation here: two states over different basins
    /// mean the calling code is wrong, and silently combining the overlap
    /// would be exactly the silent clamping CODING_STANDARDS.md forbids.
    fn check_same_grid(&self, other: &Self) {
        assert!(
            self.grid == other.grid,
            "ocean states cover different basins: this grid is {:?}, the other is {:?}",
            self.grid,
            other.grid
        );
    }
}

impl StateVector for OceanState {
    fn assign(&mut self, source: &Self) {
        self.check_same_grid(source);
        for (target, values) in [
            (&mut self.h, source.h.as_slice()),
            (&mut self.u, source.u.as_slice()),
            (&mut self.v, source.v.as_slice()),
        ] {
            target.as_mut_slice().copy_from_slice(values);
        }
    }

    fn add_scaled(&mut self, factor: f64, other: &Self) {
        self.check_same_grid(other);
        for (target, values) in [
            (&mut self.h, other.h.as_slice()),
            (&mut self.u, other.u.as_slice()),
            (&mut self.v, other.v.as_slice()),
        ] {
            for (value, term) in target.as_mut_slice().iter_mut().zip(values) {
                *value += factor * term;
            }
        }
    }
}
