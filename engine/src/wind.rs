//! The surface wind stress forcing the upper layer.
//!
//! [`WindStress`] is the *stub* Epic 02 is allowed to carry (see
//! `docs/planning/epics/EPIC-02-ocean-dynamics.md`, § Out of scope): two
//! staggered fields of stress in pascals, each sitting on the face of the
//! velocity component it accelerates, so the momentum equations pick up
//! `τ/(ρ₀·H)` without an interpolation. The real forcing — the `WindStress`
//! trait, the trade-wind and seasonal scenarios, the westerly wind burst — is
//! Epic 03's, and it replaces this type rather than extending it.
//!
//! Sign convention comes from `CONTEXT.md`: easterly trade-wind stress is
//! `τx < 0`.

use termocline_grid::{Field2D, Grid, U_STAGGERING, V_STAGGERING};

/// Stress of a calm ocean surface, in Pa.
const CALM: f64 = 0.0;

/// A surface wind stress field over one basin, in pascals.
///
/// `τx` sits on the east/west faces with the zonal current anomaly `u`, and
/// `τy` on the north/south faces with `v`.
#[derive(Debug, Clone, PartialEq)]
pub struct WindStress {
    /// Shape of the basin the two fields cover.
    grid: Grid,
    /// Zonal wind stress `τx`, in Pa, at east/west faces. Negative is
    /// easterly — the direction the trade winds blow.
    tau_x_pa: Field2D<f64>,
    /// Meridional wind stress `τy`, in Pa, at north/south faces.
    tau_y_pa: Field2D<f64>,
}

impl WindStress {
    /// No wind at all over `grid`: both components exactly zero.
    ///
    /// The unforced limit the wave tests of Epic 07 run in.
    #[must_use]
    pub fn calm(grid: Grid) -> Self {
        Self::uniform(grid, CALM, CALM)
    }

    /// A stress of `tau_x_pa` by `tau_y_pa` pascals, the same at every point.
    ///
    /// The constant test forcing Epic 02 is scoped to; a stress that varies
    /// with position and time arrives with the scenarios of Epic 03.
    #[must_use]
    pub fn uniform(grid: Grid, tau_x_pa: f64, tau_y_pa: f64) -> Self {
        Self {
            grid,
            tau_x_pa: grid.allocate(U_STAGGERING, tau_x_pa),
            tau_y_pa: grid.allocate(V_STAGGERING, tau_y_pa),
        }
    }

    /// Shape of the basin this stress covers.
    #[must_use]
    pub const fn grid(&self) -> Grid {
        self.grid
    }

    /// Zonal wind stress `τx`, in Pa, at east/west faces.
    #[must_use]
    pub const fn tau_x_pa(&self) -> &Field2D<f64> {
        &self.tau_x_pa
    }

    /// Meridional wind stress `τy`, in Pa, at north/south faces.
    #[must_use]
    pub const fn tau_y_pa(&self) -> &Field2D<f64> {
        &self.tau_y_pa
    }
}
