//! Where a basin sits, and therefore where each of its C-grid points is.
//!
//! The grid knows the *shape* of a basin in cells and the spacing knows how
//! wide a cell is; neither knows where the basin's southwest corner lies. That
//! origin is scenario input (`CONTEXT.md`, *Basin*), and anything stated as a
//! function of position — the wind forcing of Epic 03, the real Pacific
//! truncation of T-04.1 — needs it before it can be evaluated on a grid at all.
//!
//! [`row_position_m`] and [`column_position_m`] are the whole of the index-to-
//! metres arithmetic, and they are `pub(crate)` so that
//! [`BetaPlane`](crate::BetaPlane) — which carries the meridional half of the
//! same geometry, because `f = β·y` needs it and predates this module — reads
//! rows from one definition rather than a second copy of it. Folding the two
//! types into one is T-04.1's, the ticket that makes basin geometry a
//! configuration parameter; until then the shared arithmetic is what keeps
//! them from drifting.

use std::fmt;

use termocline_grid::{Grid, Staggering};
use termocline_numerics::Spacing;

/// Meridional position of the row `j` of a field at `staggering`, in metres
/// north of a boundary at `southern_edge_y_m`.
///
/// The half-cell offset that separates a cell-center row from a
/// north/south-face row comes from [`Staggering::offset_in_cells`] rather than
/// from a literal here, per CODING_STANDARDS.md § Scope guards: the grid knows
/// about staggering, the physics does not.
pub(crate) fn row_position_m(
    southern_edge_y_m: f64,
    dy_m: f64,
    staggering: Staggering,
    j: usize,
) -> f64 {
    let (_, offset_in_cells) = staggering.offset_in_cells();
    southern_edge_y_m + (j as f64 + offset_in_cells) * dy_m
}

/// Zonal position of the column `i` of a field at `staggering`, in metres east
/// of a boundary at `western_edge_x_m`. The zonal twin of [`row_position_m`].
pub(crate) fn column_position_m(
    western_edge_x_m: f64,
    dx_m: f64,
    staggering: Staggering,
    i: usize,
) -> f64 {
    let (offset_in_cells, _) = staggering.offset_in_cells();
    western_edge_x_m + (i as f64 + offset_in_cells) * dx_m
}

/// Why a basin could not be placed.
///
/// This describes invalid *scenario input* — a run asking for a basin at a
/// position that is not a position — so it is returned rather than panicked,
/// and it names the offending value (CODING_STANDARDS.md § Correctness and
/// failure).
#[derive(Debug, Clone, PartialEq)]
pub enum BasinError {
    /// One of the basin's two edges was not a finite position.
    NotFinite {
        /// Name of the parameter, matching its accessor.
        parameter: &'static str,
        /// The value supplied, in metres.
        value_m: f64,
    },
}

impl fmt::Display for BasinError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFinite { parameter, value_m } => {
                write!(f, "{parameter} is {value_m}; it must be a finite position")
            }
        }
    }
}

impl std::error::Error for BasinError {}

/// One basin's shape, cell size and position: everything needed to say where
/// a given C-grid point is in metres.
///
/// `y` is measured north from the equator, so that [`Basin::y_of_row_m`] and
/// [`BetaPlane::y_of_row_m`](crate::BetaPlane::y_of_row_m) agree on a basin
/// built from the same southern edge: the forcing and the rotation must not
/// disagree about which row is the equator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Basin {
    /// Shape of the basin, in cells.
    grid: Grid,
    /// Cell width and height, in metres.
    spacing: Spacing,
    /// Position of the basin's western boundary, in metres.
    western_edge_x_m: f64,
    /// Position of the basin's southern boundary, in metres north of the
    /// equator — negative for a basin whose southern edge is in the southern
    /// hemisphere, which is the usual case.
    southern_edge_y_m: f64,
}

impl Basin {
    /// A basin of `grid` cells at `spacing`, with its southwest corner at
    /// `(western_edge_x_m, southern_edge_y_m)`.
    ///
    /// # Errors
    /// [`BasinError::NotFinite`] if either edge is not a finite position.
    pub fn new(
        grid: Grid,
        spacing: Spacing,
        western_edge_x_m: f64,
        southern_edge_y_m: f64,
    ) -> Result<Self, BasinError> {
        check_finite("western_edge_x_m", western_edge_x_m)?;
        check_finite("southern_edge_y_m", southern_edge_y_m)?;
        Ok(Self {
            grid,
            spacing,
            western_edge_x_m,
            southern_edge_y_m,
        })
    }

    /// A basin straddling the equator symmetrically, with its western
    /// boundary at `x = 0`.
    ///
    /// The idealized configuration the Epic 02 and 03 tests run in, and the
    /// one [`BetaPlane::centered_on_equator`](crate::BetaPlane::centered_on_equator)
    /// assumes: the equatorial waveguide is centred on the basin, so no wave
    /// is trapped against a wall. The real Pacific's truncation arrives with
    /// T-04.1.
    #[must_use]
    pub fn centered_on_equator(grid: Grid, spacing: Spacing) -> Self {
        Self {
            grid,
            spacing,
            western_edge_x_m: 0.0,
            southern_edge_y_m: -(grid.ny() as f64 * spacing.dy_m()) / 2.0,
        }
    }

    /// Shape of this basin, in cells.
    #[must_use]
    pub const fn grid(self) -> Grid {
        self.grid
    }

    /// Cell width and height of this basin, in metres.
    #[must_use]
    pub const fn spacing(self) -> Spacing {
        self.spacing
    }

    /// Position of the basin's western boundary, in metres.
    #[must_use]
    pub const fn western_edge_x_m(self) -> f64 {
        self.western_edge_x_m
    }

    /// Position of the basin's southern boundary, in metres north of the
    /// equator.
    #[must_use]
    pub const fn southern_edge_y_m(self) -> f64 {
        self.southern_edge_y_m
    }

    /// Width of the basin, in metres — the `L` of the analytic wind-driven
    /// tilt.
    #[must_use]
    pub fn zonal_extent_m(self) -> f64 {
        self.grid.nx() as f64 * self.spacing.dx_m()
    }

    /// Zonal position of the column `i` of a field at `staggering`, in metres.
    #[must_use]
    pub fn x_of_column_m(self, staggering: Staggering, i: usize) -> f64 {
        column_position_m(self.western_edge_x_m, self.spacing.dx_m(), staggering, i)
    }

    /// Meridional position of the row `j` of a field at `staggering`, in
    /// metres north of the equator.
    #[must_use]
    pub fn y_of_row_m(self, staggering: Staggering, j: usize) -> f64 {
        row_position_m(self.southern_edge_y_m, self.spacing.dy_m(), staggering, j)
    }
}

fn check_finite(parameter: &'static str, value_m: f64) -> Result<(), BasinError> {
    if value_m.is_finite() {
        return Ok(());
    }
    Err(BasinError::NotFinite { parameter, value_m })
}
