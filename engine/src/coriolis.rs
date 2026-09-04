//! The equatorial beta-plane and the Coriolis terms of the momentum equations.
//!
//! The model rotates on a beta-plane: `f = β·y`, exactly zero at the equator
//! and linear in the distance from it (`CONTEXT.md`, *Beta-plane*). That is
//! the whole of the equatorial waveguide — it is what traps Kelvin and Rossby
//! waves near `y = 0` — so the sign change at the equator has to be exact
//! rather than nearly exact, which is why [`BetaPlane`] evaluates `β·y` from
//! the row's position rather than from a tabulated field.
//!
//! From the momentum equations of `docs/planning/01-scientific-model.md`,
//!
//! ```text
//! ∂u/∂t − f·v = …        ∂v/∂t + f·u = …
//! ```
//!
//! the two contributions [`CoriolisTerm`] adds to a tendency are `+f·v` in the
//! u equation and `−f·u` in the v equation. Neither product is collocated on
//! the Arakawa C-grid of [ADR-0003]: `u` lives on east/west faces and `v` on
//! north/south faces, so each has to be interpolated onto the other's points.
//! The centred four-point average of the four opposite-face values surrounding
//! a point lands exactly on that point, and is therefore second-order accurate
//! with no one-sided guess anywhere in the interior.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use std::fmt;

use termocline_grid::{Field2D, Grid, Staggering, U_STAGGERING, V_STAGGERING};
use termocline_numerics::Spacing;

use crate::params::PhysicalParams;
use crate::state::OceanState;

/// Weight each of the four surrounding points carries in the centred average
/// that moves a velocity from its own faces onto the other component's.
const FOUR_POINT_AVERAGE_WEIGHT: f64 = 0.25;

/// Tendency written on a face that lies on the closed basin's wall.
///
/// A stated contract, not a quiet substitution: the basin is closed on all
/// four sides (`docs/planning/01-scientific-model.md`), so the velocity normal
/// to a wall is zero and does not evolve. Writing zero rather than skipping
/// the point also means a reused tendency buffer cannot leak a previous
/// stage's value into a wall. Epic 04 owns the boundary conditions proper;
/// until then this matches what `termocline-numerics` writes on the same
/// faces.
const WALL_TENDENCY: f64 = 0.0;

/// Why a beta-plane could not be built.
#[derive(Debug, Clone, PartialEq)]
pub enum BetaPlaneError {
    /// The basin's southern edge was not a finite position.
    NotFinite {
        /// Name of the parameter, matching its accessor.
        parameter: &'static str,
        /// The value supplied, in metres.
        value_m: f64,
    },
}

impl fmt::Display for BetaPlaneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFinite { parameter, value_m } => {
                write!(f, "{parameter} is {value_m}; it must be a finite position")
            }
        }
    }
}

impl std::error::Error for BetaPlaneError {}

/// Where a basin sits on the equatorial beta-plane, and the `f = β·y` it
/// implies at each row of the C-grid.
///
/// Only the meridional geometry is here, because `f` depends on `y` alone:
/// how far south of the equator the basin's southern edge lies, and how tall
/// a cell is. Both are scenario input; `β` comes from
/// [`PhysicalParams`](crate::PhysicalParams).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BetaPlane {
    /// Meridional gradient of the Coriolis parameter `β`, in m⁻¹s⁻¹.
    beta_per_m_per_s: f64,
    /// Cell height, in metres.
    dy_m: f64,
    /// Position of the basin's southern boundary, in metres north of the
    /// equator — negative for a basin whose southern edge is in the southern
    /// hemisphere, which is the usual case.
    southern_edge_y_m: f64,
}

impl BetaPlane {
    /// The beta-plane of a basin whose southern boundary lies
    /// `southern_edge_y_m` metres north of the equator (so a negative value
    /// for a southern-hemisphere edge).
    ///
    /// # Errors
    /// [`BetaPlaneError::NotFinite`] if the edge is not a finite position.
    pub fn new(
        params: PhysicalParams,
        spacing: Spacing,
        southern_edge_y_m: f64,
    ) -> Result<Self, BetaPlaneError> {
        if !southern_edge_y_m.is_finite() {
            return Err(BetaPlaneError::NotFinite {
                parameter: "southern_edge_y_m",
                value_m: southern_edge_y_m,
            });
        }
        Ok(Self {
            beta_per_m_per_s: params.beta_per_m_per_s(),
            dy_m: spacing.dy_m(),
            southern_edge_y_m,
        })
    }

    /// The beta-plane of a basin straddling the equator symmetrically: half
    /// its meridional extent to the south, half to the north.
    ///
    /// The idealized Epic 02 configuration, and the one the equatorial-wave
    /// validation of Epic 07 assumes — the waveguide is centred on the
    /// equator, so a basin that is not would trap waves against a wall.
    #[must_use]
    pub fn centered_on_equator(params: PhysicalParams, spacing: Spacing, grid: Grid) -> Self {
        let half_extent_m = 0.5 * grid.ny() as f64 * spacing.dy_m();
        Self {
            beta_per_m_per_s: params.beta_per_m_per_s(),
            dy_m: spacing.dy_m(),
            southern_edge_y_m: -half_extent_m,
        }
    }

    /// Position of the basin's southern boundary, in metres north of the
    /// equator.
    #[must_use]
    pub const fn southern_edge_y_m(self) -> f64 {
        self.southern_edge_y_m
    }

    /// Meridional position of the row `j` of a field at `staggering`, in
    /// metres north of the equator.
    ///
    /// The half-cell offset that separates a cell-center row from a
    /// north/south-face row comes from [`Staggering::offset_in_cells`] rather
    /// than from a literal here, per CODING_STANDARDS.md § Scope guards: the
    /// grid knows about staggering, the physics does not.
    #[must_use]
    pub fn y_of_row_m(self, staggering: Staggering, j: usize) -> f64 {
        let (_, offset_in_cells) = staggering.offset_in_cells();
        self.southern_edge_y_m + (j as f64 + offset_in_cells) * self.dy_m
    }

    /// The Coriolis parameter `f = β·y`, in s⁻¹, on the row `j` of a field at
    /// `staggering`.
    ///
    /// Exactly zero on a row that lies on the equator, and equal and opposite
    /// on rows equidistant either side of it.
    #[must_use]
    pub fn coriolis_at_row_per_s(self, staggering: Staggering, j: usize) -> f64 {
        self.beta_per_m_per_s * self.y_of_row_m(staggering, j)
    }
}

/// The Coriolis contribution to the momentum equations over one basin.
///
/// Built once per run and applied at every right-hand-side evaluation, so it
/// allocates nothing: [`CoriolisTerm::add_to_tendency`] reads a state and adds
/// into a caller-owned tendency buffer (CODING_STANDARDS.md § Performance).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoriolisTerm {
    grid: Grid,
    plane: BetaPlane,
}

impl CoriolisTerm {
    /// The Coriolis term for `grid` on `plane`.
    #[must_use]
    pub const fn new(grid: Grid, plane: BetaPlane) -> Self {
        Self { grid, plane }
    }

    /// The beta-plane this term rotates on.
    #[must_use]
    pub const fn plane(&self) -> BetaPlane {
        self.plane
    }

    /// Add `+f·v` to the zonal tendency and `−f·u` to the meridional one.
    ///
    /// Both products are formed at the point they are needed: `v` is averaged
    /// over the four north/south faces around each east/west face, and `u`
    /// over the four east/west faces around each north/south face. The four
    /// walls of the closed basin are set to [`WALL_TENDENCY`] rather than
    /// interpolated, since a wall has faces on one side only.
    ///
    /// `tendency` is added to, not overwritten: the Coriolis term is one
    /// contribution among the pressure gradient, the wind stress and the
    /// Rayleigh damping. The thermocline tendency is untouched — rotation does
    /// not appear in the continuity equation.
    ///
    /// # Panics
    /// If either state's fields do not have the shapes this term's grid asks
    /// for. A mis-shaped buffer means the calling code is wrong, which is what
    /// panics are for (CODING_STANDARDS.md § Correctness and failure).
    pub fn add_to_tendency(&self, state: &OceanState, tendency: &mut OceanState) {
        self.check_grid("state", state);
        self.check_grid("tendency", tendency);
        self.add_zonal_tendency(state.v(), tendency.u_mut());
        self.add_meridional_tendency(state.u(), tendency.v_mut());
    }

    /// `∂u/∂t += +f·v̄`, with `v̄` the four north/south faces around each
    /// east/west face.
    fn add_zonal_tendency(&self, v: &Field2D<f64>, du: &mut Field2D<f64>) {
        let (nx, ny) = (self.grid.nx(), self.grid.ny());
        for j in 0..ny {
            let f_per_s = self.plane.coriolis_at_row_per_s(U_STAGGERING, j);
            *at_mut(du, 0, j) = WALL_TENDENCY;
            for i in 1..nx {
                // The east/west face (i, j) is flanked by the cells i−1 and i,
                // each of which carries a southern face at j and a northern
                // one at j+1.
                let mean_v_m_per_s = FOUR_POINT_AVERAGE_WEIGHT
                    * (at(v, i - 1, j) + at(v, i, j) + at(v, i - 1, j + 1) + at(v, i, j + 1));
                *at_mut(du, i, j) += f_per_s * mean_v_m_per_s;
            }
            *at_mut(du, nx, j) = WALL_TENDENCY;
        }
    }

    /// `∂v/∂t += −f·ū`, with `ū` the four east/west faces around each
    /// north/south face.
    fn add_meridional_tendency(&self, u: &Field2D<f64>, dv: &mut Field2D<f64>) {
        let (nx, ny) = (self.grid.nx(), self.grid.ny());
        for i in 0..nx {
            *at_mut(dv, i, 0) = WALL_TENDENCY;
            *at_mut(dv, i, ny) = WALL_TENDENCY;
        }
        for j in 1..ny {
            let f_per_s = self.plane.coriolis_at_row_per_s(V_STAGGERING, j);
            for i in 0..nx {
                // The north/south face (i, j) is flanked by the cells j−1 and
                // j, each of which carries a western face at i and an eastern
                // one at i+1.
                let mean_u_m_per_s = FOUR_POINT_AVERAGE_WEIGHT
                    * (at(u, i, j - 1) + at(u, i + 1, j - 1) + at(u, i, j) + at(u, i + 1, j));
                *at_mut(dv, i, j) -= f_per_s * mean_u_m_per_s;
            }
        }
    }

    /// Panic unless `state` covers the same basin as this term.
    fn check_grid(&self, role: &str, state: &OceanState) {
        assert!(
            state.grid() == self.grid,
            "{role} covers {:?}, but this Coriolis term is built for {:?}",
            state.grid(),
            self.grid
        );
    }
}

/// The value at `(i, j)`, which the grid check has already proved is present.
fn at(field: &Field2D<f64>, i: usize, j: usize) -> f64 {
    *field
        .get(i, j)
        .expect("grid checked on entry, so this point exists")
}

/// Mutable access to `(i, j)`, likewise already proved present.
fn at_mut(field: &mut Field2D<f64>, i: usize, j: usize) -> &mut f64 {
    field
        .get_mut(i, j)
        .expect("grid checked on entry, so this point exists")
}
