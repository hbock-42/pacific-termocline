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
//! That move is `termocline-numerics`' job, not this module's — the physics
//! here calls [`CGridOperators::face_y_to_face_x`] and its twin and never does
//! neighbour arithmetic of its own.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use std::fmt;

use termocline_grid::{Field2D, Grid, Staggering, U_STAGGERING, V_STAGGERING};
use termocline_numerics::{CGridOperators, Spacing};

use crate::params::PhysicalParams;
use crate::state::OceanState;

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
    /// Delegated to the `basin` module, which owns the index-to-metres
    /// arithmetic, so that this plane and the [`Basin`](crate::Basin) a
    /// forcing is sampled over cannot disagree about which row is the equator.
    #[must_use]
    pub fn y_of_row_m(self, staggering: Staggering, j: usize) -> f64 {
        crate::basin::row_position_m(self.southern_edge_y_m, self.dy_m, staggering, j)
    }

    /// The largest `|f|` any row of `grid` carries, in s⁻¹.
    ///
    /// `f = β·y` is largest in magnitude at whichever of the basin's two
    /// meridional boundaries lies further from the equator, and those are
    /// north/south-face rows — the outermost `v` rows, `0` and `ny`. It is the
    /// fastest rotation the scheme has to resolve, which is what makes it a
    /// bound on the timestep (see [`Solver`](crate::Solver)).
    #[must_use]
    pub fn largest_coriolis_magnitude_per_s(self, grid: Grid) -> f64 {
        let southern = self.coriolis_at_row_per_s(V_STAGGERING, 0).abs();
        let northern = self.coriolis_at_row_per_s(V_STAGGERING, grid.ny()).abs();
        southern.max(northern)
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
/// allocates nothing per call: the two interpolation buffers are allocated
/// here and reused across steps (CODING_STANDARDS.md § Performance).
#[derive(Debug, Clone, PartialEq)]
pub struct CoriolisTerm {
    grid: Grid,
    plane: BetaPlane,
    operators: CGridOperators,
    /// `v` interpolated onto the east/west faces, where the u equation needs
    /// it. Reused every call; never read before it is written.
    v_on_u_faces: Field2D<f64>,
    /// `u` interpolated onto the north/south faces, where the v equation needs
    /// it.
    u_on_v_faces: Field2D<f64>,
}

impl CoriolisTerm {
    /// The Coriolis term for `grid` at `spacing`, on `plane`.
    #[must_use]
    pub fn new(grid: Grid, spacing: Spacing, plane: BetaPlane) -> Self {
        Self {
            grid,
            plane,
            operators: CGridOperators::new(grid, spacing),
            v_on_u_faces: grid.allocate(U_STAGGERING, 0.0),
            u_on_v_faces: grid.allocate(V_STAGGERING, 0.0),
        }
    }

    /// Add `+f·v` to the zonal tendency and `−f·u` to the meridional one.
    ///
    /// Each velocity is first moved onto the other component's faces by the
    /// C-grid four-point average, then multiplied by the `f` of the row it
    /// landed on — cell-center rows for the u equation, north/south-face rows
    /// for the v equation.
    ///
    /// `tendency` is added to, never overwritten, at every point including the
    /// closed basin's four walls: a wall carries no interpolated velocity (the
    /// operators leave it at zero, since it has cells on one side only), so it
    /// receives an exactly-zero contribution rather than an assignment that
    /// would discard what the pressure-gradient, wind-stress or damping terms
    /// had already written there. The thermocline tendency is untouched —
    /// rotation does not appear in the continuity equation.
    ///
    /// # Panics
    /// If either state's fields do not have the shapes this term's grid asks
    /// for. A mis-shaped buffer means the calling code is wrong, which is what
    /// panics are for (CODING_STANDARDS.md § Correctness and failure).
    pub fn add_to_tendency(&mut self, state: &OceanState, tendency: &mut OceanState) {
        self.check_grid("state", state);
        self.check_grid("tendency", tendency);

        self.operators
            .face_y_to_face_x(state.v(), &mut self.v_on_u_faces);
        self.operators
            .face_x_to_face_y(state.u(), &mut self.u_on_v_faces);

        accumulate_rows(tendency.u_mut(), &self.v_on_u_faces, |j| {
            self.plane.coriolis_at_row_per_s(U_STAGGERING, j)
        });
        accumulate_rows(tendency.v_mut(), &self.u_on_v_faces, |j| {
            -self.plane.coriolis_at_row_per_s(V_STAGGERING, j)
        });
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

/// `tendency += gain(j) · velocity` at every point, with `gain` evaluated once
/// per row because `f` depends on `y` alone.
///
/// Both fields are at the same staggering, so this is a pointwise loop: the
/// staggered neighbour arithmetic all happened in `termocline-numerics`.
fn accumulate_rows(
    tendency: &mut Field2D<f64>,
    velocity: &Field2D<f64>,
    gain_per_s: impl Fn(usize) -> f64,
) {
    let points_per_row = tendency.nx();
    for (j, row) in tendency
        .as_mut_slice()
        .chunks_exact_mut(points_per_row)
        .enumerate()
    {
        let gain = gain_per_s(j);
        let velocities = &velocity.as_slice()[j * points_per_row..][..points_per_row];
        for (rate, velocity_m_per_s) in row.iter_mut().zip(velocities) {
            *rate += gain * velocity_m_per_s;
        }
    }
}
