//! Finite-difference operators on the Arakawa C-grid.
//!
//! `termocline-grid` says *where* each variable sits; this crate says how to
//! differentiate and interpolate between those positions. It is the metric
//! half of the numerical core: it is the first place a cell width in metres
//! appears, and it is still physics-free — nothing here knows that `h` is a
//! thermocline depth anomaly, only that `h` lives at cell centers.
//!
//! The staggering is the one fixed in [ADR-0003]: `h` at cell centers, `u` at
//! east/west faces, `v` at north/south faces. The two useful facts that
//! follow are that a difference of neighbouring center values lands exactly on
//! the face between them, and a difference of neighbouring face values lands
//! exactly on the center between them. Both are therefore centred differences
//! over one cell width and second-order accurate, with no averaging needed.
//!
//! CODING_STANDARDS.md tells solver code to reach for the named C-grid offsets
//! rather than raw `+1`/`-1` index arithmetic. This crate is the one place that
//! rule does not bite: the neighbour arithmetic here *is* the definition the
//! rule points at, and it is confined to the four private writers below so
//! that no caller has to repeat it.
//!
//! The same metric information answers a second question — how long a step the
//! scheme may take before it goes unstable — so the CFL bound lives here too,
//! in [`cfl`], next to the [`Spacing`] it is derived from.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

pub mod cfl;

pub use cfl::{
    check_timestep, max_stable_dt, CflError, WaveSpeed, CFL_SAFETY_FACTOR, RK4_IMAGINARY_AXIS_LIMIT,
};

use std::fmt;
use termocline_grid::{Axis, Field2D, Grid, Staggering};

/// Why a spacing could not be built.
#[derive(Debug, Clone, PartialEq)]
pub enum SpacingError {
    /// A cell width was zero, negative, or not a finite number.
    NotPositive {
        /// The axis whose spacing was rejected.
        axis: Axis,
        /// The value supplied, in metres.
        value_m: f64,
    },
}

impl fmt::Display for SpacingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotPositive { axis, value_m } => write!(
                f,
                "d{axis} is {value_m} m; cell spacing must be finite and greater than 0"
            ),
        }
    }
}

impl std::error::Error for SpacingError {}

/// Uniform cell spacing of a C-grid, in metres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spacing {
    dx_m: f64,
    dy_m: f64,
}

impl Spacing {
    /// A spacing of `dx_m` by `dy_m` metres per cell.
    ///
    /// # Errors
    /// [`SpacingError::NotPositive`] if either width is zero, negative or not
    /// finite.
    pub fn new(dx_m: f64, dy_m: f64) -> Result<Self, SpacingError> {
        check_positive(Axis::X, dx_m)?;
        check_positive(Axis::Y, dy_m)?;
        Ok(Self { dx_m, dy_m })
    }

    /// Cell width along x, in metres.
    #[must_use]
    pub const fn dx_m(self) -> f64 {
        self.dx_m
    }

    /// Cell height along y, in metres.
    #[must_use]
    pub const fn dy_m(self) -> f64 {
        self.dy_m
    }
}

/// The C-grid derivative and interpolation operators for one basin shape and
/// one cell spacing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CGridOperators {
    grid: Grid,
    spacing: Spacing,
}

impl CGridOperators {
    /// Operators for `grid` at `spacing`.
    #[must_use]
    pub const fn new(grid: Grid, spacing: Spacing) -> Self {
        Self { grid, spacing }
    }

    /// `∂/∂x` of a cell-centered field, written to an east/west-face field.
    ///
    /// The face between two cells is exactly the midpoint of their centers, so
    /// `(center[i] − center[i−1]) / dx` is a centred difference landing on the
    /// face: second-order accurate, with no averaging.
    ///
    /// The basin's western and eastern boundary faces have a cell on one side
    /// only and are set to zero rather than to a one-sided guess. They are the
    /// closed basin's walls, where boundary conditions (Epic 04) have the last
    /// word.
    ///
    /// # Panics
    /// If either field's shape does not match [`CGridOperators::grid`].
    pub fn ddx_center_to_face(&self, center: &Field2D<f64>, face: &mut Field2D<f64>) {
        self.check_shape("center input", center, Staggering::CellCenter);
        self.check_shape("east/west-face output", face, Staggering::EastWestFace);
        let inverse_dx = 1.0 / self.spacing.dx_m;
        self.write_x_faces(center, face, |west, east| (east - west) * inverse_dx);
    }

    /// `∂/∂y` of a cell-centered field, written to a north/south-face field.
    ///
    /// The meridional twin of [`CGridOperators::ddx_center_to_face`], including
    /// its treatment of the southern and northern boundary faces.
    ///
    /// # Panics
    /// If either field's shape does not match [`CGridOperators::grid`].
    pub fn ddy_center_to_face(&self, center: &Field2D<f64>, face: &mut Field2D<f64>) {
        self.check_shape("center input", center, Staggering::CellCenter);
        self.check_shape("north/south-face output", face, Staggering::NorthSouthFace);
        let inverse_dy = 1.0 / self.spacing.dy_m;
        self.write_y_faces(center, face, |south, north| (north - south) * inverse_dy);
    }

    /// `∂/∂x` of an east/west-face field, written to a cell-centered field.
    ///
    /// A cell center is exactly the midpoint of its two east/west faces, so
    /// `(face[i+1] − face[i]) / dx` is again a centred difference — and this
    /// direction is defined at every center, with no boundary gap.
    ///
    /// # Panics
    /// If either field's shape does not match [`CGridOperators::grid`].
    pub fn ddx_face_to_center(&self, face: &Field2D<f64>, center: &mut Field2D<f64>) {
        self.check_shape("east/west-face input", face, Staggering::EastWestFace);
        self.check_shape("center output", center, Staggering::CellCenter);
        let inverse_dx = 1.0 / self.spacing.dx_m;
        self.write_centers_from_x_faces(face, center, |west, east| (east - west) * inverse_dx);
    }

    /// `∂/∂y` of a north/south-face field, written to a cell-centered field.
    ///
    /// The meridional twin of [`CGridOperators::ddx_face_to_center`].
    ///
    /// # Panics
    /// If either field's shape does not match [`CGridOperators::grid`].
    pub fn ddy_face_to_center(&self, face: &Field2D<f64>, center: &mut Field2D<f64>) {
        self.check_shape("north/south-face input", face, Staggering::NorthSouthFace);
        self.check_shape("center output", center, Staggering::CellCenter);
        let inverse_dy = 1.0 / self.spacing.dy_m;
        self.write_centers_from_y_faces(face, center, |south, north| (north - south) * inverse_dy);
    }

    /// A cell-centered field interpolated onto the east/west faces.
    ///
    /// The two-point average of the neighbouring centers, which is the value
    /// at their midpoint to second order.
    ///
    /// The two boundary faces are **not** interpolated — they have a cell on
    /// one side only — and are set to zero. Read a boundary face
    /// of the output as "not computed here", never as an interpolated value.
    ///
    /// # Panics
    /// If either field's shape does not match [`CGridOperators::grid`].
    pub fn center_to_face_x(&self, center: &Field2D<f64>, face: &mut Field2D<f64>) {
        self.check_shape("center input", center, Staggering::CellCenter);
        self.check_shape("east/west-face output", face, Staggering::EastWestFace);
        self.write_x_faces(center, face, |west, east| AVERAGE_WEIGHT * (west + east));
    }

    /// A cell-centered field interpolated onto the north/south faces.
    ///
    /// The meridional twin of [`CGridOperators::center_to_face_x`], including
    /// its treatment of the two boundary faces.
    ///
    /// # Panics
    /// If either field's shape does not match [`CGridOperators::grid`].
    pub fn center_to_face_y(&self, center: &Field2D<f64>, face: &mut Field2D<f64>) {
        self.check_shape("center input", center, Staggering::CellCenter);
        self.check_shape("north/south-face output", face, Staggering::NorthSouthFace);
        self.write_y_faces(center, face, |south, north| {
            AVERAGE_WEIGHT * (south + north)
        });
    }

    /// An east/west-face field interpolated onto the cell centers.
    ///
    /// Defined at every center, since a center always has a face on each side.
    ///
    /// # Panics
    /// If either field's shape does not match [`CGridOperators::grid`].
    pub fn face_to_center_x(&self, face: &Field2D<f64>, center: &mut Field2D<f64>) {
        self.check_shape("east/west-face input", face, Staggering::EastWestFace);
        self.check_shape("center output", center, Staggering::CellCenter);
        self.write_centers_from_x_faces(face, center, |west, east| AVERAGE_WEIGHT * (west + east));
    }

    /// A north/south-face field interpolated onto the cell centers.
    ///
    /// # Panics
    /// If either field's shape does not match [`CGridOperators::grid`].
    pub fn face_to_center_y(&self, face: &Field2D<f64>, center: &mut Field2D<f64>) {
        self.check_shape("north/south-face input", face, Staggering::NorthSouthFace);
        self.check_shape("center output", center, Staggering::CellCenter);
        self.write_centers_from_y_faces(face, center, |south, north| {
            AVERAGE_WEIGHT * (south + north)
        });
    }

    /// Fill an east/west-face field from the pair of centers flanking each
    /// interior face, zeroing the two boundary lines.
    fn write_x_faces(
        &self,
        center: &Field2D<f64>,
        face: &mut Field2D<f64>,
        combine: impl Fn(f64, f64) -> f64,
    ) {
        let (nx, ny) = (self.grid.nx(), self.grid.ny());
        for j in 0..ny {
            *at_mut(face, 0, j) = BOUNDARY_FACE;
            for i in 1..nx {
                *at_mut(face, i, j) = combine(at(center, i - 1, j), at(center, i, j));
            }
            *at_mut(face, nx, j) = BOUNDARY_FACE;
        }
    }

    /// Fill a north/south-face field from the pair of centers flanking each
    /// interior face, zeroing the two boundary lines.
    fn write_y_faces(
        &self,
        center: &Field2D<f64>,
        face: &mut Field2D<f64>,
        combine: impl Fn(f64, f64) -> f64,
    ) {
        let (nx, ny) = (self.grid.nx(), self.grid.ny());
        for i in 0..nx {
            *at_mut(face, i, 0) = BOUNDARY_FACE;
            *at_mut(face, i, ny) = BOUNDARY_FACE;
        }
        for j in 1..ny {
            for i in 0..nx {
                *at_mut(face, i, j) = combine(at(center, i, j - 1), at(center, i, j));
            }
        }
    }

    /// Fill a cell-centered field from the pair of east/west faces flanking
    /// each center.
    fn write_centers_from_x_faces(
        &self,
        face: &Field2D<f64>,
        center: &mut Field2D<f64>,
        combine: impl Fn(f64, f64) -> f64,
    ) {
        for j in 0..self.grid.ny() {
            for i in 0..self.grid.nx() {
                *at_mut(center, i, j) = combine(at(face, i, j), at(face, i + 1, j));
            }
        }
    }

    /// Fill a cell-centered field from the pair of north/south faces flanking
    /// each center.
    fn write_centers_from_y_faces(
        &self,
        face: &Field2D<f64>,
        center: &mut Field2D<f64>,
        combine: impl Fn(f64, f64) -> f64,
    ) {
        for j in 0..self.grid.ny() {
            for i in 0..self.grid.nx() {
                *at_mut(center, i, j) = combine(at(face, i, j), at(face, i, j + 1));
            }
        }
    }

    /// Panic unless `field` has the shape this grid asks for at `staggering`.
    ///
    /// A mis-shaped buffer means the calling code is wrong, which is what
    /// panics are for (CODING_STANDARDS.md); a scenario cannot produce one.
    fn check_shape(&self, role: &str, field: &Field2D<f64>, staggering: Staggering) {
        let expected = self.grid.field_shape(staggering);
        let actual = (field.nx(), field.ny());
        assert!(
            actual == expected,
            "{role}: shape {actual:?} does not match the {expected:?} this grid asks for at \
             {staggering:?}"
        );
    }
}

/// Weight each of the two neighbours carries in a centred two-point average.
const AVERAGE_WEIGHT: f64 = 0.5;
/// What a center→face operator writes on the basin's boundary faces, which
/// have a cell on one side only.
///
/// This is a stated contract, not a quiet substitution: these operators define
/// the interior faces and reset the two boundary lines, so a reused buffer
/// cannot leak a previous step's values into a face nobody computed. Zero is
/// the value a closed wall carries until Epic 04 gives the boundary a
/// condition of its own, and
/// [`CGridOperators::center_to_face_x`] and its siblings say so on the way in.
/// It is deliberately not a NaN: an RK4 stage multiplies the wall by a normal
/// velocity that is itself zero there, and a NaN would poison the whole field.
const BOUNDARY_FACE: f64 = 0.0;

/// The value at `(i, j)`, which the shape check has already proved is present.
fn at(field: &Field2D<f64>, i: usize, j: usize) -> f64 {
    *field
        .get(i, j)
        .expect("shape checked on entry, so this point exists")
}

/// Mutable access to `(i, j)`, likewise already proved present.
fn at_mut(field: &mut Field2D<f64>, i: usize, j: usize) -> &mut f64 {
    field
        .get_mut(i, j)
        .expect("shape checked on entry, so this point exists")
}

fn check_positive(axis: Axis, value_m: f64) -> Result<(), SpacingError> {
    if !value_m.is_finite() || value_m <= 0.0 {
        return Err(SpacingError::NotPositive { axis, value_m });
    }
    Ok(())
}
