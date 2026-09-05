//! Where a basin sits, and therefore where each of its C-grid points is.
//!
//! Two types, one geometry, told twice because the two audiences count in
//! different units. [`BasinBounds`] is the scenario's truncation of the ocean
//! in *degrees* — "roughly 120°E–80°W by 25°S–25°N" (`CONTEXT.md`, *Basin*),
//! which is how a basin is named — and it is where every bound is checked.
//! [`Basin`] is the same truncation in *metres*, which is the only unit the
//! equations of `docs/planning/01-scientific-model.md` are written in; it is
//! what the forcing, the rotation and the solver read.
//!
//! The projection between them is one multiplication, [`METRES_PER_DEGREE_OF_ARC`],
//! because the model is an equatorial beta-plane: `f = β·y` is a linearization
//! about `φ = 0`, and on that plane a degree of longitude and a degree of
//! latitude are the same degree of arc. Anything more faithful — converging
//! meridians, an ellipsoid — would place the grid on a geometry the equations
//! are not solved on.
//!
//! [`row_position_m`] and [`column_position_m`] are the whole of the index-to-
//! metres arithmetic, and they are `pub(crate)` so that
//! [`BetaPlane`](crate::BetaPlane) — which carries the meridional half of the
//! same geometry, because `f = β·y` needs it and predates this module — reads
//! rows from one definition rather than a second copy of it. The two types
//! still want folding into one; until then the shared arithmetic is what keeps
//! them from drifting.

use std::fmt;

use termocline_grid::{Axis, Grid, Staggering};
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
    /// is trapped against a wall. A scenario file names a real truncation
    /// instead ([`BasinBounds`]), which is symmetric about the equator only if
    /// its two latitudes are; this constructor is for the tests that want the
    /// symmetry without stating a geography.
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

    /// Height of the basin, in metres — the meridional twin of
    /// [`Basin::zonal_extent_m`].
    #[must_use]
    pub fn meridional_extent_m(self) -> f64 {
        self.grid.ny() as f64 * self.spacing.dy_m()
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

/// Earth's mean radius `R`, in metres.
///
/// The IUGG arithmetic mean radius `R₁ = (2a + b)/3` of the WGS-84 ellipsoid.
/// It is the radius the model's `β = 2Ω·cos(φ)/R` is quoted from
/// ([`EQUATORIAL_BETA_PER_M_PER_S`](crate::EQUATORIAL_BETA_PER_M_PER_S)), so
/// the geometry and the rotation describe the same planet.
pub const EARTH_MEAN_RADIUS_M: f64 = 6_371_008.8;

/// Metres per degree of arc at Earth's mean radius: `R·π/180`.
///
/// The model is an equatorial beta-plane — a linearization about `φ = 0`
/// (`CONTEXT.md`, *Beta-plane*) — so a degree of longitude and a degree of
/// latitude are the same distance, and the projection from degrees to metres
/// is this one multiplication. The `cos(φ)` convergence of the meridians is
/// exactly the term the beta-plane approximation drops; reintroducing it here
/// would place the grid on a geometry the equations are not solved on.
pub const METRES_PER_DEGREE_OF_ARC: f64 = EARTH_MEAN_RADIUS_M * std::f64::consts::PI / 180.0;

/// Western boundary of the default basin, in degrees east: the Maritime
/// Continent edge of the Pacific (`CONTEXT.md`, *Basin*).
pub const PACIFIC_WESTERN_LONGITUDE_DEG: f64 = 120.0;
/// Eastern boundary of the default basin, in degrees east: 80°W, the South
/// American coast (`CONTEXT.md`, *Basin*).
pub const PACIFIC_EASTERN_LONGITUDE_DEG: f64 = -80.0;
/// Southern boundary of the default basin, in degrees north (`CONTEXT.md`,
/// *Basin*).
pub const PACIFIC_SOUTHERN_LATITUDE_DEG: f64 = -25.0;
/// Northern boundary of the default basin, in degrees north (`CONTEXT.md`,
/// *Basin*).
pub const PACIFIC_NORTHERN_LATITUDE_DEG: f64 = 25.0;
/// Cell size of the default basin, in degrees.
///
/// Half a degree is ≈ 55.6 km, so the equatorial deformation radius
/// `Le = √(c/β) ≈ 361 km` (`CONTEXT.md`) spans about six and a half cells:
/// enough to resolve the meridional structure of the waveguide the whole
/// model is about, at 320 × 100 cells for the Pacific bounds above.
pub const PACIFIC_RESOLUTION_DEG: f64 = 0.5;

/// A full turn of longitude, in degrees. Named because it is the modulus the
/// zonal span is measured in, not a magic number.
const FULL_TURN_DEG: f64 = 360.0;

/// The pole, in degrees of latitude: the bound a latitude may not exceed.
const POLE_LATITUDE_DEG: f64 = 90.0;

/// The most cells an axis may hold: the largest count a `usize` index can
/// address, as the `f64` the cell count is computed as. Beyond it the cast to
/// `usize` would saturate, turning an absurd resolution into a plausible-
/// looking grid rather than into an error.
const MAX_CELLS_PER_AXIS: f64 = usize::MAX as f64;

/// How far a cell count may sit from a whole number and still be one, as a
/// fraction of the count itself.
///
/// Bounds are written in decimal degrees, which are not exact in binary, so
/// `span/resolution` for a basin that *is* a whole number of cells lands a few
/// ulp off it: a relative error of order the machine epsilon, 2.2e-16, however
/// many cells there are. This bound is seven orders of magnitude looser than
/// that slack — so a basin stated in round degrees is never refused, at any
/// resolution — and still far tighter than any mis-specification worth
/// catching: 1e-9 of a cell of half a degree is 56 µm. Relative rather than
/// absolute precisely so the two ends of that sentence stay true together as
/// the cell count grows.
const WHOLE_CELL_TOLERANCE: f64 = 1e-9;

/// Why a set of basin bounds is not a basin.
///
/// Bounds are scenario input, so every variant is returned rather than
/// panicked and names the offending value and the bound it violated
/// (CODING_STANDARDS.md § *Correctness and failure*).
#[derive(Debug, Clone, PartialEq)]
pub enum BasinBoundsError {
    /// A boundary or the resolution was not a finite number of degrees.
    NotFinite {
        /// Name of the parameter, matching its accessor.
        parameter: &'static str,
        /// The value supplied, in degrees.
        value_deg: f64,
    },
    /// A latitude was off the planet.
    LatitudeOffThePlanet {
        /// Name of the parameter, matching its accessor.
        parameter: &'static str,
        /// The value supplied, in degrees north.
        value_deg: f64,
    },
    /// The northern boundary was not north of the southern one.
    LatitudesNotOrdered {
        /// The southern boundary supplied, in degrees north.
        southern_latitude_deg: f64,
        /// The northern boundary supplied, in degrees north.
        northern_latitude_deg: f64,
    },
    /// The resolution was not a positive number of degrees.
    ResolutionNotPositive {
        /// The value supplied, in degrees.
        value_deg: f64,
    },
    /// An axis of the basin was shorter than one cell — including the zero
    /// span of two boundaries at the same place.
    AxisShorterThanACell {
        /// The axis that is too short.
        axis: Axis,
        /// The span supplied, in degrees.
        span_deg: f64,
        /// The resolution supplied, in degrees.
        resolution_deg: f64,
    },
    /// An axis of the basin held more cells than a machine can index — a
    /// resolution so fine that the count does not fit in a `usize`.
    MoreCellsThanCanBeIndexed {
        /// The axis with too many cells.
        axis: Axis,
        /// The span supplied, in degrees.
        span_deg: f64,
        /// The resolution supplied, in degrees.
        resolution_deg: f64,
        /// The number of cells the two of them ask for.
        cells: f64,
    },
    /// An axis of the basin was not a whole number of cells. Refused rather
    /// than rounded: rounding it would silently run a basin nobody asked for
    /// (CODING_STANDARDS.md § *No silent clamping*).
    SpanNotAWholeNumberOfCells {
        /// The axis that does not divide.
        axis: Axis,
        /// The span supplied, in degrees.
        span_deg: f64,
        /// The resolution supplied, in degrees.
        resolution_deg: f64,
    },
}

impl fmt::Display for BasinBoundsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFinite {
                parameter,
                value_deg,
            } => write!(f, "{parameter} is {value_deg}; it must be a finite number of degrees"),
            Self::LatitudeOffThePlanet {
                parameter,
                value_deg,
            } => write!(
                f,
                "{parameter} is {value_deg}; it must be between {} and {POLE_LATITUDE_DEG} degrees north",
                -POLE_LATITUDE_DEG
            ),
            Self::LatitudesNotOrdered {
                southern_latitude_deg,
                northern_latitude_deg,
            } => write!(
                f,
                "northern_latitude_deg is {northern_latitude_deg}, which is not north of \
                 southern_latitude_deg {southern_latitude_deg}"
            ),
            Self::ResolutionNotPositive { value_deg } => write!(
                f,
                "resolution_deg is {value_deg}; it must be finite and greater than 0"
            ),
            Self::MoreCellsThanCanBeIndexed {
                axis,
                span_deg,
                resolution_deg,
                cells,
            } => write!(
                f,
                "the basin spans {span_deg} degrees of {} in cells of resolution_deg \
                 {resolution_deg}, which is {cells} cells: more than can be indexed",
                degrees_of(*axis)
            ),
            Self::AxisShorterThanACell {
                axis,
                span_deg,
                resolution_deg,
            } => write!(
                f,
                "the basin spans {span_deg} degrees of {}, which is less than one cell of \
                 resolution_deg {resolution_deg}",
                degrees_of(*axis)
            ),
            Self::SpanNotAWholeNumberOfCells {
                axis,
                span_deg,
                resolution_deg,
            } => write!(
                f,
                "the basin spans {span_deg} degrees of {}, which is not a whole number of cells \
                 of resolution_deg {resolution_deg}",
                degrees_of(*axis)
            ),
        }
    }
}

impl std::error::Error for BasinBoundsError {}

/// What a degree along `axis` is a degree of, for an error message: the file
/// states the zonal axis in longitude and the meridional one in latitude, and
/// the message has to use the word the file does.
fn degrees_of(axis: Axis) -> &'static str {
    match axis {
        Axis::X => "longitude",
        Axis::Y => "latitude",
    }
}

/// The truncation of the ocean a scenario runs on, in degrees: two meridians,
/// two parallels and a cell size.
///
/// This is the human-facing half of [`Basin`], and the reason it exists is
/// that the Pacific is stated in degrees — "roughly 120°E–80°W by 25°S–25°N"
/// (`CONTEXT.md`, *Basin*) — while the equations are solved in metres. Every
/// bound is checked here, once, so that [`BasinBounds::basin`] is total: a
/// value of this type is by construction a whole number of cells of a positive
/// size between two boundaries that are the right way round.
///
/// Longitude is counted eastward and wraps, so a basin may cross the
/// dateline — the Pacific has to. `-80.0` and `280.0` name the same eastern
/// boundary, and the span is always measured east from the western boundary:
/// two equal longitudes are a basin of zero width, not a basin around the
/// whole planet.
///
/// The output format states the same four degrees, as
/// `termocline_format::BasinExtent`, so that a written run records the stretch
/// of ocean it covers. That one is the *record* — a plain wire struct in the
/// crate ADR-0004 reserves for the file format, with no validation of its own;
/// this one is the *input*, and the place the bounds are checked. Neither can
/// stand in for the other without moving simulation logic into the format
/// crate, which ADR-0004 forbids.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BasinBounds {
    western_longitude_deg: f64,
    eastern_longitude_deg: f64,
    southern_latitude_deg: f64,
    northern_latitude_deg: f64,
    resolution_deg: f64,
}

impl BasinBounds {
    /// The bounds of a basin from its four boundaries and its cell size, in
    /// degrees.
    ///
    /// # Errors
    /// A [`BasinBoundsError`] naming the boundary that is not one: a value
    /// that is not finite, a latitude off the planet or south of the basin's
    /// own southern edge, a resolution that is not a positive size, or a span
    /// that is not a whole number of cells of it.
    pub fn new(
        western_longitude_deg: f64,
        eastern_longitude_deg: f64,
        southern_latitude_deg: f64,
        northern_latitude_deg: f64,
        resolution_deg: f64,
    ) -> Result<Self, BasinBoundsError> {
        check_finite_deg("western_longitude_deg", western_longitude_deg)?;
        check_finite_deg("eastern_longitude_deg", eastern_longitude_deg)?;
        check_latitude("southern_latitude_deg", southern_latitude_deg)?;
        check_latitude("northern_latitude_deg", northern_latitude_deg)?;
        if !(resolution_deg.is_finite() && resolution_deg > 0.0) {
            return Err(BasinBoundsError::ResolutionNotPositive {
                value_deg: resolution_deg,
            });
        }
        if northern_latitude_deg <= southern_latitude_deg {
            return Err(BasinBoundsError::LatitudesNotOrdered {
                southern_latitude_deg,
                northern_latitude_deg,
            });
        }

        let bounds = Self {
            western_longitude_deg,
            eastern_longitude_deg,
            southern_latitude_deg,
            northern_latitude_deg,
            resolution_deg,
        };
        // Both counts are computed and discarded here purely to reject the
        // spans that have no count; `nx` and `ny` recompute them once the
        // bounds are known to have one.
        cells_along(Axis::X, bounds.zonal_span_deg(), resolution_deg)?;
        cells_along(Axis::Y, bounds.meridional_span_deg(), resolution_deg)?;
        Ok(bounds)
    }

    /// The default basin: the equatorial Pacific of `CONTEXT.md`, 120°E–80°W
    /// by 25°S–25°N at [`PACIFIC_RESOLUTION_DEG`].
    ///
    /// # Panics
    /// Never: the constants are a valid basin, and a build that breaks that is
    /// a broken build rather than bad input.
    #[must_use]
    pub fn pacific() -> Self {
        Self::new(
            PACIFIC_WESTERN_LONGITUDE_DEG,
            PACIFIC_EASTERN_LONGITUDE_DEG,
            PACIFIC_SOUTHERN_LATITUDE_DEG,
            PACIFIC_NORTHERN_LATITUDE_DEG,
            PACIFIC_RESOLUTION_DEG,
        )
        .expect("the Pacific constants are a valid basin")
    }

    /// Western boundary, in degrees east.
    #[must_use]
    pub const fn western_longitude_deg(self) -> f64 {
        self.western_longitude_deg
    }

    /// Eastern boundary, in degrees east, as it was written.
    #[must_use]
    pub const fn eastern_longitude_deg(self) -> f64 {
        self.eastern_longitude_deg
    }

    /// Southern boundary, in degrees north.
    #[must_use]
    pub const fn southern_latitude_deg(self) -> f64 {
        self.southern_latitude_deg
    }

    /// Northern boundary, in degrees north.
    #[must_use]
    pub const fn northern_latitude_deg(self) -> f64 {
        self.northern_latitude_deg
    }

    /// Cell size, in degrees.
    #[must_use]
    pub const fn resolution_deg(self) -> f64 {
        self.resolution_deg
    }

    /// Longitude spanned by the basin, in degrees, counted eastward from the
    /// western boundary — so a basin crossing the dateline has the span one
    /// would draw on a map rather than a negative one.
    #[must_use]
    pub fn zonal_span_deg(self) -> f64 {
        (self.eastern_longitude_deg - self.western_longitude_deg).rem_euclid(FULL_TURN_DEG)
    }

    /// Latitude spanned by the basin, in degrees.
    #[must_use]
    pub fn meridional_span_deg(self) -> f64 {
        self.northern_latitude_deg - self.southern_latitude_deg
    }

    /// The [`Basin`] these bounds describe: the same truncation in metres, on
    /// the equatorial beta-plane the equations are solved on.
    ///
    /// `x` is measured east from the western boundary and `y` north from the
    /// equator, which is what lets [`Basin::y_of_row_m`] and
    /// [`BetaPlane::y_of_row_m`](crate::BetaPlane::y_of_row_m) agree about
    /// which row has `f = 0`.
    ///
    /// # Panics
    /// Never: [`BasinBounds::new`] has already refused every input a grid, a
    /// spacing or a basin would refuse, so a panic here is a broken invariant
    /// rather than bad input.
    #[must_use]
    pub fn basin(self) -> Basin {
        let grid =
            Grid::new(self.nx(), self.ny()).expect("validated bounds have cells on both axes");
        let cell_size_m = self.resolution_deg * METRES_PER_DEGREE_OF_ARC;
        let spacing = Spacing::new(cell_size_m, cell_size_m)
            .expect("a positive resolution is a positive cell");
        Basin::new(
            grid,
            spacing,
            WESTERN_BOUNDARY_X_M,
            self.southern_latitude_deg * METRES_PER_DEGREE_OF_ARC,
        )
        .expect("finite bounds place the basin at a finite position")
    }

    /// Cells east–west.
    ///
    /// # Panics
    /// Never, for the reason [`BasinBounds::basin`] gives.
    #[must_use]
    pub fn nx(self) -> usize {
        cells_along(Axis::X, self.zonal_span_deg(), self.resolution_deg)
            .expect("validated bounds are a whole number of cells")
    }

    /// Cells north–south.
    ///
    /// # Panics
    /// Never, for the reason [`BasinBounds::basin`] gives.
    #[must_use]
    pub fn ny(self) -> usize {
        cells_along(Axis::Y, self.meridional_span_deg(), self.resolution_deg)
            .expect("validated bounds are a whole number of cells")
    }
}

impl Default for BasinBounds {
    fn default() -> Self {
        Self::pacific()
    }
}

/// Position of the western boundary of a basin built from bounds, in metres:
/// `x` is measured east from that wall, so it is the origin. Named rather than
/// written as a bare `0.0` because it is the convention every `x_m` in a
/// scenario — a wind burst's centre, a diagnostic's station — is stated in.
const WESTERN_BOUNDARY_X_M: f64 = 0.0;

/// The whole number of cells of `resolution_deg` in `span_deg`.
///
/// # Errors
/// [`BasinBoundsError::AxisShorterThanACell`] if there is no cell in the span
/// at all, and [`BasinBoundsError::SpanNotAWholeNumberOfCells`] if the count is
/// not whole to within [`WHOLE_CELL_TOLERANCE`].
fn cells_along(axis: Axis, span_deg: f64, resolution_deg: f64) -> Result<usize, BasinBoundsError> {
    let cells = span_deg / resolution_deg;
    let whole = cells.round();
    if (cells - whole).abs() > WHOLE_CELL_TOLERANCE * whole.max(1.0) {
        return Err(BasinBoundsError::SpanNotAWholeNumberOfCells {
            axis,
            span_deg,
            resolution_deg,
        });
    }
    if whole < 1.0 {
        return Err(BasinBoundsError::AxisShorterThanACell {
            axis,
            span_deg,
            resolution_deg,
        });
    }
    if whole > MAX_CELLS_PER_AXIS {
        return Err(BasinBoundsError::MoreCellsThanCanBeIndexed {
            axis,
            span_deg,
            resolution_deg,
            cells: whole,
        });
    }
    // The two bounds above put `whole` in `1..=MAX_CELLS_PER_AXIS`, so the cast
    // is exact rather than saturating.
    Ok(whole as usize)
}

fn check_finite_deg(parameter: &'static str, value_deg: f64) -> Result<(), BasinBoundsError> {
    if value_deg.is_finite() {
        return Ok(());
    }
    Err(BasinBoundsError::NotFinite {
        parameter,
        value_deg,
    })
}

fn check_latitude(parameter: &'static str, value_deg: f64) -> Result<(), BasinBoundsError> {
    check_finite_deg(parameter, value_deg)?;
    if value_deg.abs() <= POLE_LATITUDE_DEG {
        return Ok(());
    }
    Err(BasinBoundsError::LatitudeOffThePlanet {
        parameter,
        value_deg,
    })
}
