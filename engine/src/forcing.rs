//! The wind forcing: what the alizés are, and how they land on the C-grid.
//!
//! Forcing is a *function of position and time*, not a field: `CONTEXT.md`
//! defines wind stress as `τx(x, y, t)`, `τy(x, y, t)`, and the scientific
//! model doc asks for it to be pluggable so a new scenario arrives as a new
//! implementation rather than an edit to the solver. [`WindStress`] is that
//! plug — one method, `stress(x, y, t) -> (τx, τy)` in pascals — and
//! [`SteadyTradeWinds`] is the control scenario: a steady easterly stress,
//! `τx < 0`, optionally decaying away from the equator.
//!
//! The solver cannot integrate a function, though. It needs the stress at the
//! points its momentum equations live on, which on the Arakawa C-grid of
//! [ADR-0003] are the east/west faces for `τx` and the north/south faces for
//! `τy`. [`WindStressField`] is that discretisation — the trait sampled onto
//! one basin at one instant — and [`Basin`] is what turns a `(staggering, i,
//! j)` into the `(x, y)` in metres the trait is asked about.
//!
//! # Why the basin's walls carry no stress
//!
//! [`WindStressField::sample`] leaves the four wall faces at exactly zero: the
//! `τx` columns `i = 0` and `i = nx`, and the `τy` rows `j = 0` and `j = ny`.
//! This is the same rule the C-grid derivative operators of T-01.1 already
//! apply to the pressure gradient, and for the same reason — a wall face has
//! water on one side only, so a stress there would accelerate a velocity that
//! is not a degree of freedom of a closed basin.
//!
//! It matters more than it looks. Nothing in the engine yet holds the wall
//! velocities to zero: no-normal-flow is T-04.2's, and until it lands the only
//! thing keeping `u` at the coast at rest is that no term forces it there. The
//! pressure gradient does not (the operators write zero), the Coriolis term
//! does not (it interpolates a velocity that is itself zero at the wall), and
//! damping cannot start a flow. A stress applied at the wall would, and a
//! basin whose coasts pass water does not tilt at all — the wind simply
//! accelerates the whole layer westward until damping balances it, which is
//! the open-channel solution and not the equatorial Pacific's. Zeroing the
//! wall stress is what makes the closed-basin steady state of
//! `tests/wind_forcing.rs` the one the physics predicts.
//!
//! It is a *sampling* rule, deliberately, and not a boundary condition: it
//! lives in this module because it says where a prescribed field is defined,
//! and it is one line for T-04.2 to subsume once the boundary owns the wall
//! velocities outright. [`WindStressField::uniform`] deliberately does *not*
//! apply it — the Epic 02 term tests use it to probe what the right-hand side
//! does with a stress at the wall, and that question stays askable.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use std::fmt;

use termocline_grid::{Field2D, Grid, Staggering, U_STAGGERING, V_STAGGERING};
use termocline_numerics::Spacing;

/// Stress of a calm ocean surface, in Pa.
const CALM: f64 = 0.0;

/// Half of something, in the "midpoint of an extent" sense.
const HALF: f64 = 0.5;

/// Why a basin could not be placed on the plane.
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

/// Where a basin sits in the plane, and therefore where each C-grid point of
/// it is.
///
/// The grid knows the *shape* of a basin in cells and the spacing knows how
/// wide a cell is; neither knows where the basin's southwest corner lies, and
/// a forcing stated as `τ(x, y, t)` cannot be sampled without it. `Basin` is
/// that missing origin, and the map from `(staggering, i, j)` to metres that
/// it implies.
///
/// `y` is measured north from the equator, so that
/// [`Basin::y_of_row_m`] and
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
            southern_edge_y_m: -HALF * grid.ny() as f64 * spacing.dy_m(),
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

    /// Width of the basin, in metres — the `L` of the analytic tilt.
    #[must_use]
    pub fn zonal_extent_m(self) -> f64 {
        self.grid.nx() as f64 * self.spacing.dx_m()
    }

    /// Height of the basin, in metres.
    #[must_use]
    pub fn meridional_extent_m(self) -> f64 {
        self.grid.ny() as f64 * self.spacing.dy_m()
    }

    /// Zonal position of the column `i` of a field at `staggering`, in metres.
    ///
    /// The half-cell offset that separates a cell-center column from an
    /// east/west-face column comes from [`Staggering::offset_in_cells`] rather
    /// than from a literal here, per CODING_STANDARDS.md § Scope guards: the
    /// grid knows about staggering, the physics does not.
    #[must_use]
    pub fn x_of_column_m(self, staggering: Staggering, i: usize) -> f64 {
        let (offset_in_cells, _) = staggering.offset_in_cells();
        self.western_edge_x_m + (i as f64 + offset_in_cells) * self.spacing.dx_m()
    }

    /// Meridional position of the row `j` of a field at `staggering`, in
    /// metres north of the equator.
    #[must_use]
    pub fn y_of_row_m(self, staggering: Staggering, j: usize) -> f64 {
        let (_, offset_in_cells) = staggering.offset_in_cells();
        self.southern_edge_y_m + (j as f64 + offset_in_cells) * self.spacing.dy_m()
    }
}

fn check_finite(parameter: &'static str, value_m: f64) -> Result<(), BasinError> {
    if value_m.is_finite() {
        return Ok(());
    }
    Err(BasinError::NotFinite { parameter, value_m })
}

/// A prescribed surface wind stress, as a function of position and time.
///
/// One method, because that is the whole contract: given a point of the basin
/// and an instant, say what the stress there is. Everything a scenario varies
/// — how strong the alizés are, how far off the equator they reach, whether
/// they breathe with the season or carry a burst — is a different
/// implementation of this trait, and none of it reaches the solver, which only
/// ever sees the [`WindStressField`] the trait was sampled into.
///
/// Implementations must be pure functions of `(x, y, t)`: the same arguments
/// give the same stress, every time and in any order. Runs are deterministic
/// (CODING_STANDARDS.md § Correctness and failure), and RK4 samples the
/// forcing four times per step at three distinct times, so a stateful
/// implementation would make the integration depend on stage order.
pub trait WindStress {
    /// The stress `(τx, τy)` in pascals at `x_m` metres east, `y_m` metres
    /// north of the equator, `t_s` seconds into the run.
    ///
    /// Easterly stress — the alizés — is `τx < 0` (`CONTEXT.md`, *Wind
    /// stress*).
    fn stress(&self, x_m: f64, y_m: f64, t_s: f64) -> (f64, f64);
}

/// Why a wind-stress scenario could not be built.
///
/// Both variants describe invalid *scenario input*, so they are returned
/// rather than panicked, and each names the value it rejected and the bound it
/// violated (CODING_STANDARDS.md § Correctness and failure).
#[derive(Debug, Clone, PartialEq)]
pub enum WindStressError {
    /// The zonal stress was not easterly. The alizés blow from the east, which
    /// is `τx < 0`; a positive or zero value describes some other wind, and
    /// naming it `SteadyTradeWinds` would be a lie.
    NotEasterly {
        /// The zonal stress supplied, in Pa.
        value_pa: f64,
    },
    /// A length scale was not a finite, strictly positive distance.
    ScaleNotPositive {
        /// Name of the parameter, matching its accessor.
        parameter: &'static str,
        /// The value supplied, in metres.
        value_m: f64,
    },
}

impl fmt::Display for WindStressError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotEasterly { value_pa } => write!(
                f,
                "the zonal trade-wind stress is {value_pa} Pa; the alizés blow from the east, \
                 so it must be strictly negative"
            ),
            Self::ScaleNotPositive { parameter, value_m } => write!(
                f,
                "{parameter} is {value_m} m; it must be a finite, strictly positive distance"
            ),
        }
    }
}

impl std::error::Error for WindStressError {}

/// The steady easterly trade winds — the control scenario of
/// `docs/planning/01-scientific-model.md`.
///
/// A zonal stress that does not vary with `x` or with `t`, and decays away
/// from the equator as a Gaussian in `y`:
///
/// ```text
/// τx(x, y, t) = τ₀ · exp(−(y / Ly)²)        τy(x, y, t) = 0
/// ```
///
/// with `τ₀ < 0` the stress on the equator and `Ly` the meridional decay
/// scale. [`SteadyTradeWinds::uniform`] is the `Ly → ∞` limit, the profile
/// that admits a closed-form steady state and therefore the one the analytic
/// tilt check runs in.
///
/// There is no meridional stress: the alizés are zonal to the accuracy this
/// model cares about, and a `τy` would drive an Ekman response the linear
/// core has nothing to say about yet.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteadyTradeWinds {
    /// Zonal stress `τ₀` on the equator, in Pa. Strictly negative.
    equatorial_zonal_stress_pa: f64,
    /// Meridional decay scale `Ly`, in metres — the `y` at which the stress
    /// has fallen to `1/e` of its equatorial value. Infinite for a uniform
    /// field.
    meridional_decay_scale_m: f64,
}

impl SteadyTradeWinds {
    /// Trade winds of `equatorial_zonal_stress_pa` everywhere, with no
    /// meridional structure at all.
    ///
    /// The `Ly → ∞` limit of [`SteadyTradeWinds::with_meridional_decay`], and
    /// the one case where the closed basin has a closed-form steady state: a
    /// stress independent of `y` is balanced by a thermocline tilt
    /// independent of `y`, so the whole problem collapses to the
    /// one-dimensional zonal balance the acceptance test checks against.
    ///
    /// # Errors
    /// [`WindStressError::NotEasterly`] unless the stress is strictly
    /// negative.
    pub fn uniform(equatorial_zonal_stress_pa: f64) -> Result<Self, WindStressError> {
        check_easterly(equatorial_zonal_stress_pa)?;
        Ok(Self {
            equatorial_zonal_stress_pa,
            meridional_decay_scale_m: f64::INFINITY,
        })
    }

    /// Trade winds of `equatorial_zonal_stress_pa` on the equator, falling to
    /// `1/e` of that at `meridional_decay_scale_m` either side of it.
    ///
    /// The realistic choice of scale is the equatorial deformation radius
    /// `Le = √(c/β)` (`CONTEXT.md`), the width of the waveguide the stress is
    /// meant to drive, but the scale is scenario input rather than a constant
    /// here: it is exactly the knob the forcing sensitivity of Epic 07 varies.
    ///
    /// # Errors
    /// [`WindStressError::NotEasterly`] unless the stress is strictly
    /// negative, or [`WindStressError::ScaleNotPositive`] unless the decay
    /// scale is a finite, strictly positive distance.
    pub fn with_meridional_decay(
        equatorial_zonal_stress_pa: f64,
        meridional_decay_scale_m: f64,
    ) -> Result<Self, WindStressError> {
        check_easterly(equatorial_zonal_stress_pa)?;
        if !meridional_decay_scale_m.is_finite() || meridional_decay_scale_m <= 0.0 {
            return Err(WindStressError::ScaleNotPositive {
                parameter: "meridional_decay_scale_m",
                value_m: meridional_decay_scale_m,
            });
        }
        Ok(Self {
            equatorial_zonal_stress_pa,
            meridional_decay_scale_m,
        })
    }

    /// Zonal stress `τ₀` on the equator, in Pa.
    #[must_use]
    pub const fn equatorial_zonal_stress_pa(self) -> f64 {
        self.equatorial_zonal_stress_pa
    }

    /// Meridional decay scale `Ly`, in metres; infinite for a uniform field.
    #[must_use]
    pub const fn meridional_decay_scale_m(self) -> f64 {
        self.meridional_decay_scale_m
    }
}

impl WindStress for SteadyTradeWinds {
    fn stress(&self, _x_m: f64, y_m: f64, _t_s: f64) -> (f64, f64) {
        // `y / ∞ = 0` and `exp(−0) = 1`, so the uniform case falls out of the
        // same expression rather than needing a branch of its own.
        let scaled = y_m / self.meridional_decay_scale_m;
        (
            self.equatorial_zonal_stress_pa * (-scaled * scaled).exp(),
            CALM,
        )
    }
}

fn check_easterly(value_pa: f64) -> Result<(), WindStressError> {
    if value_pa.is_finite() && value_pa < 0.0 {
        return Ok(());
    }
    Err(WindStressError::NotEasterly { value_pa })
}

/// A surface wind stress field over one basin, in pascals: a [`WindStress`]
/// evaluated at the C-grid points the momentum equations need it on.
///
/// `τx` sits on the east/west faces with the zonal current anomaly `u`, and
/// `τy` on the north/south faces with `v`, so the momentum equations pick up
/// `τ/(ρ₀·H)` without an interpolation.
#[derive(Debug, Clone, PartialEq)]
pub struct WindStressField {
    /// Shape of the basin the two fields cover.
    grid: Grid,
    /// Zonal wind stress `τx`, in Pa, at east/west faces. Negative is
    /// easterly — the direction the trade winds blow.
    tau_x_pa: Field2D<f64>,
    /// Meridional wind stress `τy`, in Pa, at north/south faces.
    tau_y_pa: Field2D<f64>,
}

impl WindStressField {
    /// No wind at all over `grid`: both components exactly zero.
    ///
    /// The unforced limit the wave tests of Epic 07 run in.
    #[must_use]
    pub fn calm(grid: Grid) -> Self {
        Self::uniform(grid, CALM, CALM)
    }

    /// A stress of `tau_x_pa` by `tau_y_pa` pascals at *every* face of `grid`,
    /// the basin's walls included.
    ///
    /// The raw constructor, and the one the Epic 02 right-hand-side tests use
    /// to ask what the momentum equations do with a stress at a wall. A field
    /// sampled from a [`WindStress`] leaves the walls at zero instead — see
    /// [`WindStressField::sample`] and this module's header for why the
    /// difference matters.
    #[must_use]
    pub fn uniform(grid: Grid, tau_x_pa: f64, tau_y_pa: f64) -> Self {
        Self {
            grid,
            tau_x_pa: grid.allocate(U_STAGGERING, tau_x_pa),
            tau_y_pa: grid.allocate(V_STAGGERING, tau_y_pa),
        }
    }

    /// `wind` sampled over `basin` at `t_s` seconds, as a freshly allocated
    /// field.
    ///
    /// The convenient form for a test or a steady scenario, which samples once
    /// and reuses the result for the whole run. A time-varying scenario steps
    /// a field it already owns through [`WindStressField::sample`] instead, so
    /// that a run allocates its forcing exactly once
    /// (CODING_STANDARDS.md § Performance).
    #[must_use]
    pub fn sampled<W: WindStress + ?Sized>(basin: Basin, wind: &W, t_s: f64) -> Self {
        let mut field = Self::calm(basin.grid());
        field.sample(basin, wind, t_s);
        field
    }

    /// Overwrite this field with `wind` sampled over `basin` at `t_s` seconds.
    ///
    /// Every interior face is written, so the same buffer can be re-sampled at
    /// each RK4 stage without carrying a stage's values into the next. The
    /// basin's wall faces — the `τx` columns `i = 0` and `i = nx`, and the
    /// `τy` rows `j = 0` and `j = ny` — are set to exactly zero rather than to
    /// `wind`'s value there: a wall face has water on one side only, and this
    /// module's header explains at length why forcing it would open the
    /// closed basin.
    ///
    /// # Panics
    /// If `basin` covers a different grid from the one this field was built
    /// for. A shape mismatch means the calling code is wrong, which is what
    /// panics are for (CODING_STANDARDS.md § Correctness and failure).
    pub fn sample<W: WindStress + ?Sized>(&mut self, basin: Basin, wind: &W, t_s: f64) {
        assert!(
            basin.grid() == self.grid,
            "basin covers {:?}, but this wind stress field was built for {:?}",
            basin.grid(),
            self.grid
        );
        let nx = self.grid.nx();
        let ny = self.grid.ny();
        write_component(
            &mut self.tau_x_pa,
            basin,
            U_STAGGERING,
            |i, _j| i == 0 || i == nx,
            |stress| stress.0,
            wind,
            t_s,
        );
        write_component(
            &mut self.tau_y_pa,
            basin,
            V_STAGGERING,
            |_i, j| j == 0 || j == ny,
            |stress| stress.1,
            wind,
            t_s,
        );
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

/// Write one component of `wind` into `component`, zeroing the faces
/// `is_wall` names.
///
/// Shared by the two halves of [`WindStressField::sample`], which differ only
/// in where their points sit, which wall faces they have, and which half of
/// the returned pair they keep.
fn write_component<W, Wall, Pick>(
    component: &mut Field2D<f64>,
    basin: Basin,
    staggering: Staggering,
    is_wall: Wall,
    pick: Pick,
    wind: &W,
    t_s: f64,
) where
    W: WindStress + ?Sized,
    Wall: Fn(usize, usize) -> bool,
    Pick: Fn((f64, f64)) -> f64,
{
    for j in 0..component.ny() {
        let y_m = basin.y_of_row_m(staggering, j);
        for i in 0..component.nx() {
            let value = if is_wall(i, j) {
                CALM
            } else {
                pick(wind.stress(basin.x_of_column_m(staggering, i), y_m, t_s))
            };
            *component
                .get_mut(i, j)
                .expect("the loop bounds are the field's own extents") = value;
        }
    }
}
