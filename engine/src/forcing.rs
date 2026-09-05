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
//! [`SeasonalTradeWinds`] is the second scenario of the scientific model doc:
//! the same field breathing with the year, `1 + a·cos(2π(t − t_peak)/T_year)`.
//! It is a wrapper rather than a variant of [`SteadyTradeWinds`] because the
//! two answer different questions — what the alizés look like in space, and
//! how strong they are this month — and keeping them apart is what lets the
//! composable burst of T-03.3 stack on either.
//!
//! The solver cannot integrate a function, though. It needs the stress at the
//! points its momentum equations live on, which on the Arakawa C-grid of
//! [ADR-0003] are the east/west faces for `τx` and the north/south faces for
//! `τy`. [`WindStressField`] is that discretisation — the trait sampled onto
//! one basin at one instant — and [`Basin`](crate::Basin), from the `basin`
//! module, is what turns a `(staggering, i, j)` into the `(x, y)` in metres
//! the trait is asked about.
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

use std::f64::consts::TAU;
use std::fmt;

use termocline_grid::{Field2D, Grid, Staggering, U_STAGGERING, V_STAGGERING};

use crate::basin::Basin;

/// Stress of a calm ocean surface, in Pa.
const CALM: f64 = 0.0;

/// One mean solar day, in seconds.
const SOLAR_DAY_S: f64 = 86_400.0;

/// The tropical year `T_year`, in seconds — the period of the seasonal cycle.
///
/// 365.2422 mean solar days, the mean tropical year of the *Astronomical
/// Almanac*: the equinox-to-equinox year the seasons follow, rather than the
/// sidereal year or the calendar's 365. It is a `const` and not scenario input
/// because the year is not a knob — T-03.2 asks for the *amplitude* and
/// *phase* of the annual harmonic to be configurable, and a scenario wanting
/// some other period is asking for a different forcing, not for a differently
/// tuned season.
pub const TROPICAL_YEAR_S: f64 = 365.2422 * SOLAR_DAY_S;

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
    /// The seasonal modulation amplitude was not a fraction of the steady
    /// field. Outside `[0, 1]` the annual harmonic `1 + a·cos(…)` turns
    /// negative somewhere in the year, which flips the stress westerly; a
    /// westerly stress is the wind burst of T-03.3, not a season, and a
    /// scenario named for the alizés must not quietly become one.
    ModulationNotAFraction {
        /// The value supplied, dimensionless.
        relative_amplitude: f64,
    },
    /// The seasonal phase was not a finite instant.
    PhaseNotFinite {
        /// The value supplied, in seconds.
        peak_time_s: f64,
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
            Self::ModulationNotAFraction { relative_amplitude } => write!(
                f,
                "relative_amplitude is {relative_amplitude}; it must be a fraction between 0 \
                 and 1, or the modulated alizés would turn westerly within the year"
            ),
            Self::PhaseNotFinite { peak_time_s } => write!(
                f,
                "peak_time_s is {peak_time_s} s; it must be a finite instant"
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
    /// has fallen to `1/e` of its equatorial value. `None` for a field with no
    /// meridional structure at all, which is not the same thing as a scale of
    /// any particular size.
    meridional_decay_scale_m: Option<f64>,
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
            meridional_decay_scale_m: None,
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
            meridional_decay_scale_m: Some(meridional_decay_scale_m),
        })
    }

    /// Zonal stress `τ₀` on the equator, in Pa.
    #[must_use]
    pub const fn equatorial_zonal_stress_pa(self) -> f64 {
        self.equatorial_zonal_stress_pa
    }

    /// Meridional decay scale `Ly`, in metres, or `None` for a field with no
    /// meridional structure.
    #[must_use]
    pub const fn meridional_decay_scale_m(self) -> Option<f64> {
        self.meridional_decay_scale_m
    }
}

impl WindStress for SteadyTradeWinds {
    fn stress(&self, _x_m: f64, y_m: f64, _t_s: f64) -> (f64, f64) {
        let decay = match self.meridional_decay_scale_m {
            None => 1.0,
            Some(scale_m) => {
                let scaled = y_m / scale_m;
                (-scaled * scaled).exp()
            }
        };
        (self.equatorial_zonal_stress_pa * decay, CALM)
    }
}

/// The trade winds breathing with the year — the seasonal-cycle scenario of
/// `docs/planning/01-scientific-model.md`.
///
/// A [`SteadyTradeWinds`] field scaled by an annual harmonic, the same factor
/// everywhere in the basin at a given instant:
///
/// ```text
/// τ(x, y, t) = τ_steady(x, y) · (1 + a·cos(2π·(t − t_peak)/T_year))
/// ```
///
/// with `a` the relative amplitude, `t_peak` the phase — written as the
/// instant the alizés are *strongest* rather than as an angle, so that it
/// carries a unit like every other quantity here — and `T_year` the
/// [`TROPICAL_YEAR_S`].
///
/// The modulation is a pure scaling, so it does not move the wind's structure
/// in `y`: the alizés of March and of September have the same shape and
/// different strength. Meridional migration of the wind belt is a real feature
/// of the seasonal cycle and deliberately not modelled here — the ticket asks
/// for an annual harmonic on the steady field, and a migrating belt is a
/// different scenario.
///
/// `a` is required to be a fraction, so the harmonic never turns negative and
/// the stress never reverses. At `a = 1` the basin goes momentarily calm once
/// a year, which is the strongest season this scenario can describe; anything
/// beyond it is a westerly wind burst wearing a season's name, and bursts are
/// T-03.3's, superimposed rather than substituted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeasonalTradeWinds {
    /// The field the harmonic modulates.
    steady: SteadyTradeWinds,
    /// Relative amplitude `a` of the annual harmonic, dimensionless, in
    /// `[0, 1]`.
    relative_amplitude: f64,
    /// Phase, in seconds: the instant at which the harmonic peaks and the
    /// alizés are therefore strongest.
    peak_time_s: f64,
}

impl SeasonalTradeWinds {
    /// `steady` modulated by an annual harmonic of relative amplitude
    /// `relative_amplitude`, peaking `peak_time_s` seconds into the run.
    ///
    /// # Errors
    /// [`WindStressError::ModulationNotAFraction`] unless the amplitude is a
    /// number in `[0, 1]`, and [`WindStressError::PhaseNotFinite`] unless the
    /// phase is a finite instant.
    pub fn new(
        steady: SteadyTradeWinds,
        relative_amplitude: f64,
        peak_time_s: f64,
    ) -> Result<Self, WindStressError> {
        if !(0.0..=1.0).contains(&relative_amplitude) {
            return Err(WindStressError::ModulationNotAFraction { relative_amplitude });
        }
        if !peak_time_s.is_finite() {
            return Err(WindStressError::PhaseNotFinite { peak_time_s });
        }
        Ok(Self {
            steady,
            relative_amplitude,
            peak_time_s,
        })
    }

    /// The steady field this season modulates.
    #[must_use]
    pub const fn steady(self) -> SteadyTradeWinds {
        self.steady
    }

    /// Relative amplitude `a` of the annual harmonic, dimensionless.
    #[must_use]
    pub const fn relative_amplitude(self) -> f64 {
        self.relative_amplitude
    }

    /// Phase, in seconds: the instant the alizés are strongest.
    #[must_use]
    pub const fn peak_time_s(self) -> f64 {
        self.peak_time_s
    }

    /// The harmonic `1 + a·cos(2π(t − t_peak)/T_year)` at `t_s`, dimensionless.
    ///
    /// Never negative, because `a ∈ [0, 1]` and `cos ≥ −1`: that is what keeps
    /// the modulated alizés easterly.
    fn modulation(self, t_s: f64) -> f64 {
        let phase_rad = TAU * (t_s - self.peak_time_s) / TROPICAL_YEAR_S;
        self.relative_amplitude.mul_add(phase_rad.cos(), 1.0)
    }
}

impl WindStress for SeasonalTradeWinds {
    fn stress(&self, x_m: f64, y_m: f64, t_s: f64) -> (f64, f64) {
        let modulation = self.modulation(t_s);
        let (tau_x_pa, tau_y_pa) = self.steady.stress(x_m, y_m, t_s);
        (tau_x_pa * modulation, tau_y_pa * modulation)
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
        Self::uniform_including_walls(grid, CALM, CALM)
    }

    /// A stress of `tau_x_pa` by `tau_y_pa` pascals at *every* face of `grid`,
    /// the basin's walls included.
    ///
    ///
    /// The raw constructor, and the one the Epic 02 right-hand-side tests use
    /// to ask what the momentum equations do with a stress at a wall. A field
    /// sampled from a [`WindStress`] leaves the walls at zero instead — see
    /// [`WindStressField::sample`] and this module's header for why the
    /// difference matters.
    #[must_use]
    pub fn uniform_including_walls(grid: Grid, tau_x_pa: f64, tau_y_pa: f64) -> Self {
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
