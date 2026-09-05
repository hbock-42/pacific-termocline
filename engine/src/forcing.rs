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
//! That burst is [`WindBurstAnomaly`], and stacking is [`CompositeWind`].
//! Scenarios add rather than exclude one another: a westerly wind burst is an
//! anomaly *superimposed on* the trades (`CONTEXT.md`), and the equations are
//! linear in the stress, so "superimposed on" is a sum of two scenarios that
//! each remain a `WindStress` in their own right.
//!
//! The solver cannot integrate a function, though. It needs the stress at the
//! points its momentum equations live on, which on the Arakawa C-grid of
//! [ADR-0003] are the east/west faces for `τx` and the north/south faces for
//! `τy`. [`WindStressField`] is that discretisation — the trait sampled onto
//! one basin at one instant — and [`Basin`](crate::Basin), from the `basin`
//! module, is what turns a `(staggering, i, j)` into the `(x, y)` in metres
//! the trait is asked about.
//!
//! # The wind blows over the coast too
//!
//! A sampled field carries the wind's stress at every face it has, the
//! basin's four wall lines included. That is what the trait says the stress
//! there is, and this module's job is to report it, not to edit it.
//!
//! It was not always so. T-03.1 found that a stress at the wall opens the
//! closed basin — the coasts pass water, and instead of tilting, the wind
//! accelerates the whole layer westward until damping balances it, which is
//! the open-channel solution and not the equatorial Pacific's — and, with no
//! boundary condition in the engine yet, worked around it by zeroing the wall
//! faces here as a deliberately interim *sampling* rule. T-04.2 replaced that
//! with the thing it stood in for: [`NoNormalFlow`](crate::NoNormalFlow) holds
//! the wall velocities at rest at every RK4 stage, so the `τ/(ρ₀·H)` a wall
//! face receives is discarded where the integration happens. One invariant,
//! one owner; a forcing field is no longer part of what keeps the basin
//! closed.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use std::f64::consts::TAU;
use std::fmt;

use termocline_grid::{Field2D, Grid, Staggering, U_STAGGERING, V_STAGGERING};

use crate::basin::Basin;
use crate::state::OceanState;

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

    /// Whether this wind depends on `t` at all — what lets
    /// [`WindForcing`] reuse a field it has already sampled.
    ///
    /// The default is [`TimeDependence::Varying`], which is the answer that is
    /// always safe: an implementation that says nothing is re-sampled at every
    /// instant, exactly as before this method existed. Overriding it with
    /// [`TimeDependence::Steady`] is a *promise* that `stress(x, y, t)` does
    /// not depend on `t`, and breaking that promise silently freezes the
    /// forcing of a whole run.
    fn time_dependence(&self) -> TimeDependence {
        TimeDependence::Varying
    }
}

/// Whether a [`WindStress`] varies with time.
///
/// The declaration T-10.5 turns on. `docs/performance-notes.md` measured
/// re-sampling the wind at 71% of a timestep, all of it recomputing what the
/// previous RK4 stage had already computed — but only for a forcing that does
/// not actually move, and [`SeasonalTradeWinds`] and [`WindBurstAnomaly`]
/// genuinely do. So the property is asked of the wind rather than assumed of
/// the scenario, and the variant that costs nothing to be wrong about is the
/// default. [ADR-0009] is why it is a declaration rather than something the
/// engine works out, and what a wrong one costs.
///
/// [ADR-0009]: ../../docs/planning/adr/0009-wind-declares-its-time-dependence.md
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeDependence {
    /// `stress(x, y, t)` is the same for every `t`: one sampling is the field
    /// for the whole run.
    Steady,
    /// The stress moves with time, so a field is only the field of the instant
    /// it was sampled at.
    Varying,
}

impl TimeDependence {
    /// The dependence of a sum of two winds: steady only if both are, because
    /// a sum is constant in time exactly when every term is.
    #[must_use]
    const fn and(self, other: Self) -> Self {
        match (self, other) {
            (Self::Steady, Self::Steady) => Self::Steady,
            _ => Self::Varying,
        }
    }
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
    /// The zonal stress of a burst was not westerly. A westerly wind burst is
    /// an anomaly *against* the alizés (`CONTEXT.md`), which is `τx > 0`; a
    /// negative or zero value describes a strengthening of the trades or no
    /// burst at all, and either is a different scenario.
    NotWesterly {
        /// The zonal stress supplied, in Pa.
        value_pa: f64,
    },
    /// A burst's duration was not a finite, strictly positive time.
    DurationNotPositive {
        /// The value supplied, in seconds.
        duration_s: f64,
    },
    /// A burst's zonal centre was not a finite position.
    CenterNotAPosition {
        /// The value supplied, in metres.
        center_x_m: f64,
    },
    /// A burst's peak time was not a finite instant.
    ///
    /// The burst's own twin of [`WindStressError::PhaseNotFinite`], which
    /// rejects the same thing for the seasonal cycle: the two carry different
    /// scenario parameters, and an error that named neither would leave a
    /// scenario author guessing which one they got wrong.
    PeakTimeNotFinite {
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
            Self::NotWesterly { value_pa } => write!(
                f,
                "the peak zonal stress of the burst is {value_pa} Pa; a westerly wind burst \
                 blows against the alizés, so it must be strictly positive"
            ),
            Self::DurationNotPositive { duration_s } => write!(
                f,
                "duration_s is {duration_s} s; it must be a finite, strictly positive time"
            ),
            Self::CenterNotAPosition { center_x_m } => write!(
                f,
                "center_x_m is {center_x_m} m; it must be a finite position"
            ),
            Self::PeakTimeNotFinite { peak_time_s } => write!(
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
/// scale. [`SteadyTradeWinds::uniform`] is the `Ly → ∞` limit, the profile the
/// analytic tilt check of T-07.4 runs in.
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
    /// the forcing whose steady state has a closed form: a stress independent
    /// of `y` drives a thermocline tilt independent of `y`, so the problem
    /// collapses to the one-dimensional zonal balance
    /// `tests/steady_wind_tilt.rs` checks against.
    ///
    /// That collapse needs the meridional momentum equation to be satisfied by
    /// the collapsed solution, and `f·u + g'·∂h/∂y = 0` with a `y`-independent
    /// `h` says `f·u = 0`. So it is exact in a basin one cell tall — where
    /// [`NoNormalFlow`](crate::NoNormalFlow) pins `v` at zero and the rotation
    /// terms drop out — and only approximate in a rotating one, where the
    /// steady current the damping requires is turned by the beta-plane. T-07.4
    /// runs the closed-form comparison in the first and the exact modal
    /// balances of the equatorial waveguide in the second.
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
        check_scale("meridional_decay_scale_m", meridional_decay_scale_m)?;
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
            Some(scale_m) => gaussian(y_m, scale_m),
        };
        (self.equatorial_zonal_stress_pa * decay, CALM)
    }

    /// Steady is the scenario's name and its whole content: `stress` above
    /// does not read `t_s`.
    fn time_dependence(&self) -> TimeDependence {
        TimeDependence::Steady
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

    /// A season varies with time: that is what it is for.
    ///
    /// Stated of the type rather than of the instance, even though a season of
    /// zero amplitude is a steady field in exact arithmetic and in `f64`. The
    /// ticket asks for the winds that genuinely vary to stay correct, not for
    /// the ones that could be talked out of varying to go faster, and a second
    /// claimant to `Steady` is a second thing a reader has to check.
    fn time_dependence(&self) -> TimeDependence {
        TimeDependence::Varying
    }
}

fn check_easterly(value_pa: f64) -> Result<(), WindStressError> {
    if value_pa.is_finite() && value_pa < 0.0 {
        return Ok(());
    }
    Err(WindStressError::NotEasterly { value_pa })
}

fn check_westerly(value_pa: f64) -> Result<(), WindStressError> {
    if value_pa.is_finite() && value_pa > 0.0 {
        return Ok(());
    }
    Err(WindStressError::NotWesterly { value_pa })
}

/// A length scale is a finite, strictly positive distance, whatever it scales.
fn check_scale(parameter: &'static str, value_m: f64) -> Result<(), WindStressError> {
    if value_m.is_finite() && value_m > 0.0 {
        return Ok(());
    }
    Err(WindStressError::ScaleNotPositive { parameter, value_m })
}

/// `exp(−(offset / scale)²)`, the Gaussian factor every profile in this module
/// is built from. `scale` is checked strictly positive at construction, so the
/// division is safe.
fn gaussian(offset: f64, scale: f64) -> f64 {
    let scaled = offset / scale;
    (-scaled * scaled).exp()
}

/// An idealized westerly wind burst: a positive-`τx` anomaly, Gaussian in `x`,
/// in `y` about the equator, and in `t` (`CONTEXT.md`, *Westerly wind burst*).
///
/// ```text
///                       ⎛   ⎛x − x₀⎞²⎞     ⎛   ⎛ y⎞²⎞     ⎛   ⎛t − t₀⎞²⎞
/// τx(x, y, t) = τ_burst·exp⎜− ⎜──────⎟ ⎟·exp⎜− ⎜──⎟ ⎟·exp⎜− ⎜──────⎟ ⎟
///                       ⎝   ⎝  Lx  ⎠ ⎠     ⎝   ⎝Ly⎠ ⎠     ⎝   ⎝  Lt  ⎠ ⎠
/// ```
///
/// with `τ_burst > 0`, the opposite sign to the alizés: the burst blows
/// *against* the trades, which is the perturbation known to trigger El Niño
/// onset. It is meant to be added to a base scenario rather than to replace
/// one — a burst on its own is not a state of the equatorial Pacific — and
/// [`CompositeWind`] is what performs that addition.
///
/// The meridional Gaussian is centred on the equator, like
/// [`SteadyTradeWinds`]'s: the bursts this models are equatorial, and the
/// waveguide they are meant to excite is centred there too. The three factors
/// multiply rather than combining into one distance, so each scale is
/// independent — a burst can be broad and brief, or narrow and long-lived.
///
/// There is no meridional stress, for [`SteadyTradeWinds`]'s reason: a `τy`
/// would drive an Ekman response the linear core has nothing to say about yet.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindBurstAnomaly {
    /// Peak zonal stress `τ_burst`, in Pa, at the centre of all three
    /// Gaussians. Strictly positive — westerly.
    peak_zonal_stress_pa: f64,
    /// Zonal centre `x₀` of the burst, in metres.
    center_x_m: f64,
    /// Zonal `e`-folding scale `Lx`, in metres.
    zonal_scale_m: f64,
    /// Meridional `e`-folding scale `Ly`, in metres, about the equator.
    meridional_scale_m: f64,
    /// Instant `t₀` of the burst's peak, in seconds since the start of the run.
    peak_time_s: f64,
    /// Temporal `e`-folding scale `Lt`, in seconds — how long the burst lasts.
    duration_s: f64,
}

impl WindBurstAnomaly {
    /// A burst peaking at `peak_zonal_stress_pa` at `center_x_m` on the
    /// equator, `peak_time_s` seconds into the run, and falling to `1/e` of
    /// that at `zonal_scale_m` east or west of its centre, at
    /// `meridional_scale_m` either side of the equator, and `duration_s`
    /// before or after its peak.
    ///
    /// The realistic choice of `meridional_scale_m` is the equatorial
    /// deformation radius `Le = √(c/β)` (`CONTEXT.md`) — the width of the
    /// waveguide the burst is meant to excite — but every scale here is
    /// scenario input, because they are exactly the knobs the forcing
    /// sensitivity of Epic 07 varies.
    ///
    /// # Errors
    /// [`WindStressError::NotWesterly`] unless the peak stress is strictly
    /// positive, [`WindStressError::ScaleNotPositive`] unless both length
    /// scales are finite, strictly positive distances,
    /// [`WindStressError::DurationNotPositive`] unless the duration is a
    /// finite, strictly positive time,
    /// [`WindStressError::CenterNotAPosition`] unless the zonal centre is a
    /// finite position, and [`WindStressError::PeakTimeNotFinite`] unless the
    /// peak time is a finite instant.
    pub fn new(
        peak_zonal_stress_pa: f64,
        center_x_m: f64,
        zonal_scale_m: f64,
        meridional_scale_m: f64,
        peak_time_s: f64,
        duration_s: f64,
    ) -> Result<Self, WindStressError> {
        check_westerly(peak_zonal_stress_pa)?;
        check_scale("zonal_scale_m", zonal_scale_m)?;
        check_scale("meridional_scale_m", meridional_scale_m)?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(WindStressError::DurationNotPositive { duration_s });
        }
        if !center_x_m.is_finite() {
            return Err(WindStressError::CenterNotAPosition { center_x_m });
        }
        if !peak_time_s.is_finite() {
            return Err(WindStressError::PeakTimeNotFinite { peak_time_s });
        }
        Ok(Self {
            peak_zonal_stress_pa,
            center_x_m,
            zonal_scale_m,
            meridional_scale_m,
            peak_time_s,
            duration_s,
        })
    }

    /// Peak zonal stress `τ_burst`, in Pa.
    #[must_use]
    pub const fn peak_zonal_stress_pa(self) -> f64 {
        self.peak_zonal_stress_pa
    }

    /// Zonal centre `x₀` of the burst, in metres.
    #[must_use]
    pub const fn center_x_m(self) -> f64 {
        self.center_x_m
    }

    /// Zonal `e`-folding scale `Lx`, in metres.
    #[must_use]
    pub const fn zonal_scale_m(self) -> f64 {
        self.zonal_scale_m
    }

    /// Meridional `e`-folding scale `Ly`, in metres.
    #[must_use]
    pub const fn meridional_scale_m(self) -> f64 {
        self.meridional_scale_m
    }

    /// Instant `t₀` of the burst's peak, in seconds since the start of the run.
    #[must_use]
    pub const fn peak_time_s(self) -> f64 {
        self.peak_time_s
    }

    /// Temporal `e`-folding scale `Lt` of the burst, in seconds.
    #[must_use]
    pub const fn duration_s(self) -> f64 {
        self.duration_s
    }
}

impl WindStress for WindBurstAnomaly {
    fn stress(&self, x_m: f64, y_m: f64, t_s: f64) -> (f64, f64) {
        let envelope = gaussian(x_m - self.center_x_m, self.zonal_scale_m)
            * gaussian(y_m, self.meridional_scale_m)
            * gaussian(t_s - self.peak_time_s, self.duration_s);
        (self.peak_zonal_stress_pa * envelope, CALM)
    }

    /// A burst is a Gaussian in `t` with a strictly positive duration, so it
    /// varies at every instant: there is no configuration of one that is
    /// steady.
    fn time_dependence(&self) -> TimeDependence {
        TimeDependence::Varying
    }
}

/// Several wind scenarios blowing at once: the pointwise sum of its
/// components.
///
/// The combinator the burst needs. `CONTEXT.md` calls a westerly wind burst an
/// anomaly *superimposed on* the trades, and the equations are linear in the
/// stress, so "superimposed on" is addition — a burst is stacked on a base
/// scenario rather than replacing it, and the same holds for the seasonal
/// modulation of T-03.2 or any pair of them together.
///
/// An empty composite is calm, so a scenario with no forcing at all is the
/// zero of this combinator rather than a special case. Components are summed
/// in the order they were added, which fixes the floating-point result and so
/// keeps runs deterministic (CODING_STANDARDS.md § Correctness and failure).
#[derive(Default)]
pub struct CompositeWind {
    /// The components, in the order they will be summed.
    components: Vec<Box<dyn WindStress>>,
}

impl CompositeWind {
    /// A composite of nothing at all: a calm ocean surface.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// This composite with `wind` added to it, by value — the builder form,
    /// for assembling a scenario in one expression.
    #[must_use]
    pub fn with<W: WindStress + 'static>(mut self, wind: W) -> Self {
        self.push(wind);
        self
    }

    /// Add `wind` to this composite, in place.
    pub fn push<W: WindStress + 'static>(&mut self, wind: W) {
        self.components.push(Box::new(wind));
    }

    /// How many components this composite sums.
    #[must_use]
    pub fn len(&self) -> usize {
        self.components.len()
    }

    /// Whether this composite has no components at all, and is therefore calm.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.components.is_empty()
    }
}

impl WindStress for CompositeWind {
    fn stress(&self, x_m: f64, y_m: f64, t_s: f64) -> (f64, f64) {
        self.components
            .iter()
            .fold((CALM, CALM), |(tau_x_pa, tau_y_pa), component| {
                let (component_x_pa, component_y_pa) = component.stress(x_m, y_m, t_s);
                (tau_x_pa + component_x_pa, tau_y_pa + component_y_pa)
            })
    }

    /// Steady only if every component is, which is superposition read
    /// backwards: the sum stops being constant in time as soon as one term
    /// does. An empty composite is calm, and calm does not change — the same
    /// reason it is the zero of this combinator.
    fn time_dependence(&self) -> TimeDependence {
        self.components
            .iter()
            .fold(TimeDependence::Steady, |dependence, component| {
                dependence.and(component.time_dependence())
            })
    }
}

impl fmt::Debug for CompositeWind {
    /// A `WindStress` is a function rather than data, so a composite can only
    /// report how many of them it holds.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CompositeWind")
            .field("components", &self.components.len())
            .finish()
    }
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

    /// A stress of `tau_x_pa` by `tau_y_pa` pascals at every face of `grid`,
    /// the basin's walls included.
    ///
    /// The raw constructor, and the one the Epic 02 right-hand-side tests use
    /// to ask what the momentum equations do with a stress at a wall — a
    /// question the boundary condition of T-04.2 answers at the solver rather
    /// than by leaving the stress unstated.
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
    /// Every face is written, the basin's wall lines included, so the same
    /// buffer can be re-sampled at each RK4 stage without carrying a stage's
    /// values into the next. What the solver does with a stress at a wall is
    /// the boundary condition's business, not this module's — see the header.
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
        write_component(
            &mut self.tau_x_pa,
            basin,
            U_STAGGERING,
            |stress| stress.0,
            wind,
            t_s,
        );
        write_component(
            &mut self.tau_y_pa,
            basin,
            V_STAGGERING,
            |stress| stress.1,
            wind,
            t_s,
        );
    }

    /// Overwrite this field with `other`, component for component.
    ///
    /// What a forcing that *adds* two winds needs: the sum has to start from
    /// one of them, and the one it starts from is a buffer somebody else owns
    /// and re-samples on its own schedule
    /// ([`CoupledWind`](crate::CoupledWind)).
    ///
    /// # Panics
    /// If `other` covers a different basin. A shape mismatch means the calling
    /// code is wrong, which is what panics are for (CODING_STANDARDS.md
    /// § Correctness and failure).
    pub fn assign(&mut self, other: &Self) {
        assert!(
            other.grid == self.grid,
            "that stress covers {:?}, but this one was built for {:?}",
            other.grid,
            self.grid
        );
        self.tau_x_pa
            .as_mut_slice()
            .copy_from_slice(other.tau_x_pa.as_slice());
        self.tau_y_pa
            .as_mut_slice()
            .copy_from_slice(other.tau_y_pa.as_slice());
    }

    /// Add `wind` sampled over `basin` at `t_s` seconds to what this field
    /// already holds.
    ///
    /// [`WindStressField::sample`] at the field level: the superposition of
    /// T-03.3, performed on a field rather than inside a
    /// [`CompositeWind`]. The two are the same sum — the equations are linear
    /// in the stress — and which one a forcing uses is decided by whether its
    /// components can share a cache, not by the physics.
    ///
    /// # Panics
    /// If `basin` covers a different grid from the one this field was built
    /// for.
    pub fn add_sampled<W: WindStress + ?Sized>(&mut self, basin: Basin, wind: &W, t_s: f64) {
        assert!(
            basin.grid() == self.grid,
            "basin covers {:?}, but this wind stress field was built for {:?}",
            basin.grid(),
            self.grid
        );
        add_component(
            &mut self.tau_x_pa,
            basin,
            U_STAGGERING,
            |stress| stress.0,
            wind,
            t_s,
        );
        add_component(
            &mut self.tau_y_pa,
            basin,
            V_STAGGERING,
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

/// Which instant a sampled field is the field of, and whether that can go
/// stale.
///
/// The whole invalidation rule of T-10.5, in one place because two forcings
/// need it: [`WindForcing`], which keeps a field across the steps of a run,
/// and [`BorrowedForcing`], which borrows the solver's buffer for the stages
/// of one step.
///
/// A held field may be reused on exactly two grounds, both of them properties
/// of the [`WindStress`] contract rather than of any scenario:
///
/// - the instant asked for is the instant held, bit for bit — a `WindStress`
///   is a pure function of `(x, y, t)`, so the same `t` gives the same field.
///   RK4's four stages ask about three instants (`t`, `t + dt/2` twice,
///   `t + dt`), and a step's last stage asks about what the next step's first
///   stage asks about, so this alone halves the sampling of a wind that varies
///   in time as fast as it likes;
/// - the wind declared itself [`TimeDependence::Steady`], in which case every
///   instant is the instant held.
///
/// Anything else re-samples. There is deliberately no third ground.
#[derive(Debug, Clone, Copy)]
struct HeldInstant {
    /// What the wind said about its own dependence on time. Captured once: a
    /// `WindStress` is a pure function, so its answer cannot change under an
    /// immutable wind.
    dependence: TimeDependence,
    /// The instant the field holds, in seconds, or `None` before anything has
    /// been sampled into it.
    sampled_t_s: Option<f64>,
}

impl HeldInstant {
    /// A cache holding nothing yet, for a wind of `dependence`.
    const fn nothing(dependence: TimeDependence) -> Self {
        Self {
            dependence,
            sampled_t_s: None,
        }
    }

    /// Whether the field this describes is already the field at `t_s`.
    fn holds(self, t_s: f64) -> bool {
        match self.sampled_t_s {
            None => false,
            Some(held_t_s) => matches!(self.dependence, TimeDependence::Steady) || held_t_s == t_s,
        }
    }

    /// `field` as the stress of `wind` over `basin` at `t_s`, sampling it only
    /// if what it holds is not already that.
    ///
    /// The one place the rule above is acted on, so that the two forcings
    /// below differ in what they own rather than in when they invalidate.
    fn ensure<'f, W: WindStress + ?Sized>(
        &mut self,
        field: &'f mut WindStressField,
        basin: Basin,
        wind: &W,
        t_s: f64,
    ) -> &'f WindStressField {
        if !self.holds(t_s) {
            field.sample(basin, wind, t_s);
            self.sampled_t_s = Some(t_s);
        }
        field
    }
}

/// Where one RK4 stage gets the surface stress at its own instant.
///
/// Implemented by the two forcings of this module, which differ only in what
/// they own: [`WindForcing`] holds a wind and a field for a whole run,
/// [`BorrowedForcing`] borrows a wind and the solver's buffer for one step.
/// The solver integrates against this, so it has one stage body rather than
/// one per way of holding a forcing.
/// [`CoupledWind`](crate::CoupledWind) is the third, and it is why `at` is
/// handed the state as well as the instant: the atmospheric response of Epic 12
/// answers the SST anomaly of the stage being evaluated, not the clock. The
/// two prescribed forcings ignore it, which is exactly what makes them
/// prescribed.
pub trait StageForcing {
    /// The basin this forcing samples its wind over.
    fn basin(&self) -> Basin;

    /// The stress field at `t_s` seconds, for `state`.
    fn at(&mut self, t_s: f64, state: &OceanState) -> &WindStressField;
}

/// A wind and a field the caller already owns, for the length of one step.
///
/// What [`Solver::step_forced_by`](crate::Solver::step_forced_by) is left with:
/// it is handed the wind as an argument, so nothing it caches can be trusted
/// past the end of the call, and the field it samples into is the solver's own
/// buffer rather than one of this type's. Within the call it is the same rule
/// as [`WindForcing`]'s — which is worth one sampling per step instead of four
/// for a steady wind, and three instead of four for a wind that varies.
pub(crate) struct BorrowedForcing<'a, W: WindStress + ?Sized> {
    /// The basin the wind is sampled over.
    basin: Basin,
    /// The wind this forcing is of.
    wind: &'a W,
    /// The caller's buffer, which this samples into.
    field: &'a mut WindStressField,
    /// Which instant `field` is the field of, and whether it can go stale.
    held: HeldInstant,
}

impl<'a, W: WindStress + ?Sized> BorrowedForcing<'a, W> {
    /// `wind` over `basin`, sampling into `field`, holding nothing yet.
    pub(crate) fn new(basin: Basin, wind: &'a W, field: &'a mut WindStressField) -> Self {
        Self {
            basin,
            wind,
            field,
            held: HeldInstant::nothing(wind.time_dependence()),
        }
    }
}

impl<W: WindStress + ?Sized> StageForcing for BorrowedForcing<'_, W> {
    fn basin(&self) -> Basin {
        self.basin
    }

    /// A prescribed wind is a pure function of `(x, y, t)`, so the state the
    /// stage carries is none of its business.
    fn at(&mut self, t_s: f64, _state: &OceanState) -> &WindStressField {
        self.held.ensure(self.field, self.basin, self.wind, t_s)
    }
}

/// A wind and the field it was last sampled into: the forcing of a run.
///
/// [`WindStressField`] is a stress at one instant; this is the pair a time
/// loop actually holds — the [`WindStress`] itself, the basin it is sampled
/// over, and the field, re-sampled only when the wind can have changed since
/// it was written.
///
/// # Why it owns the wind
///
/// T-10.5's cache is only correct as long as the field belongs to the wind
/// that produced it, and a cache that is *told* which wind it belongs to on
/// every call has to trust the caller not to change its mind. Owning the wind
/// removes the question: a `WindForcing` cannot outlive, be re-pointed at, or
/// disagree with the wind it was built from, so "is this field this wind's?"
/// is answered by the type rather than by a comment. A different wind is a
/// different `WindForcing`.
///
/// Sampling writes into a buffer this owns, so a whole run allocates its
/// forcing exactly once (CODING_STANDARDS.md § Performance), as
/// [`Solver::step_forced_by`](crate::Solver::step_forced_by) already did — what
/// is new is that it does not *rewrite* it four times a step for a wind that
/// is not going anywhere.
#[derive(Debug)]
pub struct WindForcing<W: WindStress> {
    /// The basin the wind is sampled over.
    basin: Basin,
    /// The wind this forcing is of.
    wind: W,
    /// The stress at the instant `held` names.
    field: WindStressField,
    /// Which instant `field` is the field of, and whether it can go stale.
    held: HeldInstant,
}

impl<W: WindStress> StageForcing for WindForcing<W> {
    fn basin(&self) -> Basin {
        Self::basin(self)
    }

    /// A prescribed wind is a pure function of `(x, y, t)`, so the state the
    /// stage carries is none of its business.
    fn at(&mut self, t_s: f64, _state: &OceanState) -> &WindStressField {
        Self::at(self, t_s)
    }
}

impl<W: WindStress> WindForcing<W> {
    /// `wind` over `basin`, with nothing sampled yet.
    ///
    /// The field is allocated here and reused for the life of the forcing; the
    /// first call to [`WindForcing::at`] fills it.
    #[must_use]
    pub fn new(basin: Basin, wind: W) -> Self {
        let dependence = wind.time_dependence();
        Self {
            basin,
            field: WindStressField::calm(basin.grid()),
            wind,
            held: HeldInstant::nothing(dependence),
        }
    }

    /// The stress field at `t_s` seconds, sampling the wind only if the field
    /// in hand is not already that field.
    pub fn at(&mut self, t_s: f64) -> &WindStressField {
        self.held
            .ensure(&mut self.field, self.basin, &self.wind, t_s)
    }

    /// The wind this forcing is of.
    pub const fn wind(&self) -> &W {
        &self.wind
    }

    /// The basin the wind is sampled over.
    pub const fn basin(&self) -> Basin {
        self.basin
    }
}

/// Add one component of `wind` to `component`, at every face it has.
///
/// [`write_component`] with `+=` for `=`: the superposition
/// [`WindStressField::add_sampled`] performs.
fn add_component<W, Pick>(
    component: &mut Field2D<f64>,
    basin: Basin,
    staggering: Staggering,
    pick: Pick,
    wind: &W,
    t_s: f64,
) where
    W: WindStress + ?Sized,
    Pick: Fn((f64, f64)) -> f64,
{
    for j in 0..component.ny() {
        let y_m = basin.y_of_row_m(staggering, j);
        for i in 0..component.nx() {
            let value = pick(wind.stress(basin.x_of_column_m(staggering, i), y_m, t_s));
            *component
                .get_mut(i, j)
                .expect("the loop bounds are the field's own extents") += value;
        }
    }
}

/// Write one component of `wind` into `component`, at every face it has.
///
/// Shared by the two halves of [`WindStressField::sample`], which differ only
/// in where their points sit and which half of the returned pair they keep.
fn write_component<W, Pick>(
    component: &mut Field2D<f64>,
    basin: Basin,
    staggering: Staggering,
    pick: Pick,
    wind: &W,
    t_s: f64,
) where
    W: WindStress + ?Sized,
    Pick: Fn((f64, f64)) -> f64,
{
    for j in 0..component.ny() {
        let y_m = basin.y_of_row_m(staggering, j);
        for i in 0..component.nx() {
            let value = pick(wind.stress(basin.x_of_column_m(staggering, i), y_m, t_s));
            *component
                .get_mut(i, j)
                .expect("the loop bounds are the field's own extents") = value;
        }
    }
}
