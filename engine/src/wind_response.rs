//! The atmosphere's answer to the SST anomaly: the wind half of the Bjerknes
//! feedback, and the piece that closes the loop.
//!
//! T-12.1 made `T'` respond to the ocean. Nothing yet made the ocean respond
//! to `T'`, so the coupling ran one way and the wind was still whatever the
//! scenario prescribed. This module is the other direction — `CONTEXT.md`,
//! *Bjerknes feedback*: **weaker trade winds → flatter thermocline → warmer
//! eastern SST → weaker trade winds**, and it is the last arrow of that chain
//! that lives here.
//!
//! # The response
//!
//! ```text
//! τx'(x, y, t) = μ · ⟨T'⟩(t) · exp(−(y/L_a)²)        τy'(x, y, t) = 0
//! ```
//!
//! - `μ` is the **feedback strength**, in Pa/K, and it is a scenario
//!   parameter rather than a constant: `μ = 0` is the prescribed-wind model of
//!   Epics 01–07 with T-12.1's SST equation riding along, and turning it up is
//!   what T-12.3 does to find the oscillation. It is positive because the
//!   feedback is: a *warm* anomaly *weakens* the easterly alizés, so it adds a
//!   westerly (`τx > 0`) anomaly to them.
//! - `⟨T'⟩` is the SST anomaly projected onto the same equatorial Gaussian —
//!   one number for the whole basin, recomputed from the state at every
//!   right-hand-side evaluation.
//! - `L_a` is the meridional scale of the atmospheric response.
//!
//! # Why this shape, and what it is not
//!
//! It is the *statistical* atmosphere of `docs/planning/01-scientific-model.md`
//! § *Phase 2*, with a Gill-type spatial pattern — the shape the delayed
//! oscillator literature runs on (Suarez & Schopf, *J. Atmos. Sci.* 45, 1988;
//! Battisti & Hirst, *J. Atmos. Sci.* 46, 1989), where the wind anomaly is a
//! regression of a fixed pattern onto an SST index. It is **not** a solution of
//! the atmospheric equations: no Gill model is integrated here, and nothing
//! answers what the wind is doing away from the equatorial waveguide.
//!
//! Two properties of the real thing justify the shape:
//!
//! - **The atmosphere is fast.** A tropical atmospheric adjustment takes days;
//!   the ocean's basin-crossing waves take months. So the wind anomaly is
//!   *diagnostic* — a function of the SST of the instant, with no memory —
//!   which is what lets it be a [`WindStress`] at all rather than a fifth
//!   prognostic variable.
//! - **The response is equatorially trapped.** Gill (*Q. J. R. Meteorol. Soc.*
//!   106, 1980) § 2: the zonal wind of the Kelvin part of the heating response
//!   falls off as `exp(−βy²/(2·c_a))`, which is `exp(−(y/L_a)²)` with
//!   `L_a = √(2·c_a/β)` — the atmospheric equatorial Rossby radius, and the
//!   same Gaussian-about-the-equator form
//!   [`SteadyTradeWinds`](crate::SteadyTradeWinds) and
//!   [`WindBurstAnomaly`](crate::WindBurstAnomaly) already use.
//!
//! The response carries no zonal structure. Gill's solution does — the
//! westerly anomaly sits west of the heating and the easterly east of it — but
//! reproducing that needs the atmospheric model this deliberately does not
//! solve, and the delayed oscillator does not turn on it: the delay that makes
//! the loop oscillate rather than run away is the ocean's, an off-equatorial
//! Rossby wave crossing to the western wall and reflecting back as a Kelvin
//! wave, and a zonally uniform stress anomaly excites that just as a patch
//! does. A scenario wanting a localised anomaly already has
//! [`WindBurstAnomaly`](crate::WindBurstAnomaly) to add.
//!
//! # Why the index, and not the local anomaly
//!
//! `τx` is a function of the whole SST field through one scalar, rather than
//! of `T'` at the same point. That is the *statistical* half of "statistical
//! or Gill-type": the atmosphere integrates its heating over the basin, so a
//! warm patch anywhere on the equator moves the wind everywhere along it.
//! Making the stress local would also make it grid-dependent — `τx` lives on
//! the east/west faces and `T'` at the cell centers, so a pointwise
//! [`WindStress`] would have to interpolate a field it is not given — and a
//! wind that reads the grid is no longer the pure function of `(x, y, t)` the
//! trait is.
//!
//! # Purity, and where the state comes in
//!
//! [`WindStress`] requires a pure function of `(x, y, t)`, and it means it:
//! RK4 samples the forcing four times a step and the integration must not
//! depend on the order. [`SstWindResponse`] keeps that contract. Between two
//! calls to [`SstWindResponse::observe`] it *is* a pure function — the index is
//! a number it holds — and `observe` is called exactly once per stage, from
//! [`CoupledWind::at`], with the state of *that* stage. So each stage's
//! right-hand side is a function of that stage's state and time, which is what
//! an ODE right-hand side is; no stage sees another stage's wind.
//!
//! That is also why the response is not put inside a
//! [`WindForcing`](crate::WindForcing) beside the prescribed winds: that cache
//! reuses a field whenever the *instant* repeats, and RK4 asks about
//! `t + dt/2` twice with two different states. [`CoupledWind`] keeps the two
//! apart — the prescribed forcing cached on time as T-10.5 left it, the
//! response re-sampled every stage — and adds them, which is the same
//! superposition [`CompositeWind`](crate::CompositeWind) performs and the same
//! one T-03.3 established.
//!
//! [ADR-0010] records that decision, the alternatives weighed against it, and
//! why the feedback strength is a parameter of `[sst]` rather than a
//! `[[wind]]` entry.
//!
//! [ADR-0010]: ../../docs/planning/adr/0010-wind-response-is-diagnosed-per-stage.md

use std::fmt;

use termocline_grid::{Field2D, Grid, H_STAGGERING};

use crate::basin::Basin;
use crate::forcing::{
    gaussian, StageForcing, TimeDependence, WindForcing, WindStress, WindStressField,
};
use crate::state::OceanState;

/// Value the meridional stress of this response takes everywhere.
///
/// The statistical response is zonal only: what forces the equatorial ocean,
/// and what the Bjerknes loop runs on, is `τx`.
const NO_MERIDIONAL_STRESS: f64 = 0.0;

/// Equivalent gravity-wave speed `c_a` of the first baroclinic mode of the
/// tropical atmosphere, in m/s.
///
/// Gill (*Q. J. R. Meteorol. Soc.* 106, 1980), whose heating-response solutions
/// are scaled by it; 60 m/s is the value for an equivalent depth of about
/// 370 m, the first internal mode of the tropical troposphere.
pub const ATMOSPHERIC_GRAVITY_WAVE_SPEED_M_PER_S: f64 = 60.0;

/// Meridional scale `L_a` of the atmospheric response, in metres, when a
/// scenario does not state one.
///
/// The atmospheric equatorial Rossby radius `√(2·c_a/β)`, from
/// [`ATMOSPHERIC_GRAVITY_WAVE_SPEED_M_PER_S`] and the `β` of `CONTEXT.md`
/// (2.3 × 10⁻¹¹ m⁻¹s⁻¹): 2.3 × 10⁶ m, quoted to two significant figures.
///
/// It is far wider than the *ocean's* deformation radius — the atmosphere's
/// gravity waves are two orders of magnitude faster than the ocean's — so the
/// wind anomaly this model puts on the basin is broad, covering the whole
/// equatorial waveguide rather than a band inside it. That is the physics
/// rather than a simplification: an equatorial heating anomaly moves the
/// tropical winds over tens of degrees of latitude.
pub const DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M: f64 = 2.3e6;

/// Why a set of wind-response parameters was rejected.
///
/// These describe invalid *scenario input* — an `[sst]` section asking for an
/// atmosphere that does not exist — so they are returned rather than panicked,
/// and each names the offending parameter and the value it carried
/// (CODING_STANDARDS.md § *Correctness and failure*).
#[derive(Debug, Clone, PartialEq)]
pub enum WindResponseError {
    /// A parameter that must be non-negative and finite was negative, or not a
    /// number.
    Negative {
        /// Name of the parameter, matching the config key.
        parameter: &'static str,
        /// The value supplied, in the unit the parameter's name states.
        value: f64,
    },
    /// A parameter that must be strictly positive and finite was not.
    NotPositive {
        /// Name of the parameter, matching the config key.
        parameter: &'static str,
        /// The value supplied, in the unit the parameter's name states.
        value: f64,
    },
}

impl fmt::Display for WindResponseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Negative { parameter, value } => write!(
                f,
                "{parameter} is {value}; it must be finite and at least 0"
            ),
            Self::NotPositive { parameter, value } => write!(
                f,
                "{parameter} is {value}; it must be finite and greater than 0"
            ),
        }
    }
}

impl std::error::Error for WindResponseError {}

/// The constants of one scenario's atmospheric response: how hard the wind
/// answers the SST anomaly, and over how much latitude.
///
/// Constructed once per run and read at every stage, so it is `Copy` and
/// validated at the boundary rather than at each use — the same shape as
/// [`SstParams`](crate::SstParams).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindResponseParams {
    /// Feedback strength `μ`, in Pa/K.
    feedback_strength_pa_per_k: f64,
    /// Meridional scale `L_a` of the response, in metres.
    meridional_scale_m: f64,
}

impl WindResponseParams {
    /// The atmospheric response's parameter set, in SI units.
    ///
    /// `μ` may be zero, which is the switch this ticket's acceptance criterion
    /// turns on: a run at zero feedback strength is the prescribed-wind run,
    /// bit for bit. It may never be negative — a negative `μ` would make a
    /// warm anomaly *strengthen* the alizés, which is the Bjerknes feedback
    /// run backwards and not a state of the equatorial Pacific. `L_a` is
    /// divided by, and a zero-width atmosphere is not a scenario, so it must
    /// be strictly positive.
    ///
    /// # Errors
    /// A [`WindResponseError`] naming the first parameter that failed and the
    /// value it carried.
    pub fn new(
        feedback_strength_pa_per_k: f64,
        meridional_scale_m: f64,
    ) -> Result<Self, WindResponseError> {
        if !feedback_strength_pa_per_k.is_finite() || feedback_strength_pa_per_k < 0.0 {
            return Err(WindResponseError::Negative {
                parameter: "wind_feedback_strength_pa_per_k",
                value: feedback_strength_pa_per_k,
            });
        }
        if !meridional_scale_m.is_finite() || meridional_scale_m <= 0.0 {
            return Err(WindResponseError::NotPositive {
                parameter: "wind_response_meridional_scale_m",
                value: meridional_scale_m,
            });
        }
        Ok(Self {
            feedback_strength_pa_per_k,
            meridional_scale_m,
        })
    }

    /// Feedback strength `μ`, in Pa/K. Zero is a prescribed-wind run.
    #[must_use]
    pub const fn feedback_strength_pa_per_k(self) -> f64 {
        self.feedback_strength_pa_per_k
    }

    /// Meridional scale `L_a` of the response, in metres.
    #[must_use]
    pub const fn meridional_scale_m(self) -> f64 {
        self.meridional_scale_m
    }
}

/// The statistical, Gill-type atmospheric wind response to the SST anomaly:
/// the wind half of the Bjerknes feedback.
///
/// A [`WindStress`] like any other — it stacks in a
/// [`CompositeWind`](crate::CompositeWind) beside the trades exactly as the
/// burst of T-03.3 does — with one thing the prescribed scenarios do not have:
/// it must be shown the SST anomaly of the stage being evaluated, through
/// [`SstWindResponse::observe`], before it is sampled. Between two such calls
/// it is the pure function of `(x, y, t)` the trait requires; the module header
/// is why that is enough, and [`CoupledWind`] is what does the showing.
///
/// **Whoever holds one owes it that call.** Dropped into a
/// [`CompositeWind`](crate::CompositeWind) it sums like any other forcing, but
/// a composite cannot reach inside a boxed component to refresh it, so a
/// response nobody observes serves one frozen index for a whole run. A run
/// therefore holds it through [`CoupledWind`], which owns both the response and
/// the refresh; ADR-0010 records why that is the shape rather than a
/// `[[wind]]` entry.
///
/// Built once per run: the row weights of its index are allocated here and
/// reused at every stage (CODING_STANDARDS.md § *Performance*).
#[derive(Debug, Clone, PartialEq)]
pub struct SstWindResponse {
    /// How hard, and over how much latitude, the atmosphere answers.
    params: WindResponseParams,
    /// The basin the SST anomaly this reads is expected to cover.
    grid: Grid,
    /// `exp(−(yⱼ/L_a)²)` at each cell-center row — the weight that row's SST
    /// anomaly carries in the index. One per row, because the weight is a
    /// function of `y` alone.
    row_weights: Vec<f64>,
    /// `Σᵢⱼ exp(−(yⱼ/L_a)²)` over every cell, the divisor that makes the index
    /// a weighted *mean* and so a temperature rather than a sum of them.
    weight_sum: f64,
    /// The index `⟨T'⟩`, in kelvin, as of the last
    /// [`SstWindResponse::observe`]. Zero before the first, which is the
    /// climatology — the state a coupled run starts from.
    index_k: f64,
}

impl SstWindResponse {
    /// The response of the atmosphere over `basin`, with `params`.
    ///
    /// The index weights are computed here, from the latitudes of `basin`'s
    /// cell-center rows, and never again.
    #[must_use]
    pub fn new(basin: Basin, params: WindResponseParams) -> Self {
        let grid = basin.grid();
        let row_weights: Vec<f64> = (0..grid.ny())
            .map(|j| {
                gaussian(
                    basin.y_of_row_m(H_STAGGERING, j),
                    params.meridional_scale_m(),
                )
            })
            .collect();
        // Every row carries `nx` cells, so the total weight is the row weights
        // summed once and scaled, rather than `nx · ny` additions of the same
        // few numbers.
        let weight_sum = row_weights.iter().sum::<f64>() * grid.nx() as f64;
        Self {
            params,
            grid,
            row_weights,
            weight_sum,
            index_k: 0.0,
        }
    }

    /// Recompute the index `⟨T'⟩` from `sst_anomaly_k`.
    ///
    /// Called once per right-hand-side evaluation, with the SST anomaly of the
    /// stage about to be evaluated: what this reads is what that stage's wind
    /// answers.
    ///
    /// # Panics
    /// If `sst_anomaly_k` covers a different basin from the one this response
    /// was built for. A shape mismatch means the calling code is wrong, which
    /// is what panics are for (CODING_STANDARDS.md § *Correctness and
    /// failure*).
    pub fn observe(&mut self, sst_anomaly_k: &Field2D<f64>) {
        let (nx, ny) = (self.grid.nx(), self.grid.ny());
        assert!(
            sst_anomaly_k.nx() == nx && sst_anomaly_k.ny() == ny,
            "the SST anomaly is {} by {}, but this wind response was built for {:?}",
            sst_anomaly_k.nx(),
            sst_anomaly_k.ny(),
            self.grid
        );
        // Rows summed first and weighted once, in row order, so the sum's
        // floating-point result is fixed by the grid rather than by an
        // iteration order (CODING_STANDARDS.md § *Correctness and failure*).
        let weighted_k: f64 = sst_anomaly_k
            .as_slice()
            .chunks_exact(nx)
            .zip(&self.row_weights)
            .map(|(row, weight)| weight * row.iter().sum::<f64>())
            .sum();
        self.index_k = weighted_k / self.weight_sum;
    }

    /// The index `⟨T'⟩`, in kelvin, as of the last
    /// [`SstWindResponse::observe`].
    ///
    /// The SST anomaly projected onto the atmosphere's equatorial mode: the
    /// one number the whole basin's wind anomaly is scaled by, and what
    /// T-12.3's oscillation is looked for in.
    #[must_use]
    pub const fn index_k(&self) -> f64 {
        self.index_k
    }

    /// The basin whose SST anomaly this response reads.
    #[must_use]
    pub const fn grid(&self) -> Grid {
        self.grid
    }
}

impl WindStress for SstWindResponse {
    /// `μ · ⟨T'⟩ · exp(−(y/L_a)²)`, and no meridional stress.
    ///
    /// Independent of `x`, because a statistical atmosphere answers the basin's
    /// heating rather than the heating underneath it, and of `t`, because the
    /// only time in it is the index — which belongs to the stage, not to the
    /// clock.
    fn stress(&self, _x_m: f64, y_m: f64, _t_s: f64) -> (f64, f64) {
        (
            self.params.feedback_strength_pa_per_k()
                * self.index_k
                * gaussian(y_m, self.params.meridional_scale_m()),
            NO_MERIDIONAL_STRESS,
        )
    }

    /// [`TimeDependence::Varying`], which is the answer that is always safe.
    ///
    /// The field this response produces genuinely does change from one stage
    /// to the next — the index does — even though `t` is not what changes it,
    /// so a cache keyed on the instant must never hold one of its fields.
    /// [`CoupledWind`] does not put it in one; this declaration is what makes a
    /// scenario that stacks the response into a
    /// [`CompositeWind`](crate::CompositeWind) safe too (ADR-0009).
    fn time_dependence(&self) -> TimeDependence {
        TimeDependence::Varying
    }
}

/// The forcing of a coupled run: the prescribed winds plus the atmosphere's
/// answer to the SST anomaly.
///
/// The two are summed, which is the superposition T-03.3 established — a
/// response is an anomaly *superimposed on* the trades, and the equations are
/// linear in the stress. They are held apart rather than in one
/// [`CompositeWind`](crate::CompositeWind) because they invalidate differently:
/// the prescribed half is a pure function of time and keeps T-10.5's cache, so
/// a steady scenario samples it once for the run; the response half is
/// re-observed and re-sampled at every stage, because that is what makes it
/// answer the state of *that* stage.
///
/// Everything it writes into is allocated here and reused, so a coupled run
/// allocates its forcing exactly once (CODING_STANDARDS.md § *Performance*).
#[derive(Debug)]
pub struct CoupledWind<W: WindStress> {
    /// The scenario's prescribed winds, cached on the instant as T-10.5 left
    /// them.
    prescribed: WindForcing<W>,
    /// The atmosphere's answer to the SST anomaly.
    response: SstWindResponse,
    /// The sum of the two: the stress a stage actually reads.
    total: WindStressField,
}

impl<W: WindStress> CoupledWind<W> {
    /// `wind` over `basin`, with `response` added to it.
    ///
    /// # Panics
    /// If `response` was built for a different basin from `basin`. A shape
    /// mismatch means the calling code is wrong (CODING_STANDARDS.md
    /// § *Correctness and failure*).
    #[must_use]
    pub fn new(basin: Basin, wind: W, response: SstWindResponse) -> Self {
        assert!(
            response.grid() == basin.grid(),
            "the wind response covers {:?}, but this forcing is over {:?}",
            response.grid(),
            basin.grid()
        );
        Self {
            prescribed: WindForcing::new(basin, wind),
            response,
            total: WindStressField::calm(basin.grid()),
        }
    }

    /// The stress field at `t_s` seconds for `state`: the prescribed winds of
    /// that instant plus the atmosphere's answer to that state's SST anomaly.
    ///
    /// # Panics
    /// If `state` is not a coupled state. A response with nothing to respond
    /// to means the run was assembled wrong: the `[sst]` section is what
    /// switches both the SST anomaly and this forcing on, so they arrive
    /// together or not at all.
    pub fn at(&mut self, t_s: f64, state: &OceanState) -> &WindStressField {
        let Self {
            prescribed,
            response,
            total,
        } = self;
        response.observe(
            state
                .sst_anomaly_k()
                .expect("a coupled run's wind response reads a coupled state"),
        );
        let basin = prescribed.basin();
        total.assign(prescribed.at(t_s));
        total.add_sampled(basin, response, t_s);
        total
    }

    /// The basin the two halves are sampled over.
    #[must_use]
    pub const fn basin(&self) -> Basin {
        self.prescribed.basin()
    }
}

impl<W: WindStress> StageForcing for CoupledWind<W> {
    fn basin(&self) -> Basin {
        Self::basin(self)
    }

    fn at(&mut self, t_s: f64, state: &OceanState) -> &WindStressField {
        Self::at(self, t_s, state)
    }
}
