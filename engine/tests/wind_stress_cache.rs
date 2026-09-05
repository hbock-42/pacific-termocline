//! Acceptance tests for T-10.5 — caching the sampled wind-stress field.
//!
//! `docs/performance-notes.md` measured re-sampling the wind at 71% of a
//! timestep: every C-grid face, at every one of RK4's four stages, through two
//! virtual calls ending in a libm `exp`. For the control scenario every one of
//! those evaluations returns what the one before it returned, because
//! [`SteadyTradeWinds`] does not depend on `t` at all.
//!
//! The ticket's own words are the reason this file exists: **a cache that is
//! correct only for steady winds is a bug, not an optimisation.** So the
//! tests here are in two halves, and the second is the load-bearing one.
//!
//! The first half is that the cache is *taken*: a steady wind is sampled once
//! for a whole run, and a time-varying one is sampled once per *distinct*
//! stage time rather than once per stage. Both are counted at the
//! [`WindStress`] itself, by a wrapper that tallies the calls, so what is
//! asserted is the number of evaluations the ticket is about rather than a
//! duration — a wall-clock threshold in the gate is the flaky measurement
//! `docs/benchmarks.md` exists to avoid, and the speed-up is cited from
//! `cargo bench` in the note instead.
//!
//! The second half is that the cache is *correct*, for the winds that
//! genuinely vary in time. [`SeasonalTradeWinds`] breathes with the year and
//! [`WindBurstAnomaly`] is a Gaussian in `t`; a `CompositeWind` stacks either
//! on the trades. Against those, the cached path is compared with an
//! uncached reference stepper written out in this file — the same four
//! evaluators in the same order, sampling the wind afresh at every stage —
//! and the comparison is **bit for bit**, because both perform the same
//! operations in the same order over the same `f64` fields. A tolerance here
//! would be a place for a real divergence to hide, which is
//! `tests/step_profile.rs`'s reasoning for comparing profiler to solver the
//! same way.
//!
//! The stage times are the RK4 tableau's own, from `integrator.rs`: `t`,
//! `t + dt/2`, `t + dt/2` and `t + dt`. Four stages, three distinct instants —
//! which is where the ceiling on a time-varying wind's saving comes from, and
//! it is derived here from the tableau rather than read off a run.

use std::cell::Cell;

use engine::{
    Basin, BetaPlane, CompositeWind, CoriolisTerm, Grid, NoNormalFlow, OceanState, PhysicalParams,
    Rk4, SeasonalTradeWinds, ShallowWaterRhs, Solver, Spacing, SteadyTradeWinds, TimeDependence,
    WindBurstAnomaly, WindForcing, WindStress, WindStressField, TROPICAL_YEAR_S,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981).
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;
/// Rayleigh damping rate `r`, in s⁻¹ — a 100-day decay, the scale
/// `docs/planning/01-scientific-model.md` gives for the linear drag.
const PACIFIC_DAMPING_PER_S: f64 = 1.0 / (100.0 * DAY_S);

/// One solar day, in seconds.
const DAY_S: f64 = 86_400.0;

/// Zonal stress `τ₀` of the alizés, in Pa — easterly, and the observed scale
/// of the equatorial Pacific's mean zonal stress.
const TRADE_WIND_STRESS_PA: f64 = -0.05;
/// Meridional decay scale of the alizés, in metres — a few degrees of
/// latitude, the width of the belt the trades occupy.
const TRADE_WIND_DECAY_M: f64 = 5.0e5;

/// Peak zonal stress of the burst, in Pa — westerly, so positive, and the
/// observed scale of an equatorial westerly wind burst.
const BURST_STRESS_PA: f64 = 0.04;
/// Zonal centre of the burst, in metres east of the western boundary.
const BURST_CENTER_X_M: f64 = 2.0e6;
/// Zonal `e`-folding scale of the burst, in metres — about 5° of longitude.
const BURST_ZONAL_SCALE_M: f64 = 5.0e5;
/// Meridional `e`-folding scale of the burst, in metres.
const BURST_MERIDIONAL_SCALE_M: f64 = 3.45e5;
/// Instant of the burst's peak, in seconds — 3 days in, so a short comparison
/// run climbs the leading flank of the Gaussian rather than sitting on a
/// plateau where a broken cache could pass unnoticed.
const BURST_PEAK_TIME_S: f64 = 3.0 * DAY_S;
/// Temporal `e`-folding scale of the burst, in seconds.
const BURST_DURATION_S: f64 = 5.0 * DAY_S;

/// Relative amplitude `a` of the seasonal harmonic, dimensionless — a strong
/// but physical season, and strictly inside `[0, 1]` so the modulation is
/// never stationary at an extremum for the length of a comparison run.
const SEASONAL_AMPLITUDE: f64 = 0.4;
/// Instant the modulated alizés are strongest, in seconds.
const SEASONAL_PEAK_TIME_S: f64 = 60.0 * DAY_S;

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*).
const BASIN_LX_M: f64 = 1.0e7;
/// Meridional extent of the test basin, in metres.
const BASIN_LY_M: f64 = 2.0e6;
/// Zonal cells of the test basin. Small: these tests are about how many times
/// the wind is evaluated and about whether two steppers agree, neither of
/// which needs resolution.
const BASIN_NX: usize = 24;
/// Meridional cells of the test basin.
const BASIN_NY: usize = 10;

/// Timestep of the comparison runs, in seconds. Well inside the gravity-wave
/// CFL bound of a basin this coarse, and short beside every scale of the
/// forcing, so the stage times of successive steps are distinct instants of
/// the burst's Gaussian rather than one instant repeated.
const DT_S: f64 = 900.0;

/// Steps a comparison run takes.
///
/// Enough that a stale cached stress would have propagated: the shallow-water
/// terms couple `h`, `u` and `v` at every step, so a wrong stress at any stage
/// reaches all three within a handful of steps and grows. Sixteen steps span
/// four hours, over which the burst's envelope changes by about 1%, so a cache
/// that never invalidated would be wrong by far more than the last bit.
const COMPARED_STEPS: u64 = 16;

/// Node of the two midpoint stages of the RK4 tableau, as a fraction of `dt`
/// (`integrator.rs`). The stage times are `t`, `t + dt/2`, `t + dt/2`,
/// `t + dt`.
const MIDPOINT_NODE: f64 = 0.5;

/// Relative slack allowed where a check is exact in exact arithmetic: a few
/// tens of ulps of `f64` (ε ≈ 2.2×10⁻¹⁶) for the handful of operations per
/// point the closed form costs.
const ROUNDING_TOLERANCE: f64 = 1.0e-14;

// --- Fixtures. ---

/// The equatorial-Pacific parameter set of the scenarios above.
fn pacific_params() -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        PACIFIC_DAMPING_PER_S,
        engine::EQUATORIAL_BETA_PER_M_PER_S,
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

/// A basin [`BASIN_LX_M`] by [`BASIN_LY_M`] metres centred on the equator.
fn test_basin() -> Basin {
    let grid = Grid::new(BASIN_NX, BASIN_NY).expect("extents are non-zero");
    let spacing = Spacing::new(BASIN_LX_M / BASIN_NX as f64, BASIN_LY_M / BASIN_NY as f64)
        .expect("a basin spanned by whole cells has positive spacing");
    Basin::centered_on_equator(grid, spacing)
}

/// The control scenario's forcing: steady alizés decaying away from the
/// equator.
fn trade_winds() -> SteadyTradeWinds {
    SteadyTradeWinds::with_meridional_decay(TRADE_WIND_STRESS_PA, TRADE_WIND_DECAY_M)
        .expect("an easterly stress with a positive decay scale is a trade wind")
}

/// The seasonal cycle: the same alizés breathing with the year.
fn seasonal_winds() -> SeasonalTradeWinds {
    SeasonalTradeWinds::new(trade_winds(), SEASONAL_AMPLITUDE, SEASONAL_PEAK_TIME_S)
        .expect("a fractional amplitude and a finite phase are a season")
}

/// The westerly wind burst: a Gaussian in `x`, in `y` and in `t`.
fn burst() -> WindBurstAnomaly {
    WindBurstAnomaly::new(
        BURST_STRESS_PA,
        BURST_CENTER_X_M,
        BURST_ZONAL_SCALE_M,
        BURST_MERIDIONAL_SCALE_M,
        BURST_PEAK_TIME_S,
        BURST_DURATION_S,
    )
    .expect("a westerly burst with positive scales and a positive duration")
}

/// The hardest case the ticket names: a burst superimposed on a seasonal
/// cycle, so both a scalar harmonic and a travelling-in-time Gaussian are in
/// the same field.
fn seasonal_winds_with_a_burst() -> CompositeWind {
    CompositeWind::new().with(seasonal_winds()).with(burst())
}

/// A [`WindStress`] that counts how many times it is asked for a stress.
///
/// The instrument of the first half of this file. It forwards both the stress
/// and the *declared* time dependence, so wrapping a wind changes how often it
/// is called and nothing else.
struct CountingWind<W: WindStress> {
    /// The wind being counted.
    inner: W,
    /// Calls to [`WindStress::stress`] so far. A `Cell` because the trait
    /// takes `&self` — it must, since a stateful wind would make RK4 depend on
    /// stage order (`forcing.rs`).
    calls: Cell<u64>,
}

impl<W: WindStress> CountingWind<W> {
    /// `inner`, with a call counter at zero.
    fn new(inner: W) -> Self {
        Self {
            inner,
            calls: Cell::new(0),
        }
    }

    /// How many times the wrapped wind has been evaluated.
    fn calls(&self) -> u64 {
        self.calls.get()
    }
}

impl<W: WindStress> WindStress for CountingWind<W> {
    fn stress(&self, x_m: f64, y_m: f64, t_s: f64) -> (f64, f64) {
        self.calls.set(self.calls.get() + 1);
        self.inner.stress(x_m, y_m, t_s)
    }

    fn time_dependence(&self) -> TimeDependence {
        self.inner.time_dependence()
    }
}

/// A solver for [`test_basin`] at [`DT_S`].
fn test_solver(basin: Basin) -> Solver {
    let params = pacific_params();
    Solver::new(
        basin.grid(),
        basin.spacing(),
        params,
        BetaPlane::of_basin(params, basin),
        DT_S,
    )
    .expect("a coarse basin at a quarter-hour timestep clears both bounds")
}

/// How many wind evaluations one full sampling pass over `basin` costs: every
/// east/west face for `τx`, every north/south face for `τy`.
///
/// Counted from the field's own extents rather than from `(nx+1)·ny +
/// nx·(ny+1)` written out here, so the test does not keep a second copy of
/// where the C-grid puts a stress.
fn evaluations_per_sampling(basin: Basin) -> u64 {
    let field = WindStressField::calm(basin.grid());
    let tau_x = field.tau_x_pa();
    let tau_y = field.tau_y_pa();
    (tau_x.nx() * tau_x.ny() + tau_y.nx() * tau_y.ny()) as u64
}

/// The stage times a run of `steps` steps from `t = 0` asks the wind about, in
/// order, with the repeats the RK4 tableau produces.
///
/// The tableau is `integrator.rs`'s: `t`, `t + dt/2`, `t + dt/2`, `t + dt` for
/// each step, and a step starts where the schedule says rather than where the
/// previous step's arithmetic landed — which is what
/// [`run_scenario`](engine::run_scenario) does, and it is reproduced here so
/// that the expected count below is derived from the tableau rather than from
/// a run.
fn stage_times(steps: u64) -> Vec<f64> {
    let mut times = Vec::with_capacity(steps as usize * 4);
    for step in 0..steps {
        let t_s = step as f64 * DT_S;
        times.push(t_s);
        times.push(t_s + MIDPOINT_NODE * DT_S);
        times.push(t_s + MIDPOINT_NODE * DT_S);
        times.push(t_s + DT_S);
    }
    times
}

/// How many of `times` differ from the one before them — the number of
/// samplings a cache holding exactly one instant can be asked for.
fn distinct_consecutive(times: &[f64]) -> u64 {
    let mut held: Option<f64> = None;
    let mut samplings = 0;
    for &t_s in times {
        if held != Some(t_s) {
            samplings += 1;
            held = Some(t_s);
        }
    }
    samplings
}

/// One step of the engine's own scheme, with the wind sampled afresh at every
/// stage: the uncached reference the cached path is held to.
///
/// The evaluators, their order and the boundary condition are
/// [`Solver::step_forced_by`]'s, written out here so that the comparison is
/// against a stepper that has no cache at all rather than against another
/// caller of the same cache. `tests/step_profile.rs` holds
/// `engine::profiling`'s open-coded step to the solver in the same way and for
/// the same reason.
struct UncachedStepper {
    /// The basin the wind is sampled over.
    basin: Basin,
    /// The pressure-gradient, continuity, surface-stress and damping terms.
    rhs: ShallowWaterRhs,
    /// The beta-plane rotation terms.
    coriolis: CoriolisTerm,
    /// The integrator and its stage buffers.
    integrator: Rk4<OceanState>,
}

impl UncachedStepper {
    /// A stepper for `basin` with the parameters of [`pacific_params`].
    fn new(basin: Basin) -> Self {
        let params = pacific_params();
        Self {
            basin,
            rhs: ShallowWaterRhs::new(basin.grid(), basin.spacing(), params),
            coriolis: CoriolisTerm::new(
                basin.grid(),
                basin.spacing(),
                BetaPlane::of_basin(params, basin),
            ),
            integrator: Rk4::new(&OceanState::at_rest(basin.grid())),
        }
    }

    /// Advance `state` from `t_s` by [`DT_S`] seconds under `wind`.
    fn step<W: WindStress + ?Sized>(&mut self, state: &mut OceanState, t_s: f64, wind: &W) {
        let Self {
            basin,
            rhs,
            coriolis,
            integrator,
        } = self;
        NoNormalFlow::apply_to_state(state);
        integrator.step(
            state,
            t_s,
            DT_S,
            &mut |now: &OceanState, stage_t_s: f64, tendency: &mut OceanState| {
                let stress = WindStressField::sampled(*basin, wind, stage_t_s);
                rhs.evaluate(now, &stress, tendency);
                coriolis.add_to_tendency(now, tendency);
                NoNormalFlow::apply_to_tendency(tendency);
            },
        );
    }
}

// --- What a wind says about its own dependence on time. ---

#[test]
fn steady_trade_winds_declare_themselves_steady() {
    // The whole basis of the cache: `SteadyTradeWinds::stress` ignores `t_s`,
    // so a field sampled from it at one instant is the field at every instant.
    assert_eq!(trade_winds().time_dependence(), TimeDependence::Steady);
}

#[test]
fn a_season_and_a_burst_declare_themselves_time_varying() {
    // The two scenarios the ticket names as the ones a naive cache would
    // break. Both are genuinely functions of `t`, and both must say so.
    assert_eq!(seasonal_winds().time_dependence(), TimeDependence::Varying);
    assert_eq!(burst().time_dependence(), TimeDependence::Varying);
}

#[test]
fn a_season_of_no_amplitude_at_all_is_steady() {
    // `1 + a·cos(…)` with `a = 0` is exactly 1 for every `t`, in exact
    // arithmetic and in `f64` — `0·cos + 1` is `1.0` whatever the cosine is.
    // The declaration is a property of the instance, not of the type.
    let unmodulated = SeasonalTradeWinds::new(trade_winds(), 0.0, SEASONAL_PEAK_TIME_S)
        .expect("a zero amplitude is a fraction");

    assert_eq!(unmodulated.time_dependence(), TimeDependence::Steady);
}

#[test]
fn a_composite_is_steady_only_when_every_component_is() {
    // Superposition: a sum is constant in time exactly when every term is.
    // The empty composite is calm, and calm does not change.
    assert_eq!(
        CompositeWind::new().time_dependence(),
        TimeDependence::Steady
    );
    assert_eq!(
        CompositeWind::new()
            .with(trade_winds())
            .with(trade_winds())
            .time_dependence(),
        TimeDependence::Steady
    );
    assert_eq!(
        CompositeWind::new()
            .with(trade_winds())
            .with(burst())
            .time_dependence(),
        TimeDependence::Varying
    );
    assert_eq!(
        seasonal_winds_with_a_burst().time_dependence(),
        TimeDependence::Varying
    );
}

#[test]
fn a_wind_that_says_nothing_about_time_is_assumed_to_vary() {
    // The default has to be the conservative one: a `WindStress` written
    // outside this crate, by someone who never read this ticket, must be
    // re-sampled rather than silently frozen.
    struct UndeclaredWind;
    impl WindStress for UndeclaredWind {
        fn stress(&self, _x_m: f64, _y_m: f64, t_s: f64) -> (f64, f64) {
            (-t_s, 0.0)
        }
    }

    assert_eq!(UndeclaredWind.time_dependence(), TimeDependence::Varying);
}

// --- That the cache is taken. ---

#[test]
fn a_steady_wind_is_sampled_once_for_a_whole_run() {
    // The finding of `docs/performance-notes.md`, stated as a test: the
    // control scenario's 4 samplings per step, every step, are one sampling
    // for the run.
    let basin = test_basin();
    let mut solver = test_solver(basin);
    let mut forcing = WindForcing::new(basin, CountingWind::new(trade_winds()));
    let mut state = OceanState::at_rest(basin.grid());

    for step in 0..COMPARED_STEPS {
        solver.step_with_forcing(&mut state, step as f64 * DT_S, &mut forcing);
    }

    assert_eq!(
        forcing.wind().calls(),
        evaluations_per_sampling(basin),
        "a steady wind was re-evaluated after the field was already in hand"
    );
}

#[test]
fn a_time_varying_wind_is_sampled_once_per_distinct_stage_time() {
    // The saving a wind that genuinely varies is entitled to, and no more:
    // RK4's four stages ask about three instants, and the last stage of a step
    // asks about the instant the next step's first stage asks about. What is
    // *not* allowed is holding a stale field across a change of `t`, which is
    // what the equivalence tests below check the values of.
    let basin = test_basin();
    let mut solver = test_solver(basin);
    let mut forcing = WindForcing::new(basin, CountingWind::new(burst()));
    let mut state = OceanState::at_rest(basin.grid());

    for step in 0..COMPARED_STEPS {
        solver.step_with_forcing(&mut state, step as f64 * DT_S, &mut forcing);
    }

    let expected = distinct_consecutive(&stage_times(COMPARED_STEPS));
    assert_eq!(
        forcing.wind().calls(),
        expected * evaluations_per_sampling(basin),
        "a time-varying wind was sampled a different number of times than the \
         RK4 tableau has distinct stage times"
    );
}

#[test]
fn the_solvers_wind_taking_step_samples_a_steady_wind_once_per_step() {
    // `Solver::step_forced_by` takes the wind itself rather than a forcing it
    // can keep, so its cache cannot outlive the call: a steady wind costs one
    // sampling per step instead of four, and a time-varying one three.
    let basin = test_basin();
    let mut solver = test_solver(basin);
    let steady = CountingWind::new(trade_winds());
    let varying = CountingWind::new(burst());
    let mut state = OceanState::at_rest(basin.grid());
    let mut other_state = OceanState::at_rest(basin.grid());

    for step in 0..COMPARED_STEPS {
        let t_s = step as f64 * DT_S;
        solver.step_forced_by(&mut state, t_s, basin, &steady);
        solver.step_forced_by(&mut other_state, t_s, basin, &varying);
    }

    let per_sampling = evaluations_per_sampling(basin);
    assert_eq!(steady.calls(), COMPARED_STEPS * per_sampling);
    // Three distinct stage times per step: `t`, `t + dt/2` and `t + dt`.
    assert_eq!(varying.calls(), 3 * COMPARED_STEPS * per_sampling);
}

// --- That the cache is correct. ---

#[test]
fn a_cached_field_is_the_field_a_fresh_sample_gives_at_that_instant() {
    // Asked for instants in an order no run would use — forwards, backwards,
    // repeated — the forcing must return the field of the instant asked for,
    // not the one it happens to hold.
    let basin = test_basin();
    let wind = seasonal_winds_with_a_burst();
    let mut forcing = WindForcing::new(basin, seasonal_winds_with_a_burst());

    for &t_s in &[
        0.0,
        0.0,
        BURST_PEAK_TIME_S,
        BURST_PEAK_TIME_S,
        DT_S,
        0.0,
        BURST_PEAK_TIME_S + BURST_DURATION_S,
        TROPICAL_YEAR_S / 4.0,
        DT_S,
    ] {
        assert_eq!(
            forcing.at(t_s),
            &WindStressField::sampled(basin, &wind, t_s),
            "the field held at t = {t_s} s is not the field of that instant"
        );
    }
}

#[test]
fn a_cached_burst_carries_its_closed_form_in_space_and_time() {
    // Tied to the physics rather than to another code path: `WindBurstAnomaly`
    // is `τ₀·exp(−((x−x₀)/Lx)²)·exp(−(y/Ly)²)·exp(−((t−t₀)/Lt)²)`
    // (`CONTEXT.md`, *Westerly wind burst*), evaluated here independently of
    // the engine's own arithmetic.
    let basin = test_basin();
    let mut forcing = WindForcing::new(basin, burst());
    let closed_form = |x_m: f64, y_m: f64, t_s: f64| {
        BURST_STRESS_PA
            * (-(((x_m - BURST_CENTER_X_M) / BURST_ZONAL_SCALE_M).powi(2))).exp()
            * (-((y_m / BURST_MERIDIONAL_SCALE_M).powi(2))).exp()
            * (-(((t_s - BURST_PEAK_TIME_S) / BURST_DURATION_S).powi(2))).exp()
    };

    for &t_s in &[0.0, BURST_PEAK_TIME_S, BURST_PEAK_TIME_S + BURST_DURATION_S] {
        let field = forcing.at(t_s);
        for (i, j) in [(0_usize, 0_usize), (3, 4), (BASIN_NX, BASIN_NY - 1)] {
            let x_m = basin.x_of_column_m(engine::U_STAGGERING, i);
            let y_m = basin.y_of_row_m(engine::U_STAGGERING, j);
            let expected_pa = closed_form(x_m, y_m, t_s);
            let sampled_pa = *field
                .tau_x_pa()
                .get(i, j)
                .expect("the probe is inside the field");

            assert!(
                (sampled_pa - expected_pa).abs()
                    <= ROUNDING_TOLERANCE * expected_pa.abs().max(BURST_STRESS_PA),
                "τx at ({x_m} m, {y_m} m, {t_s} s) is {sampled_pa} Pa, not the \
                 closed form's {expected_pa} Pa"
            );
        }
    }
}

#[test]
fn a_seasonal_run_reaches_the_state_the_uncached_stepper_reaches() {
    // The ticket's headline: a season genuinely varies, so a cache that never
    // invalidated would freeze the modulation at its `t = 0` value and the two
    // runs would part company within a step.
    assert_cached_run_matches_uncached(seasonal_winds());
}

#[test]
fn a_burst_run_reaches_the_state_the_uncached_stepper_reaches() {
    // `engine/tests/wind_burst.rs` is the physics of this forcing; this is the
    // arithmetic being unchanged by the cache, at the point of the burst's
    // steepest growth in `t`.
    assert_cached_run_matches_uncached(burst());
}

#[test]
fn a_seasonal_run_with_a_burst_reaches_the_state_the_uncached_stepper_reaches() {
    // The composite the ticket calls out: two time-varying components stacked,
    // so neither the sum nor either term may be frozen.
    assert_cached_run_matches_uncached(seasonal_winds_with_a_burst());
}

#[test]
fn a_steady_run_reaches_the_state_the_uncached_stepper_reaches() {
    // The case the cache is for. It would be a poor optimisation that changed
    // the control scenario's answer while making it faster, and the Epic 07
    // validations run on exactly this forcing.
    assert_cached_run_matches_uncached(trade_winds());
}

/// Run `wind` through both steppers for [`COMPARED_STEPS`] steps and require
/// the states to agree bit for bit at every one of them.
///
/// Bit-for-bit rather than within a tolerance: the cached and uncached paths
/// perform the same operations in the same order over the same `f64` fields,
/// so the only thing a difference can mean is that they were given different
/// stresses.
fn assert_cached_run_matches_uncached<W: WindStress>(wind: W) {
    let basin = test_basin();
    let mut solver = test_solver(basin);
    let mut cached_state = OceanState::at_rest(basin.grid());
    let mut forced_state = OceanState::at_rest(basin.grid());
    let mut reference_state = OceanState::at_rest(basin.grid());
    let mut reference = UncachedStepper::new(basin);
    let mut forcing = WindForcing::new(basin, wind);

    for step in 0..COMPARED_STEPS {
        let t_s = step as f64 * DT_S;
        solver.step_with_forcing(&mut cached_state, t_s, &mut forcing);
        solver.step_forced_by(&mut forced_state, t_s, basin, forcing.wind());
        reference.step(&mut reference_state, t_s, forcing.wind());

        assert_eq!(
            &cached_state,
            &reference_state,
            "the cached run diverged from the uncached one at step {}",
            step + 1
        );
        assert_eq!(
            &forced_state,
            &reference_state,
            "`step_forced_by` diverged from the uncached stepper at step {}",
            step + 1
        );
    }
}
