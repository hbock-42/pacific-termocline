//! Acceptance tests for T-03.3 — the idealized westerly wind-burst anomaly
//! and the `CompositeWind` combinator that stacks it on a base scenario.
//!
//! Three things are under test.
//!
//! The first is the burst field itself: a `WindBurstAnomaly` is a *westerly*
//! stress — `τx > 0`, the opposite sign to the alizés — that is Gaussian in
//! `x`, in `y` about the equator, and in `t` (`CONTEXT.md`, *Westerly wind
//! burst*). Every expected value below is the closed form
//! `τ₀·exp(−((x−x₀)/Lx)²)·exp(−(y/Ly)²)·exp(−((t−t₀)/Lt)²)` written out from
//! the definition and evaluated independently in the test.
//!
//! The second is composability, which is the half of the ticket that says the
//! burst is "addable on top of steady or seasonal winds, not exclusive with
//! them": a `CompositeWind` is the pointwise sum of its components, and
//! sampling one onto the C-grid gives the sum of the sampled components.
//!
//! The third is the ticket's acceptance criterion: *injecting a burst on top
//! of steady trade winds and running forward shows a visible eastward-
//! propagating thermocline-depth signal consistent with Kelvin wave
//! behaviour.*
//!
//! # How the burst's signal is isolated
//!
//! The v1 core is linear (CODING_STANDARDS.md § Scope guards), so the state of
//! a run forced by `trades + burst` is exactly the state of a run forced by
//! `trades` alone plus the state of a run forced by the burst alone. The
//! signal the criterion is about is therefore the *difference* between the
//! composite run and the trade-wind run — the burst's own response, with the
//! basin's much larger wind-driven adjustment removed.
//! [`the_burst_response_is_the_composite_run_minus_the_trade_wind_run`] is the
//! test that this differencing is legitimate rather than a convenient fiction.
//!
//! # What "consistent with Kelvin wave behaviour" is checked to mean
//!
//! A Kelvin wave travels **eastward only**, non-dispersively, at
//! `c = √(g'·H)`, and is trapped within a few equatorial deformation radii
//! `Le = √(c/β)` of the equator (`CONTEXT.md`). Each clause is a test:
//!
//! - the differenced signal is a *deepening* (`h > 0`) that arrives at an
//!   eastern station later than at a station west of it;
//! - the speed implied by the two arrival times is `c`;
//! - nothing arrives ahead of the `c` front;
//! - the signal decays away from the equator on the `Le` scale.
//!
//! The runs are undamped (`r = 0`, which
//! [`PhysicalParams`](engine::PhysicalParams) admits as the validation limit
//! of `01-scientific-model.md`): a Rayleigh decay `exp(−r·t)` would shift a
//! pulse's peak time and blur the very measurement these tests make.
//!
//! The issue calls this a qualitative smoke test and defers the rigorous wave
//! check to Epic 07; the phase-speed measurement here is nevertheless made at
//! two resolutions, so that the error is shown to shrink under refinement
//! rather than merely to sit under one threshold
//! (CODING_STANDARDS.md § Tests).

use engine::{
    max_stable_dt, Basin, BetaPlane, CompositeWind, Grid, OceanState, PhysicalParams, Solver,
    Spacing, SteadyTradeWinds, WaveSpeed, WindBurstAnomaly, WindStress, WindStressError,
    WindStressField, H_STAGGERING, U_STAGGERING,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981).
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;

/// One solar day, in seconds.
const DAY_S: f64 = 86_400.0;

/// Zonal wind stress `τ₀` of the background alizés, in Pa — easterly, and the
/// observed scale of the equatorial Pacific's mean zonal stress.
const TRADE_WIND_STRESS_PA: f64 = -0.05;

/// Peak zonal stress of the burst, in Pa. Westerly, so positive, and 0.04 Pa
/// is the observed scale of an equatorial westerly wind burst — comparable to
/// the mean trades it briefly overwhelms.
const BURST_STRESS_PA: f64 = 0.04;
/// Zonal centre of the burst, in metres east of the western boundary: in the
/// western basin, leaving 8×10⁶ m of open water to its east for the Kelvin
/// signal to cross.
const BURST_CENTER_X_M: f64 = 2.0e6;
/// Zonal `e`-folding scale of the burst, in metres — about 5° of longitude,
/// the scale of an observed burst, and five cells of the coarse test basin.
const BURST_ZONAL_SCALE_M: f64 = 5.0e5;
/// Time of the burst's peak, in seconds — 15 days in, three duration scales
/// after the run starts, so the run begins with the burst switched off to
/// within `exp(−9) ≈ 10⁻⁴` of its peak.
const BURST_PEAK_TIME_S: f64 = 15.0 * DAY_S;
/// Temporal `e`-folding scale of the burst, in seconds — 5 days, the duration
/// of an observed westerly wind burst.
const BURST_DURATION_S: f64 = 5.0 * DAY_S;

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*).
const BASIN_LX_M: f64 = 1.0e7;
/// Meridional extent of the test basin, in metres: ±10⁶ m, which is 2.9
/// equatorial deformation radii either side of the equator, so the Kelvin
/// structure `exp(−y²/(2·Le²))` is down to 1.5% of its peak at the walls.
const BASIN_LY_M: f64 = 2.0e6;

/// Length of the propagation runs, in seconds — 50 days.
///
/// Long enough for the burst's Kelvin signal to reach the eastern measurement
/// station ([`EAST_STATION_X_M`], at 40 days), and short enough that the
/// pulse's reflection off the eastern wall — which would arrive back there at
/// 57 days — has not yet contaminated it.
const RUN_S: f64 = 50.0 * DAY_S;

/// Western measurement station, in metres east of the western boundary.
const WEST_STATION_X_M: f64 = 4.0e6;
/// Eastern measurement station, in metres east of the western boundary. Four
/// million metres — 24 hours' short of 17 days of Kelvin travel — east of
/// [`WEST_STATION_X_M`], which is the baseline the phase speed is measured
/// over.
const EAST_STATION_X_M: f64 = 8.0e6;

/// Meridional positions the burst profile is probed at, in metres.
const PROBE_LATITUDES_M: [f64; 4] = [0.0, 1.0e5, -5.0e5, 2.0e6];

/// Relative slack allowed where a check is exact in exact arithmetic: a few
/// tens of ulps of `f64` (ε ≈ 2.2×10⁻¹⁶) for the handful of operations per
/// point the expression costs.
const ROUNDING_TOLERANCE: f64 = 1.0e-14;

/// The equatorial-Pacific parameter set, undamped.
///
/// `r = 0` is the validation limit of `docs/planning/01-scientific-model.md`
/// and the one in which a Kelvin pulse keeps its shape, so that a peak time is
/// a propagation time and not a propagation time biased by decay.
fn undamped_pacific_params() -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        0.0,
        engine::EQUATORIAL_BETA_PER_M_PER_S,
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

/// The equatorial deformation radius `Le = √(c/β)`, in metres (`CONTEXT.md`).
/// About 3.45×10⁵ m for the parameters above.
fn equatorial_deformation_radius_m(params: PhysicalParams) -> f64 {
    (params.kelvin_wave_speed_m_per_s() / params.beta_per_m_per_s()).sqrt()
}

/// A basin [`BASIN_LX_M`] by [`BASIN_LY_M`] metres, resolved by `nx` by `ny`
/// cells and centred on the equator.
///
/// `ny` must be even, so that the equator falls exactly between two rows of
/// `h` and [`equatorial_h_m`] can average the pair straddling it.
fn equatorial_basin(nx: usize, ny: usize) -> Basin {
    assert!(
        ny.is_multiple_of(2),
        "the equator must fall between two rows of h"
    );
    let grid = Grid::new(nx, ny).expect("extents are non-zero");
    let spacing = Spacing::new(BASIN_LX_M / nx as f64, BASIN_LY_M / ny as f64)
        .expect("a basin spanned by whole cells has positive spacing");
    Basin::centered_on_equator(grid, spacing)
}

/// The burst of the propagation tests: the constants above, as a scenario.
fn pacific_burst() -> WindBurstAnomaly {
    WindBurstAnomaly::new(
        BURST_STRESS_PA,
        BURST_CENTER_X_M,
        BURST_ZONAL_SCALE_M,
        equatorial_deformation_radius_m(undamped_pacific_params()),
        BURST_PEAK_TIME_S,
        BURST_DURATION_S,
    )
    .expect("a westerly burst with positive scales and a positive duration")
}

/// The steady alizés the burst is superimposed on.
fn pacific_trade_winds() -> SteadyTradeWinds {
    SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA).expect("an easterly stress is a trade wind")
}

// --- The burst field itself. ---

#[test]
fn a_westerly_burst_is_strongest_at_its_own_centre_in_space_and_time() {
    // "Westerly" is the whole name of the scenario: `τx > 0`, the opposite
    // sign to the alizés of `CONTEXT.md`. At the centre of the three Gaussians
    // every factor is exactly one, so the stress is the configured peak.
    let burst = pacific_burst();

    let (tau_x_pa, tau_y_pa) = burst.stress(BURST_CENTER_X_M, 0.0, BURST_PEAK_TIME_S);

    assert!(
        tau_x_pa > 0.0,
        "a burst must be westerly, got {tau_x_pa} Pa"
    );
    assert_eq!(tau_x_pa, BURST_STRESS_PA);
    assert_eq!(
        tau_y_pa, 0.0,
        "an idealized burst carries no meridional stress"
    );
}

#[test]
fn the_burst_decays_as_a_gaussian_in_x_in_y_and_in_t() {
    // Checked against `exp(−s²)` at `s` scales off centre, evaluated here
    // rather than read back from the module: one scale away the stress is
    // `τ₀/e`, two away `τ₀/e⁴`, and each direction is symmetric about its
    // centre because the offset enters squared.
    let burst = pacific_burst();
    let meridional_scale_m = burst.meridional_scale_m();

    for (scales_off_centre, expected_fraction) in [
        (0.0, 1.0),
        (1.0, 1.0 / std::f64::consts::E),
        (2.0, (-4.0_f64).exp()),
        (3.0, (-9.0_f64).exp()),
    ] {
        let expected_pa = BURST_STRESS_PA * expected_fraction;
        for sign in [-1.0, 1.0] {
            let offset = sign * scales_off_centre;
            let probes = [
                (
                    BURST_CENTER_X_M + offset * BURST_ZONAL_SCALE_M,
                    0.0,
                    BURST_PEAK_TIME_S,
                ),
                (
                    BURST_CENTER_X_M,
                    offset * meridional_scale_m,
                    BURST_PEAK_TIME_S,
                ),
                (
                    BURST_CENTER_X_M,
                    0.0,
                    BURST_PEAK_TIME_S + offset * BURST_DURATION_S,
                ),
            ];
            for (x_m, y_m, t_s) in probes {
                let (tau_x_pa, _) = burst.stress(x_m, y_m, t_s);
                assert!(
                    (tau_x_pa - expected_pa).abs() <= ROUNDING_TOLERANCE * expected_pa.abs(),
                    "at (x = {x_m} m, y = {y_m} m, t = {t_s} s) the stress is {tau_x_pa} Pa, \
                     expected {expected_pa} Pa"
                );
            }
        }
    }
}

#[test]
fn the_burst_is_the_product_of_its_three_gaussians() {
    // Displaced in all three directions at once, the stress is the product of
    // the three factors — the separable form the ticket asks for, rather than
    // a Gaussian in one combined distance.
    let burst = pacific_burst();
    let meridional_scale_m = burst.meridional_scale_m();
    let (x_offset, y_offset, t_offset) = (1.5, 0.75, 2.25);

    let (tau_x_pa, _) = burst.stress(
        BURST_CENTER_X_M + x_offset * BURST_ZONAL_SCALE_M,
        y_offset * meridional_scale_m,
        BURST_PEAK_TIME_S + t_offset * BURST_DURATION_S,
    );

    let expected_pa = BURST_STRESS_PA
        * (-(x_offset * x_offset + y_offset * y_offset + t_offset * t_offset)).exp();
    assert!(
        (tau_x_pa - expected_pa).abs() <= ROUNDING_TOLERANCE * expected_pa.abs(),
        "the stress is {tau_x_pa} Pa, expected the product {expected_pa} Pa"
    );
}

#[test]
fn a_burst_that_is_not_westerly_is_refused_by_name() {
    // Invalid scenario input is a `Result` naming the offending value, not a
    // panic and not a silently flipped sign (CODING_STANDARDS.md
    // § Correctness and failure).
    for value_pa in [0.0, TRADE_WIND_STRESS_PA, f64::NAN, f64::INFINITY] {
        let error = WindBurstAnomaly::new(
            value_pa,
            BURST_CENTER_X_M,
            BURST_ZONAL_SCALE_M,
            BURST_ZONAL_SCALE_M,
            BURST_PEAK_TIME_S,
            BURST_DURATION_S,
        )
        .expect_err("only a strictly positive stress is a westerly burst");
        let WindStressError::NotWesterly {
            value_pa: rejected_pa,
        } = error
        else {
            panic!("expected the stress itself to be rejected, got {error}");
        };
        // Compared bitwise rather than with `==`, so that the NaN case checks
        // the value was carried through rather than trivially passing.
        assert_eq!(rejected_pa.to_bits(), value_pa.to_bits());
        let message = error.to_string();
        assert!(message.contains("strictly positive"), "{message}");
    }
}

#[test]
fn a_burst_scale_that_is_not_a_distance_is_refused_by_name() {
    for parameter in ["zonal_scale_m", "meridional_scale_m"] {
        for value_m in [0.0, -BURST_ZONAL_SCALE_M, f64::NAN, f64::INFINITY] {
            let (zonal_scale_m, meridional_scale_m) = if parameter == "zonal_scale_m" {
                (value_m, BURST_ZONAL_SCALE_M)
            } else {
                (BURST_ZONAL_SCALE_M, value_m)
            };
            let error = WindBurstAnomaly::new(
                BURST_STRESS_PA,
                BURST_CENTER_X_M,
                zonal_scale_m,
                meridional_scale_m,
                BURST_PEAK_TIME_S,
                BURST_DURATION_S,
            )
            .expect_err("a burst scale must be a finite, positive distance");
            let WindStressError::ScaleNotPositive {
                parameter: rejected_parameter,
                value_m: rejected_m,
            } = error
            else {
                panic!("expected the scale to be rejected, got {error}");
            };
            assert_eq!(rejected_parameter, parameter);
            assert_eq!(rejected_m.to_bits(), value_m.to_bits());
            assert!(error.to_string().contains(parameter), "{error}");
        }
    }
}

#[test]
fn a_burst_duration_that_is_not_a_duration_is_refused_by_name() {
    for value_s in [0.0, -BURST_DURATION_S, f64::NAN, f64::INFINITY] {
        let error = WindBurstAnomaly::new(
            BURST_STRESS_PA,
            BURST_CENTER_X_M,
            BURST_ZONAL_SCALE_M,
            BURST_ZONAL_SCALE_M,
            BURST_PEAK_TIME_S,
            value_s,
        )
        .expect_err("a burst duration must be a finite, positive time");
        let WindStressError::DurationNotPositive {
            value_s: rejected_s,
        } = error
        else {
            panic!("expected the duration to be rejected, got {error}");
        };
        assert_eq!(rejected_s.to_bits(), value_s.to_bits());
        assert!(error.to_string().contains("duration_s"), "{error}");
    }
}

#[test]
fn a_burst_centred_nowhere_is_refused_by_name() {
    for value_m in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let error = WindBurstAnomaly::new(
            BURST_STRESS_PA,
            value_m,
            BURST_ZONAL_SCALE_M,
            BURST_ZONAL_SCALE_M,
            BURST_PEAK_TIME_S,
            BURST_DURATION_S,
        )
        .expect_err("a burst centre must be a finite position");
        let WindStressError::CenterNotAPosition {
            value_m: rejected_m,
        } = error
        else {
            panic!("expected the centre to be rejected, got {error}");
        };
        assert_eq!(rejected_m.to_bits(), value_m.to_bits());
        assert!(error.to_string().contains("center_x_m"), "{error}");
    }
}

#[test]
fn a_burst_that_peaks_at_no_time_is_refused_by_name() {
    for value_s in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let error = WindBurstAnomaly::new(
            BURST_STRESS_PA,
            BURST_CENTER_X_M,
            BURST_ZONAL_SCALE_M,
            BURST_ZONAL_SCALE_M,
            value_s,
            BURST_DURATION_S,
        )
        .expect_err("a burst must peak at a finite time");
        let WindStressError::PeakTimeNotFinite {
            value_s: rejected_s,
        } = error
        else {
            panic!("expected the peak time to be rejected, got {error}");
        };
        assert_eq!(rejected_s.to_bits(), value_s.to_bits());
        assert!(error.to_string().contains("peak_time_s"), "{error}");
    }
}

// --- Composability: stacking scenarios. ---

#[test]
fn a_composite_wind_is_the_sum_of_its_components() {
    // The point of the combinator: a burst is *added to* the alizés rather
    // than replacing them, so at every point the composite is the arithmetic
    // sum, computed here from the two components separately.
    let trades = pacific_trade_winds();
    let burst = pacific_burst();
    let composite = CompositeWind::new().with(trades).with(burst);

    for x_m in [0.0, BURST_CENTER_X_M, BASIN_LX_M] {
        for y_m in PROBE_LATITUDES_M {
            for t_s in [0.0, BURST_PEAK_TIME_S, RUN_S] {
                let (trade_x_pa, trade_y_pa) = trades.stress(x_m, y_m, t_s);
                let (burst_x_pa, burst_y_pa) = burst.stress(x_m, y_m, t_s);
                let (tau_x_pa, tau_y_pa) = composite.stress(x_m, y_m, t_s);
                assert_eq!(
                    tau_x_pa,
                    trade_x_pa + burst_x_pa,
                    "τx at ({x_m}, {y_m}, {t_s})"
                );
                assert_eq!(
                    tau_y_pa,
                    trade_y_pa + burst_y_pa,
                    "τy at ({x_m}, {y_m}, {t_s})"
                );
            }
        }
    }
}

#[test]
fn a_composite_wind_of_nothing_is_calm() {
    // The identity of the combinator, and the case a scenario with no forcing
    // at all lands on.
    let composite = CompositeWind::new();

    assert!(composite.is_empty());
    for y_m in PROBE_LATITUDES_M {
        assert_eq!(
            composite.stress(BURST_CENTER_X_M, y_m, BURST_PEAK_TIME_S),
            (0.0, 0.0)
        );
    }
}

#[test]
fn a_composite_wind_of_one_component_is_that_component() {
    let trades = pacific_trade_winds();
    let composite = CompositeWind::new().with(trades);

    assert_eq!(composite.len(), 1);
    for y_m in PROBE_LATITUDES_M {
        assert_eq!(
            composite.stress(BURST_CENTER_X_M, y_m, BURST_PEAK_TIME_S),
            trades.stress(BURST_CENTER_X_M, y_m, BURST_PEAK_TIME_S)
        );
    }
}

#[test]
fn a_composite_wind_stacks_a_burst_on_a_burst() {
    // "Composable" has to mean composable with anything implementing the
    // trait, not just with the trades: two bursts at different places stack
    // into one field with two maxima.
    let first = pacific_burst();
    let second = WindBurstAnomaly::new(
        BURST_STRESS_PA,
        BURST_CENTER_X_M + 4.0 * BURST_ZONAL_SCALE_M,
        BURST_ZONAL_SCALE_M,
        first.meridional_scale_m(),
        BURST_PEAK_TIME_S,
        BURST_DURATION_S,
    )
    .expect("a westerly burst with positive scales and a positive duration");
    let composite = CompositeWind::new().with(first).with(second);

    assert_eq!(composite.len(), 2);
    for centre_x_m in [
        BURST_CENTER_X_M,
        BURST_CENTER_X_M + 4.0 * BURST_ZONAL_SCALE_M,
    ] {
        let (tau_x_pa, _) = composite.stress(centre_x_m, 0.0, BURST_PEAK_TIME_S);
        // Each centre carries its own peak plus the far tail of the other,
        // `exp(−16) ≈ 10⁻⁷` of a peak: bounded well inside 1% of the peak.
        let excess = (tau_x_pa - BURST_STRESS_PA).abs();
        assert!(
            excess > 0.0 && excess <= 0.01 * BURST_STRESS_PA,
            "at x = {centre_x_m} m the composite is {tau_x_pa} Pa, expected one peak plus the \
             other's tail"
        );
    }
}

#[test]
fn sampling_a_composite_is_sampling_the_components_and_adding_them() {
    // The trait composes; the discretisation must too. Interior faces carry
    // the sum, and the basin's walls stay at exactly zero — the sampling rule
    // of the `forcing` module header applies to a composite like any other
    // scenario.
    let basin = equatorial_basin(20, 4);
    let trades = pacific_trade_winds();
    let burst = pacific_burst();
    let composite = CompositeWind::new().with(trades).with(burst);
    let t_s = BURST_PEAK_TIME_S;

    let field = WindStressField::sampled(basin, &composite, t_s);

    let trade_field = WindStressField::sampled(basin, &trades, t_s);
    let burst_field = WindStressField::sampled(basin, &burst, t_s);
    let nx = basin.grid().nx();
    for j in 0..field.tau_x_pa().ny() {
        for i in 0..=nx {
            let sampled_pa = *field.tau_x_pa().get(i, j).expect("an in-bounds face");
            let expected_pa = *trade_field.tau_x_pa().get(i, j).expect("an in-bounds face")
                + *burst_field.tau_x_pa().get(i, j).expect("an in-bounds face");
            assert!(
                (sampled_pa - expected_pa).abs() <= ROUNDING_TOLERANCE * expected_pa.abs().max(1.0),
                "τx at face ({i}, {j}) is {sampled_pa} Pa, expected {expected_pa} Pa"
            );
            if i == 0 || i == nx {
                assert_eq!(
                    sampled_pa, 0.0,
                    "the wall face ({i}, {j}) must carry no stress"
                );
            }
        }
    }
    for j in 0..field.tau_y_pa().ny() {
        for i in 0..field.tau_y_pa().nx() {
            assert_eq!(
                *field.tau_y_pa().get(i, j).expect("an in-bounds face"),
                0.0,
                "neither the alizés nor an idealized burst carries a meridional stress"
            );
        }
    }
}

// --- The acceptance criterion: an eastward-propagating thermocline signal. ---

/// The equatorial thermocline anomaly `h` of the column `i`, in metres: the
/// mean of the two rows straddling the equator.
///
/// The basin has an even meridional cell count, so no row of `h` sits *on* the
/// equator; averaging the pair either side of it is the symmetric reading of
/// the equatorial value, and it is where a Kelvin wave's amplitude peaks.
fn equatorial_h_m(state: &OceanState, i: usize) -> f64 {
    let ny = state.grid().ny();
    let south = *state.h().get(i, ny / 2 - 1).expect("an in-bounds cell");
    let north = *state.h().get(i, ny / 2).expect("an in-bounds cell");
    (south + north) / 2.0
}

/// One run's equatorial thermocline anomaly, sampled at every step.
struct EquatorialHistory {
    /// Length of one step, in seconds.
    dt_s: f64,
    /// `h` in metres, indexed by step then by column: `h_m[step][i]`.
    h_m: Vec<Vec<f64>>,
}

impl EquatorialHistory {
    /// Time of the step `index`, in seconds since the start of the run.
    fn time_of_step_s(&self, index: usize) -> f64 {
        index as f64 * self.dt_s
    }

    /// This history minus `other`, step by step and column by column.
    ///
    /// The two runs must share a basin and a timestep, which they do because
    /// only the forcing differs between them.
    fn difference(&self, other: &Self) -> Self {
        assert_eq!(self.dt_s, other.dt_s);
        assert_eq!(self.h_m.len(), other.h_m.len());
        Self {
            dt_s: self.dt_s,
            h_m: self
                .h_m
                .iter()
                .zip(&other.h_m)
                .map(|(mine, theirs)| {
                    mine.iter()
                        .zip(theirs)
                        .map(|(mine, theirs)| mine - theirs)
                        .collect()
                })
                .collect(),
        }
    }

    /// The largest `h` reached in the column `i` over the whole run, in metres.
    fn peak_h_m(&self, i: usize) -> f64 {
        self.h_m
            .iter()
            .map(|row| row[i])
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// The time, in seconds, at which `h` in the column `i` is largest.
    ///
    /// Refined below the sampling interval by fitting a parabola through the
    /// largest sample and its two neighbours — the standard sub-sample peak
    /// estimate, which is what makes the arrival time second-order accurate in
    /// `dt` rather than quantised by it.
    fn peak_time_s(&self, i: usize) -> f64 {
        let steps = self.h_m.len();
        let argmax = (0..steps)
            .max_by(|&a, &b| {
                self.h_m[a][i]
                    .partial_cmp(&self.h_m[b][i])
                    .expect("a finite run carries no NaN")
            })
            .expect("a run has at least one step");
        assert!(
            argmax > 0 && argmax + 1 < steps,
            "the peak in column {i} is at step {argmax} of {steps}, on the edge of the run: the \
             run is too short, or too long and contaminated by a wall reflection"
        );
        let (before, at, after) = (
            self.h_m[argmax - 1][i],
            self.h_m[argmax][i],
            self.h_m[argmax + 1][i],
        );
        let curvature = before - 2.0 * at + after;
        let offset_in_steps = if curvature == 0.0 {
            0.0
        } else {
            0.5 * (before - after) / curvature
        };
        self.time_of_step_s(argmax) + offset_in_steps * self.dt_s
    }
}

/// A solver for `basin` at the longest timestep its CFL bound admits, and that
/// timestep.
///
/// The rig every run in this file starts from: the two ways a run is made —
/// [`run_recording_equatorial_h`], which keeps the equatorial history, and
/// [`state_after`], which keeps only the final state — differ in what they
/// record, not in how they are set up.
fn solver_for(basin: Basin, params: PhysicalParams) -> (Solver, f64) {
    let wave_speed =
        WaveSpeed::new(params.kelvin_wave_speed_m_per_s()).expect("a positive wave speed");
    let dt_s = max_stable_dt(basin.spacing(), wave_speed);
    let plane = BetaPlane::centered_on_equator(params, basin.spacing(), basin.grid());
    let solver = Solver::new(basin.grid(), basin.spacing(), params, plane, dt_s)
        .unwrap_or_else(|error| panic!("the test's own timestep must be admissible: {error}"));
    (solver, dt_s)
}

/// Run `basin` from rest under `wind` for [`RUN_S`] seconds, recording the
/// equatorial `h` of every column at every step.
fn run_recording_equatorial_h(
    basin: Basin,
    params: PhysicalParams,
    wind: &dyn WindStress,
) -> EquatorialHistory {
    let (mut solver, dt_s) = solver_for(basin, params);

    let mut state = OceanState::at_rest(basin.grid());
    let steps = (RUN_S / dt_s).ceil() as usize;
    let nx = basin.grid().nx();
    let mut h_m = Vec::with_capacity(steps + 1);
    for step in 0..=steps {
        h_m.push((0..nx).map(|i| equatorial_h_m(&state, i)).collect());
        if step < steps {
            solver.step_forced_by(&mut state, step as f64 * dt_s, basin, wind);
        }
    }
    EquatorialHistory { dt_s, h_m }
}

/// The burst's own signal on `basin`: the composite run minus the run forced
/// by the trade winds alone.
fn burst_signal(basin: Basin, params: PhysicalParams) -> EquatorialHistory {
    let trades = pacific_trade_winds();
    let composite = CompositeWind::new().with(trades).with(pacific_burst());
    let with_burst = run_recording_equatorial_h(basin, params, &composite);
    let without_burst = run_recording_equatorial_h(basin, params, &trades);
    with_burst.difference(&without_burst)
}

/// The index in `0..count` whose `position_m` is nearest `target_m`.
fn index_nearest_m(count: usize, target_m: f64, position_m: impl Fn(usize) -> f64) -> usize {
    (0..count)
        .min_by(|&a, &b| {
            let distance = |index: usize| (position_m(index) - target_m).abs();
            distance(a)
                .partial_cmp(&distance(b))
                .expect("positions are finite")
        })
        .expect("a basin has at least one cell in each direction")
}

/// The column of `basin` whose centre is nearest `x_m` metres east of the
/// western boundary.
fn column_nearest_x(basin: Basin, x_m: f64) -> usize {
    index_nearest_m(basin.grid().nx(), x_m, |i| {
        basin.x_of_column_m(H_STAGGERING, i)
    })
}

#[test]
fn a_westerly_burst_deepens_the_thermocline_and_the_signal_travels_east() {
    // The acceptance criterion, in its two halves. A westerly burst drives an
    // eastward acceleration that converges the layer to its east, so the
    // signal is a *deepening* — `h > 0`, the El Niño-onset sense of
    // `CONTEXT.md`, *Westerly wind burst* — and it must reach the eastern
    // station after the western one, never before.
    //
    // "Visible in raw field data" is a statement about magnitude. The scale
    // analysis `h ~ (τ/(ρ₀·H))·Lt·(H/Lx)·Lt ≈ 10 m` bounds the whole response
    // of the burst above, of which the Kelvin mode carries a fraction — it is
    // an upper bound, not a prediction — so the threshold is set two orders of
    // magnitude below it. A tenth of a metre still excludes a numerically
    // negligible wobble.
    const VISIBLE_SIGNAL_M: f64 = 0.1;
    let basin = equatorial_basin(100, 20);
    let signal = burst_signal(basin, undamped_pacific_params());

    let west = column_nearest_x(basin, WEST_STATION_X_M);
    let east = column_nearest_x(basin, EAST_STATION_X_M);
    let west_peak_m = signal.peak_h_m(west);
    let east_peak_m = signal.peak_h_m(east);

    assert!(
        west_peak_m > VISIBLE_SIGNAL_M && east_peak_m > VISIBLE_SIGNAL_M,
        "a westerly burst must leave a visible deepening at both stations, got {west_peak_m} m \
         in the west and {east_peak_m} m in the east"
    );
    assert!(
        signal.peak_time_s(east) > signal.peak_time_s(west),
        "the signal peaks at {} s in the east and {} s in the west: it is not travelling eastward",
        signal.peak_time_s(east),
        signal.peak_time_s(west)
    );
}

/// The zonal phase speed of the burst's signal on `basin`, in m/s: the
/// distance between the two measurement stations divided by the difference of
/// the times at which each sees its peak.
fn measured_phase_speed_m_per_s(basin: Basin, params: PhysicalParams) -> f64 {
    let signal = burst_signal(basin, params);
    let west = column_nearest_x(basin, WEST_STATION_X_M);
    let east = column_nearest_x(basin, EAST_STATION_X_M);
    let distance_m =
        basin.x_of_column_m(H_STAGGERING, east) - basin.x_of_column_m(H_STAGGERING, west);
    distance_m / (signal.peak_time_s(east) - signal.peak_time_s(west))
}

#[test]
fn the_signal_travels_at_the_kelvin_wave_speed() {
    // The independent expected value is `c = √(g'·H)` from the published
    // parameters — 2.74 m/s — not anything read out of the model.
    //
    // The tolerance is the centred difference's numerical dispersion. Its
    // phase speed is `c·sin(k·Δx)/(k·Δx)`, low by `(k·Δx)²/6`; the pulse
    // carries its power around `k ≈ 1/Lx` and `Lx` is five cells of this
    // basin, so `k·Δx ≈ 0.2` and the error is `(0.2)²/6 ≈ 0.7%`. One percent
    // is that rounded up. The convergence test below is what shows the
    // discrepancy really is this discretisation error and not a coincidence.
    const DISPERSION_TOLERANCE: f64 = 0.01;
    let params = undamped_pacific_params();
    let basin = equatorial_basin(100, 20);

    let measured_m_per_s = measured_phase_speed_m_per_s(basin, params);

    let kelvin_m_per_s = params.kelvin_wave_speed_m_per_s();
    let relative_error = (measured_m_per_s - kelvin_m_per_s).abs() / kelvin_m_per_s;
    assert!(
        relative_error <= DISPERSION_TOLERANCE,
        "the signal travels at {measured_m_per_s} m/s against a Kelvin speed of \
         {kelvin_m_per_s} m/s, a relative error of {relative_error}"
    );
}

#[test]
fn the_measured_speed_approaches_the_kelvin_speed_under_refinement() {
    // Convergence rather than a point check (CODING_STANDARDS.md § Tests). The
    // scheme is second-order in space, so halving both cell dimensions — which
    // halves the timestep with them — must cut the phase-speed error by four.
    // The bound is 0.35 rather than 0.25 because the measurement is not exact:
    // each arrival time is a sub-sample parabolic estimate, itself only
    // second-order accurate in `dt`, so the measured ratio carries a little of
    // its own error. A modelling error — one that survived refinement — would
    // leave the ratio at one.
    const REQUIRED_ERROR_REDUCTION: f64 = 0.35;
    let params = undamped_pacific_params();
    let kelvin_m_per_s = params.kelvin_wave_speed_m_per_s();

    let error_at = |nx: usize, ny: usize| {
        let measured = measured_phase_speed_m_per_s(equatorial_basin(nx, ny), params);
        (measured - kelvin_m_per_s).abs() / kelvin_m_per_s
    };
    let coarse = error_at(100, 20);
    let fine = error_at(200, 40);

    assert!(
        fine <= REQUIRED_ERROR_REDUCTION * coarse,
        "halving the cell size took the phase-speed error from {coarse} to {fine}, which is not \
         the reduction a second-order discretisation must show"
    );
}

#[test]
fn nothing_arrives_ahead_of_the_kelvin_front() {
    // The Kelvin speed is the fastest signal in the model (`CONTEXT.md`), so
    // the eastern station must still be quiet three duration scales before the
    // burst's own crest could have reached it. The burst's Gaussian tail is
    // `exp(−9) ≈ 10⁻⁴` of its peak by then, and the second-order scheme's
    // dispersive ripples trail the front rather than lead it, so 1% of the
    // station's peak is that bound with two orders of magnitude of headroom.
    const AHEAD_OF_FRONT_TOLERANCE: f64 = 0.01;
    let params = undamped_pacific_params();
    let basin = equatorial_basin(100, 20);
    let signal = burst_signal(basin, params);

    let east = column_nearest_x(basin, EAST_STATION_X_M);
    let travel_s = (basin.x_of_column_m(H_STAGGERING, east) - BURST_CENTER_X_M)
        / params.kelvin_wave_speed_m_per_s();
    let quiet_until_s = BURST_PEAK_TIME_S + travel_s - 3.0 * BURST_DURATION_S;
    let peak_m = signal.peak_h_m(east);

    for (step, row) in signal.h_m.iter().enumerate() {
        let t_s = signal.time_of_step_s(step);
        if t_s > quiet_until_s {
            break;
        }
        assert!(
            row[east].abs() <= AHEAD_OF_FRONT_TOLERANCE * peak_m,
            "at t = {t_s} s the eastern station already reads {} m against a peak of {peak_m} m, \
             {} s before the Kelvin front can reach it",
            row[east],
            quiet_until_s - t_s
        );
    }
}

#[test]
fn the_signal_is_trapped_near_the_equator() {
    // The last clause of "consistent with Kelvin wave behaviour": the wave's
    // meridional structure is `exp(−y²/(2·Le²))`, so two deformation radii off
    // the equator its amplitude is down to `exp(−2) ≈ 14%` of its equatorial
    // value. Thirty percent is that bound roughly doubled, for the meridional
    // resolution of the test basin — `Le` is only 3.5 cells — and for the
    // westward Rossby response, which is not trapped as tightly.
    const OFF_EQUATOR_FRACTION: f64 = 0.3;
    let params = undamped_pacific_params();
    let basin = equatorial_basin(100, 20);
    let trades = pacific_trade_winds();
    let composite = CompositeWind::new().with(trades).with(pacific_burst());

    // Re-run rather than reuse `burst_signal`, which keeps only the equatorial
    // rows: this test needs a whole meridional column. The runs stop when the
    // crest is *at* the eastern station — `t₀ + (x − x₀)/c`, the arrival time
    // the speed test measures independently — because that is when the wave's
    // meridional structure is there to be read; at the end of the full run the
    // pulse has long passed the station and only its wake remains.
    let east = column_nearest_x(basin, EAST_STATION_X_M);
    let arrival_s = BURST_PEAK_TIME_S
        + (basin.x_of_column_m(H_STAGGERING, east) - BURST_CENTER_X_M)
            / params.kelvin_wave_speed_m_per_s();
    let with_burst = state_after(basin, params, &composite, arrival_s);
    let without_burst = state_after(basin, params, &trades, arrival_s);

    let anomaly_m = |j: usize| {
        *with_burst.h().get(east, j).expect("an in-bounds cell")
            - *without_burst.h().get(east, j).expect("an in-bounds cell")
    };
    let ny = basin.grid().ny();
    let equatorial_m = (anomaly_m(ny / 2 - 1) + anomaly_m(ny / 2)) / 2.0;
    let two_radii_m = 2.0 * equatorial_deformation_radius_m(params);
    let far = index_nearest_m(ny, two_radii_m, |j| basin.y_of_row_m(H_STAGGERING, j));

    assert!(
        equatorial_m.abs() > 0.0,
        "there is no equatorial signal to be trapped"
    );
    assert!(
        anomaly_m(far).abs() <= OFF_EQUATOR_FRACTION * equatorial_m.abs(),
        "the anomaly is {} m at y = {} m against {equatorial_m} m on the equator: it is not \
         equatorially trapped",
        anomaly_m(far),
        basin.y_of_row_m(H_STAGGERING, far)
    );
}

/// The state of `basin` after `run_s` seconds under `wind`, from rest.
fn state_after(
    basin: Basin,
    params: PhysicalParams,
    wind: &dyn WindStress,
    run_s: f64,
) -> OceanState {
    let (mut solver, dt_s) = solver_for(basin, params);
    let mut state = OceanState::at_rest(basin.grid());
    let steps = (run_s / dt_s).ceil() as usize;
    for step in 0..steps {
        solver.step_forced_by(&mut state, step as f64 * dt_s, basin, wind);
    }
    state
}

#[test]
fn the_burst_response_is_the_composite_run_minus_the_trade_wind_run() {
    // The v1 core is linear, so a run forced by `trades + burst` must be the
    // sum of a run forced by each alone. This is what licenses every test
    // above to read the burst's signal off a difference of two runs — and it
    // is a check on the combinator itself, since a `CompositeWind` that did
    // anything other than add would break it.
    //
    // The bound is `f64` rounding accumulated over the run: ε ≈ 2.2×10⁻¹⁶ over
    // ~150 steps of four RK4 stages, relative to the largest `h` in the basin.
    // 10⁻¹⁰ is that with four orders of magnitude of headroom for growth of
    // the rounding through the integration.
    const SUPERPOSITION_TOLERANCE: f64 = 1.0e-10;
    let params = undamped_pacific_params();
    let basin = equatorial_basin(100, 20);
    let trades = pacific_trade_winds();
    let burst = pacific_burst();
    let composite = CompositeWind::new().with(trades).with(burst);

    let with_both = state_after(basin, params, &composite, RUN_S);
    let trades_only = state_after(basin, params, &trades, RUN_S);
    let burst_only = state_after(basin, params, &burst, RUN_S);

    let scale_m = with_both
        .h()
        .as_slice()
        .iter()
        .fold(0.0_f64, |largest, &h_m| largest.max(h_m.abs()));
    assert!(scale_m > 0.0, "the composite run left the basin at rest");
    for j in 0..basin.grid().ny() {
        for i in 0..basin.grid().nx() {
            let both_m = *with_both.h().get(i, j).expect("an in-bounds cell");
            let sum_m = *trades_only.h().get(i, j).expect("an in-bounds cell")
                + *burst_only.h().get(i, j).expect("an in-bounds cell");
            assert!(
                (both_m - sum_m).abs() <= SUPERPOSITION_TOLERANCE * scale_m,
                "at cell ({i}, {j}) the composite run gives {both_m} m against {sum_m} m for the \
                 sum of the two runs"
            );
        }
    }
}

#[test]
fn a_burst_is_sampled_onto_the_zonal_faces_like_any_other_scenario() {
    // The burst reaches the solver through the same discretisation as the
    // trades: `τx` on the east/west faces at those faces' own positions.
    let basin = equatorial_basin(20, 4);
    let burst = pacific_burst();
    let t_s = BURST_PEAK_TIME_S;

    let field = WindStressField::sampled(basin, &burst, t_s);

    let nx = basin.grid().nx();
    for j in 0..field.tau_x_pa().ny() {
        let y_m = basin.y_of_row_m(U_STAGGERING, j);
        for i in 1..nx {
            let x_m = basin.x_of_column_m(U_STAGGERING, i);
            let expected_pa = burst.stress(x_m, y_m, t_s).0;
            let sampled_pa = *field.tau_x_pa().get(i, j).expect("an interior face");
            assert!(
                (sampled_pa - expected_pa).abs() <= ROUNDING_TOLERANCE * expected_pa.abs().max(1.0),
                "τx at face ({i}, {j}) is {sampled_pa} Pa, expected {expected_pa} Pa"
            );
        }
    }
}
