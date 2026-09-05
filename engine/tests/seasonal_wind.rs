//! Acceptance tests for T-03.2 — the seasonal trade-wind scenario.
//!
//! `SeasonalTradeWinds` is the steady field of T-03.1 multiplied by an annual
//! harmonic, `1 + a·cos(2π·(t − t_peak)/T_year)`, with the relative amplitude
//! `a` and the phase — expressed as the instant `t_peak` at which the alizés
//! are strongest — both scenario input.
//!
//! The ticket's acceptance criterion is that *the output field's time series
//! at a fixed point shows the expected annual periodicity, via FFT or
//! peak-detection on a short run*. It is checked three ways, in increasing
//! order of how much of the engine has to be right:
//!
//! - [`the_sampled_time_series_peaks_once_a_tropical_year`] — peak detection
//!   on `τx` at one interior C-grid face, sampled over three years.
//! - [`the_sampled_time_series_carries_a_single_annual_spectral_line`] — the
//!   discrete Fourier transform of that same series. For
//!   `τ₀·(1 + a·cos(2π(t − t_peak)/T))` sampled at `N` equally spaced instants
//!   spanning a whole number `Y` of periods, the transform is known in closed
//!   form: the mean is `τ₀`, the component at `Y` cycles per window has
//!   amplitude `a·|τ₀|`, and *every* other component is exactly zero. That is
//!   the expected value, written out here rather than read back from the code.
//! - [`a_seasonally_forced_basin_oscillates_at_the_annual_period`] — the same
//!   transform applied to `h` at one cell of a basin actually integrated
//!   forward under the scenario. The v1 core is linear
//!   (CODING_STANDARDS.md § Scope guards), so a forcing carrying only the
//!   frequencies `{0, 1/T}` drives a response carrying only `{0, 1/T}`: a
//!   linear system generates no harmonics. Any power the run shows at `2/T` or
//!   `3/T` is numerical, and the test bounds it.
//!
//! No tolerance below was obtained by running the code; each is derived in the
//! comment beside it.

use std::f64::consts::PI;

use engine::{
    max_stable_dt, Basin, BetaPlane, Grid, OceanState, PhysicalParams, SeasonalTradeWinds, Solver,
    Spacing, SteadyTradeWinds, WaveSpeed, WindStress, WindStressError, WindStressField,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981).
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;
/// The equatorial beta-plane gradient, in m⁻¹s⁻¹ — `CONTEXT.md`, *Beta-plane*.
const BETA_PER_M_PER_S: f64 = engine::EQUATORIAL_BETA_PER_M_PER_S;
/// Reference seawater density `ρ₀`, in kg/m³ — `CONTEXT.md` and Gill, appendix 3.
const REFERENCE_DENSITY_KG_PER_M3: f64 = engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3;

/// Rayleigh damping `r` the forced run uses, in s⁻¹: an `e`-folding time of
/// about 11.6 days. Far stronger than the equatorial Pacific's own damping,
/// for the reason `rayleigh_damping.rs` spells out — the basin has to settle
/// into its periodic state inside a run of CFL-admissible steps.
const STRONG_DAMPING_PER_S: f64 = 1.0e-6;

/// Zonal wind stress `τ₀` of the trade-wind scenarios, in Pa. Easterly
/// trade-wind stress is `τx < 0` (`CONTEXT.md`), and 0.05 Pa is the observed
/// scale of the equatorial Pacific's mean zonal stress.
const TRADE_WIND_STRESS_PA: f64 = -0.05;

/// One solar day, in seconds.
const DAY_S: f64 = 86_400.0;
/// The tropical year, in seconds: 365.2422 mean solar days, the *Astronomical
/// Almanac*'s mean tropical year. Every expected value below is written in
/// terms of this rather than of the engine's own constant, so a change to the
/// period would have to be made here too before the suite agreed with it.
const YEAR_S: f64 = 365.2422 * DAY_S;

/// Relative amplitude `a` of the seasonal modulation. The equatorial Pacific's
/// zonal stress varies by a few tens of percent over the year; the value is
/// scenario input, and nothing below depends on it beyond its being a strict
/// fraction — the closed forms are written in terms of `a` itself.
const SEASONAL_RELATIVE_AMPLITUDE: f64 = 0.4;
/// Phase of the modulation, in seconds: the instant the alizés are strongest,
/// a quarter of a year into the run. Chosen off both the start and the middle
/// of the sampling windows below so that a peak-detection or transform test
/// could not pass by landing on a boundary.
const PEAK_TIME_S: f64 = 0.25 * YEAR_S;

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*).
const BASIN_LX_M: f64 = 1.0e7;
/// Cell height of the test basin, in metres.
const BASIN_DY_M: f64 = 1.0e5;
/// Zonal cell count of the test basin.
const BASIN_NX: usize = 40;
/// Meridional cell count of the test basin. Narrow, so the run stays close to
/// the equatorial waveguide the scenario is about and well inside the
/// rotation bound on the timestep.
const BASIN_NY: usize = 4;

/// Relative slack allowed where a check is exact in exact arithmetic: a few
/// tens of ulps of `f64` (ε ≈ 2.2×10⁻¹⁶) for the handful of operations per
/// point the expression costs.
const ROUNDING_TOLERANCE: f64 = 1.0e-14;

/// Samples per year of the sampled-field time series. A power of two so that
/// a whole number of years is a whole number of samples, which is what makes
/// the discrete transform below land on exact frequencies.
const FIELD_SAMPLES_PER_YEAR: usize = 512;
/// Years the sampled-field time series spans.
const FIELD_SAMPLE_YEARS: usize = 3;

/// Steps per year of the forced run, and therefore samples per year of its
/// `h` time series. `T_year/1024 = 3.08×10⁴ s` is inside the CFL bound of
/// `3.8×10⁴ s` for this basin, and dividing the year exactly is what makes the
/// discrete solution exactly periodic over 1024 steps and the transform's
/// frequencies exact.
const RUN_STEPS_PER_YEAR: usize = 1024;
/// Years of spin-up discarded before the `h` time series is recorded. The
/// slowest transient decays like `exp(−r·t)`, so after three years it is down
/// by `exp(−r·3·T) = exp(−94)`, far below every tolerance here.
const RUN_SPIN_UP_YEARS: usize = 3;
/// Years of `h` recorded, and therefore the number of cycles the fundamental
/// occupies in the transform window.
const RUN_SAMPLE_YEARS: usize = 2;

/// The equatorial-Pacific parameter set at a given Rayleigh damping.
fn pacific_params(rayleigh_damping_per_s: f64) -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        rayleigh_damping_per_s,
        BETA_PER_M_PER_S,
        REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

/// The equatorial deformation radius `Le = √(c/β)`, in metres — the meridional
/// scale over which equatorial waves decay away from the equator
/// (`CONTEXT.md`), and the natural width for a wind field meant to drive that
/// waveguide.
fn equatorial_deformation_radius_m(params: PhysicalParams) -> f64 {
    (params.kelvin_wave_speed_m_per_s() / params.beta_per_m_per_s()).sqrt()
}

/// A basin [`BASIN_LX_M`] wide and [`BASIN_NY`] cells tall, centred on the
/// equator.
fn equatorial_basin() -> Basin {
    let grid = Grid::new(BASIN_NX, BASIN_NY).expect("extents are non-zero");
    let spacing = Spacing::new(BASIN_LX_M / BASIN_NX as f64, BASIN_DY_M)
        .expect("a basin spanned by whole cells has positive spacing");
    Basin::centered_on_equator(grid, spacing)
}

/// The uniform trade winds the seasonal scenarios below modulate.
fn steady_uniform() -> SteadyTradeWinds {
    SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA).expect("an easterly stress is a trade wind")
}

/// The seasonal scenario under test: [`steady_uniform`] breathing at
/// [`SEASONAL_RELATIVE_AMPLITUDE`], strongest at [`PEAK_TIME_S`].
fn seasonal_uniform() -> SeasonalTradeWinds {
    SeasonalTradeWinds::new(steady_uniform(), SEASONAL_RELATIVE_AMPLITUDE, PEAK_TIME_S)
        .expect("a fractional amplitude at a finite phase is a season")
}

/// The same season over trade winds that decay away from the equator at
/// `decay_scale_m`: the profile that shows whether the modulation preserves
/// the field's shape in `y` while changing its strength.
fn seasonal_decaying(decay_scale_m: f64) -> SeasonalTradeWinds {
    let steady = SteadyTradeWinds::with_meridional_decay(TRADE_WIND_STRESS_PA, decay_scale_m)
        .expect("an easterly stress with a positive decay scale");
    SeasonalTradeWinds::new(steady, SEASONAL_RELATIVE_AMPLITUDE, PEAK_TIME_S)
        .expect("a fractional amplitude at a finite phase is a season")
}

/// The annual harmonic `1 + a·cos(2π(t − t_peak)/T_year)`, evaluated here from
/// the ticket's formula rather than read back from the engine.
fn expected_modulation(t_s: f64) -> f64 {
    1.0 + SEASONAL_RELATIVE_AMPLITUDE * (2.0 * PI * (t_s - PEAK_TIME_S) / YEAR_S).cos()
}

// --- The harmonic itself. ---

#[test]
fn the_period_of_the_season_is_the_tropical_year() {
    // Which year the season follows is a modelling decision, not an
    // implementation detail: it fixes what "annual" means for every scenario
    // and for T-03.4's config files. This pins it to the *Astronomical
    // Almanac*'s mean tropical year of 365.2422 days, so that adopting the
    // sidereal year or the calendar's 365 has to be a deliberate change to a
    // stated value rather than a quiet one.
    assert_eq!(engine::TROPICAL_YEAR_S, YEAR_S);
}

#[test]
fn the_season_modulates_the_steady_field_by_the_annual_harmonic() {
    // The ticket's formula, point by point: the seasonal stress is the steady
    // stress at the same place times `1 + a·cos(2π(t − t_peak)/T)`, with the
    // harmonic evaluated in this file and the steady profile — a Gaussian in
    // `y` — evaluated by the field it wraps.
    let decay_scale_m = equatorial_deformation_radius_m(pacific_params(STRONG_DAMPING_PER_S));
    let seasonal = seasonal_decaying(decay_scale_m);

    for y_m in [0.0, decay_scale_m, -2.0 * decay_scale_m] {
        let scaled = y_m / decay_scale_m;
        let steady_pa = TRADE_WIND_STRESS_PA * (-scaled * scaled).exp();
        for eighths in 0..16 {
            let t_s = eighths as f64 * YEAR_S / 8.0;
            let expected_pa = steady_pa * expected_modulation(t_s);
            let (tau_x_pa, tau_y_pa) = seasonal.stress(BASIN_LX_M / 3.0, y_m, t_s);
            assert!(
                (tau_x_pa - expected_pa).abs() <= ROUNDING_TOLERANCE * expected_pa.abs(),
                "at y = {y_m} m, t = {t_s} s the stress is {tau_x_pa} Pa, expected {expected_pa} Pa"
            );
            assert_eq!(tau_y_pa, 0.0, "the alizés carry no meridional stress");
        }
    }
}

#[test]
fn the_alizes_are_strongest_at_the_configured_phase_and_weakest_half_a_year_later() {
    // What "amplitude and phase configurable" has to mean. At `t_peak` the
    // harmonic is `1 + a`, half a year later `1 − a`, and at the two quarter
    // points it is exactly one — the steady field's own value.
    let seasonal = seasonal_uniform();

    for (offset_s, expected_factor) in [
        (0.0, 1.0 + SEASONAL_RELATIVE_AMPLITUDE),
        (YEAR_S / 2.0, 1.0 - SEASONAL_RELATIVE_AMPLITUDE),
        (YEAR_S / 4.0, 1.0),
        (-YEAR_S / 4.0, 1.0),
    ] {
        let expected_pa = TRADE_WIND_STRESS_PA * expected_factor;
        let (tau_x_pa, _) = seasonal.stress(0.0, 0.0, PEAK_TIME_S + offset_s);
        assert!(
            (tau_x_pa - expected_pa).abs() <= ROUNDING_TOLERANCE * expected_pa.abs(),
            "{offset_s} s after the peak the stress is {tau_x_pa} Pa, expected {expected_pa} Pa"
        );
    }
    // "Strongest" is about easterly strength, and `τx < 0`: the peak is the
    // most negative value the year takes.
    let peak_pa = seasonal.stress(0.0, 0.0, PEAK_TIME_S).0;
    for twelfths in 1..12 {
        let t_s = PEAK_TIME_S + twelfths as f64 * YEAR_S / 12.0;
        assert!(
            seasonal.stress(0.0, 0.0, t_s).0 > peak_pa,
            "the alizés must be strongest at the configured phase, but {t_s} s beats it"
        );
    }
}

#[test]
fn the_season_repeats_every_tropical_year() {
    // Periodicity as an identity on the trait, before any sampling or
    // transform: the stress at `t` and at `t + n·T` is the same field.
    //
    // Not an exact equality: the cosine's argument is formed by dividing an
    // elapsed time by the period, so `n` years on, the argument carries a
    // relative rounding error of order `n·ε`. With `|d(cos)/dθ| ≤ 1` and
    // `θ ≈ 2πn`, that is an absolute error below `2π·n·ε ≈ 5×10⁻¹⁵` on the
    // harmonic for the `n ≤ 4` here, which [`ROUNDING_TOLERANCE`] covers.
    let seasonal = seasonal_uniform();

    for twelfths in 0..12 {
        let t_s = twelfths as f64 * YEAR_S / 12.0;
        let reference_pa = seasonal.stress(0.0, 0.0, t_s).0;
        for years in 1..=4 {
            let later_pa = seasonal.stress(0.0, 0.0, t_s + years as f64 * YEAR_S).0;
            assert!(
                (later_pa - reference_pa).abs() <= ROUNDING_TOLERANCE * reference_pa.abs(),
                "at t = {t_s} s the stress is {reference_pa} Pa but {years} years later it is \
                 {later_pa} Pa"
            );
        }
    }
}

#[test]
fn a_season_of_no_amplitude_is_the_steady_field() {
    // The `a → 0` limit, exactly: a scenario configured with no seasonal cycle
    // must reproduce T-03.1's field bit for bit, at every instant.
    let steady = steady_uniform();
    let seasonal = SeasonalTradeWinds::new(steady, 0.0, PEAK_TIME_S)
        .expect("no modulation at all is a degenerate but valid season");

    for t_s in [0.0, DAY_S, PEAK_TIME_S, YEAR_S, 10.0 * YEAR_S] {
        assert_eq!(seasonal.stress(0.0, 0.0, t_s), steady.stress(0.0, 0.0, t_s));
    }
}

#[test]
fn the_modulated_alizes_never_turn_westerly() {
    // Why the amplitude is required to be a fraction: at `a = 1` the harmonic
    // touches zero once a year and the basin goes calm, and beyond that the
    // "trade winds" would reverse. A westerly stress is T-03.3's wind burst,
    // not a season, and a scenario named for the alizés must not quietly
    // become one.
    for amplitude in [0.0, 0.25, SEASONAL_RELATIVE_AMPLITUDE, 1.0] {
        let seasonal = SeasonalTradeWinds::new(steady_uniform(), amplitude, PEAK_TIME_S)
            .expect("a fraction is a valid amplitude");
        for sample in 0..=FIELD_SAMPLES_PER_YEAR {
            let t_s = sample as f64 * YEAR_S / FIELD_SAMPLES_PER_YEAR as f64;
            let tau_x_pa = seasonal.stress(0.0, 0.0, t_s).0;
            assert!(
                tau_x_pa <= 0.0,
                "at a = {amplitude}, t = {t_s} s the stress is {tau_x_pa} Pa, which is westerly"
            );
        }
    }
}

#[test]
fn a_season_reports_the_steady_field_it_modulates() {
    // The scenario is the steady field plus a harmonic; a config writer, and
    // T-03.4's loader, must be able to read back what they asked for.
    let seasonal = seasonal_uniform();

    assert_eq!(seasonal.steady(), steady_uniform());
    assert_eq!(seasonal.relative_amplitude(), SEASONAL_RELATIVE_AMPLITUDE);
    assert_eq!(seasonal.peak_time_s(), PEAK_TIME_S);
}

#[test]
fn a_seasonal_amplitude_that_is_not_a_fraction_is_refused_by_name() {
    // Invalid scenario input is a `Result` naming the offending value, not a
    // panic and not a silently clamped amplitude (CODING_STANDARDS.md
    // § Correctness and failure).
    for value in [-0.1, 1.000_001, 2.0, f64::NAN, f64::INFINITY] {
        let error = SeasonalTradeWinds::new(steady_uniform(), value, PEAK_TIME_S)
            .expect_err("only a fraction of the steady field is a seasonal modulation");
        let WindStressError::ModulationNotAFraction {
            relative_amplitude: rejected,
        } = error
        else {
            panic!("expected the amplitude itself to be rejected, got {error}");
        };
        // Compared bitwise rather than with `==`, so that the NaN case checks
        // the value was carried through rather than trivially passing.
        assert_eq!(rejected.to_bits(), value.to_bits());
        let message = error.to_string();
        assert!(message.contains("relative_amplitude"), "{message}");
    }
}

#[test]
fn a_phase_that_is_not_an_instant_is_refused_by_name() {
    for value_s in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let error = SeasonalTradeWinds::new(steady_uniform(), SEASONAL_RELATIVE_AMPLITUDE, value_s)
            .expect_err("a phase must be a finite instant");
        let WindStressError::PhaseNotFinite {
            peak_time_s: rejected_s,
        } = error
        else {
            panic!("expected the phase itself to be rejected, got {error}");
        };
        assert_eq!(rejected_s.to_bits(), value_s.to_bits());
        let message = error.to_string();
        assert!(message.contains("peak_time_s"), "{message}");
    }
}

// --- The sampled field's time series: the acceptance criterion. ---

/// `τx` in pascals at one fixed interior east/west face of `basin`, sampled
/// from `wind` at `samples` instants spaced `interval_s` apart from zero.
///
/// The face is an interior one mid-basin, so the series is the stress the
/// solver would actually read at a degree of freedom rather than at a wall,
/// where a sampled field is zero by the rule of T-03.1 and would carry no
/// season at all.
fn sampled_tau_x_series(
    basin: Basin,
    wind: &impl WindStress,
    samples: usize,
    interval_s: f64,
) -> Vec<f64> {
    let (probe_i, probe_j) = (BASIN_NX / 2, BASIN_NY / 2);
    let mut field = WindStressField::calm(basin.grid());
    (0..samples)
        .map(|sample| {
            let t_s = sample as f64 * interval_s;
            field.sample(basin, wind, t_s);
            *field
                .tau_x_pa()
                .get(probe_i, probe_j)
                .expect("an interior east/west face")
        })
        .collect()
}

/// The indices at which `samples` has a strict interior local maximum.
fn local_maxima(samples: &[f64]) -> Vec<usize> {
    (1..samples.len().saturating_sub(1))
        .filter(|&index| samples[index] > samples[index - 1] && samples[index] > samples[index + 1])
        .collect()
}

/// Amplitude of the component of `samples` completing exactly `cycles` cycles
/// over the window, in the samples' own units.
///
/// The plain definition, `(2/N)·|Σ xₙ·e^(−2πi·k·n/N)|`, which for a window
/// spanning a whole number of periods returns the amplitude of a cosine at
/// that frequency and zero at every other integer `k ≥ 1`. Written out rather
/// than pulled from an FFT crate: `N` is a thousand-odd samples, so the O(N²)
/// sum costs nothing, and the engine gains no dependency for a test.
fn fourier_amplitude(samples: &[f64], cycles: usize) -> f64 {
    let n = samples.len();
    let (mut real, mut imaginary) = (0.0, 0.0);
    for (index, &value) in samples.iter().enumerate() {
        // Reduced modulo `n` before scaling, so the angle stays inside one
        // turn and the cosine is not asked to cancel a large argument.
        let angle = -2.0 * PI * ((cycles * index) % n) as f64 / n as f64;
        real += value * angle.cos();
        imaginary += value * angle.sin();
    }
    2.0 * real.hypot(imaginary) / n as f64
}

/// The mean of `samples`, the zero-frequency component.
fn mean(samples: &[f64]) -> f64 {
    samples.iter().sum::<f64>() / samples.len() as f64
}

#[test]
fn the_sampled_time_series_peaks_once_a_tropical_year() {
    // Peak detection, the first form the acceptance criterion offers. Three
    // years of `τx` at one face, sampled 512 times a year: the easterly
    // strength `−τx` must have exactly one maximum per year, the first at the
    // configured phase and each a tropical year after the last.
    let basin = equatorial_basin();
    let interval_s = YEAR_S / FIELD_SAMPLES_PER_YEAR as f64;
    let samples = FIELD_SAMPLE_YEARS * FIELD_SAMPLES_PER_YEAR + 1;
    let easterly_strength_pa: Vec<f64> =
        sampled_tau_x_series(basin, &seasonal_uniform(), samples, interval_s)
            .into_iter()
            .map(|tau_x_pa| -tau_x_pa)
            .collect();

    let peaks = local_maxima(&easterly_strength_pa);

    assert_eq!(
        peaks.len(),
        FIELD_SAMPLE_YEARS,
        "a {FIELD_SAMPLE_YEARS}-year series of an annual harmonic has {FIELD_SAMPLE_YEARS} \
         peaks, found {peaks:?}"
    );
    // A peak is resolved no better than the sampling interval, so each is
    // required to land within half an interval of where the harmonic puts it.
    // (It lands exactly on a sample here — `t_peak` is a quarter year and the
    // year is 512 samples — but the bound is the honest one.)
    for (year, &peak) in peaks.iter().enumerate() {
        let peak_time_s = peak as f64 * interval_s;
        let expected_s = PEAK_TIME_S + year as f64 * YEAR_S;
        assert!(
            (peak_time_s - expected_s).abs() <= interval_s / 2.0,
            "peak {year} is at {peak_time_s} s, expected {expected_s} s"
        );
    }
    for pair in peaks.windows(2) {
        let spacing_s = (pair[1] - pair[0]) as f64 * interval_s;
        assert!(
            (spacing_s - YEAR_S).abs() <= interval_s,
            "consecutive peaks are {spacing_s} s apart, expected a tropical year ({YEAR_S} s)"
        );
    }
}

#[test]
fn the_sampled_time_series_carries_a_single_annual_spectral_line() {
    // The transform, the criterion's other form, over a window of exactly
    // three years. For `τ₀·(1 + a·cos(2π(t − t_peak)/T))` the answer is known
    // in closed form: mean `τ₀`, amplitude `a·|τ₀|` at three cycles per window
    // — one per year — and exactly zero at every other frequency.
    let basin = equatorial_basin();
    let interval_s = YEAR_S / FIELD_SAMPLES_PER_YEAR as f64;
    let samples = FIELD_SAMPLE_YEARS * FIELD_SAMPLES_PER_YEAR;
    let series = sampled_tau_x_series(basin, &seasonal_uniform(), samples, interval_s);

    // The transform sums `N ≈ 1.5×10³` terms of magnitude at most
    // `|τ₀|(1 + a)`, so its rounding error is bounded by `N·ε·|τ₀|(1 + a)`,
    // about `5×10⁻¹⁶·|τ₀|`. A relative bound of 10⁻¹² is that rounded up by
    // three orders of magnitude.
    const SPECTRAL_ROUNDING_TOLERANCE: f64 = 1.0e-12;
    let scale_pa = TRADE_WIND_STRESS_PA.abs();

    let mean_pa = mean(&series);
    assert!(
        (mean_pa - TRADE_WIND_STRESS_PA).abs() <= SPECTRAL_ROUNDING_TOLERANCE * scale_pa,
        "the annual mean stress is {mean_pa} Pa, expected the steady {TRADE_WIND_STRESS_PA} Pa"
    );

    let annual_pa = fourier_amplitude(&series, FIELD_SAMPLE_YEARS);
    let expected_pa = SEASONAL_RELATIVE_AMPLITUDE * scale_pa;
    assert!(
        (annual_pa - expected_pa).abs() <= SPECTRAL_ROUNDING_TOLERANCE * scale_pa,
        "the annual line has amplitude {annual_pa} Pa, expected {expected_pa} Pa"
    );

    for cycles in 1..=(3 * FIELD_SAMPLE_YEARS) {
        if cycles == FIELD_SAMPLE_YEARS {
            continue;
        }
        let amplitude_pa = fourier_amplitude(&series, cycles);
        assert!(
            amplitude_pa <= SPECTRAL_ROUNDING_TOLERANCE * scale_pa,
            "the series carries {amplitude_pa} Pa at {cycles} cycles over the window; only the \
             annual line may be non-zero"
        );
    }
}

#[test]
fn the_sampled_field_breathes_over_the_whole_basin() {
    // The season multiplies the field everywhere, not just at the one face the
    // series above watches: every face must scale by the harmonic, the basin's
    // wall lines included. What the solver does with a stress at the coast is
    // the no-normal-flow condition of T-04.2's business, not the sampler's.
    let basin = equatorial_basin();
    let decay_scale_m = equatorial_deformation_radius_m(pacific_params(STRONG_DAMPING_PER_S));
    let seasonal = seasonal_decaying(decay_scale_m);
    let (nx, ny) = (basin.grid().nx(), basin.grid().ny());

    let reference = WindStressField::sampled(basin, &seasonal.steady(), 0.0);
    for eighths in 0..8 {
        let t_s = eighths as f64 * YEAR_S / 8.0;
        let modulation = expected_modulation(t_s);
        let field = WindStressField::sampled(basin, &seasonal, t_s);

        for j in 0..field.tau_x_pa().ny() {
            for i in 0..=nx {
                let expected_pa =
                    modulation * *reference.tau_x_pa().get(i, j).expect("an in-bounds face");
                let sampled_pa = *field.tau_x_pa().get(i, j).expect("an in-bounds face");
                assert!(
                    (sampled_pa - expected_pa).abs() <= ROUNDING_TOLERANCE * expected_pa.abs(),
                    "τx at face ({i}, {j}) at t = {t_s} s is {sampled_pa} Pa, expected \
                     {expected_pa} Pa"
                );
            }
        }
        for i in 0..field.tau_y_pa().nx() {
            for j in 0..=ny {
                assert_eq!(
                    *field.tau_y_pa().get(i, j).expect("an in-bounds face"),
                    0.0,
                    "the alizés carry no meridional stress"
                );
            }
        }
    }
}

// --- The basin's response: annual periodicity in the model's own output. ---

#[test]
fn a_seasonally_forced_basin_oscillates_at_the_annual_period() {
    // The acceptance criterion applied to a short run rather than to the
    // forcing alone: three years of spin-up, then two years of `h` recorded at
    // one cell, one sample per step.
    //
    // The v1 core is linear, so the response to a forcing carrying only the
    // frequencies `{0, 1/T}` carries only `{0, 1/T}`: a mean tilt and one
    // annual oscillation, with no harmonics generated. Three statements of
    // that are checked — the series repeats exactly one year on, the annual
    // line is a substantial fraction of the mean tilt, and the harmonics are
    // numerical noise.
    let basin = equatorial_basin();
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let seasonal = seasonal_uniform();
    let dt_s = YEAR_S / RUN_STEPS_PER_YEAR as f64;

    let wave_speed =
        WaveSpeed::new(params.kelvin_wave_speed_m_per_s()).expect("a positive wave speed");
    assert!(
        dt_s < max_stable_dt(basin.spacing(), wave_speed),
        "the run's timestep must divide the year and still be CFL-admissible"
    );
    let plane = BetaPlane::centered_on_equator(params, basin.spacing(), basin.grid());
    let mut solver = Solver::new(basin.grid(), basin.spacing(), params, plane, dt_s)
        .unwrap_or_else(|error| panic!("the test's own timestep must be admissible: {error}"));

    // The westernmost column on the equatorward side of the basin: where the
    // wind-driven tilt is deepest, and so where its seasonal breathing is
    // largest.
    let (probe_i, probe_j) = (0, BASIN_NY / 2);
    let mut state = OceanState::at_rest(basin.grid());
    let mut step = 0usize;
    for _ in 0..RUN_SPIN_UP_YEARS * RUN_STEPS_PER_YEAR {
        solver.step_forced_by(&mut state, step as f64 * dt_s, basin, &seasonal);
        step += 1;
    }
    let mut h_m = Vec::with_capacity(RUN_SAMPLE_YEARS * RUN_STEPS_PER_YEAR);
    for _ in 0..RUN_SAMPLE_YEARS * RUN_STEPS_PER_YEAR {
        h_m.push(*state.h().get(probe_i, probe_j).expect("an in-bounds cell"));
        solver.step_forced_by(&mut state, step as f64 * dt_s, basin, &seasonal);
        step += 1;
    }

    let mean_h_m = mean(&h_m);
    assert!(
        mean_h_m > 0.0,
        "the easterly alizés must leave the western thermocline deeper than the mean on \
         average, got {mean_h_m} m"
    );

    let annual_m = fourier_amplitude(&h_m, RUN_SAMPLE_YEARS);
    // The forcing's annual line is a fraction `a` of its mean, and the basin's
    // adjustment time `L/c = 3.7×10⁶ s` is a tenth of a year, so the response
    // is close to quasi-steady and its annual line must be a comparable
    // fraction of the mean tilt. Requiring it to exceed `a/2` of the mean is
    // the weak form of that: it says the season is visible in `h`, without
    // asserting a transfer function Epic 07 has not yet pinned down.
    assert!(
        annual_m >= 0.5 * SEASONAL_RELATIVE_AMPLITUDE * mean_h_m,
        "the annual oscillation of h is {annual_m} m against a mean tilt of {mean_h_m} m; the \
         season has barely reached the ocean"
    );

    // A linear system generates no harmonics, so every other line is
    // numerical. The residual spin-up transient is `exp(−r·3·T) ≈ 10⁻⁴¹` of
    // its initial size, and the transform's own rounding is `N·ε ≈ 5×10⁻¹³`
    // relative, so 10⁻⁹ of the annual line is that bounded well from above.
    const HARMONIC_TOLERANCE: f64 = 1.0e-9;
    for cycles in 1..=(4 * RUN_SAMPLE_YEARS) {
        if cycles == RUN_SAMPLE_YEARS {
            continue;
        }
        let amplitude_m = fourier_amplitude(&h_m, cycles);
        assert!(
            amplitude_m <= HARMONIC_TOLERANCE * annual_m,
            "h carries {amplitude_m} m at {cycles} cycles over the window against {annual_m} m \
             at the annual line; the linear core must generate no harmonics"
        );
    }

    // Periodicity in the run's own output, stated directly: a year on, to the
    // step, the thermocline is where it was. The bound is the same 10⁻⁹ of the
    // oscillation's amplitude, for the same reasons.
    for (index, &value_m) in h_m.iter().take(RUN_STEPS_PER_YEAR).enumerate() {
        let a_year_later_m = h_m[index + RUN_STEPS_PER_YEAR];
        assert!(
            (a_year_later_m - value_m).abs() <= HARMONIC_TOLERANCE * annual_m,
            "at step {index} h is {value_m} m, but a year later it is {a_year_later_m} m"
        );
    }
}
