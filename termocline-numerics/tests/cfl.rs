//! Acceptance tests for T-01.3 — CFL-based timestep selection.
//!
//! Every expected number here is derived on paper from the stability analysis
//! of the scheme ADR-0003 fixes — centred differences on the C-grid advanced
//! by classic RK4 — and never from running the code.
//!
//! The derivation, once, so each test below can refer to it:
//!
//! A plane wave `exp(i(kx + ly))` sampled on the C-grid turns the centred
//! difference `(f[i+1] − f[i]) / dx` into a multiplication by
//! `2i·sin(k·dx/2)/dx`, so the linear gravity-wave operator of the
//! shallow-water system has purely imaginary eigenvalues `±i·c·κ` with
//! `κ = √( (2 sin(k dx/2)/dx)² + (2 sin(l dy/2)/dy)² )`. The grid-scale mode
//! (`k dx = l dy = π`) maximises it at `κ_max = 2·√(1/dx² + 1/dy²)`.
//!
//! Classic RK4 is stable on the imaginary axis for `|λ·dt| ≤ 2√2` (the
//! standard result, re-derived independently in
//! `the_stability_limit_is_the_rk4_imaginary_axis_bound` below). Combining
//! the two:
//!
//! ```text
//! dt_stable = 2√2 / (2·c·√(1/dx² + 1/dy²))
//! ```
//!
//! and `max_stable_dt` returns `CFL_SAFETY_FACTOR` times that.

use termocline_numerics::{
    max_stable_dt, CflError, Spacing, WaveSpeed, CFL_SAFETY_FACTOR, RK4_IMAGINARY_AXIS_LIMIT,
};

/// Relative tolerance for comparing a returned `dt` against a value worked out
/// by hand.
///
/// Both sides evaluate the same closed-form expression in `f64`; they differ
/// only by the handful of ulps the two orderings of square roots and divisions
/// accumulate. `1e-12` is roughly 4500 ulps at these magnitudes — orders of
/// magnitude tighter than any change of formula or safety factor, which is
/// what the test exists to catch.
const RELATIVE_TOLERANCE: f64 = 1e-12;

fn assert_close(actual: f64, expected: f64, what: &str) {
    let relative = (actual - expected).abs() / expected.abs();
    assert!(
        relative <= RELATIVE_TOLERANCE,
        "{what}: got {actual}, expected {expected} (relative error {relative:e})"
    );
}

fn spacing(dx_m: f64, dy_m: f64) -> Spacing {
    Spacing::new(dx_m, dy_m).expect("positive spacing")
}

fn wave_speed(m_per_s: f64) -> WaveSpeed {
    WaveSpeed::new(m_per_s).expect("positive wave speed")
}

// --- Acceptance criterion: unit tests pin the formula and the safety factor.

#[test]
fn the_stability_limit_is_the_rk4_imaginary_axis_bound() {
    // Independent check of the constant the formula rests on. The RK4
    // amplification factor for `y' = λy` is the degree-4 truncation of the
    // exponential, `R(z) = 1 + z + z²/2 + z³/6 + z⁴/24`. On the imaginary axis
    // `z = iθ` that is
    //     Re = 1 − θ²/2 + θ⁴/24,  Im = θ − θ³/6.
    // `|R| = 1` exactly at θ = 2√2 — the endpoint of RK4's imaginary-axis
    // stability interval, and the number the CFL bound divides by.
    let amplification = |theta: f64| {
        let real = 1.0 - theta * theta / 2.0 + theta.powi(4) / 24.0;
        let imaginary = theta - theta.powi(3) / 6.0;
        real.hypot(imaginary)
    };

    assert_close(
        RK4_IMAGINARY_AXIS_LIMIT,
        2.0 * std::f64::consts::SQRT_2,
        "RK4 imaginary-axis limit",
    );
    // Machine-epsilon-scale bound: |R(i·2√2)| = 1 is an identity, so the only
    // departure is roundoff in evaluating a handful of powers.
    assert!(
        (amplification(RK4_IMAGINARY_AXIS_LIMIT) - 1.0).abs() < 1e-14,
        "|R(i·2√2)| = {}, expected 1",
        amplification(RK4_IMAGINARY_AXIS_LIMIT)
    );
    assert!(
        amplification(RK4_IMAGINARY_AXIS_LIMIT * 0.99) < 1.0,
        "RK4 must be stable just inside the limit"
    );
    assert!(
        amplification(RK4_IMAGINARY_AXIS_LIMIT * 1.01) > 1.0,
        "RK4 must be unstable just outside the limit"
    );
}

#[test]
fn the_safety_factor_is_the_documented_one() {
    // Pinned so that changing the margin is a deliberate, reviewed edit rather
    // than a side effect. It is dimensionless and strictly inside (0, 1]: a
    // factor above 1 would put the returned `dt` past the stability boundary.
    assert!(
        (0.0..=1.0).contains(&CFL_SAFETY_FACTOR),
        "safety factor {CFL_SAFETY_FACTOR} must be a margin, not an amplifier"
    );
    assert_eq!(CFL_SAFETY_FACTOR, 0.8);
}

#[test]
fn the_formula_matches_the_hand_worked_isotropic_case() {
    // dx = dy = 100 km, c = 2 m/s (a plausible Kelvin wave speed for
    // g' ≈ 0.03 m/s², H ≈ 150 m). Then κ_max = 2·√2/10⁵ m⁻¹ and
    //     dt = 0.8 · 2√2 / (2 · 2 · √2/10⁵) = 0.8 · 10⁵/2 = 40 000 s.
    let dt_s = max_stable_dt(spacing(100_000.0, 100_000.0), wave_speed(2.0));
    assert_close(dt_s, 40_000.0, "isotropic 100 km grid at c = 2 m/s");
}

#[test]
fn the_formula_uses_both_axes_when_the_spacing_is_anisotropic() {
    // dx = 100 km, dy = 50 km, c = 2 m/s. Now
    //     √(1/dx² + 1/dy²) = √(5)/10⁵ m⁻¹
    // so dt = 0.8 · 2√2 · 10⁵ / (4·√5) = 4·10⁴·√(2/5) = 25 298.221 281 347 04 s.
    // A formula that looked only at the smaller spacing would return
    // 20 000 s; one that looked only at the larger, 40 000 s. Both are wrong,
    // and this case separates them.
    let dt_s = max_stable_dt(spacing(100_000.0, 50_000.0), wave_speed(2.0));
    assert_close(dt_s, 25_298.221_281_347_04, "100 km by 50 km grid");
}

#[test]
fn the_bound_is_inversely_proportional_to_the_wave_speed() {
    // The CFL number c·dt/dx is what the bound holds fixed, so doubling the
    // fastest wave speed must halve the timestep exactly.
    let slow = max_stable_dt(spacing(80_000.0, 80_000.0), wave_speed(1.5));
    let fast = max_stable_dt(spacing(80_000.0, 80_000.0), wave_speed(3.0));
    assert_close(fast, slow / 2.0, "doubling c halves dt");
}

#[test]
fn the_bound_shrinks_linearly_as_the_grid_is_refined() {
    // A CFL limit is first order in the spacing: over a sequence of
    // resolutions the ratio dt/dx must be constant, not merely small. Three
    // resolutions, each a halving of the last.
    let resolutions_m = [200_000.0, 100_000.0, 50_000.0];
    let courant: Vec<f64> = resolutions_m
        .iter()
        .map(|&dx_m| {
            let dt_s = max_stable_dt(spacing(dx_m, dx_m), wave_speed(2.5));
            dt_s / dx_m
        })
        .collect();
    for ratio in &courant[1..] {
        assert_close(*ratio, courant[0], "dt/dx across a halving of the grid");
    }
}

#[test]
fn the_returned_timestep_is_the_safety_factor_times_the_raw_stability_bound() {
    // The raw bound for this case, from the derivation in the module comment:
    //     2√2 / (2 · 2 · √2/10⁵) = 50 000 s.
    let raw_bound_s = 50_000.0;
    let dt_s = max_stable_dt(spacing(100_000.0, 100_000.0), wave_speed(2.0));
    assert_close(
        dt_s,
        CFL_SAFETY_FACTOR * raw_bound_s,
        "safety factor applied to the raw bound",
    );
}

// --- Acceptance criterion: a deliberately-too-large `dt` produces a clear,
// actionable error rather than running and producing garbage.

#[test]
fn a_timestep_past_the_cfl_bound_is_refused_with_the_value_and_the_bound() {
    let grid_spacing = spacing(100_000.0, 100_000.0);
    let c = wave_speed(2.0);
    let max_s = max_stable_dt(grid_spacing, c); // 40 000 s, worked out above.
    let requested_s = 90_000.0;

    let err = termocline_numerics::check_timestep(requested_s, grid_spacing, c)
        .expect_err("90 000 s is more than twice the CFL-stable maximum");
    assert_eq!(
        err,
        CflError::TimestepExceedsCfl {
            requested_s,
            max_stable_s: max_s,
        }
    );

    // Actionable per CODING_STANDARDS.md: the message names the offending
    // value and the bound it violated, and there is no silent clamping — the
    // call fails rather than substituting `max_s`.
    let message = err.to_string();
    assert!(message.contains("90000"), "{message}");
    assert!(message.contains("40000"), "{message}");
}

#[test]
fn a_timestep_inside_the_cfl_bound_is_accepted_up_to_the_bound_itself() {
    let grid_spacing = spacing(100_000.0, 100_000.0);
    let c = wave_speed(2.0);
    let max_s = max_stable_dt(grid_spacing, c);

    assert_eq!(
        termocline_numerics::check_timestep(max_s, grid_spacing, c),
        Ok(())
    );
    assert_eq!(
        termocline_numerics::check_timestep(max_s / 10.0, grid_spacing, c),
        Ok(())
    );
}

#[test]
fn a_timestep_that_is_not_a_positive_duration_is_refused() {
    let grid_spacing = spacing(100_000.0, 100_000.0);
    let c = wave_speed(2.0);
    for bad_s in [0.0, -60.0, f64::NAN, f64::INFINITY] {
        let err = termocline_numerics::check_timestep(bad_s, grid_spacing, c)
            .expect_err("a run needs a finite, positive dt");
        assert!(
            matches!(err, CflError::TimestepNotPositive { .. }),
            "dt = {bad_s} gave {err:?}"
        );
    }
}

#[test]
fn a_wave_speed_that_is_not_a_positive_speed_is_refused() {
    // c = √(g'H) is positive for any physical scenario; zero would make the
    // CFL bound infinite and silently disable the check.
    for bad in [0.0, -2.0, f64::NAN, f64::INFINITY] {
        let err = WaveSpeed::new(bad).expect_err("the fastest wave speed must be positive");
        // `bad` is compared by pattern rather than by `==` because one of the
        // cases is NaN, which is unequal to itself.
        assert!(
            matches!(err, CflError::WaveSpeedNotPositive { value_m_per_s } if value_m_per_s.to_bits() == bad.to_bits()),
            "c = {bad} gave {err:?}"
        );
    }
    let message = WaveSpeed::new(-2.0)
        .expect_err("negative speed")
        .to_string();
    assert!(message.contains("-2"), "{message}");
}
