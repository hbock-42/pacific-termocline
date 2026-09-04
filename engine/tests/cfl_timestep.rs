//! Acceptance tests for T-01.3 from the engine's side: the runtime check a run
//! passes through before it starts, and empirical confirmation that the
//! timestep the check accepts really does keep RK4 stable.
//!
//! The stability facts here are the analytic ones derived in
//! `termocline-numerics/tests/cfl.rs`; nothing is calibrated by running the
//! solver.

use engine::{check_timestep, max_stable_dt, CflError, Rk4, Spacing, WaveSpeed, CFL_SAFETY_FACTOR};

/// Cell spacing and wave speed used throughout: a 100 km square grid and a
/// 2 m/s Kelvin wave speed, for which the CFL-stable maximum is 40 000 s.
fn basin() -> (Spacing, WaveSpeed) {
    (
        Spacing::new(100_000.0, 100_000.0).expect("positive spacing"),
        WaveSpeed::new(2.0).expect("positive wave speed"),
    )
}

/// Angular frequency of the fastest mode the C-grid can carry, in rad/s.
///
/// The grid-scale mode of the centred-difference gravity-wave operator, from
/// the derivation in the numerics acceptance tests:
/// `ω = c · 2 · √(1/dx² + 1/dy²)`.
///
/// Deliberately written out here rather than taken from `termocline-numerics`:
/// this is the analytic system the test integrates, and a test that asked the
/// code under review for it could not catch the code getting it wrong.
fn fastest_mode_rad_per_s(spacing: Spacing, wave_speed: WaveSpeed) -> f64 {
    let dx_m = spacing.dx_m();
    let dy_m = spacing.dy_m();
    wave_speed.m_per_s() * 2.0 * (1.0 / (dx_m * dx_m) + 1.0 / (dy_m * dy_m)).sqrt()
}

/// Amplitude of that mode after `steps` RK4 steps of `dt_s`, starting from
/// unit amplitude.
///
/// The mode obeys `dy/dt = iωy`; written over the reals with
/// `y = a + ib` that is `[a, b]' = ω[−b, a]`, whose exact solution is a
/// rotation — amplitude 1 for all time. Any growth is the scheme's, not the
/// physics'.
fn amplitude_after(steps: usize, dt_s: f64, omega_rad_per_s: f64) -> f64 {
    let mut integrator = Rk4::new(&[0.0_f64, 0.0]);
    let mut state = [1.0_f64, 0.0];
    for n in 0..steps {
        integrator.step(
            &mut state,
            n as f64 * dt_s,
            dt_s,
            &mut |y: &[f64; 2], _t, out: &mut [f64; 2]| {
                out[0] = -omega_rad_per_s * y[1];
                out[1] = omega_rad_per_s * y[0];
            },
        );
    }
    state[0].hypot(state[1])
}

// --- Acceptance criterion: a deliberately-too-large `dt` produces a clear,
// actionable error rather than running and producing garbage.

#[test]
fn the_engine_refuses_an_unstable_timestep_before_starting_a_run() {
    let (spacing, wave_speed) = basin();
    let err = check_timestep(86_400.0, spacing, wave_speed)
        .expect_err("a one-day step is more than twice the 40 000 s CFL bound");

    let message = err.to_string();
    assert!(message.contains("86400"), "{message}");
    assert!(message.contains("40000"), "{message}");
    assert!(
        matches!(err, CflError::TimestepExceedsCfl { .. }),
        "{err:?}"
    );
}

#[test]
fn the_rejected_timestep_really_does_blow_the_run_up() {
    // The point of the refusal: had the run been allowed to proceed, it would
    // have produced garbage rather than a wrong-but-bounded answer. The
    // amplification per step at this timestep is |R(i·4.887)| ≈ 24, so twenty
    // steps of a mode whose amplitude should stay at 1 is already past 1e6 —
    // and few enough steps that the growth does not overflow `f64` into a NaN
    // the assertion could not read.
    let (spacing, wave_speed) = basin();
    let omega = fastest_mode_rad_per_s(spacing, wave_speed);
    let amplitude = amplitude_after(20, 86_400.0, omega);
    assert!(
        amplitude > 1e6,
        "an unstable step should diverge, got amplitude {amplitude}"
    );
}

#[test]
fn the_accepted_timestep_keeps_the_fastest_mode_bounded() {
    // The converse: at exactly the timestep the check accepts, the grid-scale
    // mode must not grow. `|R(iθ)| ≤ 1` for `θ ≤ 2√2` is the analytic
    // statement; 10 000 steps is long enough that any amplification per step
    // above 1 + 1e-6 would show up as growth well past the bound below.
    let (spacing, wave_speed) = basin();
    let omega = fastest_mode_rad_per_s(spacing, wave_speed);
    let dt_s = max_stable_dt(spacing, wave_speed);
    assert_eq!(check_timestep(dt_s, spacing, wave_speed), Ok(()));

    let amplitude = amplitude_after(10_000, dt_s, omega);
    assert!(
        amplitude <= 1.0,
        "the accepted timestep must not amplify the grid-scale mode, got {amplitude}"
    );
}

#[test]
fn the_safety_factor_leaves_real_margin_below_the_stability_boundary() {
    // The accepted timestep sits strictly inside the stability region, so
    // scaling it back up by more than the margin must cross the boundary.
    // With a safety factor of 0.8 the raw bound is dt/0.8, and 1.2× the raw
    // bound is unambiguously unstable: |R(i·1.2·2√2)| ≈ 3.2 per step, so forty
    // steps grow the mode by twenty orders of magnitude while staying inside
    // `f64`'s range.
    let (spacing, wave_speed) = basin();
    let omega = fastest_mode_rad_per_s(spacing, wave_speed);
    let raw_bound_s = max_stable_dt(spacing, wave_speed) / CFL_SAFETY_FACTOR;

    assert!(
        amplitude_after(40, raw_bound_s * 1.2, omega) > 10.0,
        "past the raw stability bound the mode must grow"
    );
    assert!(
        check_timestep(raw_bound_s * 1.2, spacing, wave_speed).is_err(),
        "and the check must refuse it"
    );
}

#[test]
fn the_refusal_message_reads_as_something_a_user_can_act_on() {
    // The whole point of the check is a message that says what to do next, so
    // the exact wording is pinned rather than left to drift.
    let (spacing, wave_speed) = basin();
    let message = check_timestep(86_400.0, spacing, wave_speed)
        .expect_err("a one-day step is unstable here")
        .to_string();
    assert_eq!(
        message,
        "dt is 86400 s, past the CFL-stable maximum of 40000 s for this grid spacing and wave \
         speed; the run would go unstable. Set dt to at most 40000 s, or coarsen the grid"
    );
}
