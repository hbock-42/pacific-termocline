//! Acceptance tests for T-01.2 — the generic RK4 integrator.
//!
//! Every expected value here comes from an analytic solution of the test ODE
//! or from the Runge-Kutta tableau worked out by hand, never from running the
//! integrator. The central criterion is the *order* of the scheme: RK4 has a
//! global truncation error `O(dt^4)`, so halving `dt` must divide the error by
//! roughly 16. Each convergence test measures that ratio across three
//! resolutions rather than checking a single error against a fixed threshold.

use engine::{Rk4, StateVector};

/// Acceptable window on the measured convergence order.
///
/// The exact value is 4 (RK4's global order). The window allows for the fact
/// that a finite `dt` is only asymptotically in the `O(dt^4)` regime — the
/// leading error term carries `O(dt^5)` corrections — and for floating-point
/// roundoff. It is far too tight to admit a 3rd- or 5th-order scheme, which is
/// what it exists to rule out.
const ORDER_WINDOW: std::ops::RangeInclusive<f64> = 3.8..=4.2;

/// Estimated order of accuracy from two errors measured at `dt` and `dt / 2`.
///
/// If `error ∝ dt^p` then `error(dt) / error(dt/2) = 2^p`, so
/// `p = log2(coarse / fine)`.
fn observed_order(coarse_error: f64, fine_error: f64) -> f64 {
    (coarse_error / fine_error).log2()
}

/// Integrate `rhs` from `t = 0` to `t = t_end_s` seconds in `steps` uniform
/// steps.
fn integrate<S, F>(initial: &S, t_end_s: f64, steps: usize, mut rhs: F) -> S
where
    S: StateVector,
    F: FnMut(&S, f64, &mut S),
{
    let dt = t_end_s / steps as f64;
    let mut integrator = Rk4::new(initial);
    let mut state = initial.clone();
    for n in 0..steps {
        integrator.step(&mut state, n as f64 * dt, dt, &mut rhs);
    }
    state
}

fn assert_fourth_order(errors: [f64; 3], what: &str) {
    for pair in errors.windows(2) {
        let order = observed_order(pair[0], pair[1]);
        assert!(
            ORDER_WINDOW.contains(&order),
            "{what}: measured order {order} outside {:?} (errors {:e} then {:e})",
            ORDER_WINDOW,
            pair[0],
            pair[1]
        );
    }
}

// --- Acceptance criterion: 4th-order convergence on ODEs with analytic
// solutions, measured as `dt` is halved. ---

/// Exponential decay `y' = -λy`, whose analytic solution is `y = y0 e^{-λt}`.
#[test]
fn exponential_decay_converges_at_fourth_order() {
    // Rate and horizon chosen so that λ·t_end = O(1): the solution decays by a
    // factor e, far from both the trivial regime and stiffness.
    const DECAY_RATE_PER_S: f64 = 1.0;
    const T_END_S: f64 = 1.0;
    const Y0: f64 = 2.0;

    let exact = Y0 * (-DECAY_RATE_PER_S * T_END_S).exp();

    let errors = [10_usize, 20, 40].map(|steps| {
        let y = integrate(&Y0, T_END_S, steps, |y: &f64, _t, out: &mut f64| {
            *out = -DECAY_RATE_PER_S * y;
        });
        (y - exact).abs()
    });

    assert_fourth_order(errors, "exponential decay");
}

/// Simple harmonic oscillator `x'' = -ω²x`, written as the first-order system
/// `[x, v]' = [v, -ω²x]`, with analytic solution `x = cos(ωt)`,
/// `v = -ω sin(ωt)`.
#[test]
fn harmonic_oscillator_converges_at_fourth_order() {
    const ANGULAR_FREQUENCY_PER_S: f64 = 1.0;
    // One full period: the phase error of a time integrator is what shows up
    // over a whole oscillation, and it is the term that must vanish at 4th
    // order.
    const T_END_S: f64 = std::f64::consts::TAU;

    let initial = [1.0_f64, 0.0_f64];
    let exact = [
        (ANGULAR_FREQUENCY_PER_S * T_END_S).cos(),
        -ANGULAR_FREQUENCY_PER_S * (ANGULAR_FREQUENCY_PER_S * T_END_S).sin(),
    ];

    let errors = [20_usize, 40, 80].map(|steps| {
        let state = integrate(&initial, T_END_S, steps, |s: &[f64; 2], _t, out| {
            out[0] = s[1];
            out[1] = -ANGULAR_FREQUENCY_PER_S * ANGULAR_FREQUENCY_PER_S * s[0];
        });
        ((state[0] - exact[0]).powi(2) + (state[1] - exact[1]).powi(2)).sqrt()
    });

    assert_fourth_order(errors, "harmonic oscillator");
}

/// A time-dependent right-hand side `y' = cos(t)`, analytic solution
/// `y = sin(t)`.
///
/// The forced shallow-water problem is non-autonomous — wind stress is a
/// function of `t` — so the stage times `t`, `t + dt/2`, `t + dt` matter. An
/// integrator that evaluated every stage at the step's start time would still
/// converge, but only at 1st order, which this test rules out.
#[test]
fn time_dependent_forcing_converges_at_fourth_order() {
    const T_END_S: f64 = 2.0;

    let exact = T_END_S.sin();

    let errors = [8_usize, 16, 32].map(|steps| {
        let y = integrate(&0.0_f64, T_END_S, steps, |_y: &f64, t, out: &mut f64| {
            *out = t.cos();
        });
        (y - exact).abs()
    });

    assert_fourth_order(errors, "time-dependent forcing");
}

// --- Pinning the tableau itself, so a scheme with the right order but the
// wrong coefficients cannot pass. ---

/// One step of classic RK4 on `y' = y` from `y = 1` with `dt = 1` reproduces
/// the exponential series truncated after the `dt⁴` term:
/// `1 + 1 + 1/2 + 1/6 + 1/24`. Derived from the tableau by hand, not measured.
#[test]
fn one_step_matches_the_hand_worked_tableau() {
    let mut integrator = Rk4::new(&1.0_f64);
    let mut y = 1.0_f64;
    integrator.step(&mut y, 0.0, 1.0, &mut |y: &f64, _t, out: &mut f64| {
        *out = *y;
    });

    let expected = 1.0 + 1.0 + 1.0 / 2.0 + 1.0 / 6.0 + 1.0 / 24.0;
    // Both sides are short sums of order-one values; the bound is a few
    // machine epsilons, i.e. roundoff only.
    assert!(
        (y - expected).abs() < 8.0 * f64::EPSILON,
        "one RK4 step gave {y}, hand-worked tableau gives {expected}"
    );
}

/// RK4 applied to `y' = g(t)` reduces to Simpson's rule, which integrates
/// cubics exactly. So `y' = 3t²` from `y(0) = 0` must reproduce `y = t³` to
/// roundoff at any resolution — a property no lower-order scheme has.
#[test]
fn cubic_quadrature_is_exact() {
    const T_END_S: f64 = 3.0;
    let exact = T_END_S.powi(3);

    let y = integrate(&0.0_f64, T_END_S, 4, |_y: &f64, t, out: &mut f64| {
        *out = 3.0 * t * t;
    });

    // Exact in exact arithmetic; the bound is accumulated roundoff over four
    // steps, relative to a value of order 27.
    assert!(
        (y - exact).abs() < 32.0 * f64::EPSILON * exact,
        "Simpson-exact cubic gave {y}, analytic value is {exact}"
    );
}

// --- The integrator is generic over the state type, per the ticket. ---

/// The same integrator advances a scalar and a two-component system with no
/// change to the integrator itself.
#[test]
fn the_same_integrator_serves_any_state_vector() {
    // `y' = 0` leaves any state untouched, whatever its shape.
    let scalar = integrate(&5.0_f64, 1.0, 3, |_y: &f64, _t, out: &mut f64| *out = 0.0);
    assert_eq!(scalar, 5.0);

    let vector = integrate(&[1.0_f64, -2.0], 1.0, 3, |_s: &[f64; 2], _t, out| {
        *out = [0.0, 0.0];
    });
    assert_eq!(vector, [1.0, -2.0]);
}

/// Determinism: the same problem integrated twice gives bit-identical output
/// (CODING_STANDARDS.md § Correctness and failure).
#[test]
fn repeated_runs_are_bit_identical() {
    let run = || {
        integrate(&1.0_f64, 1.0, 17, |y: &f64, t, out: &mut f64| {
            *out = -y * t;
        })
    };
    assert_eq!(run().to_bits(), run().to_bits());
}
