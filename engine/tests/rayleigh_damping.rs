//! Acceptance tests for T-02.4 — the Rayleigh damping terms of the 1.5-layer
//! reduced-gravity shallow-water right-hand side.
//!
//! The ticket folds `−r·u`, `−r·v` and `−r·h` into `shallow_water_rhs`, so the
//! equations under test are those of T-02.3 with one linear term added to each:
//!
//! ```text
//! ∂u/∂t = −g'·∂h/∂x + τx/(ρ₀·H) − r·u
//! ∂v/∂t = −g'·∂h/∂y + τy/(ρ₀·H) − r·v
//! ∂h/∂t = −H·(∂u/∂x + ∂v/∂y)    − r·h
//! ```
//!
//! Two analytic facts drive every quantitative check here, and neither comes
//! from running this code:
//!
//! - A spatially uniform state has no gradient and no divergence, so it is an
//!   eigenvector of the whole right-hand side with eigenvalue `−r`: each of
//!   `h`, `u` and `v` obeys `ẏ = −r·y` exactly and therefore decays as
//!   `y(0)·exp(−r·t)`. That is the single-mode test case the acceptance
//!   criteria ask for, and RK4 approximates it to fourth order in `dt`
//!   ([ADR-0003]). It is an eigenvector of *this* right-hand side, which does
//!   not yet carry the Coriolis term of T-02.2: `−f·v` with `f = β·y` varies
//!   down the basin, so folding [`CoriolisTerm`](engine::CoriolisTerm) into
//!   the same evaluation will cost the uniform state its eigenvector status
//!   and this file its single-mode case. The energy identity below survives
//!   that — rotation is skew too — and is the check to lean on when it
//!   happens.
//! - Damping `h`, `u` and `v` at the same rate makes the system
//!   `ẋ = (L − r)·x`, where `L` is the conservative pressure-gradient and
//!   continuity pair. `L` is skew in the discrete energy
//!   `E = (g'/2)·Σh² + (H/2)·Σ(u² + v²)` — summation by parts on the C-grid,
//!   whose boundary term vanishes because the operators leave the wall faces
//!   at zero (T-01.1). So `Ė = −2·r·E` exactly, whatever `L` does, and
//!   `E(t) = E(0)·exp(−2·r·t)`.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use std::f64::consts::PI;

use engine::{
    check_timestep, shallow_water_rhs, Field2D, Grid, OceanState, PhysicalParams, Rk4,
    ShallowWaterRhs, Spacing, Staggering, WaveSpeed, WindStressField, H_STAGGERING, U_STAGGERING,
    V_STAGGERING,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981).
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;

/// Rayleigh damping `r` every time-integration test here runs at, in s⁻¹: an
/// `e`-folding time of about 11.6 days.
///
/// Far stronger than the equatorial Pacific's own damping (order 1/(2 yr)),
/// and deliberately so. A run has to cover several `e`-folding times to say
/// anything about a decay rate, and it has to do so in timesteps the CFL bound
/// of T-01.3 admits; at the physical `r` that is millions of steps, and the
/// truncation error being measured would sit below `f64` roundoff. Nothing in
/// the terms under test depends on the magnitude of `r`, which is exactly what
/// [`damping_adds_exactly_minus_r_times_the_state`] pins down.
const STRONG_DAMPING_PER_S: f64 = 1.0e-6;

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*). Deliberately different from
/// [`BASIN_LY_M`] so an x/y swap cannot pass.
const BASIN_LX_M: f64 = 1.0e7;
/// Meridional extent of the test basin, in metres.
const BASIN_LY_M: f64 = 5.0e6;

/// Amplitude of the test thermocline depth anomaly, in metres. A 20 m
/// departure is the scale of an observed equatorial Pacific anomaly.
const H_AMPLITUDE_M: f64 = 20.0;
/// Amplitude of the test zonal current anomaly, in m/s.
const U_AMPLITUDE_M_PER_S: f64 = 0.3;
/// Amplitude of the test meridional current anomaly, in m/s. Different from
/// [`U_AMPLITUDE_M_PER_S`] so a u/v swap cannot pass.
const V_AMPLITUDE_M_PER_S: f64 = 0.1;

/// Zonal wind stress of the forced test case, in Pa. Easterly trade-wind
/// stress is `τx < 0` (`CONTEXT.md`).
const TRADE_WIND_STRESS_X_PA: f64 = -0.05;
/// Meridional wind stress of the forced test case, in Pa. Different in
/// magnitude from [`TRADE_WIND_STRESS_X_PA`] so an x/y swap cannot pass.
const TRADE_WIND_STRESS_Y_PA: f64 = 0.02;

/// Relative slack allowed where a check is exact in exact arithmetic: a few
/// ulps of `f64` (ε ≈ 2.2e-16) for the handful of operations each point costs.
const ROUNDING_TOLERANCE: f64 = 1.0e-14;

/// The equatorial-Pacific parameter set at a given Rayleigh damping.
fn pacific_params(rayleigh_damping_per_s: f64) -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        rayleigh_damping_per_s,
        engine::EQUATORIAL_BETA_PER_M_PER_S,
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

/// A basin of `n` by `n` cells spanning [`BASIN_LX_M`] by [`BASIN_LY_M`].
/// `dx ≠ dy`, since the basin is not square, so an x/y swap cannot pass.
fn basin(n: usize) -> (Grid, Spacing) {
    let grid = Grid::new(n, n).expect("extents are non-zero");
    let spacing = Spacing::new(BASIN_LX_M / n as f64, BASIN_LY_M / n as f64)
        .expect("a basin spanned by whole cells has positive spacing");
    (grid, spacing)
}

/// Position of a field point, in metres from the basin's southwest corner.
fn position_m(spacing: Spacing, staggering: Staggering, i: usize, j: usize) -> (f64, f64) {
    let (offset_x, offset_y) = staggering.offset_in_cells();
    (
        (i as f64 + offset_x) * spacing.dx_m(),
        (j as f64 + offset_y) * spacing.dy_m(),
    )
}

/// A field at `staggering` holding `g` sampled at each of its points.
fn sample(
    grid: Grid,
    spacing: Spacing,
    staggering: Staggering,
    g: impl Fn(f64, f64) -> f64,
) -> Field2D<f64> {
    let mut field = grid.allocate(staggering, 0.0_f64);
    for j in 0..field.ny() {
        for i in 0..field.nx() {
            let (x_m, y_m) = position_m(spacing, staggering, i, j);
            *field.get_mut(i, j).expect("in-bounds point") = g(x_m, y_m);
        }
    }
    field
}

// --- The term itself: `−r` times the state, at every point. ---

#[test]
fn damping_adds_exactly_minus_r_times_the_state() {
    // Rayleigh damping is linear and pointwise, so the whole of its effect on
    // the right-hand side is the difference between the same evaluation at
    // `r > 0` and at `r = 0`: exactly `−r·h`, `−r·u` and `−r·v`, at every
    // point of every field. Nothing else about the state, the wind or the
    // basin may change with `r`.
    //
    // Every point is checked, the four walls included: damping is a property
    // of the water at a point, not of a difference stencil, so unlike the
    // pressure gradient it has no boundary gap (T-02.3).
    //
    // The `r = 0` evaluation is a baseline rather than an expected value: what
    // it holds is T-02.3's right-hand side, whose values are pinned against
    // analytic derivatives in `shallow_water_rhs.rs`. The analytic claim being
    // made here is about the *difference* of the two, which is `−r·x` on
    // paper — the baseline cancels out of it.
    let (grid, spacing) = basin(SMALL_BASIN_CELLS);
    let damping_per_s = STRONG_DAMPING_PER_S;
    let state = lopsided_state(grid, spacing);
    let wind = WindStressField::uniform(grid, TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);

    let undamped = shallow_water_rhs(&state, pacific_params(0.0), spacing, &wind);
    let damped = shallow_water_rhs(&state, pacific_params(damping_per_s), spacing, &wind);

    for (name, value, undamped_rate, damped_rate) in [
        ("h", state.h(), undamped.h(), damped.h()),
        ("u", state.u(), undamped.u(), damped.u()),
        ("v", state.v(), undamped.v(), damped.v()),
    ] {
        for j in 0..value.ny() {
            for i in 0..value.nx() {
                let anomaly = *value.get(i, j).expect("in-bounds");
                let before = *undamped_rate.get(i, j).expect("in-bounds");
                let after = *damped_rate.get(i, j).expect("in-bounds");
                let expected = before - damping_per_s * anomaly;
                // Relative to the size of the two terms being added rather
                // than to their sum, which can cancel to near zero at a point.
                // Where both terms are exactly zero the slack is zero too, and
                // rightly so: nothing was rounded.
                let slack = ROUNDING_TOLERANCE * (before.abs() + damping_per_s * anomaly.abs());
                assert!(
                    (after - expected).abs() <= slack,
                    "{name} at ({i}, {j}): expected {expected}, got {after}"
                );
            }
        }
    }
}

/// Cells across the small test basin used by the pointwise checks. Odd, so no
/// face falls on the crest of the test bump.
const SMALL_BASIN_CELLS: usize = 9;

/// A state whose three fields are all non-zero and all different: a sinusoid
/// in `h`, and two sinusoids of different amplitude and phase in `u` and `v`.
///
/// The point is only that no field is zero and no two are equal, so that a
/// damping term applied to the wrong variable, or with the wrong sign, cannot
/// pass unnoticed.
fn lopsided_state(grid: Grid, spacing: Spacing) -> OceanState {
    let kx = 2.0 * PI / BASIN_LX_M;
    let ky = 2.0 * PI / BASIN_LY_M;
    let mut state = OceanState::at_rest(grid);
    *state.h_mut() = sample(grid, spacing, H_STAGGERING, |x_m, y_m| {
        H_AMPLITUDE_M * (kx * x_m).sin() * (ky * y_m).sin()
    });
    *state.u_mut() = sample(grid, spacing, U_STAGGERING, |x_m, y_m| {
        U_AMPLITUDE_M_PER_S * (kx * x_m).cos() * (ky * y_m).sin()
    });
    *state.v_mut() = sample(grid, spacing, V_STAGGERING, |x_m, y_m| {
        V_AMPLITUDE_M_PER_S * (kx * x_m).sin() * (ky * y_m).cos()
    });
    state
}

// --- The single-mode test case: a uniform state decays as exp(−r·t). ---

/// Cells across the basin the uniform-mode test runs on. Coarse on purpose:
/// the mode has no spatial structure, so resolution buys nothing, while the
/// wide cells raise the CFL bound of T-01.3 far enough that the timesteps
/// below are admissible for the wave speed the basin also carries.
const COARSE_BASIN_CELLS: usize = 4;

/// Step counts of the uniform-mode convergence check: the same run length in
/// halved timesteps, so each pair gives an independent order estimate.
const STEP_COUNTS: [usize; 3] = [8, 16, 32];

/// How far a measured RK4 order may sit from the theoretical 4.
///
/// For `ẏ = λ·y` the RK4 error per step is `−(λ·dt)⁵/120` with a relative
/// correction of `−(5/6)·λ·dt` at the next order, so a run of fixed length
/// gives an error `∝ dt⁴·(1 + (5/6)·r·dt)`. At the coarsest pair here
/// (`r·dt = 1/8` and `1/16`) that shifts the measured ratio from 16 to about
/// 16.8, i.e. the order from 4 to 4.07. A tolerance of 0.15 absorbs that with
/// room to spare while still excluding third- or fifth-order behaviour.
const ORDER_TOLERANCE: f64 = 0.15;

#[test]
fn a_uniform_state_decays_as_exp_minus_r_t() {
    // The single-mode case of the acceptance criteria. A state that is the
    // same at every point has no pressure gradient (`∂h/∂x = ∂h/∂y = 0`) and
    // no divergence (`∂u/∂x = ∂v/∂y = 0`), so with `τ = 0` the right-hand side
    // collapses to `ẏ = −r·y` for each of `h`, `u` and `v` separately. The
    // analytic solution is `y(0)·exp(−r·t)`, and the error of RK4 against it
    // must fall like `dt⁴` (ADR-0003).
    let (grid, spacing) = basin(COARSE_BASIN_CELLS);
    let params = pacific_params(STRONG_DAMPING_PER_S);
    // One `e`-folding time, which is where `exp(−r·t)` says the most.
    let run_length_s = 1.0 / STRONG_DAMPING_PER_S;

    let errors: Vec<f64> = STEP_COUNTS
        .iter()
        .map(|&steps| {
            let dt_s = run_length_s / steps as f64;
            assert_admissible_timestep(dt_s, spacing, params);

            let mut state = OceanState::at_rest(grid);
            state.h_mut().as_mut_slice().fill(H_AMPLITUDE_M);
            state.u_mut().as_mut_slice().fill(U_AMPLITUDE_M_PER_S);
            state.v_mut().as_mut_slice().fill(V_AMPLITUDE_M_PER_S);

            let mut worst_relative_error = 0.0_f64;
            let mut previous = state.clone();
            for_each_step(
                grid,
                spacing,
                params,
                &mut state,
                dt_s,
                steps,
                |step, now| {
                    let elapsed_s = step as f64 * dt_s;
                    let expected_factor = (-STRONG_DAMPING_PER_S * elapsed_s).exp();
                    for (amplitude, field) in [
                        (H_AMPLITUDE_M, now.h()),
                        (U_AMPLITUDE_M_PER_S, now.u()),
                        (V_AMPLITUDE_M_PER_S, now.v()),
                    ] {
                        let expected = amplitude * expected_factor;
                        for value in field.as_slice() {
                            worst_relative_error =
                                worst_relative_error.max((value - expected).abs() / expected);
                        }
                    }
                    // Monotone decay towards rest, the other half of the
                    // acceptance criterion: every anomaly shrinks in magnitude at
                    // every step, and none of them changes sign.
                    for (before, after) in [
                        (previous.h(), now.h()),
                        (previous.u(), now.u()),
                        (previous.v(), now.v()),
                    ] {
                        for (was, is) in before.as_slice().iter().zip(after.as_slice()) {
                            assert!(
                                is.abs() < was.abs() && is.signum() == was.signum(),
                                "a damped anomaly went from {was} to {is}"
                            );
                        }
                    }
                    previous = now.clone();
                },
            );
            worst_relative_error
        })
        .collect();

    for window in errors.windows(2) {
        let order = (window[0] / window[1]).log2();
        assert!(
            (order - 4.0).abs() < ORDER_TOLERANCE,
            "expected fourth-order convergence on exp(-r·t), measured order {order} from errors {errors:?}"
        );
    }
}

/// Panic unless `dt_s` is a timestep this basin's CFL bound admits.
///
/// The decay tests are not about stability, but a run on a timestep the engine
/// would refuse (T-01.3) would not be a run of this model.
fn assert_admissible_timestep(dt_s: f64, spacing: Spacing, params: PhysicalParams) {
    let wave_speed =
        WaveSpeed::new(params.kelvin_wave_speed_m_per_s()).expect("a positive wave speed");
    check_timestep(dt_s, spacing, wave_speed).unwrap_or_else(|error| {
        panic!("the test's own timestep must be admissible: {error}");
    });
}

// --- Energy: an unforced perturbation decays monotonically to rest. ---

/// Cells across the basin the energy tests run on. Fine enough that the
/// gravest zonal mode below is smooth on the grid, coarse enough to stay fast.
const WAVE_BASIN_CELLS: usize = 16;
/// Steps of an energy run.
const ENERGY_RUN_STEPS: usize = 256;
/// Length of an energy run, in `e`-folding times of the current anomaly. Five
/// leaves `exp(−2·r·t) = 4.5e-5` of the initial energy, which is "the rest
/// state" to any resolution this test can see.
const ENERGY_RUN_EFOLDINGS: f64 = 5.0;

/// Relative tolerance on the energy decay `E(t) = E(0)·exp(−2·r·t)`.
///
/// The identity is exact for the semi-discrete system, so the only error is
/// RK4's. The initial condition below is a single discrete normal mode, whose
/// eigenvalues are `−r ± i·ω` with `ω = c·(2/dx)·sin(k·dx/2) ≈ 8.6e-7 s⁻¹`;
/// at this run's `dt ≈ 2.0e4 s` that is `|λ·dt| ≈ 0.026`, and RK4's global
/// relative error over `N = 256` steps is about `N·|λ·dt|⁵/120 ≈ 2e-8`, so
/// twice that in the energy. A bound of 1e-6 sits an order of magnitude above
/// it and eight below the effect being measured.
const ENERGY_TOLERANCE: f64 = 1.0e-6;

#[test]
fn an_unforced_perturbation_decays_monotonically_to_the_rest_state() {
    // The acceptance criterion, in the norm the question is well posed in: an
    // individual anomaly oscillates as the perturbation radiates away as
    // waves, but the energy of the whole basin can only fall, because the
    // wave terms merely move energy around and the damping removes it at
    // `Ė = −2·r·E`.
    let (grid, spacing) = basin(WAVE_BASIN_CELLS);
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let (dt_s, steps) = energy_run_schedule(spacing, params);
    let initial = gravest_zonal_mode(grid, spacing);

    let energies = energy_history(grid, spacing, params, initial, dt_s, steps);

    for (step, pair) in energies.windows(2).enumerate() {
        assert!(
            pair[1] < pair[0],
            "energy rose from {} to {} at step {step}",
            pair[0],
            pair[1]
        );
    }

    // …and it ends at rest. The rate it got there at is
    // `perturbation_energy_decays_as_exp_minus_2_r_t`'s business; what this
    // test adds is that the destination really is the rest state.
    let final_ratio = energies.last().expect("a non-empty history") / energies[0];
    assert!(
        final_ratio < AT_REST_ENERGY_FRACTION,
        "after {ENERGY_RUN_EFOLDINGS} e-folding times the basin should be at rest, but {final_ratio} of its energy is left"
    );
}

/// Fraction of its initial energy a basin may still hold and count as being
/// back at the rest state.
///
/// Not a tolerance on a measurement but a reading of "at rest": after
/// [`ENERGY_RUN_EFOLDINGS`] `e`-folding times the analytic energy ratio is
/// `exp(−2·5) = 4.5e-5`, so 1e-4 is the next round number above it. A run that
/// damped at even half the required rate would leave `exp(−5) = 6.7e-3` and
/// fail.
const AT_REST_ENERGY_FRACTION: f64 = 1.0e-4;

#[test]
fn perturbation_energy_decays_as_exp_minus_2_r_t() {
    // The rate, checked the whole way down rather than only at the end:
    // damping `h`, `u` and `v` at the same `r` makes `Ė = −2·r·E` an identity
    // of the discrete system, independent of what the wave terms do.
    let (grid, spacing) = basin(WAVE_BASIN_CELLS);
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let (dt_s, steps) = energy_run_schedule(spacing, params);
    let initial = gravest_zonal_mode(grid, spacing);

    let energies = energy_history(grid, spacing, params, initial, dt_s, steps);

    for (step, energy) in energies.iter().enumerate() {
        let elapsed_s = step as f64 * dt_s;
        let expected = energies[0] * (-2.0 * STRONG_DAMPING_PER_S * elapsed_s).exp();
        assert!(
            (energy - expected).abs() <= ENERGY_TOLERANCE * expected,
            "at step {step} ({elapsed_s} s): expected energy {expected}, got {energy}"
        );
    }
}

#[test]
fn an_undamped_perturbation_keeps_its_energy() {
    // The control on the two tests above: at `r = 0` the same run must
    // conserve energy to the same tolerance. Without this, a damping term
    // applied unconditionally — ignoring `r` — would still pass a decay test
    // as long as the rate happened to match.
    let (grid, spacing) = basin(WAVE_BASIN_CELLS);
    let params = pacific_params(0.0);
    let (dt_s, steps) = energy_run_schedule(spacing, params);
    let initial = gravest_zonal_mode(grid, spacing);

    let energies = energy_history(grid, spacing, params, initial, dt_s, steps);

    for (step, energy) in energies.iter().enumerate() {
        assert!(
            (energy - energies[0]).abs() <= ENERGY_TOLERANCE * energies[0],
            "at step {step}: energy drifted from {} to {energy}",
            energies[0]
        );
    }
}

/// Timestep and step count of an energy run: [`ENERGY_RUN_EFOLDINGS`]
/// `e`-folding times of [`STRONG_DAMPING_PER_S`] in [`ENERGY_RUN_STEPS`]
/// steps, checked against the CFL bound.
fn energy_run_schedule(spacing: Spacing, params: PhysicalParams) -> (f64, usize) {
    let run_length_s = ENERGY_RUN_EFOLDINGS / STRONG_DAMPING_PER_S;
    let dt_s = run_length_s / ENERGY_RUN_STEPS as f64;
    assert_admissible_timestep(dt_s, spacing, params);
    (dt_s, ENERGY_RUN_STEPS)
}

/// The gravest zonal standing mode of the closed basin: `h = A·cos(π·x/Lx)`,
/// at rest, uniform in `y`.
///
/// Two properties earn it its place as the initial condition of the energy
/// tests. It is an exact eigenvector of the discrete wave operator — the
/// C-grid difference of `cos(k·x)` at cell centers is `sin(k·x)` at faces,
/// which vanishes on both walls — so the run excites one frequency, low enough
/// that RK4's time-truncation error stays far below the decay being measured.
/// And its velocities start and stay exactly zero on the four walls, which is
/// the condition under which the discrete energy's boundary term vanishes.
fn gravest_zonal_mode(grid: Grid, spacing: Spacing) -> OceanState {
    let wavenumber_per_m = PI / BASIN_LX_M;
    let mut state = OceanState::at_rest(grid);
    *state.h_mut() = sample(grid, spacing, H_STAGGERING, |x_m, _y_m| {
        H_AMPLITUDE_M * (wavenumber_per_m * x_m).cos()
    });
    state
}

/// The discrete energy `E = (g'/2)·Σh² + (H/2)·Σ(u² + v²)`, summed over grid
/// points, in m³/s² — energy per unit reference density and per unit cell
/// area, which is constant and so drops out of every ratio taken here.
///
/// The weights are the ones that make the pressure-gradient and continuity
/// pair skew: `g'·h·ḣ + H·u·u̇` telescopes to a boundary term under summation
/// by parts, and the boundary term is zero because the wall faces carry no
/// velocity. Hence `Ė = −2·r·E`.
fn wave_energy(state: &OceanState, params: PhysicalParams) -> f64 {
    let sum_of_squares =
        |field: &Field2D<f64>| -> f64 { field.as_slice().iter().map(|value| value * value).sum() };
    0.5 * params.reduced_gravity_m_per_s2() * sum_of_squares(state.h())
        + 0.5
            * params.mean_thermocline_depth_m()
            * (sum_of_squares(state.u()) + sum_of_squares(state.v()))
}

/// [`wave_energy`] after each of `steps` unforced steps, starting with the
/// initial state's own energy.
fn energy_history(
    grid: Grid,
    spacing: Spacing,
    params: PhysicalParams,
    mut state: OceanState,
    dt_s: f64,
    steps: usize,
) -> Vec<f64> {
    let mut energies = vec![wave_energy(&state, params)];
    for_each_step(grid, spacing, params, &mut state, dt_s, steps, |_, now| {
        energies.push(wave_energy(now, params));
    });
    energies
}

/// Advance `state` through `steps` unforced RK4 steps, calling `observe` with
/// the step number (1-based) and the state after each one.
fn for_each_step(
    grid: Grid,
    spacing: Spacing,
    params: PhysicalParams,
    state: &mut OceanState,
    dt_s: f64,
    steps: usize,
    mut observe: impl FnMut(usize, &OceanState),
) {
    let mut evaluator = ShallowWaterRhs::new(grid, spacing, params);
    let calm = WindStressField::calm(grid);
    let mut integrator = Rk4::new(state);
    for step in 0..steps {
        let t_s = step as f64 * dt_s;
        integrator.step(state, t_s, dt_s, &mut |now: &OceanState, _t, out| {
            evaluator.evaluate(now, &calm, out);
        });
        observe(step + 1, state);
    }
}
