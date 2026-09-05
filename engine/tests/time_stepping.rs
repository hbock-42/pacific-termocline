//! Acceptance tests for T-02.5 — the whole shallow-water right-hand side wired
//! through the RK4 stepper of T-01.3, run end to end on a test grid.
//!
//! The step under test is the complete one Epic 02 set out to build: the
//! pressure gradient and continuity terms of T-02.3, the surface stress stub,
//! the Rayleigh damping of T-02.4 *and* the beta-plane Coriolis term of T-02.2,
//! all evaluated together at each of RK4's four stages.
//!
//! ```text
//! ∂u/∂t = +f·v − g'·∂h/∂x + τx/(ρ₀·H) − r·u
//! ∂v/∂t = −f·u − g'·∂h/∂y + τy/(ρ₀·H) − r·v
//! ∂h/∂t =       −H·(∂u/∂x + ∂v/∂y)    − r·h
//! ```
//!
//! Every expected value below comes from an independent source — an analytic
//! solution of the equations above, the RK4 amplification polynomial quoted in
//! `termocline-numerics`' CFL derivation, or the beta-plane definition `f = β·y`
//! written out from `CONTEXT.md` — never from running this code.
//!
//! # The energy budget of the discrete system
//!
//! The acceptance criterion asks an undamped, unforced run to conserve the
//! total energy `E = (g'/2)·Σh² + (H/2)·Σ(u² + v²)` "to within
//! numerical-diffusion-scale tolerance". Two facts about this discretisation
//! say what that scale is, and both are properties of the scheme rather than
//! measurements of it:
//!
//! - The pressure-gradient/continuity pair is **exactly** skew in that energy.
//!   Summation by parts on the C-grid leaves a boundary term `h·u` at the
//!   walls, and the operators of T-01.1 leave every wall face at zero, so a
//!   basin that starts with quiescent walls and is never forced keeps them
//!   quiescent and the boundary term vanishes for all time. This pair
//!   therefore contributes nothing to `dE/dt`, at any resolution.
//! - The Coriolis pair is skew only to `O(Δy²)`. Writing the four-point
//!   averages out, `Σ_u u·f̄·v̄` and `Σ_v v·f·ū` pair up term by term over the
//!   same `(u, v)` neighbour pairs, so they would cancel exactly if the two
//!   sums weighted each pair with the same `f`. They do not: the u equation
//!   evaluates `f` on the cell-center rows and the v equation on the
//!   north/south-face rows, half a cell apart, so each pair survives with a
//!   residual weight `±β·Δy/2` and
//!
//!   ```text
//!   dE/dt = −H·(β·Δy²/4)·Σ_u u·∂v/∂y + O(Δy⁴)
//!   ```
//!
//!   a truncation-scale term that vanishes at second order as the grid is
//!   refined. It is the discretisation T-02.2 delivered and the one its
//!   hand-computed example pins; making it exactly skew would need `f` at the
//!   cell corners, which is a change to that ticket's contract, not this one's.
//!
//! So the energy is conserved to truncation scale, and the honest test of that
//! is a convergence test rather than a fixed threshold
//! (CODING_STANDARDS.md § Tests): the drift must *shrink at the scheme's
//! order* as the basin is resolved. RK4's own contribution shrinks faster
//! still — the timestep falls with the cell width under the CFL bound, so its
//! share of the drift over a fixed run length falls like `dt³` — which leaves
//! second order as the rate to hold the scheme to.

use std::cell::RefCell;

use engine::{
    max_stable_dt, step, BetaPlane, CflError, Field2D, Grid, OceanState, PhysicalParams, Solver,
    SolverError, Spacing, Staggering, WaveSpeed, WindStressField, H_STAGGERING, U_STAGGERING,
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

/// Rayleigh damping `r` the decay checks run at, in s⁻¹: an `e`-folding time of
/// about 11.6 days. Far stronger than the equatorial Pacific's own damping, for
/// the reason spelled out in `rayleigh_damping.rs` — a decay has to be visible
/// inside a run of CFL-admissible steps.
const STRONG_DAMPING_PER_S: f64 = 1.0e-6;

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*). Deliberately different from
/// [`BASIN_LY_M`] so an x/y swap cannot pass.
const BASIN_LX_M: f64 = 1.0e7;
/// Meridional extent of the test basin, in metres: an equatorial channel
/// reaching ±500 km, or about 1.4 equatorial deformation radii
/// (`Le = √(c/β) = 3.45×10⁵ m` for the parameters above), either side of the
/// equator.
///
/// Narrower than the ±2500 km basin of `CONTEXT.md` on purpose, and the reason
/// is the second timestep bound `Solver` enforces. `|f| = β·|y|` grows away
/// from the equator, and RK4 can only follow the rotation while `|f|·dt` stays
/// inside its stability region; at ±2500 km on a grid coarse enough to test
/// quickly, the gravity-wave CFL bound would admit a timestep twice the
/// inertial period at the walls, which `Solver::new` rightly refuses (see
/// `a_timestep_that_cannot_resolve_the_basins_rotation_is_refused`). A channel
/// a few deformation radii wide is where the equatorial wave physics lives, is
/// what "small idealized grid" means in this ticket, and leaves the
/// gravity-wave bound the binding one — which is the configuration the
/// acceptance criteria are about.
const BASIN_LY_M: f64 = 1.0e6;

/// Amplitude of the test thermocline depth anomaly, in metres. A 20 m
/// departure is the scale of an observed equatorial Pacific anomaly.
const H_AMPLITUDE_M: f64 = 20.0;
/// Amplitude of the test meridional current anomaly, in m/s.
const V_AMPLITUDE_M_PER_S: f64 = 0.1;

/// Zonal wind stress of the forced test cases, in Pa. Easterly trade-wind
/// stress is `τx < 0` (`CONTEXT.md`).
const TRADE_WIND_STRESS_X_PA: f64 = -0.05;
/// Meridional wind stress of the forced test cases, in Pa. Different in
/// magnitude from [`TRADE_WIND_STRESS_X_PA`] so an x/y swap cannot pass.
const TRADE_WIND_STRESS_Y_PA: f64 = 0.02;

/// Relative slack allowed where a check is exact in exact arithmetic: a few
/// tens of ulps of `f64` (ε ≈ 2.2e-16) for the handful of operations per point
/// a step costs.
const ROUNDING_TOLERANCE: f64 = 1.0e-14;

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

/// A basin of `nx` by `ny` cells spanning [`BASIN_LX_M`] by [`BASIN_LY_M`].
/// `dx ≠ dy`, since the basin is not square, so an x/y swap cannot pass.
fn basin(nx: usize, ny: usize) -> (Grid, Spacing) {
    let grid = Grid::new(nx, ny).expect("extents are non-zero");
    let spacing = Spacing::new(BASIN_LX_M / nx as f64, BASIN_LY_M / ny as f64)
        .expect("a basin spanned by whole cells has positive spacing");
    (grid, spacing)
}

/// The beta-plane of a test basin: the equator through its middle.
fn equatorial_plane(params: PhysicalParams, spacing: Spacing, grid: Grid) -> BetaPlane {
    BetaPlane::centered_on_equator(params, spacing, grid)
}

/// The CFL-stable maximum timestep for this basin and parameter set, in
/// seconds — the "CFL-safe `dt` from Epic 01" the acceptance criteria name.
fn cfl_safe_dt_s(spacing: Spacing, params: PhysicalParams) -> f64 {
    let wave_speed =
        WaveSpeed::new(params.kelvin_wave_speed_m_per_s()).expect("a positive wave speed");
    max_stable_dt(spacing, wave_speed)
}

/// A solver for this basin, or a panic naming the timestep it refused.
fn solver_for(grid: Grid, spacing: Spacing, params: PhysicalParams, dt_s: f64) -> Solver {
    Solver::new(
        grid,
        spacing,
        params,
        equatorial_plane(params, spacing, grid),
        dt_s,
    )
    .unwrap_or_else(|error| panic!("the test's own timestep must be admissible: {error}"))
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

// --- A single step, against the RK4 amplification polynomial. ---

/// One step of `ẏ = λ·y` under classic RK4 multiplies `y` by
/// `R(z) = 1 + z + z²/2 + z³/6 + z⁴/24` with `z = λ·dt`: the degree-4
/// truncation of `exp(z)`, and the standard amplification factor of the method
/// (Hairer & Wanner, *Solving Ordinary Differential Equations I*, § II.2 — the
/// same result `termocline-numerics`' CFL bound is derived from).
fn rk4_amplification(z: f64) -> f64 {
    1.0 + z + z * z / 2.0 + z * z * z / 6.0 + z * z * z * z / 24.0
}

#[test]
fn one_step_of_a_uniform_thermocline_anomaly_is_the_rk4_amplification_polynomial() {
    // The single-timestep half of the deliverable, checked against a closed
    // form. A basin whose `h` is the same everywhere and whose currents are at
    // rest has no pressure gradient, no divergence and no Coriolis coupling
    // (both products carry a velocity, and both velocities are zero), so with
    // `τ = 0` the whole right-hand side collapses to `ḣ = −r·h`, `u̇ = 0`,
    // `v̇ = 0`. RK4 turns that into a multiplication by `R(−r·dt)` exactly, and
    // leaves the currents untouched.
    let (grid, spacing) = basin(UNIFORM_BASIN_CELLS, UNIFORM_BASIN_CELLS);
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let dt_s = cfl_safe_dt_s(spacing, params);

    let mut state = OceanState::at_rest(grid);
    state.h_mut().as_mut_slice().fill(H_AMPLITUDE_M);
    let calm = WindStressField::calm(grid);

    let advanced = step(
        &state,
        dt_s,
        params,
        spacing,
        equatorial_plane(params, spacing, grid),
        |_t_s| &calm,
    )
    .expect("the CFL-safe timestep is admissible");

    let expected_m = H_AMPLITUDE_M * rk4_amplification(-params.rayleigh_damping_per_s() * dt_s);
    for value in advanced.h().as_slice() {
        assert!(
            (value - expected_m).abs() <= ROUNDING_TOLERANCE * expected_m.abs(),
            "expected h = {expected_m} m after one step, got {value} m"
        );
    }
    // Nothing set the currents in motion, so they are still exactly at rest —
    // not nearly at rest.
    assert!(advanced.u().as_slice().iter().all(|value| *value == 0.0));
    assert!(advanced.v().as_slice().iter().all(|value| *value == 0.0));
}

/// Cells across the basin the uniform-state checks run on. Coarse on purpose:
/// the state has no spatial structure, so resolution buys nothing.
const UNIFORM_BASIN_CELLS: usize = 4;

/// Timestep of the pressure-gradient check, in seconds; see
/// [`PRESSURE_GRADIENT_TOLERANCE`].
const PRESSURE_GRADIENT_STEP_S: f64 = 100.0;

/// Relative tolerance on the leading-order acceleration `u(dt) = −g'·∂h/∂x·dt`.
///
/// From a basin at rest in `u` and `v`, the next term is `(dt²/2)·ü` with
/// `ü = −g'·∂ḣ/∂x = g'·H·∂²(∂u/∂x + ∂v/∂y)/∂x²`, i.e. a relative correction of
/// order `(c·k·dt)²/2` for the mode below. Its wavenumber is `k = π/Lx`, so
/// `c·k·dt = 2.74 × 3.14×10⁻⁷ × 100 = 8.6×10⁻⁵` and the correction is
/// `3.7×10⁻⁹`. A bound of 1e-6 clears that by two orders of magnitude and
/// still fails by six if the term is missing or mis-signed.
const PRESSURE_GRADIENT_TOLERANCE: f64 = 1.0e-6;

#[test]
fn a_step_accelerates_a_thermocline_slope_down_its_own_gradient() {
    // The other end-to-end single step: a basin whose only anomaly is a slope
    // in `h`. All three of the momentum equation's other terms are zero at
    // rest, so the whole of the initial acceleration is `−g'·∂h/∂x`, and the
    // C-grid difference of `A·cos(k·x)` at cell centers is exactly
    // `−A·(2/dx)·sin(k·dx/2)·sin(k·x_face)` at the face between them — the
    // discrete derivative of the mode, written out from the stencil rather
    // than measured.
    let (grid, spacing) = basin(CORIOLIS_BASIN_CELLS, CORIOLIS_BASIN_CELLS);
    let params = pacific_params(0.0);
    let mut state = gravest_zonal_mode(grid, spacing);

    let mut solver = solver_for(grid, spacing, params, PRESSURE_GRADIENT_STEP_S);
    let calm = WindStressField::calm(grid);
    solver.step(&mut state, 0.0, |_t_s| &calm);

    let wavenumber_per_m = std::f64::consts::PI / BASIN_LX_M;
    let dx_m = spacing.dx_m();
    let slope_amplitude_per_m = H_AMPLITUDE_M * 2.0 / dx_m * (0.5 * wavenumber_per_m * dx_m).sin();
    for j in 0..state.u().ny() {
        // The two wall faces have a cell on one side only, so the operators
        // leave the gradient at zero there (T-01.1) and the wall stays at rest.
        for i in 1..state.u().nx() - 1 {
            let x_m = i as f64 * dx_m;
            let expected_m_per_s = params.reduced_gravity_m_per_s2()
                * slope_amplitude_per_m
                * (wavenumber_per_m * x_m).sin()
                * PRESSURE_GRADIENT_STEP_S;
            assert!(
                (state.u().get(i, j).expect("in-bounds point") - expected_m_per_s).abs()
                    <= PRESSURE_GRADIENT_TOLERANCE * expected_m_per_s.abs(),
                "u at ({i}, {j}): expected {expected_m_per_s} m/s, got {} m/s",
                state.u().get(i, j).expect("in-bounds point")
            );
        }
    }
    // The walls really are still at rest, exactly.
    for j in 0..state.u().ny() {
        assert_eq!(state.u().get(0, j), Some(&0.0));
        assert_eq!(state.u().get(state.u().nx() - 1, j), Some(&0.0));
    }
}

// --- The Coriolis term reaches the step. ---

/// Cells across the basin the Coriolis check runs on. `ny` is even, so no `u`
/// row lies on the equator and every expected value below is non-zero, which
/// is what lets the check be a relative one.
const CORIOLIS_BASIN_CELLS: usize = 8;

/// Timestep of the Coriolis check, in seconds. Short compared with the
/// rotation period `1/f ≳ 1.7×10⁴ s` so that the leading-order rotation is
/// what the check sees; see [`CORIOLIS_TOLERANCE`].
const CORIOLIS_STEP_S: f64 = 100.0;

/// Relative tolerance on the leading-order rotation `u(dt) = f·v·dt`.
///
/// From a basin at rest in `u` and `h`, the next term of the Taylor expansion
/// of `u` is third order: `ü(0) = f·v̇ = f·(−f·ū) = 0` because `u` starts at
/// zero, so the error is `(dt³/6)·u'''`, i.e. a relative `(f·dt)²/6`. The
/// largest `|f|` in this basin is `β·Ly/2 = 1.15×10⁻⁵ s⁻¹`, so at
/// [`CORIOLIS_STEP_S`] that is `(1.15×10⁻³)²/6 ≈ 2.2×10⁻⁷`. A bound of 1e-4
/// clears it by nearly three orders of magnitude and still fails by a factor
/// of 10⁴ if the term is missing, mis-signed or read at the wrong row.
const CORIOLIS_TOLERANCE: f64 = 1.0e-4;

#[test]
fn a_step_rotates_a_meridional_current_into_a_zonal_one() {
    // The Coriolis term of T-02.2 has to be part of the step, not just part of
    // the crate. A basin with a uniform northward current and nothing else has
    // `u̇ = +f·v̄ = f·v` at every interior face — `f = β·y` read on the
    // cell-center rows, written out here from `CONTEXT.md` rather than asked of
    // the code under test.
    let (grid, spacing) = basin(CORIOLIS_BASIN_CELLS, CORIOLIS_BASIN_CELLS);
    let params = pacific_params(0.0);
    let mut state = OceanState::at_rest(grid);
    state.v_mut().as_mut_slice().fill(V_AMPLITUDE_M_PER_S);

    let mut solver = solver_for(grid, spacing, params, CORIOLIS_STEP_S);
    let calm = WindStressField::calm(grid);
    solver.step(&mut state, 0.0, |_t_s| &calm);

    for j in 0..state.u().ny() {
        // Measured from the equator, which the basin straddles symmetrically,
        // rather than from its southwest corner. Where the `u` row sits within
        // its cell is the grid's business, not this test's
        // (CODING_STANDARDS.md § Scope guards).
        let (_, y_from_corner_m) = position_m(spacing, U_STAGGERING, 0, j);
        let y_m = y_from_corner_m - 0.5 * BASIN_LY_M;
        // The four-point average that carries `v` onto a `u` face takes the
        // two `v` rows flanking it. On the outermost `u` rows one of those two
        // is the coast, which the no-normal-flow condition of T-04.2 holds at
        // rest, so `v̄` there is half the uniform current — an exact factor,
        // not an approximation, and one worth pinning rather than skipping.
        let on_the_coast = j == 0 || j == state.u().ny() - 1;
        let interpolated_v_m_per_s = if on_the_coast {
            0.5 * V_AMPLITUDE_M_PER_S
        } else {
            V_AMPLITUDE_M_PER_S
        };
        let expected_m_per_s = BETA_PER_M_PER_S * y_m * interpolated_v_m_per_s * CORIOLIS_STEP_S;
        // The two wall faces have a cell on one side only, so they carry no
        // interpolated `v` and stay at rest (T-01.1, and now T-04.2); the
        // interior is what this check is about.
        for i in 1..state.u().nx() - 1 {
            let value = *state.u().get(i, j).expect("in-bounds point");
            assert!(
                (value - expected_m_per_s).abs() <= CORIOLIS_TOLERANCE * expected_m_per_s.abs(),
                "u at ({i}, {j}): expected {expected_m_per_s} m/s, got {value} m/s"
            );
        }
    }
}

// --- The wind-stress function reaches the step, at the stage times. ---

/// Timestep of the wind-forcing check, in seconds; see [`WIND_TOLERANCE`].
const WIND_STEP_S: f64 = 10.0;

/// Relative tolerance on the leading-order acceleration `u(dt) = τx·dt/(ρ₀·H)`.
///
/// From rest the next term is `(dt²/2)·ü` with `ü(0) = f·v̇ = f·τy/(ρ₀·H)`, a
/// relative correction of `(|f|·dt/2)·(τy/τx) = 2.3×10⁻⁵` at the basin's
/// largest `|f| = β·Ly/2 = 1.15×10⁻⁵ s⁻¹` and [`WIND_STEP_S`]. A bound of 1e-3
/// clears that and still fails by three orders of magnitude if the stress
/// never reaches the step.
const WIND_TOLERANCE: f64 = 1.0e-3;

#[test]
fn a_constant_wind_stress_accelerates_a_basin_at_rest() {
    // The stubbed constant forcing of the deliverable. A basin at rest has no
    // gradient, no divergence and no rotation to feel, so the whole of its
    // initial acceleration is `τ/(ρ₀·H)` — the surface stress divided by the
    // mass of the upper layer per unit area.
    let (grid, spacing) = basin(CORIOLIS_BASIN_CELLS, CORIOLIS_BASIN_CELLS);
    let params = pacific_params(0.0);
    let mut state = OceanState::at_rest(grid);

    let mut solver = solver_for(grid, spacing, params, WIND_STEP_S);
    let trade_winds =
        WindStressField::uniform(grid, TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);
    solver.step(&mut state, 0.0, |_t_s| &trade_winds);

    // Away from the coast, that is. The stress is applied at every face
    // including the walls, and the no-normal-flow condition of T-04.2 discards
    // the wall faces' share of it: the normal velocity at a coast is not a
    // degree of freedom, so it stays at exactly rest however hard the wind
    // blows there. Both halves are checked.
    let layer_mass_kg_per_m2 = REFERENCE_DENSITY_KG_PER_M3 * PACIFIC_MEAN_DEPTH_M;
    for (name, stress_pa, field, walls_are_columns) in [
        ("u", TRADE_WIND_STRESS_X_PA, state.u(), true),
        ("v", TRADE_WIND_STRESS_Y_PA, state.v(), false),
    ] {
        let expected_m_per_s = stress_pa / layer_mass_kg_per_m2 * WIND_STEP_S;
        for j in 0..field.ny() {
            for i in 0..field.nx() {
                let value = *field.get(i, j).expect("an in-bounds face");
                // `u`'s walls are its first and last columns, `v`'s its first
                // and last rows: each field's coast is the line it is
                // staggered across.
                let on_the_coast = if walls_are_columns {
                    i == 0 || i == field.nx() - 1
                } else {
                    j == 0 || j == field.ny() - 1
                };
                if on_the_coast {
                    assert_eq!(
                        value, 0.0,
                        "{name} at the wall face ({i}, {j}) must stay at rest"
                    );
                } else {
                    assert!(
                        (value - expected_m_per_s).abs() <= WIND_TOLERANCE * expected_m_per_s.abs(),
                        "{name}: expected {expected_m_per_s} m/s after one step at ({i}, {j}), \
                         got {value} m/s"
                    );
                }
            }
        }
    }
    // Easterly trade winds push the water westward, so `u < 0` — the sign
    // convention of `CONTEXT.md`, which a stress applied with the wrong sign
    // would invert.
    for j in 0..state.u().ny() {
        for i in 1..state.u().nx() - 1 {
            assert!(*state.u().get(i, j).expect("an interior face") < 0.0);
        }
    }
}

#[test]
fn the_wind_stress_function_is_sampled_at_the_four_rk4_stage_times() {
    // The forcing is a function of time, not a constant baked into the solver:
    // Epic 03's scenarios vary with `t`, and a step that sampled the stress
    // once would integrate the wrong forcing. The stage times are the nodes
    // `c = [0, 1/2, 1/2, 1]` of the classic RK4 tableau (ADR-0003), offset
    // from the step's start time.
    let (grid, spacing) = basin(UNIFORM_BASIN_CELLS, UNIFORM_BASIN_CELLS);
    let params = pacific_params(0.0);
    let dt_s = cfl_safe_dt_s(spacing, params);
    let start_s = 7.0 * dt_s;

    let mut state = OceanState::at_rest(grid);
    let mut solver = solver_for(grid, spacing, params, dt_s);
    let calm = WindStressField::calm(grid);
    let sampled_at_s = RefCell::new(Vec::new());
    solver.step(&mut state, start_s, |t_s| {
        sampled_at_s.borrow_mut().push(t_s);
        &calm
    });

    assert_eq!(
        sampled_at_s.into_inner(),
        vec![
            start_s,
            start_s + 0.5 * dt_s,
            start_s + 0.5 * dt_s,
            start_s + dt_s
        ]
    );
}

// --- Acceptance criterion: an undamped, unforced run conserves energy. ---

/// The discrete energy `E = (g'/2)·Σh² + (H/2)·Σ(u² + v²)`, summed over grid
/// points, in m³/s² — energy per unit reference density and per unit cell
/// area, which is constant across a run and so drops out of every ratio taken
/// here.
///
/// These are the weights that make the pressure-gradient and continuity pair
/// skew under summation by parts on the C-grid; see the module comment.
fn wave_energy(state: &OceanState, params: PhysicalParams) -> f64 {
    let sum_of_squares =
        |field: &Field2D<f64>| -> f64 { field.as_slice().iter().map(|value| value * value).sum() };
    0.5 * params.reduced_gravity_m_per_s2() * sum_of_squares(state.h())
        + 0.5
            * params.mean_thermocline_depth_m()
            * (sum_of_squares(state.u()) + sum_of_squares(state.v()))
}

/// The gravest zonal standing mode of the closed basin: `h = A·cos(π·x/Lx)`,
/// at rest, uniform in `y`.
///
/// It is an exact eigenvector of the discrete wave operator — the C-grid
/// difference of `cos(k·x)` at cell centers is `sin(k·x)` at faces, which
/// vanishes on both walls — so a run excites one low frequency, and its
/// velocities start and stay exactly zero on the four walls, which is the
/// condition under which the discrete energy's boundary term vanishes.
fn gravest_zonal_mode(grid: Grid, spacing: Spacing) -> OceanState {
    let wavenumber_per_m = std::f64::consts::PI / BASIN_LX_M;
    let mut state = OceanState::at_rest(grid);
    *state.h_mut() = sample(grid, spacing, H_STAGGERING, |x_m, _y_m| {
        H_AMPLITUDE_M * (wavenumber_per_m * x_m).cos()
    });
    state
}

/// Length of the energy runs, in crossings of the basin at the Kelvin wave
/// speed `c = √(g'·H)`: long enough that the perturbation has radiated,
/// reflected off both walls and interfered with itself several times over
/// rather than merely drifted.
const ENERGY_RUN_CROSSINGS: f64 = 4.0;

/// Length of an energy run, in seconds. `c = √(g'·H)` is written out from the
/// definition in `CONTEXT.md` rather than asked of the code under test.
fn energy_run_length_s() -> f64 {
    let wave_speed_m_per_s = (PACIFIC_REDUCED_GRAVITY_M_PER_S2 * PACIFIC_MEAN_DEPTH_M).sqrt();
    ENERGY_RUN_CROSSINGS * BASIN_LX_M / wave_speed_m_per_s
}

/// Cell counts of the energy convergence study: the same basin resolved at
/// three resolutions, each a halving of the last.
const ENERGY_RESOLUTIONS: [usize; 3] = [16, 32, 64];

/// How far below second order the measured convergence of the energy drift may
/// sit.
///
/// The drift is dominated by the `O(Δy²)` skewness defect of the C-grid
/// Coriolis pair derived in the module comment; RK4's own share falls faster
/// still, since the CFL bound ties `dt` to the cell width. Second order is
/// therefore the rate the drift falls at, and the check is one-sided: falling
/// faster is not a defect, falling slower is.
///
/// The slack is for the sub-leading `O(Δy⁴)` correction, which is what keeps a
/// measurement at finite resolution off the asymptote. Writing the drift as
/// `A·Δy²·(1 + c·Δy²)`, a halving of `Δy` gives a ratio `4·(1 + (3/4)·c·Δy²)`,
/// so a correction worth a fraction `ε` of the leading term at the coarser of
/// the two resolutions moves the measured order by about `(3/4)·ε/ln 2`. A
/// tolerance of 0.2 therefore admits a sub-leading term up to about 18% of the
/// leading one — generous for the resolutions studied here, and still far from
/// the first-order or resolution-independent behaviour a mis-wired step would
/// show.
const ENERGY_ORDER_TOLERANCE: f64 = 0.2;

/// Largest relative energy drift the finest resolution of the convergence
/// study may show.
///
/// Not a fitted number: it is the statement that at 32×32 the drift is a small
/// correction rather than the leading behaviour, so that the orders measured
/// across the study are orders of a truncation error. One part in a thousand
/// over sixteen basin crossings is two orders of magnitude below the energy
/// itself.
///
/// It is a ceiling on *this* four-crossing run and nothing more. The bound the
/// drift is actually held to — derived from the C-grid Coriolis pair's
/// skewness defect and RK4's amplification polynomial, resolution-dependent,
/// and checked over a run eight times as long — is T-07.5's, in
/// `conservation.rs`.
const ENERGY_DRIFT_CEILING: f64 = 1.0e-3;

/// The largest relative departure of the energy from its initial value over an
/// undamped, unforced run of `cells`×`cells` at the CFL-safe timestep.
fn energy_drift(cells: usize) -> f64 {
    let (grid, spacing) = basin(cells, cells);
    let params = pacific_params(0.0);
    let dt_s = cfl_safe_dt_s(spacing, params);
    let steps = (energy_run_length_s() / dt_s).round() as usize;

    let mut state = gravest_zonal_mode(grid, spacing);
    let initial = wave_energy(&state, params);
    let mut solver = solver_for(grid, spacing, params, dt_s);
    let calm = WindStressField::calm(grid);

    let mut worst = 0.0_f64;
    for n in 0..steps {
        solver.step(&mut state, n as f64 * dt_s, |_t_s| &calm);
        worst = worst.max((wave_energy(&state, params) / initial - 1.0).abs());
    }
    worst
}

#[test]
fn an_undamped_unforced_run_conserves_energy_to_truncation_scale() {
    // The acceptance criterion. `r = 0` and `τ = 0` leave a system whose only
    // departure from exact energy conservation is its own truncation error, so
    // the drift must fall at the scheme's order as the basin is resolved — see
    // the module comment for why that order is two.
    let drifts: Vec<f64> = ENERGY_RESOLUTIONS
        .iter()
        .map(|&n| energy_drift(n))
        .collect();

    for (coarse, pair) in drifts.windows(2).enumerate() {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            order >= 2.0 - ENERGY_ORDER_TOLERANCE,
            "energy drift must fall at least at second order under refinement, but \
             {}x{} to {}x{} gives order {order} (drifts {drifts:?})",
            ENERGY_RESOLUTIONS[coarse],
            ENERGY_RESOLUTIONS[coarse],
            ENERGY_RESOLUTIONS[coarse + 1],
            ENERGY_RESOLUTIONS[coarse + 1],
        );
    }

    let finest = *drifts.last().expect("a non-empty study");
    assert!(
        finest < ENERGY_DRIFT_CEILING,
        "at {}x{} the energy drifted by a relative {finest} over the run",
        ENERGY_RESOLUTIONS[ENERGY_RESOLUTIONS.len() - 1],
        ENERGY_RESOLUTIONS[ENERGY_RESOLUTIONS.len() - 1],
    );
}

// --- Acceptance criterion: a long run at the CFL-safe dt stays finite. ---

/// Cells across the basin the stability runs use.
const STABILITY_BASIN_CELLS: usize = 16;
/// Steps of a stability run. At the CFL-safe timestep of this basin
/// (≈1.15×10⁵ s) that is about seven years of simulated time, and around
/// thirty crossings of the basin at the Kelvin wave speed.
const STABILITY_RUN_STEPS: usize = 2_000;

/// Factor by which the energy of an undamped run may exceed its initial value
/// before the run counts as blowing up.
///
/// Not a tolerance on a measurement: the energy of this system is conserved to
/// truncation scale (see the convergence test above), so a doubling could only
/// come from an instability amplifying a mode step after step. A scheme run
/// past its CFL bound grows by orders of magnitude over this many steps, not
/// by a factor of two.
const BLOW_UP_FACTOR: f64 = 2.0;

#[test]
fn a_multi_step_run_at_the_cfl_safe_timestep_stays_finite() {
    // The second acceptance criterion, in both the unforced and the forced
    // configuration the epic is scoped to: the run neither produces a NaN nor
    // amplifies.
    let (grid, spacing) = basin(STABILITY_BASIN_CELLS, STABILITY_BASIN_CELLS);

    let undamped = pacific_params(0.0);
    let dt_s = cfl_safe_dt_s(spacing, undamped);
    let mut state = gravest_zonal_mode(grid, spacing);
    let initial = wave_energy(&state, undamped);
    let mut solver = solver_for(grid, spacing, undamped, dt_s);
    let calm = WindStressField::calm(grid);
    for n in 0..STABILITY_RUN_STEPS {
        solver.step(&mut state, n as f64 * dt_s, |_t_s| &calm);
        assert_finite(&state, n);
        let energy = wave_energy(&state, undamped);
        assert!(
            energy < BLOW_UP_FACTOR * initial,
            "the unforced run amplified: energy went from {initial} to {energy} by step {n}"
        );
    }

    // The forced, damped configuration: a basin started at rest and pushed by
    // constant trade winds. Its energy rises — the wind is doing work — so
    // finiteness is what there is to assert here; where it settles is the
    // thermocline-tilt validation of Epic 07.
    let damped = pacific_params(STRONG_DAMPING_PER_S);
    let dt_s = cfl_safe_dt_s(spacing, damped);
    let mut state = OceanState::at_rest(grid);
    let mut solver = solver_for(grid, spacing, damped, dt_s);
    let trade_winds =
        WindStressField::uniform(grid, TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);
    for n in 0..STABILITY_RUN_STEPS {
        solver.step(&mut state, n as f64 * dt_s, |_t_s| &trade_winds);
        assert_finite(&state, n);
    }
}

/// Panic if any point of `state` has stopped being a finite number.
fn assert_finite(state: &OceanState, step: usize) {
    for (name, field) in [("h", state.h()), ("u", state.u()), ("v", state.v())] {
        for (flat, value) in field.as_slice().iter().enumerate() {
            assert!(
                value.is_finite(),
                "{name} went to {value} at flat index {flat} on step {step}"
            );
        }
    }
}

// --- The timestep is checked, not clamped. ---

#[test]
fn a_timestep_past_the_cfl_bound_is_refused_rather_than_shortened() {
    // CODING_STANDARDS.md § No silent clamping: a solver built on an unstable
    // timestep is a scenario error to report, not something to quietly fix.
    let (grid, spacing) = basin(STABILITY_BASIN_CELLS, STABILITY_BASIN_CELLS);
    let params = pacific_params(0.0);
    let bound_s = cfl_safe_dt_s(spacing, params);
    let plane = equatorial_plane(params, spacing, grid);

    let too_long_s = 2.0 * bound_s;
    let error = Solver::new(grid, spacing, params, plane, too_long_s)
        .expect_err("a timestep twice the CFL bound must be refused");
    assert_eq!(
        error,
        SolverError::Cfl(CflError::TimestepExceedsCfl {
            requested_s: too_long_s,
            max_stable_s: bound_s,
        })
    );

    // …and the bound itself is admissible, so the check is a bound rather than
    // an unreachable one.
    assert!(Solver::new(grid, spacing, params, plane, bound_s).is_ok());
}

/// Meridional extent of the basin the rotation bound is checked on, in metres
/// — the ±2500 km of `CONTEXT.md`'s Pacific basin.
const WIDE_BASIN_LY_M: f64 = 5.0e6;

#[test]
fn a_timestep_that_cannot_resolve_the_basins_rotation_is_refused() {
    // The gravity-wave CFL bound of T-01.3 is not the only limit on the step:
    // the rotation pair `u̇ = +f·v`, `v̇ = −f·u` oscillates at `|f| = β·|y|`,
    // and RK4 follows it only while `|f|·dt` stays inside the same stability
    // region. On the full-width Pacific basin at 16 cells across, the two
    // bounds disagree — the CFL bound admits a step longer than the inertial
    // period at the walls — and the solver must refuse rather than run.
    let grid = Grid::new(STABILITY_BASIN_CELLS, STABILITY_BASIN_CELLS).expect("non-zero extents");
    let spacing = Spacing::new(
        BASIN_LX_M / STABILITY_BASIN_CELLS as f64,
        WIDE_BASIN_LY_M / STABILITY_BASIN_CELLS as f64,
    )
    .expect("a basin spanned by whole cells has positive spacing");
    let params = pacific_params(0.0);
    let plane = BetaPlane::centered_on_equator(params, spacing, grid);
    let cfl_bound_s = cfl_safe_dt_s(spacing, params);

    // `|f|` at the walls is `β·Ly/2 = 5.75×10⁻⁵ s⁻¹`, so RK4's imaginary-axis
    // limit `2√2` — with the same safety factor the CFL bound holds back —
    // allows at most 0.8·2√2/5.75×10⁻⁵ = 3.94×10⁴ s, well under the CFL
    // bound's 1.15×10⁵ s.
    //
    // `CFL_SAFETY_FACTOR` is read from the crate rather than written out
    // because the claim under test is that the *same* margin governs both
    // bounds; its value is pinned independently by T-01.3's own tests. `2√2`
    // is RK4's imaginary-axis limit, written out from Hairer & Wanner.
    let largest_coriolis_per_s = BETA_PER_M_PER_S * WIDE_BASIN_LY_M / 2.0;
    let rotation_bound_s =
        engine::CFL_SAFETY_FACTOR * 2.0 * std::f64::consts::SQRT_2 / largest_coriolis_per_s;
    assert!(
        rotation_bound_s < cfl_bound_s,
        "this basin is meant to be one the rotation bound binds on, but {rotation_bound_s} s is \
         not shorter than the CFL bound of {cfl_bound_s} s"
    );

    let error = Solver::new(grid, spacing, params, plane, cfl_bound_s)
        .expect_err("a step longer than the rotation allows must be refused");
    assert_eq!(
        error,
        SolverError::TimestepExceedsRotationLimit {
            requested_s: cfl_bound_s,
            max_stable_s: rotation_bound_s,
            largest_coriolis_per_s,
        }
    );
    // Actionable, per CODING_STANDARDS.md § Correctness and failure: the
    // message names the value it rejected and the bound it violated.
    let message = error.to_string();
    assert!(message.contains(&rotation_bound_s.to_string()), "{message}");

    // …and the rotation bound itself is admissible, so it is a bound rather
    // than an unreachable one.
    assert!(Solver::new(grid, spacing, params, plane, rotation_bound_s).is_ok());
}

// --- Determinism. ---

#[test]
fn two_identical_runs_produce_identical_states() {
    // CODING_STANDARDS.md § Determinism: identical scenario in, identical
    // state out, to the last bit — no iteration-order dependence and no
    // unseeded state anywhere in the step.
    let (grid, spacing) = basin(STABILITY_BASIN_CELLS, STABILITY_BASIN_CELLS);
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let dt_s = cfl_safe_dt_s(spacing, params);
    let trade_winds =
        WindStressField::uniform(grid, TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);

    let run = || {
        let mut state = gravest_zonal_mode(grid, spacing);
        let mut solver = solver_for(grid, spacing, params, dt_s);
        for n in 0..DETERMINISM_RUN_STEPS {
            solver.step(&mut state, n as f64 * dt_s, |_t_s| &trade_winds);
        }
        state
    };

    assert_eq!(run(), run());
}

/// Steps of the determinism run: enough for a difference to have somewhere to
/// grow from, cheap enough to run twice.
const DETERMINISM_RUN_STEPS: usize = 64;

// --- The allocating wrapper and the reusable solver agree. ---

#[test]
fn the_convenience_step_is_the_same_computation_as_the_reusable_solver() {
    // `step` exists for tests and one-off evaluation; a time loop uses
    // `Solver` so that a run allocates its buffers once
    // (CODING_STANDARDS.md § Performance). The two must not drift apart, so
    // the wrapper's result is asserted bit-identical to the solver's.
    let (grid, spacing) = basin(CORIOLIS_BASIN_CELLS, CORIOLIS_BASIN_CELLS);
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let dt_s = cfl_safe_dt_s(spacing, params);
    let plane = equatorial_plane(params, spacing, grid);
    let trade_winds =
        WindStressField::uniform(grid, TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);
    let initial = gravest_zonal_mode(grid, spacing);

    let wrapped = step(&initial, dt_s, params, spacing, plane, |_t_s| &trade_winds)
        .expect("the CFL-safe timestep is admissible");

    let mut stepped = initial.clone();
    solver_for(grid, spacing, params, dt_s).step(&mut stepped, 0.0, |_t_s| &trade_winds);

    assert_eq!(wrapped, stepped);
}
