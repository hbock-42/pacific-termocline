//! Acceptance tests for T-02.3 — the pressure-gradient and continuity terms of
//! the 1.5-layer reduced-gravity shallow-water right-hand side.
//!
//! The equations under test are the linear ones of
//! `docs/planning/01-scientific-model.md`, restricted to the terms this ticket
//! owns:
//!
//! ```text
//! ∂u/∂t = −g'·∂h/∂x + τx/(ρ₀·H)
//! ∂v/∂t = −g'·∂h/∂y + τy/(ρ₀·H)
//! ∂h/∂t = −H·(∂u/∂x + ∂v/∂y)
//! ```
//!
//! The Rayleigh damping T-02.4 later folded into the same evaluation is
//! `engine/tests/rayleigh_damping.rs`' business, not this file's; it shows up
//! here only where a test's state makes it non-zero.
//!
//! Every expected value below is calculus done on paper — the derivatives of
//! `sin(kx·x)·sin(ky·y)` and of a Gaussian bump — never the output of running
//! this code. The spatial operators are the centred C-grid differences of
//! T-01.1, second order in the cell width (see `termocline-numerics`), so the
//! quantitative checks assert that order across three resolutions rather than
//! a single fixed threshold.

use std::f64::consts::PI;

use engine::{
    shallow_water_rhs, Field2D, Grid, OceanState, PhysicalParams, ShallowWaterRhs, Spacing,
    Staggering, WindStress, H_STAGGERING, U_STAGGERING, V_STAGGERING,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981).
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;
/// Rayleigh damping `r`, in s⁻¹: a damping timescale of about two years, the
/// order quoted for the equatorial Pacific. Only the tests whose state is
/// non-zero in the variable being checked see it; `engine/tests/
/// rayleigh_damping.rs` is where the damping term itself is pinned down.
const PACIFIC_DAMPING_PER_S: f64 = 1.0 / (2.0 * 365.0 * 86_400.0);

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

/// The parameter set every test here runs on.
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

/// Zonal wavenumber of the test fields: one full wave across the basin.
fn wavenumber_x() -> f64 {
    2.0 * PI / BASIN_LX_M
}

/// Meridional wavenumber of the test fields: one full wave across the basin.
fn wavenumber_y() -> f64 {
    2.0 * PI / BASIN_LY_M
}

/// `s(x, y) = sin(kx·x)·sin(ky·y)`, the shape every smooth test field takes.
fn s(x_m: f64, y_m: f64) -> f64 {
    (wavenumber_x() * x_m).sin() * (wavenumber_y() * y_m).sin()
}

/// `∂s/∂x = kx·cos(kx·x)·sin(ky·y)`, by hand.
fn dsdx(x_m: f64, y_m: f64) -> f64 {
    wavenumber_x() * (wavenumber_x() * x_m).cos() * (wavenumber_y() * y_m).sin()
}

/// `∂s/∂y = ky·sin(kx·x)·cos(ky·y)`, by hand.
fn dsdy(x_m: f64, y_m: f64) -> f64 {
    wavenumber_y() * (wavenumber_x() * x_m).sin() * (wavenumber_y() * y_m).cos()
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

/// Largest `|actual − expected|` over the points `keep` accepts.
fn max_error(
    actual: &Field2D<f64>,
    expected: &Field2D<f64>,
    keep: impl Fn(usize, usize) -> bool,
) -> f64 {
    assert_eq!((actual.nx(), actual.ny()), (expected.nx(), expected.ny()));
    let mut worst = 0.0_f64;
    for j in 0..actual.ny() {
        for i in 0..actual.nx() {
            if !keep(i, j) {
                continue;
            }
            let error = (actual.get(i, j).expect("in-bounds")
                - expected.get(i, j).expect("in-bounds"))
            .abs();
            worst = worst.max(error);
        }
    }
    worst
}

/// The east/west faces the pressure-gradient operator actually computes: those
/// with a cell on either side. The two boundary lines are the closed basin's
/// walls, which belong to Epic 04.
fn interior_x_faces(nx: usize) -> impl Fn(usize, usize) -> bool {
    move |i, _j| i > 0 && i < nx
}

/// As [`interior_x_faces`], for the north/south faces.
fn interior_y_faces(ny: usize) -> impl Fn(usize, usize) -> bool {
    move |_i, j| j > 0 && j < ny
}

/// Every point.
fn everywhere(_i: usize, _j: usize) -> bool {
    true
}

/// Order of accuracy implied by two errors measured at spacings differing by a
/// factor of two: `log2(coarse / fine)`.
fn observed_order(error_coarse: f64, error_fine: f64) -> f64 {
    (error_coarse / error_fine).log2()
}

/// Resolutions used for every convergence check. Three points, each a halving
/// of the last, so two independent order estimates are available.
const RESOLUTIONS: [usize; 3] = [16, 32, 64];

/// How far a measured order may sit from the theoretical 2.
///
/// The bound is the one T-01.1 derived for the same operators on the same test
/// function (`termocline-numerics/tests/operators.rs`): the next Taylor term
/// shifts the ratio by under 0.002 at these resolutions, and sampling the
/// remainder's peak more closely as the grid refines shifts it by about 0.021.
/// Wide enough to absorb both, narrow enough that a first- or third-order
/// scheme cannot pass. Scaling a field by a constant — `−g'` or `−H` — changes
/// neither the order nor the ratio.
const ORDER_TOLERANCE: f64 = 0.1;

/// Assert that `measure(n)` shrinks like `Δ²` across [`RESOLUTIONS`].
fn assert_second_order(measure: impl Fn(usize) -> f64) {
    let errors: Vec<f64> = RESOLUTIONS.iter().map(|&n| measure(n)).collect();
    for window in errors.windows(2) {
        let order = observed_order(window[0], window[1]);
        assert!(
            (order - 2.0).abs() < ORDER_TOLERANCE,
            "expected second-order convergence, measured order {order} from errors {errors:?}"
        );
    }
}

// --- Acceptance criterion: a Gaussian bump accelerates the flow outward. ---

/// Cells across the small test grid the acceptance criterion asks for, and the
/// basin every non-convergence test below runs on. Odd, so the Gaussian bump
/// can sit on the center of the middle cell and no face falls exactly on its
/// crest.
const SMALL_BASIN_CELLS: usize = 9;
/// `e`-folding width of the bump, in cell widths. Two cells is wide enough for
/// the centred difference to see a smooth hill and narrow enough that the bump
/// decays well inside the basin.
const BUMP_WIDTH_CELLS: f64 = 2.0;

#[test]
fn a_gaussian_bump_accelerates_the_current_outward() {
    // A deeper-than-average thermocline (`h > 0`) is a high-pressure hill in
    // the reduced-gravity system: `∂u/∂t = −g'·∂h/∂x` pushes water down the
    // hill, so east of the crest (`∂h/∂x < 0`) the flow accelerates eastward
    // and west of it westward. Both signs are read from the analytic gradient
    // of the bump, `∂h/∂x = −2(x−x₀)/L²·h`, not from the code.
    let (grid, spacing) = basin(SMALL_BASIN_CELLS);
    let crest = crest_position_m(spacing);
    let mut state = OceanState::at_rest(grid);
    *state.h_mut() = sample(grid, spacing, H_STAGGERING, |x_m, y_m| {
        gaussian_bump_m(spacing, crest, x_m, y_m)
    });

    let tendency = shallow_water_rhs(&state, pacific_params(), spacing, &WindStress::calm(grid));

    for j in 0..tendency.u().ny() {
        for i in 1..tendency.u().nx() - 1 {
            let (x_m, _) = position_m(spacing, U_STAGGERING, i, j);
            let acceleration = *tendency.u().get(i, j).expect("in-bounds");
            assert!(
                (x_m - crest.0).signum() == acceleration.signum(),
                "u face ({i}, {j}) at x = {x_m} m sits {} the crest but accelerates by {acceleration} m/s²",
                if x_m > crest.0 { "east of" } else { "west of" }
            );
        }
    }

    for j in 1..tendency.v().ny() - 1 {
        for i in 0..tendency.v().nx() {
            let (_, y_m) = position_m(spacing, V_STAGGERING, i, j);
            let acceleration = *tendency.v().get(i, j).expect("in-bounds");
            assert!(
                (y_m - crest.1).signum() == acceleration.signum(),
                "v face ({i}, {j}) at y = {y_m} m sits {} the crest but accelerates by {acceleration} m/s²",
                if y_m > crest.1 { "north of" } else { "south of" }
            );
        }
    }

    // The bump has not started moving yet: with `u = v = 0` the divergence is
    // identically zero — a sum of zeros, with no rounding to absorb — so the
    // whole of `∂h/∂t` is the Rayleigh damping T-02.4 folded in, `−r·h`,
    // exact to the one multiplication it costs.
    let params = pacific_params();
    for j in 0..tendency.h().ny() {
        for i in 0..tendency.h().nx() {
            let (x_m, y_m) = position_m(spacing, H_STAGGERING, i, j);
            let expected =
                -params.rayleigh_damping_per_s() * gaussian_bump_m(spacing, crest, x_m, y_m);
            let rate = *tendency.h().get(i, j).expect("in-bounds");
            assert!(
                (rate - expected).abs() <= ROUNDING_TOLERANCE * expected.abs(),
                "cell ({i}, {j}): a motionless ocean's thickness tendency is the damping alone, {expected} m/s, but it is {rate} m/s"
            );
        }
    }
}

/// Crest of the test bump: the center of the middle cell of the basin.
fn crest_position_m(spacing: Spacing) -> (f64, f64) {
    position_m(
        spacing,
        H_STAGGERING,
        SMALL_BASIN_CELLS / 2,
        SMALL_BASIN_CELLS / 2,
    )
}

/// `h(x, y) = A·exp(−((x−x₀)² + (y−y₀)²)/L²)`, in metres.
fn gaussian_bump_m(spacing: Spacing, crest: (f64, f64), x_m: f64, y_m: f64) -> f64 {
    let width_x_m = BUMP_WIDTH_CELLS * spacing.dx_m();
    let width_y_m = BUMP_WIDTH_CELLS * spacing.dy_m();
    let east_m = (x_m - crest.0) / width_x_m;
    let north_m = (y_m - crest.1) / width_y_m;
    H_AMPLITUDE_M * (-(east_m * east_m + north_m * north_m)).exp()
}

// --- The three terms against their analytic values. ---

#[test]
fn the_zonal_pressure_gradient_converges_on_minus_g_prime_dhdx() {
    let params = pacific_params();
    assert_second_order(|n| {
        let (grid, spacing) = basin(n);
        let mut state = OceanState::at_rest(grid);
        *state.h_mut() = sample(grid, spacing, H_STAGGERING, |x_m, y_m| {
            H_AMPLITUDE_M * s(x_m, y_m)
        });
        let expected = sample(grid, spacing, U_STAGGERING, |x_m, y_m| {
            -params.reduced_gravity_m_per_s2() * H_AMPLITUDE_M * dsdx(x_m, y_m)
        });

        let tendency = shallow_water_rhs(&state, params, spacing, &WindStress::calm(grid));
        max_error(tendency.u(), &expected, interior_x_faces(n))
    });
}

#[test]
fn the_meridional_pressure_gradient_converges_on_minus_g_prime_dhdy() {
    let params = pacific_params();
    assert_second_order(|n| {
        let (grid, spacing) = basin(n);
        let mut state = OceanState::at_rest(grid);
        *state.h_mut() = sample(grid, spacing, H_STAGGERING, |x_m, y_m| {
            H_AMPLITUDE_M * s(x_m, y_m)
        });
        let expected = sample(grid, spacing, V_STAGGERING, |x_m, y_m| {
            -params.reduced_gravity_m_per_s2() * H_AMPLITUDE_M * dsdy(x_m, y_m)
        });

        let tendency = shallow_water_rhs(&state, params, spacing, &WindStress::calm(grid));
        max_error(tendency.v(), &expected, interior_y_faces(n))
    });
}

#[test]
fn the_continuity_term_converges_on_minus_big_h_times_the_divergence() {
    let params = pacific_params();
    assert_second_order(|n| {
        let (grid, spacing) = basin(n);
        let mut state = OceanState::at_rest(grid);
        *state.u_mut() = sample(grid, spacing, U_STAGGERING, |x_m, y_m| {
            U_AMPLITUDE_M_PER_S * s(x_m, y_m)
        });
        *state.v_mut() = sample(grid, spacing, V_STAGGERING, |x_m, y_m| {
            V_AMPLITUDE_M_PER_S * s(x_m, y_m)
        });
        // `∂h/∂t = −H·(∂u/∂x + ∂v/∂y)`, differentiated on paper. A face→center
        // difference is centred at every cell, so there is no boundary gap
        // here and the check runs over the whole field.
        let expected = sample(grid, spacing, H_STAGGERING, |x_m, y_m| {
            -params.mean_thermocline_depth_m()
                * (U_AMPLITUDE_M_PER_S * dsdx(x_m, y_m) + V_AMPLITUDE_M_PER_S * dsdy(x_m, y_m))
        });

        let tendency = shallow_water_rhs(&state, params, spacing, &WindStress::calm(grid));
        max_error(tendency.h(), &expected, everywhere)
    });
}

/// Relative slack allowed where a check is exact in exact arithmetic: a few
/// ulps of `f64` (ε ≈ 2.2e-16) for the handful of operations each point costs.
const ROUNDING_TOLERANCE: f64 = 1.0e-14;

/// Zonal wind stress of the test forcing, in Pa. Easterly trade winds are
/// `τx < 0` (`CONTEXT.md`), so a negative value is the physical sign and the
/// one that must drive the current westward.
const TRADE_WIND_STRESS_X_PA: f64 = -0.05;
/// Meridional wind stress of the test forcing, in Pa. Different in magnitude
/// from [`TRADE_WIND_STRESS_X_PA`] so an x/y swap cannot pass.
const TRADE_WIND_STRESS_Y_PA: f64 = 0.02;

#[test]
fn a_uniform_wind_stress_accelerates_the_current_by_tau_over_rho_h() {
    // With a flat thermocline the pressure gradient vanishes and the momentum
    // equations reduce to the body force a surface stress exerts on a layer of
    // thickness `H`: `∂u/∂t = τx/(ρ₀·H)` (Gill, ch. 9). Easterly stress
    // (`τx < 0`) must therefore drive a westward current.
    let (grid, spacing) = basin(SMALL_BASIN_CELLS);
    let params = pacific_params();
    let state = OceanState::at_rest(grid);
    let wind = WindStress::uniform(grid, TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);

    let tendency = shallow_water_rhs(&state, params, spacing, &wind);

    let layer_mass_kg_per_m2 =
        params.reference_density_kg_per_m3() * params.mean_thermocline_depth_m();
    let expected_u_m_per_s2 = TRADE_WIND_STRESS_X_PA / layer_mass_kg_per_m2;
    let expected_v_m_per_s2 = TRADE_WIND_STRESS_Y_PA / layer_mass_kg_per_m2;
    assert!(expected_u_m_per_s2 < 0.0, "easterly stress drives westward");

    for rate in tendency.u().as_slice() {
        assert!(
            (rate - expected_u_m_per_s2).abs() <= ROUNDING_TOLERANCE * expected_u_m_per_s2.abs(),
            "expected {expected_u_m_per_s2} m/s², got {rate} m/s²"
        );
    }
    for rate in tendency.v().as_slice() {
        assert!(
            (rate - expected_v_m_per_s2).abs() <= ROUNDING_TOLERANCE * expected_v_m_per_s2.abs(),
            "expected {expected_v_m_per_s2} m/s², got {rate} m/s²"
        );
    }
}

#[test]
fn a_calm_ocean_at_rest_stays_at_rest() {
    // Every term of this ticket's right-hand side is linear and homogeneous in
    // the state except the wind stress, so the rest state with no wind is a
    // fixed point — exactly, since every product involved has a zero factor.
    let (grid, spacing) = basin(SMALL_BASIN_CELLS);
    let state = OceanState::at_rest(grid);

    let tendency = shallow_water_rhs(&state, pacific_params(), spacing, &WindStress::calm(grid));

    assert_eq!(tendency, OceanState::at_rest(grid));
}

#[test]
fn the_tendency_carries_the_staggering_of_the_state_it_came_from() {
    // A tendency is an `OceanState` of per-second units, so `∂h/∂t` sits at
    // cell centers with `h`, and each acceleration on the face its velocity
    // lives on — which is what lets RK4 combine the two (T-01.2).
    let (grid, spacing) = basin(SMALL_BASIN_CELLS);
    let state = OceanState::at_rest(grid);

    let tendency = shallow_water_rhs(&state, pacific_params(), spacing, &WindStress::calm(grid));

    assert_eq!(tendency.grid(), grid);
    assert_eq!(
        (tendency.h().nx(), tendency.h().ny()),
        grid.field_shape(H_STAGGERING)
    );
    assert_eq!(
        (tendency.u().nx(), tendency.u().ny()),
        grid.field_shape(U_STAGGERING)
    );
    assert_eq!(
        (tendency.v().nx(), tendency.v().ny()),
        grid.field_shape(V_STAGGERING)
    );
}

#[test]
fn a_reused_evaluator_writes_every_point_of_its_output() {
    // The time loop reuses one output buffer across steps and stages
    // (CODING_STANDARDS.md § Performance), so a right-hand side that left any
    // point untouched would leak the previous stage's tendency into this one.
    // Seeding the buffer with NaN catches that: an unwritten point stays NaN,
    // and NaN survives every later arithmetic operation.
    //
    // This is an equivalence check between the two entry points and a check
    // that every point is written — not a check of the values, which the
    // analytic tests above own. The finiteness assertion is what makes it fail
    // on a missed point, since `NaN != NaN` would also make the comparison
    // below fail for the wrong reason.
    let (grid, spacing) = basin(SMALL_BASIN_CELLS);
    let params = pacific_params();
    let crest = crest_position_m(spacing);
    let mut state = OceanState::at_rest(grid);
    *state.h_mut() = sample(grid, spacing, H_STAGGERING, |x_m, y_m| {
        gaussian_bump_m(spacing, crest, x_m, y_m)
    });
    let wind = WindStress::uniform(grid, TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);
    let expected = shallow_water_rhs(&state, params, spacing, &wind);

    let mut evaluator = ShallowWaterRhs::new(grid, spacing, params);
    let mut tendency = OceanState::at_rest(grid);
    tendency.h_mut().as_mut_slice().fill(f64::NAN);
    tendency.u_mut().as_mut_slice().fill(f64::NAN);
    tendency.v_mut().as_mut_slice().fill(f64::NAN);

    evaluator.evaluate(&state, &wind, &mut tendency);

    for rates in [
        tendency.h().as_slice(),
        tendency.u().as_slice(),
        tendency.v().as_slice(),
    ] {
        assert!(
            rates.iter().all(|rate| rate.is_finite()),
            "a point of the seeded buffer was left unwritten"
        );
    }
    assert_eq!(tendency, expected);
}

#[test]
fn the_basin_walls_carry_the_wind_stress_alone() {
    // A center→face difference needs a cell on either side, so the four
    // boundary faces have no pressure gradient to speak of and the operators
    // write zero there (T-01.1). The acceleration left on a wall is therefore
    // exactly the surface stress term, until Epic 04 gives the boundary a
    // condition of its own — a documented contract of the right-hand side, so
    // it is asserted rather than left to prose.
    let (grid, spacing) = basin(SMALL_BASIN_CELLS);
    let params = pacific_params();
    let crest = crest_position_m(spacing);
    let mut state = OceanState::at_rest(grid);
    *state.h_mut() = sample(grid, spacing, H_STAGGERING, |x_m, y_m| {
        gaussian_bump_m(spacing, crest, x_m, y_m)
    });
    let wind = WindStress::uniform(grid, TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);

    let tendency = shallow_water_rhs(&state, params, spacing, &wind);

    let layer_mass_kg_per_m2 =
        params.reference_density_kg_per_m3() * params.mean_thermocline_depth_m();
    let expected_u_m_per_s2 = TRADE_WIND_STRESS_X_PA / layer_mass_kg_per_m2;
    let expected_v_m_per_s2 = TRADE_WIND_STRESS_Y_PA / layer_mass_kg_per_m2;

    for j in 0..tendency.u().ny() {
        for i in [0, tendency.u().nx() - 1] {
            let rate = *tendency.u().get(i, j).expect("in-bounds");
            assert!(
                (rate - expected_u_m_per_s2).abs() <= ROUNDING_TOLERANCE * expected_u_m_per_s2.abs(),
                "west/east wall face ({i}, {j}): expected {expected_u_m_per_s2} m/s², got {rate} m/s²"
            );
        }
    }
    for j in [0, tendency.v().ny() - 1] {
        for i in 0..tendency.v().nx() {
            let rate = *tendency.v().get(i, j).expect("in-bounds");
            assert!(
                (rate - expected_v_m_per_s2).abs() <= ROUNDING_TOLERANCE * expected_v_m_per_s2.abs(),
                "south/north wall face ({i}, {j}): expected {expected_v_m_per_s2} m/s², got {rate} m/s²"
            );
        }
    }
}
