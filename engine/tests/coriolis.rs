//! Acceptance tests for T-02.2 — the beta-plane Coriolis parameter and the
//! Coriolis terms of the momentum equations.
//!
//! Every expected value here comes from an independent source: the beta-plane
//! definition `f = β·y` in `CONTEXT.md`, the momentum equations in
//! `docs/planning/01-scientific-model.md`, arithmetic done on paper for the
//! small-grid example, or the Taylor expansion of a centred average. None of
//! them was produced by running this code.
//!
//! The two Coriolis contributions the model asks for are, from
//! `∂u/∂t − f·v = …` and `∂v/∂t + f·u = …`:
//!
//! ```text
//! ∂u/∂t += +f·v      ∂v/∂t += −f·u
//! ```
//!
//! On the Arakawa C-grid of [ADR-0003] neither product is collocated: `v` has
//! to reach an east/west face and `u` a north/south face. The centred
//! four-point average of the surrounding opposite-face values lands exactly on
//! the target point, so it is second order — a centred two-point average over
//! `Δ` satisfies `(g(x+Δ/2) + g(x−Δ/2))/2 = g(x) + (Δ²/8)·g''(x) + O(Δ⁴)`, and
//! averaging in x and then in y just adds one such remainder per axis.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use std::f64::consts::PI;

use engine::{
    BetaPlane, CoriolisTerm, Field2D, Grid, OceanState, PhysicalParams, Spacing, H_STAGGERING,
    U_STAGGERING, V_STAGGERING,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s² (Gill, *Atmosphere–Ocean Dynamics*, ch. 11; Cane & Sarachik 1981).
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;
/// Rayleigh damping `r`, in s⁻¹: a damping timescale of about two years. The
/// Coriolis term does not read it, but `PhysicalParams` is validated as a set.
const PACIFIC_DAMPING_PER_S: f64 = 1.0 / (2.0 * 365.0 * 86_400.0);

/// The equatorial beta-plane gradient, in m⁻¹s⁻¹ — `CONTEXT.md`, *Beta-plane*.
const BETA_PER_M_PER_S: f64 = engine::EQUATORIAL_BETA_PER_M_PER_S;

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*). Deliberately different from
/// `BASIN_LY_M` so an x/y swap cannot pass.
const BASIN_LX_M: f64 = 1.0e7;
/// Meridional extent of the test basin, in metres.
const BASIN_LY_M: f64 = 5.0e6;

/// Relative tolerance for a result that is a handful of `f64` operations away
/// from the hand-computed value: four additions, one multiplication by 0.25
/// and one by `f`, so a few units in the last place. `8·ε ≈ 1.8e-15`.
const FEW_ULPS: f64 = 8.0 * f64::EPSILON;

fn pacific_params() -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        PACIFIC_DAMPING_PER_S,
        BETA_PER_M_PER_S,
        engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

/// Assert `actual` matches `expected` to within `FEW_ULPS` relative error.
fn assert_close(actual: f64, expected: f64, what: &str) {
    let tolerance = FEW_ULPS * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: got {actual:e}, expected {expected:e} (tolerance {tolerance:e})"
    );
}

#[test]
fn the_coriolis_parameter_is_beta_times_the_distance_from_the_equator() {
    // `f = β·y` (CONTEXT.md, *Beta-plane*). One cell height is 10⁶ m here, so
    // the u row 10⁶ m north of the equator carries
    // f = 2.3×10⁻¹¹ × 10⁶ = 2.3×10⁻⁵ s⁻¹ — the textbook mid-basin value.
    let grid = Grid::new(3, 3).expect("a 3x3 basin is non-degenerate");
    let spacing = Spacing::new(1.0e6, 1.0e6).expect("a 1000 km cell is positive");
    let plane = BetaPlane::centered_on_equator(pacific_params(), spacing, grid);

    // ny = 3 is odd, so the middle *cell-center* row sits on the equator and
    // the rows either side of it are one full cell height away.
    assert_close(
        plane.y_of_row_m(U_STAGGERING, 2),
        1.0e6,
        "y of the u row one cell north of the equator",
    );
    assert_close(
        plane.coriolis_at_row_per_s(U_STAGGERING, 2),
        2.3e-5,
        "f one cell north of the equator",
    );
    // The v rows are offset half a cell, so the nearest one north of the
    // equator is at 5×10⁵ m: f = 2.3×10⁻¹¹ × 5×10⁵ = 1.15×10⁻⁵ s⁻¹.
    assert_close(
        plane.coriolis_at_row_per_s(V_STAGGERING, 2),
        1.15e-5,
        "f half a cell north of the equator",
    );
}

#[test]
fn the_coriolis_parameter_vanishes_exactly_on_the_equator() {
    let spacing = Spacing::new(1.0e6, 1.0e6).expect("a 1000 km cell is positive");
    let params = pacific_params();

    // An odd number of rows puts a cell-center line — where `h` lives and
    // where the u rows sit — exactly on the equator.
    let odd = Grid::new(3, 3).expect("a 3x3 basin is non-degenerate");
    let odd_plane = BetaPlane::centered_on_equator(params, spacing, odd);
    assert_eq!(odd_plane.y_of_row_m(H_STAGGERING, 1), 0.0);
    assert_eq!(odd_plane.coriolis_at_row_per_s(H_STAGGERING, 1), 0.0);
    assert_eq!(odd_plane.coriolis_at_row_per_s(U_STAGGERING, 1), 0.0);

    // An even number puts a north/south face line — where `v` lives — there
    // instead. `f = 0` exactly at the equator in both staggerings: the whole
    // point of the beta-plane is that the equator is a zero, not a small
    // number.
    let even = Grid::new(3, 4).expect("a 3x4 basin is non-degenerate");
    let even_plane = BetaPlane::centered_on_equator(params, spacing, even);
    assert_eq!(even_plane.y_of_row_m(V_STAGGERING, 2), 0.0);
    assert_eq!(even_plane.coriolis_at_row_per_s(V_STAGGERING, 2), 0.0);
}

#[test]
fn the_coriolis_parameter_changes_sign_across_the_equator() {
    // `f = β·y` is odd in y: rows equidistant from the equator carry equal and
    // opposite values, positive in the northern hemisphere.
    let grid = Grid::new(3, 6).expect("a 3x6 basin is non-degenerate");
    let spacing = Spacing::new(1.0e6, 1.0e6).expect("a 1000 km cell is positive");
    let plane = BetaPlane::centered_on_equator(pacific_params(), spacing, grid);

    // ny = 6: the equator is the v row j = 3, so j = 3 ± k are mirror rows.
    for offset in 1..=3 {
        let north = plane.coriolis_at_row_per_s(V_STAGGERING, 3 + offset);
        let south = plane.coriolis_at_row_per_s(V_STAGGERING, 3 - offset);
        assert!(
            north > 0.0,
            "f must be positive north of the equator: {north:e}"
        );
        assert!(
            south < 0.0,
            "f must be negative south of the equator: {south:e}"
        );
        assert_close(north, -south, "f is odd in y");
    }
}

/// A 3×3 basin of 1000 km cells with the equator through its middle row, and
/// the Coriolis term over it.
fn small_basin() -> (Grid, CoriolisTerm) {
    let grid = Grid::new(3, 3).expect("a 3x3 basin is non-degenerate");
    let spacing = Spacing::new(1.0e6, 1.0e6).expect("a 1000 km cell is positive");
    let plane = BetaPlane::centered_on_equator(pacific_params(), spacing, grid);
    (grid, CoriolisTerm::new(grid, plane))
}

fn set(field: &mut Field2D<f64>, i: usize, j: usize, value: f64) {
    *field.get_mut(i, j).expect("point inside the field") = value;
}

fn at(field: &Field2D<f64>, i: usize, j: usize) -> f64 {
    *field.get(i, j).expect("point inside the field")
}

#[test]
fn the_zonal_tendency_is_a_hand_computed_four_point_average_of_v() {
    // Worked by hand. The u point (i=1, j=2) sits on the cell-center row one
    // cell north of the equator, so f = β·10⁶ = 2.3×10⁻⁵ s⁻¹. The four v
    // points around it are (0,2), (1,2), (0,3) and (1,3); set them to
    // 1, 2, 3 and 4 m/s, so their average is (1+2+3+4)/4 = 2.5 m/s and
    //     ∂u/∂t = +f·v̄ = 2.3×10⁻⁵ × 2.5 = 5.75×10⁻⁵ m/s².
    let (grid, coriolis) = small_basin();
    let mut state = OceanState::at_rest(grid);
    set(state.v_mut(), 0, 2, 1.0);
    set(state.v_mut(), 1, 2, 2.0);
    set(state.v_mut(), 0, 3, 3.0);
    set(state.v_mut(), 1, 3, 4.0);

    let mut tendency = OceanState::at_rest(grid);
    coriolis.add_to_tendency(&state, &mut tendency);

    assert_close(at(tendency.u(), 1, 2), 5.75e-5, "∂u/∂t at (1, 2)");
    // The v equation reads `u`, which is still at rest, so it stays at rest.
    assert!(tendency.v().as_slice().iter().all(|dv| *dv == 0.0));
    // Coriolis touches neither the thermocline nor the continuity equation.
    assert!(tendency.h().as_slice().iter().all(|dh| *dh == 0.0));
}

#[test]
fn the_meridional_tendency_is_a_hand_computed_four_point_average_of_u() {
    // Worked by hand. The v point (i=1, j=2) sits half a cell north of the
    // equator, so f = β·5×10⁵ = 1.15×10⁻⁵ s⁻¹. The four u points around it
    // are (1,1), (2,1), (1,2) and (2,2); set them to 1, 3, 5 and 7 m/s, whose
    // average is 4 m/s, and
    //     ∂v/∂t = −f·ū = −1.15×10⁻⁵ × 4 = −4.6×10⁻⁵ m/s².
    // The sign is the physical one: eastward flow north of the equator is
    // deflected to the right, i.e. southward.
    let (grid, coriolis) = small_basin();
    let mut state = OceanState::at_rest(grid);
    set(state.u_mut(), 1, 1, 1.0);
    set(state.u_mut(), 2, 1, 3.0);
    set(state.u_mut(), 1, 2, 5.0);
    set(state.u_mut(), 2, 2, 7.0);

    let mut tendency = OceanState::at_rest(grid);
    coriolis.add_to_tendency(&state, &mut tendency);

    assert_close(at(tendency.v(), 1, 2), -4.6e-5, "∂v/∂t at (1, 2)");
    assert!(tendency.u().as_slice().iter().all(|du| *du == 0.0));
    assert!(tendency.h().as_slice().iter().all(|dh| *dh == 0.0));
}

#[test]
fn eastward_flow_is_deflected_right_in_the_north_left_in_the_south_and_not_at_all_on_the_equator() {
    // The sign change the acceptance criteria ask for, seen through the term
    // rather than through `f` alone: a uniform eastward current is turned
    // equatorward in both hemispheres and left alone exactly on the equator.
    let grid = Grid::new(3, 4).expect("a 3x4 basin is non-degenerate");
    let spacing = Spacing::new(1.0e6, 1.0e6).expect("a 1000 km cell is positive");
    let plane = BetaPlane::centered_on_equator(pacific_params(), spacing, grid);
    let coriolis = CoriolisTerm::new(grid, plane);

    let mut state = OceanState::at_rest(grid);
    state.u_mut().as_mut_slice().fill(1.0);
    let mut tendency = OceanState::at_rest(grid);
    coriolis.add_to_tendency(&state, &mut tendency);

    // ny = 4 puts the equator on the v row j = 2; j = 1 is south of it and
    // j = 3 north. i = 1 is an interior column in every case.
    assert!(
        at(tendency.v(), 1, 3) < 0.0,
        "northern flow turns southward"
    );
    assert!(
        at(tendency.v(), 1, 1) > 0.0,
        "southern flow turns northward"
    );
    assert_eq!(
        at(tendency.v(), 1, 2),
        0.0,
        "f = 0 on the equator, so eastward flow is not deflected there"
    );
    assert_close(
        at(tendency.v(), 1, 3),
        -at(tendency.v(), 1, 1),
        "the deflection is antisymmetric about the equator",
    );
}

#[test]
fn the_rest_state_has_no_coriolis_tendency() {
    // `f·0 = 0`: an ocean at rest stays at rest under rotation alone.
    let (grid, coriolis) = small_basin();
    let state = OceanState::at_rest(grid);
    let mut tendency = OceanState::at_rest(grid);
    coriolis.add_to_tendency(&state, &mut tendency);
    assert_eq!(tendency, OceanState::at_rest(grid));
}

#[test]
fn the_closed_basin_walls_carry_no_coriolis_tendency() {
    // The basin is closed on all four sides (01-scientific-model.md), so the
    // normal velocity on a wall does not evolve: the western and eastern u
    // faces and the southern and northern v faces are left untouched.
    let (grid, coriolis) = small_basin();
    let mut state = OceanState::at_rest(grid);
    state.u_mut().as_mut_slice().fill(1.0);
    state.v_mut().as_mut_slice().fill(1.0);

    let mut tendency = OceanState::at_rest(grid);
    coriolis.add_to_tendency(&state, &mut tendency);

    for j in 0..tendency.u().ny() {
        assert_eq!(at(tendency.u(), 0, j), 0.0, "western wall at j = {j}");
        assert_eq!(
            at(tendency.u(), grid.nx(), j),
            0.0,
            "eastern wall at j = {j}"
        );
    }
    for i in 0..tendency.v().nx() {
        assert_eq!(at(tendency.v(), i, 0), 0.0, "southern wall at i = {i}");
        assert_eq!(
            at(tendency.v(), i, grid.ny()),
            0.0,
            "northern wall at i = {i}"
        );
    }
}

#[test]
fn the_term_accumulates_into_the_tendency_rather_than_overwriting_it() {
    // The Coriolis term is one contribution among several (pressure gradient,
    // friction, wind stress), so applying it twice must double it rather than
    // replace what a previous contribution wrote.
    let (grid, coriolis) = small_basin();
    let mut state = OceanState::at_rest(grid);
    set(state.v_mut(), 0, 2, 1.0);
    set(state.v_mut(), 1, 2, 2.0);
    set(state.v_mut(), 0, 3, 3.0);
    set(state.v_mut(), 1, 3, 4.0);

    let mut tendency = OceanState::at_rest(grid);
    coriolis.add_to_tendency(&state, &mut tendency);
    coriolis.add_to_tendency(&state, &mut tendency);

    assert_close(
        at(tendency.u(), 1, 2),
        2.0 * 5.75e-5,
        "twice ∂u/∂t at (1, 2)",
    );
}

/// Zonal wavenumber of the analytic test flow: one full wave across the basin.
fn wavenumber_x() -> f64 {
    2.0 * PI / BASIN_LX_M
}

/// Meridional wavenumber of the analytic test flow: one full wave across the
/// basin.
fn wavenumber_y() -> f64 {
    2.0 * PI / BASIN_LY_M
}

/// The analytic velocity field both components are sampled from,
/// `w(x, y) = sin(kx·x)·cos(ky·y)`. Smooth and non-separable in the sense that
/// matters here: its second derivatives are non-zero everywhere, so the
/// leading truncation term of a centred average does not accidentally vanish.
fn analytic_velocity_m_per_s(x_m: f64, y_m: f64) -> f64 {
    (wavenumber_x() * x_m).sin() * (wavenumber_y() * y_m).cos()
}

/// Largest error, over the interior points of an `n`-by-`n`-cell basin
/// spanning the fixed domain, between the discrete Coriolis tendencies and the
/// analytic `+f·v` and `−f·u` evaluated at the same points.
fn largest_coriolis_error(cells: usize) -> f64 {
    let grid = Grid::new(cells, cells).expect("a non-degenerate basin");
    let (dx_m, dy_m) = (BASIN_LX_M / cells as f64, BASIN_LY_M / cells as f64);
    let spacing = Spacing::new(dx_m, dy_m).expect("a positive cell size");
    // The basin straddles the equator, so `f` changes sign inside it.
    let southern_edge_y_m = -0.5 * BASIN_LY_M;
    let plane = BetaPlane::new(pacific_params(), spacing, southern_edge_y_m)
        .expect("a finite southern edge");
    let coriolis = CoriolisTerm::new(grid, plane);

    // Sample the same analytic field at the u points and at the v points.
    let mut state = OceanState::at_rest(grid);
    for j in 0..state.u().ny() {
        for i in 0..state.u().nx() {
            let y_m = plane.y_of_row_m(U_STAGGERING, j);
            set(
                state.u_mut(),
                i,
                j,
                analytic_velocity_m_per_s(i as f64 * dx_m, y_m),
            );
        }
    }
    for j in 0..state.v().ny() {
        for i in 0..state.v().nx() {
            let y_m = plane.y_of_row_m(V_STAGGERING, j);
            set(
                state.v_mut(),
                i,
                j,
                analytic_velocity_m_per_s((i as f64 + 0.5) * dx_m, y_m),
            );
        }
    }

    let mut tendency = OceanState::at_rest(grid);
    coriolis.add_to_tendency(&state, &mut tendency);

    let mut worst: f64 = 0.0;
    for j in 0..grid.ny() {
        for i in 1..grid.nx() {
            let y_m = plane.y_of_row_m(U_STAGGERING, j);
            let exact = plane.coriolis_at_row_per_s(U_STAGGERING, j)
                * analytic_velocity_m_per_s(i as f64 * dx_m, y_m);
            worst = worst.max((at(tendency.u(), i, j) - exact).abs());
        }
    }
    for j in 1..grid.ny() {
        for i in 0..grid.nx() {
            let y_m = plane.y_of_row_m(V_STAGGERING, j);
            let exact = -plane.coriolis_at_row_per_s(V_STAGGERING, j)
                * analytic_velocity_m_per_s((i as f64 + 0.5) * dx_m, y_m);
            worst = worst.max((at(tendency.v(), i, j) - exact).abs());
        }
    }
    worst
}

#[test]
fn the_staggered_interpolation_is_second_order_accurate() {
    // The four-point average is a centred average in x followed by one in y,
    // and a centred average over `Δ` carries a `(Δ²/8)·g''` remainder. So the
    // error is O(Δ²) and halving the spacing must quarter it. Asserting the
    // *order* rather than a single threshold is CODING_STANDARDS.md
    // § Tests, *Convergence over point checks*.
    let coarse = largest_coriolis_error(16);
    let fine = largest_coriolis_error(32);
    let finer = largest_coriolis_error(64);

    // 3.6 rather than 4.0 leaves room for the O(Δ⁴) term that is still
    // visible at 16 cells; it is far above the 2.0 a first-order scheme would
    // reach, which is what the assertion is really excluding.
    let first_refinement = coarse / fine;
    let second_refinement = fine / finer;
    assert!(
        first_refinement > 3.6,
        "16 → 32 cells reduced the error by {first_refinement:.3}×, not the ~4× of a \
         second-order scheme"
    );
    assert!(
        second_refinement > 3.6,
        "32 → 64 cells reduced the error by {second_refinement:.3}×, not the ~4× of a \
         second-order scheme"
    );
}

#[test]
fn a_basin_centred_on_the_equator_is_symmetric_about_it() {
    // `centered_on_equator` is the idealized configuration of Epic 02: half
    // the basin's meridional extent on each side of y = 0.
    let grid = Grid::new(3, 4).expect("a 3x4 basin is non-degenerate");
    let spacing = Spacing::new(1.0e6, 2.0e6).expect("positive cell sizes");
    let plane = BetaPlane::centered_on_equator(pacific_params(), spacing, grid);

    // 4 cells of 2000 km spans 8000 km, so the southern edge is 4000 km south.
    assert_close(plane.southern_edge_y_m(), -4.0e6, "southern edge");
    assert_close(plane.y_of_row_m(V_STAGGERING, 0), -4.0e6, "southern v row");
    assert_close(plane.y_of_row_m(V_STAGGERING, 4), 4.0e6, "northern v row");
}

#[test]
fn a_non_finite_southern_edge_is_rejected_with_the_value_it_carried() {
    // Where the basin sits is scenario input, so it comes back as a `Result`
    // naming the offending value (CODING_STANDARDS.md § Correctness).
    let spacing = Spacing::new(1.0e6, 1.0e6).expect("a 1000 km cell is positive");
    let error =
        BetaPlane::new(pacific_params(), spacing, f64::NAN).expect_err("NaN is not a position");
    let message = error.to_string();
    assert!(message.contains("southern_edge_y_m"), "{message}");
    assert!(message.contains("NaN"), "{message}");
}
