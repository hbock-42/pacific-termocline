//! Acceptance tests for T-01.1 — C-grid spatial derivative operators.
//!
//! Every expected value here comes from calculus done on paper, never from
//! running the operators: the test fields are samples of
//! `f(x, y) = sin(kx·x)·sin(ky·y)`, whose derivatives and Taylor remainders
//! are known in closed form.
//!
//! Two facts about the Arakawa C-grid ([ADR-0003]) drive the expected orders:
//!
//! - The midpoint of two neighbouring cell centers is exactly the face between
//!   them, and the midpoint of two neighbouring faces is exactly the center
//!   between them. So every operator below is a *centred* difference or a
//!   *centred* two-point average over one cell width.
//! - A centred difference over `Δ` satisfies
//!   `(f(x+Δ/2) − f(x−Δ/2))/Δ = f'(x) + (Δ²/24)·f'''(x) + O(Δ⁴)`, and a
//!   centred average satisfies
//!   `(f(x+Δ/2) + f(x−Δ/2))/2 = f(x) + (Δ²/8)·f''(x) + O(Δ⁴)`.
//!
//! Both are second order, so halving the spacing must quarter the error.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use std::f64::consts::PI;
use termocline_grid::Axis;
use termocline_grid::{Field2D, Grid, Staggering, H_STAGGERING, U_STAGGERING, V_STAGGERING};
use termocline_numerics::{CGridOperators, Spacing, SpacingError};

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (see CONTEXT.md, *Basin*). Nothing depends on the value;
/// it is deliberately different from `BASIN_LY_M` so an x/y swap cannot pass.
const BASIN_LX_M: f64 = 1.0e7;
/// Meridional extent of the test basin, in metres.
const BASIN_LY_M: f64 = 5.0e6;

/// Zonal wavenumber of the test field: one full wave across the basin.
fn wavenumber_x() -> f64 {
    2.0 * PI / BASIN_LX_M
}

/// Meridional wavenumber of the test field: one full wave across the basin.
fn wavenumber_y() -> f64 {
    2.0 * PI / BASIN_LY_M
}

/// `f(x, y) = sin(kx·x)·sin(ky·y)`.
fn f(x_m: f64, y_m: f64) -> f64 {
    (wavenumber_x() * x_m).sin() * (wavenumber_y() * y_m).sin()
}

/// `∂f/∂x = kx·cos(kx·x)·sin(ky·y)`, by hand.
fn dfdx(x_m: f64, y_m: f64) -> f64 {
    wavenumber_x() * (wavenumber_x() * x_m).cos() * (wavenumber_y() * y_m).sin()
}

/// `∂f/∂y = ky·sin(kx·x)·cos(ky·y)`, by hand.
fn dfdy(x_m: f64, y_m: f64) -> f64 {
    wavenumber_y() * (wavenumber_x() * x_m).sin() * (wavenumber_y() * y_m).cos()
}

/// A basin of `nx` by `ny` cells spanning `BASIN_LX_M` by `BASIN_LY_M`, and
/// the operators over it.
fn rectangular_basin(nx: usize, ny: usize) -> (Grid, Spacing, CGridOperators) {
    let grid = Grid::new(nx, ny).expect("extents are non-zero");
    let spacing = Spacing::new(BASIN_LX_M / nx as f64, BASIN_LY_M / ny as f64)
        .expect("a basin spanned by whole cells has positive spacing");
    (grid, spacing, CGridOperators::new(grid, spacing))
}

/// The square-celled case of [`rectangular_basin`]. `dx ≠ dy` still, since the
/// basin itself is not square, so an x/y swap cannot pass.
fn basin(n: usize) -> (Grid, Spacing, CGridOperators) {
    rectangular_basin(n, n)
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

/// Points a center→face operator actually computes: the faces with a cell on
/// both sides. The two boundary lines are the basin's edges, which belong to
/// Epic 04.
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
/// of the last, so two independent order estimates are available — the
/// acceptance criterion asks for at least two resolutions.
const RESOLUTIONS: [usize; 3] = [16, 32, 64];

/// How far a measured order may sit from the theoretical 2.
///
/// Two things move the measured ratio off 2, and the smaller one is the
/// obvious one: the next Taylor term is down by `(kΔ)²/80 ≤ 0.002` at these
/// resolutions. The larger is that the error is a maximum over *sampled*
/// points, and the remainder's shape `cos(kx·x)·sin(ky·y)` is sampled closer
/// to its true peak as the grid refines — 0.9808 at n = 16 against 0.9952 at
/// n = 32, which alone shifts `log2` of the ratio by about 0.021. Ten times
/// the Taylor effect, and still twenty times inside this band, which is set
/// wide enough to absorb both and narrow enough that a first- or third-order
/// scheme cannot pass.
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

// --- Acceptance criterion 1: analytic functions, at the right order. ---

#[test]
fn ddx_center_to_face_converges_at_second_order() {
    assert_second_order(|n| {
        let (grid, spacing, ops) = basin(n);
        let center = sample(grid, spacing, H_STAGGERING, f);
        let expected = sample(grid, spacing, U_STAGGERING, dfdx);
        let mut face = grid.allocate(U_STAGGERING, 0.0_f64);
        ops.ddx_center_to_face(&center, &mut face);
        max_error(&face, &expected, interior_x_faces(n))
    });
}

#[test]
fn ddy_center_to_face_converges_at_second_order() {
    assert_second_order(|n| {
        let (grid, spacing, ops) = basin(n);
        let center = sample(grid, spacing, H_STAGGERING, f);
        let expected = sample(grid, spacing, V_STAGGERING, dfdy);
        let mut face = grid.allocate(V_STAGGERING, 0.0_f64);
        ops.ddy_center_to_face(&center, &mut face);
        max_error(&face, &expected, interior_y_faces(n))
    });
}

#[test]
fn ddx_face_to_center_converges_at_second_order() {
    assert_second_order(|n| {
        let (grid, spacing, ops) = basin(n);
        let face = sample(grid, spacing, U_STAGGERING, f);
        let expected = sample(grid, spacing, H_STAGGERING, dfdx);
        let mut center = grid.allocate(H_STAGGERING, 0.0_f64);
        ops.ddx_face_to_center(&face, &mut center);
        max_error(&center, &expected, everywhere)
    });
}

#[test]
fn ddy_face_to_center_converges_at_second_order() {
    assert_second_order(|n| {
        let (grid, spacing, ops) = basin(n);
        let face = sample(grid, spacing, V_STAGGERING, f);
        let expected = sample(grid, spacing, H_STAGGERING, dfdy);
        let mut center = grid.allocate(H_STAGGERING, 0.0_f64);
        ops.ddy_face_to_center(&face, &mut center);
        max_error(&center, &expected, everywhere)
    });
}

#[test]
fn center_to_face_interpolation_converges_at_second_order() {
    assert_second_order(|n| {
        let (grid, spacing, ops) = basin(n);
        let center = sample(grid, spacing, H_STAGGERING, f);

        let expected_u = sample(grid, spacing, U_STAGGERING, f);
        let mut on_u = grid.allocate(U_STAGGERING, 0.0_f64);
        ops.center_to_face_x(&center, &mut on_u);

        let expected_v = sample(grid, spacing, V_STAGGERING, f);
        let mut on_v = grid.allocate(V_STAGGERING, 0.0_f64);
        ops.center_to_face_y(&center, &mut on_v);

        max_error(&on_u, &expected_u, interior_x_faces(n)).max(max_error(
            &on_v,
            &expected_v,
            interior_y_faces(n),
        ))
    });
}

#[test]
fn face_to_center_interpolation_converges_at_second_order() {
    assert_second_order(|n| {
        let (grid, spacing, ops) = basin(n);
        let expected = sample(grid, spacing, H_STAGGERING, f);

        let u = sample(grid, spacing, U_STAGGERING, f);
        let mut from_u = grid.allocate(H_STAGGERING, 0.0_f64);
        ops.face_to_center_x(&u, &mut from_u);

        let v = sample(grid, spacing, V_STAGGERING, f);
        let mut from_v = grid.allocate(H_STAGGERING, 0.0_f64);
        ops.face_to_center_y(&v, &mut from_v);

        max_error(&from_u, &expected, everywhere).max(max_error(&from_v, &expected, everywhere))
    });
}

/// Fraction by which a measured error may exceed the leading Taylor remainder
/// before the test fails. The next term in the series is smaller by
/// `(kΔ)²/80 ≤ 0.002` at these resolutions, so 5% is slack, not licence.
const TRUNCATION_SLACK: f64 = 0.05;
/// How much of the leading remainder the measured error must reach. The
/// remainder's shape here is `cos(kx·x)·sin(ky·y)`, whose magnitude comes
/// within a fraction of a percent of 1 somewhere on a 32×32 sample, so an
/// error far below the bound would mean the operator is not the second-order
/// difference it claims to be.
const TRUNCATION_FLOOR: f64 = 0.8;

/// Resolution the remainder checks run at. One value is enough: convergence
/// *order* is what needs several resolutions, and the tests above cover it.
const REMAINDER_RESOLUTION: usize = 32;

/// Leading remainder of a centred difference along x, `(dx²/24)·max|∂³f/∂x³|`,
/// with `|∂³f/∂x³| = kx³·|cos(kx·x)·sin(ky·y)| ≤ kx³`.
fn difference_remainder_x(spacing: Spacing) -> f64 {
    spacing.dx_m().powi(2) * wavenumber_x().powi(3) / 24.0
}

/// Leading remainder of a centred difference along y.
fn difference_remainder_y(spacing: Spacing) -> f64 {
    spacing.dy_m().powi(2) * wavenumber_y().powi(3) / 24.0
}

/// Leading remainder of a centred two-point average along x,
/// `(dx²/8)·max|∂²f/∂x²|`, with `|∂²f/∂x²| = kx²·|f| ≤ kx²`.
fn average_remainder_x(spacing: Spacing) -> f64 {
    spacing.dx_m().powi(2) * wavenumber_x().powi(2) / 8.0
}

/// Leading remainder of a centred two-point average along y.
fn average_remainder_y(spacing: Spacing) -> f64 {
    spacing.dy_m().powi(2) * wavenumber_y().powi(2) / 8.0
}

/// Assert that one operator's error sits *on* the analytic Taylor remainder,
/// not merely under some threshold — an operator can converge at the right
/// order and still carry the wrong constant.
fn assert_error_matches_remainder(
    operator: &str,
    input_staggering: Staggering,
    output_staggering: Staggering,
    input: fn(f64, f64) -> f64,
    expected: fn(f64, f64) -> f64,
    remainder: fn(Spacing) -> f64,
    apply: fn(&CGridOperators, &Field2D<f64>, &mut Field2D<f64>),
) {
    let n = REMAINDER_RESOLUTION;
    let (grid, spacing, ops) = basin(n);
    let input_field = sample(grid, spacing, input_staggering, input);
    let expected_field = sample(grid, spacing, output_staggering, expected);
    let mut output = grid.allocate(output_staggering, f64::NAN);
    apply(&ops, &input_field, &mut output);

    // A face output is only defined where a cell sits on both sides.
    let measured = match output_staggering {
        Staggering::CellCenter => max_error(&output, &expected_field, everywhere),
        Staggering::EastWestFace => max_error(&output, &expected_field, interior_x_faces(n)),
        Staggering::NorthSouthFace => max_error(&output, &expected_field, interior_y_faces(n)),
    };
    let bound = remainder(spacing);
    assert!(
        measured <= bound * (1.0 + TRUNCATION_SLACK),
        "{operator}: error {measured} exceeds the analytic remainder {bound}"
    );
    assert!(
        measured >= bound * TRUNCATION_FLOOR,
        "{operator}: error {measured} is implausibly far below the analytic remainder {bound}"
    );
}

#[test]
fn every_derivative_error_matches_the_analytic_truncation_remainder() {
    assert_error_matches_remainder(
        "ddx_center_to_face",
        H_STAGGERING,
        U_STAGGERING,
        f,
        dfdx,
        difference_remainder_x,
        CGridOperators::ddx_center_to_face,
    );
    assert_error_matches_remainder(
        "ddy_center_to_face",
        H_STAGGERING,
        V_STAGGERING,
        f,
        dfdy,
        difference_remainder_y,
        CGridOperators::ddy_center_to_face,
    );
    assert_error_matches_remainder(
        "ddx_face_to_center",
        U_STAGGERING,
        H_STAGGERING,
        f,
        dfdx,
        difference_remainder_x,
        CGridOperators::ddx_face_to_center,
    );
    assert_error_matches_remainder(
        "ddy_face_to_center",
        V_STAGGERING,
        H_STAGGERING,
        f,
        dfdy,
        difference_remainder_y,
        CGridOperators::ddy_face_to_center,
    );
}

#[test]
fn every_interpolation_error_matches_the_analytic_truncation_remainder() {
    assert_error_matches_remainder(
        "center_to_face_x",
        H_STAGGERING,
        U_STAGGERING,
        f,
        f,
        average_remainder_x,
        CGridOperators::center_to_face_x,
    );
    assert_error_matches_remainder(
        "center_to_face_y",
        H_STAGGERING,
        V_STAGGERING,
        f,
        f,
        average_remainder_y,
        CGridOperators::center_to_face_y,
    );
    assert_error_matches_remainder(
        "face_to_center_x",
        U_STAGGERING,
        H_STAGGERING,
        f,
        f,
        average_remainder_x,
        CGridOperators::face_to_center_x,
    );
    assert_error_matches_remainder(
        "face_to_center_y",
        V_STAGGERING,
        H_STAGGERING,
        f,
        f,
        average_remainder_y,
        CGridOperators::face_to_center_y,
    );
}

// --- Acceptance criterion 2: generic over grid size. ---

/// Deliberately awkward basin shapes: non-square, prime, and one-cell-thin, so
/// nothing can pass by assuming a square or an even dimension.
const AWKWARD_SHAPES: [(usize, usize); 5] = [(1, 1), (1, 7), (7, 1), (5, 3), (23, 4)];

/// A centred difference of a linear field is exact, and a centred average of a
/// linear field is exact, at any spacing — so the only error left is
/// floating-point rounding. `16·f64::EPSILON` scaled by the magnitude of the
/// values involved covers the handful of operations each point costs.
const ROUNDING_SLACK: f64 = 16.0 * f64::EPSILON;

#[test]
fn operators_are_exact_on_a_linear_field_at_any_basin_shape() {
    // A linear field has zero third and second derivative, so both the
    // difference and the average remainders vanish identically.
    let slope_x_per_m = 3.0e-6;
    let slope_y_per_m = -7.0e-6;
    let linear = |x_m: f64, y_m: f64| slope_x_per_m * x_m + slope_y_per_m * y_m;

    for (nx, ny) in AWKWARD_SHAPES {
        let grid = Grid::new(nx, ny).expect("non-zero extents");
        let spacing =
            Spacing::new(BASIN_LX_M / nx as f64, BASIN_LY_M / ny as f64).expect("positive spacing");
        let ops = CGridOperators::new(grid, spacing);
        let scale = BASIN_LX_M * slope_x_per_m.abs() + BASIN_LY_M * slope_y_per_m.abs();
        let tolerance = scale * ROUNDING_SLACK / spacing.dx_m().min(spacing.dy_m());

        let center = sample(grid, spacing, H_STAGGERING, linear);
        let mut on_u = grid.allocate(U_STAGGERING, 0.0_f64);
        ops.ddx_center_to_face(&center, &mut on_u);
        for j in 0..ny {
            for i in 1..nx {
                let got = *on_u.get(i, j).expect("in-bounds");
                assert!(
                    (got - slope_x_per_m).abs() <= tolerance,
                    "{nx}x{ny}: d/dx at face ({i}, {j}) was {got}, expected {slope_x_per_m}"
                );
            }
        }

        let mut on_v = grid.allocate(V_STAGGERING, 0.0_f64);
        ops.ddy_center_to_face(&center, &mut on_v);
        for j in 1..ny {
            for i in 0..nx {
                let got = *on_v.get(i, j).expect("in-bounds");
                assert!(
                    (got - slope_y_per_m).abs() <= tolerance,
                    "{nx}x{ny}: d/dy at face ({i}, {j}) was {got}, expected {slope_y_per_m}"
                );
            }
        }

        let u = sample(grid, spacing, U_STAGGERING, linear);
        let mut from_u = grid.allocate(H_STAGGERING, 0.0_f64);
        ops.ddx_face_to_center(&u, &mut from_u);
        let v = sample(grid, spacing, V_STAGGERING, linear);
        let mut from_v = grid.allocate(H_STAGGERING, 0.0_f64);
        ops.ddy_face_to_center(&v, &mut from_v);
        for j in 0..ny {
            for i in 0..nx {
                let dx = *from_u.get(i, j).expect("in-bounds");
                let dy = *from_v.get(i, j).expect("in-bounds");
                assert!(
                    (dx - slope_x_per_m).abs() <= tolerance,
                    "{nx}x{ny}: d/dx at center ({i}, {j}) was {dx}"
                );
                assert!(
                    (dy - slope_y_per_m).abs() <= tolerance,
                    "{nx}x{ny}: d/dy at center ({i}, {j}) was {dy}"
                );
            }
        }

        // Interpolating a linear field is exact too: the average of the two
        // neighbours is the value at their midpoint.
        let mut interpolated = grid.allocate(H_STAGGERING, 0.0_f64);
        ops.face_to_center_x(&u, &mut interpolated);
        for j in 0..ny {
            for i in 0..nx {
                let (x_m, y_m) = position_m(spacing, H_STAGGERING, i, j);
                let got = *interpolated.get(i, j).expect("in-bounds");
                assert!(
                    (got - linear(x_m, y_m)).abs() <= scale * ROUNDING_SLACK,
                    "{nx}x{ny}: interpolation at center ({i}, {j}) was {got}"
                );
            }
        }
    }
}

#[test]
fn a_constant_field_has_zero_derivative_at_any_basin_shape() {
    for (nx, ny) in AWKWARD_SHAPES {
        let grid = Grid::new(nx, ny).expect("non-zero extents");
        let spacing =
            Spacing::new(BASIN_LX_M / nx as f64, BASIN_LY_M / ny as f64).expect("positive spacing");
        let ops = CGridOperators::new(grid, spacing);

        let center = grid.allocate(H_STAGGERING, 12.5_f64);
        let mut on_u = grid.allocate(U_STAGGERING, f64::NAN);
        let mut on_v = grid.allocate(V_STAGGERING, f64::NAN);
        ops.ddx_center_to_face(&center, &mut on_u);
        ops.ddy_center_to_face(&center, &mut on_v);
        // Exactly zero: `(a − a)/Δ` is exact in floating point.
        assert!(on_u.as_slice().iter().all(|&value| value == 0.0));
        assert!(on_v.as_slice().iter().all(|&value| value == 0.0));

        let u = grid.allocate(U_STAGGERING, 12.5_f64);
        let v = grid.allocate(V_STAGGERING, 12.5_f64);
        let mut from_u = grid.allocate(H_STAGGERING, f64::NAN);
        let mut from_v = grid.allocate(H_STAGGERING, f64::NAN);
        ops.ddx_face_to_center(&u, &mut from_u);
        ops.ddy_face_to_center(&v, &mut from_v);
        assert!(from_u.as_slice().iter().all(|&value| value == 0.0));
        assert!(from_v.as_slice().iter().all(|&value| value == 0.0));

        // And a constant interpolates to itself, everywhere it is defined.
        let mut interpolated = grid.allocate(H_STAGGERING, f64::NAN);
        ops.face_to_center_y(&v, &mut interpolated);
        assert!(interpolated.as_slice().iter().all(|&value| value == 12.5));
    }
}

#[test]
fn operators_reuse_the_output_buffer_rather_than_allocating() {
    // Time stepping must not allocate per step (CODING_STANDARDS.md), so the
    // operators write through a caller-owned buffer. Running twice into the
    // same buffer must overwrite it completely, not accumulate.
    let (grid, spacing, ops) = basin(8);
    let center = sample(grid, spacing, H_STAGGERING, f);
    let expected = sample(grid, spacing, U_STAGGERING, dfdx);

    let mut face = grid.allocate(U_STAGGERING, 0.0_f64);
    ops.ddx_center_to_face(&center, &mut face);
    let first = face.as_slice().to_vec();
    ops.ddx_center_to_face(&center, &mut face);
    assert_eq!(face.as_slice(), first.as_slice());

    // Sanity: the buffer really did receive the derivative.
    assert!(max_error(&face, &expected, interior_x_faces(8)) < expected.as_slice()[1].abs());
}

// --- Boundaries, shapes and invalid input. ---

#[test]
fn center_to_face_operators_leave_the_basin_boundary_faces_at_zero() {
    // A center→face stencil needs a cell on each side, which the outermost
    // faces do not have. Those faces are the closed basin's walls, where
    // Epic 04 sets the boundary condition; the operators write zero there
    // rather than inventing a one-sided value.
    let (grid, spacing, ops) = basin(8);
    let center = sample(grid, spacing, H_STAGGERING, f);

    let mut on_u = grid.allocate(U_STAGGERING, f64::NAN);
    ops.ddx_center_to_face(&center, &mut on_u);
    let mut interpolated_u = grid.allocate(U_STAGGERING, f64::NAN);
    ops.center_to_face_x(&center, &mut interpolated_u);
    for j in 0..grid.ny() {
        assert_eq!(on_u.get(0, j), Some(&0.0));
        assert_eq!(on_u.get(grid.nx(), j), Some(&0.0));
        assert_eq!(interpolated_u.get(0, j), Some(&0.0));
        assert_eq!(interpolated_u.get(grid.nx(), j), Some(&0.0));
    }

    let mut on_v = grid.allocate(V_STAGGERING, f64::NAN);
    ops.ddy_center_to_face(&center, &mut on_v);
    let mut interpolated_v = grid.allocate(V_STAGGERING, f64::NAN);
    ops.center_to_face_y(&center, &mut interpolated_v);
    for i in 0..grid.nx() {
        assert_eq!(on_v.get(i, 0), Some(&0.0));
        assert_eq!(on_v.get(i, grid.ny()), Some(&0.0));
        assert_eq!(interpolated_v.get(i, 0), Some(&0.0));
        assert_eq!(interpolated_v.get(i, grid.ny()), Some(&0.0));
    }
}

#[test]
fn spacing_rejects_non_positive_and_non_finite_cell_widths() {
    // A cell width is scenario input, so a bad one is a Result naming the
    // offending value and the bound it violated (CODING_STANDARDS.md).
    assert_eq!(
        Spacing::new(0.0, 1.0),
        Err(SpacingError::NotPositive {
            axis: Axis::X,
            value_m: 0.0
        })
    );
    assert!(matches!(
        Spacing::new(1.0, -3.0),
        Err(SpacingError::NotPositive { axis: Axis::Y, .. })
    ));
    assert!(Spacing::new(f64::NAN, 1.0).is_err());
    assert!(Spacing::new(1.0, f64::INFINITY).is_err());

    let message = Spacing::new(0.0, 1.0)
        .expect_err("zero is rejected")
        .to_string();
    assert!(message.contains("dx is 0"), "{message}");

    let spacing = Spacing::new(2.0, 5.0).expect("positive widths");
    assert_eq!((spacing.dx_m(), spacing.dy_m()), (2.0, 5.0));
}

#[test]
#[should_panic(expected = "shape")]
fn an_output_buffer_of_the_wrong_shape_is_a_bug_and_panics() {
    // A mis-shaped buffer means the calling code is wrong, not that the user
    // asked for something impossible — so it panics rather than returning.
    let (grid, spacing, ops) = basin(8);
    let center = sample(grid, spacing, H_STAGGERING, f);
    let mut wrong = grid.allocate(H_STAGGERING, 0.0_f64);
    ops.ddx_center_to_face(&center, &mut wrong);
}

#[test]
#[should_panic(expected = "shape")]
fn an_input_field_of_the_wrong_shape_is_a_bug_and_panics() {
    let (grid, spacing, ops) = basin(8);
    let not_a_center_field = sample(grid, spacing, U_STAGGERING, f);
    let mut face = grid.allocate(U_STAGGERING, 0.0_f64);
    ops.ddx_center_to_face(&not_a_center_field, &mut face);
}

// --- The face↔face averages the Coriolis term needs (T-02.2). ---

/// Leading remainder of the four-point average that moves a value between the
/// two face staggerings. It is a centred average along x composed with one
/// along y, and for `f = sin(kx·x)·sin(ky·y)` both remainders are
/// proportional to `−f`, so they add rather than partly cancelling.
fn four_point_average_remainder(spacing: Spacing) -> f64 {
    average_remainder_x(spacing) + average_remainder_y(spacing)
}

#[test]
fn face_to_face_interpolation_converges_at_second_order() {
    assert_second_order(|n| {
        let (grid, spacing, ops) = basin(n);

        let expected_u = sample(grid, spacing, U_STAGGERING, f);
        let v = sample(grid, spacing, V_STAGGERING, f);
        let mut on_u = grid.allocate(U_STAGGERING, 0.0_f64);
        ops.face_y_to_face_x(&v, &mut on_u);

        let expected_v = sample(grid, spacing, V_STAGGERING, f);
        let u = sample(grid, spacing, U_STAGGERING, f);
        let mut on_v = grid.allocate(V_STAGGERING, 0.0_f64);
        ops.face_x_to_face_y(&u, &mut on_v);

        max_error(&on_u, &expected_u, interior_x_faces(n)).max(max_error(
            &on_v,
            &expected_v,
            interior_y_faces(n),
        ))
    });
}

#[test]
fn every_face_to_face_error_matches_the_analytic_truncation_remainder() {
    assert_error_matches_remainder(
        "face_y_to_face_x",
        V_STAGGERING,
        U_STAGGERING,
        f,
        f,
        four_point_average_remainder,
        CGridOperators::face_y_to_face_x,
    );
    assert_error_matches_remainder(
        "face_x_to_face_y",
        U_STAGGERING,
        V_STAGGERING,
        f,
        f,
        four_point_average_remainder,
        CGridOperators::face_x_to_face_y,
    );
}

#[test]
fn the_face_to_face_averages_reset_the_boundary_lines_they_do_not_compute() {
    // Same contract as the center→face operators: a boundary face has a cell
    // on one side only, so it carries zero rather than a one-sided guess, and
    // a reused buffer cannot leak a previous value into it.
    let (grid, _spacing, ops) = basin(4);

    let mut on_u = grid.allocate(U_STAGGERING, f64::NAN);
    let v = grid.allocate(V_STAGGERING, 1.0_f64);
    ops.face_y_to_face_x(&v, &mut on_u);
    for j in 0..on_u.ny() {
        assert_eq!(on_u.get(0, j), Some(&0.0), "western boundary at j = {j}");
        assert_eq!(
            on_u.get(grid.nx(), j),
            Some(&0.0),
            "eastern boundary at j = {j}"
        );
        for i in 1..grid.nx() {
            // Averaging a constant field must reproduce the constant.
            assert_eq!(on_u.get(i, j), Some(&1.0), "interior face at ({i}, {j})");
        }
    }

    let mut on_v = grid.allocate(V_STAGGERING, f64::NAN);
    let u = grid.allocate(U_STAGGERING, 1.0_f64);
    ops.face_x_to_face_y(&u, &mut on_v);
    for i in 0..on_v.nx() {
        assert_eq!(on_v.get(i, 0), Some(&0.0), "southern boundary at i = {i}");
        assert_eq!(
            on_v.get(i, grid.ny()),
            Some(&0.0),
            "northern boundary at i = {i}"
        );
        for j in 1..grid.ny() {
            assert_eq!(on_v.get(i, j), Some(&1.0), "interior face at ({i}, {j})");
        }
    }
}
