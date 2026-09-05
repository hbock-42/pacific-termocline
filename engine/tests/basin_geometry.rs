//! T-04.1 — the basin's extent is scenario input, in degrees of longitude and
//! latitude, and it is what places the grid in metres.
//!
//! The acceptance criterion is one sentence: *changing the config's basin
//! bounds changes the resulting grid's physical extent, verified by checking
//! physical (not just index) coordinates of grid corners.* So every assertion
//! below is about a corner's position in metres — never about `nx`, `ny` or a
//! cell index alone.
//!
//! Expected values come from the projection the equatorial beta-plane is
//! defined on, not from running the loader: the model linearizes about the
//! equator, so a degree of longitude and a degree of latitude are both one
//! degree of arc on a sphere of Earth's mean radius,
//!
//! ```text
//! R = 6 371 008.8 m          (IUGG mean radius R₁ = (2a + b)/3)
//! 1° of arc = R·π/180 = 111 194.926 644… m
//! ```
//!
//! and that is the number this file computes its expectations from. It is the
//! same `R` the model's `β = 2Ω/R = 2.3×10⁻¹¹ m⁻¹s⁻¹` is quoted from
//! (`CONTEXT.md`, *Beta-plane*), which is what keeps the geometry and the
//! rotation talking about the same planet.

use engine::basin::{BasinBounds, BasinBoundsError};
use engine::scenario::{Scenario, ScenarioError};
use termocline_grid::{Staggering, H_STAGGERING, U_STAGGERING, V_STAGGERING};

/// Earth's mean radius, in metres: IUGG R₁ = (2a + b)/3 for the WGS-84
/// ellipsoid. Written here independently of the engine's constant so that a
/// change to that constant fails these tests rather than silently moving the
/// basin.
const EARTH_MEAN_RADIUS_M: f64 = 6_371_008.8;

/// One degree of arc at Earth's mean radius, in metres.
fn metres_per_degree() -> f64 {
    EARTH_MEAN_RADIUS_M * std::f64::consts::PI / 180.0
}

/// Relative tolerance for a position compared against `R·Δ°·π/180` computed
/// here. Both sides are the same two or three `f64` multiplications in a
/// different order, so the gap is a handful of ulp; 1e-12 is roughly 4500 ulp
/// at basin scale (1e7 m) and still under 10 µm, which is nine orders of
/// magnitude finer than the 55 km cell it positions.
const POSITION_RELATIVE_TOLERANCE: f64 = 1e-12;

fn assert_close(actual_m: f64, expected_m: f64, what: &str) {
    let tolerance = expected_m.abs() * POSITION_RELATIVE_TOLERANCE;
    assert!(
        (actual_m - expected_m).abs() <= tolerance,
        "{what}: expected {expected_m} m, got {actual_m} m (tolerance {tolerance} m)"
    );
}

/// The `[basin]` section of a scenario stating every bound, wrapped in the
/// rest of a scenario that is known good so that only the basin varies.
fn scenario_with_basin(basin_section: &str) -> Result<Scenario, ScenarioError> {
    Scenario::from_toml(&format!(
        "{basin_section}
[physics]
reduced_gravity_m_per_s2 = 0.06
mean_thermocline_depth_m = 150.0
rayleigh_damping_per_s = 1.0e-7

[run]
dt_s = 3600.0
total_steps = 24
output_every_n_steps = 24
"
    ))
}

fn basin_of(basin_section: &str) -> engine::Basin {
    scenario_with_basin(basin_section)
        .unwrap_or_else(|error| panic!("this scenario should load: {error}"))
        .basin()
}

fn error_from(basin_section: &str) -> ScenarioError {
    scenario_with_basin(basin_section).expect_err("this basin should be rejected")
}

// ---------------------------------------------------------------------------
// The bounds place the grid's corners in metres.
// ---------------------------------------------------------------------------

/// The four corners of a basin, in metres: the walls themselves, addressed
/// through the C-grid faces that sit on them.
///
/// `u` lives on east/west faces indexed `0..=nx`, so column `0` is the western
/// wall and column `nx` the eastern one; `v` lives on north/south faces
/// indexed `0..=ny`, so row `0` is the southern wall and row `ny` the northern
/// one (`termocline_grid::Staggering`).
fn walls_m(basin: engine::Basin) -> (f64, f64, f64, f64) {
    (
        basin.x_of_column_m(U_STAGGERING, 0),
        basin.x_of_column_m(U_STAGGERING, basin.grid().nx()),
        basin.y_of_row_m(V_STAGGERING, 0),
        basin.y_of_row_m(V_STAGGERING, basin.grid().ny()),
    )
}

#[test]
fn a_basin_stated_in_degrees_puts_its_corners_where_the_projection_says() {
    // 120°E–160°E by 10°S–10°N at half a degree: 40° of longitude and 20° of
    // latitude, both a whole number of cells.
    let basin = basin_of(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = -10.0
northern_latitude_deg = 10.0
resolution_deg = 0.5
",
    );

    let degree_m = metres_per_degree();
    let (west_m, east_m, south_m, north_m) = walls_m(basin);
    // x is measured east from the western wall, so the western wall is the
    // origin and the eastern wall is the basin's zonal extent.
    assert_eq!(west_m, 0.0, "the western wall is the origin of x");
    assert_close(east_m, 40.0 * degree_m, "the eastern wall");
    // y is measured north from the equator, because `f = β·y` needs it to be
    // (`CONTEXT.md`, *Beta-plane*).
    assert_close(south_m, -10.0 * degree_m, "the southern wall");
    assert_close(north_m, 10.0 * degree_m, "the northern wall");
}

#[test]
fn changing_the_bounds_changes_the_physical_extent_of_the_grid() {
    let narrow = basin_of(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = -10.0
northern_latitude_deg = 10.0
resolution_deg = 0.5
",
    );
    // The same resolution over twice the longitude and twice the latitude.
    let wide = basin_of(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 200.0
southern_latitude_deg = -20.0
northern_latitude_deg = 20.0
resolution_deg = 0.5
",
    );

    let degree_m = metres_per_degree();
    let (_, narrow_east_m, narrow_south_m, narrow_north_m) = walls_m(narrow);
    let (_, wide_east_m, wide_south_m, wide_north_m) = walls_m(wide);

    assert_close(
        narrow_east_m,
        40.0 * degree_m,
        "the narrow basin's east wall",
    );
    assert_close(wide_east_m, 80.0 * degree_m, "the wide basin's east wall");
    assert_close(
        wide_south_m,
        -20.0 * degree_m,
        "the wide basin's south wall",
    );
    assert_close(wide_north_m, 20.0 * degree_m, "the wide basin's north wall");
    assert!(
        wide_east_m > narrow_east_m && wide_north_m > narrow_north_m,
        "widening the bounds has to widen the basin, not just renumber its cells"
    );
    assert!(narrow_south_m > wide_south_m);
}

#[test]
fn the_resolution_sets_the_cell_size_and_leaves_the_extent_alone() {
    const BOUNDS: &str = "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = -10.0
northern_latitude_deg = 10.0
";
    let fine = basin_of(&format!("{BOUNDS}resolution_deg = 0.5\n"));
    let coarse = basin_of(&format!("{BOUNDS}resolution_deg = 1.0\n"));

    let degree_m = metres_per_degree();
    // A cell is one resolution-degree of arc on each side: the beta-plane is a
    // linearization about the equator, so a degree of longitude and a degree
    // of latitude are the same distance.
    assert_close(fine.spacing().dx_m(), 0.5 * degree_m, "the fine cell width");
    assert_close(
        fine.spacing().dy_m(),
        0.5 * degree_m,
        "the fine cell height",
    );
    assert_close(coarse.spacing().dx_m(), degree_m, "the coarse cell width");
    assert_close(coarse.spacing().dy_m(), degree_m, "the coarse cell height");

    // Halving the resolution halves the cell count on each axis and leaves
    // every wall exactly where it was.
    assert_eq!(fine.grid().nx(), 2 * coarse.grid().nx());
    assert_eq!(fine.grid().ny(), 2 * coarse.grid().ny());
    let (fine_west_m, fine_east_m, fine_south_m, fine_north_m) = walls_m(fine);
    let (coarse_west_m, coarse_east_m, coarse_south_m, coarse_north_m) = walls_m(coarse);
    assert_close(
        coarse_east_m,
        fine_east_m,
        "the eastern wall at half the resolution",
    );
    assert_close(
        coarse_north_m,
        fine_north_m,
        "the northern wall at half the resolution",
    );
    assert_close(
        coarse_south_m,
        fine_south_m,
        "the southern wall at half the resolution",
    );
    assert_eq!(coarse_west_m, fine_west_m);
}

#[test]
fn the_cell_centres_of_the_corner_cells_sit_half_a_cell_inside_the_walls() {
    // `h` lives at cell centres, so the southwest-most value of the field the
    // simulation is about is half a cell east and half a cell north of the
    // corner — index arithmetic that is only right if the staggering is
    // applied to the basin's own origin.
    let basin = basin_of(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = -10.0
northern_latitude_deg = 10.0
resolution_deg = 0.5
",
    );
    let degree_m = metres_per_degree();
    let half_cell_m = 0.25 * degree_m;

    assert_close(
        basin.x_of_column_m(H_STAGGERING, 0),
        half_cell_m,
        "the westernmost h column",
    );
    assert_close(
        basin.y_of_row_m(H_STAGGERING, 0),
        -10.0 * degree_m + half_cell_m,
        "the southernmost h row",
    );
    assert_close(
        basin.x_of_column_m(H_STAGGERING, basin.grid().nx() - 1),
        40.0 * degree_m - half_cell_m,
        "the easternmost h column",
    );
    assert_close(
        basin.y_of_row_m(H_STAGGERING, basin.grid().ny() - 1),
        10.0 * degree_m - half_cell_m,
        "the northernmost h row",
    );
}

#[test]
fn an_equatorially_symmetric_basin_has_a_row_of_faces_exactly_on_the_equator() {
    // 10°S–10°N at 0.5° is 40 cells, so the equator falls on a north/south
    // face rather than between two of them, and `f = β·y` is exactly zero
    // there. The rotation and the geometry must agree about which row that is.
    let basin = basin_of(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = -10.0
northern_latitude_deg = 10.0
resolution_deg = 0.5
",
    );
    let equator_row = basin.grid().ny() / 2;
    assert_eq!(
        basin.y_of_row_m(V_STAGGERING, equator_row),
        0.0,
        "the equator has to be exactly y = 0, not nearly"
    );
}

#[test]
fn an_omitted_basin_section_is_the_pacific() {
    // "Sensible Pacific-basin defaults": 120°E–80°W by 25°S–25°N
    // (`CONTEXT.md`, *Basin*). A scenario that says nothing about its basin
    // gets the basin this project is about.
    let stated = basin_of(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = -80.0
southern_latitude_deg = -25.0
northern_latitude_deg = 25.0
resolution_deg = 0.5
",
    );
    let defaulted = basin_of("");

    assert_eq!(defaulted, stated, "the default basin is the Pacific");

    let degree_m = metres_per_degree();
    let (west_m, east_m, south_m, north_m) = walls_m(defaulted);
    assert_eq!(west_m, 0.0);
    // 120°E to 80°W eastward is 160° of longitude: the width of the Pacific.
    assert_close(east_m, 160.0 * degree_m, "the eastern wall of the Pacific");
    assert_close(south_m, -25.0 * degree_m, "25°S");
    assert_close(north_m, 25.0 * degree_m, "25°N");
}

#[test]
fn an_eastern_longitude_west_of_the_western_one_crosses_the_dateline() {
    // The Pacific is the basin that makes this necessary: its eastern boundary
    // at 80°W is east of its western boundary at 120°E only if longitude is
    // counted eastward across the dateline. Writing 280°E must mean the same
    // basin as writing 80°W.
    let signed = basin_of(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = -80.0
southern_latitude_deg = -25.0
northern_latitude_deg = 25.0
resolution_deg = 0.5
",
    );
    let unwrapped = basin_of(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 280.0
southern_latitude_deg = -25.0
northern_latitude_deg = 25.0
resolution_deg = 0.5
",
    );
    assert_eq!(signed, unwrapped);
}

// ---------------------------------------------------------------------------
// A basin that is not a basin is refused by name, not clamped or panicked.
// ---------------------------------------------------------------------------

#[test]
fn a_northern_boundary_south_of_the_southern_one_is_rejected_by_name() {
    let message = error_from(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = 10.0
northern_latitude_deg = -10.0
resolution_deg = 0.5
",
    )
    .to_string();
    assert!(
        message.contains("northern_latitude_deg") && message.contains("southern_latitude_deg"),
        "the message should name both boundaries, got: {message}"
    );
}

#[test]
fn a_latitude_off_the_planet_is_rejected_by_name() {
    let message = error_from(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = -95.0
northern_latitude_deg = 25.0
resolution_deg = 0.5
",
    )
    .to_string();
    assert!(
        message.contains("southern_latitude_deg") && message.contains("-95"),
        "the message should name the latitude it rejected, got: {message}"
    );
}

#[test]
fn a_non_positive_resolution_is_rejected_by_name() {
    let message = error_from(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = -10.0
northern_latitude_deg = 10.0
resolution_deg = 0.0
",
    )
    .to_string();
    assert!(
        message.contains("resolution_deg") && message.contains('0'),
        "the message should name the resolution it rejected, got: {message}"
    );
}

#[test]
fn bounds_that_are_not_a_whole_number_of_cells_are_refused_rather_than_rounded() {
    // 20° of latitude at 0.3° is 66.67 cells. Rounding it silently would run a
    // basin nobody asked for (CODING_STANDARDS.md § *No silent clamping*).
    let message = error_from(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.2
southern_latitude_deg = -10.0
northern_latitude_deg = 10.0
resolution_deg = 0.3
",
    )
    .to_string();
    assert!(
        message.contains("0.3"),
        "the message should name the resolution that does not divide the span, got: {message}"
    );
    assert!(
        message.contains("20"),
        "the message should name the span it could not divide, got: {message}"
    );
    assert!(
        message.contains("latitude"),
        "the message should say which axis is wrong, got: {message}"
    );
}

#[test]
fn a_basin_narrower_than_one_cell_is_rejected() {
    let message = error_from(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 120.0
southern_latitude_deg = -10.0
northern_latitude_deg = 10.0
resolution_deg = 0.5
",
    )
    .to_string();
    assert!(
        !message.is_empty(),
        "a zero-width basin still has to say something"
    );
}

#[test]
fn a_misspelled_basin_key_is_rejected_rather_than_ignored() {
    let message = error_from(
        "[basin]
western_longitude_deg = 120.0
eastern_longitude_deg = 160.0
southern_latitude_deg = -10.0
northern_latitude_deg = 10.0
resolution_degrees = 0.5
",
    )
    .to_string();
    assert!(
        message.contains("resolution_degrees"),
        "the message should name the key it did not recognise, got: {message}"
    );
}

// ---------------------------------------------------------------------------
// The bounds are a type, not just a section: the same validation is available
// to anything constructing a basin in code.
// ---------------------------------------------------------------------------

#[test]
fn the_pacific_bounds_are_the_bounds_context_md_states() {
    let pacific = BasinBounds::pacific();
    assert_eq!(pacific.western_longitude_deg(), 120.0);
    assert_eq!(pacific.eastern_longitude_deg(), -80.0);
    assert_eq!(pacific.southern_latitude_deg(), -25.0);
    assert_eq!(pacific.northern_latitude_deg(), 25.0);
    // 160° by 50° of arc, in cells of the default resolution.
    let basin = pacific.basin();
    assert_eq!(
        basin.grid().nx() as f64 * pacific.resolution_deg(),
        160.0,
        "the Pacific is 160° wide"
    );
    assert_eq!(
        basin.grid().ny() as f64 * pacific.resolution_deg(),
        50.0,
        "the Pacific of this model is 50° tall"
    );
}

#[test]
fn bounds_that_are_not_finite_are_refused_by_name() {
    let error = BasinBounds::new(f64::NAN, 160.0, -10.0, 10.0, 0.5)
        .expect_err("a basin cannot start at an unknown longitude");
    assert!(matches!(error, BasinBoundsError::NotFinite { .. }));
    assert!(
        error.to_string().contains("western_longitude_deg"),
        "the message should name the parameter, got: {error}"
    );
}

#[test]
fn the_extent_of_the_basin_matches_its_walls() {
    let bounds =
        BasinBounds::new(120.0, 160.0, -10.0, 10.0, 0.5).expect("these bounds are a basin");
    let basin = bounds.basin();
    let degree_m = metres_per_degree();
    assert_close(basin.zonal_extent_m(), 40.0 * degree_m, "the zonal extent");
    assert_close(
        basin.meridional_extent_m(),
        20.0 * degree_m,
        "the meridional extent",
    );
    // And the extents are the walls: the same numbers reached two ways.
    let (west_m, east_m, south_m, north_m) = walls_m(basin);
    assert_close(
        basin.zonal_extent_m(),
        east_m - west_m,
        "extent versus walls, zonally",
    );
    assert_close(
        basin.meridional_extent_m(),
        north_m - south_m,
        "extent versus walls, meridionally",
    );
    // Staggering is a position on the cell, not a shift of the basin.
    assert_eq!(
        basin.x_of_column_m(Staggering::EastWestFace, 0),
        basin.western_edge_x_m()
    );
}
