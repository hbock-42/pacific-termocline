//! T-09.1 acceptance criteria: the wind-stress overlay shows easterly arrows
//! along the equator for the steady trade-wind scenario, and toggling it on
//! and off does not affect the heatmap underneath.
//!
//! Like the heatmap tests, these assert the layer instead of looking at it. A
//! [`WindOverlay`] is a list of arrows in the map's own cell coordinates, built
//! without a window or a device, so "arrows point west along the equator" is a
//! sign check on a tip rather than a screenshot — and the same code draws the
//! overlay in a browser and natively (ADR-0006).
//!
//! # Where the expected values come from
//!
//! The stress field under test is the control scenario's, written out from
//! `engine/scenarios/steady-trades.toml` and the `SteadyTradeWinds` definition
//! in `docs/planning/01-scientific-model.md`:
//!
//! ```text
//! τx(x, y) = τ₀ · exp(−(y / Ly)²),   τy = 0
//! ```
//!
//! with `τ₀ = −0.05 Pa` and `Ly = 3.61 × 10⁵ m`, the equatorial deformation
//! radius `Le = √(c/β)` (`CONTEXT.md`). That formula is the *input* these tests
//! draw and also what they check the sampled arrows against; it is never read
//! back off the visualizer's own output. The sign is what the criterion turns
//! on: easterly trade-wind stress is `τx < 0` (`CONTEXT.md`, *Wind stress*), and
//! an easterly stress pushes the ocean westward, so its arrow points west.

mod common;

use common::{encoded_frames_with_fields, steady_trades_header, FrameFields, NX, NY, PACIFIC};
use termocline_format::{RunHeader, Variable};
use visualizer::{
    Heatmap, LoadedRun, RunBytes, WindOverlay, ARROW_SPACING_CELLS, MAX_ARROW_LENGTH_CELLS,
    MIN_ARROW_LENGTH_CELLS,
};

/// `τ₀`, the equatorial zonal stress of `steady-trades.toml`, in pascals.
/// Negative because the alizés are easterly.
const EQUATORIAL_ZONAL_STRESS_PA: f64 = -0.05;

/// `Ly`, the meridional decay scale of `steady-trades.toml`, in metres: the
/// equatorial deformation radius `Le = √(c/β)` for `c = 3.0 m s⁻¹` and
/// `β = 2.3 × 10⁻¹¹ m⁻¹ s⁻¹`.
const MERIDIONAL_DECAY_SCALE_M: f64 = 361_000.0;

/// Metres per degree of arc on a sphere of the Earth's mean radius,
/// 6 371 008.8 m (IUGG mean radius `R₁`). The visualizer must not link the
/// engine (ADR-0001), so the projection the scenario's `Ly` is expressed
/// against is restated here rather than imported.
const METRES_PER_DEGREE_OF_ARC: f64 = 6_371_008.8 * std::f64::consts::PI / 180.0;

/// A meridional decay scale wide enough that the waveguide covers most of the
/// basin, in metres: ten times the control scenario's.
///
/// The shape of the forcing off the equator is what
/// [`the_arrows_shorten_where_the_stress_weakens`] is about, and under the
/// control scenario's own `Ly` the stress a few degrees off the equator is
/// already too weak to draw at all. A wider waveguide is a scenario the engine
/// takes — the decay scale is scenario input — and it puts several arrow rows
/// on the falling part of the same Gaussian.
const WIDE_DECAY_SCALE_M: f64 = 10.0 * MERIDIONAL_DECAY_SCALE_M;

/// A relative tolerance of a few ULP.
///
/// The overlay reports `τx` at a cell centre as the mean of the two east/west
/// faces the cell sits between. Where the field does not vary with `x` — which
/// is every trade-wind scenario — those two values are equal and their mean is
/// exact but for the rounding of one addition and one halving. Four ULP is the
/// next round number clear of that; nothing here is a physical tolerance.
const FEW_ULP: f64 = 4.0 * f64::EPSILON;

/// The trade-wind `τx` of `steady-trades.toml` at latitude `latitude_deg`.
///
/// The scenario's analytic profile, evaluated independently of the code under
/// test.
fn trade_wind_stress_pa(latitude_deg: f64, decay_scale_m: f64) -> f64 {
    let y_m = latitude_deg * METRES_PER_DEGREE_OF_ARC;
    let scaled = y_m / decay_scale_m;
    EQUATORIAL_ZONAL_STRESS_PA * (-scaled * scaled).exp()
}

/// The strongest trade-wind stress anywhere on `header`'s grid, in pascals.
///
/// Not `τ₀`: the equator falls on a cell *edge* of this basin, so the two rows
/// either side of it sit a quarter degree off and carry very slightly less than
/// the equatorial stress. The profile says how much less.
fn strongest_trade_wind_stress_pa(header: &RunHeader, decay_scale_m: f64) -> f64 {
    (0..header.grid.ny())
        .map(|j| trade_wind_stress_pa(cell_centre_latitude_deg(header, j), decay_scale_m).abs())
        .fold(0.0_f64, f64::max)
}

/// The latitude of the centre of cell row `j`, in degrees north.
fn cell_centre_latitude_deg(header: &RunHeader, j: usize) -> f64 {
    let extent = header.grid.extent();
    #[allow(clippy::cast_precision_loss)]
    let resolution_deg =
        (extent.north_deg_north - extent.south_deg_north) / header.grid.ny() as f64;
    #[allow(clippy::cast_precision_loss)]
    let offset = j as f64 + 0.5;
    extent.south_deg_north + offset * resolution_deg
}

/// Trade winds of decay scale `decay_scale_m` sampled onto `header`'s C-grid:
/// `τx` at the east/west faces, `τy` calm.
///
/// The scale is a parameter because it is scenario input — exactly the knob the
/// forcing sensitivity of Epic 07 varies — not because anything here needs it
/// to vary for its own convenience.
///
/// `τx` does not vary with `x`, so every face of a row carries the row's own
/// value; the row's `y` is the cell centre's, which is where `u` — and so `τx`
/// — sits on the meridional axis (`Staggering::EastWestFace`).
fn trade_wind_fields(header: &RunHeader, decay_scale_m: f64) -> FrameFields {
    let (faces_x, faces_y) = header
        .grid
        .grid()
        .field_shape(Variable::ZonalWindStress.staggering());
    let mut tau_x_pa = Vec::with_capacity(faces_x * faces_y);
    for j in 0..faces_y {
        let stress_pa = trade_wind_stress_pa(cell_centre_latitude_deg(header, j), decay_scale_m);
        tau_x_pa.extend(std::iter::repeat_n(stress_pa, faces_x));
    }
    FrameFields {
        tau_x_pa,
        ..FrameFields::calm(header)
    }
}

/// A uniform easterly stress of `tau_x_pa` over the whole basin.
///
/// Not a scenario — the trades decay away from the equator — but the field that
/// asks a layout question without a magnitude question mixed in: every cell
/// carries the same stress, so every cell the layout picks gets an arrow.
fn uniform_stress_fields(header: &RunHeader, tau_x_pa: f64) -> FrameFields {
    let calm = FrameFields::calm(header);
    FrameFields {
        tau_x_pa: vec![tau_x_pa; calm.tau_x_pa.len()],
        ..calm
    }
}

/// The image row of the arrow line closest to the equator.
fn nearest_arrow_row(header: &RunHeader) -> usize {
    arrow_line(header.grid.ny())
        .min_by(|a, b| {
            latitude_of_row(header, *a as f64)
                .abs()
                .total_cmp(&latitude_of_row(header, *b as f64).abs())
        })
        .expect("the basin carries arrows")
}

/// A one-frame run of `header`'s shape whose only frame carries `fields`.
fn run_of_one_frame(header: &RunHeader, fields: FrameFields) -> LoadedRun {
    let bytes = RunBytes {
        header: serde_json::to_vec(header).expect("a header serializes"),
        frames: encoded_frames_with_fields(header, 1, |_| FrameFields {
            h_m: fields.h_m.clone(),
            tau_x_pa: fields.tau_x_pa.clone(),
            tau_y_pa: fields.tau_y_pa.clone(),
        }),
    };
    LoadedRun::from_bytes("run-steady-trades", bytes).expect("the run loads")
}

/// The overlay of the only frame of a one-frame run carrying `fields`.
fn overlay_of(header: &RunHeader, fields: FrameFields) -> WindOverlay {
    let run = run_of_one_frame(header, fields);
    let frame = run.frame(0).expect("a one-frame run has a frame 0");
    WindOverlay::of_frame(run.header().grid, &frame, run.wind_stress_scale())
        .expect("the frame fits its own grid")
}

#[test]
fn the_steady_trades_are_drawn_as_easterly_arrows_along_the_equator() {
    // The acceptance criterion itself. Every arrow inside the waveguide the
    // scenario forces must carry an easterly stress and be drawn pointing
    // west: on the image, west is the −x direction.
    let header = steady_trades_header(1);
    let overlay = overlay_of(
        &header,
        trade_wind_fields(&header, MERIDIONAL_DECAY_SCALE_M),
    );

    let waveguide_deg = MERIDIONAL_DECAY_SCALE_M / METRES_PER_DEGREE_OF_ARC;
    let mut equatorial = 0;
    for arrow in overlay.arrows() {
        let (tail_x, tail_y) = arrow.tail_cells();
        let latitude_deg = latitude_of_row(&header, tail_y);
        if latitude_deg.abs() > waveguide_deg {
            continue;
        }
        equatorial += 1;
        assert!(
            arrow.tau_x_pa() < 0.0,
            "the alizés are easterly, so τx < 0 at {latitude_deg}°: {}",
            arrow.tau_x_pa()
        );
        let (tip_x, tip_y) = arrow.tip_cells();
        assert!(
            tip_x < tail_x,
            "an easterly stress at {latitude_deg}° must be drawn pointing west: \
             tail {tail_x}, tip {tip_x}"
        );
        // τy is calm in this scenario, so the arrow is purely zonal.
        assert!(
            (tip_y - tail_y).abs() < FEW_ULP * MAX_ARROW_LENGTH_CELLS,
            "τy is zero here, so the arrow must not tilt: tail {tail_y}, tip {tip_y}"
        );
    }
    assert!(
        equatorial >= arrow_line(NX).count(),
        "a whole row of the lattice falls inside the waveguide; only {equatorial} arrows did"
    );
}

#[test]
fn each_arrow_carries_the_stress_the_scenario_puts_at_its_latitude() {
    // Pointing the right way is not enough: an overlay that read the wrong row,
    // or averaged across the wrong pair of faces, would still point west
    // everywhere. The analytic profile says what the value must be.
    let header = steady_trades_header(1);
    let overlay = overlay_of(
        &header,
        trade_wind_fields(&header, MERIDIONAL_DECAY_SCALE_M),
    );

    for arrow in overlay.arrows() {
        let (_, tail_y) = arrow.tail_cells();
        let latitude_deg = latitude_of_row(&header, tail_y);
        let expected_pa = trade_wind_stress_pa(latitude_deg, MERIDIONAL_DECAY_SCALE_M);
        assert!(
            (arrow.tau_x_pa() - expected_pa).abs() <= FEW_ULP * expected_pa.abs(),
            "τx at {latitude_deg}°: expected {expected_pa}, got {}",
            arrow.tau_x_pa()
        );
        assert_eq!(
            arrow.tau_y_pa(),
            0.0,
            "the trades have no meridional stress"
        );
    }
}

#[test]
fn the_arrows_shorten_where_the_stress_weakens() {
    // The Gaussian is the whole meridional structure of the forcing. An overlay
    // that drew every arrow the same length would hide it, and hiding it is
    // hiding why the response is trapped near the equator.
    let header = steady_trades_header(1);
    let overlay = overlay_of(&header, trade_wind_fields(&header, WIDE_DECAY_SCALE_M));

    let mut by_latitude: Vec<(f64, f64)> = overlay
        .arrows()
        .iter()
        // One meridional profile: the arrows of the westernmost column.
        .filter(|arrow| arrow.tail_cells().0 == overlay.arrows()[0].tail_cells().0)
        .map(|arrow| {
            (
                latitude_of_row(&header, arrow.tail_cells().1),
                arrow.length_fraction(),
            )
        })
        .collect();
    by_latitude.sort_by(|a, b| a.0.abs().total_cmp(&b.0.abs()));
    assert!(
        by_latitude.len() >= 4,
        "a profile of {} arrows says nothing about a shape",
        by_latitude.len()
    );

    for pair in by_latitude.windows(2) {
        let ((near_deg, near), (far_deg, far)) = (pair[0], pair[1]);
        assert!(
            near >= far,
            "the stress falls away from the equator, so the arrow at {near_deg}° \
             cannot be shorter than the one at {far_deg}°: {near} vs {far}"
        );
    }

    // And the lengths are the profile itself, not merely decreasing: each is
    // the stress the scenario puts at that latitude over the run's scale.
    let scale_pa = overlay.scale().max_magnitude_pa();
    for (latitude_deg, fraction) in by_latitude {
        let expected = trade_wind_stress_pa(latitude_deg, WIDE_DECAY_SCALE_M).abs() / scale_pa;
        assert!(
            (fraction - expected).abs() <= FEW_ULP,
            "the arrow at {latitude_deg}° is {fraction} of the scale; the profile says {expected}"
        );
    }
}

#[test]
fn a_northward_stress_is_drawn_pointing_up_the_image() {
    // `y` increases northward in the field and downward in the image, exactly
    // as it does for the heatmap. An overlay that forgot the flip would draw a
    // southerly wind as a northerly one, which reads as perfectly plausible.
    let header = steady_trades_header(1);
    let mut fields = FrameFields::calm(&header);
    fields.tau_y_pa = vec![EQUATORIAL_ZONAL_STRESS_PA.abs(); fields.tau_y_pa.len()];
    let overlay = overlay_of(&header, fields);

    for arrow in overlay.arrows() {
        let (tail_x, tail_y) = arrow.tail_cells();
        let (tip_x, tip_y) = arrow.tip_cells();
        assert!(
            tip_y < tail_y,
            "a northward stress belongs pointing up: tail {tail_y}, tip {tip_y}"
        );
        assert!(
            (tip_x - tail_x).abs() < FEW_ULP * MAX_ARROW_LENGTH_CELLS,
            "τx is zero here, so the arrow must not tilt: tail {tail_x}, tip {tip_x}"
        );
    }
}

#[test]
fn an_arrow_reports_the_stress_at_the_cell_centre_between_the_faces_it_spans() {
    // `τx` lives on east/west faces and `τy` on north/south faces (ADR-0003),
    // so neither is defined where an arrow is drawn. The overlay averages each
    // to the cell centre; reading one face instead would put the arrow half a
    // cell from the stress it draws.
    let header = steady_trades_header(1);
    let mut fields = FrameFields::calm(&header);
    let (faces_x, faces_y) = header
        .grid
        .grid()
        .field_shape(Variable::ZonalWindStress.staggering());
    // Linear in `i`, so the mean of two neighbouring faces is the value at the
    // centre between them and nothing else is.
    for j in 0..faces_y {
        for i in 0..faces_x {
            #[allow(clippy::cast_precision_loss)]
            let value = i as f64;
            fields.tau_x_pa[j * faces_x + i] = -value;
        }
    }
    let overlay = overlay_of(&header, fields);

    for arrow in overlay.arrows() {
        let (tail_x, _) = arrow.tail_cells();
        // The centre of cell `i` is at `i + 0.5` cells from the western edge,
        // and the two faces around it carry −i and −(i + 1).
        let expected = -tail_x;
        assert!(
            (arrow.tau_x_pa() - expected).abs() <= FEW_ULP * expected.abs(),
            "at x = {tail_x} cells τx should be {expected}, got {}",
            arrow.tau_x_pa()
        );
    }
}

#[test]
fn the_arrows_sit_at_cell_centres_spaced_across_the_whole_basin() {
    // The overlay is drawn over the map in the map's own coordinates, so an
    // arrow's position is a cell-centre position: an arrow that landed on a
    // cell corner would be drawn half a cell from the stress it reports. The
    // stress is uniform here so that every cell the layout picks carries an
    // arrow, and the layout is the only thing under test.
    let header = steady_trades_header(1);
    let overlay = overlay_of(
        &header,
        uniform_stress_fields(&header, EQUATORIAL_ZONAL_STRESS_PA),
    );

    let columns: Vec<usize> = arrow_line(NX).collect();
    let rows: Vec<usize> = arrow_line(NY).collect();
    assert_eq!(overlay.arrows().len(), columns.len() * rows.len());
    for arrow in overlay.arrows() {
        let (x_cells, y_cells) = arrow.tail_cells();
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let (column, row) = (x_cells.floor() as usize, y_cells.floor() as usize);
        assert!(
            columns.contains(&column) && rows.contains(&row),
            "an arrow at ({x_cells}, {y_cells}) is off the layout's lattice"
        );
        assert!(
            (x_cells.fract() - 0.5).abs() < FEW_ULP && (y_cells.fract() - 0.5).abs() < FEW_ULP,
            "({x_cells}, {y_cells}) is not a cell centre"
        );
    }
}

#[test]
fn a_row_of_arrows_is_drawn_along_the_equator_itself() {
    // Where the arrows fall is what decides whether the criterion can be read
    // off the screen at all. The basin is symmetric about the equator, and the
    // equator is where the trades are strongest and where the response they
    // drive lives, so a layout that left it between two rows would draw
    // everything except the thing being looked for.
    let header = steady_trades_header(1);
    let overlay = overlay_of(
        &header,
        trade_wind_fields(&header, MERIDIONAL_DECAY_SCALE_M),
    );

    // Half a cell is as close as a cell centre gets: the equator falls on a
    // cell edge of this basin, not on a centre.
    #[allow(clippy::cast_precision_loss)]
    let half_cell_deg = (PACIFIC.north_deg_north - PACIFIC.south_deg_north) / (2.0 * NY as f64);
    let on_the_equator: Vec<&visualizer::WindArrow> = overlay
        .arrows()
        .iter()
        .filter(|arrow| latitude_of_row(&header, arrow.tail_cells().1).abs() <= half_cell_deg)
        .collect();
    assert_eq!(
        on_the_equator.len(),
        arrow_line(NX).count(),
        "the equatorial row must carry an arrow in every column of the lattice"
    );
    for arrow in on_the_equator {
        assert!(
            arrow.tau_x_pa() < 0.0,
            "the equatorial arrows are the easterly ones: {}",
            arrow.tau_x_pa()
        );
    }
}

#[test]
fn a_stress_too_weak_to_draw_an_arrow_is_not_drawn_as_a_dot() {
    // The trades fall off as a Gaussian, so most of the basin carries a stress
    // orders of magnitude below the strongest in the run. Drawn, those arrows
    // are shorter than the cell they sit on — a stipple of marks over the map
    // standing for a wind that is not there.
    let header = steady_trades_header(1);
    // Short enough that its arrow cannot reach the minimum length, by the
    // overlay's own two constants.
    let too_weak_pa =
        EQUATORIAL_ZONAL_STRESS_PA * MIN_ARROW_LENGTH_CELLS / MAX_ARROW_LENGTH_CELLS / 2.0;
    let mut fields = uniform_stress_fields(&header, too_weak_pa);
    let (faces_x, _) = header
        .grid
        .grid()
        .field_shape(Variable::ZonalWindStress.staggering());
    // One band at full strength, so the run has a scale to be weak against.
    let strong_row = NY - 1 - nearest_arrow_row(&header);
    for i in 0..faces_x {
        fields.tau_x_pa[strong_row * faces_x + i] = EQUATORIAL_ZONAL_STRESS_PA;
    }
    let overlay = overlay_of(&header, fields);

    assert_eq!(
        overlay.arrows().len(),
        arrow_line(NX).count(),
        "only the band at full strength is long enough to draw"
    );
    for arrow in overlay.arrows() {
        assert!(
            arrow.length_cells() >= MIN_ARROW_LENGTH_CELLS,
            "an arrow shorter than a cell must not be drawn: {} cells",
            arrow.length_cells()
        );
    }
}

#[test]
fn one_stress_scale_covers_the_whole_run_so_a_weakening_wind_is_seen_to_weaken() {
    // The same reason the colour scale is the run's: rescaling per frame would
    // redraw every frame at full length, and a wind burst relaxing back to the
    // trades — the thing this overlay exists to show beside its effect — would
    // look like a wind that never moved.
    let header = steady_trades_header(2);
    let strong = trade_wind_fields(&header, MERIDIONAL_DECAY_SCALE_M);
    let weak: Vec<f64> = strong.tau_x_pa.iter().map(|tau| tau / 4.0).collect();
    let bytes = RunBytes {
        header: serde_json::to_vec(&header).expect("a header serializes"),
        frames: encoded_frames_with_fields(&header, 2, |index| FrameFields {
            tau_x_pa: if index == 0 {
                strong.tau_x_pa.clone()
            } else {
                weak.clone()
            },
            ..FrameFields::calm(&header)
        }),
    };
    let run = LoadedRun::from_bytes("run", bytes).expect("the run loads");
    let strongest_pa = strongest_trade_wind_stress_pa(&header, MERIDIONAL_DECAY_SCALE_M);
    assert!(
        (run.wind_stress_scale().max_magnitude_pa() - strongest_pa).abs() <= FEW_ULP * strongest_pa,
        "the run's scale is the strongest stress anywhere in it, {strongest_pa} Pa: {}",
        run.wind_stress_scale().max_magnitude_pa()
    );

    let overlay_of_frame = |index| {
        let frame = run.frame(index).expect("the run has two frames");
        WindOverlay::of_frame(run.header().grid, &frame, run.wind_stress_scale())
            .expect("the frame fits its own grid")
    };
    let longest = |overlay: &WindOverlay| {
        overlay
            .arrows()
            .iter()
            .map(visualizer::WindArrow::length_fraction)
            .fold(0.0_f64, f64::max)
    };
    // A quarter of the stress is a quarter of the length: the reader compares
    // arrows between frames by looking at them. The two frames are drawn
    // against the one scale, so the ratio of the lengths is the ratio of the
    // stresses and nothing else.
    let (strong_cells, weak_cells) = (longest(&overlay_of_frame(0)), longest(&overlay_of_frame(1)));
    assert!(
        (weak_cells - strong_cells / 4.0).abs() <= FEW_ULP,
        "a quarter-strength wind must be drawn a quarter as long: {strong_cells} then {weak_cells}"
    );
    // And the lengths are the stresses against the run's scale, not against
    // each frame's own strongest arrow — which would make both frames equal.
    let strongest_drawn_pa = trade_wind_stress_pa(
        nearest_arrow_latitude_deg(&header),
        MERIDIONAL_DECAY_SCALE_M,
    )
    .abs();
    assert!(
        (strong_cells - strongest_drawn_pa / strongest_pa).abs() <= FEW_ULP,
        "the longest arrow of the strong frame is its own stress over the run's scale"
    );
}

#[test]
fn a_calm_ocean_is_drawn_with_no_arrows_at_all() {
    // A zero-length arrow is a dot, and a scatter of dots over the map reads as
    // data. Nothing to draw is drawn as nothing.
    let header = steady_trades_header(1);
    let overlay = overlay_of(&header, FrameFields::calm(&header));
    assert_eq!(overlay.scale().max_magnitude_pa(), 0.0);
    assert!(overlay.arrows().is_empty());
}

#[test]
fn a_stress_that_has_blown_up_is_drawn_as_missing_rather_than_as_a_gale() {
    // A non-finite `τ` has no direction to point in. Drawing it at full length
    // would put the strongest wind in the basin wherever the arithmetic broke,
    // and rescale every honest arrow against it.
    let header = steady_trades_header(1);
    let clean = overlay_of(
        &header,
        trade_wind_fields(&header, MERIDIONAL_DECAY_SCALE_M),
    );
    let mut fields = trade_wind_fields(&header, MERIDIONAL_DECAY_SCALE_M);
    let (faces_x, _) = header
        .grid
        .grid()
        .field_shape(Variable::ZonalWindStress.staggering());
    // The western face of the cell the northwesternmost arrow samples, so that
    // the broken value is one an arrow would otherwise have drawn.
    // The western face of a cell on the equatorial arrow row: the broken value
    // is one an arrow would otherwise have drawn.
    let column = arrow_line(NX).next().expect("the basin carries arrows");
    let broken = (NY - 1 - nearest_arrow_row(&header)) * faces_x + column;
    fields.tau_x_pa[broken] = f64::NAN;
    let overlay = overlay_of(&header, fields);

    assert_eq!(
        overlay.scale().max_magnitude_pa(),
        clean.scale().max_magnitude_pa(),
        "the scale ignores the non-finite value rather than being destroyed by it"
    );
    assert_eq!(
        overlay.arrows().len(),
        clean.arrows().len() - 1,
        "the one arrow that would have drawn the broken stress is not drawn"
    );
    for arrow in overlay.arrows() {
        assert!(
            arrow.tau_x_pa().is_finite() && arrow.tau_y_pa().is_finite(),
            "a non-finite stress must not reach a drawn arrow"
        );
    }
}

#[test]
fn toggling_the_overlay_does_not_change_the_heatmap_under_it() {
    // The second acceptance criterion. The overlay is a separate layer built
    // from the same frame: it takes the frame by reference and returns
    // geometry, so there is no path by which showing or hiding it could reach
    // the colour-mapped image the shell uploads.
    let header = steady_trades_header(1);
    let mut fields = trade_wind_fields(&header, MERIDIONAL_DECAY_SCALE_M);
    fields.h_m = tilt_field(NX, NY);
    let run = run_of_one_frame(&header, fields);
    let frame = run.frame(0).expect("a one-frame run has a frame 0");

    let without = Heatmap::of_frame(run.header().grid, &frame, run.anomaly_scale())
        .expect("the frame fits its own grid");
    let overlay = WindOverlay::of_frame(run.header().grid, &frame, run.wind_stress_scale())
        .expect("the frame fits its own grid");
    let with = Heatmap::of_frame(run.header().grid, &frame, run.anomaly_scale())
        .expect("the frame fits its own grid");

    assert!(
        !overlay.arrows().is_empty(),
        "the overlay must have something to hide for this to mean anything"
    );
    assert_eq!(
        without.rgb(),
        with.rgb(),
        "the map is the same image whether or not the overlay is drawn over it"
    );
    assert_eq!(without.scale(), with.scale());
}

/// The latitude of the arrow row closest to the equator, in degrees north.
///
/// Which rows carry an arrow is a property of the overlay's layout rather than
/// of the stress field, so it is restated here from that rule.
fn nearest_arrow_latitude_deg(header: &RunHeader) -> f64 {
    arrow_line(header.grid.ny())
        .map(|row| cell_centre_latitude_deg(header, header.grid.ny() - 1 - row))
        .min_by(|a, b| a.abs().total_cmp(&b.abs()))
        .expect("the basin is taller than half a spacing")
}

/// The latitude of the cell row an arrow at `y_cells` down the image sits in.
///
/// Row 0 of the image is the northernmost, and the field's row 0 is the
/// southernmost — the same flip the heatmap makes.
fn latitude_of_row(header: &RunHeader, y_cells: f64) -> f64 {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let row = y_cells.floor() as usize;
    cell_centre_latitude_deg(header, header.grid.ny() - 1 - row)
}

/// The lines of cells that carry arrows along an axis `extent` cells long.
///
/// The overlay anchors its pattern on the middle of each axis rather than on
/// its start, so that a basin laid out symmetrically about the equator gets a
/// row of arrows along the equator itself.
fn arrow_line(extent: usize) -> impl Iterator<Item = usize> {
    (extent / 2 % ARROW_SPACING_CELLS..extent).step_by(ARROW_SPACING_CELLS)
}

/// A thermocline tilt to draw under the arrows: deep in the west, shallow in
/// the east, linear between. T-07.4's equilibrium of this scenario, as in
/// `tests/heatmap.rs`.
fn tilt_field(nx: usize, ny: usize) -> Vec<f64> {
    use common::{EASTERN_WALL_H_M, WESTERN_WALL_H_M};
    let mut field = Vec::with_capacity(nx * ny);
    for _ in 0..ny {
        for i in 0..nx {
            #[allow(clippy::cast_precision_loss)]
            let fraction = i as f64 / (nx - 1) as f64;
            field.push(WESTERN_WALL_H_M + (EASTERN_WALL_H_M - WESTERN_WALL_H_M) * fraction);
        }
    }
    field
}
