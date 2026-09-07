//! T-13.5 acceptance criteria: the shell opens on the story in a browser and
//! on the instrument on a desktop, and one control moves between them without
//! either losing anything.
//!
//! Epic 13 built a second reading of the same run — the ocean from the side
//! (T-13.1), in real longitudes, metres and months (T-13.2), started from
//! scenarios that are stories rather than files (T-13.3), captioned from what
//! was measured along the equator (T-13.4). This ticket is only *which of the
//! two is on screen*, so what is asserted here is the rule and not a picture:
//! which mode a platform opens on, which of the seven Epics 08–09 views each
//! mode draws, and what the one control between them says.
//!
//! The shell's own behaviour — that the browser build opens with a scenario
//! already computing, that switching modes rebuilds nothing, that every cost
//! property still holds — is asserted in `src/app.rs`'s test module, which can
//! see the state those claims are about.
//!
//! # Where the expected values come from
//!
//! The mode rules are the ticket's, restated: the toggle reaches every view
//! Epics 08–09 built, the native default is deliberate, and the control is
//! discoverable without being noisy.
//!
//! The frame arithmetic in the last two tests is the format's, done here
//! rather than read back out of the crate. The presets run on an 80 × 25 basin
//! (`crate::presets`, `BROWSER_BASIN`: 120°E–80°W by 25°S–25°N at 2°), and a
//! frame carries one `f64` per point of each variable at that variable's
//! staggering — `h` at cell centres, `u` and `τx` on east–west faces, `v` and
//! `τy` on north–south faces (`termocline_grid::Staggering`). The SST anomaly
//! `T'` sits at cell centres with `h`.

use termocline_format::{BasinExtent, GridSpec, Variable};
use visualizer::{FrameBudget, Mode, Platform, ScenarioPreset, ScientificView};

/// Cells east–west and north–south of the basin every preset runs on.
const NX: u64 = 80;
const NY: u64 = 25;

/// Bytes an `f64` occupies in an encoded frame: the frame encoding is
/// `bincode`'s fixed-width standard configuration.
const BYTES_PER_VALUE: u64 = 8;

/// Values one frame of the linear core carries over that basin: `h` at cell
/// centres (80 × 25), `u` and `τx` on east–west faces (81 × 25), `v` and `τy`
/// on north–south faces (80 × 26).
const LINEAR_CORE_VALUES_PER_FRAME: u64 = NX * NY + 2 * (NX + 1) * NY + 2 * NX * (NY + 1);

/// Values `T'` would add to each of those frames: one per cell.
const SST_VALUES_PER_FRAME: u64 = NX * NY;

/// The basin every preset runs on, as the format describes it.
fn browser_grid() -> GridSpec {
    GridSpec::new(
        NX as usize,
        NY as usize,
        BasinExtent::new(120.0, -80.0, -25.0, 25.0),
    )
    .expect("the presets' basin is a grid the format can describe")
}

/// The preset whose frames cost the most: the wind burst buys a daily frame so
/// the pulse can be followed across the basin (`crate::presets`).
///
/// Found by the name on its own button rather than by its position in the
/// list, so that reordering the presets moves this test with them instead of
/// silently retargeting it at a cheaper scenario.
fn wind_burst() -> ScenarioPreset {
    *ScenarioPreset::ALL
        .iter()
        .find(|preset| preset.name() == "A westerly wind burst")
        .expect("the teaching panel offers the wind burst")
}

#[test]
fn a_browser_opens_on_the_story_and_a_desktop_on_the_instrument() {
    assert_eq!(Mode::opening_on(Platform::Browser), Mode::Teaching);
    assert_eq!(Mode::opening_on(Platform::Desktop), Mode::Scientific);
}

/// The native binary is the desktop one, so the rule above is the rule this
/// build actually gets. Compiled only where it is true: the other half of the
/// rule is asserted above against the [`Platform`] value itself, because a
/// test compiled for the host can only ever answer for the host.
#[test]
#[cfg(not(target_arch = "wasm32"))]
fn the_native_build_is_a_desktop_and_so_opens_on_the_instrument() {
    assert_eq!(Platform::of_target(), Platform::Desktop);
    assert_eq!(
        Mode::opening_on(Platform::of_target()),
        Mode::Scientific,
        "a scientist opening the desktop app is not shown the teaching view first"
    );
}

/// Nothing built so far is removed: every view has a mode that draws it, and
/// the one control reaches that mode from wherever the reader is.
#[test]
fn the_toggle_reaches_every_view_epics_08_and_09_built() {
    for view in ScientificView::ALL {
        assert!(
            [Mode::Teaching, Mode::Scientific]
                .into_iter()
                .any(|mode| mode.draws(view)),
            "no mode draws {}, so the toggle cannot reach it",
            view.label()
        );
        assert!(
            Mode::Scientific.draws(view),
            "and the mode that draws all of them is the scientific one, so one \
             switch reaches every view rather than some of them each: {}",
            view.label()
        );
    }
}

/// And the seven are seven: a view listed twice would let the test above pass
/// while a view of Epics 08-09 went unnamed.
#[test]
fn the_seven_views_are_seven_different_views() {
    for (index, view) in ScientificView::ALL.into_iter().enumerate() {
        for other in ScientificView::ALL.into_iter().skip(index + 1) {
            assert_ne!(view, other, "{} is listed twice", view.label());
            assert_ne!(view.label(), other.label(), "two views share a name");
        }
    }
}

/// The two the story keeps are the two that are not readings of the ocean but
/// ways of moving through time; the other five say in a colour, a chart or a
/// number what the captions say in words.
#[test]
fn the_story_keeps_only_the_ways_of_moving_through_time() {
    for view in ScientificView::ALL {
        let kept = matches!(view, ScientificView::Scrubber | ScientificView::Playback);
        assert_eq!(
            Mode::Teaching.draws(view),
            kept,
            "the teaching mode's answer for {} is wrong",
            view.label()
        );
    }
}

#[test]
fn one_control_moves_between_the_two_and_back() {
    assert_eq!(Mode::Teaching.other(), Mode::Scientific);
    assert_eq!(Mode::Scientific.other(), Mode::Teaching);
    for mode in [Mode::Teaching, Mode::Scientific] {
        assert_eq!(mode.other().other(), mode);
    }
}

/// Discoverable: the control names where it goes, so a reader looking for the
/// instrument reads the word *scientific* on it rather than having to work out
/// which half of a switch they are standing on.
#[test]
fn the_control_names_where_it_goes_rather_than_where_it_is() {
    assert!(
        Mode::Teaching.switch_label().contains("Scientific"),
        "from the story, the control offers the instrument: {}",
        Mode::Teaching.switch_label()
    );
    assert!(
        Mode::Scientific.switch_label().contains("Teaching"),
        "from the instrument, the control offers the story: {}",
        Mode::Scientific.switch_label()
    );
}

/// Not noisy: the list of what is on the other side is in the hover text
/// rather than on screen, and it names the views so the offer is legible.
#[test]
fn the_control_says_what_is_on_the_other_side_of_it() {
    let offered = Mode::Teaching.switch_hover().to_lowercase();
    for fragment in [
        "basin map",
        "wind stress",
        "cross-section",
        "time series",
        "side by side",
    ] {
        assert!(
            offered.contains(fragment),
            "the offer of the instrument does not mention {fragment}: {offered}"
        );
    }
    let offered = Mode::Scientific.switch_hover().to_lowercase();
    for fragment in ["side", "months", "captions"] {
        assert!(
            offered.contains(fragment),
            "the offer of the story does not mention {fragment}: {offered}"
        );
    }
}

/// The teaching scenarios are the validated linear ocean of Epics 01–07, with
/// no `[sst]` section and so no `T'`.
///
/// The alternative was tempting: El Niño is a warm event and sea-surface
/// temperature is what a general audience has intuition for. It was not taken,
/// and the test below is the half of the reason that is arithmetic. The other
/// half is that every teaching view Epic 13 built — the side view, the coast
/// depths, all four caption topics — is a function of `h`, so a coupled run
/// would put a second headline quantity on screen that nothing in the teaching
/// mode can describe.
#[test]
fn every_teaching_preset_is_an_uncoupled_run() {
    for preset in ScenarioPreset::ALL {
        assert!(
            preset.config().sst.is_none(),
            "{} carries an [sst] section",
            preset.name()
        );
    }
}

/// And the arithmetic half: `T'` would put the most expensive preset over the
/// budget ADR-0012 holds a tab to, so "make the teaching scenarios coupled" is
/// not a change of one field.
#[test]
fn coupling_the_wind_burst_would_break_the_budget_a_tab_is_held_to() {
    let cost = wind_burst().cost();
    let budget = FrameBudget::browser().max_bytes();
    assert_eq!(
        cost.frame_bytes,
        cost.frame_count * LINEAR_CORE_VALUES_PER_FRAME * BYTES_PER_VALUE,
        "the preset's stated size is the format's arithmetic over its own frames"
    );
    assert!(
        cost.frame_bytes <= budget,
        "the uncoupled wind burst fits: {} of {budget}",
        cost.frame_bytes
    );
    assert_eq!(
        browser_grid().field_len(Variable::SstAnomaly) as u64,
        SST_VALUES_PER_FRAME,
        "T' sits with h at cell centres, so it costs one f64 per cell per frame"
    );
    let coupled_bytes =
        cost.frame_bytes + cost.frame_count * SST_VALUES_PER_FRAME * BYTES_PER_VALUE;
    assert!(
        coupled_bytes > budget,
        "a coupled wind burst would cost {coupled_bytes} against a budget of {budget}, so it \
         would be refused before its first step"
    );
}
