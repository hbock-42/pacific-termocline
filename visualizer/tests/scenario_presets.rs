//! T-13.3 acceptance criteria: the teaching panel's scenarios are presets
//! built in code, they say what they cost, they fit the tab, switching between
//! them starts clean — and the one that tells the El Niño story says out loud
//! what it is not.
//!
//! Nothing here is measured out of a run. The preset parameters are read off
//! `visualizer/src/presets.rs` by eye; the sizes are the arithmetic of the
//! format, one `f64` per point of each field per frame; and the two physical
//! assertions come from the linear model's own steady state, not from what
//! this code happens to produce.
//!
//! # Where the physical expectations come from
//!
//! The equilibrium of the linear shallow-water core under a steady zonal
//! stress is a thermocline tilted so that the pressure-gradient term balances
//! the wind: `g'·∂h/∂x = τx/(ρ₀·H)` (`docs/planning/01-scientific-model.md`,
//! and the tilt `docs/validation-report.md` measures for T-07.4). Two
//! consequences are used below, both of them qualitative on purpose, because
//! what is being checked is that the *captions are true* rather than that the
//! solver is accurate — the solver is Epic 07's business and is validated
//! there, natively, at half a degree.
//!
//! 1. The tilt is proportional to the stress. Weaken the wind over the basin
//!    and the tilt the ocean relaxes towards is weaker in the same proportion,
//!    so the west-minus-east depth difference shrinks and the eastern
//!    thermocline deepens.
//! 2. Under a stress that does not change, the tilt stops changing. The
//!    approach is on the Rayleigh timescale `r⁻¹` ≈ 116 days and the basin
//!    crossing time `L/c` ≈ 69 days, so a year in, a two-year run is on its
//!    equilibrium and stays there.

use engine::{ScenarioWind, WindStress};
use termocline_format::Variable;
use visualizer::{ComputedRun, FrameBudget, ScenarioPreset};

/// Steps taken per call while running a preset to the end. Large enough that a
/// two-year run is a handful of calls, and it changes nothing about the run
/// (`tests/computed_run.rs` holds chunking to that).
const CHUNK_STEPS: u64 = 4_096;

/// Seconds in a day.
const DAY_S: f64 = 86_400.0;

/// The longest a preset may take to compute natively, in seconds.
///
/// The ticket asks that a visitor "sees something within seconds". The
/// measured reference is 0.72 s for the 17 520 steps of the control preset
/// (release build, Apple M1 Pro; `visualizer/src/compute.rs`), and a preset
/// twice that long would be past what a visitor waits for on a machine some
/// unmeasured factor slower than the reference one. Two seconds is that
/// reference doubled, and it is a bound on the *whole* run: frames are drawn
/// as they arrive, so what a visitor waits for is the first of them.
const LONGEST_NATIVE_COMPUTE_S: f64 = 2.0;

/// The preset whose caption is the El Niño story, by name.
const RELAXATION: &str = "The winds relax";

/// The Kelvin wave speed of every preset's ocean, in m/s: `c = √(g'H)` with
/// `g'` = 0.06 m/s² and `H` = 150 m (`CONTEXT.md`).
const WAVE_SPEED_M_PER_S: f64 = 3.0;

/// The width of the basin every preset runs on, in metres.
///
/// 160° of longitude at the equatorial radius the engine uses for a degree of
/// arc, 111 320 m (`CONTEXT.md`, *Basin*): 1.78 × 10⁷ m.
const BASIN_WIDTH_M: f64 = 160.0 * 111_320.0;

/// The preset named `name`.
fn preset(name: &str) -> ScenarioPreset {
    ScenarioPreset::ALL
        .into_iter()
        .find(|preset| preset.name() == name)
        .unwrap_or_else(|| panic!("no preset is called {name}"))
}

/// `preset` computed to the end of its schedule.
fn computed(preset: ScenarioPreset) -> ComputedRun {
    let mut run = ComputedRun::of_preset(preset, FrameBudget::browser())
        .unwrap_or_else(|error| panic!("{} does not start: {error}", preset.name()));
    while !run.is_finished() {
        run.advance_steps(CHUNK_STEPS)
            .unwrap_or_else(|error| panic!("{} does not compute: {error}", preset.name()));
    }
    run
}

/// The thermocline depth anomaly `h`, in metres, at the western and eastern
/// ends of the equator in frame `index` of `run`.
///
/// The basin runs 25°S–25°N in 25 rows of two degrees, so the middle row is
/// centred on the equator; `h` is stored row-major on cell centres.
fn equatorial_ends_m(run: &visualizer::LoadedRun, index: u64) -> (f64, f64) {
    let grid = run.header().grid;
    let (nx, ny) = (grid.nx(), grid.ny());
    let frame = run.frame(index).expect("a frame the run holds");
    let row = (ny / 2) * nx;
    (frame.h()[row], frame.h()[row + nx - 1])
}

#[test]
fn every_preset_is_a_scenario_built_from_the_engines_own_types() {
    // The first acceptance criterion. There is no file: a preset assembles the
    // engine's `ScenarioConfig` and hands it to the same validation a file
    // would have gone through, and what comes out names the engine's own
    // forcings.
    assert!(
        ScenarioPreset::ALL.len() >= 3,
        "the ticket names three stories at minimum"
    );
    for preset in ScenarioPreset::ALL {
        let scenario = preset
            .scenario()
            .unwrap_or_else(|error| panic!("{} is not a scenario: {error}", preset.name()));
        assert!(
            !scenario.winds().is_empty(),
            "{} drives the ocean with nothing",
            preset.name()
        );
        // And the config round-trips the sections it was built from, so what
        // the engine validated is what the preset states.
        assert_eq!(preset.config().build().expect("it builds"), scenario);
    }

    // The three stories, each with the forcing that tells it.
    let control = preset("Normal trade winds").scenario().expect("it builds");
    assert!(
        matches!(control.winds(), [ScenarioWind::Steady(_)]),
        "the control is the steady trade winds and nothing else"
    );

    let relaxed = preset(RELAXATION).scenario().expect("it builds");
    let [ScenarioWind::Steady(trades), ScenarioWind::Burst(relaxation)] = relaxed.winds() else {
        panic!("the relaxation is the trades plus a westerly anomaly");
    };
    // The anomaly weakens the trades without reversing them: −0.05 Pa of
    // easterly plus +0.04 Pa of westerly leaves the equator easterly, which is
    // the observed regime.
    let (net_pa, _) = relaxed
        .wind()
        .stress(relaxation.center_x_m(), 0.0, relaxation.peak_time_s());
    assert!(
        net_pa < 0.0,
        "the relaxed trades reverse: {net_pa} Pa at the centre of the anomaly"
    );
    let (full_pa, _) = trades.stress(relaxation.center_x_m(), 0.0, 0.0);
    assert!(
        net_pa > full_pa / 2.0,
        "the trades are not appreciably weakened: {net_pa} Pa against {full_pa} Pa"
    );
    // A season rather than a gust, which is what makes it a relaxation: the
    // ocean adjusts across the basin in L/c ≈ 69 days, so an anomaly that
    // lasted days would be the burst below rather than a weakened trade wind.
    assert!(
        relaxation.duration_s() > BASIN_WIDTH_M / WAVE_SPEED_M_PER_S,
        "the winds relax for less than the ocean takes to answer"
    );

    let burst = preset("A westerly wind burst")
        .scenario()
        .expect("it builds");
    let [ScenarioWind::Steady(_), ScenarioWind::Burst(gust)] = burst.winds() else {
        panic!("the burst is the trades plus T-03.3's Gaussian anomaly");
    };
    // T-03.3's burst: Lt = 10 days, the duration of an observed one.
    assert_eq!(gust.duration_s(), 10.0 * DAY_S);

    let seasonal = preset("The year's cycle").scenario().expect("it builds");
    assert!(
        matches!(seasonal.winds(), [ScenarioWind::Seasonal(_)]),
        "the seasonal cycle is the trades modulated by the year"
    );
}

#[test]
fn every_preset_runs_its_first_steps() {
    // The rest of the first criterion: a preset that builds but does not step
    // is a button that fails when pressed. Two frames' worth of steps is
    // enough to have integrated the initial state and saved the first frames.
    for preset in ScenarioPreset::ALL {
        let mut run = ComputedRun::of_preset(preset, FrameBudget::browser())
            .unwrap_or_else(|error| panic!("{} does not start: {error}", preset.name()));
        assert_eq!(run.run().frame_count(), 0, "{}", preset.name());

        let cost = preset.cost();
        let steps_per_frame = cost.steps / (cost.frame_count - 1);
        run.advance_steps(2 * steps_per_frame)
            .unwrap_or_else(|error| panic!("{} does not step: {error}", preset.name()));

        assert!(
            run.run().frame_count() >= 2,
            "{} produced no frames in two frames' worth of steps",
            preset.name()
        );
        let frame = run.run().frame(1).expect("the second frame");
        assert!(
            frame.h().iter().all(|h_m| h_m.is_finite()),
            "{} produced a thermocline that is not a number",
            preset.name()
        );
        // The first saved frame is the initial state, which is at rest.
        let first = run.run().frame(0).expect("the first frame");
        assert!(first.h().iter().all(|h_m| *h_m == 0.0));
        assert_eq!(first.t_s(), 0.0);
    }
}

#[test]
fn every_preset_states_the_grid_duration_and_cost_the_run_turns_out_to_have() {
    // The second criterion. A preset states its grid, its length and what it
    // will cost *before* it is pressed — so what it states is held here
    // against the run it actually produces, header field by header field.
    for preset in ScenarioPreset::ALL {
        let cost = preset.cost();
        let run = ComputedRun::of_preset(preset, FrameBudget::browser())
            .unwrap_or_else(|error| panic!("{} does not start: {error}", preset.name()));
        let header = run.run().header();

        assert_eq!((cost.nx, cost.ny), (header.grid.nx(), header.grid.ny()));
        assert_eq!(cost.frame_count, header.output.frame_count);
        assert_eq!(cost.frame_interval_s, header.output.interval_s);
        assert_eq!(
            cost.frame_bytes,
            FrameBudget::bytes_of(header),
            "{} misstates what it will hold",
            preset.name()
        );
        // The duration is the schedule's own: steps times the timestep.
        let last_frame_s = cost.frame_interval_s * (cost.frame_count - 1) as f64;
        assert!(cost.duration_days * DAY_S >= last_frame_s);

        // And it is all in the line the panel prints.
        let line = cost.line();
        for expected in [
            format!("{} × {} cells", cost.nx, cost.ny),
            format!("{} frames", cost.frame_count),
        ] {
            assert!(
                line.contains(&expected),
                "{} does not state its {expected}: {line}",
                preset.name()
            );
        }
        assert!(
            line.contains(" MB") && line.contains(" s to compute"),
            "{} states neither its size nor its compute time: {line}",
            preset.name()
        );
    }
}

#[test]
fn every_preset_fits_the_frame_budget_a_tab_is_held_to() {
    // The second criterion's other half, checked rather than assumed. The
    // budget is 33.6 MB of frames (ADR-0012), enforced against the header
    // before the first step — so a preset past it is a button that refuses
    // when pressed.
    let budget = FrameBudget::browser();
    for preset in ScenarioPreset::ALL {
        let cost = preset.cost();
        assert!(
            cost.fits(budget),
            "{} holds {} bytes of frames against a budget of {}",
            preset.name(),
            cost.frame_bytes,
            budget.max_bytes()
        );
        // The estimate the panel prints is what the run is admitted on: the
        // budget is checked against the header the engine writes, and a preset
        // whose estimate disagreed would be admitted or refused on a number
        // the reader was never shown.
        ComputedRun::of_preset(preset, budget)
            .unwrap_or_else(|error| panic!("{} is refused: {error}", preset.name()));

        // A frame is one `f64` per point of each of the five variables of the
        // linear core: the arithmetic of the format, done here rather than
        // read back out of the crate.
        let scenario = preset.scenario().expect("it builds");
        let bounds = scenario.bounds();
        let per_frame: u64 = [
            bounds.nx() * bounds.ny(),
            (bounds.nx() + 1) * bounds.ny(),
            bounds.nx() * (bounds.ny() + 1),
            (bounds.nx() + 1) * bounds.ny(),
            bounds.nx() * (bounds.ny() + 1),
        ]
        .iter()
        .map(|len| *len as u64 * 8)
        .sum();
        assert_eq!(Variable::LINEAR_CORE.len(), 5);
        assert_eq!(cost.frame_bytes, per_frame * cost.frame_count);
    }
}

#[test]
fn a_visitor_waits_seconds_rather_than_minutes() {
    // The rest of the second criterion: "fast enough that a visitor sees
    // something within seconds".
    for preset in ScenarioPreset::ALL {
        let cost = preset.cost();
        assert!(
            cost.native_compute_s <= LONGEST_NATIVE_COMPUTE_S,
            "{} takes {:.1} s to compute natively, past the {LONGEST_NATIVE_COMPUTE_S} s a \
             visitor waits for",
            preset.name(),
            cost.native_compute_s
        );
        // The estimate is the measurement scaled by cells and steps, and the
        // measurement was taken on this grid — so every preset is on it, and
        // the scaling is exercised only by the step count.
        assert_eq!(
            (cost.nx, cost.ny),
            (80, 25),
            "{} is not on the grid the compute-time measurement was taken on",
            preset.name()
        );
    }
}

#[test]
fn switching_scenario_starts_a_run_with_nothing_of_the_last_one_in_it() {
    // The third criterion. A preset produces a new `ComputedRun` — a new
    // engine loop over a new run — so the frames, the header, the scales and
    // the progress of whatever was being computed are dropped rather than
    // reset. What is asserted is the observable consequence: a run started
    // after another one is frame-for-frame the run started on its own.
    let first = preset("Normal trade winds");
    let second = preset(RELAXATION);

    let mut abandoned = ComputedRun::of_preset(first, FrameBudget::browser()).expect("it starts");
    abandoned.advance_steps(CHUNK_STEPS).expect("it computes");
    assert!(abandoned.run().frame_count() > 1);

    let mut switched = ComputedRun::of_preset(second, FrameBudget::browser()).expect("it starts");
    assert_eq!(switched.run().frame_count(), 0, "the new run starts empty");
    assert_eq!(
        switched.progress(),
        (0, second.cost().frame_count),
        "the new run's progress is its own"
    );
    assert_eq!(switched.run().source(), second.name());
    switched.advance_steps(CHUNK_STEPS).expect("it computes");

    let mut alone = ComputedRun::of_preset(second, FrameBudget::browser()).expect("it starts");
    alone.advance_steps(CHUNK_STEPS).expect("it computes");

    assert_eq!(switched.run().frame_count(), alone.run().frame_count());
    assert_eq!(switched.run().anomaly_scale(), alone.run().anomaly_scale());
    for index in 0..alone.run().frame_count() {
        let after = switched.run().frame(index).expect("a frame");
        let fresh = alone.run().frame(index).expect("the same frame");
        // Byte equality: identical scenario in, identical output out
        // (CODING_STANDARDS.md § *Correctness and failure*).
        assert_eq!(after.t_s(), fresh.t_s());
        assert_eq!(after.h(), fresh.h());
        assert_eq!(after.tau_x(), fresh.tau_x());
    }
}

#[test]
fn the_relaxed_winds_preset_is_not_offered_as_the_model_producing_el_nino() {
    // The fourth criterion, and the one that matters most. The relaxation is a
    // *forced* response: the wind change is prescribed and there is no
    // feedback from the ocean to the wind, so nothing in this run is the model
    // producing a cycle of its own. `docs/enso-oscillation-report.md` records
    // that the coupled model's period is 1.03 years, outside the observed 2–7
    // year band, with that acceptance criterion explicitly unmet — a teaching
    // panel claiming the cycle would contradict the project's own validation
    // ledger.
    let relaxation = preset(RELAXATION);
    let caveat = relaxation
        .caveat()
        .expect("the El Niño story is the one preset that has something to disclaim");

    assert!(
        caveat.contains("1.03 years"),
        "the caveat does not give the period the model actually oscillates at: {caveat}"
    );
    assert!(
        caveat.contains("2 to 7 years") || caveat.contains("2–7 years"),
        "the caveat does not give the band that period misses: {caveat}"
    );
    assert!(
        caveat.contains("docs/enso-oscillation-report.md"),
        "the caveat does not say where the measurement is recorded: {caveat}"
    );
    assert!(
        caveat.contains("unmet") || caveat.contains("not met"),
        "the caveat does not say the criterion was missed: {caveat}"
    );
    assert!(
        caveat.contains("prescribed"),
        "the caveat does not say the wind change was imposed rather than produced: {caveat}"
    );

    // And nothing anywhere in the preset's own prose claims the cycle.
    for preset in ScenarioPreset::ALL {
        let prose = format!(
            "{} {} {} {}",
            preset.name(),
            preset.summary(),
            preset.story(),
            preset.caveat().unwrap_or_default()
        )
        .to_lowercase();
        for claim in [
            "reproduces el niño",
            "el niño cycle",
            "enso cycle",
            "the model oscillates",
            "emergent",
        ] {
            assert!(
                !prose.contains(claim),
                "{} claims \"{claim}\", which the model does not do",
                preset.name()
            );
        }
        // Any mention of El Niño at all is a mention that has to be qualified.
        if prose.contains("el niño") {
            assert!(
                preset
                    .caveat()
                    .is_some_and(|caveat| caveat.contains("1.03 years")),
                "{} names El Niño without saying what this model does not do",
                preset.name()
            );
        }
    }
}

#[test]
fn the_control_preset_establishes_a_tilt_and_holds_it() {
    // The control's caption says the tilt builds and then holds "because the
    // wind holds". That is a claim about the run, so it is checked: under a
    // steady stress the linear core relaxes onto the equilibrium tilt
    // `g'·∂h/∂x = τx/(ρ₀H)` on the Rayleigh timescale r⁻¹ ≈ 116 days, and a
    // year in, a two-year run is on it.
    let control = preset("Normal trade winds");
    let run = computed(control);
    let cost = control.cost();
    let frames_per_year = (365.0 * DAY_S / cost.frame_interval_s) as u64;

    let (west_m, east_m) = equatorial_ends_m(run.run(), cost.frame_count - 1);
    assert!(
        west_m > 0.0 && east_m < 0.0,
        "the easterly trades did not pile warm water in the west: {west_m} m against {east_m} m"
    );

    let tilt_of = |index| {
        let (west_m, east_m) = equatorial_ends_m(run.run(), index);
        west_m - east_m
    };
    let after_a_year = tilt_of(frames_per_year);
    let at_the_end = tilt_of(cost.frame_count - 1);
    // A tenth: the remaining approach to equilibrium after a year is
    // exp(−365/116) ≈ 4 %, and the reflected waves that ring through the
    // adjustment are smaller still by then. A tilt that had not settled would
    // be out by far more than this, and one that had drifted away would be out
    // in the other direction.
    assert!(
        (at_the_end - after_a_year).abs() / after_a_year < 0.1,
        "the tilt did not hold: {after_a_year} m a year in against {at_the_end} m at the end"
    );
}

#[test]
fn the_relaxed_winds_preset_flattens_the_tilt_and_deepens_the_east() {
    // The relaxation's caption says the tilt flattens, the west thins and the
    // east deepens. The equilibrium tilt is proportional to the stress that
    // drives it, so weakening the wind over the basin weakens the tilt the
    // ocean relaxes towards, and the warm water that was held in the west
    // moves east.
    //
    // How much: the tilt scales with the zonal integral of the stress. The
    // anomaly integrates to τ·Lx·√π ≈ 0.04 · 6 000 km · 1.7725 over the part
    // of it that lies inside the basin (96 %), which is 46 % of the trades'
    // own 0.05 · 17 800 km — so a quasi-steady ocean would flatten by about
    // 46 %. It is not quasi-steady: the anomaly lasts 90 days against a basin
    // adjustment of L/c ≈ 69 days, so only part of that is reached. A tenth is
    // asserted, far below the bound the forcing sets and far above the
    // nothing a mislabelled preset would show.
    let relaxation = preset(RELAXATION);
    let run = computed(relaxation);
    let cost = relaxation.cost();
    let frame_of_day = |day: f64| (day * DAY_S / cost.frame_interval_s) as u64;

    // Day 180: the tilt is established and the anomaly, peaking at day 365
    // with `Lt` = 90 days, is exp(−(185/90)²) = 1.5 % of itself.
    let (west_before_m, east_before_m) = equatorial_ends_m(run.run(), frame_of_day(180.0));
    // Day 420: the anomaly has been at or near its peak for a season and the
    // ocean has had a basin crossing to answer it.
    let (west_after_m, east_after_m) = equatorial_ends_m(run.run(), frame_of_day(420.0));

    let before_m = west_before_m - east_before_m;
    let after_m = west_after_m - east_after_m;
    assert!(
        before_m > 0.0,
        "the trades had not tilted the thermocline before the winds relaxed: {before_m} m"
    );
    assert!(
        after_m < 0.9 * before_m,
        "the tilt did not flatten when the winds relaxed: {before_m} m before, {after_m} m after"
    );
    assert!(
        after_m > 0.0,
        "the tilt reversed, which a weakened easterly cannot do: {after_m} m"
    );
    assert!(
        east_after_m > east_before_m,
        "the warm water did not move back east: {east_before_m} m before, {east_after_m} m after"
    );
    assert!(
        west_after_m < west_before_m,
        "the west did not give up the water the east gained: {west_before_m} m before, \
         {west_after_m} m after"
    );
}

#[test]
fn the_burst_preset_runs_long_enough_to_watch_the_wave_cross() {
    // The burst's caption says the pulse can be followed the whole way across.
    // A Kelvin wave travels at c = √(g'H) = 3.0 m/s, so it crosses the
    // 1.78 × 10⁷ m basin in 6.9 × 10⁵ s ≈ 69 days: the run has to continue
    // that long past the burst, and sample it often enough to be followed.
    let burst = preset("A westerly wind burst");
    let scenario = burst.scenario().expect("it builds");
    let [_, ScenarioWind::Burst(gust)] = scenario.winds() else {
        panic!("the burst preset carries a burst");
    };
    let cost = burst.cost();
    let crossing_s = BASIN_WIDTH_M / WAVE_SPEED_M_PER_S;

    let after_the_burst_s = cost.duration_days * DAY_S - gust.peak_time_s();
    assert!(
        after_the_burst_s > crossing_s,
        "the run ends {after_the_burst_s} s after the burst, before a wave crossing of {crossing_s} s"
    );
    // And the frames are close enough together to see it move: at least fifty
    // of them across the basin, which at 80 cells is a wave that never jumps
    // more than a couple of cells between frames.
    let frames_across = crossing_s / cost.frame_interval_s;
    assert!(
        frames_across >= 50.0,
        "the burst is sampled {frames_across} times as it crosses, which is not a wave being \
         watched"
    );
}
