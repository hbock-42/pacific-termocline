//! The scenarios the teaching mode offers, built in code rather than read from
//! files.
//!
//! A visitor picking "the trade winds relax" is not choosing a configuration;
//! they are choosing a story to watch. So a preset is not a path to a TOML
//! file, and there is no TOML file: each one assembles the engine's own
//! [`ScenarioConfig`] out of [`BasinSection`], [`PhysicsSection`],
//! [`RunSection`] and the [`WindSection`] entries that name
//! `SteadyTradeWinds`, `SeasonalTradeWinds` and `WindBurstAnomaly`, and hands
//! it to the same [`ScenarioConfig::build`] a file would have gone through.
//! Everything the engine checks about a scenario is still checked; what is
//! gone is the file, and the reader's need to know one exists.
//!
//! # A preset says what it costs
//!
//! ADR-0012 makes memory the binding constraint on a run in a tab, and
//! [`crate::FrameBudget`] enforces it — but enforcement is a refusal, and a
//! refusal is the wrong place for a visitor to learn how big a scenario is. So
//! every preset states its grid, its length and what it will cost *before* it
//! is pressed ([`PresetCost`]). Three of the four are facts about the run and
//! are held against the run itself by `tests/scenario_presets.rs` — the grid,
//! the frame count and the bytes those frames occupy, checked against the
//! header the engine writes. The fourth, the compute time, is an *estimate*:
//! one measurement scaled by cells and steps
//! ([`NATIVE_COST_S_PER_CELL_STEP`]), bounded by a test rather than timed by
//! one, because a wall clock in a test suite measures the machine CI happened
//! to run on. It is labelled as an estimate on screen for the same reason.
//! The budget is then a backstop rather than the first thing a reader hears
//! about.
//!
//! # A preset says what it is not
//!
//! One of these is the El Niño story — the alizés weaken and the warm water
//! slides back east — and it would be very easy, and quite wrong, to caption
//! it as this model reproducing ENSO. It is not: the wind change is
//! *prescribed*, and the run is the validated linear ocean of Epics 01–07
//! answering it. The project's own coupled experiment, which does oscillate on
//! its own, oscillates with a period of 1.03 years, outside the observed 2–7
//! year band, and `docs/enso-oscillation-report.md` records that acceptance
//! criterion as not met. A teaching panel that quietly claimed otherwise would
//! contradict the project's own validation ledger, so [`ScenarioPreset::caveat`]
//! carries that sentence to the screen and `tests/scenario_presets.rs` checks
//! that no preset's prose claims the cycle.

use engine::{
    BasinSection, PhysicsSection, RunSection, Scenario, ScenarioConfig, ScenarioError, WindSection,
};
use termocline_format::{GridSpec, OutputTiming};

use crate::compute::FrameBudget;
use crate::run::SECONDS_PER_DAY as DAY_S;
use crate::RunClock;

/// Zonal stress `τ₀` the steady trade winds put on the equator, in Pa.
///
/// −0.05 Pa is the mean easterly stress over the equatorial Pacific
/// (`docs/planning/01-scientific-model.md`, *Wind stress*), and the value
/// every scenario in `engine/scenarios/` is written with.
const TRADE_WIND_STRESS_PA: f64 = -0.05;

/// Meridional decay scale `Ly` of the trade winds, in metres.
///
/// The equatorial deformation radius `Le = √(c/β)` = 3.61 × 10⁵ m for
/// `c` = 3.0 m/s and `β` = 2.3 × 10⁻¹¹ m⁻¹s⁻¹ (`CONTEXT.md`): the width of the
/// waveguide, and so of the wind that drives it.
const WAVEGUIDE_SCALE_M: f64 = 361_000.0;

/// Reduced gravity `g'`, in m/s².
///
/// With [`MEAN_THERMOCLINE_DEPTH_M`] it gives `c = √(g'H)` = 3.0 m/s, the
/// observed first-baroclinic Kelvin speed of the equatorial Pacific
/// (`CONTEXT.md`).
const REDUCED_GRAVITY_M_PER_S2: f64 = 0.06;

/// Mean thermocline depth `H`, in metres.
const MEAN_THERMOCLINE_DEPTH_M: f64 = 150.0;

/// Rayleigh damping coefficient `r`, in s⁻¹: `r⁻¹` ≈ 116 days, the weak
/// damping of `docs/planning/01-scientific-model.md`.
const RAYLEIGH_DAMPING_PER_S: f64 = 1.0e-7;

/// Solver timestep, in seconds.
///
/// One hour, far inside the CFL bound `0.8·dx/c` ≈ 59 200 s of the browser
/// grid below.
const DT_S: f64 = 3600.0;

/// Solver steps in a day of model time.
const STEPS_PER_DAY: u64 = 24;

/// The basin every preset runs on: the equatorial Pacific of `CONTEXT.md`,
/// *Basin* — 120°E–80°W by 25°S–25°N — in cells of two degrees, which is
/// 80 × 25.
///
/// Two degrees rather than the half-degree of `engine/scenarios/`: ADR-0012
/// makes memory the constraint, and a half-degree frame is 1.29 MB against
/// this grid's 82 kB. A 2° cell is ~222 km of arc against a deformation radius
/// of 361 km, so the waveguide these scenarios are about is still resolved —
/// but a *validated* run is a native run of the engine's own scenarios
/// (`docs/validation-report.md`), never one of these.
const BROWSER_BASIN: BasinSection = BasinSection {
    western_longitude_deg: 120.0,
    eastern_longitude_deg: -80.0,
    southern_latitude_deg: -25.0,
    northern_latitude_deg: 25.0,
    resolution_deg: 2.0,
};

/// The ocean every preset runs in.
const BROWSER_PHYSICS: PhysicsSection = PhysicsSection {
    reduced_gravity_m_per_s2: REDUCED_GRAVITY_M_PER_S2,
    mean_thermocline_depth_m: MEAN_THERMOCLINE_DEPTH_M,
    rayleigh_damping_per_s: RAYLEIGH_DAMPING_PER_S,
    beta_per_m_per_s: None,
    reference_density_kg_per_m3: None,
};

/// Two years of model time, saved every three days: 730 days in 17 520 hourly
/// steps, and 244 frames of them.
///
/// The length the slow stories need — a tilt that establishes and holds, a
/// wind that relaxes over a season, a year that comes round twice — and the
/// cadence that keeps them inside the frame budget: 244 frames is 19.9 MB of
/// the 33.6 MB a tab is allowed.
const TWO_YEARS_EVERY_THREE_DAYS: RunSection = RunSection {
    dt_s: DT_S,
    total_steps: 730 * STEPS_PER_DAY,
    output_every_n_steps: 3 * STEPS_PER_DAY,
};

/// One year of model time, saved every day: 365 days in 8 760 hourly steps,
/// and 366 frames of them.
///
/// The burst is the one story whose object *moves*, and a frame every three
/// days would sample a basin crossing 23 times. A daily frame samples it 69
/// times, at 29.9 MB of the 33.6 MB budget — the one preset that spends its
/// budget on cadence rather than on length.
const ONE_YEAR_DAILY: RunSection = RunSection {
    dt_s: DT_S,
    total_steps: 365 * STEPS_PER_DAY,
    output_every_n_steps: STEPS_PER_DAY,
};

/// The steady alizés, as every preset carries them.
const STEADY_TRADES: WindSection = WindSection::SteadyTradeWinds {
    equatorial_zonal_stress_pa: TRADE_WIND_STRESS_PA,
    meridional_decay_scale_m: Some(WAVEGUIDE_SCALE_M),
};

/// The day the wind changes, in the two presets where it does.
///
/// A year into a two-year run and half a year into a one-year run: in both
/// cases long enough after the start for the wind-driven tilt to have
/// established (it settles on the Rayleigh timescale `r⁻¹` ≈ 116 days), and
/// long enough before the end for the answer to be watched.
const RELAXATION_PEAK_DAY: f64 = 365.0;
/// The same, for the burst's one-year run.
const BURST_PEAK_DAY: f64 = 180.0;

/// Peak westerly stress of the relaxation anomaly, in Pa.
///
/// +0.04 Pa against −0.05 Pa of alizés: at the centre of the patch, at its
/// peak, 80 % of the easterly stress is gone and the remaining 20 % still
/// blows east to west. The trades weaken; they do not reverse, which is the
/// observed regime (`docs/planning/01-scientific-model.md`, *Wind burst*).
const RELAXATION_STRESS_PA: f64 = 0.04;

/// Zonal `e`-folding scale of the relaxation patch, in metres.
///
/// 6 000 km of a basin 17 800 km wide, centred on it: broad enough that what a
/// reader sees is the alizés giving way across the middle of the Pacific
/// rather than a local gust, and narrow enough to leave the eastern end of the
/// basin under its own wind.
const RELAXATION_ZONAL_SCALE_M: f64 = 6_000_000.0;

/// Where that patch is centred, in metres east of the western wall: the middle
/// of a basin 160° ≈ 17 800 km wide.
const RELAXATION_CENTER_X_M: f64 = 8_900_000.0;

/// Temporal `e`-folding scale of the relaxation, in seconds: 90 days.
///
/// A season, not a gust. It is what separates this preset from the burst
/// below, which is the same Gaussian with `Lt` = 10 days: the ocean has time
/// to come most of the way to the equilibrium the weakened wind implies, which
/// is why the tilt is seen to flatten rather than merely to wobble.
const RELAXATION_DURATION_S: f64 = 90.0 * DAY_S;

/// Peak westerly stress of the wind burst, in Pa (T-03.3).
const BURST_STRESS_PA: f64 = 0.04;

/// Zonal `e`-folding scale of the burst, in metres (T-03.3).
const BURST_ZONAL_SCALE_M: f64 = 1_000_000.0;

/// Where the burst sits, in metres east of the western wall: 2 000 km, about
/// 138°E, the warm pool where observed westerly bursts sit.
const BURST_CENTER_X_M: f64 = 2_000_000.0;

/// Temporal `e`-folding scale of the burst, in seconds: `Lt` = 10 days, the
/// duration of an observed burst (T-03.3).
const BURST_DURATION_S: f64 = 10.0 * DAY_S;

/// Relative amplitude of the seasonal harmonic, dimensionless.
///
/// The alizés vary by roughly ±20 % over the year in the equatorial Pacific
/// (`docs/planning/01-scientific-model.md`, *Seasonal cycle*).
const SEASONAL_RELATIVE_AMPLITUDE: f64 = 0.2;

/// The day the alizés are strongest: 210 days into a run started on 1 January,
/// the boreal summer.
const SEASONAL_PEAK_DAY: f64 = 210.0;

/// Seconds one cell costs for one step, natively.
///
/// Measured, not modelled: 17 520 steps of the 80 × 25 browser grid take
/// 0.72 s in a release build on an Apple M1 Pro — the 41 µs a step
/// [`crate::compute`] states — and 0.72 / (17 520 · 2 000) = 2.05 × 10⁻⁸ s.
/// Scaling it by cells and steps assumes the solver's cost is linear in both,
/// which is what a fixed-stencil explicit scheme over a rectangular grid is;
/// every preset here is on the grid the measurement was taken on, so the
/// scaling is exercised only by the step count.
///
/// It is a *native* figure. A browser is slower by a factor nothing in this
/// project has profiled, which is why [`PresetCost::line`] says "natively" out
/// loud rather than quoting a number for a machine it has not measured.
const NATIVE_COST_S_PER_CELL_STEP: f64 = 2.05e-8;

/// What one preset will cost to compute and to hold.
///
/// Every field is arithmetic on the preset's own sections — no run is started
/// to produce it, which is the point: a visitor is told what a scenario costs
/// before they press it. `tests/scenario_presets.rs` holds each field against
/// the run the preset actually produces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PresetCost {
    /// Cells east–west.
    pub nx: usize,
    /// Cells north–south.
    pub ny: usize,
    /// Solver steps the run takes.
    pub steps: u64,
    /// Model time the run integrates, in days.
    ///
    /// Up to one frame interval longer than the span of the frames it saves,
    /// which is what [`PresetCost::line`] states and what a reader can scrub
    /// through: the last step is saved only when the cadence divides it. The
    /// two are held within an interval of each other by
    /// `tests/scenario_presets.rs`.
    pub duration_days: f64,
    /// Frames it saves.
    pub frame_count: u64,
    /// Model time between two frames, in seconds.
    pub frame_interval_s: f64,
    /// Bytes of encoded frames it will hold.
    pub frame_bytes: u64,
    /// Wall-clock seconds it takes to compute, natively.
    pub native_compute_s: f64,
}

impl PresetCost {
    /// Whether the frames of this run fit in `budget`.
    #[must_use]
    pub const fn fits(&self, budget: FrameBudget) -> bool {
        self.frame_bytes <= budget.max_bytes()
    }

    /// The cadence of the run, as the run's own header will declare it.
    #[must_use]
    pub const fn timing(&self) -> OutputTiming {
        OutputTiming {
            frame_count: self.frame_count,
            interval_s: self.frame_interval_s,
        }
    }

    /// The whole cost in one line: `80 × 25 cells · 729 days (2 years) · 244
    /// frames, 19.9 MB of the 33.6 MB a tab holds · about 0.7 s to compute
    /// natively`.
    ///
    /// The three things the ticket asks a preset to state — grid, duration,
    /// expected compute time — plus the memory, because that is the limit a
    /// tab actually dies on. The duration is [`RunClock`]'s phrase over the
    /// *frames*, which is the span a reader can scrub through and the same
    /// phrase the caption under the ocean uses, rather than a second way of
    /// writing a run length here. The compute time says "natively" out loud:
    /// it is [`NATIVE_COST_S_PER_CELL_STEP`] scaled, and a browser is slower
    /// by a factor nothing here has measured.
    #[must_use]
    pub fn line(&self) -> String {
        format!(
            "{nx} × {ny} cells · {span} · {frames} frames, {held} of the {budget} a tab holds · \
             about {compute:.1} s to compute natively",
            nx = self.nx,
            ny = self.ny,
            span = RunClock::of_run(self.timing()).span_phrase(),
            frames = self.frame_count,
            held = crate::InMegabytes(self.frame_bytes),
            budget = crate::InMegabytes(FrameBudget::browser().max_bytes()),
            compute = self.native_compute_s,
        )
    }
}

/// One of the scenarios the teaching panel offers as a button.
///
/// The prose and the physics travel together deliberately: `story` is what the
/// run will show and the sections below are what produces it, so a preset
/// whose caption stopped matching its wind is one file to notice it in rather
/// than two.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScenarioPreset {
    /// What the button says.
    name: &'static str,
    /// One line for the button's tooltip: what this run shows.
    summary: &'static str,
    /// The paragraph under the picker: what to watch for, in plain words.
    story: &'static str,
    /// What this run is *not*, where a reader could reasonably conclude
    /// something the model does not support. `None` where there is nothing to
    /// disclaim.
    caveat: Option<&'static str>,
    /// The basin, the ocean, the schedule and the winds — the four sections a
    /// scenario file would have carried.
    basin: BasinSection,
    /// The constants of the ocean.
    physics: PhysicsSection,
    /// How long the run is and how often it is saved.
    run: RunSection,
    /// The wind forcings, in the order they are summed.
    winds: &'static [WindSection],
}

impl ScenarioPreset {
    /// The presets the teaching panel offers, in the order it offers them.
    ///
    /// The first three are the ticket's own minimum — the control, the
    /// relaxation, the burst — in the order a reader should meet them: the
    /// tilt has to exist before its collapse means anything. The seasonal
    /// cycle follows as the one that needs no event at all.
    pub const ALL: [Self; 4] = [
        Self {
            name: "Normal trade winds",
            summary: "The control: steady winds from the east, and the tilt they hold",
            story: "The trade winds blow steadily from east to west along the equator. They push \
                    warm surface water into the west, so the warm layer grows deep off Indonesia \
                    and thin off South America. Watch the tilt build over the first few months \
                    and then hold, because the wind holds.",
            caveat: None,
            basin: BROWSER_BASIN,
            physics: BROWSER_PHYSICS,
            run: TWO_YEARS_EVERY_THREE_DAYS,
            winds: &[STEADY_TRADES],
        },
        Self {
            name: "The winds relax",
            summary: "The easterly winds fade for a season, and the warm water slides back east",
            story: "A year in, the easterly winds fade over the middle of the Pacific and stay \
                    weak for a season. With less wind holding it there, the warm water that was \
                    piled up in the west slides back east: the tilt flattens, the warm layer \
                    thins off Indonesia and deepens off South America. This is the ocean half of \
                    the El Niño story.",
            caveat: Some(
                "This is the ocean answering a wind change we imposed by hand. Nothing here is \
                 the model producing El Niño on its own: the wind is prescribed, so there is no \
                 feedback from the ocean back to it. The coupled version of this model, which \
                 does oscillate by itself, swings with a period near 1.03 years — well short of \
                 the 2 to 7 years real El Niño events take, an acceptance criterion this project \
                 records as unmet in docs/enso-oscillation-report.md.",
            ),
            basin: BROWSER_BASIN,
            physics: BROWSER_PHYSICS,
            run: TWO_YEARS_EVERY_THREE_DAYS,
            winds: &[
                STEADY_TRADES,
                WindSection::WindBurstAnomaly {
                    peak_zonal_stress_pa: RELAXATION_STRESS_PA,
                    center_x_m: RELAXATION_CENTER_X_M,
                    zonal_scale_m: RELAXATION_ZONAL_SCALE_M,
                    meridional_scale_m: WAVEGUIDE_SCALE_M,
                    peak_time_s: RELAXATION_PEAK_DAY * DAY_S,
                    duration_s: RELAXATION_DURATION_S,
                },
            ],
        },
        Self {
            name: "A westerly wind burst",
            summary: "Ten days of wind the wrong way, and the pulse it sends east",
            story: "Half a year in, a patch of wind over the warm pool blows the other way — west \
                    to east — for about ten days. The pulse it launches travels east along the \
                    equator at about 3 metres a second and takes some ten weeks to cross the \
                    ocean. Frames here are a day apart so it can be followed the whole way.",
            caveat: None,
            basin: BROWSER_BASIN,
            physics: BROWSER_PHYSICS,
            run: ONE_YEAR_DAILY,
            winds: &[
                STEADY_TRADES,
                WindSection::WindBurstAnomaly {
                    peak_zonal_stress_pa: BURST_STRESS_PA,
                    center_x_m: BURST_CENTER_X_M,
                    zonal_scale_m: BURST_ZONAL_SCALE_M,
                    meridional_scale_m: WAVEGUIDE_SCALE_M,
                    peak_time_s: BURST_PEAK_DAY * DAY_S,
                    duration_s: BURST_DURATION_S,
                },
            ],
        },
        Self {
            name: "The year's cycle",
            summary: "The same winds breathing with the seasons, give or take a fifth",
            story: "The trade winds are not the same all year: they strengthen through the middle \
                    of the year and ease off again, by about a fifth either way. The ocean \
                    breathes with them — the tilt steepens and slackens once a year, twice over \
                    this run.",
            caveat: None,
            basin: BROWSER_BASIN,
            physics: BROWSER_PHYSICS,
            run: TWO_YEARS_EVERY_THREE_DAYS,
            winds: &[WindSection::SeasonalTradeWinds {
                equatorial_zonal_stress_pa: TRADE_WIND_STRESS_PA,
                meridional_decay_scale_m: Some(WAVEGUIDE_SCALE_M),
                relative_amplitude: SEASONAL_RELATIVE_AMPLITUDE,
                peak_time_s: SEASONAL_PEAK_DAY * DAY_S,
            }],
        },
    ];

    /// The preset a panel starts on: the control, because every other story
    /// here is a departure from it.
    #[must_use]
    pub const fn default_preset() -> Self {
        Self::ALL[0]
    }

    /// What the button says.
    #[must_use]
    pub const fn name(&self) -> &'static str {
        self.name
    }

    /// One line on what this run shows.
    #[must_use]
    pub const fn summary(&self) -> &'static str {
        self.summary
    }

    /// What to watch for, in plain words.
    #[must_use]
    pub const fn story(&self) -> &'static str {
        self.story
    }

    /// What this run is not, where a reader could reasonably conclude
    /// something the model does not support.
    #[must_use]
    pub const fn caveat(&self) -> Option<&'static str> {
        self.caveat
    }

    /// This preset as the record a scenario file would have held.
    ///
    /// The same [`ScenarioConfig`] `Scenario::from_toml` would have parsed, so
    /// nothing downstream can tell that no file was involved — and so the
    /// engine's own validation is what accepts or refuses a preset.
    #[must_use]
    pub fn config(&self) -> ScenarioConfig {
        ScenarioConfig {
            basin: self.basin,
            physics: self.physics,
            run: self.run,
            wind: self.winds.to_vec(),
            sst: None,
        }
    }

    /// The runnable scenario, every value through the constructor that checks
    /// it.
    ///
    /// # Errors
    /// Whatever [`ScenarioConfig::build`] objected to. For the shipped presets
    /// that is nothing — `tests/scenario_presets.rs` builds every one of them
    /// — so a failure here means a preset was edited into something the engine
    /// will not run.
    pub fn scenario(&self) -> Result<Scenario, ScenarioError> {
        self.config().build()
    }

    /// What this preset will cost to compute and to hold.
    ///
    /// # Panics
    /// If the preset's basin is not a grid the run format can describe. The
    /// presets are compiled-in constants held to that by
    /// `tests/scenario_presets.rs`, so this is a broken invariant of this
    /// crate rather than anything a reader did (CODING_STANDARDS.md §
    /// *Correctness and failure*), and the panel that asks what a scenario
    /// costs has nowhere to report it to.
    #[must_use]
    pub fn cost(&self) -> PresetCost {
        let bounds = self
            .basin
            .bounds()
            .expect("a preset's basin is a basin the engine accepts");
        let grid = GridSpec::new(bounds.nx(), bounds.ny(), bounds.into())
            .expect("a basin the engine accepts is a grid the format can describe");
        // The engine's own schedule rather than this module's arithmetic on
        // the same three numbers: the frame count is what the budget is spent
        // on, and two places to compute it is one place for them to disagree.
        let schedule = self
            .run
            .build()
            .expect("a preset's run section is a schedule the engine accepts");
        let timing = schedule.timing();
        let steps = self.run.total_steps;
        #[allow(clippy::cast_precision_loss)]
        let cells = (bounds.nx() * bounds.ny()) as f64;
        #[allow(clippy::cast_precision_loss)]
        let native_compute_s = NATIVE_COST_S_PER_CELL_STEP * cells * steps as f64;
        PresetCost {
            nx: bounds.nx(),
            ny: bounds.ny(),
            steps,
            duration_days: schedule.model_time_at_step(steps) / DAY_S,
            frame_count: timing.frame_count,
            frame_interval_s: timing.interval_s,
            frame_bytes: FrameBudget::bytes_of_linear_core(grid, timing.frame_count),
            native_compute_s,
        }
    }
}
