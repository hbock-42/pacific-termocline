//! What the ocean on screen is doing, in sentences — read off the run rather
//! than off a clock.
//!
//! A visitor watching the basin needs to be told what is happening, and the
//! easy way to tell them is a script: at frame 40 say the winds are piling
//! water west, at frame 120 say it is sliding back. That script is a lie
//! waiting to happen. The scenario is a button since T-13.3, the grid, the
//! physics and the length of the run are all a reader's to change, and a
//! sentence keyed to a frame index survives none of it — it would go on
//! narrating a collapse that is not happening, on a run that never had one.
//!
//! So nothing here knows what frame it is on. An [`EquatorialReading`] is what
//! was *measured* along the equator — how deep the thermocline is at each end,
//! how that difference has moved since a measured number of days ago, what the
//! wind stress is and how it compares with the strongest this run reaches,
//! and how far east a recent change has reached — and every caption is a
//! function of that value. `tests/captions.rs` fabricates readings that should
//! and should not produce each one.
//!
//! # The sentence this module may not write
//!
//! It would be very easy to caption the relaxation story as El Niño arriving
//! on schedule, and quite wrong. `docs/enso-oscillation-report.md` records the
//! coupled model's measured period as 1.03 years against an observed 2–7 year
//! band, and records that acceptance criterion as **not met**; T-12.3 refused
//! to tune the model until it agreed. A caption that spoke of a cycle, of an
//! ENSO period, or of what the ocean is about to do next would quietly
//! contradict the project's own validation ledger — which is the most valuable
//! thing in this repository — so the captions describe the present state and
//! nothing else, and `tests/captions.rs` sweeps every state this module can be
//! in for that vocabulary.
//!
//! Every number a caption states is the run's: the metres come off the two
//! ends of this run's equator, the pascals off its wind stress, the days off
//! its own cadence ([`crate::RunClock`]) and the longitude off the basin its
//! header declares. There are no numbers about the ocean in this file.
//!
//! Device-free like every view beside it ([ADR-0006]): a caption is a string
//! and a topic, so what the panel says is asserted rather than looked at.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use std::error::Error;
use std::fmt;

use termocline_format::FormatError;

use crate::cross_section::CrossSection;
use crate::geography::longitude_text;
use crate::wind::CellCentreStress;
use crate::wording::{PlainTerm, THERMOCLINE, TRADE_WINDS};
use crate::{LoadedRun, RunClock};

/// How far back a change is measured, in days.
///
/// Thirty days is short against both timescales this ocean moves on — an
/// equatorial Kelvin wave crosses the 17 800 km basin in 17 800 km / 3.0 m/s ≈
/// 69 days, and the Rayleigh damping adjusts the tilt on `r⁻¹` ≈ 116 days
/// (`CONTEXT.md`) — so a difference taken over it is a change happening *now*
/// rather than the whole history of the run. It is also long against the one-
/// to three-day cadence the presets write at, so it spans enough frames that a
/// single frame's value cannot set the trend.
///
/// It is a *requested* window. How far back the run can actually be read is
/// [`EquatorialReading::lookback_days`], which is what the captions state:
/// near the start of a run there is less history than this, and a caption that
/// said "30 days" anyway would be quoting a constant instead of the run.
const LOOKBACK_DAYS: f64 = 30.0;

/// A tilt smaller than this share of the run's full anomaly range counts as no
/// tilt at all.
///
/// The scale every view of the run is drawn on reaches the same distance
/// either side of zero ([`crate::DivergingScale`]), so an equilibrium profile
/// — deepest in the west, shallowest in the east — spans very nearly the whole
/// of it, and the full range is the tilt this run is *about*. A tenth of that
/// is a slope a reader cannot pick out of the picture, and calling it one
/// would be describing something the screen does not show.
///
/// A share of the run's own range rather than a count of metres, because the
/// metres are a scenario's: the same threshold has to work for a run whose
/// tilt reaches 66 m and one whose wind is a tenth as strong.
const LEVEL_FRACTION: f64 = 0.1;

/// Below this share of the strongest easterly its equator carries, the run's
/// wind counts as slack.
///
/// What the threshold has to do is separate a wind merely breathing with the
/// year from one that has actually let go, and the two presets that do those
/// things say where it goes. The seasonal preset swings the alizés by ±20 %
/// about their mean (`crate::presets`), so its equatorial stress runs between
/// 0.04 and 0.06 Pa and never falls below 0.04/0.06 = 0.67 of the strongest
/// the run sees. The relaxation preset lays a westerly Gaussian of +0.04 Pa
/// and zonal scale `Lx` = 6 000 km over −0.05 Pa of alizés across a basin
/// 17 810 km wide, and the mean of that Gaussian along the equator is
/// `Lx·√π·erf(L/2Lx)/L` = 0.575 of its peak — so at the height of the
/// relaxation the equatorial wind is 0.027 Pa, which is 0.54 of the strongest.
///
/// Six tenths sits between the two, near enough the middle of them: above
/// every point of a seasonal cycle and below the peak of a relaxation. Both
/// figures are ratios of the *same* quantity the comparison is made on — the
/// mean zonal stress along the equator, against the strongest that mean
/// reaches in the run — and it is a share rather than a count of pascals, so a
/// scenario forced ten times harder is judged by its own scale.
const SLACK_FRACTION: f64 = 0.6;

/// A change smaller than this share of the run's full anomaly range, taken
/// over the window, is not a change worth a sentence.
///
/// Two per cent. The tilt relaxes towards whatever the wind implies on the
/// Rayleigh timescale `r⁻¹` ≈ 116 days (`CONTEXT.md`), so a wind that has
/// genuinely changed drags the tilt a good fraction of the way to its new
/// equilibrium over the thirty days read back over — tens of per cent of the
/// range — while a wind that has not changed leaves it where it is. Two per
/// cent is an order of magnitude below the first and well above the drift of
/// the second, so what it separates is motion from stillness rather than one
/// speed from another.
///
/// It is deliberately *not* [`LEVEL_FRACTION`]. Whether there is a tilt at all
/// and whether the tilt is moving are different questions, and a difference
/// far too small to see can still be closing fast enough to be the thing on
/// screen.
const MOVING_FRACTION: f64 = 0.02;

/// The share of its own peak at which a change counts as having reached a
/// place.
///
/// Half, the convention a pulse's width is quoted at — full width at half
/// maximum. The leading edge of a departure is then the easternmost point
/// still carrying half of the largest change anywhere, which is a place a
/// reader can be pointed at rather than the tail of a signal that decays
/// smoothly to nothing.
const FRONT_FRACTION: f64 = 0.5;

/// The most captions one reading can produce.
///
/// One per topic of [`CaptionTopic`]: the panel reserves room for this many
/// rows under the ocean so that turning captions on shrinks the picture rather
/// than pushing anything off the window, and `tests/captions.rs` holds the
/// module to it.
pub const MAX_CAPTIONS: usize = 5;

/// What one caption is about.
///
/// Carried so that a test can ask whether *this* caption fired without
/// matching on its prose, and so that the panel could one day order or filter
/// them. The wording is the module's; the topic is the contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CaptionTopic {
    /// The wind measured along the equator.
    Wind,
    /// How much deeper the thermocline is at one end of the basin than at the
    /// other.
    Tilt,
    /// Whether that difference is growing, shrinking or holding.
    TiltTrend,
    /// The warm water piled up in the west running back east under a slack
    /// wind — the one compound state in the set.
    WarmWaterEast,
    /// How far east a recent change has reached — which is a statement about
    /// where the ocean has and has not moved, and deliberately not one about
    /// something travelling: what is measured is a profile against an earlier
    /// profile, and the model's own wave speeds are Epic 04's to assert.
    Front,
}

/// One sentence about the run on screen.
///
/// Text and a topic, and no colour: a caption is read, and a reader who cannot
/// see the picture's colours — or has the picture switched off — loses nothing
/// of it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Caption {
    /// What it is about.
    topic: CaptionTopic,
    /// The sentence itself.
    text: String,
}

impl Caption {
    /// What this caption is about.
    #[must_use]
    pub const fn topic(&self) -> CaptionTopic {
        self.topic
    }

    /// The sentence.
    #[must_use]
    pub fn text(&self) -> &str {
        &self.text
    }
}

/// How far east a recent change has reached, as measured against the run's
/// own earlier state.
///
/// The difference between the equatorial profile now and the profile
/// [`EquatorialReading::lookback_days`] ago, reduced to the three facts a
/// caption needs: how big the largest change is, where its leading edge is,
/// and whether that edge is still inside the basin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Departure {
    /// The largest change in thermocline depth anywhere along the equator over
    /// the window, in metres, with its sign: positive where the thermocline
    /// deepened, negative where it shoaled.
    pub peak_change_m: f64,
    /// The longitude of the easternmost point that changed by at least half
    /// that peak *in the same direction*, in degrees east of the prime
    /// meridian.
    pub half_peak_longitude_deg_east: f64,
    /// Whether that point is the easternmost column of the basin — in which
    /// case the change has run out of ocean and there is no leading edge left
    /// to point at.
    pub reaches_eastern_wall: bool,
}

/// What the run is doing along the equator at one frame, measured.
///
/// Public fields rather than accessors, like [`crate::PresetCost`]: this is a
/// record of measurements, and the point of it is that a test can fabricate
/// one directly and ask what the captions say about it. Nothing is derived
/// here that a caption then re-derives — the thresholds are
/// [`EquatorialReading::captions`]'s, because they are judgements about what is
/// worth saying rather than facts about the ocean.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquatorialReading {
    /// How far the thermocline sits below the sea surface at the western end
    /// of the equator, in metres — the total depth `H + h`, never the anomaly
    /// (`CONTEXT.md`). `None` where the model has broken down there and there
    /// is no depth to report ([`crate::SideViewColumn::interface_fraction`]).
    pub west_depth_m: Option<f64>,
    /// The same at the eastern end.
    pub east_depth_m: Option<f64>,
    /// The difference between those two depths as it stood
    /// [`EquatorialReading::lookback_days`] ago, in metres, west minus east.
    /// `None` where the run does not reach back that far, or where the model
    /// had broken down at one end of it.
    pub earlier_tilt_m: Option<f64>,
    /// How far back that earlier tilt was read, in days of model time. The gap
    /// the run could actually supply, which near the start of a run is shorter
    /// than the window that was asked for.
    pub lookback_days: f64,
    /// The mean zonal wind stress along the equator, in pascals. Negative is
    /// easterly — the trade winds, which blow east to west (`CONTEXT.md`).
    /// `None` where the run gave no wind to read, which is a different thing
    /// from a wind of zero and is captioned as nothing at all.
    pub equatorial_stress_pa: Option<f64>,
    /// The strongest easterly wind this run's equator carries, in pascals, as
    /// a magnitude ([`crate::LoadedRun::strongest_equatorial_easterly_pa`]).
    ///
    /// What "slack" is measured against, so that a run driven by a weak wind
    /// is not permanently described as slack. It is the same quantity
    /// [`EquatorialReading::equatorial_stress_pa`] is — the mean zonal stress
    /// along the equator — taken in the frame where it is strongest, so the
    /// two can honestly be quoted against each other in one sentence.
    pub strongest_easterly_pa: f64,
    /// Half the run's anomaly range, in metres: the same
    /// [`crate::DivergingScale`] every view of the run is drawn on. What a
    /// difference is judged large or small against.
    pub anomaly_half_range_m: f64,
    /// How far east a recent change has travelled, where one was measured.
    pub departure: Option<Departure>,
}

impl EquatorialReading {
    /// Read the run at frame `index`: the equator now, the equator as it stood
    /// a window ago, and the wind between them.
    ///
    /// The window is [`LOOKBACK_DAYS`] rounded to whole frames of this run's
    /// own cadence, and at least one frame, so a run written every three days
    /// and one written every hour both look back over about a month rather
    /// than over a fixed number of frames.
    ///
    /// # Errors
    /// [`ReadingError::NoSuchFrame`] where the run holds no frame `index`, and
    /// [`ReadingError::Frame`] where a frame does not fit the grid its own
    /// header declares.
    pub fn of_run(run: &LoadedRun, index: u64) -> Result<Self, ReadingError> {
        let header = run.header();
        let grid = header.grid;
        let scale = run.anomaly_scale();
        let clock = RunClock::of_run(header.output);
        let mean_depth_m = header.physical_params.mean_depth_m;

        let frame = run.frame(index).ok_or(ReadingError::NoSuchFrame {
            index,
            frame_count: run.frame_count(),
        })?;
        let now = CrossSection::of_frame(grid, &frame, scale)?;
        let equatorial_stress_pa =
            CellCentreStress::of_frame(grid, &frame)?.equatorial_zonal_mean_pa(grid.extent());

        let earlier_index = index.saturating_sub(frames_in_lookback(header.output.interval_s));
        let earlier = if earlier_index == index {
            None
        } else {
            run.frame(earlier_index)
                .map(|frame| CrossSection::of_frame(grid, &frame, scale))
                .transpose()?
        };

        Ok(Self {
            west_depth_m: end_depth(&now, End::West, mean_depth_m),
            east_depth_m: end_depth(&now, End::East, mean_depth_m),
            earlier_tilt_m: earlier.as_ref().and_then(|earlier| {
                Some(
                    end_depth(earlier, End::West, mean_depth_m)?
                        - end_depth(earlier, End::East, mean_depth_m)?,
                )
            }),
            lookback_days: clock.day_of_frame(index) - clock.day_of_frame(earlier_index),
            equatorial_stress_pa,
            strongest_easterly_pa: run.strongest_equatorial_easterly_pa(),
            anomaly_half_range_m: scale.half_range_m(),
            departure: earlier
                .as_ref()
                .and_then(|earlier| Departure::between(earlier, &now)),
        })
    }

    /// How much deeper the thermocline is in the west than in the east, in
    /// metres, or `None` where either end is not a depth the run produced.
    #[must_use]
    pub fn tilt_m(&self) -> Option<f64> {
        Some(self.west_depth_m? - self.east_depth_m?)
    }

    /// What this reading says, in sentences, in the order a panel shows them.
    ///
    /// Every one of them is a branch on the measurements above and on nothing
    /// else. Where a measurement is missing the sentence that would have used
    /// it is missing too: silence is the only honest thing to say about a
    /// quantity the run did not produce.
    ///
    /// A term of the project's vocabulary is glossed the first time the set
    /// uses it and named plainly after that (`crate::wording`), so a reader
    /// meets the explanation once instead of in every line.
    #[must_use]
    pub fn captions(&self) -> Vec<Caption> {
        let mut words = Glossary::new();
        let mut captions = Vec::with_capacity(MAX_CAPTIONS);
        let mut push = |topic, text| captions.push(Caption { topic, text });

        if let Some(text) = self.wind_sentence(&mut words) {
            push(CaptionTopic::Wind, text);
        }
        if let Some(tilt_m) = self.tilt_m() {
            push(CaptionTopic::Tilt, self.tilt_sentence(tilt_m, &mut words));
            if let Some(spread_change_m) = self.spread_change_m() {
                push(
                    CaptionTopic::TiltTrend,
                    self.trend_sentence(spread_change_m),
                );
                if let Some(text) = self.running_back_east_sentence(tilt_m, spread_change_m) {
                    push(CaptionTopic::WarmWaterEast, text);
                }
            }
        }
        if let Some(text) = self.front_sentence(&mut words) {
            push(CaptionTopic::Front, text);
        }
        captions
    }

    /// `fraction` of the run's full anomaly range, in metres — twice the
    /// half-range the run's scale carries.
    ///
    /// The two thresholds that judge a number of metres large or small are
    /// both shares of this, so they move with the run rather than with the
    /// scenario that happened to be written first.
    fn fraction_of_range_m(&self, fraction: f64) -> f64 {
        fraction * 2.0 * self.anomaly_half_range_m
    }

    /// The metres below which a difference between the two ends of the basin
    /// is no difference at all.
    fn level_tilt_m(&self) -> f64 {
        self.fraction_of_range_m(LEVEL_FRACTION)
    }

    /// The metres a difference has to move over the window before anything is
    /// said to have changed.
    fn moving_m(&self) -> f64 {
        self.fraction_of_range_m(MOVING_FRACTION)
    }

    /// What the wind along the equator is doing, or `None` where the run gave
    /// no wind to read.
    ///
    /// Classified once. Two sentences below turn on which of the four states
    /// this is, and deciding it twice is two places for them to come to
    /// disagree about the same run.
    fn wind_state(&self) -> Option<WindState> {
        let stress_pa = self.equatorial_stress_pa?;
        if !stress_pa.is_finite() {
            return None;
        }
        if stress_pa > 0.0 {
            return Some(WindState::Reversed);
        }
        if stress_pa == 0.0 {
            return Some(WindState::Calm);
        }
        if stress_pa.abs() < SLACK_FRACTION * self.strongest_easterly_pa {
            return Some(WindState::Slack);
        }
        Some(WindState::Blowing)
    }

    /// How much the *size* of the difference between the two ends has changed
    /// over the window, in metres: positive where the two ends have moved
    /// apart, negative where they have come together.
    ///
    /// The size rather than the signed tilt, so a tilt that collapses through
    /// zero and builds the other way reads as shrinking and then growing,
    /// which is what the picture does.
    fn spread_change_m(&self) -> Option<f64> {
        Some(self.tilt_m()?.abs() - self.earlier_tilt_m?.abs())
    }

    /// The wind, in one sentence, or `None` where no stress was measured.
    ///
    /// Each arm glosses the trade winds for itself rather than ahead of the
    /// match, because the one that names no wind names no trade winds either,
    /// and a term glossed into a sentence it does not appear in is a gloss
    /// spent.
    fn wind_sentence(&self, words: &mut Glossary) -> Option<String> {
        let stress_pa = self.equatorial_stress_pa?;
        Some(match self.wind_state()? {
            WindState::Calm => "No wind is blowing on the equator in this frame.".to_owned(),
            WindState::Reversed => format!(
                "The wind over the equator has turned around: {stress_pa:.2} Pa blowing west to \
                 east, the opposite of the {}.",
                words.say(TRADE_WINDS),
            ),
            WindState::Slack => format!(
                "The {} have gone slack along the equator: {:.2} Pa, against {:.2} Pa at their \
                 strongest along this run's equator.",
                words.say(TRADE_WINDS),
                stress_pa.abs(),
                self.strongest_easterly_pa,
            ),
            WindState::Blowing => format!(
                "The {} are blowing east to west along the equator at {:.2} Pa, piling warm \
                 water into the west.",
                words.say(TRADE_WINDS),
                stress_pa.abs(),
            ),
        })
    }

    /// The difference between the two ends, in one sentence.
    fn tilt_sentence(&self, tilt_m: f64, words: &mut Glossary) -> String {
        let thermocline = words.say(THERMOCLINE);
        let level_m = self.level_tilt_m();
        if tilt_m.abs() <= level_m {
            return format!(
                "The {thermocline} sits at nearly the same depth right across the basin: {:.0} m \
                 between the two ends.",
                tilt_m.abs(),
            );
        }
        let (deep, shallow) = if tilt_m > 0.0 {
            ("west", "east")
        } else {
            ("east", "west")
        };
        format!(
            "The {thermocline} is now {:.0} m deeper in the {deep} than in the {shallow}.",
            tilt_m.abs(),
        )
    }

    /// Whether that difference is growing, shrinking or holding, in one
    /// sentence.
    fn trend_sentence(&self, spread_change_m: f64) -> String {
        let days = self.lookback_days;
        let moving_m = self.moving_m();
        if spread_change_m > moving_m {
            return format!(
                "That difference is still growing: {spread_change_m:.0} m more over the last \
                 {days:.0} days."
            );
        }
        if spread_change_m < -moving_m {
            return format!(
                "That difference is shrinking: {:.0} m less over the last {days:.0} days.",
                spread_change_m.abs(),
            );
        }
        format!(
            "That difference is holding: it has moved {:.0} m in the last {days:.0} days.",
            spread_change_m.abs(),
        )
    }

    /// The one compound state in the set, in a sentence, or `None` where it
    /// does not hold.
    ///
    /// Three measurements have to agree before this can be said: the wind has
    /// let go, the difference between the two ends is closing, and the west is
    /// still the deep end — so there is a pile of warm water there for the
    /// ocean to run back *from*. Any one of them missing and the sentence
    /// would be describing something else's run.
    ///
    /// How the wind let go is part of what is said, because "slack" and
    /// "blowing the other way" are different states and the run measured which
    /// one this is.
    fn running_back_east_sentence(&self, tilt_m: f64, spread_change_m: f64) -> Option<String> {
        let wind = self.wind_state()?.no_longer_holding_the_west()?;
        if spread_change_m >= -self.moving_m() || tilt_m <= 0.0 {
            return None;
        }
        Some(format!(
            "{wind}, the warm water piled up in the west is sliding back east."
        ))
    }

    /// How far east the recent change has reached, in one sentence, or `None`
    /// where nothing with a leading edge was measured.
    ///
    /// It states where the ocean has moved and where it has not, which is what
    /// was measured. It does not say that anything is travelling: that would
    /// be a claim about a wave speed, and the reading holds two profiles
    /// rather than a trajectory.
    fn front_sentence(&self, words: &mut Glossary) -> Option<String> {
        let departure = self.departure?;
        if departure.reaches_eastern_wall || departure.peak_change_m <= self.moving_m() {
            return None;
        }
        let thermocline = words.say(THERMOCLINE);
        Some(format!(
            "The change of the last {:.0} days reaches as far east as {}; east of that the \
             {thermocline} has barely moved.",
            self.lookback_days,
            longitude_text(departure.half_peak_longitude_deg_east),
        ))
    }
}

/// What the wind along the equator is doing, as the measurements classify it.
///
/// Four states rather than a number, because the sentences turn on which of
/// them it is: whether the trades are described as piling water west, as
/// having gone slack against the strongest this run reaches, as having turned
/// around, or as absent. Two captions read it, and this is where the reading
/// happens.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WindState {
    /// Easterly, and at a good share of the strongest stress the run reaches.
    Blowing,
    /// Easterly, but well below that — [`SLACK_FRACTION`] of it.
    Slack,
    /// No stress on the equator at all.
    Calm,
    /// Westerly: blowing the opposite way to the trades.
    Reversed,
}

impl WindState {
    /// How to open a sentence about water running back east, or `None` where
    /// the wind is still holding it in the west.
    ///
    /// A wind that has let go has let go in one of three ways, and which one
    /// was measured is part of what the caption says: "slack" and "blowing the
    /// other way" are different states of the ocean's forcing.
    const fn no_longer_holding_the_west(self) -> Option<&'static str> {
        match self {
            Self::Blowing => None,
            Self::Slack => Some("With the winds slack"),
            Self::Calm => Some("With no wind holding it there"),
            Self::Reversed => Some("With the wind pushing the other way"),
        }
    }
}

impl Departure {
    /// The change between two equatorial sections of the same run, reduced to
    /// what a caption can say about it, or `None` where nothing finite was
    /// measured.
    ///
    /// The two sections come from the same run and so have the same columns;
    /// a caller that mixed two runs gets `None` rather than a change read
    /// across mismatched basins.
    fn between(earlier: &CrossSection, now: &CrossSection) -> Option<Self> {
        let (before, after) = (earlier.points(), now.points());
        if before.len() != after.len() || after.is_empty() {
            return None;
        }
        let change_m: Vec<f64> = before
            .iter()
            .zip(after)
            .map(|(before, after)| after.h_m() - before.h_m())
            .collect();
        // The peak keeps its sign, and the leading edge is looked for in that
        // same sign: a relaxation deepens the east while it shoals the west,
        // and an edge picked on magnitude alone would point at whichever of
        // the two happened to reach furthest — a longitude belonging to the
        // other half of the signal.
        let peak_change_m =
            change_m
                .iter()
                .filter(|metres| metres.is_finite())
                .fold(0.0_f64, |peak, metres| {
                    if metres.abs() > peak.abs() {
                        *metres
                    } else {
                        peak
                    }
                });
        if peak_change_m == 0.0 {
            return None;
        }
        let edge = FRONT_FRACTION * peak_change_m;
        let reached = |metres: f64| {
            metres.is_finite()
                && if peak_change_m > 0.0 {
                    metres >= edge
                } else {
                    metres <= edge
                }
        };
        let index = change_m.iter().rposition(|&metres| reached(metres))?;
        Some(Self {
            peak_change_m,
            half_peak_longitude_deg_east: after[index].longitude_deg_east(),
            reaches_eastern_wall: index + 1 == after.len(),
        })
    }
}

/// Which end of the basin a depth was read at.
#[derive(Debug, Clone, Copy)]
enum End {
    West,
    East,
}

/// The thermocline depth at one end of `section`, in metres below the sea
/// surface, or `None` where the model has broken down there.
///
/// The same three refusals the side view makes ([`crate::SideViewColumn`]): a
/// non-finite anomaly means the integration diverged, and an interface at or
/// above the sea surface means the upper layer has no thickness left and the
/// linear model has stopped describing anything.
fn end_depth(section: &CrossSection, end: End, mean_depth_m: f64) -> Option<f64> {
    let point = match end {
        End::West => section.points().first(),
        End::East => section.points().last(),
    }?;
    let depth_m = mean_depth_m + point.h_m();
    (depth_m.is_finite() && depth_m > 0.0).then_some(depth_m)
}

/// How many frames of a run written every `interval_s` seconds span
/// [`LOOKBACK_DAYS`], and never fewer than one.
fn frames_in_lookback(interval_s: f64) -> u64 {
    let window_s = LOOKBACK_DAYS * crate::run::SECONDS_PER_DAY;
    if !interval_s.is_finite() || interval_s <= 0.0 {
        return 1;
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let frames = (window_s / interval_s).round() as u64;
    frames.max(1)
}

/// The terms a set of captions has already explained.
///
/// A gloss is worth reading once. The first caption to use a term carries its
/// plain-language meaning with it, and every caption after that in the same set
/// uses the bare term — so a reader meets "thermocline (the boundary between
/// …)" once and then reads sentences rather than parentheses.
struct Glossary {
    /// The terms already glossed, in the order they were used.
    said: Vec<&'static str>,
}

impl Glossary {
    /// A set of captions that has explained nothing yet.
    const fn new() -> Self {
        Self { said: Vec::new() }
    }

    /// `term` as this set should write it here: glossed on its first use,
    /// bare after that.
    fn say(&mut self, term: PlainTerm) -> String {
        if self.said.contains(&term.term()) {
            return term.term().to_owned();
        }
        self.said.push(term.term());
        term.glossed()
    }
}

/// What stopped a reading being taken.
#[derive(Debug, Clone, PartialEq)]
pub enum ReadingError {
    /// The run holds no frame at that index.
    NoSuchFrame {
        /// The index that was asked for.
        index: u64,
        /// How many frames the run holds.
        frame_count: u64,
    },
    /// A frame of the run did not fit the grid its own header declares.
    Frame(FormatError),
}

impl From<FormatError> for ReadingError {
    fn from(error: FormatError) -> Self {
        Self::Frame(error)
    }
}

impl fmt::Display for ReadingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoSuchFrame { index, frame_count } => write!(
                f,
                "this run holds {frame_count} frames and was asked for frame {index}"
            ),
            Self::Frame(error) => error.fmt(f),
        }
    }
}

impl Error for ReadingError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::NoSuchFrame { .. } => None,
            Self::Frame(error) => Some(error),
        }
    }
}
