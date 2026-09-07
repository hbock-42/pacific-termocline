//! Which of the two screens the shell is on, and which screen it opens on.
//!
//! Epics 08–09 built an instrument: a basin map of `h` on a diverging scale,
//! the wind stress over it, a scrubber and a clock, a cross-section, a point
//! time series, and two runs held against each other. Epic 13 built a
//! *reading* of the same run: the ocean in cross-section at real depths
//! ([`crate::side_view`]), longitudes and coasts and months
//! ([`crate::geography`], [`crate::clock`]), scenarios as stories
//! ([`crate::presets`]), and sentences measured off the equator
//! ([`crate::captions`]).
//!
//! Both are drawn from the same [`crate::LoadedRun`] and neither is removed.
//! [`Mode`] is only which one is on screen.
//!
//! # Why the two targets open differently
//!
//! [`Platform`] is what decides that, and it decides it for a reason that is
//! about *who arrived*, not about what the machine can do.
//!
//! A visitor to the web build asked for a URL. Per [ADR-0012] there is nothing
//! for them to open — the file format is not served to the browser at all — so
//! the shell computes a scenario and the only question left is what it draws
//! of it. Someone who has not met the word *thermocline* is the expected
//! visitor there, and the instrument's first row is a diverging colour scale
//! in metres of anomaly. So the browser opens on [`Mode::Teaching`].
//!
//! Someone starting the native binary asked for something else. The desktop
//! build takes run directories on the command line — `termocline-viz
//! /tmp/run-control /tmp/run-burst` opens a comparison — reads runs the engine
//! wrote, and is the only build that can. Nobody reaches it by accident, and a
//! reader who has just named two runs on a command line has said what they
//! want: the instrument. Opening them on a captioned ocean would be a screen
//! they have to dismiss before doing the thing they launched the program for.
//! So the desktop opens on [`Mode::Scientific`], and the story is one click
//! away for the times it is what is wanted — a lecture given off a laptop,
//! say.
//!
//! Neither default is a claim that one screen is for experts. The toggle is
//! the same toggle on both, and it starts where the arriving reader most
//! likely wanted it.
//!
//! [ADR-0012]: ../../docs/planning/adr/0012-the-browser-runs-the-engine.md

/// Which build of the shell this is.
///
/// A value rather than a `cfg!` read at the point of use, so that "the browser
/// opens on the story and the desktop opens on the instrument" is a total
/// function this crate's tests can evaluate both branches of. A test compiled
/// for the host can otherwise only ever assert the host's own answer, which is
/// half a rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Platform {
    /// The `wasm32` build, served as a page (ADR-0006, ADR-0012).
    Browser,
    /// The native build: a window, a command line, and a filesystem.
    Desktop,
}

impl Platform {
    /// The platform this binary was compiled for.
    #[must_use]
    pub const fn of_target() -> Self {
        if cfg!(target_arch = "wasm32") {
            Self::Browser
        } else {
            Self::Desktop
        }
    }
}

/// One of the views Epics 08–09 built, by name.
///
/// The list the ticket names, written down so that "the toggle reaches every
/// scientific view" is a statement about an enumeration rather than about
/// whatever happens to be drawn: [`Mode::draws`] is a total function of it, so
/// a view added here has to be answered for in both modes before this crate
/// compiles, and `tests/teaching_mode.rs` reads that answer view by view.
///
/// Two of the seven are also what [`crate::app`] asks about at runtime, and
/// they are the two whose answer changes the *shape* of the shell rather than
/// what one panel draws: [`ScientificView::Heatmap`] chooses which of the two
/// readings of a run is drawn at all, and [`ScientificView::RunComparison`]
/// chooses whether the second panel is on screen. The other five ride on the
/// first — they are drawn under, over and beside the basin map, so a shell
/// that is not drawing the map is not drawing them either.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScientificView {
    /// The basin map of `h`, on the run's diverging colour scale and under its
    /// colour bar (T-08.1, T-08.2).
    Heatmap,
    /// The wind-stress arrows over the map (T-09.1).
    WindOverlay,
    /// The frame chooser: the slider and its keys (T-08.3).
    Scrubber,
    /// The clock that chooses frames on the reader's behalf (T-09.2).
    Playback,
    /// The equatorial cross-section of `h` against longitude (T-09.3).
    CrossSection,
    /// The time series of one picked cell over the whole run (T-09.4).
    PointTimeSeries,
    /// Two runs side by side on one index and one scale (T-09.5).
    RunComparison,
}

impl ScientificView {
    /// Every view the toggle has to reach, in the order the shell draws them.
    pub const ALL: [Self; 7] = [
        Self::Heatmap,
        Self::WindOverlay,
        Self::Scrubber,
        Self::Playback,
        Self::CrossSection,
        Self::PointTimeSeries,
        Self::RunComparison,
    ];

    /// What this view is called on screen.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Heatmap => "Basin map",
            Self::WindOverlay => "Wind stress τ",
            Self::Scrubber => "Frame chooser",
            Self::Playback => "Playback",
            Self::CrossSection => "Equatorial cross-section",
            Self::PointTimeSeries => "Point time series",
            Self::RunComparison => "Compare two runs",
        }
    }
}

/// Which screen the shell is on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// The story: the ocean in cross-section, in metres and months, with the
    /// scenario's own account of itself and captions read off the run.
    Teaching,
    /// The instrument: every view of Epics 08–09.
    Scientific,
}

impl Mode {
    /// The mode a shell on `platform` opens in — see the module's own note on
    /// why the two differ.
    #[must_use]
    pub const fn opening_on(platform: Platform) -> Self {
        match platform {
            Platform::Browser => Self::Teaching,
            Platform::Desktop => Self::Scientific,
        }
    }

    /// The other mode: what the one control in the shell switches to.
    #[must_use]
    pub const fn other(self) -> Self {
        match self {
            Self::Teaching => Self::Scientific,
            Self::Scientific => Self::Teaching,
        }
    }

    /// Whether this mode draws `view`.
    ///
    /// [`Mode::Scientific`] draws all seven, which is the whole of what the
    /// toggle has to reach. [`Mode::Teaching`] keeps two of them, and the two
    /// are the ones that are not readings of the ocean but ways of moving
    /// through time: a story that cannot be replayed is one a reader sees once
    /// and only if they were already looking. The other five are readings —
    /// an anomaly in metres on a diverging scale, a stress in pascals, a
    /// chart, a series, a second run — and each is a thing the story says in
    /// words instead.
    #[must_use]
    pub const fn draws(self, view: ScientificView) -> bool {
        match self {
            Self::Scientific => true,
            Self::Teaching => matches!(view, ScientificView::Scrubber | ScientificView::Playback),
        }
    }

    /// What the control that leaves this mode says.
    ///
    /// One button rather than a pair of tabs or a labelled switch, and it
    /// names the *destination*: a reader who wants the instrument is looking
    /// for the word, and a reader who wants the story never has to work out
    /// which half of a switch they are on.
    #[must_use]
    pub const fn switch_label(self) -> &'static str {
        match self {
            Self::Teaching => "Scientific views",
            Self::Scientific => "Teaching view",
        }
    }

    /// What that control says on hover: what is on the other side of it, named
    /// view by view, so the button is an offer rather than a mystery.
    #[must_use]
    pub const fn switch_hover(self) -> &'static str {
        match self {
            Self::Teaching => {
                "The instrument underneath: the basin map of the depth anomaly and its colour \
                 scale, the wind stress over it, the equatorial cross-section, the time series of \
                 a cell you pick, and two runs side by side."
            }
            Self::Scientific => {
                "The same run as a story: the ocean seen from the side at real depths, in \
                 longitudes, metres and months, with captions read off this run."
            }
        }
    }
}
