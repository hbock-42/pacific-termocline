//! Two runs side by side: when that says something true, and when it does not.
//!
//! The view this backs is the one a reader reaches for to answer "what did the
//! wind burst do?" — the control run in one panel, the perturbed run in the
//! other, both on the same frame. Everything with a value in it is here, and
//! none of it knows what a GPU is: two loaded runs in, a scale and a verdict
//! out, so what the ticket asks for can be asserted in `tests/comparison.rs`
//! without a device (ADR-0006).
//!
//! # One scale, not two
//!
//! Every view in this crate draws a run on a scale that covers the *whole
//! run* rather than the frame in front of it, so that a tilt which collapses
//! is seen to collapse instead of being renormalized back to full saturation
//! ([`crate::heatmap`]). A comparison is that argument again, one level up: two
//! runs each drawn on its own scale would look identical however far apart
//! they are, because each would reach the ends of the same ramp. So a
//! comparison has exactly one scale — [`DivergingScale::widened`] over both
//! runs — and the quieter run occupies the middle of it, which is the true
//! reading. The wind overlay's scale is shared for the same reason.
//!
//! When the two runs have genuinely different ranges, the louder one sets the
//! scale and the difference is stated in words as well as drawn
//! ([`Difference::AnomalyRange`]), so a reader is told why one panel is pale
//! rather than left to infer that its run did nothing.
//!
//! # What comparison requires, and what it merely notes
//!
//! Index-synced panels claim two things: that a cell in one panel is the same
//! patch of ocean as the cell beside it, and that a frame in one panel is the
//! same model time as the frame beside it. Two runs that cannot support those
//! claims are **refused** ([`Mismatch`]) rather than drawn misleadingly.
//! Everything else two runs may differ in — length, forcing, parameters,
//! whether they couple SST — is what a comparison is *for*, so it is carried
//! and stated ([`Difference`]).

use std::fmt;

use termocline_format::{GridSpec, PhysicalParams, Variable};

use crate::run::{grid_description, SECONDS_PER_DAY};
use crate::{DivergingScale, LoadedRun, StressScale};

/// Which of the two panels of a comparison something belongs to.
///
/// Two rather than a list of runs: a side-by-side is a pair, and it is the
/// pair that lets one frame index mean one moment in both panels. A third run
/// would need a third meaning for the index, and there is none.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    /// The panel a single run is shown in, and the left of a comparison.
    Left,
    /// The panel the run being compared against is shown in.
    Right,
}

impl Side {
    /// Both sides, in the order they are drawn.
    pub const BOTH: [Self; 2] = [Self::Left, Self::Right];

    /// Where this side's state sits in the shell's pairs.
    #[must_use]
    pub const fn index(self) -> usize {
        match self {
            Self::Left => 0,
            Self::Right => 1,
        }
    }

    /// What the reader calls this panel.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Left => "A",
            Self::Right => "B",
        }
    }
}

/// Why two runs cannot be shown side by side.
///
/// Both variants are claims the side-by-side would make silently and could not
/// support. Neither is bad input — each run loaded fine on its own — so this
/// is a refusal to draw, reported to the reader in the words below.
#[derive(Debug, Clone, PartialEq)]
pub enum Mismatch {
    /// The runs are not over the same cells: a different resolution, a
    /// different basin, or both.
    ///
    /// The panels are drawn side by side at the same size, so the same place
    /// on screen is the same cell of the grid. Two grids that are not the same
    /// grid make that false: at two resolutions a pixel covers two different
    /// patches of ocean, and over two basins it covers two different oceans.
    Grid {
        /// The left panel's grid.
        left: GridSpec,
        /// The right panel's grid.
        right: GridSpec,
    },
    /// The runs wrote frames at different cadences, so one index is not one
    /// model time.
    ///
    /// The frame index is what the two panels share (`crate::app`). Where the
    /// output intervals differ, sharing it syncs the panels to the same
    /// *number* while showing two different moments, which is precisely the
    /// comparison a reader would believe and should not.
    FrameInterval {
        /// The left run's output interval, in seconds.
        left_s: f64,
        /// The right run's output interval, in seconds.
        right_s: f64,
    },
}

impl fmt::Display for Mismatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Grid { left, right } => write!(
                f,
                "these runs are not over the same cells — {} against {} — so the same place \
                 in the two panels would not be the same place in the ocean",
                grid_description(*left),
                grid_description(*right)
            ),
            Self::FrameInterval { left_s, right_s } => write!(
                f,
                "these runs wrote frames at different cadences — every {left_s} s against every \
                 {right_s} s — so one frame index is not one model time",
            ),
        }
    }
}

/// Something the two runs differ in that a comparison carries rather than
/// refuses — and says out loud, because it changes how the panels read.
#[derive(Debug, Clone, PartialEq)]
pub enum Difference {
    /// The runs hold different numbers of frames.
    Length {
        /// Frames in the left run.
        left_frames: u64,
        /// Frames in the right run.
        right_frames: u64,
        /// Frames the comparison covers: the shorter of the two.
        shared_frames: u64,
    },
    /// One run couples SST and the other does not (T-05.4).
    SstAnomaly {
        /// Whether the left run's frames carry `T'`.
        left: bool,
        /// Whether the right run's frames carry `T'`.
        right: bool,
    },
    /// The runs were integrated with a different value of one physical
    /// parameter — often the reason they are being compared at all.
    PhysicalParam {
        /// The parameter, labelled as the metadata panel labels it.
        name: &'static str,
        /// Its SI unit.
        unit: &'static str,
        /// The left run's value.
        left: f64,
        /// The right run's value.
        right: f64,
    },
    /// The runs reach different distances from zero, so the shared scale is
    /// wider than one of them needs.
    AnomalyRange {
        /// How far the left run's own scale would have reached, in metres.
        left_half_range_m: f64,
        /// How far the right run's own scale would have reached, in metres.
        right_half_range_m: f64,
    },
}

impl fmt::Display for Difference {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Length {
                left_frames,
                right_frames,
                shared_frames,
            } => write!(
                f,
                "Different lengths: {left_frames} frames against {right_frames}. The panels \
                 compare the {shared_frames} frames both runs reach; past that the shorter run \
                 has nothing to show.",
            ),
            Self::SstAnomaly { left, right } => {
                let coupled = if *left { "left" } else { "right" };
                debug_assert_ne!(left, right, "an SST difference is a difference");
                write!(
                    f,
                    "Only the {coupled} run couples SST: its frames carry T' and the other's do \
                     not. Both panels draw the thermocline depth anomaly h, which every run \
                     carries, so the two are comparable in what is on screen.",
                )
            }
            Self::PhysicalParam {
                name,
                unit,
                left,
                right,
            } => write!(f, "{name}: {left} {unit} against {right} {unit}."),
            Self::AnomalyRange {
                left_half_range_m,
                right_half_range_m,
            } => write!(
                f,
                "The runs reach ±{left_half_range_m:.1} m and ±{right_half_range_m:.1} m. Both \
                 panels are drawn on the wider of the two, so the quieter run reads as quieter \
                 rather than being redrawn at full saturation.",
            ),
        }
    }
}

/// Two runs held against each other: the scales both panels draw on, and how
/// far the shared frame index reaches.
///
/// Borrowed rather than owning, and built afresh each repaint: everything in
/// it is a handful of comparisons and two `f64` maxima, so a panel can ask
/// what the shared scale is on every frame of the display without that being
/// work. Nothing here decodes a frame.
#[derive(Debug, Clone, Copy)]
pub struct Comparison<'a> {
    /// The run in the left panel.
    left: &'a LoadedRun,
    /// The run in the right panel.
    right: &'a LoadedRun,
    /// The one colour scale both panels are drawn on.
    scale: DivergingScale,
    /// The one stress scale both overlays are drawn on.
    stress_scale: StressScale,
    /// Frames the comparison covers: the shorter run's count.
    frame_count: u64,
}

impl<'a> Comparison<'a> {
    /// Hold `left` and `right` against each other, or say why they cannot be.
    ///
    /// # Errors
    /// [`Mismatch`] when the side-by-side would claim something the two runs
    /// do not support: a shared grid, or a shared meaning for the frame index.
    pub fn of(left: &'a LoadedRun, right: &'a LoadedRun) -> Result<Self, Mismatch> {
        let (left_header, right_header) = (left.header(), right.header());
        if left_header.grid != right_header.grid {
            return Err(Mismatch::Grid {
                left: left_header.grid,
                right: right_header.grid,
            });
        }
        if left_header.output.interval_s != right_header.output.interval_s {
            return Err(Mismatch::FrameInterval {
                left_s: left_header.output.interval_s,
                right_s: right_header.output.interval_s,
            });
        }
        Ok(Self {
            left,
            right,
            scale: left.anomaly_scale().widened(right.anomaly_scale()),
            stress_scale: left.wind_stress_scale().widened(right.wind_stress_scale()),
            frame_count: left_header
                .output
                .frame_count
                .min(right_header.output.frame_count),
        })
    }

    /// The run in the left panel.
    #[must_use]
    pub const fn left(&self) -> &'a LoadedRun {
        self.left
    }

    /// The run in the right panel.
    #[must_use]
    pub const fn right(&self) -> &'a LoadedRun {
        self.right
    }

    /// The colour scale **both** panels are drawn on: the one that covers
    /// every frame of both runs.
    #[must_use]
    pub const fn scale(&self) -> DivergingScale {
        self.scale
    }

    /// The length scale both wind overlays are drawn on, for the reason
    /// [`Comparison::scale`] is shared.
    #[must_use]
    pub const fn wind_stress_scale(&self) -> StressScale {
        self.stress_scale
    }

    /// Frames the comparison covers, counting from the start of both runs.
    ///
    /// The shorter run's count. A shared index that ran past it would leave one
    /// panel holding its last frame while the other moved on, which draws a
    /// steady ocean neither run produced.
    #[must_use]
    pub const fn frame_count(&self) -> u64 {
        self.frame_count
    }

    /// Model time at frame `index`, in days, counting from the start of the
    /// runs.
    ///
    /// One number for both panels: the cadences match, or there would be no
    /// comparison ([`Mismatch::FrameInterval`]).
    #[must_use]
    pub fn model_time_days(&self, index: u64) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let index = index as f64;
        index * self.left.header().output.interval_s / SECONDS_PER_DAY
    }

    /// What the two runs differ in, in the order the panel states it.
    ///
    /// Empty for two runs of one scenario. Everything here is drawn anyway —
    /// what is refused is in [`Mismatch`] — but a reader comparing two panels
    /// is entitled to know what besides the forcing is not the same.
    #[must_use]
    pub fn differences(&self) -> Vec<Difference> {
        let (left, right) = (self.left.header(), self.right.header());
        let mut differences = Vec::new();
        if left.output.frame_count != right.output.frame_count {
            differences.push(Difference::Length {
                left_frames: left.output.frame_count,
                right_frames: right.output.frame_count,
                shared_frames: self.frame_count,
            });
        }
        differences.extend(differing_params(
            left.physical_params,
            right.physical_params,
        ));
        let (carries_left, carries_right) = (
            left.carries(Variable::SstAnomaly),
            right.carries(Variable::SstAnomaly),
        );
        if carries_left != carries_right {
            differences.push(Difference::SstAnomaly {
                left: carries_left,
                right: carries_right,
            });
        }
        let (left_half_range_m, right_half_range_m) = (
            self.left.anomaly_scale().half_range_m(),
            self.right.anomaly_scale().half_range_m(),
        );
        if left_half_range_m != right_half_range_m {
            differences.push(Difference::AnomalyRange {
                left_half_range_m,
                right_half_range_m,
            });
        }
        differences
    }
}

/// The parameters `left` and `right` were integrated with that are not the
/// same, labelled as [`LoadedRun::metadata`] labels them.
///
/// Compared exactly rather than within a tolerance: these are the numbers the
/// two scenarios stated, not quantities either run computed, so "the same" here
/// means the same value was written down.
fn differing_params(left: PhysicalParams, right: PhysicalParams) -> Vec<Difference> {
    [
        ("Mean depth H", "m", left.mean_depth_m, right.mean_depth_m),
        (
            "Reduced gravity g'",
            "m s^-2",
            left.reduced_gravity_m_per_s2,
            right.reduced_gravity_m_per_s2,
        ),
        (
            "Coriolis gradient β",
            "m^-1 s^-1",
            left.beta_per_m_per_s,
            right.beta_per_m_per_s,
        ),
        (
            "Rayleigh damping r",
            "s^-1",
            left.rayleigh_damping_per_s,
            right.rayleigh_damping_per_s,
        ),
        (
            "Reference density ρ₀",
            "kg m^-3",
            left.reference_density_kg_per_m3,
            right.reference_density_kg_per_m3,
        ),
    ]
    .into_iter()
    .filter(|(_, _, left, right)| left != right)
    .map(|(name, unit, left, right)| Difference::PhysicalParam {
        name,
        unit,
        left,
        right,
    })
    .collect()
}
