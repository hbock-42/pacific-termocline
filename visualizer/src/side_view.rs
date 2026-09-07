//! The equatorial side view: the ocean along the equator drawn side-on, warm
//! water above the thermocline and cold water below it.
//!
//! Every scientific view of this run plots `h` — a depth *anomaly*
//! (`CONTEXT.md`) — and plots it looking down: the basin map colours it, the
//! cross-section draws it as a line about zero. Both show the number the model
//! computes. Neither shows the thing the number describes, and reading either
//! one takes a reader who already knows the thermocline is an interface, that
//! red means deeper, and that a deeper warm layer means warm water piled up.
//!
//! This view turns the picture on its side. The vertical axis is depth below
//! the sea surface, and the interface is drawn at the depth the model actually
//! puts it at — `H + h`, the thickness of the upper layer — so the western
//! wall's +38.2 m of anomaly is an interface 188 m down rather than a red
//! patch. Above it is the warm upper layer; below it the cold abyss. That is
//! the diagram in every textbook of the equatorial Pacific, and it is the
//! same data [`crate::CrossSection`] already extracts — which is why this
//! takes a section rather than a frame. One extraction of the equator, drawn
//! two ways.
//!
//! # What is *not* drawn, and why
//!
//! No sea floor. The 1.5-layer model sits its active layer on a motionless
//! abyss of infinite depth (`CONTEXT.md`, *Upper layer*): the ocean has no
//! bottom in this model, so drawing one would be drawing something the run
//! does not contain. The panel simply ends, with cold water running off it.
//!
//! No relief in the sea surface either. The free surface deviation is not part
//! of the state the model integrates, so the surface is flat across the top of
//! the panel at 0 m — the depth every interface here is measured from.
//!
//! And no water at all where `h` is not a number: the integration diverged
//! there, and an ocean drawn across the gap would claim a depth the run never
//! produced. That is the same refusal the cross-section makes by breaking its
//! line.
//!
//! # Device-free, like every view beside it
//!
//! A [`SideView`] is a list of columns and the coloured bands each one is
//! drawn as, on a unit rectangle with `y` down — no textures, no device, no
//! window. So the acceptance criteria are asserted in `tests/side_view.rs`
//! rather than looked at, and the same code draws the ocean in a browser and
//! natively ([ADR-0006]), exactly as [`crate::heatmap`], [`crate::wind`] and
//! [`crate::cross_section`] do.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use termocline_format::BasinExtent;

use crate::cross_section::CrossSection;
use crate::geography::{coasts, latitude_phrase, longitude_ticks, Coast, LongitudeTick};
use crate::wording::{THERMOCLINE, TRADE_WINDS};

/// The colour of the warm upper layer.
///
/// ColorBrewer's 11-class `RdYlBu`, class 3 (<https://colorbrewer2.org>) — the
/// warm end of the same published, colour-blind-safe family the basin map's
/// ramp comes from, so the two views do not disagree about which end of the
/// ocean is warm.
///
/// It is the *lighter* of the two layers on purpose. Two blocks of colour are
/// told apart by their contrast, and contrast that survives every dichromacy
/// is contrast in luminance rather than in hue; sunlit water above a dark
/// abyss is also the way round a reader already expects. `tests/side_view.rs`
/// holds the pair to WCAG 2.1's 3:1 non-text contrast ratio.
pub const WARM_LAYER_RGB: [u8; 3] = [253, 174, 97];

/// The colour of the cold water beneath the thermocline.
///
/// ColorBrewer's 9-class `Blues`, class 9 (<https://colorbrewer2.org>): the
/// darkest blue of the published family, for the contrast reason
/// [`WARM_LAYER_RGB`] gives.
pub const COLD_LAYER_RGB: [u8; 3] = [8, 48, 107];

/// The colour of the line marking the sea surface across the top of the panel.
///
/// A near-white, so that the top edge of the warm layer reads as the surface
/// of the ocean rather than as the edge of the drawing.
pub const SEA_SURFACE_RGB: [u8; 3] = [247, 247, 247];

/// How far below the deepest interface the run reaches the panel is drawn, as
/// a fraction of that depth.
///
/// The bottom of the panel is a drawing decision and nothing else: the model's
/// abyss is unbounded, so there is no depth the picture *has* to stop at. Half
/// again below the deepest interface any frame of the run reaches leaves cold
/// water plainly visible under the warm pool at every longitude and in every
/// frame, which is what makes the picture read as two layers rather than as
/// one layer with a floor.
const ABYSS_MARGIN_FRACTION: f64 = 0.5;

/// One column of the side view: the ocean at one cell of the equatorial row,
/// and where on the panel it is drawn.
///
/// It carries a depth and a place on the panel and nothing else. Where the
/// column *is* — its longitude, its fraction of the way east — is the
/// cross-section's ([`crate::CrossSectionPoint`]), and the columns here are
/// that section's own points in that section's own order, so restating either
/// would be two records of one fact.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SideViewColumn {
    /// How far the thermocline sits below the sea surface here, in metres:
    /// the total depth `H + h`, not the anomaly (`CONTEXT.md`).
    interface_depth_m: f64,
    /// The same depth as a fraction of the panel's, measured down from the sea
    /// surface — or `None` where there is nothing to draw. Settled when the
    /// view is built, so a column carries its own place and cannot be asked
    /// for it against a panel that is not the one it belongs to.
    interface_fraction: Option<f64>,
}

impl SideViewColumn {
    /// How far the thermocline sits below the sea surface here, in metres.
    ///
    /// The total depth of the upper layer, `H + h` — the depth the model puts
    /// the interface at, and never the anomaly about `H`.
    #[must_use]
    pub const fn interface_depth_m(&self) -> f64 {
        self.interface_depth_m
    }

    /// Where the thermocline sits, as a fraction of the panel's depth measured
    /// down from the sea surface — or `None` where there is nothing to draw.
    ///
    /// `None` in three cases, and they are one case: the run gives no ocean
    /// here that can honestly be drawn. A `NaN` in `h` means the integration
    /// diverged. An interface at or above the sea surface means `h ≤ −H`, so
    /// the upper layer has no thickness left — the linear model has broken
    /// down, and a column drawn full of cold water would report a state the
    /// model cannot be read as describing. And a panel with no depth to it
    /// means a scenario whose mean thermocline depth is not a depth.
    #[must_use]
    pub const fn interface_fraction(&self) -> Option<f64> {
        self.interface_fraction
    }

    /// The two layers of ocean drawn here: the warm upper layer from the sea
    /// surface down to the thermocline, then the cold water filling the panel
    /// beneath it.
    ///
    /// `None` where there is no interface to draw, in which case no water is
    /// drawn in this column at all — see
    /// [`SideViewColumn::interface_fraction`].
    ///
    /// This is the whole of what the picture *is*: a caller turns each band
    /// into a rectangle of whatever size it drew the panel at, and decides
    /// nothing else.
    #[must_use]
    pub fn bands(&self) -> Option<[LayerBand; 2]> {
        let interface = self.interface_fraction?;
        Some([
            LayerBand {
                rgb: WARM_LAYER_RGB,
                top_fraction: 0.0,
                bottom_fraction: interface,
            },
            LayerBand {
                rgb: COLD_LAYER_RGB,
                top_fraction: interface,
                bottom_fraction: 1.0,
            },
        ])
    }
}

/// One coloured band of one column: a layer of the ocean, and the stretch of
/// the panel it fills.
///
/// The fractions run *down* from the sea surface, the convention every view in
/// this crate places things in, so `top_fraction` is always the smaller of the
/// two.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LayerBand {
    /// The colour this layer is drawn in.
    pub rgb: [u8; 3],
    /// Where the band starts, as a fraction of the panel's depth, measured
    /// down from the sea surface.
    pub top_fraction: f64,
    /// Where the band ends, measured the same way.
    pub bottom_fraction: f64,
}

/// The tick spacings the depth axis may use, in metres, coarsest last.
///
/// Round depths a reader counts in. The first that fits the panel in
/// [`MAX_DEPTH_TICKS`] is the one used, so a shallow scenario gets a fine axis
/// and a deep one a coarse axis without either being told which.
const DEPTH_STEPS_M: [f64; 5] = [25.0, 50.0, 100.0, 200.0, 500.0];

/// The most ticks the depth axis is drawn with, counting the sea surface.
///
/// Six across the panel's height leaves the labels a clear gap at the height a
/// window gives the view, and six labels is already more than a reader of a
/// teaching picture needs to place the interface.
const MAX_DEPTH_TICKS: usize = 6;

/// One labelled depth of the vertical axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DepthTick {
    /// The depth below the sea surface, in metres.
    depth_m: f64,
    /// The same depth as a fraction of the panel's, measured down from the sea
    /// surface — the convention every band and interface here is placed in.
    fraction: f64,
}

impl DepthTick {
    /// The depth below the sea surface, in metres.
    #[must_use]
    pub const fn depth_m(&self) -> f64 {
        self.depth_m
    }

    /// Where it sits on the panel, as a fraction of the panel's depth measured
    /// down from the sea surface.
    #[must_use]
    pub const fn fraction(&self) -> f64 {
        self.fraction
    }

    /// The depth as a label: metres, with the unit on it.
    #[must_use]
    pub fn label(&self) -> String {
        format!("{:.0} m", self.depth_m)
    }
}

/// One frame's ocean along the equator, drawn side-on.
#[derive(Debug, Clone)]
pub struct SideView {
    /// One column per cell of the basin's zonal axis, west to east.
    columns: Vec<SideViewColumn>,
    /// The mean thermocline depth `H` the interfaces are measured against, in
    /// metres — a scenario parameter, never a field (`CONTEXT.md`).
    mean_depth_m: f64,
    /// The depth at the bottom edge of the panel, in metres. The run's rather
    /// than this frame's, so scrubbing does not move the ocean.
    deepest_drawn_depth_m: f64,
    /// The latitude the section this was built from was read at, in degrees
    /// north.
    latitude_deg_north: f64,
    /// The basin the section was read across, as the header declares it: what
    /// the longitude axis and the coast names are derived from.
    extent: BasinExtent,
}

impl SideView {
    /// The ocean `section` describes, with its interface at total depth
    /// against a mean thermocline depth of `mean_depth_m`.
    ///
    /// A section rather than a frame because this is not a second extraction
    /// of the equator: it is the cross-section's own data — the same rows
    /// averaged the same way — drawn as an ocean instead of as a line.
    ///
    /// The depth axis comes from the section's scale, which is the *run's*
    /// (`crate::cross_section`), for the reason the basin map's colours are
    /// the run's: a panel rescaled frame by frame would draw every frame's
    /// interface in the same place, and a tilt that collapses over a run —
    /// which is El Niño (`CONTEXT.md`, *ENSO*) — would look like a tilt that
    /// never moved.
    #[must_use]
    pub fn of_section(section: &CrossSection, mean_depth_m: f64) -> Self {
        // The deepest the interface reaches anywhere in the run is `H` plus
        // the far end of the run's symmetric anomaly scale, and the panel goes
        // further still so there is cold water to see beneath it.
        let deepest_interface_m = mean_depth_m + section.scale().half_range_m();
        let deepest_drawn_depth_m = deepest_interface_m * (1.0 + ABYSS_MARGIN_FRACTION);
        let columns = section
            .points()
            .iter()
            .map(|point| {
                let interface_depth_m = mean_depth_m + point.h_m();
                SideViewColumn {
                    interface_depth_m,
                    interface_fraction: interface_fraction(
                        interface_depth_m,
                        deepest_drawn_depth_m,
                    ),
                }
            })
            .collect();
        Self {
            columns,
            mean_depth_m,
            deepest_drawn_depth_m,
            latitude_deg_north: section.latitude_deg_north(),
            extent: section.extent(),
        }
    }

    /// The columns of the ocean, west to east.
    #[must_use]
    pub fn columns(&self) -> &[SideViewColumn] {
        &self.columns
    }

    /// The mean thermocline depth `H` the interfaces are measured against, in
    /// metres.
    #[must_use]
    pub const fn mean_depth_m(&self) -> f64 {
        self.mean_depth_m
    }

    /// The depth at the bottom edge of the panel, in metres.
    ///
    /// Not a sea floor: the model's abyss is unbounded (`CONTEXT.md`, *Upper
    /// layer*), so this is where the picture stops and not where the ocean
    /// does.
    #[must_use]
    pub const fn deepest_drawn_depth_m(&self) -> f64 {
        self.deepest_drawn_depth_m
    }

    /// The latitude the equatorial section behind this view was read at, in
    /// degrees north. Zero for a basin laid out symmetrically about the
    /// equator, which every scenario's is (`CONTEXT.md`, *Basin*).
    #[must_use]
    pub const fn latitude_deg_north(&self) -> f64 {
        self.latitude_deg_north
    }

    /// The basin this ocean was read across, as the run's header declares it.
    #[must_use]
    pub const fn extent(&self) -> BasinExtent {
        self.extent
    }

    /// The two ends of the basin, west first: which coast each is, and where.
    ///
    /// Derived from the header's basin rather than named here, because the
    /// basin is a scenario parameter (`CONTEXT.md`, *Basin*) — see
    /// [`crate::geography`].
    #[must_use]
    pub fn coasts(&self) -> [Coast; 2] {
        coasts(self.extent)
    }

    /// The meridians of the longitude axis under the ocean, west to east.
    #[must_use]
    pub fn longitude_ticks(&self) -> Vec<LongitudeTick> {
        longitude_ticks(self.extent)
    }

    /// The labelled depths of the vertical axis, from the sea surface down.
    ///
    /// Round depths rather than the interface's own: the axis exists so a
    /// reader can read a depth *off* it, which needs numbers they can
    /// interpolate between.
    #[must_use]
    pub fn depth_ticks(&self) -> Vec<DepthTick> {
        let deepest_m = self.deepest_drawn_depth_m;
        if !deepest_m.is_finite() || deepest_m <= 0.0 {
            return Vec::new();
        }
        let Some(step_m) = DEPTH_STEPS_M
            .iter()
            .copied()
            .find(|&step| depth_tick_count(deepest_m, step) <= MAX_DEPTH_TICKS)
        else {
            return Vec::new();
        };
        (0..depth_tick_count(deepest_m, step_m))
            .map(|index| {
                #[allow(clippy::cast_precision_loss)]
                let depth_m = index as f64 * step_m;
                DepthTick {
                    depth_m,
                    fraction: depth_m / deepest_m,
                }
            })
            .collect()
    }

    /// What the picture is, in one line a visitor can read without knowing any
    /// of this project's vocabulary.
    ///
    /// It names the thermocline and explains it in the same breath: the word
    /// is what the whole view is about, so it is glossed rather than avoided
    /// (`crate::wording`).
    #[must_use]
    pub fn caption(&self) -> String {
        format!(
            "The ocean {}, seen from the side: warm water on top, cold water below, and the {} \
             between them at the depth the model puts it",
            latitude_phrase(self.latitude_deg_north),
            THERMOCLINE.glossed()
        )
    }

    /// What the vertical axis measures, and the one thing about it a reader
    /// could otherwise get wrong: the foot of the panel is where the picture
    /// stops, not where the ocean does.
    ///
    /// Stated per view because in a comparison the two runs may declare
    /// different mean depths, and so be drawn to different axes.
    #[must_use]
    pub fn depth_axis_note(&self) -> String {
        format!(
            "Depth in metres below the sea surface: 0 m at the surface, {:.0} m at the foot of \
             the panel — the model's deep ocean has no floor",
            self.deepest_drawn_depth_m
        )
    }

    /// What the wind has done to the ocean in this picture, in plain words and
    /// with the number it is read off the picture with.
    ///
    /// Measured rather than asserted: the tilt is the difference between the
    /// two ends of *this frame's* drawn interface, so the note cannot claim a
    /// warm pool piled in the west while the panel shows a level ocean. The
    /// easterly trade winds are what piles it (`CONTEXT.md`, *Thermocline
    /// tilt*), and they are named — and explained — wherever the ocean is
    /// tilted the way they tilt it.
    ///
    /// A frame whose ends are not both drawn — the model has broken down at
    /// one of them (`SideViewColumn::interface_fraction`) — gets no note at
    /// all rather than a tilt measured across a gap.
    #[must_use]
    pub fn tilt_note(&self) -> Option<String> {
        let (west, east) = (self.columns.first()?, self.columns.last()?);
        // Both ends have to be drawn: a tilt measured across a column where
        // the model has broken down is not a tilt of anything.
        west.interface_fraction()?;
        east.interface_fraction()?;
        let drop_m = west.interface_depth_m() - east.interface_depth_m();
        // Rounded to the metre first, because the metre is what the note
        // states: a difference that prints as "0 m deeper" is a level ocean as
        // far as this sentence is concerned.
        let stated_m = drop_m.round();
        if stated_m > 0.0 {
            return Some(format!(
                "The {} push the warm water westward: here it is piled up against the western \
                 coast, where the boundary is {stated_m:.0} m deeper than at the eastern coast",
                TRADE_WINDS.glossed()
            ));
        }
        if stated_m < 0.0 {
            return Some(format!(
                "Here the boundary is {:.0} m deeper at the eastern coast than at the western — \
                 the warm water is not piled up in the west",
                stated_m.abs()
            ));
        }
        Some(format!(
            "The {} push the warm water westward; in this frame the boundary is level from coast \
             to coast",
            TRADE_WINDS.glossed()
        ))
    }

    /// The width of one column, as a fraction of the panel: the columns tile
    /// the panel edge to edge, so this is one over their count.
    #[must_use]
    pub fn column_width_fraction(&self) -> f64 {
        if self.columns.is_empty() {
            return 0.0;
        }
        #[allow(clippy::cast_precision_loss)]
        let count = self.columns.len() as f64;
        1.0 / count
    }
}

/// How many round depths of `step_m` fit on a panel `deepest_m` deep,
/// counting the sea surface at zero.
fn depth_tick_count(deepest_m: f64, step_m: f64) -> usize {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let count = (deepest_m / step_m).floor() as usize + 1;
    count
}

/// Where an interface `interface_depth_m` below the surface sits on a panel
/// `deepest_drawn_depth_m` deep, or `None` where there is nothing honest to
/// draw — the three refusals [`SideViewColumn::interface_fraction`] states.
///
/// A finite interface deeper than the panel is clamped to its foot. That cannot
/// happen for a panel built over the run being drawn — the panel is half again
/// the deepest interface the run reaches — but it can for a caller that mixed
/// two runs, and an interface at the edge of the picture is a better answer
/// than one drawn beyond it.
fn interface_fraction(interface_depth_m: f64, deepest_drawn_depth_m: f64) -> Option<f64> {
    if !interface_depth_m.is_finite()
        || interface_depth_m <= 0.0
        || !deepest_drawn_depth_m.is_finite()
        || deepest_drawn_depth_m <= 0.0
    {
        return None;
    }
    Some((interface_depth_m / deepest_drawn_depth_m).min(1.0))
}
