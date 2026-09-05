//! The shell itself: a window (or a canvas) that loads a run, says what it is,
//! and draws one frame of it.
//!
//! It draws the header first — the grid, the scenario, the frame count —
//! because that is what tells a reader the run they think they opened is the
//! run they opened, and under it the basin map of one chosen frame.
//!
//! Everything with a value in it lives in [`crate::run`], [`crate::heatmap`],
//! [`crate::wind`], [`crate::cross_section`], [`crate::time_series`] and
//! [`crate::pending`]; this module is the part that needs a
//! GPU, and so is deliberately thin. What it adds on top of them is a texture
//! cache and a layout, and neither is where a wrong basin map would come from.

use egui::{Color32, RichText};

use termocline_format::GridSpec;

use crate::comparison::Side;
use crate::loading::Loaded;
use crate::run::SECONDS_PER_DAY;
use crate::{
    BasinPoint, Comparison, CrossSection, DivergingScale, Heatmap, LoadedRun, Loader, Mismatch,
    PendingRun, Playback, PointSeries, Scrubber, StressScale, WindOverlay,
};

/// What one panel is showing.
enum Shown {
    /// No run yet: the panel explains how to load one.
    Nothing,
    /// A load is in flight from the named source.
    Loading(String),
    /// A run, with its metadata.
    Run(Box<LoadedRun>),
    /// A source that did not yield a run.
    Failed {
        /// Where the run was to come from.
        source: String,
        /// What went wrong, in the words of whatever refused it.
        message: String,
    },
}

impl Shown {
    /// The run this panel is showing, if it is showing one.
    const fn run(&self) -> Option<&LoadedRun> {
        match self {
            Self::Run(run) => Some(run),
            _ => None,
        }
    }
}

/// The visualizer's application state.
///
/// Everything a panel needs is held per [`Side`], and a single run is the left
/// panel with the right one closed (T-09.5). One state rather than two — a
/// "single run" mode beside a "comparison" mode — because the two would drift:
/// the frame chooser, the drop handling and the loader would each have to be
/// taught the difference, and each is a place the two panels could come to
/// disagree.
pub struct VisualizerApp {
    /// What each panel is showing.
    shown: [Shown; 2],
    /// Dropped files seen so far, waiting for their pair, per panel.
    pending: [PendingRun; 2],
    /// The one channel every source of run bytes posts to, whichever panel
    /// asked for it.
    loader: Loader,
    /// The URL a run is served under, as typed or as passed in the query, per
    /// panel.
    run_url: [String; 2],
    /// Whether the second panel is open.
    comparing: bool,
    /// Which panel a dropped run is loaded into.
    ///
    /// Chosen by the reader rather than inferred from which panel is empty: a
    /// run arrives as two files, sometimes in two drops, and a target that
    /// moved between them would split one run across both panels.
    drop_side: Side,
    /// The frame chooser both panels share, and what each of them last drew.
    basin_map: BasinMap,
}

impl Default for VisualizerApp {
    fn default() -> Self {
        Self {
            shown: [Shown::Nothing, Shown::Nothing],
            pending: [PendingRun::default(), PendingRun::default()],
            loader: Loader::default(),
            run_url: [String::new(), String::new()],
            comparing: false,
            drop_side: Side::Left,
            basin_map: BasinMap::default(),
        }
    }
}

impl VisualizerApp {
    /// A shell with nothing loaded.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Start fetching the run served under `base_url` into the first panel,
    /// showing it as loading until it lands.
    pub fn fetch_run(&mut self, base_url: &str, ctx: &egui::Context) {
        self.fetch_into(Side::Left, base_url, ctx);
    }

    /// Start fetching the run served under `base_url` into the second panel,
    /// opening the comparison.
    pub fn fetch_run_to_compare(&mut self, base_url: &str, ctx: &egui::Context) {
        self.comparing = true;
        self.drop_side = Side::Right;
        self.fetch_into(Side::Right, base_url, ctx);
    }

    /// Start fetching the run served under `base_url` into `side`.
    fn fetch_into(&mut self, side: Side, base_url: &str, ctx: &egui::Context) {
        self.run_url[side.index()] = base_url.to_owned();
        self.shown[side.index()] = Shown::Loading(base_url.to_owned());
        let ctx = ctx.clone();
        self.loader
            .fetch(side, base_url, move || ctx.request_repaint());
    }

    /// Load the run in `directory` into the first panel. Native only: a
    /// browser has no directories to open (ADR-0006).
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load_directory(&mut self, directory: &std::path::Path) {
        self.load_directory_into(Side::Left, directory);
    }

    /// Load the run in `directory` into the second panel, opening the
    /// comparison. Native only, for the reason [`VisualizerApp::load_directory`]
    /// is.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load_directory_to_compare(&mut self, directory: &std::path::Path) {
        self.comparing = true;
        self.drop_side = Side::Right;
        self.load_directory_into(Side::Right, directory);
    }

    /// Load the run in `directory` into `side`.
    #[cfg(not(target_arch = "wasm32"))]
    fn load_directory_into(&mut self, side: Side, directory: &std::path::Path) {
        let source = directory.display().to_string();
        self.shown[side.index()] = Shown::Loading(source.clone());
        self.loader.deliver(
            side,
            source,
            crate::loading::native::read_run_directory(directory),
        );
    }

    /// Ask the user for a run directory for `side`, on a thread of its own so
    /// the frame loop keeps running while the dialog is open.
    #[cfg(not(target_arch = "wasm32"))]
    fn pick_directory(&self, side: Side, ctx: &egui::Context) {
        let sender = self.loader.sender();
        let ctx = ctx.clone();
        let title = format!("Open a run directory for panel {}", side.label());
        std::thread::spawn(move || {
            let Some(directory) = rfd::FileDialog::new().set_title(title).pick_folder() else {
                return;
            };
            let _ = sender.send(Loaded {
                side,
                source: directory.display().to_string(),
                bytes: crate::loading::native::read_run_directory(&directory),
            });
            ctx.request_repaint();
        });
    }

    /// Take whatever finished loading since the last frame and show it in the
    /// panel that asked for it.
    fn absorb_finished_loads(&mut self) {
        while let Some(Loaded {
            side,
            source,
            bytes,
        }) = self.loader.poll()
        {
            // Whatever arrives, the map of the last run is no longer the map
            // of the run on screen — and in a comparison the scale is both
            // runs', so the other panel's map is no longer its map either.
            self.basin_map.forget();
            self.shown[side.index()] = match bytes.and_then(|bytes| {
                LoadedRun::from_bytes(source.clone(), bytes).map_err(|error| error.to_string())
            }) {
                Ok(run) => Shown::Run(Box::new(run)),
                Err(message) => Shown::Failed { source, message },
            };
        }
    }

    /// Take the files dropped on the window this frame, and load the run into
    /// the chosen panel once both of them have arrived.
    fn absorb_dropped_files(&mut self, ctx: &egui::Context) {
        let dropped = ctx.input(|input| input.raw.dropped_files.clone());
        let side = if self.comparing {
            self.drop_side
        } else {
            Side::Left
        };
        for file in &dropped {
            #[cfg(not(target_arch = "wasm32"))]
            if let Some(path) = file.path.as_ref().filter(|path| path.is_dir()) {
                self.load_directory_into(side, path);
                continue;
            }
            match dropped_file_contents(file) {
                Ok((name, bytes)) => {
                    if !self.pending[side.index()].offer(&name, bytes) {
                        self.shown[side.index()] = Shown::Failed {
                            source: name.clone(),
                            message: format!(
                                "{name} is not part of a run; drop {} and {}",
                                termocline_format::HEADER_FILE_NAME,
                                termocline_format::FRAME_FILE_NAME
                            ),
                        };
                    }
                }
                Err(message) => {
                    self.shown[side.index()] = Shown::Failed {
                        source: "dropped file".to_owned(),
                        message,
                    };
                }
            }
        }
        if let Some(bytes) = self.pending[side.index()].take_run() {
            self.loader.deliver(side, "dropped files", Ok(bytes));
        }
    }

    /// The bar of run-loading affordances: one row per open panel, and the
    /// switch that opens the second.
    ///
    /// Returns whether a run-URL field has the keyboard. They are the one
    /// thing in the shell that a keystroke means something different to, so
    /// they are what decide whether the scrubber's keys are the scrubber's to
    /// take.
    fn draw_controls(&mut self, ui: &mut egui::Ui) -> bool {
        let mut url_has_keyboard = false;
        for side in Side::BOTH {
            if side == Side::Right && !self.comparing {
                continue;
            }
            url_has_keyboard |= self.draw_side_controls(ui, side);
        }
        ui.checkbox(&mut self.comparing, "Compare two runs")
            .on_hover_text(
                "Show a second run beside this one, on one frame index and one colour scale",
            );
        url_has_keyboard
    }

    /// One panel's row of loading affordances.
    fn draw_side_controls(&mut self, ui: &mut egui::Ui, side: Side) -> bool {
        ui.horizontal(|ui| {
            if self.comparing {
                ui.radio_value(&mut self.drop_side, side, format!("Run {}", side.label()))
                    .on_hover_text("Dropped files load into this panel");
            }
            #[cfg(not(target_arch = "wasm32"))]
            if ui.button("Open run directory…").clicked() {
                self.pick_directory(side, ui.ctx());
            }
            ui.label("Run URL:");
            let url = ui.add(
                egui::TextEdit::singleline(&mut self.run_url[side.index()])
                    .hint_text("https://…/run-demo/")
                    .desired_width(260.0),
            );
            let submitted =
                url.lost_focus() && ui.input(|input| input.key_pressed(egui::Key::Enter));
            let fetch = ui.button("Fetch").clicked();
            if (submitted || fetch) && !self.run_url[side.index()].trim().is_empty() {
                let (url, ctx) = (self.run_url[side.index()].clone(), ui.ctx().clone());
                self.fetch_into(side, &url, &ctx);
            }
            url.has_focus()
        })
        .inner
    }
}

/// What to show when a panel has no run yet, or its last one failed.
fn draw_instructions(ui: &mut egui::Ui, pending: &PendingRun) {
    ui.label(format!(
        "Drop a run's {} and {} onto this window{}.",
        termocline_format::HEADER_FILE_NAME,
        termocline_format::FRAME_FILE_NAME,
        if cfg!(target_arch = "wasm32") {
            ""
        } else {
            ", or drop the run directory itself"
        }
    ));
    let still_needed = pending.still_needed();
    if still_needed.len() == 1 {
        ui.label(format!("Waiting for {}.", still_needed[0]));
    }
}

impl eframe::App for VisualizerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.absorb_finished_loads();
        self.absorb_dropped_files(ctx);

        let url_has_keyboard = egui::TopBottomPanel::top("controls")
            .show(ctx, |ui| {
                ui.add_space(4.0);
                ui.heading(crate::APP_NAME);
                let url_has_keyboard = self.draw_controls(ui);
                ui.add_space(4.0);
                url_has_keyboard
            })
            .inner;
        let keyboard_free = !url_has_keyboard;

        // Disjoint borrows: the panels read their runs while the maps they
        // draw cache a texture of one frame each.
        let Self {
            shown,
            pending,
            basin_map,
            comparing,
            ..
        } = self;
        egui::CentralPanel::default().show(ctx, |ui| {
            if *comparing {
                draw_comparison(ui, shown, pending, basin_map, keyboard_free);
            } else {
                draw_single(ui, &shown[0], &pending[0], basin_map, keyboard_free);
            }
        });

        if ctx.input(|input| !input.raw.hovered_files.is_empty()) {
            // Dropping is the web's only local affordance, so say the window
            // will take the file before the user lets go of it.
            egui::Area::new(egui::Id::new("drop-hint"))
                .anchor(egui::Align2::CENTER_BOTTOM, egui::vec2(0.0, -24.0))
                .show(ctx, |ui| {
                    ui.label(RichText::new("Drop to load").strong());
                });
        }
    }
}

/// The one open panel: its metadata, and under it the basin map of the chosen
/// frame.
fn draw_single(
    ui: &mut egui::Ui,
    shown: &Shown,
    pending: &PendingRun,
    basin_map: &mut BasinMap,
    keyboard_free: bool,
) {
    match shown {
        Shown::Run(run) => {
            draw_run_header(ui, Side::Left, run);
            ui.add_space(12.0);
            ui.separator();
            basin_map.draw(ui, run, keyboard_free);
        }
        other => draw_waiting(ui, other, pending),
    }
}

/// Two runs held against each other, or why they are not.
///
/// The two panels are only worth drawing together once both runs are there and
/// [`Comparison`] has accepted them; short of that this says what each panel is
/// waiting for, side by side, so a reader can see which half is missing.
fn draw_comparison(
    ui: &mut egui::Ui,
    shown: &[Shown; 2],
    pending: &[PendingRun; 2],
    basin_map: &mut BasinMap,
    keyboard_free: bool,
) {
    let (Some(left), Some(right)) = (shown[0].run(), shown[1].run()) else {
        ui.columns(2, |columns| {
            for side in Side::BOTH {
                let column = &mut columns[side.index()];
                column.label(RichText::new(format!("Run {}", side.label())).strong());
                match &shown[side.index()] {
                    Shown::Run(run) => draw_run_header(column, side, run),
                    other => draw_waiting(column, other, &pending[side.index()]),
                }
            }
        });
        return;
    };
    match Comparison::of(left, right) {
        Ok(comparison) => basin_map.draw_comparison(ui, &comparison, keyboard_free),
        Err(mismatch) => draw_refusal(ui, [left, right], &mismatch),
    }
}

/// Why these two runs are not being drawn side by side, and what each of them
/// is.
///
/// A refusal rather than a picture: the alternative is two panels claiming a
/// correspondence — same place, same moment — that the runs do not support,
/// and a reader has no way to see that from the picture itself
/// (`crate::comparison`).
fn draw_refusal(ui: &mut egui::Ui, runs: [&LoadedRun; 2], mismatch: &Mismatch) {
    ui.label(
        RichText::new("These two runs cannot be compared")
            .color(Color32::LIGHT_RED)
            .strong(),
    );
    ui.label(mismatch.to_string());
    ui.label("Load another run into one of the panels, or turn the comparison off to look at them one at a time.");
    ui.add_space(12.0);
    ui.separator();
    ui.columns(2, |columns| {
        for side in Side::BOTH {
            draw_run_header(&mut columns[side.index()], side, runs[side.index()]);
        }
    });
}

/// A panel with no run in it yet: what it is doing, or how to give it one.
fn draw_waiting(ui: &mut egui::Ui, shown: &Shown, pending: &PendingRun) {
    match shown {
        Shown::Nothing => draw_instructions(ui, pending),
        Shown::Loading(source) => {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label(format!("Loading {source}…"));
            });
        }
        Shown::Failed { source, message } => {
            ui.label(
                RichText::new(format!("{source} could not be loaded"))
                    .color(Color32::LIGHT_RED)
                    .strong(),
            );
            ui.label(message);
            ui.separator();
            draw_instructions(ui, pending);
        }
        Shown::Run(_) => {}
    }
}

/// Where a run came from, and the metadata panel of what it says it is.
///
/// The grid is named after the panel rather than after the run: two panels may
/// hold the same run, and two grids under one name would lay out as one.
fn draw_run_header(ui: &mut egui::Ui, side: Side, run: &LoadedRun) {
    ui.label(RichText::new(run.source()).strong());
    ui.add_space(6.0);
    egui::Grid::new(format!("run-metadata-{}", side.label()))
        .num_columns(2)
        .striped(true)
        .spacing([24.0, 4.0])
        .show(ui, |ui| {
            for row in run.metadata() {
                ui.label(row.label);
                ui.label(row.value);
                ui.end_row();
            }
        });
}

/// The name and bytes of a dropped file, from whichever of the two egui fills
/// in: a path natively, the bytes themselves in a browser.
fn dropped_file_contents(file: &egui::DroppedFile) -> Result<(String, Vec<u8>), String> {
    #[cfg(not(target_arch = "wasm32"))]
    if let Some(path) = &file.path {
        let bytes = std::fs::read(path)
            .map_err(|error| format!("{} could not be read: {error}", path.display()))?;
        return Ok((path.display().to_string(), bytes));
    }
    let bytes = file
        .bytes
        .as_ref()
        .ok_or_else(|| format!("{} arrived with no contents", file.name))?;
    Ok((file.name.clone(), bytes.to_vec()))
}

/// Height of the colour bar, in points.
const COLOR_BAR_HEIGHT: f32 = 14.0;

/// Samples across the colour bar. More than the eye can separate at the width
/// a panel gives it, and few enough that building it is not worth caching
/// beyond the frame it belongs to.
const COLOR_BAR_SAMPLES: usize = 256;

/// Width of the arrow line, in points.
const ARROW_WIDTH_PT: f32 = 1.2;

/// Anything drawn *over* the basin map — the wind arrows, and the ring round
/// the cell the time series is of: near-black, so it reads over the pale
/// middle of the colour scale and over the casing at the scale's two dark ends.
const OVERLAY_COLOR: Color32 = Color32::from_rgb(16, 16, 16);

/// The casing under it: opaque white, so the near-black mark stays separable
/// from the dark blue and dark red the scale ends on.
const OVERLAY_CASING_COLOR: Color32 = Color32::from_rgb(245, 245, 245);

/// Height of the cross-section chart, in points.
///
/// Tall enough that the ±`half_range` axis has room to show a tilt changing by
/// a few per cent, and short enough that the basin map above it still gets the
/// bulk of a window: the chart says how much, the map says where.
const CROSS_SECTION_HEIGHT_PT: f32 = 120.0;

/// Width of a line on either chart, in points.
const CHART_LINE_WIDTH_PT: f32 = 1.6;

/// The line either chart draws `h` in: the deep end of the basin map's colour
/// scale, so the two charts and the map are read as one picture of the same
/// field.
const H_LINE_COLOR: Color32 = Color32::from_rgb(178, 24, 43);

/// A chart's frame and its zero line: grey, so the line drawn over them is
/// what the eye lands on. Shared by the cross-section and the time series, so
/// the two read as the same kind of picture.
const CHART_AXIS_COLOR: Color32 = Color32::from_gray(128);

/// Height of the point time-series chart, in points.
///
/// The same height as the cross-section: the two charts answer the same
/// question on different axes — how much anomaly, against distance and against
/// time — and drawing one taller than the other would say one mattered more.
const TIME_SERIES_HEIGHT_PT: f32 = 120.0;

/// The `T'` line: ColorBrewer `PRGn`'s dark purple
/// (<https://colorbrewer2.org>), chosen because it is off the red-blue axis the
/// map and the `h` line already use, so the second line cannot be mistaken for
/// another reading of the first.
const TIME_SERIES_SST_COLOR: Color32 = Color32::from_rgb(118, 42, 131);

/// The marker naming the frame the map is showing: darker than the axis, so it
/// reads as an annotation on the chart rather than part of its frame.
const TIME_SERIES_MARKER_COLOR: Color32 = Color32::from_gray(90);

/// Width of the ring drawn round the cell the time series is of, in points,
/// and of the pale casing under it.
const SELECTED_CELL_WIDTH_PT: f32 = 1.5;
const SELECTED_CELL_CASING_WIDTH_PT: f32 = 3.5;

/// The smallest the ring is drawn, in points.
///
/// A cell of the scenario basin is a fraction of a point wide on any window a
/// browser tab is likely to give it, and a ring that small would be invisible.
/// Nine points is about the size of a text caret: findable without hiding the
/// field around it.
const SELECTED_CELL_MIN_PT: f32 = 9.0;

/// Width of the pale casing drawn under each arrow, in points.
///
/// The map underneath runs from a dark blue through a near-white to a dark red,
/// so no single colour reads over all of it. Each arrow is drawn twice — a
/// wider pale stroke, then a narrow dark one — which is the cartographic casing
/// that keeps a line legible over any ground.
const ARROW_CASING_WIDTH_PT: f32 = 3.0;

/// The frame chooser, and the basin map each open panel draws of the frame it
/// names.
///
/// The map is a texture rather than a mesh: `h` is one value per cell, and the
/// cheapest honest way to show a cell grid is one pixel per cell, magnified
/// without interpolation. It is also the way that costs the same on both
/// targets, which is what ADR-0006 asks of anything drawn here.
///
/// # One chooser, two panels
///
/// A comparison (T-09.5) does not *keep* the two panels synced; there is
/// nothing to keep. The chooser, the clock, the layer toggles, the colour bar
/// and the picked cell are held here, once, and each panel holds only what it
/// built for itself — so scrubbing or playing writes one `u64` that both
/// panels read, and there is no state in which they could disagree. That is
/// the same argument `crate::playback` makes for playback owning no index, one
/// level up.
///
/// # What a drag costs
///
/// The scrubber is dragged, so everything under it is on a path that runs once
/// per frame of the *display*, not once per frame of the run. Three things
/// were made not to happen there (T-08.3): reaching the run's frame walks no
/// other frames ([`LoadedRun::frame`]); a repaint that lands on the frame
/// already drawn rebuilds nothing ([`Attempt`]); and the colour bar, which is
/// the run's and not the frame's, is uploaded once per run ([`ColorBar`]).
/// What is left is one frame decoded, colour-mapped and uploaded per frame the
/// reader actually asks for — one per open panel — which is the work the drag
/// is for.
///
/// Playback (T-09.2) adds nothing to that: it writes the scrubber's index and
/// nothing else, so a played frame and a dragged one cost the same, and a
/// paused run — which is how a run loads — asks for no repaint at all.
///
/// The wind overlay rides on that same [`Attempt`]: it is built with the frame
/// and drawn from geometry, so neither showing it nor hiding it is a reason to
/// rebuild anything, and it adds nothing per frame a drag passes through. The
/// equatorial cross-section (T-09.3) rides on it the same way, and is drawn
/// from the same scale as the colour bar, so it too costs nothing per repaint
/// and nothing per toggle.
///
/// # What the time series costs, and why it is not on that path
///
/// The point time series (T-09.4) is the one view whose cost is shaped
/// differently, because it is the one view that is not of a frame: it is one
/// cell of *every* frame, so the indexed lookup that makes the others cheap
/// buys it nothing and a rebuild walks the whole run — 731 frames of the
/// scenario run, against the one frame everything above reads.
///
/// It is therefore held beside the frame cache rather than inside it, and
/// rebuilt only when the reader picks a **different** cell. Scrubbing,
/// playing, toggling any chart and re-picking the selected cell all leave it
/// alone; [`SeriesCache::walks`] counts the rebuilds so the tests can say so
/// by name. In a comparison the picked cell is the pair's — the two runs are
/// over one grid, or they would not be comparable — so one click asks each
/// panel for the series of the same place in its own run, and each panel walks
/// its own run once for it.
struct BasinMap {
    /// The frame the reader has chosen, and the ways they choose another.
    /// Shared by both panels of a comparison.
    scrubber: Scrubber,
    /// The clock that chooses frames on the reader's behalf, for both panels.
    playback: Playback,
    /// Which layers are drawn over and under the maps.
    layers: Layers,
    /// The cell the reader picked, in the coordinates of the grid both panels
    /// are over.
    selected: Option<BasinPoint>,
    /// What each panel last drew, and the series it holds. The second is
    /// untouched until a comparison is open.
    panels: [FramePanel; 2],
    /// The colour bar of whatever is being drawn, built the first time a frame
    /// of it is drawn: the run's scale on its own, and both runs' shared scale
    /// in a comparison.
    bar: Option<ColorBar>,
}

/// Which layers the reader has asked for, over and under every open map.
///
/// One set for the pair rather than one each: they are a question about what a
/// reader wants to see, not about a run, and two panels showing different
/// layers would be two pictures rather than a comparison.
#[derive(Debug, Clone, Copy)]
struct Layers {
    /// Whether the wind-stress overlay is drawn over the maps.
    wind: bool,
    /// Whether the equatorial cross-section is drawn under the maps.
    section: bool,
    /// Whether the point time series is drawn under the maps.
    series: bool,
}

impl Default for BasinMap {
    fn default() -> Self {
        Self {
            scrubber: Scrubber::default(),
            // Paused, as `crate::playback` says: a run that started moving on
            // load would be off its first frame before the header had been
            // read, and it is also what leaves an idle repaint with nothing to
            // rebuild.
            playback: Playback::new(),
            layers: Layers {
                // On by default: the forcing is why the map looks the way it
                // does, and a reader who does not know the overlay exists
                // cannot ask for it.
                wind: true,
                // On by default: the tilt is what the run is about, and the
                // section is the view that states it as a number rather than
                // as a colour (T-09.3).
                section: true,
                // On by default, with nothing selected: the chart says how to
                // pick a point, and a reader who does not know the map is
                // clickable never finds out (T-09.4).
                series: true,
            },
            selected: None,
            panels: [FramePanel::default(), FramePanel::default()],
            bar: None,
        }
    }
}

/// One panel's caches: the last frame it drew, and the series it holds.
#[derive(Default)]
struct FramePanel {
    /// The last attempt at drawing a frame, if there has been one.
    attempt: Option<Attempt>,
    /// The series at the picked cell of this panel's run, and what building it
    /// has cost.
    series: SeriesCache,
}

/// The point time series on screen, and how many times one has been built.
///
/// The two travel together everywhere, because the count is only meaningful
/// against the series it counts: this is the one path in the shell that walks a
/// whole run ([`crate::time_series`]), and the count is how a test says so.
#[derive(Default)]
struct SeriesCache {
    /// The series at the cell the reader picked, if they have picked one.
    series: Option<PointSeries>,
    /// How many times a series has been built since this map was created.
    ///
    /// Not shown anywhere. It exists because "the run is walked once per cell
    /// the reader picks" is the cost property this view lives or dies by, and a
    /// property no test can name is a property the next change quietly breaks.
    walks: usize,
}

impl SeriesCache {
    /// The series on screen, if the reader has picked a cell.
    const fn shown(&self) -> Option<&PointSeries> {
        self.series.as_ref()
    }

    /// Forget the cell, which was a place in a run that has gone. The count is
    /// the panel's own tally and survives, as it must to mean anything.
    fn forget(&mut self) {
        self.series = None;
    }

    /// Show the series at `point`, walking `run` for it unless the series
    /// already in hand is that cell's.
    ///
    /// The guard is on the cell rather than on the click: a reader clicking
    /// twice in the same cell — or dragging across one — has asked for the
    /// series they already have.
    fn select(&mut self, run: &LoadedRun, point: BasinPoint) {
        if self
            .series
            .as_ref()
            .is_some_and(|held| held.point() == point)
        {
            return;
        }
        self.series = Some(
            PointSeries::at_point(run, point)
                .expect("a point picked off this run's own map is a cell of its basin"),
        );
        self.walks += 1;
    }
}

/// What came of trying to draw one frame.
///
/// Kept because a panel repaints many times per second while the chosen frame
/// does not change, and each repaint would otherwise decode and upload the
/// frame again. The failure is kept for the same reason: retrying it every
/// repaint would fail every repaint.
struct Attempt {
    /// The frame this was an attempt at.
    index: u64,
    /// The colour scale it was drawn on. Part of what identifies the attempt
    /// because in a comparison the scale is both runs' (`crate::comparison`):
    /// loading a louder run beside this one changes what this panel's own
    /// frame should look like, and a cache that ignored that would leave one
    /// panel drawn on a scale the colour bar no longer states.
    scale: DivergingScale,
    /// The stress scale its overlay was drawn on, kept for the same reason.
    stress_scale: StressScale,
    /// The map of it, or what refused to build one, in the words of whatever
    /// refused it.
    outcome: Result<DrawnFrame, String>,
}

/// A frame already colour-mapped and uploaded, and the overlay that goes over
/// it.
struct DrawnFrame {
    /// Its model time, in seconds since the start of the run.
    t_s: f64,
    /// The map itself, one texel per cell.
    map: egui::TextureHandle,
    /// The frame's wind stress, as arrows in the map's own cell coordinates.
    ///
    /// Built whether or not it is currently shown, so that toggling the overlay
    /// neither rebuilds the texture under it nor re-decodes the frame — and so
    /// that there is no path at all from the toggle to the map.
    wind: WindOverlay,
    /// The frame's `h` along the equator, in the same unit rectangle the chart
    /// is drawn into.
    ///
    /// Built with the frame, like the overlay and for the same reason: the
    /// section is the frame's, so a repaint that lands on the frame already
    /// drawn rebuilds nothing, and toggling the chart cannot reach the map.
    section: CrossSection,
}

/// The colour bar of what is on screen, and the scale it was sampled from.
///
/// The scale is the run's rather than the frame's (`crate::heatmap`), and in a
/// comparison it is both runs' (`crate::comparison`) — so the bar is the same
/// in every frame, and there is one of it for two panels rather than one each.
struct ColorBar {
    /// The scale the bar was sampled from, so a scale that is not this one
    /// gets its own bar.
    scale: DivergingScale,
    /// The bar itself, sampled across the scale.
    texture: egui::TextureHandle,
}

impl FramePanel {
    /// Forget the run this was a panel of.
    fn forget(&mut self) {
        self.attempt = None;
        self.series.forget();
    }

    /// The frame `index` of `run` drawn on `scale`, from the cache where the
    /// last attempt was at the same frame on the same scales.
    fn drawn(
        &mut self,
        ui: &egui::Ui,
        run: &LoadedRun,
        index: u64,
        scale: DivergingScale,
        stress_scale: StressScale,
    ) -> &Result<DrawnFrame, String> {
        drawn_in(&mut self.attempt, ui, run, index, scale, stress_scale)
    }

    /// Draw this panel: the map of frame `index` of `run`, and whichever
    /// layers the reader has asked for.
    ///
    /// The picked cell is the caller's rather than this panel's, because in a
    /// comparison it is the pair's: a click here picks the same place in both
    /// panels, and each panel then holds that place's series in its own run.
    fn draw(
        &mut self,
        ui: &mut egui::Ui,
        run: &LoadedRun,
        index: u64,
        scale: DivergingScale,
        stress_scale: StressScale,
        layers: Layers,
        selected: &mut Option<BasinPoint>,
    ) {
        // Disjoint borrows: the frame already drawn is read while the series —
        // which is of the whole run rather than of that frame — is picked and
        // built beside it.
        let Self { attempt, series } = self;
        let drawn = match drawn_in(attempt, ui, run, index, scale, stress_scale) {
            Ok(drawn) => drawn,
            Err(message) => {
                ui.label(
                    RichText::new(format!("Frame {} could not be drawn", index + 1))
                        .color(Color32::LIGHT_RED)
                        .strong(),
                );
                ui.label(message);
                return;
            }
        };
        ui.horizontal(|ui| {
            ui.label("west");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label("east");
            });
        });
        let reserved_pt = reserved_below_map_pt(ui, layers.section, layers.series);
        let map_area = draw_texture_fitted(ui, &drawn.map, reserved_pt);
        let map = map_area.rect;
        if layers.wind {
            draw_wind_arrows(ui, map, drawn);
        }
        // The click that picks a place, and the only path in the shell that
        // walks a whole run. It is taken before the charts are drawn, so the
        // series on screen is the one the click just asked for. A repaint is
        // asked for with it because the other panel of a comparison may
        // already have been drawn this pass, and it holds the same cell of its
        // own run.
        if let Some(point) = clicked_point(&map_area, map, run.header().grid) {
            *selected = Some(point);
            ui.ctx().request_repaint();
        }
        if let Some(point) = *selected {
            series.select(run, point);
        }
        // Where the series came from, marked on the map it was picked off: a
        // chart labelled with a longitude says where in words, and the ring
        // says where by pointing at it.
        if let Some(shown) = series.shown() {
            draw_selected_cell(ui, map, shown.point());
        }
        // Directly under the map and across exactly its width, so a longitude
        // on the chart sits under the column of the map it came from. The
        // colour bar goes below both: the scale is the same one.
        if layers.section {
            draw_cross_section(ui, &drawn.section, map);
        }
        // The time series goes under the section, and only its horizontal
        // extent is the map's: its axis is time, so nothing on it lines up with
        // a column of the map.
        if layers.series {
            draw_time_series(ui, series.shown(), drawn.t_s, map);
        }
    }
}

impl BasinMap {
    /// Forget the runs this was a map of. The layer toggles are the reader's
    /// choice rather than the run's, so they survive.
    fn forget(&mut self) {
        self.scrubber = Scrubber::new();
        // The clock was playing the run that has just gone; the speed it was
        // playing at is the reader's choice, so that stays.
        self.playback.pause();
        // The cell was a place in the basin of the run that has just gone; the
        // next run's cell of that index is a different place.
        self.selected = None;
        for panel in &mut self.panels {
            panel.forget();
        }
        self.bar = None;
    }

    /// Draw the frame chooser and the map of the chosen frame of one run.
    fn draw(&mut self, ui: &mut egui::Ui, run: &LoadedRun, keyboard_free: bool) {
        let frame_count = run.header().output.frame_count;
        let Some(index) = self.chooser(ui, frame_count, keyboard_free) else {
            ui.label("This run holds no frames to draw.");
            return;
        };
        let (scale, stress_scale) = (run.anomaly_scale(), run.wind_stress_scale());
        // Before the frame is read back, because that read borrows the panel
        // for as long as the frame is drawn and this is the one thing left
        // that writes to it.
        self.ensure_color_bar(ui, scale);
        let Self {
            panels,
            bar,
            layers,
            selected,
            ..
        } = self;
        let panel = &mut panels[Side::Left.index()];
        // The frame is built before the line above it is written, because the
        // time in that line is the frame's. Building it here costs nothing:
        // the panel keeps what it built, and drawing it reads that back.
        let Some(t_s) = panel
            .drawn(ui, run, index, scale, stress_scale)
            .as_ref()
            .ok()
            .map(|drawn| drawn.t_s)
        else {
            // The panel says what stopped it, where the map would have been.
            panel.draw(ui, run, index, scale, stress_scale, *layers, selected);
            return;
        };
        // Counted from one, as the metadata panel counts the run's frames.
        // The scrubber's own index starts at zero, and it shows no number.
        ui.label(format!(
            "Frame {} of {frame_count} — thermocline depth anomaly h at {:.2} days",
            index + 1,
            t_s / SECONDS_PER_DAY
        ));
        panel.draw(ui, run, index, scale, stress_scale, *layers, selected);
        draw_color_bar(ui, bar.as_ref().expect("a colour bar was just built"));
        if layers.wind {
            ui.label(format!(
                "Wind stress τ: the longest arrow is {:.3} N m^-2, the strongest in the run",
                stress_scale.max_magnitude_pa()
            ));
        }
    }

    /// Draw the frame chooser and, side by side, the map each of the two runs
    /// has of the frame it names.
    ///
    /// One chooser, one clock, one colour bar and one picked cell for the pair:
    /// the two panels are drawn from the same index and the same scale, which
    /// is what makes them comparable at all (`crate::comparison`).
    fn draw_comparison(
        &mut self,
        ui: &mut egui::Ui,
        comparison: &Comparison<'_>,
        keyboard_free: bool,
    ) {
        let frame_count = comparison.frame_count();
        let Some(index) = self.chooser(ui, frame_count, keyboard_free) else {
            ui.label("These runs share no frames to draw.");
            return;
        };
        let (scale, stress_scale) = (comparison.scale(), comparison.wind_stress_scale());
        self.ensure_color_bar(ui, scale);
        let Self {
            panels,
            bar,
            layers,
            selected,
            ..
        } = self;
        let runs = [comparison.left(), comparison.right()];
        // The frames are built before the line above them is written, because
        // what that line says about the time is the frames' to say and not the
        // cadence's. Building them here costs nothing: each panel keeps what it
        // built, and the column below reads it back rather than building it
        // again.
        let mut times_s = [None, None];
        for side in Side::BOTH {
            times_s[side.index()] = panels[side.index()]
                .drawn(ui, runs[side.index()], index, scale, stress_scale)
                .as_ref()
                .ok()
                .map(|drawn| drawn.t_s);
        }
        ui.label(frame_line(index, frame_count, times_s));
        ui.columns(2, |columns| {
            for (side, (panel, run)) in Side::BOTH.into_iter().zip(panels.iter_mut().zip(runs)) {
                let column = &mut columns[side.index()];
                column.label(RichText::new(run.source()).strong());
                panel.draw(column, run, index, scale, stress_scale, *layers, selected);
            }
        });
        // Under both panels and across the whole width, because it is the one
        // scale both of them were drawn on.
        draw_color_bar(ui, bar.as_ref().expect("a colour bar was just built"));
        if layers.wind {
            ui.label(format!(
                "Wind stress τ: the longest arrow is {:.3} N m^-2, the strongest in either run",
                stress_scale.max_magnitude_pa()
            ));
        }
        for difference in comparison.differences() {
            ui.label(difference.to_string());
        }
    }

    /// Draw the frame chooser over a run of `frame_count` frames and say which
    /// frame it names, or `None` when there are none to choose between.
    fn chooser(&mut self, ui: &mut egui::Ui, frame_count: u64, keyboard_free: bool) -> Option<u64> {
        self.scrubber.fit_to(frame_count);
        self.scrubber.last()?;
        // The clock before the controls: the frame this repaint draws is the
        // one whatever time has passed since the last one has bought, and
        // every affordance below writes the same index it does.
        self.playback
            .advance(&mut self.scrubber, f64::from(ui.input(|i| i.stable_dt)));
        self.scrubber.draw(ui, keyboard_free);
        self.playback.draw(ui, &mut self.scrubber, keyboard_free);
        // The layers are drawn over and under the map, never into it: nothing
        // the map is built from depends on these.
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.layers.wind, "Wind stress τ");
            ui.checkbox(&mut self.layers.section, "Equatorial cross-section");
            ui.checkbox(&mut self.layers.series, "Point time series");
        });
        Some(self.scrubber.index())
    }

    /// Sample the colour bar of `scale`, the first time a frame is drawn on it.
    fn ensure_color_bar(&mut self, ui: &egui::Ui, scale: DivergingScale) {
        if self.bar.as_ref().is_none_or(|bar| bar.scale != scale) {
            self.bar = Some(ColorBar {
                scale,
                texture: ui.ctx().load_texture(
                    "basin-map-scale",
                    egui::ColorImage::from_rgb(
                        [COLOR_BAR_SAMPLES, 1],
                        &scale.bar_rgb(COLOR_BAR_SAMPLES),
                    ),
                    egui::TextureOptions::LINEAR,
                ),
            });
        }
    }
}

/// Frame `index` of `run` on `scale`, from `attempt` where that is what it
/// last held, and built into it where it is not.
///
/// Free rather than a method because a panel draws its series while holding
/// the frame it drew, and the two are different fields of it.
fn drawn_in<'a>(
    attempt: &'a mut Option<Attempt>,
    ui: &egui::Ui,
    run: &LoadedRun,
    index: u64,
    scale: DivergingScale,
    stress_scale: StressScale,
) -> &'a Result<DrawnFrame, String> {
    let stale = attempt.as_ref().is_none_or(|last| {
        last.index != index || last.scale != scale || last.stress_scale != stress_scale
    });
    if stale {
        *attempt = Some(Attempt {
            index,
            scale,
            stress_scale,
            outcome: build(ui, run, index, scale, stress_scale),
        });
    }
    &attempt.as_ref().expect("an attempt was just made").outcome
}

/// The line above the two panels: which frame they are both on, and what the
/// runs themselves date it.
///
/// The time is read off the frames rather than counted from the output
/// cadence, so two runs that date the same index differently say so instead of
/// being labelled with one time neither of them wrote. Compared exactly: two
/// runs put the same number on the same moment or they do not, and no
/// tolerance would make a disagreement of a second more honest than one of a
/// day.
fn frame_line(index: u64, frame_count: u64, times_s: [Option<f64>; 2]) -> String {
    let frame = format!("Frame {} of {frame_count} in both panels", index + 1);
    let days = |t_s: f64| t_s / SECONDS_PER_DAY;
    match times_s {
        [Some(left_s), Some(right_s)] if left_s == right_s => format!(
            "{frame} — thermocline depth anomaly h at {:.2} days",
            days(left_s)
        ),
        [Some(left_s), Some(right_s)] => format!(
            "{frame} — thermocline depth anomaly h, which run {} dates {:.2} days and run {} \
             dates {:.2} days",
            Side::Left.label(),
            days(left_s),
            Side::Right.label(),
            days(right_s),
        ),
        // A panel that could not draw its frame has no time to report; the
        // panel itself says what stopped it.
        _ => frame,
    }
}

/// Colour-map frame `index` of `run` on `scale` and upload it, or say what
/// stopped that.
fn build(
    ui: &egui::Ui,
    run: &LoadedRun,
    index: u64,
    scale: DivergingScale,
    stress_scale: StressScale,
) -> Result<DrawnFrame, String> {
    let frame = run
        .frame(index)
        .ok_or_else(|| format!("this run holds no frame {index}"))?;
    let heatmap =
        Heatmap::of_frame(run.header().grid, &frame, scale).map_err(|error| error.to_string())?;
    let wind = WindOverlay::of_frame(run.header().grid, &frame, stress_scale)
        .map_err(|error| error.to_string())?;
    let section = CrossSection::of_frame(run.header().grid, &frame, scale)
        .map_err(|error| error.to_string())?;
    let image = egui::ColorImage::from_rgb([heatmap.width(), heatmap.height()], heatmap.rgb());
    Ok(DrawnFrame {
        t_s: frame.t_s(),
        wind,
        section,
        // Nearest, not linear: a texel is a cell of the model, and
        // smoothing between them would draw an anomaly the run never
        // produced.
        map: ui
            .ctx()
            .load_texture("basin-map", image, egui::TextureOptions::NEAREST),
    })
}

/// Draw the colour bar and the anomalies its ends stand for.
fn draw_color_bar(ui: &mut egui::Ui, bar: &ColorBar) {
    let width = ui.available_width();
    ui.add(egui::Image::new(egui::load::SizedTexture::new(
        bar.texture.id(),
        egui::vec2(width, COLOR_BAR_HEIGHT),
    )));
    let half_range_m = bar.scale.half_range_m();
    // Negating zero would label a run at rest "-0.0 m".
    let shallow_m = if half_range_m == 0.0 {
        0.0
    } else {
        -half_range_m
    };
    // Three equal columns, so the middle label sits under the middle of the
    // bar — where the neutral colour is — rather than wherever the two end
    // labels happen to leave room.
    ui.columns(3, |columns| {
        columns[0].label(format!("{shallow_m:+.1} m (shallower)"));
        columns[1].with_layout(egui::Layout::top_down(egui::Align::Center), |ui| {
            ui.label("0 m")
        });
        columns[2].with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(format!("{half_range_m:+.1} m (deeper)"));
        });
    });
}

/// Draw the equatorial cross-section: `h` along the equator against longitude.
///
/// The chart's coordinates are the section's own — a unit rectangle, `y` down
/// — so placing a point is one multiplication by however large the chart was
/// drawn, and the module that knows what `h` is never learns what a panel is
/// ([`crate::cross_section`]).
///
/// It is drawn across `map`, the rectangle the basin map landed in, rather
/// than across the panel: the map is fitted to its own shape and so is
/// narrower than the panel on a short window, and a chart wider than the map
/// would put a longitude under the wrong column of it. The two views are only
/// worth having together if a place on one is the same place on the other.
///
/// The line breaks where a point has no position: a `NaN` in `h` means the
/// integration diverged there, and a segment drawn across the gap would claim
/// a value the run never produced.
fn draw_cross_section(ui: &mut egui::Ui, section: &CrossSection, map: egui::Rect) {
    ui.label(format!(
        "Thermocline depth anomaly h {}",
        section_latitude(section)
    ));
    let (row, _response) = ui.allocate_exact_size(
        egui::vec2(ui.available_width(), CROSS_SECTION_HEIGHT_PT),
        egui::Sense::hover(),
    );
    let chart = egui::Rect::from_min_max(
        egui::pos2(map.left(), row.top()),
        egui::pos2(map.right(), row.bottom()),
    );
    let painter = ui.painter().with_clip_rect(chart);
    let axis = egui::Stroke::new(1.0_f32, CHART_AXIS_COLOR);
    painter.rect_stroke(chart, 0.0, axis, egui::StrokeKind::Inside);
    // The zero line, at the middle of a scale that is symmetric about zero: it
    // is where the thermocline sits at its mean depth, and where the tilt
    // changes sign.
    let zero_y = chart.center().y;
    painter.line_segment(
        [
            egui::pos2(chart.left(), zero_y),
            egui::pos2(chart.right(), zero_y),
        ],
        axis,
    );

    draw_chart_line(
        &painter,
        chart,
        section
            .points()
            .iter()
            .map(|point| section.plot_position(point)),
        H_LINE_COLOR,
    );

    let half_range_m = section.scale().half_range_m();
    ui.columns(3, |columns| {
        columns[0].label(format!(
            "{} (west) — the axis reaches ±{half_range_m:.1} m",
            longitude_label(section.points().first())
        ));
        columns[1].with_layout(egui::Layout::top_down(egui::Align::Center), |ui| {
            ui.label("0 m")
        });
        columns[2].with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(format!(
                "{} (east)",
                longitude_label(section.points().last())
            ));
        });
    });
}

/// Draw the point time series: `h` at the chosen cell against model time, with
/// the run's SST anomaly beside it where the run has one.
///
/// The chart's coordinates are the series' own — a unit rectangle, `y` down —
/// so placing a sample is one multiplication by however large the chart was
/// drawn, and the module that knows what `h` is never learns what a panel is
/// ([`crate::time_series`]).
///
/// It spans the width of `map` for the same reason the cross-section does:
/// the two charts and the map are read as one picture, and a chart wider than
/// the map would sit under nothing. Unlike the section, though, only the width
/// is shared — the axis here is time, so no point on it lines up with a column
/// of the map.
///
/// The line breaks where a sample has no position, exactly as the section's
/// does: a `NaN` means the integration diverged there, and a segment drawn
/// across the gap would claim a value the run never produced.
fn draw_time_series(
    ui: &mut egui::Ui,
    series: Option<&PointSeries>,
    shown_t_s: f64,
    map: egui::Rect,
) {
    ui.label(match series {
        Some(series) => format!(
            "Thermocline depth anomaly h at {} through the run — the axis reaches ±{:.1} m",
            place_label(series.point()),
            series.scale().half_range_m()
        ),
        None => "Point time series — click the basin map to pick a place".to_owned(),
    });
    let (row, _response) = ui.allocate_exact_size(
        egui::vec2(ui.available_width(), TIME_SERIES_HEIGHT_PT),
        egui::Sense::hover(),
    );
    let chart = egui::Rect::from_min_max(
        egui::pos2(map.left(), row.top()),
        egui::pos2(map.right(), row.bottom()),
    );
    let painter = ui.painter().with_clip_rect(chart);
    let axis = egui::Stroke::new(1.0_f32, CHART_AXIS_COLOR);
    painter.rect_stroke(chart, 0.0, axis, egui::StrokeKind::Inside);
    // The zero line, at the middle of an axis symmetric about zero: it is where
    // the thermocline sits at its mean depth, and — for the second line — where
    // the mixed layer sits at its climatological temperature.
    let zero_y = chart.center().y;
    painter.line_segment(
        [
            egui::pos2(chart.left(), zero_y),
            egui::pos2(chart.right(), zero_y),
        ],
        axis,
    );

    let Some(series) = series else {
        ui.label("No place chosen yet.");
        return;
    };
    draw_chart_line(
        &painter,
        chart,
        series.samples().iter().map(|s| series.plot_position(s)),
        H_LINE_COLOR,
    );
    if series.carries_sst_anomaly() {
        draw_chart_line(
            &painter,
            chart,
            series.samples().iter().map(|s| series.sst_plot_position(s)),
            TIME_SERIES_SST_COLOR,
        );
    }
    // Where on the series the frame the map is showing sits, so the two views
    // say which instant they have in common.
    #[allow(clippy::cast_possible_truncation)]
    let shown_x = series
        .time_fraction(shown_t_s)
        .mul_add(f64::from(chart.width()), f64::from(chart.left())) as f32;
    painter.line_segment(
        [
            egui::pos2(shown_x, chart.top()),
            egui::pos2(shown_x, chart.bottom()),
        ],
        egui::Stroke::new(1.0_f32, TIME_SERIES_MARKER_COLOR),
    );

    let first_days = series.samples().first().map_or(0.0, |s| s.t_s()) / SECONDS_PER_DAY;
    let last_days = series.samples().last().map_or(0.0, |s| s.t_s()) / SECONDS_PER_DAY;
    ui.columns(3, |columns| {
        columns[0].label(format!("{first_days:.0} days"));
        columns[1].with_layout(egui::Layout::top_down(egui::Align::Center), |ui| {
            ui.label(sst_label(series))
        });
        columns[2].with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(format!("{last_days:.0} days"));
        });
    });
}

/// Draw one line of a chart from `positions` — points on a unit rectangle, `y`
/// down — breaking it wherever a point has no position rather than drawing
/// across the gap.
///
/// A `NaN` in a field means the integration diverged there, and a segment drawn
/// across the gap would claim a value the run never produced. Shared by the
/// cross-section and the time series: the two differ in what their horizontal
/// axis means, not in how a line is drawn.
fn draw_chart_line(
    painter: &egui::Painter,
    chart: egui::Rect,
    positions: impl Iterator<Item = Option<(f64, f64)>>,
    color: Color32,
) {
    #[allow(clippy::cast_possible_truncation)]
    let to_screen = |(x, y): (f64, f64)| {
        chart.min + egui::vec2(x as f32 * chart.width(), y as f32 * chart.height())
    };
    let flush = |points: &mut Vec<egui::Pos2>| {
        if points.len() > 1 {
            painter.add(egui::Shape::line(
                std::mem::take(points),
                egui::Stroke::new(CHART_LINE_WIDTH_PT, color),
            ));
        } else {
            points.clear();
        }
    };
    let mut segment: Vec<egui::Pos2> = Vec::new();
    for position in positions {
        match position {
            Some(position) => segment.push(to_screen(position)),
            None => flush(&mut segment),
        }
    }
    flush(&mut segment);
}

/// What the chart says about the run's SST anomaly.
///
/// A run that never coupled SST is told so in those words. It is not drawn as a
/// flat line at zero: zero kelvin of anomaly is a claim that the mixed layer
/// sat at its climatological temperature all run, and an uncoupled run made no
/// such claim (`termocline_format::Frame`).
fn sst_label(series: &PointSeries) -> String {
    match series.sst_scale() {
        Some(scale) => format!(
            "SST anomaly T' also drawn — its axis reaches ±{:.2} K",
            scale.half_range_k()
        ),
        None => "This run carries no SST anomaly".to_owned(),
    }
}

/// A cell of the basin, as a chart labels it: a longitude and a latitude, each
/// a magnitude and a hemisphere, never a minus sign.
fn place_label(point: BasinPoint) -> String {
    let (deg_east, deg_north) = (point.longitude_deg_east(), point.latitude_deg_north());
    format!(
        "{:.1}°{}, {:.1}°{}",
        deg_east.abs(),
        if deg_east < 0.0 { 'W' } else { 'E' },
        deg_north.abs(),
        if deg_north < 0.0 { 'S' } else { 'N' }
    )
}

/// Where the section was read, in words.
///
/// Every scenario's basin is laid out symmetrically about the equator
/// (`CONTEXT.md`, *Basin*), and then the two innermost rows straddle it and
/// their mean is the equator itself. A basin that is not says so rather than
/// having its off-equator line labelled as the equator.
fn section_latitude(section: &CrossSection) -> String {
    let latitude_deg_north = section.latitude_deg_north();
    let read_at = if latitude_deg_north == 0.0 {
        "along the equator".to_owned()
    } else {
        format!(
            "along {:.2}°{}",
            latitude_deg_north.abs(),
            if latitude_deg_north < 0.0 { 'S' } else { 'N' }
        )
    };
    match section.rows_averaged() {
        1 => format!("{read_at}, on the row of cells there"),
        rows => format!("{read_at}, the mean of the {rows} rows nearest it"),
    }
}

/// A point's longitude, as a chart labels it: a magnitude and a hemisphere,
/// never a minus sign.
fn longitude_label(point: Option<&crate::CrossSectionPoint>) -> String {
    let Some(point) = point else {
        return String::new();
    };
    let deg_east = point.longitude_deg_east();
    format!(
        "{:.1}°{}",
        deg_east.abs(),
        if deg_east < 0.0 { 'W' } else { 'E' }
    )
}

/// Draw `texture` as large as the panel allows without changing its shape, and
/// hand back the response: where it landed, so a layer can be drawn over it,
/// and whether it was clicked, which is how a place is picked (T-09.4).
///
/// The basin is far wider than it is tall, so fitting it to the width alone
/// would push the colour bar off the bottom of a short window.
fn draw_texture_fitted(
    ui: &mut egui::Ui,
    texture: &egui::TextureHandle,
    reserved_pt: f32,
) -> egui::Response {
    let size = texture.size_vec2();
    let available = egui::vec2(
        ui.available_width(),
        (ui.available_height() - COLOR_BAR_HEIGHT * 3.0 - reserved_pt).max(COLOR_BAR_HEIGHT),
    );
    let scale = (available.x / size.x).min(available.y / size.y);
    ui.add(
        egui::Image::new(egui::load::SizedTexture::new(texture.id(), size * scale))
            .sense(egui::Sense::click()),
    )
}

/// Height the charts under the map need, in points: each is a chart and the
/// two rows of text around it.
///
/// It comes out of the height the map would otherwise have taken, so that
/// turning a chart on shrinks the map rather than pushing the colour bar off a
/// short window.
fn reserved_below_map_pt(ui: &egui::Ui, show_section: bool, show_series: bool) -> f32 {
    let text_pt = ui.text_style_height(&egui::TextStyle::Body) * 2.0;
    let section_pt = if show_section {
        CROSS_SECTION_HEIGHT_PT + text_pt
    } else {
        0.0
    };
    let series_pt = if show_series {
        TIME_SERIES_HEIGHT_PT + text_pt
    } else {
        0.0
    };
    section_pt + series_pt
}

/// Draw a ring round the cell `point` names on the map occupying `map`.
///
/// The map's own coordinates, like the wind arrows over it: the cell is a
/// column and a row of the field, and the field's shape is what `point` carries
/// for exactly this. Nothing here reads or writes the texture under it.
///
/// It is cased the way the arrows are, and for the same reason: the ground it
/// is drawn over runs from a dark blue through a near-white to a dark red, and
/// no single colour reads over all of it.
fn draw_selected_cell(ui: &egui::Ui, map: egui::Rect, point: BasinPoint) {
    let (width, height) = point.field_shape();
    if width == 0 || height == 0 {
        return;
    }
    #[allow(clippy::cast_precision_loss)]
    let per_cell = egui::vec2(map.width() / width as f32, map.height() / height as f32);
    // Row 0 of the map is the northernmost, which is the *last* row of the
    // field (`crate::heatmap`).
    #[allow(clippy::cast_precision_loss)]
    let centre = map.min
        + egui::vec2(
            (point.column() as f32 + 0.5) * per_cell.x,
            ((height - 1 - point.row()) as f32 + 0.5) * per_cell.y,
        );
    let half = egui::vec2(
        (per_cell.x / 2.0).max(SELECTED_CELL_MIN_PT / 2.0),
        (per_cell.y / 2.0).max(SELECTED_CELL_MIN_PT / 2.0),
    );
    let ring = egui::Rect::from_min_max(centre - half, centre + half);
    let painter = ui.painter().with_clip_rect(map);
    for (width_pt, color) in [
        (SELECTED_CELL_CASING_WIDTH_PT, OVERLAY_CASING_COLOR),
        (SELECTED_CELL_WIDTH_PT, OVERLAY_COLOR),
    ] {
        painter.rect_stroke(
            ring,
            0.0,
            egui::Stroke::new(width_pt, color),
            egui::StrokeKind::Outside,
        );
    }
}

/// The cell of `grid` a click on the basin map landed in, or `None` if this
/// was not a click, or was a click outside the map.
///
/// The whole of the click-to-select affordance that a mouse is needed for: the
/// arithmetic that turns a place on the map into a cell of the basin is
/// [`BasinPoint::at_map_fraction`], where it can be asserted without one.
fn clicked_point(response: &egui::Response, map: egui::Rect, grid: GridSpec) -> Option<BasinPoint> {
    if !response.clicked() || map.width() <= 0.0 || map.height() <= 0.0 {
        return None;
    }
    let position = response.interact_pointer_pos()?;
    BasinPoint::at_map_fraction(
        grid,
        f64::from((position.x - map.left()) / map.width()),
        f64::from((position.y - map.top()) / map.height()),
    )
}

/// Draw the frame's wind stress as arrows over the map occupying `map`.
///
/// The overlay's coordinates are the map's own — cells from its northwest
/// corner — so placing an arrow is one multiplication by however large the map
/// was drawn. Nothing here reads or writes the texture under it.
fn draw_wind_arrows(ui: &egui::Ui, map: egui::Rect, drawn: &DrawnFrame) {
    let cells = drawn.map.size_vec2();
    if cells.x <= 0.0 || cells.y <= 0.0 {
        return;
    }
    let per_cell = egui::vec2(map.width() / cells.x, map.height() / cells.y);
    #[allow(clippy::cast_possible_truncation)]
    let to_screen = |(x_cells, y_cells): (f64, f64)| {
        map.min + egui::vec2(x_cells as f32 * per_cell.x, y_cells as f32 * per_cell.y)
    };
    let painter = ui.painter().with_clip_rect(map);
    for arrow in drawn.wind.arrows() {
        let tail = to_screen(arrow.tail_cells());
        let along = to_screen(arrow.tip_cells()) - tail;
        for (width, color) in [
            (ARROW_CASING_WIDTH_PT, OVERLAY_CASING_COLOR),
            (ARROW_WIDTH_PT, OVERLAY_COLOR),
        ] {
            painter.arrow(tail, along, egui::Stroke::new(width, color));
        }
    }
}

#[cfg(test)]
mod tests {
    //! What a drag must not do per repaint, asserted on the panel itself.
    //!
    //! `egui` needs no GPU to lay a panel out and no window to run in, so the
    //! caches [`BasinMap`] keeps — the frame it already drew, and the colour
    //! bar of the run it is drawing — are checked here by texture identity: a
    //! rebuilt texture is a new handle, and a reused one is the same handle.
    //!
    //! The run is written from `termocline_format` alone, as the integration
    //! tests write theirs (`tests/common/mod.rs`).

    use termocline_format::{
        frame_encoding, BasinExtent, Frame, GridSpec, OutputTiming, PhysicalParams, RunHeader,
        Variable,
    };

    use super::{BasinMap, DrawnFrame, LoadedRun, Side};
    use crate::{BasinPoint, Comparison, RunBytes, SeriesSample};

    /// A basin small enough to build in a unit test, on the extent of
    /// `CONTEXT.md`, *Basin*.
    fn grid() -> GridSpec {
        GridSpec::new(4, 3, BasinExtent::new(120.0, -80.0, -25.0, 25.0))
            .expect("a 4x3 basin is a valid grid")
    }

    /// A run of three frames whose `h` is everywhere the frame's own index, in
    /// metres, so consecutive frames are drawn in different colours.
    fn run() -> LoadedRun {
        run_of("basin-map", 3, 1.0)
    }

    /// A run named `source` of `frame_count` frames whose `h` is everywhere
    /// `metres_per_frame` times the frame's own index, in metres.
    ///
    /// The value doubles as the frame's name: a panel drawing this run says
    /// which frame it is drawing through the numbers in it, which is what lets
    /// two panels be asserted to be on the same frame rather than merely to
    /// hold the same index.
    fn run_of(source: &str, frame_count: u64, metres_per_frame: f64) -> LoadedRun {
        let grid = grid();
        let header = RunHeader::new(
            grid,
            PhysicalParams {
                mean_depth_m: 150.0,
                reduced_gravity_m_per_s2: 0.06,
                beta_per_m_per_s: 2.3e-11,
                rayleigh_damping_per_s: 1.0e-7,
                reference_density_kg_per_m3: 1025.0,
            },
            source,
            OutputTiming {
                frame_count,
                interval_s: 86_400.0,
            },
        );
        let mut frames = Vec::new();
        for index in 0..header.output.frame_count {
            #[allow(clippy::cast_precision_loss)]
            let value = index as f64 * metres_per_frame;
            let field = |variable| vec![0.0; grid.field_len(variable)];
            #[allow(clippy::cast_precision_loss)]
            let t_s = index as f64 * header.output.interval_s;
            let frame = Frame::new(
                t_s,
                &grid,
                vec![value; grid.field_len(Variable::ThermoclineDepthAnomaly)],
                field(Variable::ZonalCurrentAnomaly),
                field(Variable::MeridionalCurrentAnomaly),
                field(Variable::ZonalWindStress),
                field(Variable::MeridionalWindStress),
            )
            .expect("fields sized from the grid fit it");
            frames.extend(
                bincode::serde::encode_to_vec(&frame, frame_encoding()).expect("a frame encodes"),
            );
        }
        LoadedRun::from_bytes(
            source,
            RunBytes {
                header: serde_json::to_vec(&header).expect("a header serializes"),
                frames,
            },
        )
        .expect("a run written from its own header loads")
    }

    /// Repaint `map` once, and say which textures it drew the run with.
    fn repaint(
        ctx: &egui::Context,
        map: &mut BasinMap,
        run: &LoadedRun,
    ) -> (egui::TextureId, egui::TextureId) {
        // The full output of the pass is what a backend would paint; what is
        // under test is which textures the panel asked for, not the pixels.
        let _painted = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| map.draw(ui, run, true));
        });
        let drawn = panel_frame(map, Side::Left);
        let bar = map.bar.as_ref().expect("a colour bar was built");
        (drawn.map.id(), bar.texture.id())
    }

    /// What panel `side` drew, on the repaint that has just happened.
    fn panel_frame(map: &BasinMap, side: Side) -> &DrawnFrame {
        map.panels[side.index()]
            .attempt
            .as_ref()
            .expect("a frame was drawn")
            .outcome
            .as_ref()
            .expect("it drew")
    }

    /// The `h` panel `side` is drawing, in metres.
    ///
    /// Every cell of a frame of these runs holds the same value, and that
    /// value names the frame ([`run_of`]) — so this is which frame the panel
    /// is showing, read off what it drew rather than off the index it was
    /// handed.
    fn panel_h_m(map: &BasinMap, side: Side) -> f64 {
        let section = &panel_frame(map, side).section;
        let first = section
            .points()
            .first()
            .expect("the section crosses the basin");
        assert!(
            section
                .points()
                .iter()
                .all(|point| point.h_m() == first.h_m()),
            "these runs are flat within a frame"
        );
        first.h_m()
    }

    /// The series panel `side` holds, and how many times it has walked a run
    /// for one.
    fn panel_series(map: &BasinMap, side: Side) -> &super::SeriesCache {
        &map.panels[side.index()].series
    }

    /// Repaint a comparison of `left` and `right` once.
    fn repaint_comparison(
        ctx: &egui::Context,
        map: &mut BasinMap,
        left: &LoadedRun,
        right: &LoadedRun,
    ) {
        let comparison = Comparison::of(left, right).expect("these runs are comparable");
        let _painted = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                map.draw_comparison(ui, &comparison, true);
            });
        });
    }

    #[test]
    fn repainting_the_same_frame_rebuilds_nothing() {
        let (ctx, run) = (egui::Context::default(), run());
        let mut map = BasinMap::default();
        let first = repaint(&ctx, &mut map, &run);
        assert_eq!(repaint(&ctx, &mut map, &run), first);
    }

    #[test]
    fn playing_carries_the_panel_to_the_end_of_the_run_and_stops_there() {
        let (ctx, run) = (egui::Context::default(), run());
        let mut map = BasinMap::default();
        let (first_map, first_bar) = repaint(&ctx, &mut map, &run);
        map.playback.set_frames_per_second(60.0);
        map.playback.play(&mut map.scrubber);
        // Sixty frames a second and a repaint loop that runs at rather less
        // than sixty hertz: ten repaints is far more than the three frames of
        // this run, so what is asserted is where it stopped, not how fast.
        let mut last = first_map;
        for _ in 0..10 {
            let (drawn, bar) = repaint(&ctx, &mut map, &run);
            last = drawn;
            // The scale is the run's, so playing through it rebuilds no bar.
            assert_eq!(bar, first_bar);
        }
        assert_eq!(
            map.scrubber.index(),
            map.scrubber.last().expect("this run has frames"),
            "playback stops at the last frame"
        );
        assert!(!map.playback.is_playing());
        assert_ne!(last, first_map, "and the panel is drawing it");
    }

    #[test]
    fn toggling_the_cross_section_rebuilds_nothing() {
        // The section is built with the frame and drawn from its own geometry,
        // so there is no path from the toggle to the texture under it — the
        // same structural claim T-09.1 makes for the wind overlay.
        let (ctx, run) = (egui::Context::default(), run());
        let mut map = BasinMap::default();
        let first = repaint(&ctx, &mut map, &run);
        map.layers.section = !map.layers.section;
        assert_eq!(repaint(&ctx, &mut map, &run), first);
        map.layers.section = !map.layers.section;
        assert_eq!(repaint(&ctx, &mut map, &run), first);
    }

    #[test]
    fn the_cross_section_shows_the_frame_the_map_shows() {
        // The one thing the ticket asks of the two views together: the section
        // is of the frame the map draws, whichever frame that is. This run's
        // `h` is everywhere the frame's own index in metres, so the section's
        // own values say which frame it was built from.
        let (ctx, run) = (egui::Context::default(), run());
        let mut map = BasinMap::default();
        for index in 0..run.header().output.frame_count {
            map.scrubber.set_index(index);
            let _ = repaint(&ctx, &mut map, &run);
            let drawn = panel_frame(&map, Side::Left);
            #[allow(clippy::cast_precision_loss)]
            let expected_m = index as f64;
            assert!(drawn
                .section
                .points()
                .iter()
                .all(|point| point.h_m() == expected_m));
        }
    }

    /// The cell a reader would pick out of the middle of this run's basin.
    fn middle_point(run: &LoadedRun) -> BasinPoint {
        BasinPoint::at_map_fraction(run.header().grid, 0.5, 0.5)
            .expect("the middle of the map is on the map")
    }

    #[test]
    fn picking_a_point_walks_the_run_once_and_picking_it_again_not_at_all() {
        // The cost property this view lives by. Every other view reads one
        // frame; a series reads all of them, so the walk must happen once per
        // cell the reader picks and not once per click, per drag step or per
        // repaint.
        let (ctx, run) = (egui::Context::default(), run());
        let mut map = BasinMap::default();
        let _ = repaint(&ctx, &mut map, &run);
        map.panels[Side::Left.index()].series.select(&run, middle_point(&run));
        assert_eq!(panel_series(&map, Side::Left).walks, 1);
        map.panels[Side::Left.index()].series.select(&run, middle_point(&run));
        assert_eq!(panel_series(&map, Side::Left).walks, 1, "the same cell is the series in hand");
        // And it is a series of this run: `h` is everywhere the frame's own
        // index in metres, so the samples say which run they came from.
        let series = panel_series(&map, Side::Left).shown().expect("a point was picked");
        #[allow(clippy::cast_precision_loss)]
        let expected: Vec<f64> = (0..run.frame_count()).map(|index| index as f64).collect();
        assert_eq!(
            series
                .samples()
                .iter()
                .map(SeriesSample::h_m)
                .collect::<Vec<f64>>(),
            expected
        );
    }

    #[test]
    fn scrubbing_and_playing_do_not_rebuild_the_time_series() {
        // The series is of the whole run, so the frame on screen is not a
        // reason to rebuild it — and playback, which changes that frame every
        // repaint, is the path where rebuilding would hurt most.
        let (ctx, run) = (egui::Context::default(), run());
        let mut map = BasinMap::default();
        let _ = repaint(&ctx, &mut map, &run);
        map.panels[Side::Left.index()].series.select(&run, middle_point(&run));
        map.scrubber.set_index(2);
        let _ = repaint(&ctx, &mut map, &run);
        map.playback.set_frames_per_second(60.0);
        map.playback.play(&mut map.scrubber);
        for _ in 0..10 {
            let _ = repaint(&ctx, &mut map, &run);
        }
        assert_eq!(panel_series(&map, Side::Left).walks, 1);
    }

    #[test]
    fn toggling_the_time_series_rebuilds_nothing() {
        // The same structural claim T-09.1 and T-09.3 make for the overlay and
        // the section: the chart is held by the panel, not the other way
        // round, so there is no path from the toggle to the frame under it or
        // to the series beside it.
        let (ctx, run) = (egui::Context::default(), run());
        let mut map = BasinMap::default();
        let first = repaint(&ctx, &mut map, &run);
        map.panels[Side::Left.index()].series.select(&run, middle_point(&run));
        map.layers.series = !map.layers.series;
        assert_eq!(repaint(&ctx, &mut map, &run), first);
        map.layers.series = !map.layers.series;
        assert_eq!(repaint(&ctx, &mut map, &run), first);
        assert_eq!(panel_series(&map, Side::Left).walks, 1);
    }

    #[test]
    fn choosing_another_frame_redraws_the_map_but_not_the_colour_bar() {
        let (ctx, run) = (egui::Context::default(), run());
        let mut map = BasinMap::default();
        let (first_map, first_bar) = repaint(&ctx, &mut map, &run);
        map.scrubber.set_index(2);
        let (second_map, second_bar) = repaint(&ctx, &mut map, &run);
        assert_ne!(
            second_map, first_map,
            "the map of another frame is another map"
        );
        // The scale is the run's, not the frame's, so the bar is the same bar.
        assert_eq!(second_bar, first_bar);
    }

    #[test]
    fn scrubbing_moves_both_panels_to_the_same_frame() {
        // The ticket's acceptance criterion, from the scrubbing half: the two
        // panels are not kept in step, they are drawn from one index, so what
        // is asserted is that each panel drew the frame the chooser names —
        // read off the values in the frame, not off the index.
        let ctx = egui::Context::default();
        let (left, right) = (
            run_of("control", 3, 1.0),
            // A different amplitude, so a panel drawing the wrong run's frame
            // would show a different number rather than a coincidentally equal
            // one.
            run_of("perturbed", 3, 10.0),
        );
        let mut map = BasinMap::default();
        for index in 0..3 {
            map.scrubber.set_index(index);
            repaint_comparison(&ctx, &mut map, &left, &right);
            #[allow(clippy::cast_precision_loss)]
            let frame = index as f64;
            assert_eq!(panel_h_m(&map, Side::Left), frame);
            assert_eq!(panel_h_m(&map, Side::Right), frame * 10.0);
        }
    }

    #[test]
    fn playing_carries_both_panels_through_the_run_together() {
        // The other half of the criterion: the clock writes the one index the
        // two panels read, so every repaint has them on the same frame, and
        // the pair stops at the last frame together.
        let ctx = egui::Context::default();
        let (left, right) = (run_of("control", 3, 1.0), run_of("perturbed", 3, 10.0));
        let mut map = BasinMap::default();
        repaint_comparison(&ctx, &mut map, &left, &right);
        map.playback.set_frames_per_second(60.0);
        map.playback.play(&mut map.scrubber);
        // Sixty frames a second against a repaint loop running rather slower:
        // ten repaints is far more than the three frames of these runs, so
        // what is asserted is where they stopped, not how fast.
        for _ in 0..10 {
            repaint_comparison(&ctx, &mut map, &left, &right);
            assert_eq!(panel_h_m(&map, Side::Left) * 10.0, panel_h_m(&map, Side::Right));
        }
        assert_eq!(map.scrubber.index(), 2, "playback stops at the last frame");
        assert!(!map.playback.is_playing());
        assert_eq!(panel_h_m(&map, Side::Left), 2.0);
        assert_eq!(panel_h_m(&map, Side::Right), 20.0);
    }

    #[test]
    fn a_shorter_run_stops_the_pair_at_the_frames_they_share() {
        // Past the shorter run's last frame there is nothing to put in the
        // second panel, so the chooser reaches only as far as both runs do.
        let ctx = egui::Context::default();
        let (left, right) = (run_of("long", 3, 1.0), run_of("short", 2, 10.0));
        let mut map = BasinMap::default();
        // The first repaint is what fits the chooser to the pair; asking for
        // frame 2 after it is asking for a frame the shorter run does not have.
        repaint_comparison(&ctx, &mut map, &left, &right);
        map.scrubber.set_index(2);
        repaint_comparison(&ctx, &mut map, &left, &right);
        assert_eq!(map.scrubber.index(), 1);
        assert_eq!(panel_h_m(&map, Side::Left), 1.0);
        assert_eq!(panel_h_m(&map, Side::Right), 10.0);
    }

    #[test]
    fn both_panels_are_drawn_on_the_one_colour_bar() {
        // One bar for the pair, sampled from the scale that covers both runs:
        // two bars would be two scales, and two scales are what would make the
        // panels incomparable (`crate::comparison`).
        let ctx = egui::Context::default();
        let (left, right) = (run_of("control", 3, 1.0), run_of("perturbed", 3, 10.0));
        let mut map = BasinMap::default();
        repaint_comparison(&ctx, &mut map, &left, &right);
        let comparison = Comparison::of(&left, &right).expect("these runs are comparable");
        let bar = map.bar.as_ref().expect("a colour bar was built");
        assert_eq!(bar.scale, comparison.scale());
        assert_eq!(bar.scale, right.anomaly_scale(), "the wider of the two");
        assert_ne!(bar.scale, left.anomaly_scale());
    }

    #[test]
    fn repainting_the_same_pair_rebuilds_nothing() {
        let ctx = egui::Context::default();
        let (left, right) = (run_of("control", 3, 1.0), run_of("perturbed", 3, 10.0));
        let mut map = BasinMap::default();
        repaint_comparison(&ctx, &mut map, &left, &right);
        let drawn = Side::BOTH.map(|side| panel_frame(&map, side).map.id());
        let bar = map
            .bar
            .as_ref()
            .expect("a colour bar was built")
            .texture
            .id();
        repaint_comparison(&ctx, &mut map, &left, &right);
        assert_eq!(
            Side::BOTH.map(|side| panel_frame(&map, side).map.id()),
            drawn
        );
        assert_eq!(
            map.bar
                .as_ref()
                .expect("a colour bar was built")
                .texture
                .id(),
            bar
        );
    }

    #[test]
    fn putting_a_louder_run_beside_a_run_redraws_it_on_the_shared_scale() {
        // A panel's map is only right for the scale it was drawn on, and in a
        // comparison the scale is both runs'. So the cache is keyed on the
        // scale as well as the frame: without that, the run drawn alone would
        // keep the colours of its own scale while the colour bar under it
        // stated another.
        let ctx = egui::Context::default();
        let (left, right) = (run_of("control", 3, 1.0), run_of("perturbed", 3, 10.0));
        let mut map = BasinMap::default();
        let _fitted = repaint(&ctx, &mut map, &left);
        map.scrubber.set_index(2);
        let (alone, _bar) = repaint(&ctx, &mut map, &left);
        repaint_comparison(&ctx, &mut map, &left, &right);
        assert_ne!(panel_frame(&map, Side::Left).map.id(), alone);
        // And it is still the same frame of the same run.
        assert_eq!(panel_h_m(&map, Side::Left), 2.0);
    }

    #[test]
    fn a_picked_cell_is_the_same_cell_in_both_panels_and_walks_each_run_once() {
        // The cell is the pair's, because the two runs are over one grid: one
        // click asks each panel for the series of the same place in its own
        // run, and each walks its own run once for it (T-09.4's cost property,
        // held per panel).
        let ctx = egui::Context::default();
        let (left, right) = (run_of("control", 3, 1.0), run_of("perturbed", 3, 10.0));
        let mut map = BasinMap::default();
        repaint_comparison(&ctx, &mut map, &left, &right);
        map.selected = Some(middle_point(&left));
        repaint_comparison(&ctx, &mut map, &left, &right);
        for side in Side::BOTH {
            let series = panel_series(&map, side);
            assert_eq!(series.walks, 1);
            assert_eq!(
                series.shown().expect("a cell was picked").point(),
                middle_point(&left)
            );
        }
        // The two series are of the same place in two runs, so they differ by
        // the amplitude the runs differ by.
        let h_m = |side| {
            panel_series(&map, side)
                .shown()
                .expect("a cell was picked")
                .samples()
                .iter()
                .map(SeriesSample::h_m)
                .collect::<Vec<f64>>()
        };
        assert_eq!(h_m(Side::Left), vec![0.0, 1.0, 2.0]);
        assert_eq!(h_m(Side::Right), vec![0.0, 10.0, 20.0]);
        // And repainting the pair walks neither run again.
        repaint_comparison(&ctx, &mut map, &left, &right);
        for side in Side::BOTH {
            assert_eq!(panel_series(&map, side).walks, 1);
        }
    }

    #[test]
    fn the_frame_line_reports_the_time_both_runs_put_the_frame_at() {
        // One time, because both runs wrote the same one: 86 400 s is a day.
        let line = super::frame_line(11, 731, [Some(1_036_800.0), Some(1_036_800.0)]);
        assert!(line.contains("Frame 12 of 731"), "{line}");
        assert!(line.contains("at 12.00 days"), "{line}");
    }

    #[test]
    fn the_frame_line_says_when_the_two_runs_date_the_frame_differently() {
        // A run continued from day 100 and a run started at zero share an
        // index and nothing else; the line reports both times rather than
        // inventing one from the cadence.
        let line = super::frame_line(0, 366, [Some(0.0), Some(8_640_000.0)]);
        assert!(line.contains("run A dates 0.00 days"), "{line}");
        assert!(line.contains("run B dates 100.00 days"), "{line}");
    }
}
