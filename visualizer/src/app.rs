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

use crate::loading::Loaded;
use crate::run::SECONDS_PER_DAY;
use crate::{
    BasinPoint, CrossSection, DivergingScale, Heatmap, LoadedRun, Loader, PendingRun, Playback,
    PointSeries, Scrubber, WindOverlay,
};

/// What the central panel is showing.
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

/// The visualizer's application state.
pub struct VisualizerApp {
    /// What the central panel is showing.
    shown: Shown,
    /// Dropped files seen so far, waiting for their pair.
    pending: PendingRun,
    /// The one channel every source of run bytes posts to.
    loader: Loader,
    /// The URL a run is served under, as typed or as passed in `?run=`.
    run_url: String,
    /// The basin map: which frame it shows, and what it last drew.
    basin_map: BasinMap,
}

impl Default for VisualizerApp {
    fn default() -> Self {
        Self {
            shown: Shown::Nothing,
            pending: PendingRun::default(),
            loader: Loader::default(),
            run_url: String::new(),
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

    /// Start fetching the run served under `base_url`, showing it as loading
    /// until it lands.
    pub fn fetch_run(&mut self, base_url: &str, ctx: &egui::Context) {
        self.run_url = base_url.to_owned();
        self.shown = Shown::Loading(base_url.to_owned());
        let ctx = ctx.clone();
        self.loader.fetch(base_url, move || ctx.request_repaint());
    }

    /// Load the run in `directory`. Native only: a browser has no directories
    /// to open (ADR-0006).
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load_directory(&mut self, directory: &std::path::Path) {
        let source = directory.display().to_string();
        self.shown = Shown::Loading(source.clone());
        self.loader.deliver(
            source,
            crate::loading::native::read_run_directory(directory),
        );
    }

    /// Ask the user for a run directory, on a thread of its own so the frame
    /// loop keeps running while the dialog is open.
    #[cfg(not(target_arch = "wasm32"))]
    fn pick_directory(&self, ctx: &egui::Context) {
        let sender = self.loader.sender();
        let ctx = ctx.clone();
        std::thread::spawn(move || {
            let Some(directory) = rfd::FileDialog::new()
                .set_title("Open a run directory")
                .pick_folder()
            else {
                return;
            };
            let _ = sender.send(Loaded {
                source: directory.display().to_string(),
                bytes: crate::loading::native::read_run_directory(&directory),
            });
            ctx.request_repaint();
        });
    }

    /// Take whatever finished loading since the last frame and show it.
    fn absorb_finished_loads(&mut self) {
        while let Some(Loaded { source, bytes }) = self.loader.poll() {
            // Whatever arrives, the map of the last run is no longer the map
            // of the run on screen.
            self.basin_map.forget();
            self.shown = match bytes.and_then(|bytes| {
                LoadedRun::from_bytes(source.clone(), bytes).map_err(|error| error.to_string())
            }) {
                Ok(run) => Shown::Run(Box::new(run)),
                Err(message) => Shown::Failed { source, message },
            };
        }
    }

    /// Take the files dropped on the window this frame, and load the run once
    /// both of them have arrived.
    fn absorb_dropped_files(&mut self, ctx: &egui::Context) {
        let dropped = ctx.input(|input| input.raw.dropped_files.clone());
        for file in &dropped {
            #[cfg(not(target_arch = "wasm32"))]
            if let Some(path) = file.path.as_ref().filter(|path| path.is_dir()) {
                self.load_directory(path);
                continue;
            }
            match dropped_file_contents(file) {
                Ok((name, bytes)) => {
                    if !self.pending.offer(&name, bytes) {
                        self.shown = Shown::Failed {
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
                    self.shown = Shown::Failed {
                        source: "dropped file".to_owned(),
                        message,
                    };
                }
            }
        }
        if let Some(bytes) = self.pending.take_run() {
            self.loader.deliver("dropped files", Ok(bytes));
        }
    }

    /// The bar of run-loading affordances.
    ///
    /// Returns whether the run-URL field has the keyboard. It is the one thing
    /// in the shell that a keystroke means something different to, so it is
    /// what decides whether the scrubber's keys are the scrubber's to take.
    fn draw_controls(&mut self, ui: &mut egui::Ui) -> bool {
        ui.horizontal(|ui| {
            #[cfg(not(target_arch = "wasm32"))]
            if ui.button("Open run directory…").clicked() {
                self.pick_directory(ui.ctx());
            }
            ui.label("Run URL:");
            let url = ui.add(
                egui::TextEdit::singleline(&mut self.run_url)
                    .hint_text("https://…/run-demo/")
                    .desired_width(260.0),
            );
            let submitted =
                url.lost_focus() && ui.input(|input| input.key_pressed(egui::Key::Enter));
            let fetch = ui.button("Fetch").clicked();
            if (submitted || fetch) && !self.run_url.trim().is_empty() {
                let (url, ctx) = (self.run_url.clone(), ui.ctx().clone());
                self.fetch_run(&url, &ctx);
            }
            url.has_focus()
        })
        .inner
    }
}

/// What to show when there is no run yet, or the last one failed.
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

        // Disjoint borrows: the panel reads the run while the map it draws
        // caches a texture of one of its frames.
        let Self {
            shown,
            pending,
            basin_map,
            ..
        } = self;
        egui::CentralPanel::default().show(ctx, |ui| match &*shown {
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
            Shown::Run(run) => draw_run(ui, run, basin_map, keyboard_free),
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

/// The metadata panel, and under it the basin map of the chosen frame.
fn draw_run(ui: &mut egui::Ui, run: &LoadedRun, basin_map: &mut BasinMap, keyboard_free: bool) {
    ui.label(RichText::new(run.source()).strong());
    ui.add_space(6.0);
    egui::Grid::new("run-metadata")
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
    ui.add_space(12.0);
    ui.separator();
    basin_map.draw(ui, run, keyboard_free);
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

/// The basin map of one chosen frame of the loaded run.
///
/// The map is a texture rather than a mesh: `h` is one value per cell, and the
/// cheapest honest way to show a cell grid is one pixel per cell, magnified
/// without interpolation. It is also the way that costs the same on both
/// targets, which is what ADR-0006 asks of anything drawn here.
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
/// reader actually asks for, which is the work the drag is for.
///
/// Playback (T-09.2) adds nothing to that: it writes the scrubber's index and
/// nothing else, so a played frame and a dragged one cost the same, and a
/// paused run — which is how a run loads — asks for no repaint at all.
///
/// The wind overlay rides on that same [`Attempt`]: it is built with the frame
/// and drawn from geometry, so neither showing it nor hiding it is a reason to
/// rebuild anything, and it adds nothing per frame a drag passes through. The
/// equatorial cross-section (T-09.3) rides on it the same way, and is drawn
/// from the same run-wide scale as the colour bar, so it too costs nothing per
/// repaint and nothing per toggle.
///
/// # What the time series costs, and why it is not on that path
///
/// The point time series (T-09.4) is the one view whose cost is shaped
/// differently, because it is the one view that is not of a frame: it is one
/// cell of *every* frame, so the indexed lookup that makes the others cheap
/// buys it nothing and a rebuild walks the whole run — 731 frames of the
/// scenario run, against the one frame everything above reads.
///
/// It is therefore held here, beside the frame cache rather than inside it,
/// and rebuilt only when the reader clicks a **different** cell. Scrubbing,
/// playing, toggling either chart and re-clicking the selected cell all leave
/// it alone; [`BasinMap::series_walks`] counts the rebuilds so the tests can
/// say so by name.
struct BasinMap {
    /// The frame the reader has chosen, and the ways they choose another.
    scrubber: Scrubber,
    /// The clock that chooses frames on the reader's behalf.
    playback: Playback,
    /// Whether the wind-stress overlay is drawn over the map.
    show_wind: bool,
    /// Whether the equatorial cross-section is drawn under the map.
    show_section: bool,
    /// Whether the point time series is drawn under the map.
    show_series: bool,
    /// The last attempt at drawing a frame, if there has been one.
    attempt: Option<Attempt>,
    /// The run's colour bar, built the first time a frame of it is drawn.
    bar: Option<ColorBar>,
    /// The series the reader picked, and what building it has cost.
    ///
    /// Beside the frame cache rather than in it: the series is of the whole
    /// run, so the frame on screen is not a reason to rebuild it.
    series: SeriesCache,
}

/// The point time series on screen, and how many times one has been built.
///
/// The two travel together everywhere, because the count is only meaningful
/// against the series it counts: this is the one path in the shell that walks a
/// whole run ([`crate::time_series`]), and the count is how a test says so.
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
    /// A cache with nothing picked, which has cost nothing.
    const fn empty() -> Self {
        Self {
            series: None,
            walks: 0,
        }
    }

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

impl Default for BasinMap {
    fn default() -> Self {
        Self {
            scrubber: Scrubber::default(),
            // Paused, as `crate::playback` says: a run that started moving on
            // load would be off its first frame before the header had been
            // read, and it is also what leaves an idle repaint with nothing to
            // rebuild.
            playback: Playback::new(),
            // On by default: the forcing is why the map looks the way it does,
            // and a reader who does not know the overlay exists cannot ask for
            // it.
            show_wind: true,
            // On by default: the tilt is what the run is about, and the
            // section is the view that states it as a number rather than as a
            // colour (T-09.3).
            show_section: true,
            // On by default, with nothing selected: the chart says how to pick
            // a point, and a reader who does not know the map is clickable
            // cannot find out by clicking it.
            show_series: true,
            attempt: None,
            bar: None,
            series: SeriesCache::empty(),
        }
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

/// The colour bar of a run, and the scale it was sampled from.
///
/// The scale is the run's rather than the frame's (`crate::heatmap`), so the
/// bar is the same in every frame and is uploaded once — not once per frame a
/// drag passes through.
struct ColorBar {
    /// The scale the bar was sampled from, so a run whose scale is not this
    /// one gets its own bar.
    scale: DivergingScale,
    /// The bar itself, sampled across the scale.
    texture: egui::TextureHandle,
}

impl BasinMap {
    /// Forget the run this was a map of. The overlay toggle is the reader's
    /// choice rather than the run's, so it survives.
    fn forget(&mut self) {
        self.scrubber = Scrubber::new();
        // The clock was playing the run that has just gone; the speed it was
        // playing at is the reader's choice, so that stays.
        self.playback.pause();
        self.attempt = None;
        self.bar = None;
        // The cell was a place in the basin of the run that has just gone; the
        // next run's cell of that index is a different place.
        self.series.forget();
    }

    /// Draw the frame chooser and the map of the chosen frame.
    fn draw(&mut self, ui: &mut egui::Ui, run: &LoadedRun, keyboard_free: bool) {
        let frame_count = run.header().output.frame_count;
        self.scrubber.fit_to(frame_count);
        if self.scrubber.last().is_none() {
            ui.label("This run holds no frames to draw.");
            return;
        }
        // The clock before the controls: the frame this repaint draws is the
        // one whatever time has passed since the last one has bought, and
        // every affordance below writes the same index it does.
        self.playback
            .advance(&mut self.scrubber, f64::from(ui.input(|i| i.stable_dt)));
        self.scrubber.draw(ui, keyboard_free);
        self.playback.draw(ui, &mut self.scrubber, keyboard_free);
        // The overlay is drawn over the map, never into it: nothing the map is
        // built from depends on this.
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_wind, "Wind stress τ");
            ui.checkbox(&mut self.show_section, "Equatorial cross-section");
            ui.checkbox(&mut self.show_series, "Point time series");
        });

        let index = self.scrubber.index();
        if self.attempt.as_ref().is_none_or(|last| last.index != index) {
            self.attempt = Some(Attempt {
                index,
                outcome: self.build(ui, run, index),
            });
        }
        // Before the frame is read back, because that read borrows the panel
        // for as long as the frame is drawn and this is the one thing left
        // that writes to it.
        self.ensure_color_bar(ui, run.anomaly_scale());
        // Disjoint borrows: the frame already drawn is read while the series —
        // which is of the whole run rather than of that frame — is picked and
        // built beside it.
        let Self {
            show_wind,
            show_section,
            show_series,
            attempt,
            bar,
            series,
            ..
        } = self;
        let outcome = &attempt.as_ref().expect("an attempt was just made").outcome;
        let drawn = match outcome {
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

        // Counted from one, as the metadata panel counts the run's frames.
        // The scrubber's own index starts at zero, and it shows no number.
        ui.label(format!(
            "Frame {} of {frame_count} — thermocline depth anomaly h at {:.2} days",
            index + 1,
            drawn.t_s / SECONDS_PER_DAY
        ));
        ui.horizontal(|ui| {
            ui.label("west");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label("east");
            });
        });
        let reserved_pt = reserved_below_map_pt(ui, *show_section, *show_series);
        let map_area = draw_texture_fitted(ui, &drawn.map, reserved_pt);
        let map = map_area.rect;
        if *show_wind {
            draw_wind_arrows(ui, map, drawn);
        }
        // The click that picks a place, and the only path in the shell that
        // walks a whole run. It is taken before the charts are drawn, so the
        // series on screen is the one the click just asked for, and it is
        // guarded on the cell rather than on the click: a reader clicking twice
        // in the same cell has asked for the series they already have.
        if let Some(point) = clicked_point(&map_area, map, run.header().grid) {
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
        if *show_section {
            draw_cross_section(ui, &drawn.section, map);
        }
        // The time series goes under the section, and only its horizontal
        // extent is the map's: its axis is time, so nothing on it lines up with
        // a column of the map.
        if *show_series {
            draw_time_series(ui, series.shown(), drawn.t_s, map);
        }
        draw_color_bar(ui, bar.as_ref().expect("a colour bar was just built"));
        if *show_wind {
            ui.label(format!(
                "Wind stress τ: the longest arrow is {:.3} N m^-2, the strongest in the run",
                drawn.wind.scale().max_magnitude_pa()
            ));
        }
    }

    /// Sample the run's colour bar, the first time a frame of it is drawn.
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

    /// Colour-map frame `index` and upload it, or say what stopped that.
    fn build(&self, ui: &egui::Ui, run: &LoadedRun, index: u64) -> Result<DrawnFrame, String> {
        let frame = run
            .frame(index)
            .ok_or_else(|| format!("this run holds no frame {index}"))?;
        let heatmap = Heatmap::of_frame(run.header().grid, &frame, run.anomaly_scale())
            .map_err(|error| error.to_string())?;
        let wind = WindOverlay::of_frame(run.header().grid, &frame, run.wind_stress_scale())
            .map_err(|error| error.to_string())?;
        let section = CrossSection::of_frame(run.header().grid, &frame, run.anomaly_scale())
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

    use super::{BasinMap, LoadedRun};
    use crate::{BasinPoint, RunBytes, SeriesSample};

    /// A basin small enough to build in a unit test, on the extent of
    /// `CONTEXT.md`, *Basin*.
    fn grid() -> GridSpec {
        GridSpec::new(4, 3, BasinExtent::new(120.0, -80.0, -25.0, 25.0))
            .expect("a 4x3 basin is a valid grid")
    }

    /// A run of three frames whose `h` is everywhere the frame's own index, in
    /// metres, so consecutive frames are drawn in different colours.
    fn run() -> LoadedRun {
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
            "basin-map",
            OutputTiming {
                frame_count: 3,
                interval_s: 86_400.0,
            },
        );
        let mut frames = Vec::new();
        for index in 0..header.output.frame_count {
            #[allow(clippy::cast_precision_loss)]
            let value = index as f64;
            let field = |variable| vec![0.0; grid.field_len(variable)];
            let frame = Frame::new(
                value * header.output.interval_s,
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
            "basin-map",
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
        let drawn = map
            .attempt
            .as_ref()
            .expect("a frame was drawn")
            .outcome
            .as_ref()
            .expect("it drew");
        let bar = map.bar.as_ref().expect("a colour bar was built");
        (drawn.map.id(), bar.texture.id())
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
        map.show_section = !map.show_section;
        assert_eq!(repaint(&ctx, &mut map, &run), first);
        map.show_section = !map.show_section;
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
            let drawn = map
                .attempt
                .as_ref()
                .expect("a frame was drawn")
                .outcome
                .as_ref()
                .expect("it drew");
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
        map.series.select(&run, middle_point(&run));
        assert_eq!(map.series.walks, 1);
        map.series.select(&run, middle_point(&run));
        assert_eq!(map.series.walks, 1, "the same cell is the series in hand");
        // And it is a series of this run: `h` is everywhere the frame's own
        // index in metres, so the samples say which run they came from.
        let series = map.series.shown().expect("a point was picked");
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
        map.series.select(&run, middle_point(&run));
        map.scrubber.set_index(2);
        let _ = repaint(&ctx, &mut map, &run);
        map.playback.set_frames_per_second(60.0);
        map.playback.play(&mut map.scrubber);
        for _ in 0..10 {
            let _ = repaint(&ctx, &mut map, &run);
        }
        assert_eq!(map.series.walks, 1);
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
        map.series.select(&run, middle_point(&run));
        map.show_series = !map.show_series;
        assert_eq!(repaint(&ctx, &mut map, &run), first);
        map.show_series = !map.show_series;
        assert_eq!(repaint(&ctx, &mut map, &run), first);
        assert_eq!(map.series.walks, 1);
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
}
