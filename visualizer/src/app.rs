//! The shell itself: a window (or a canvas) that loads a run, says what it is,
//! and draws one frame of it.
//!
//! It draws the header first — the grid, the scenario, the frame count —
//! because that is what tells a reader the run they think they opened is the
//! run they opened, and under it the basin map of one chosen frame.
//!
//! Everything with a value in it lives in [`crate::run`], [`crate::heatmap`],
//! [`crate::wind`], [`crate::cross_section`] and [`crate::pending`]; this
//! module is the part that needs a
//! GPU, and so is deliberately thin. What it adds on top of them is a texture
//! cache and a layout, and neither is where a wrong basin map would come from.

use egui::{Color32, RichText};

use crate::loading::Loaded;
use crate::run::SECONDS_PER_DAY;
use crate::{
    CrossSection, DivergingScale, Heatmap, LoadedRun, Loader, PendingRun, Playback, Scrubber,
    WindOverlay,
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

/// The arrow itself: near-black, so it reads over the pale middle of the
/// colour scale and over the casing at the scale's two dark ends.
const ARROW_COLOR: Color32 = Color32::from_rgb(16, 16, 16);

/// The casing under it: opaque white, so the near-black arrow stays separable
/// from the dark blue and dark red the scale ends on.
const ARROW_CASING_COLOR: Color32 = Color32::from_rgb(245, 245, 245);

/// Height of the cross-section chart, in points.
///
/// Tall enough that the ±`half_range` axis has room to show a tilt changing by
/// a few per cent, and short enough that the basin map above it still gets the
/// bulk of a window: the chart says how much, the map says where.
const CROSS_SECTION_HEIGHT_PT: f32 = 120.0;

/// Width of the cross-section line, in points.
const CROSS_SECTION_WIDTH_PT: f32 = 1.6;

/// The cross-section line: the deep end of the basin map's colour scale, so
/// the chart and the map are read as one picture of the same field.
const CROSS_SECTION_COLOR: Color32 = Color32::from_rgb(178, 24, 43);

/// The chart's frame and its zero line: grey, so the line drawn over them is
/// what the eye lands on.
const CROSS_SECTION_AXIS_COLOR: Color32 = Color32::from_gray(128);

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
struct BasinMap {
    /// The frame the reader has chosen, and the ways they choose another.
    scrubber: Scrubber,
    /// The clock that chooses frames on the reader's behalf.
    playback: Playback,
    /// Whether the wind-stress overlay is drawn over the map.
    show_wind: bool,
    /// Whether the equatorial cross-section is drawn under the map.
    show_section: bool,
    /// The last attempt at drawing a frame, if there has been one.
    attempt: Option<Attempt>,
    /// The run's colour bar, built the first time a frame of it is drawn.
    bar: Option<ColorBar>,
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
            attempt: None,
            bar: None,
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
        let outcome = &self
            .attempt
            .as_ref()
            .expect("an attempt was just made")
            .outcome;
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
        // The chart and the two rows of text around it come out of the height
        // the map would otherwise have taken, so that turning it on shrinks the
        // map rather than pushing the colour bar off a short window.
        let reserved_pt = if self.show_section {
            CROSS_SECTION_HEIGHT_PT + ui.text_style_height(&egui::TextStyle::Body) * 2.0
        } else {
            0.0
        };
        let map = draw_texture_fitted(ui, &drawn.map, reserved_pt);
        if self.show_wind {
            draw_wind_arrows(ui, map, drawn);
        }
        // Directly under the map and across exactly its width, so a longitude
        // on the chart sits under the column of the map it came from. The
        // colour bar goes below both: the scale is the same one.
        if self.show_section {
            draw_cross_section(ui, &drawn.section, map);
        }
        draw_color_bar(ui, self.bar.as_ref().expect("a colour bar was just built"));
        if self.show_wind {
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
    let axis = egui::Stroke::new(1.0_f32, CROSS_SECTION_AXIS_COLOR);
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

    #[allow(clippy::cast_possible_truncation)]
    let to_screen = |(x, y): (f64, f64)| {
        chart.min + egui::vec2(x as f32 * chart.width(), y as f32 * chart.height())
    };
    let mut run_of_points: Vec<egui::Pos2> = Vec::new();
    let flush = |points: &mut Vec<egui::Pos2>| {
        if points.len() > 1 {
            painter.add(egui::Shape::line(
                std::mem::take(points),
                egui::Stroke::new(CROSS_SECTION_WIDTH_PT, CROSS_SECTION_COLOR),
            ));
        } else {
            points.clear();
        }
    };
    for point in section.points() {
        match section.plot_position(point) {
            Some(position) => run_of_points.push(to_screen(position)),
            None => flush(&mut run_of_points),
        }
    }
    flush(&mut run_of_points);

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
/// say where it landed so a layer can be drawn over it.
///
/// The basin is far wider than it is tall, so fitting it to the width alone
/// would push the colour bar off the bottom of a short window.
fn draw_texture_fitted(
    ui: &mut egui::Ui,
    texture: &egui::TextureHandle,
    reserved_pt: f32,
) -> egui::Rect {
    let size = texture.size_vec2();
    let available = egui::vec2(
        ui.available_width(),
        (ui.available_height() - COLOR_BAR_HEIGHT * 3.0 - reserved_pt).max(COLOR_BAR_HEIGHT),
    );
    let scale = (available.x / size.x).min(available.y / size.y);
    ui.add(egui::Image::new(egui::load::SizedTexture::new(
        texture.id(),
        size * scale,
    )))
    .rect
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
            (ARROW_CASING_WIDTH_PT, ARROW_CASING_COLOR),
            (ARROW_WIDTH_PT, ARROW_COLOR),
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
    use crate::RunBytes;

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
