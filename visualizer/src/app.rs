//! The shell itself: a window (or a canvas) that loads a run, says what it is,
//! and draws one frame of it.
//!
//! It draws the header first — the grid, the scenario, the frame count —
//! because that is what tells a reader the run they think they opened is the
//! run they opened, and under it the basin map of one chosen frame.
//!
//! Everything with a value in it lives in [`crate::run`], [`crate::heatmap`]
//! and [`crate::pending`]; this module is the part that needs a GPU, and so is
//! deliberately thin. What it adds on top of them is a texture cache and a
//! layout, and neither is where a wrong basin map would come from.

use egui::{Color32, RichText};

use crate::loading::Loaded;
use crate::run::SECONDS_PER_DAY;
use crate::{DivergingScale, Heatmap, LoadedRun, Loader, PendingRun, Scrubber};

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
    fn draw_controls(&mut self, ui: &mut egui::Ui) {
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
        });
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

        egui::TopBottomPanel::top("controls").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.heading(crate::APP_NAME);
            self.draw_controls(ui);
            ui.add_space(4.0);
        });

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
            Shown::Run(run) => draw_run(ui, run, basin_map),
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
fn draw_run(ui: &mut egui::Ui, run: &LoadedRun, basin_map: &mut BasinMap) {
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
    basin_map.draw(ui, run);
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
#[derive(Default)]
struct BasinMap {
    /// The frame the reader has chosen, and the ways they choose another.
    scrubber: Scrubber,
    /// The last attempt at drawing a frame, if there has been one.
    attempt: Option<Attempt>,
    /// The run's colour bar, built the first time a frame of it is drawn.
    bar: Option<ColorBar>,
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

/// A frame already colour-mapped and uploaded.
struct DrawnFrame {
    /// Its model time, in seconds since the start of the run.
    t_s: f64,
    /// The map itself, one texel per cell.
    map: egui::TextureHandle,
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

/// Frames a page key moves the chooser by.
///
/// Ten, against the arrow keys' one: the scenario writes a frame a day
/// (`steady-trades.toml`), so a page is a week and a half of model time — far
/// enough to see a change, near enough to still be reading the same event.
const FRAMES_PER_PAGE: i64 = 10;

impl BasinMap {
    /// Forget the run this was a map of.
    fn forget(&mut self) {
        self.scrubber = Scrubber::new();
        self.attempt = None;
        self.bar = None;
    }

    /// Draw the frame chooser and the map of the chosen frame.
    fn draw(&mut self, ui: &mut egui::Ui, run: &LoadedRun) {
        self.scrubber.fit_to(run.header().output.frame_count);
        let Some(last) = self.scrubber.last() else {
            ui.label("This run holds no frames to draw.");
            return;
        };
        self.draw_scrubber(ui, last);

        let index = self.scrubber.index();
        if self.attempt.as_ref().is_none_or(|last| last.index != index) {
            self.attempt = Some(Attempt {
                index,
                outcome: self.build(ui, run, index),
            });
        }
        let outcome = &self
            .attempt
            .as_ref()
            .expect("an attempt was just made")
            .outcome;
        let drawn = match outcome {
            Ok(drawn) => drawn,
            Err(message) => {
                ui.label(
                    RichText::new(format!("Frame {index} could not be drawn"))
                        .color(Color32::LIGHT_RED)
                        .strong(),
                );
                ui.label(message);
                return;
            }
        };

        ui.label(format!(
            "Frame {index} of {last} — thermocline depth anomaly h at {:.2} days",
            drawn.t_s / SECONDS_PER_DAY
        ));
        ui.horizontal(|ui| {
            ui.label("west");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label("east");
            });
        });
        draw_texture_fitted(ui, &drawn.map);
        let scale = run.anomaly_scale();
        let bar = self.color_bar(ui, scale);
        draw_color_bar(ui, bar, scale);
    }

    /// The scrubber itself: the slider, the steps either side of it, and the
    /// keys that do the same thing without the mouse.
    ///
    /// Every control changes one number. What makes the frame appear is the
    /// next repaint reading that number, so a drag, an arrow key and a jump to
    /// the end of the run all cost exactly the same.
    fn draw_scrubber(&mut self, ui: &mut egui::Ui, last: u64) {
        ui.horizontal(|ui| {
            let index = self.scrubber.index();
            if step_button(ui, "⏮", "First frame (Home)", index > 0) {
                self.scrubber.to_first();
            }
            if step_button(ui, "◀", "Back one frame (left arrow)", index > 0) {
                self.scrubber.step(-1);
            }
            if step_button(ui, "▶", "On one frame (right arrow)", index < last) {
                self.scrubber.step(1);
            }
            if step_button(ui, "⏭", "Last frame (End)", index < last) {
                self.scrubber.to_last();
            }
            // The slider takes the whole rest of the row: it is dragged, and a
            // frame of a long run is worth more pixels than a short slider
            // gives it — at 731 frames a narrow one puts several frames under
            // every pixel and none of them within reach.
            let mut chosen = self.scrubber.index();
            let slider = ui.add_sized(
                egui::vec2(ui.available_width(), ui.spacing().interact_size.y),
                egui::Slider::new(&mut chosen, 0..=last)
                    .integer()
                    .show_value(false),
            );
            if slider.changed() {
                self.scrubber.set_index(chosen);
            }
        });
        self.take_scrubber_keys(ui.ctx());
    }

    /// Move the chooser by whatever the keyboard asked for this frame.
    ///
    /// Nothing is taken while a text field has the keyboard: the run URL is
    /// typed into one, and an arrow key inside it belongs to the caret.
    fn take_scrubber_keys(&mut self, ctx: &egui::Context) {
        if ctx.wants_keyboard_input() {
            return;
        }
        let pressed = |key| ctx.input(|input| input.key_pressed(key));
        if pressed(egui::Key::Home) {
            self.scrubber.to_first();
        }
        if pressed(egui::Key::End) {
            self.scrubber.to_last();
        }
        for (key, frames) in [
            (egui::Key::ArrowLeft, -1),
            (egui::Key::ArrowRight, 1),
            (egui::Key::PageUp, -FRAMES_PER_PAGE),
            (egui::Key::PageDown, FRAMES_PER_PAGE),
        ] {
            if pressed(key) {
                self.scrubber.step(frames);
            }
        }
    }

    /// The run's colour bar, sampled the first time a frame of it is drawn.
    fn color_bar(&mut self, ui: &egui::Ui, scale: DivergingScale) -> &ColorBar {
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
        self.bar.as_ref().expect("a bar was just built")
    }

    /// Colour-map frame `index` and upload it, or say what stopped that.
    fn build(&self, ui: &egui::Ui, run: &LoadedRun, index: u64) -> Result<DrawnFrame, String> {
        let frame = run
            .frame(index)
            .ok_or_else(|| format!("this run holds no frame {index}"))?;
        let heatmap = Heatmap::of_frame(run.header().grid, &frame, run.anomaly_scale())
            .map_err(|error| error.to_string())?;
        let image = egui::ColorImage::from_rgb([heatmap.width(), heatmap.height()], heatmap.rgb());
        Ok(DrawnFrame {
            t_s: frame.t_s(),
            // Nearest, not linear: a texel is a cell of the model, and
            // smoothing between them would draw an anomaly the run never
            // produced.
            map: ui
                .ctx()
                .load_texture("basin-map", image, egui::TextureOptions::NEAREST),
        })
    }
}

/// One of the scrubber's step buttons: what it says, what it does when hovered
/// over, and whether there is anywhere for it to go.
fn step_button(ui: &mut egui::Ui, label: &str, hint: &str, enabled: bool) -> bool {
    ui.add_enabled(enabled, egui::Button::new(label))
        .on_hover_text(hint)
        .clicked()
}

/// Draw the colour bar and the anomalies its ends stand for.
fn draw_color_bar(ui: &mut egui::Ui, bar: &ColorBar, scale: DivergingScale) {
    let width = ui.available_width();
    ui.add(egui::Image::new(egui::load::SizedTexture::new(
        bar.texture.id(),
        egui::vec2(width, COLOR_BAR_HEIGHT),
    )));
    let half_range_m = scale.half_range_m();
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

/// Draw `texture` as large as the panel allows without changing its shape.
///
/// The basin is far wider than it is tall, so fitting it to the width alone
/// would push the colour bar off the bottom of a short window.
fn draw_texture_fitted(ui: &mut egui::Ui, texture: &egui::TextureHandle) {
    let size = texture.size_vec2();
    let available = egui::vec2(
        ui.available_width(),
        (ui.available_height() - COLOR_BAR_HEIGHT * 3.0).max(COLOR_BAR_HEIGHT),
    );
    let scale = (available.x / size.x).min(available.y / size.y);
    ui.add(egui::Image::new(egui::load::SizedTexture::new(
        texture.id(),
        size * scale,
    )));
}
