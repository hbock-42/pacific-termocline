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
use crate::{DivergingScale, Heatmap, LoadedRun, Loader, PendingRun};

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
#[derive(Default)]
struct BasinMap {
    /// The frame the reader has chosen, by index into the run.
    index: u64,
    /// What was last uploaded to the GPU, if anything.
    drawn: Option<DrawnFrame>,
}

/// A frame already colour-mapped and uploaded.
///
/// Kept because reaching a frame decodes every frame before it
/// ([`LoadedRun::frame`]), and a panel repaints many times per second while
/// the chosen frame does not change.
struct DrawnFrame {
    /// Which frame this is a map of.
    index: u64,
    /// Its model time, in seconds since the start of the run.
    t_s: f64,
    /// The scale its colours came from, for the bar beside it.
    scale: DivergingScale,
    /// The map itself, one texel per cell.
    map: egui::TextureHandle,
    /// The colour bar, sampled across the same scale.
    bar: egui::TextureHandle,
}

impl BasinMap {
    /// Forget the run this was a map of.
    fn forget(&mut self) {
        self.index = 0;
        self.drawn = None;
    }

    /// Draw the frame chooser and the map of the chosen frame.
    fn draw(&mut self, ui: &mut egui::Ui, run: &LoadedRun) {
        let frame_count = run.header().output.frame_count;
        let Some(last) = frame_count.checked_sub(1) else {
            ui.label("This run holds no frames to draw.");
            return;
        };
        self.index = self.index.min(last);
        if last > 0 {
            ui.add(egui::Slider::new(&mut self.index, 0..=last).text("Frame"));
        }

        if self
            .drawn
            .as_ref()
            .is_none_or(|drawn| drawn.index != self.index)
        {
            self.drawn = self.build(ui, run);
        }
        let Some(drawn) = &self.drawn else {
            ui.label(
                RichText::new("This frame does not fit the grid its header describes.")
                    .color(Color32::LIGHT_RED),
            );
            return;
        };

        ui.label(format!(
            "Thermocline depth anomaly h at {:.2} days",
            drawn.t_s / SECONDS_PER_DAY
        ));
        ui.horizontal(|ui| {
            ui.label("west");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label("east");
            });
        });
        draw_texture_fitted(ui, &drawn.map);
        draw_color_bar(ui, drawn);
    }

    /// Colour-map the chosen frame and upload it, or `None` if the frame does
    /// not fit the grid.
    fn build(&self, ui: &egui::Ui, run: &LoadedRun) -> Option<DrawnFrame> {
        let frame = run.frame(self.index)?;
        let heatmap = Heatmap::of_frame(run.header().grid, &frame).ok()?;
        let image = egui::ColorImage::from_rgb([heatmap.width(), heatmap.height()], heatmap.rgb());
        Some(DrawnFrame {
            index: self.index,
            t_s: frame.t_s(),
            scale: *heatmap.scale(),
            // Nearest, not linear: a texel is a cell of the model, and
            // smoothing between them would draw an anomaly the run never
            // produced.
            map: ui
                .ctx()
                .load_texture("basin-map", image, egui::TextureOptions::NEAREST),
            bar: ui.ctx().load_texture(
                "basin-map-scale",
                color_bar_image(*heatmap.scale()),
                egui::TextureOptions::LINEAR,
            ),
        })
    }
}

/// The colour bar of `scale`, one row of samples from its shallow end to its
/// deep one.
fn color_bar_image(scale: DivergingScale) -> egui::ColorImage {
    let mut rgb = Vec::with_capacity(COLOR_BAR_SAMPLES * 3);
    for sample in 0..COLOR_BAR_SAMPLES {
        #[allow(clippy::cast_precision_loss)]
        let fraction = sample as f64 / (COLOR_BAR_SAMPLES - 1) as f64;
        rgb.extend_from_slice(&scale.color(scale.half_range_m() * (2.0 * fraction - 1.0)));
    }
    egui::ColorImage::from_rgb([COLOR_BAR_SAMPLES, 1], &rgb)
}

/// Draw the colour bar and the anomalies its ends stand for.
fn draw_color_bar(ui: &mut egui::Ui, drawn: &DrawnFrame) {
    let width = ui.available_width();
    ui.add(egui::Image::new(egui::load::SizedTexture::new(
        drawn.bar.id(),
        egui::vec2(width, COLOR_BAR_HEIGHT),
    )));
    let half_range_m = drawn.scale.half_range_m();
    ui.horizontal(|ui| {
        ui.label(format!("{:+.1} m (shallower)", -half_range_m));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(format!("{half_range_m:+.1} m (deeper)"));
            ui.with_layout(egui::Layout::top_down(egui::Align::Center), |ui| {
                ui.label("0 m");
            });
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
