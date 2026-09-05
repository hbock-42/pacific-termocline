//! The shell itself: a window (or a canvas) that loads a run, says what it is,
//! and draws one frame of it.
//!
//! It draws the header first — the grid, the scenario, the frame count —
//! because that is what tells a reader the run they think they opened is the
//! run they opened, and under it the basin map of one chosen frame.
//!
//! Everything with a value in it lives in [`crate::run`], [`crate::heatmap`],
//! [`crate::wind`] and [`crate::pending`]; this module is the part that needs a
//! GPU, and so is deliberately thin. What it adds on top of them is a texture
//! cache and a layout, and neither is where a wrong basin map would come from.

use egui::{Color32, RichText};

use crate::loading::Loaded;
use crate::run::SECONDS_PER_DAY;
use crate::{DivergingScale, Heatmap, LoadedRun, Loader, PendingRun, Scrubber, WindOverlay};

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
/// The wind overlay rides on that same [`Attempt`]: it is built with the frame
/// and drawn from geometry, so neither showing it nor hiding it is a reason to
/// rebuild anything, and it adds nothing per frame a drag passes through.
struct BasinMap {
    /// The frame the reader has chosen, and the ways they choose another.
    scrubber: Scrubber,
    /// Whether the wind-stress overlay is drawn over the map.
    show_wind: bool,
    /// The last attempt at drawing a frame, if there has been one.
    attempt: Option<Attempt>,
    /// The run's colour bar, built the first time a frame of it is drawn.
    bar: Option<ColorBar>,
}

impl Default for BasinMap {
    fn default() -> Self {
        Self {
            scrubber: Scrubber::default(),
            // On by default: the forcing is why the map looks the way it does,
            // and a reader who does not know the overlay exists cannot ask for
            // it.
            show_wind: true,
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
        self.scrubber.draw(ui, keyboard_free);
        // The overlay is drawn over the map, never into it: nothing the map is
        // built from depends on this.
        ui.checkbox(&mut self.show_wind, "Wind stress τ");

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
        let map = draw_texture_fitted(ui, &drawn.map);
        if self.show_wind {
            draw_wind_arrows(ui, map, drawn);
        }
        // Read off before the colour bar, which needs the map by mutable
        // borrow to build the run's bar the first time it is asked for.
        let strongest_pa = drawn.wind.scale().max_magnitude_pa();
        draw_color_bar(ui, self.color_bar(ui, run.anomaly_scale()));
        if self.show_wind {
            ui.label(format!(
                "Wind stress τ: the longest arrow is {strongest_pa:.3} N m^-2, the strongest in the run"
            ));
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
        let wind = WindOverlay::of_frame(run.header().grid, &frame, run.wind_stress_scale())
            .map_err(|error| error.to_string())?;
        let image = egui::ColorImage::from_rgb([heatmap.width(), heatmap.height()], heatmap.rgb());
        Ok(DrawnFrame {
            t_s: frame.t_s(),
            wind,
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

/// Draw `texture` as large as the panel allows without changing its shape, and
/// say where it landed so a layer can be drawn over it.
///
/// The basin is far wider than it is tall, so fitting it to the width alone
/// would push the colour bar off the bottom of a short window.
fn draw_texture_fitted(ui: &mut egui::Ui, texture: &egui::TextureHandle) -> egui::Rect {
    let size = texture.size_vec2();
    let available = egui::vec2(
        ui.available_width(),
        (ui.available_height() - COLOR_BAR_HEIGHT * 3.0).max(COLOR_BAR_HEIGHT),
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
