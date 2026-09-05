//! The shell itself: a window (or a canvas) that loads a run and says what it
//! is.
//!
//! T-08.1 draws no data. What it draws is the header — the grid, the scenario,
//! the frame count — because that is what tells a reader the run they think
//! they opened is the run they opened. Epic 09 fills the space below it.
//!
//! Everything with a value in it lives in [`crate::run`] and
//! [`crate::pending`]; this module is the part that needs a GPU, and so is
//! deliberately thin.

use egui::{Color32, RichText};

use crate::loading::Loaded;
use crate::{LoadedRun, Loader, PendingRun};

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
}

impl Default for VisualizerApp {
    fn default() -> Self {
        Self {
            shown: Shown::Nothing,
            pending: PendingRun::default(),
            loader: Loader::default(),
            run_url: String::new(),
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

    /// What to show when there is no run yet, or the last one failed.
    fn draw_instructions(&self, ui: &mut egui::Ui) {
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
        let still_needed = self.pending.still_needed();
        if still_needed.len() == 1 {
            ui.label(format!("Waiting for {}.", still_needed[0]));
        }
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

        egui::CentralPanel::default().show(ctx, |ui| match &self.shown {
            Shown::Nothing => self.draw_instructions(ui),
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
                self.draw_instructions(ui);
            }
            Shown::Run(run) => draw_run(ui, run),
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

/// The metadata panel: every row the run's header yields, labelled.
fn draw_run(ui: &mut egui::Ui, run: &LoadedRun) {
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
    ui.label("No frames drawn yet — rendering lands in Epic 09.");
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
