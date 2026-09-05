//! The visualizer: it reads the engine's output files through
//! [`termocline_format`] and never links against the simulation code, so
//! either side can be reimplemented without touching the other (ADR-0001).
//!
//! # One shell, two platforms
//!
//! Per [ADR-0006] this is a browser app that also runs natively, built on
//! `eframe` + `wgpu` and compiled to `wasm32-unknown-unknown` as well as to
//! the host. The split is not between two apps but between how a run reaches
//! one: a browser has no filesystem, so on the web a run arrives as dropped
//! files or an HTTP fetch, and only the native build opens a directory. Past
//! that point — [`RunBytes`] in, [`LoadedRun`] out — the two targets run the
//! same code, and the parts with a value in them are testable without a GPU.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

mod app;
mod loading;
mod pending;
mod run;

pub use app::VisualizerApp;
#[cfg(not(target_arch = "wasm32"))]
pub use loading::native::read_run_directory;
pub use loading::Loader;
pub use pending::PendingRun;
pub use run::{LoadedRun, MetadataRow, RunBytes};

/// Re-exported so the visualizer and engine agree on one format version.
pub use termocline_format::FORMAT_VERSION;

/// The window title and, on the web, the name in the tab.
pub const APP_NAME: &str = "Termocline visualizer";

#[cfg(target_arch = "wasm32")]
mod web;

#[cfg(test)]
mod tests {
    #[test]
    fn workspace_links_the_format_crate() {
        assert_eq!(crate::FORMAT_VERSION, termocline_format::FORMAT_VERSION);
    }
}
