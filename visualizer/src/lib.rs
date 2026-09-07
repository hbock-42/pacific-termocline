//! The visualizer: it draws a run — one the engine wrote to a directory, or
//! one it computed in the tab.
//!
//! # One shell, two platforms
//!
//! Per [ADR-0006] this is a browser app that also runs natively, built on
//! `eframe` + `wgpu` and compiled to `wasm32-unknown-unknown` as well as to
//! the host. The split is not between two apps but between where a run comes
//! from. Natively it is read: a directory, a pair of dropped files, or an HTTP
//! fetch, all of them [`RunBytes`] in and [`LoadedRun`] out. In a browser it
//! is *computed*: per [ADR-0012] the visualizer links the engine, holds a
//! scenario and steps it, because the file format is not served to the web at
//! all — 941 MB of control run is not a download. [`ComputedRun`] is that
//! loop.
//!
//! Both origins end in the same [`LoadedRun`], which is why the heatmap, the
//! scrubber, playback, the wind overlay, the cross-section, the point time
//! series and the comparison are unchanged by any of it: they consume a run,
//! not a file. And as before, the parts with a value in them are testable
//! without a GPU.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md
//! [ADR-0012]: ../../docs/planning/adr/0012-the-browser-runs-the-engine.md

mod app;
mod chart;
mod comparison;
mod compute;
mod cross_section;
mod heatmap;
#[cfg(not(target_arch = "wasm32"))]
mod loading;
#[cfg(not(target_arch = "wasm32"))]
mod pending;
mod playback;
mod run;
mod scrubber;
mod time_series;
mod wind;

pub use app::VisualizerApp;
pub use comparison::{Comparison, Difference, Mismatch, Side};
pub use compute::{
    BrowserScenario, BudgetExceeded, ComputeError, ComputedRun, FrameBudget, InMegabytes,
    STEP_BUDGET,
};
pub use cross_section::{CrossSection, CrossSectionPoint};
pub use heatmap::{DivergingScale, Heatmap};
/// Reading a *written* run is native-only since ADR-0012: the browser computes
/// its runs, and the file format is not served to it at all.
#[cfg(not(target_arch = "wasm32"))]
pub use loading::{native::read_run_directory, Loader};
#[cfg(not(target_arch = "wasm32"))]
pub use pending::PendingRun;
pub use playback::{Playback, MAX_STALL_S, PLAYBACK_SPEEDS_FPS};
pub use run::{FrameAppendError, LoadedRun, MetadataRow, RunBytes};
pub use scrubber::Scrubber;
pub use time_series::{BasinPoint, PointSeries, SeriesSample, SstScale};
pub use wind::{
    StressScale, WindArrow, WindOverlay, ARROW_SPACING_CELLS, MAX_ARROW_LENGTH_CELLS,
    MIN_ARROW_LENGTH_CELLS,
};

/// Re-exported so the visualizer and engine agree on one format version.
pub use termocline_format::FORMAT_VERSION;

/// The window title and, on the web, the name in the tab.
pub const APP_NAME: &str = "Termocline visualizer";

#[cfg(target_arch = "wasm32")]
mod web;
