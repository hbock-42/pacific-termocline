//! Native entry point for the visualizer.
//!
//! A run may be named on the command line — `termocline-viz /tmp/run-demo` —
//! picked from the window, or dropped onto it. In a browser there is no
//! command line and no filesystem, and [`visualizer::VisualizerApp`] is
//! started from `src/web.rs` instead (ADR-0006).

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result {
    env_logger::init();

    let directory = std::env::args_os().nth(1).map(std::path::PathBuf::from);
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([760.0, 520.0]),
        ..Default::default()
    };
    eframe::run_native(
        visualizer::APP_NAME,
        options,
        Box::new(move |_cc| {
            let mut app = visualizer::VisualizerApp::new();
            if let Some(directory) = directory {
                app.load_directory(&directory);
            }
            Ok(Box::new(app))
        }),
    )
}

/// The web build is started by `wasm-bindgen` from the library, not by a
/// `main`; this exists so the binary target still compiles for `wasm32`.
#[cfg(target_arch = "wasm32")]
fn main() {}
