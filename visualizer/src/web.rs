//! The browser entry point.
//!
//! `wasm-bindgen` calls [`start`] as soon as the module is instantiated, so
//! the page needs no JavaScript of its own beyond a canvas to draw on.
//!
//! There is no run to open. Per [ADR-0012] the file format is not served to
//! the browser at all — the control run is 941 MB against a 100 MB file cap —
//! so the app links the engine and computes a run instead, and the `?run=`
//! parameter this module used to read went with the fetch it named. What
//! first paint costs is a compile and a few steps rather than a transfer, and
//! the first scenario is started here so a visitor sees the ocean answering
//! the alizés without having pressed anything.
//!
//! [ADR-0012]: ../../docs/planning/adr/0012-the-browser-runs-the-engine.md

use wasm_bindgen::prelude::*;
use wasm_bindgen::JsCast;

use crate::VisualizerApp;

/// Id of the canvas in `index.html` that the app draws on.
const CANVAS_ID: &str = "termocline_canvas";

/// The app, plus the one measurement ADR-0012 is answerable for.
///
/// The trade the ADR made is a download for a compile and a few steps, and
/// what says whether that was a good trade is the wall-clock time between the
/// page being asked for and something being on screen. It is a property of the
/// deployed site — of the bundle's size over a real link, of the browser's
/// wasm compile, of the first slice of stepping — so it is measured there
/// rather than estimated here, and reported to the console on every load.
struct TimedFirstFrame {
    /// The shell being timed.
    app: VisualizerApp,
    /// The graphics backend `wgpu` actually got — `BrowserWebGpu`, or `Gl`
    /// where the browser has no WebGPU and the WebGL2 fallback took over
    /// (ADR-0006). Reported alongside the timing because which one is in use
    /// is otherwise invisible, and the two are not the same machine.
    backend: &'static str,
    /// Whether the first frame has been reported; the measurement is of the
    /// first one, and every frame after it costs a bool.
    reported: bool,
}

impl eframe::App for TimedFirstFrame {
    fn update(&mut self, ctx: &egui::Context, frame: &mut eframe::Frame) {
        self.app.update(ctx, frame);
        if !self.reported {
            self.reported = true;
            report_first_frame(self.backend);
        }
    }
}

/// Log how long the first frame took to arrive, in milliseconds since the
/// navigation started.
///
/// `performance.now()`'s origin is the navigation, not the module's
/// instantiation, so what this reports includes fetching the bundle and
/// compiling it — the number a visitor experiences. It is taken after
/// [`eframe::App::update`] has returned, so the frame is composed and
/// submitted; presentation is a compositor's frame later, which is below the
/// resolution of the thing being reported.
fn report_first_frame(backend: &str) {
    let Some(elapsed_ms) = web_sys::window()
        .and_then(|window| window.performance())
        .map(|performance| performance.now())
    else {
        return;
    };
    web_sys::console::log_1(
        &format!("termocline: first frame at {elapsed_ms:.0} ms since page load, on {backend}")
            .into(),
    );
}

/// Start the app on the page's canvas. Called by `wasm-bindgen` on load.
#[wasm_bindgen(start)]
pub fn start() {
    // Without this a Rust panic reaches the console as `unreachable executed`,
    // which says nothing about what failed.
    console_error_panic_hook::set_once();

    wasm_bindgen_futures::spawn_local(async {
        let canvas = canvas().expect("index.html defines the canvas the app draws on");
        let result = eframe::WebRunner::new()
            .start(
                canvas,
                eframe::WebOptions::default(),
                Box::new(|cc| {
                    let mut app = VisualizerApp::new();
                    app.compute_default_run();
                    Ok(Box::new(TimedFirstFrame {
                        app,
                        backend: backend_of(cc),
                        reported: false,
                    }))
                }),
            )
            .await;
        if let Err(error) = result {
            log::error!("the visualizer failed to start: {error:?}");
        }
    });
}

/// The name of the backend `wgpu` chose for this browser.
///
/// `"unknown"` only if `eframe` was built without its `wgpu` renderer, which
/// this crate's manifest does not allow.
fn backend_of(cc: &eframe::CreationContext<'_>) -> &'static str {
    cc.wgpu_render_state
        .as_ref()
        .map_or("unknown", |state| state.adapter.get_info().backend.to_str())
}

/// The canvas named by [`CANVAS_ID`].
fn canvas() -> Option<web_sys::HtmlCanvasElement> {
    web_sys::window()?
        .document()?
        .get_element_by_id(CANVAS_ID)?
        .dyn_into::<web_sys::HtmlCanvasElement>()
        .ok()
}
