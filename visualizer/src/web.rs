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
                Box::new(|_cc| {
                    let mut app = VisualizerApp::new();
                    app.compute_default_run();
                    Ok(Box::new(app))
                }),
            )
            .await;
        if let Err(error) = result {
            log::error!("the visualizer failed to start: {error:?}");
        }
    });
}

/// The canvas named by [`CANVAS_ID`].
fn canvas() -> Option<web_sys::HtmlCanvasElement> {
    web_sys::window()?
        .document()?
        .get_element_by_id(CANVAS_ID)?
        .dyn_into::<web_sys::HtmlCanvasElement>()
        .ok()
}
