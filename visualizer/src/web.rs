//! The browser entry point.
//!
//! `wasm-bindgen` calls [`start`] as soon as the module is instantiated, so
//! the page needs no JavaScript of its own beyond a canvas to draw on. The
//! run to open is taken from the page's `?run=` query parameter — the HTTP
//! fetch of ADR-0006, which is also what makes the browser build testable
//! without a hand on a mouse. A second run named in `?compare=` opens beside
//! it (T-09.5), so a comparison is a link a reader can be sent.

use wasm_bindgen::prelude::*;
use wasm_bindgen::JsCast;

use crate::VisualizerApp;

/// Id of the canvas in `index.html` that the app draws on.
const CANVAS_ID: &str = "termocline_canvas";

/// Query parameter naming a run served over HTTP, as a base URL holding the
/// run's two files.
const RUN_PARAM: &str = "run";

/// Query parameter naming a second run to show beside the first.
const COMPARE_PARAM: &str = "compare";

/// Start the app on the page's canvas. Called by `wasm-bindgen` on load.
#[wasm_bindgen(start)]
pub fn start() {
    // Without this a Rust panic reaches the console as `unreachable executed`,
    // which says nothing about what failed.
    console_error_panic_hook::set_once();

    wasm_bindgen_futures::spawn_local(async {
        let canvas = canvas().expect("index.html defines the canvas the app draws on");
        let (run_url, compare_url) = (url_from_query(RUN_PARAM), url_from_query(COMPARE_PARAM));
        let result = eframe::WebRunner::new()
            .start(
                canvas,
                eframe::WebOptions::default(),
                Box::new(move |cc| {
                    let mut app = VisualizerApp::new();
                    if let Some(url) = run_url {
                        app.fetch_run(&url, &cc.egui_ctx);
                    }
                    if let Some(url) = compare_url {
                        app.fetch_run_to_compare(&url, &cc.egui_ctx);
                    }
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

/// The `name` parameter of the page's URL, if it has one.
fn url_from_query(name: &str) -> Option<String> {
    let search = web_sys::window()?.location().search().ok()?;
    let params = web_sys::UrlSearchParams::new_with_str(&search).ok()?;
    params.get(name).filter(|url| !url.is_empty())
}
