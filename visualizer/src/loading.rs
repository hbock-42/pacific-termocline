//! Getting a run's bytes to the shell.
//!
//! Loading is asynchronous on both targets and for the same reason: an HTTP
//! fetch has no synchronous form in a browser, and a native folder picker
//! blocks until the user has chosen. So every source posts its result down one
//! channel, the UI drains it each frame, and nothing in [`crate::app`] has to
//! know which source it came from.

use std::sync::mpsc::{Receiver, Sender};

use termocline_format::{FRAME_FILE_NAME, HEADER_FILE_NAME};

use crate::{RunBytes, Side};

/// A run's bytes, or why they did not arrive, tagged with where they came from
/// and which panel asked for them.
pub struct Loaded {
    /// The panel the run was loaded for. A load that lands while the reader
    /// has moved on still belongs to the panel that started it, so nothing
    /// arriving late can overwrite the other one.
    pub side: Side,
    /// Where the run came from, as it is shown to the reader.
    pub source: String,
    /// The bytes, or a message naming what went wrong.
    pub bytes: Result<RunBytes, String>,
}

/// The channel every source posts to, and the UI drains.
pub struct Loader {
    /// Handed to each in-flight load.
    sender: Sender<Loaded>,
    /// Drained once per frame.
    receiver: Receiver<Loaded>,
}

impl Default for Loader {
    fn default() -> Self {
        let (sender, receiver) = std::sync::mpsc::channel();
        Self { sender, receiver }
    }
}

impl Loader {
    /// The result of a load that has finished since the last call, if any.
    pub fn poll(&self) -> Option<Loaded> {
        self.receiver.try_recv().ok()
    }

    /// Deliver a run that was assembled without an asynchronous source — a
    /// drop, which egui has already read for us.
    pub fn deliver(&self, side: Side, source: impl Into<String>, bytes: Result<RunBytes, String>) {
        // The receiver lives as long as the loader, so a send only fails if
        // the app is being torn down; there is nothing left to report to.
        let _ = self.sender.send(Loaded {
            side,
            source: source.into(),
            bytes,
        });
    }

    /// Fetch a run served under `base_url`: the two files of the format, by
    /// name, beneath it.
    ///
    /// The same code on both targets — `ehttp` is `fetch` in a browser and a
    /// request on a background thread natively. `repaint` is called when the
    /// result lands, because on the web nothing else wakes the frame loop.
    pub fn fetch(&self, side: Side, base_url: &str, repaint: impl Fn() + Send + 'static) {
        let base = base_url.trim();
        let source = base.to_owned();
        let prefix = if base.ends_with('/') {
            base.to_owned()
        } else {
            format!("{base}/")
        };
        let sender = self.sender.clone();
        let frames_url = format!("{prefix}{FRAME_FILE_NAME}");
        ehttp::fetch(
            ehttp::Request::get(format!("{prefix}{HEADER_FILE_NAME}")),
            move |header| {
                let header = match body_of(header, HEADER_FILE_NAME) {
                    Ok(bytes) => bytes,
                    Err(message) => {
                        let _ = sender.send(Loaded {
                            side,
                            source,
                            bytes: Err(message),
                        });
                        repaint();
                        return;
                    }
                };
                ehttp::fetch(ehttp::Request::get(frames_url), move |frames| {
                    let bytes =
                        body_of(frames, FRAME_FILE_NAME).map(|frames| RunBytes { header, frames });
                    let _ = sender.send(Loaded {
                        side,
                        source,
                        bytes,
                    });
                    repaint();
                });
            },
        );
    }

    /// A clone of the channel's sending half, for a source that runs on a
    /// thread of its own.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn sender(&self) -> Sender<Loaded> {
        self.sender.clone()
    }
}

/// The body of a response, or a message naming the file and what refused it.
fn body_of(response: ehttp::Result<ehttp::Response>, file_name: &str) -> Result<Vec<u8>, String> {
    let response =
        response.map_err(|error| format!("{file_name} could not be fetched: {error}"))?;
    if response.ok {
        Ok(response.bytes)
    } else {
        Err(format!(
            "{file_name} could not be fetched: {} {}",
            response.status, response.status_text
        ))
    }
}

/// Reading a run from a directory: the native convenience, and the only place
/// in the visualizer that knows what a path is (ADR-0006).
#[cfg(not(target_arch = "wasm32"))]
pub mod native {
    use std::path::Path;

    use termocline_format::{FRAME_FILE_NAME, HEADER_FILE_NAME};

    use crate::RunBytes;

    /// Read the run in `directory`: the two files of the format, by name.
    ///
    /// # Errors
    /// A message naming the file that could not be read and why.
    pub fn read_run_directory(directory: &Path) -> Result<RunBytes, String> {
        Ok(RunBytes {
            header: read_file(&directory.join(HEADER_FILE_NAME))?,
            frames: read_file(&directory.join(FRAME_FILE_NAME))?,
        })
    }

    /// One of the run's files, named in the error so a half-copied directory
    /// says which half is missing.
    fn read_file(path: &Path) -> Result<Vec<u8>, String> {
        std::fs::read(path)
            .map_err(|error| format!("{} could not be read: {error}", path.display()))
    }
}
