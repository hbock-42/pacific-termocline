//! A run arriving one dropped file at a time.
//!
//! A browser has no filesystem, so on the web a run is delivered as its two
//! files dragged onto the window (ADR-0006). They arrive in whatever order the
//! user drags them, and possibly in separate drops, so what is held here is a
//! half-assembled run: the files seen so far, matched by name against the two
//! the format defines, and handed over only once the pair is complete.

use termocline_format::{FRAME_FILE_NAME, HEADER_FILE_NAME};

use crate::RunBytes;

/// The files of a run seen so far.
#[derive(Debug, Default, Clone)]
pub struct PendingRun {
    /// The bytes of `header.json`, once dropped.
    header: Option<Vec<u8>>,
    /// The bytes of `frames.bin`, once dropped.
    frames: Option<Vec<u8>>,
}

impl PendingRun {
    /// Offer a dropped file, named by `file_name` and holding `bytes`.
    ///
    /// Returns whether the file is one of a run's two. `file_name` may be a
    /// full path — a native drop carries one — and only its final segment is
    /// matched, exactly: `Header.json` is not `header.json`, because the
    /// format names the file rather than describing it.
    ///
    /// A second drop of the same file replaces the first, so dropping the
    /// right header after the wrong one fixes the run.
    pub fn offer(&mut self, file_name: &str, bytes: Vec<u8>) -> bool {
        match base_name(file_name) {
            HEADER_FILE_NAME => self.header = Some(bytes),
            FRAME_FILE_NAME => self.frames = Some(bytes),
            _ => return false,
        }
        true
    }

    /// The files still missing, in the order they are named to the user.
    ///
    /// Empty exactly when [`PendingRun::take_run`] would return a run.
    #[must_use]
    pub fn still_needed(&self) -> Vec<&'static str> {
        let mut missing = Vec::new();
        if self.header.is_none() {
            missing.push(HEADER_FILE_NAME);
        }
        if self.frames.is_none() {
            missing.push(FRAME_FILE_NAME);
        }
        missing
    }

    /// Take the run, if both its files have arrived, leaving nothing behind.
    ///
    /// Clearing matters: files left here would pair the frames of one run with
    /// the header of the next, which reads as a truncated or overlong run
    /// rather than as the mistake it is.
    pub fn take_run(&mut self) -> Option<RunBytes> {
        let (header, frames) = (self.header.take(), self.frames.take());
        match (header, frames) {
            (Some(header), Some(frames)) => Some(RunBytes { header, frames }),
            (header, frames) => {
                self.header = header;
                self.frames = frames;
                None
            }
        }
    }
}

/// The final segment of `path`, under either platform's separator. Web drops
/// carry a bare name and native drops carry a path, and both reach here.
fn base_name(path: &str) -> &str {
    path.rsplit(['/', '\\']).next().unwrap_or(path)
}
