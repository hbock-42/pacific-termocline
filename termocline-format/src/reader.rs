//! Reading a run back: the header once, then frames one at a time.
//!
//! [`RunReader`] is the other side of the engine's writer. It takes the two
//! byte sources a run is made of, decodes the JSON [`RunHeader`] eagerly —
//! everything downstream needs it before it can do anything with a frame — and
//! then hands out [`Frame`]s lazily, one per call to [`Iterator::next`].
//!
//! # Byte sources, not paths
//!
//! Per [ADR-0006] the visualizer runs in a browser, where there is no
//! filesystem: a run arrives by file selection, by drag-and-drop, or over
//! HTTP. So the reader is defined over any [`Read`], and
//! [`RunReader::open`] — the only thing here that knows what a path is — sits
//! behind the `fs` feature as a native convenience. The bytes are the same
//! either way, so a run fetched into a browser buffer and a run on local disk
//! read back as the same run.
//!
//! Reading is strictly forward: no seek, no index, no trailer. That is a
//! stronger constraint than ADR-0006 asks for, and it is what lets a run be
//! read while it is still arriving over a network.
//!
//! # Why the frames are lazy
//!
//! A run at a realistic resolution and duration is far larger than the frame
//! being drawn. Decoding the whole series to return it would put the entire
//! run in memory at once — tolerable for a local file, not for a browser tab
//! holding a few hundred megabytes of it. So the reader keeps one frame at a
//! time and its memory does not grow with the length of the run;
//! `tests/reader_memory.rs` is that claim as an assertion.
//!
//! # What the header is trusted for
//!
//! Everything the frames do not say themselves. The bytes of a frame record
//! neither the grid nor the number of frames beside them, so the reader takes
//! both from the header and holds the file to them: a frame whose fields do
//! not cover the header's basin is [`RunReadError::Frame`], a run whose frames
//! stop early is [`RunReadError::Truncated`], and bytes past the promised
//! count are [`RunReadError::TrailingBytes`] — the mirror of the writer's
//! refusal to append past them. A truncated run that ended the iteration
//! quietly would look exactly like a complete one.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use std::fmt;
use std::io::Read;

use crate::{frame_encoding, FormatError, Frame, RunHeader, FORMAT_VERSION};

/// Why a run could not be read.
///
/// Every variant describes invalid *input* — a header from a version this
/// build does not know, a frame that does not fit the basin its header
/// describes, a stream that ended early — so they are returned rather than
/// panicked, and each names what it was expecting.
#[derive(Debug)]
pub enum RunReadError {
    /// The header could not be read, or was not the JSON a header is.
    Header(serde_json::Error),
    /// The run was written by a version of the format this build does not
    /// understand. Its frames may have any layout, so they are not decoded.
    UnsupportedVersion {
        /// Format version the run's header carries.
        found: u32,
        /// Format version this build reads.
        supported: u32,
    },
    /// A frame's bytes could not be decoded.
    Decode(bincode::error::DecodeError),
    /// A frame decoded, but did not fit the grid the header describes.
    Frame(FormatError),
    /// The frame source ended before the header's frame count was reached:
    /// the run was cut short, or its two files do not belong together.
    Truncated {
        /// Frames the header promised.
        promised: u64,
        /// Frames actually read.
        read: u64,
    },
    /// The frame source held bytes past the header's frame count. The header
    /// is what a reader counts by, so those bytes are data nothing would ever
    /// read; the run is refused rather than silently truncated to fit.
    TrailingBytes {
        /// Frames the header promised.
        promised: u64,
    },
    /// The frame source could not be read past the last complete frame.
    Io(std::io::Error),
}

impl fmt::Display for RunReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Header(error) => write!(f, "the run header could not be read: {error}"),
            Self::UnsupportedVersion { found, supported } => write!(
                f,
                "the run is in format version {found} and this build reads version {supported}"
            ),
            Self::Decode(error) => write!(f, "a frame could not be decoded: {error}"),
            Self::Frame(error) => error.fmt(f),
            Self::Truncated { promised, read } => write!(
                f,
                "the header promises {promised} frames but the frame source ended after {read}; \
                 the run was cut short, or its header belongs to another run"
            ),
            Self::TrailingBytes { promised } => write!(
                f,
                "the header promises {promised} frames and the frame source holds more; \
                 the run is longer than its header describes"
            ),
            Self::Io(error) => write!(f, "the run's frames could not be read: {error}"),
        }
    }
}

impl std::error::Error for RunReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Header(error) => Some(error),
            Self::Decode(error) => Some(error),
            Self::Frame(error) => Some(error),
            Self::Io(error) => Some(error),
            Self::UnsupportedVersion { .. }
            | Self::Truncated { .. }
            | Self::TrailingBytes { .. } => None,
        }
    }
}

impl From<serde_json::Error> for RunReadError {
    fn from(error: serde_json::Error) -> Self {
        Self::Header(error)
    }
}

impl From<FormatError> for RunReadError {
    fn from(error: FormatError) -> Self {
        Self::Frame(error)
    }
}

impl From<std::io::Error> for RunReadError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

/// An open run: the header already decoded, the frames still to come.
///
/// Built from the run's two byte sources with [`RunReader::new`] (or
/// [`RunReader::open`] on a native filesystem), asked about the run through
/// [`RunReader::header`], and iterated for its frames — which arrive in the
/// order they were written, one decode at a time.
///
/// Iteration yields `Result`s because a byte source can fail or end early. An
/// error is terminal: the stream's position is no longer known, so the reader
/// yields `None` from then on rather than trying to resynchronize on bytes it
/// cannot interpret.
#[derive(Debug)]
pub struct RunReader<R: Read> {
    /// Everything the frames do not say about themselves.
    header: RunHeader,
    /// Encoded frames, in order, with nothing between them.
    frames: R,
    /// Frames decoded so far.
    read: u64,
    /// Whether the bytes past the last promised frame have been looked at.
    tail_checked: bool,
    /// Whether an error has ended the iteration.
    failed: bool,
}

impl<R: Read> RunReader<R> {
    /// Open a run over its two byte sources, decoding `header_source` before
    /// returning.
    ///
    /// `frame_source` is not touched until the first frame is asked for, so a
    /// caller may inspect the header — and refuse the run — without reading a
    /// byte of it.
    ///
    /// # Errors
    /// [`RunReadError::Header`] if the header could not be read or is not
    /// valid JSON for a [`RunHeader`], and
    /// [`RunReadError::UnsupportedVersion`] if it was written by a format
    /// version this build does not understand.
    pub fn new<H: Read>(header_source: H, frame_source: R) -> Result<Self, RunReadError> {
        let header: RunHeader = serde_json::from_reader(header_source)?;
        if header.format_version != FORMAT_VERSION {
            return Err(RunReadError::UnsupportedVersion {
                found: header.format_version,
                supported: FORMAT_VERSION,
            });
        }
        Ok(Self {
            header,
            frames: frame_source,
            read: 0,
            tail_checked: false,
            failed: false,
        })
    }

    /// The run's header: its grid, its physical parameters, its variables and
    /// its output cadence.
    #[must_use]
    pub const fn header(&self) -> &RunHeader {
        &self.header
    }

    /// Frames the header promises that have not been read yet.
    ///
    /// Zero once the run is exhausted, and zero after an error, which ends the
    /// iteration.
    #[must_use]
    pub const fn remaining_frames(&self) -> u64 {
        if self.failed {
            0
        } else {
            self.header.output.frame_count - self.read
        }
    }

    /// Decode the next frame, whatever the counters say.
    fn decode_frame(&mut self) -> Result<Frame, RunReadError> {
        let frame: Frame = bincode::serde::decode_from_std_read(&mut self.frames, frame_encoding())
            .map_err(|error| {
                if ended_early(&error) {
                    RunReadError::Truncated {
                        promised: self.header.output.frame_count,
                        read: self.read,
                    }
                } else {
                    RunReadError::Decode(error)
                }
            })?;
        frame.validate(&self.header.grid)?;
        Ok(frame)
    }

    /// Look for bytes past the last frame the header promised, once.
    fn check_tail(&mut self) -> Option<Result<Frame, RunReadError>> {
        self.tail_checked = true;
        let mut byte = [0_u8; 1];
        match self.frames.read(&mut byte) {
            Ok(0) => None,
            Ok(_) => {
                self.failed = true;
                Some(Err(RunReadError::TrailingBytes {
                    promised: self.header.output.frame_count,
                }))
            }
            Err(error) => {
                self.failed = true;
                Some(Err(RunReadError::Io(error)))
            }
        }
    }
}

impl<R: Read> Iterator for RunReader<R> {
    type Item = Result<Frame, RunReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed {
            return None;
        }
        if self.read == self.header.output.frame_count {
            return if self.tail_checked {
                None
            } else {
                self.check_tail()
            };
        }
        match self.decode_frame() {
            Ok(frame) => {
                self.read += 1;
                Some(Ok(frame))
            }
            Err(error) => {
                self.failed = true;
                Some(Err(error))
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let Ok(remaining) = usize::try_from(self.remaining_frames()) else {
            // More frames than this target can count; the honest hint is that
            // there are at least some and no known bound.
            return (0, None);
        };
        // One more may follow: the trailing-byte check yields an error rather
        // than a frame, and it has not run yet.
        let upper = if self.tail_checked || self.failed {
            Some(remaining)
        } else {
            None
        };
        (remaining, upper)
    }
}

/// Whether a decode failure was the frame source running out, rather than the
/// bytes being wrong.
fn ended_early(error: &bincode::error::DecodeError) -> bool {
    match error {
        bincode::error::DecodeError::UnexpectedEnd { .. } => true,
        bincode::error::DecodeError::Io { inner, .. } => {
            inner.kind() == std::io::ErrorKind::UnexpectedEof
        }
        _ => false,
    }
}

#[cfg(feature = "fs")]
mod fs {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;

    use super::{RunReadError, RunReader};
    use crate::{FRAME_FILE_NAME, HEADER_FILE_NAME};

    impl RunReader<BufReader<File>> {
        /// Open the run in `directory`: [`HEADER_FILE_NAME`] beside
        /// [`FRAME_FILE_NAME`].
        ///
        /// The native convenience over [`RunReader::new`], behind the `fs`
        /// feature because a browser has no filesystem to offer it
        /// (ADR-0006). The bytes it reads are the ones any other source would
        /// supply.
        ///
        /// # Errors
        /// [`RunReadError::Io`] if either file could not be opened, and the
        /// errors of [`RunReader::new`].
        pub fn open(directory: &Path) -> Result<Self, RunReadError> {
            let header = BufReader::new(File::open(directory.join(HEADER_FILE_NAME))?);
            let frames = BufReader::new(File::open(directory.join(FRAME_FILE_NAME))?);
            Self::new(header, frames)
        }
    }
}
