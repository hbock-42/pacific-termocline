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
//! A caller that already holds the whole run in memory needs the other thing —
//! any frame, in any order, because a scrubber is dragged rather than played
//! (T-08.3). It gets it by building an index of its own: [`decode_frame`]
//! reads the frame at the start of a slice and says how many bytes it took, so
//! one forward pass learns where every frame begins and each later fetch costs
//! one decode. The index is the caller's, held for as long as it holds the
//! bytes; nothing about the file on disk changes.
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

use crate::frame::FrameV1;
use crate::{
    frame_encoding, FormatError, Frame, RunHeader, FORMAT_VERSION, OLDEST_READABLE_FORMAT_VERSION,
};

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
        /// Oldest format version this build reads.
        oldest_supported: u32,
        /// Newest format version this build reads — the one it writes.
        newest_supported: u32,
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
    /// One of the run's two files could not be opened. Only
    /// [`RunReader::open`] produces this: a byte source handed in by a caller
    /// is already open.
    Open {
        /// The file that could not be opened.
        path: std::path::PathBuf,
        /// Why the filesystem refused it.
        error: std::io::Error,
    },
}

impl fmt::Display for RunReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Header(error) => write!(f, "the run header could not be read: {error}"),
            Self::UnsupportedVersion {
                found,
                oldest_supported,
                newest_supported,
            } => write!(
                f,
                "the run is in format version {found} and this build reads versions \
                 {oldest_supported} to {newest_supported}"
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
            Self::Open { path, error } => {
                write!(f, "{} could not be opened: {error}", path.display())
            }
        }
    }
}

impl std::error::Error for RunReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Header(error) => Some(error),
            Self::Decode(error) => Some(error),
            Self::Frame(error) => Some(error),
            Self::Io(error) | Self::Open { error, .. } => Some(error),
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

/// Decode the frame that begins at the start of `bytes`, and say how many
/// bytes it occupied.
///
/// The random-access half of reading a run, for a caller holding the frames
/// whole in memory: the length it returns is what lets one forward pass note
/// where every frame begins, and an offset from that pass is what lets a later
/// frame be fetched without decoding the frames before it. Bytes past the
/// frame are ignored, so the slice may be the rest of the run.
///
/// The frame is held to `header` exactly as [`RunReader`] holds it: a frame is
/// bytes until something says what shape they should be, and that is the
/// header's job either way. `header` rather than a bare grid because the grid
/// is only half of that — the format version says which frame layout the bytes
/// are in, and version 1 has no room in it for `T'`.
///
/// # Errors
/// [`RunReadError::UnsupportedVersion`] if the header is from a version this
/// build does not read, [`RunReadError::Decode`] if the bytes are not a frame,
/// or run out part-way through one — a caller decoding a single frame is not
/// counting against a promised total, so nothing here is
/// [`RunReadError::Truncated`] — and [`RunReadError::Frame`] if the frame does
/// not fit the header's grid.
pub fn decode_frame(bytes: &[u8], header: &RunHeader) -> Result<(Frame, usize), RunReadError> {
    let layout = FrameLayout::of_version(header.format_version)?;
    let (frame, used) = layout
        .decode_from_slice(bytes)
        .map_err(RunReadError::Decode)?;
    frame.validate(&header.grid)?;
    Ok((frame, used))
}

/// Which shape a run's frames are in on disk.
///
/// The format version reduced to the only thing a frame decoder needs from it.
/// Every readable version has a variant, so a version that reached this far is
/// a version there are bytes-to-frame instructions for; there is no default
/// case that would guess.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FrameLayout {
    /// Version 1: the five fields of the linear core.
    LinearCore,
    /// Version 2: the linear core followed by the optional SST anomaly.
    WithOptionalSstAnomaly,
}

impl FrameLayout {
    /// The layout format version `version` writes its frames in.
    ///
    /// # Errors
    /// [`RunReadError::UnsupportedVersion`] for any version outside
    /// [`OLDEST_READABLE_FORMAT_VERSION`]..=[`FORMAT_VERSION`].
    fn of_version(version: u32) -> Result<Self, RunReadError> {
        match version {
            1 => Ok(Self::LinearCore),
            2 => Ok(Self::WithOptionalSstAnomaly),
            found => Err(RunReadError::UnsupportedVersion {
                found,
                oldest_supported: OLDEST_READABLE_FORMAT_VERSION,
                newest_supported: FORMAT_VERSION,
            }),
        }
    }

    /// Decode the frame at the start of `bytes`, and say how many bytes it
    /// occupied.
    fn decode_from_slice(
        self,
        bytes: &[u8],
    ) -> Result<(Frame, usize), bincode::error::DecodeError> {
        match self {
            Self::LinearCore => {
                bincode::serde::decode_from_slice::<FrameV1, _>(bytes, frame_encoding())
                    .map(|(frame, used)| (frame.into(), used))
            }
            Self::WithOptionalSstAnomaly => {
                bincode::serde::decode_from_slice(bytes, frame_encoding())
            }
        }
    }

    /// Decode the next frame from `source`.
    fn decode_from_read<R: Read>(
        self,
        source: &mut R,
    ) -> Result<Frame, bincode::error::DecodeError> {
        match self {
            Self::LinearCore => {
                bincode::serde::decode_from_std_read::<FrameV1, _, _>(source, frame_encoding())
                    .map(Frame::from)
            }
            Self::WithOptionalSstAnomaly => {
                bincode::serde::decode_from_std_read(source, frame_encoding())
            }
        }
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
    /// The shape the run's frames are in, from the header's format version.
    layout: FrameLayout,
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
    /// version outside the range this build reads.
    pub fn new<H: Read>(header_source: H, frame_source: R) -> Result<Self, RunReadError> {
        let header: RunHeader = serde_json::from_reader(header_source)?;
        let layout = FrameLayout::of_version(header.format_version)?;
        Ok(Self {
            header,
            layout,
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
    /// A promise, not a guarantee: it is the header's frame count less what
    /// has been read, and a run cut short holds fewer (which is
    /// [`RunReadError::Truncated`] when the iteration reaches the gap). Zero
    /// once the run is exhausted, and zero after an error, which ends the
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
        let frame: Frame = self
            .layout
            .decode_from_read(&mut self.frames)
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

    /// Look for bytes past the last frame the header promised, once. Returns
    /// what is wrong with them, if there are any.
    fn tail_error(&mut self) -> Option<RunReadError> {
        self.tail_checked = true;
        let mut byte = [0_u8; 1];
        let error = match self.frames.read(&mut byte) {
            Ok(0) => return None,
            Ok(_) => RunReadError::TrailingBytes {
                promised: self.header.output.frame_count,
            },
            Err(error) => RunReadError::Io(error),
        };
        self.failed = true;
        Some(error)
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
                self.tail_error().map(Err)
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
        // The lower bound stays zero however many frames the header promises.
        // `collect` reserves the lower bound up front, and the header is a
        // claim about a file that may not keep it: a truncated run yields
        // fewer frames than it promises, and a header claiming a billion
        // frames beside a two-kilobyte file would otherwise reserve gigabytes
        // before the first frame was decoded. The count a caller can act on is
        // `remaining_frames`, which says where it came from.
        //
        // The upper bound is honest in the other direction, and only once the
        // run is over: until the trailing-byte check has run, one more item —
        // an error rather than a frame — may still follow.
        let upper = if self.failed {
            Some(0)
        } else if self.tail_checked {
            usize::try_from(self.remaining_frames()).ok()
        } else {
            None
        };
        (0, upper)
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
        /// [`RunReadError::Open`], naming the file, if either could not be
        /// opened, and the errors of [`RunReader::new`].
        pub fn open(directory: &Path) -> Result<Self, RunReadError> {
            let header = open_file(&directory.join(HEADER_FILE_NAME))?;
            let frames = open_file(&directory.join(FRAME_FILE_NAME))?;
            Self::new(header, frames)
        }
    }

    /// One of the run's two files, opened by name so a half-copied run
    /// directory says which half is missing.
    fn open_file(path: &Path) -> Result<BufReader<File>, RunReadError> {
        File::open(path)
            .map(BufReader::new)
            .map_err(|error| RunReadError::Open {
                path: path.to_owned(),
                error,
            })
    }
}
