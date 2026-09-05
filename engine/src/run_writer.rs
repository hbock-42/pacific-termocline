//! Writing a run to disk: the header once, then frames at the output cadence.
//!
//! [`RunWriter`] is the engine's side of the file contract of [ADR-0001] and
//! [ADR-0004]. It is opened at the start of a run, writes the JSON
//! [`RunHeader`] immediately, and then takes one [`Frame`] per saved timestep
//! — *not* per timestep. A run long enough to be interesting is far longer
//! than a run small enough to hand to a visualizer, so the cadence is a
//! scenario parameter: [`OutputSchedule`] holds it, decides which steps are
//! written, and is the single place the frame count is derived.
//!
//! # Two files, and what is in each
//!
//! A run is a directory holding [`HEADER_FILE_NAME`] and [`FRAME_FILE_NAME`]:
//!
//! - the header is JSON, written once and never revisited, so a run that
//!   crashes mid-way still describes itself and a human can read its
//!   parameters with `cat`;
//! - the frames are `bincode`-encoded, one after another, with nothing
//!   between them and nothing after them.
//!
//! A run cut short by a crash therefore still describes its scenario — the
//! header is on disk from the first moment — and its frame file simply ends
//! early, which a reader counting the header's frames sees as the truncation
//! it is.
//!
//! Frames are encoded with [`termocline_format::frame_encoding`], which the
//! format crate owns rather than this one: nothing in the bytes says which
//! `bincode` configuration wrote them, so writer and reader have to be reading
//! the same line of the contract (ADR-0004).
//!
//! # Why the header's frame count is promised up front
//!
//! [`RunHeader`] carries the number of frames, and this writer writes the
//! header before the first of them exists. That is deliberate. Per [ADR-0006]
//! the reader is defined over a byte source — a file, an HTTP response, a
//! buffer dropped onto a browser window — and a byte source need not support
//! seeking. A reader that could not learn the frame count from the header
//! would have to seek to a trailer or decode until an error, and "decode until
//! it fails" cannot tell a finished run from a truncated one. So the count
//! comes from the [`OutputSchedule`], which knows it before the run starts,
//! and the writer holds itself to it: appending past the count is an error
//! ([`RunWriteError::TooManyFrames`]), and so is finishing below it
//! ([`RunWriteError::MissingFrames`]). A header that promised what the file
//! does not hold is exactly the silent corruption the checks exist to prevent.
//!
//! # Byte sinks, not paths
//!
//! For the same reason the reader is generic over a byte source, the writer is
//! generic over a byte sink: [`RunWriter::new`] takes any [`Write`], and
//! [`RunWriter::create`] is the convenience that turns a directory into two of
//! them. The bytes are identical either way, so a run assembled in memory and
//! a run on disk are the same run.
//!
//! [ADR-0001]: ../../docs/planning/adr/0001-engine-visualizer-split.md
//! [ADR-0004]: ../../docs/planning/adr/0004-data-interchange-format.md
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use std::fmt;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use termocline_format::{
    frame_encoding, FormatError, Frame, GridSpec, OutputTiming, RunHeader, FRAME_FILE_NAME,
    HEADER_FILE_NAME,
};

use crate::forcing::WindStressField;
use crate::params::PhysicalParams;
use crate::state::OceanState;

/// The engine's physical parameters as the format records them.
///
/// Both types are the same five SI quantities under the same names; this
/// conversion is the one place the run's parameters cross from the solver's
/// vocabulary into the file's, so a header can never disagree with the run it
/// describes.
impl From<PhysicalParams> for termocline_format::PhysicalParams {
    fn from(params: PhysicalParams) -> Self {
        Self {
            mean_depth_m: params.mean_thermocline_depth_m(),
            reduced_gravity_m_per_s2: params.reduced_gravity_m_per_s2(),
            beta_per_m_per_s: params.beta_per_m_per_s(),
            rayleigh_damping_per_s: params.rayleigh_damping_per_s(),
            reference_density_kg_per_m3: params.reference_density_kg_per_m3(),
        }
    }
}

/// Why an output schedule was rejected.
///
/// Both variants describe invalid *scenario input* — a run asking for a
/// cadence that does not describe an output series — so they are returned
/// rather than panicked, and each names the value it rejected.
#[derive(Debug, Clone, PartialEq)]
pub enum OutputScheduleError {
    /// The timestep was not a finite, positive duration in seconds.
    TimestepNotPositive {
        /// The timestep supplied, in seconds.
        dt_s: f64,
    },
    /// The output cadence was zero steps, which is not a cadence: it would
    /// ask for a frame between every pair of frames, without end.
    CadenceIsZero,
}

impl fmt::Display for OutputScheduleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TimestepNotPositive { dt_s } => {
                write!(f, "dt_s is {dt_s}; it must be finite and greater than 0")
            }
            Self::CadenceIsZero => write!(
                f,
                "every_n_steps is 0; a run writes a frame every N steps, and N must be at least 1"
            ),
        }
    }
}

impl std::error::Error for OutputScheduleError {}

/// How often a run writes a frame, and therefore how many it writes.
///
/// Output is *decimated*: a run steps at `dt_s` and saves every
/// `every_n_steps`-th step, starting from the initial state at step 0. The
/// saved series is then
///
/// ```text
/// steps 0, N, 2N, ...  up to and including the last multiple of N that the
///                      run reaches, which is total_steps / N frames plus the
///                      one at step 0
/// ```
///
/// so a run whose length is not a whole number of intervals simply stops at
/// the last interval that fits — the schedule never rounds the run up, and
/// never quietly moves a sample (CODING_STANDARDS.md § *No silent clamping*).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OutputSchedule {
    /// Length of one solver step, in seconds.
    dt_s: f64,
    /// Steps the run takes from its initial state.
    total_steps: u64,
    /// Steps between saved frames.
    every_n_steps: u64,
}

impl OutputSchedule {
    /// A run of `total_steps` steps of `dt_s` seconds, saving a frame every
    /// `every_n_steps` steps.
    ///
    /// # Errors
    /// [`OutputScheduleError::TimestepNotPositive`] if `dt_s` is not a finite,
    /// positive duration, and [`OutputScheduleError::CadenceIsZero`] if
    /// `every_n_steps` is zero.
    pub fn new(
        dt_s: f64,
        total_steps: u64,
        every_n_steps: u64,
    ) -> Result<Self, OutputScheduleError> {
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(OutputScheduleError::TimestepNotPositive { dt_s });
        }
        if every_n_steps == 0 {
            return Err(OutputScheduleError::CadenceIsZero);
        }
        Ok(Self {
            dt_s,
            total_steps,
            every_n_steps,
        })
    }

    /// Length of one solver step, in seconds.
    #[must_use]
    pub const fn dt_s(self) -> f64 {
        self.dt_s
    }

    /// Steps the run takes from its initial state.
    #[must_use]
    pub const fn total_steps(self) -> u64 {
        self.total_steps
    }

    /// Model time between consecutive frames, in seconds — the *output*
    /// interval, which is the timestep only when every step is saved.
    #[must_use]
    pub fn interval_s(self) -> f64 {
        self.every_n_steps as f64 * self.dt_s
    }

    /// Frames this run will write: one at step 0, then one per interval that
    /// fits inside the run.
    #[must_use]
    pub const fn frame_count(self) -> u64 {
        self.total_steps / self.every_n_steps + 1
    }

    /// Whether the state after `step` steps is one of the saved ones.
    #[must_use]
    pub const fn writes_at_step(self, step: u64) -> bool {
        step <= self.total_steps && step.is_multiple_of(self.every_n_steps)
    }

    /// This schedule as the header records it.
    #[must_use]
    pub fn timing(self) -> OutputTiming {
        OutputTiming {
            frame_count: self.frame_count(),
            interval_s: self.interval_s(),
        }
    }
}

/// Why a run could not be written.
///
/// Every variant describes something the *caller* got wrong or something the
/// filesystem refused — never a broken invariant of the writer — so they are
/// returned rather than panicked, and each names what it was expecting.
#[derive(Debug)]
pub enum RunWriteError {
    /// The run directory or one of its two files could not be opened.
    Io(std::io::Error),
    /// The header could not be encoded or written.
    Header(serde_json::Error),
    /// A frame did not fit the grid the header describes.
    Frame(FormatError),
    /// A frame could not be encoded or written.
    Encode(bincode::error::EncodeError),
    /// A frame was appended after the run had written every frame its header
    /// promised. A reader takes the count from the header, so the extra frame
    /// would be unreadable — the run is stopped instead of silently growing a
    /// tail nothing will read.
    TooManyFrames {
        /// Frames the header promised.
        promised: u64,
    },
    /// The run finished before writing every frame its header promised, which
    /// would leave a reader reading off the end of the frame file.
    MissingFrames {
        /// Frames the header promised.
        promised: u64,
        /// Frames actually written.
        written: u64,
    },
}

impl fmt::Display for RunWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => write!(f, "the run's files could not be written: {error}"),
            Self::Header(error) => write!(f, "the run header could not be written: {error}"),
            Self::Frame(error) => error.fmt(f),
            Self::Encode(error) => write!(f, "a frame could not be encoded: {error}"),
            Self::TooManyFrames { promised } => write!(
                f,
                "the header promises {promised} frames and {promised} have been written; \
                 a further frame would be one no reader ever sees"
            ),
            Self::MissingFrames { promised, written } => write!(
                f,
                "the header promises {promised} frames but the run wrote {written}; \
                 finish the run's output schedule, or write a header that promises {written}"
            ),
        }
    }
}

impl std::error::Error for RunWriteError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::Header(error) => Some(error),
            Self::Frame(error) => Some(error),
            Self::Encode(error) => Some(error),
            Self::TooManyFrames { .. } | Self::MissingFrames { .. } => None,
        }
    }
}

impl From<std::io::Error> for RunWriteError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<serde_json::Error> for RunWriteError {
    fn from(error: serde_json::Error) -> Self {
        Self::Header(error)
    }
}

impl From<FormatError> for RunWriteError {
    fn from(error: FormatError) -> Self {
        Self::Frame(error)
    }
}

impl From<bincode::error::EncodeError> for RunWriteError {
    fn from(error: bincode::error::EncodeError) -> Self {
        Self::Encode(error)
    }
}

/// An open run: the header already written, the frames still arriving.
///
/// Built at the start of a run from the [`RunHeader`] that describes it, fed
/// one saved timestep at a time through [`RunWriter::append`], and closed with
/// [`RunWriter::finish`], which is where a run that wrote the wrong number of
/// frames is caught.
#[derive(Debug)]
pub struct RunWriter<W: Write> {
    /// Where encoded frames go, in order, with nothing between them.
    frame_sink: W,
    /// The grid the header describes; every frame is checked against it.
    grid: GridSpec,
    /// Frames the header promised.
    promised: u64,
    /// Frames written so far.
    written: u64,
}

impl<W: Write> RunWriter<W> {
    /// Open a run over two byte sinks, writing `header` into the first of them
    /// before returning.
    ///
    /// `header_sink` is consumed: the header is written once, at the start of
    /// the run, and never revisited.
    ///
    /// # Errors
    /// [`RunWriteError::Header`] if the header could not be encoded or written.
    pub fn new<H: Write>(
        mut header_sink: H,
        frame_sink: W,
        header: &RunHeader,
    ) -> Result<Self, RunWriteError> {
        serde_json::to_writer(&mut header_sink, header)?;
        // A trailing newline so the file reads as a line of text in a
        // terminal; `serde_json` ignores it on the way back in.
        header_sink.write_all(b"\n")?;
        header_sink.flush()?;
        Ok(Self {
            frame_sink,
            grid: header.grid,
            promised: header.output.frame_count,
            written: 0,
        })
    }

    /// Append the saved timestep at model time `t_s`.
    ///
    /// `state` supplies `h`, `u` and `v`; `wind_stress` supplies the forcing
    /// `τx` and `τy` that drove the run to it, in pascals, which is the
    /// `N m^-2` the format states.
    ///
    /// A frame is built and encoded per call, so this allocates — which is why
    /// it belongs at the output cadence and not inside the time-stepping loop
    /// (CODING_STANDARDS.md § *Performance*).
    ///
    /// # Errors
    /// [`RunWriteError::Frame`] if the state or the stress does not cover the
    /// basin the header describes, [`RunWriteError::TooManyFrames`] if the
    /// header's frame count has already been written, and
    /// [`RunWriteError::Encode`] if the frame could not be encoded or written.
    pub fn append(
        &mut self,
        t_s: f64,
        state: &OceanState,
        wind_stress: &WindStressField,
    ) -> Result<(), RunWriteError> {
        if self.written == self.promised {
            return Err(RunWriteError::TooManyFrames {
                promised: self.promised,
            });
        }
        let frame = Frame::new(
            t_s,
            &self.grid,
            state.h().as_slice().to_vec(),
            state.u().as_slice().to_vec(),
            state.v().as_slice().to_vec(),
            wind_stress.tau_x_pa().as_slice().to_vec(),
            wind_stress.tau_y_pa().as_slice().to_vec(),
        )?;
        bincode::serde::encode_into_std_write(&frame, &mut self.frame_sink, frame_encoding())?;
        self.written += 1;
        Ok(())
    }

    /// Close the run, returning the frame sink.
    ///
    /// # Errors
    /// [`RunWriteError::MissingFrames`] if the run wrote fewer frames than its
    /// header promised, and [`RunWriteError::Io`] if the last of them could
    /// not be flushed.
    pub fn finish(mut self) -> Result<W, RunWriteError> {
        if self.written != self.promised {
            return Err(RunWriteError::MissingFrames {
                promised: self.promised,
                written: self.written,
            });
        }
        self.frame_sink.flush()?;
        Ok(self.frame_sink)
    }
}

impl RunWriter<BufWriter<File>> {
    /// Open a run as two files in `directory`, creating it if it does not
    /// exist: [`HEADER_FILE_NAME`] and [`FRAME_FILE_NAME`].
    ///
    /// The native convenience over [`RunWriter::new`]; the bytes it writes are
    /// the same ones any other sink would receive.
    ///
    /// # Errors
    /// [`RunWriteError::Io`] if the directory or either file could not be
    /// created, and the errors of [`RunWriter::new`].
    pub fn create(directory: &Path, header: &RunHeader) -> Result<Self, RunWriteError> {
        fs::create_dir_all(directory)?;
        let header_file = BufWriter::new(File::create(directory.join(HEADER_FILE_NAME))?);
        let frame_file = BufWriter::new(File::create(directory.join(FRAME_FILE_NAME))?);
        Self::new(header_file, frame_file, header)
    }
}
