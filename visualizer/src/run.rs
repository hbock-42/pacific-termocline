//! A run, once its bytes have arrived, and the metadata the shell shows for it.
//!
//! The shell's job in T-08.1 is to load a run and say what it is. Everything
//! about that is here, and none of it knows what a window or a GPU is: a run
//! is two byte buffers in, a labelled list of rows out. Where the bytes came
//! from — a directory, a pair of dropped files, an HTTP fetch — is the
//! caller's problem, which is what lets the same code serve the browser and
//! the native build (ADR-0006).

use std::cell::Cell;
use std::fmt;
use std::io::Read;

use termocline_format::{frame_encoding, Frame, RunHeader, RunReadError, RunReader};

use crate::DivergingScale;

/// Seconds in a day, for reporting a run's model time in the unit its output
/// cadence is chosen in (`steady-trades.toml` writes a frame a day).
pub(crate) const SECONDS_PER_DAY: f64 = 86_400.0;

/// The two byte sources a run is made of, however they arrived.
///
/// Owned buffers rather than streams: a run that arrived by drop or by fetch
/// is already whole in memory, and a native directory is read into the same
/// shape so that one loader serves both.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct RunBytes {
    /// The JSON [`RunHeader`], as read from `header.json`.
    pub header: Vec<u8>,
    /// The encoded frames, as read from `frames.bin`.
    pub frames: Vec<u8>,
}

/// One line of the metadata panel: what it is, and what this run says it is.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MetadataRow {
    /// What the value means, as shown to the reader.
    pub label: &'static str,
    /// The value, already formatted with its unit.
    pub value: String,
}

impl MetadataRow {
    /// A row labelled `label` showing `value`.
    fn new(label: &'static str, value: impl Into<String>) -> Self {
        Self {
            label,
            value: value.into(),
        }
    }
}

/// A run the shell has loaded: its header decoded, its frames still bytes.
///
/// The frames are kept rather than dropped because on the web they cannot be
/// re-read: the bytes a user dropped or a fetch returned are the only copy,
/// and any frame of the run may be the one asked for next
/// ([`LoadedRun::frame`]).
#[derive(Debug, Clone)]
pub struct LoadedRun {
    /// Where the run came from, for a reader checking they opened the one they
    /// meant to.
    source: String,
    /// Everything the frames do not say about themselves.
    header: RunHeader,
    /// The run's encoded frames, undecoded.
    frames: Vec<u8>,
    /// Where each frame begins in `frames`, by index into the run.
    ///
    /// The format is forward-only by design (`termocline_format::reader`), and
    /// this is what buys random access back without asking it to seek: the
    /// offsets are noted on the pass over the frames that load already makes,
    /// and cost eight bytes a frame — 6 KB for a run of 731 — against a frame
    /// of the scenario's basin, which is 1.3 MB. A scrubber is a random-access
    /// control, so without them dragging to the end of a run would decode
    /// every frame of it (T-08.3).
    frame_offsets: Vec<usize>,
    /// A colour scale covering every frame of the run, built on the pass that
    /// counts them.
    scale: DivergingScale,
}

impl LoadedRun {
    /// Decode `bytes`, which came from `source`.
    ///
    /// The header goes through [`RunReader`] rather than `serde_json` directly,
    /// so a run written by a format version this build does not read is
    /// refused here rather than mislabelled in the UI.
    ///
    /// Every frame is decoded once here and thrown away, and where each one
    /// began is kept. That costs a pass over the run at load, and it buys two
    /// things the panel cannot get from the header alone: the frame count it
    /// shows is a number the file keeps rather than a number the header claims
    /// — a `header.json` beside another run's `frames.bin` reads as a
    /// perfectly plausible run until the frames are counted — and the offsets
    /// that make [`LoadedRun::frame`] cost the same for every frame of the
    /// run. The pass holds one frame at a time, so it costs time and (bar the
    /// offsets, eight bytes a frame) not memory — which is the resource a
    /// browser tab is short of (ADR-0006).
    ///
    /// # Errors
    /// The errors of [`RunReader`]: a header that is not valid JSON for a
    /// [`RunHeader`], one from an unsupported [format version], frames that do
    /// not fit the basin the header describes, or a frame count the file does
    /// not keep.
    ///
    /// [format version]: termocline_format::FORMAT_VERSION
    pub fn from_bytes(source: impl Into<String>, bytes: RunBytes) -> Result<Self, RunReadError> {
        let RunBytes {
            header: header_source,
            frames,
        } = bytes;
        // Taken by value, and the frames moved rather than copied: a run in a
        // browser tab has no second copy to spare (ADR-0006).
        let consumed = Cell::new(0);
        let mut reader = RunReader::new(
            header_source.as_slice(),
            Counting {
                inner: frames.as_slice(),
                consumed: &consumed,
            },
        )?;
        let header = reader.header().clone();
        let mut scale = DivergingScale::symmetric_over(&[]);
        let mut frame_offsets = Vec::new();
        let mut offset = 0;
        for frame in reader.by_ref() {
            let frame = frame?;
            frame_offsets.push(offset);
            offset = consumed.get();
            scale = scale.widened(DivergingScale::symmetric_over(frame.h()));
        }
        Ok(Self {
            source: source.into(),
            header,
            frames,
            frame_offsets,
            scale,
        })
    }

    /// A colour scale for the run's thermocline depth anomaly, reaching as far
    /// either side of zero as the largest anomaly anywhere in the run.
    ///
    /// One scale for the whole run rather than one per frame: the same colour
    /// then means the same anomaly wherever it is seen, so a tilt that
    /// collapses over a run is seen to collapse instead of being renormalized
    /// back to full saturation frame by frame.
    #[must_use]
    pub const fn anomaly_scale(&self) -> DivergingScale {
        self.scale
    }

    /// Where the run came from: a directory, a pair of dropped files, or a URL.
    #[must_use]
    pub fn source(&self) -> &str {
        &self.source
    }

    /// The run's header.
    #[must_use]
    pub const fn header(&self) -> &RunHeader {
        &self.header
    }

    /// Frame number `index`, counting from zero, or `None` past the end of the
    /// run.
    ///
    /// One decode, wherever the frame sits in the run: the frame's bytes are
    /// found by the offset noted for it at load, and the frames before it are
    /// never touched. That is what a scrubber needs — a drag lands on frame
    /// 600 without passing through the 599 before it — and it is why dragging
    /// costs the same at both ends of the run (T-08.3). The decoded frame is
    /// not kept: a run in a browser tab has no room for a second copy of
    /// itself (ADR-0006), and the shell caches the one map it is drawing.
    ///
    /// # Panics
    /// If a frame this run already decoded at load does not decode again. That
    /// is not bad input — [`LoadedRun::from_bytes`] refused bad input — but
    /// this code disagreeing with itself.
    #[must_use]
    pub fn frame(&self, index: u64) -> Option<Frame> {
        let index = usize::try_from(index).ok()?;
        let offset = *self.frame_offsets.get(index)?;
        let (frame, _bytes) =
            bincode::serde::decode_from_slice::<Frame, _>(&self.frames[offset..], frame_encoding())
                .expect("every frame decoded at load");
        Some(frame)
    }

    /// What the shell shows about this run, in the order it shows it.
    #[must_use]
    pub fn metadata(&self) -> Vec<MetadataRow> {
        let grid = self.header.grid;
        let extent = grid.extent();
        let output = self.header.output;
        let params = self.header.physical_params;
        // The header's own symbols, not this build's: a run says what its
        // frames carry, and the panel reports what the run says.
        let variables: Vec<&str> = self
            .header
            .variables
            .iter()
            .map(|spec| spec.symbol.as_str())
            .collect();
        vec![
            MetadataRow::new("Scenario", self.header.scenario_description.clone()),
            MetadataRow::new("Grid", format!("{} × {} cells", grid.nx(), grid.ny())),
            MetadataRow::new("Frames", output.frame_count.to_string()),
            MetadataRow::new(
                "Basin",
                format!(
                    "{}, {}",
                    longitude_span(extent.west_deg_east, extent.east_deg_east),
                    latitude_span(extent.south_deg_north, extent.north_deg_north),
                ),
            ),
            MetadataRow::new(
                "Frame interval",
                format!(
                    "{} s ({:.2} days)",
                    output.interval_s,
                    output.interval_s / SECONDS_PER_DAY
                ),
            ),
            MetadataRow::new("Model time", format!("{:.2} days", self.model_time_days())),
            MetadataRow::new(
                "Mean depth H",
                format!("{} m", format_quantity(params.mean_depth_m)),
            ),
            MetadataRow::new(
                "Reduced gravity g'",
                format!(
                    "{} m s^-2",
                    format_quantity(params.reduced_gravity_m_per_s2)
                ),
            ),
            MetadataRow::new(
                "Coriolis gradient β",
                format!("{} m^-1 s^-1", format_quantity(params.beta_per_m_per_s)),
            ),
            MetadataRow::new(
                "Rayleigh damping r",
                format!("{} s^-1", format_quantity(params.rayleigh_damping_per_s)),
            ),
            MetadataRow::new(
                "Reference density ρ₀",
                format!(
                    "{} kg m^-3",
                    format_quantity(params.reference_density_kg_per_m3)
                ),
            ),
            MetadataRow::new("Variables", variables.join(", ")),
            MetadataRow::new("Format version", self.header.format_version.to_string()),
        ]
    }

    /// Model time from the first frame to the last, in days.
    ///
    /// The span of `n` frames is `n - 1` intervals, not `n`: a run of one frame
    /// spans no time at all, and a run of none spans none either.
    fn model_time_days(&self) -> f64 {
        let output = self.header.output;
        #[allow(clippy::cast_precision_loss)]
        let intervals = output.frame_count.saturating_sub(1) as f64;
        intervals * output.interval_s / SECONDS_PER_DAY
    }
}

/// A byte source that remembers how much of itself has been handed out.
///
/// The reader is forward-only and says nothing about where it has got to, so
/// this is how the offset of each frame is learned while the frames are being
/// walked at load anyway: read the counter between frames, and the difference
/// is the frame that just went past.
struct Counting<'a, R: Read> {
    /// The bytes themselves.
    inner: R,
    /// Bytes read out of `inner` so far.
    consumed: &'a Cell<usize>,
}

impl<R: Read> Read for Counting<'_, R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let read = self.inner.read(buf)?;
        self.consumed.set(self.consumed.get() + read);
        Ok(read)
    }
}

/// Below this magnitude a value is written in scientific notation.
///
/// β is 2.3e-11 m^-1 s^-1 and r is 1e-7 s^-1: written out in full they are a
/// row of zeros a reader has to count, and miscounting them is the kind of
/// error the panel exists to catch.
const SCIENTIFIC_BELOW: f64 = 1e-3;

/// A physical quantity, in the notation that keeps its magnitude readable.
fn format_quantity(value: f64) -> String {
    if value != 0.0 && value.abs() < SCIENTIFIC_BELOW {
        format!("{value:e}")
    } else {
        format!("{value}")
    }
}

/// A west-to-east longitude span, as degrees east or west of the meridian.
///
/// The basin crosses the antimeridian, so the eastern boundary is numerically
/// the smaller of the two; the pair is written in the order it is traversed
/// rather than sorted (`termocline_format::BasinExtent`).
fn longitude_span(west_deg_east: f64, east_deg_east: f64) -> String {
    format!(
        "{} – {}",
        Hemispheric::longitude(west_deg_east),
        Hemispheric::longitude(east_deg_east)
    )
}

/// A south-to-north latitude span.
fn latitude_span(south_deg_north: f64, north_deg_north: f64) -> String {
    format!(
        "{} – {}",
        Hemispheric::latitude(south_deg_north),
        Hemispheric::latitude(north_deg_north)
    )
}

/// A signed degree value written the way a chart labels it: a magnitude and a
/// hemisphere, never a minus sign.
struct Hemispheric {
    /// Degrees from the meridian or the equator, always non-negative.
    magnitude_deg: f64,
    /// `E`/`W` for a longitude, `N`/`S` for a latitude.
    hemisphere: char,
}

impl Hemispheric {
    /// `deg_east` written as degrees east or west.
    fn longitude(deg_east: f64) -> Self {
        Self {
            magnitude_deg: deg_east.abs(),
            hemisphere: if deg_east < 0.0 { 'W' } else { 'E' },
        }
    }

    /// `deg_north` written as degrees north or south.
    fn latitude(deg_north: f64) -> Self {
        Self {
            magnitude_deg: deg_north.abs(),
            hemisphere: if deg_north < 0.0 { 'S' } else { 'N' },
        }
    }
}

impl fmt::Display for Hemispheric {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:.1}°{}", self.magnitude_deg, self.hemisphere)
    }
}
