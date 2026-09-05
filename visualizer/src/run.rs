//! A run, once its bytes have arrived, and the metadata the shell shows for it.
//!
//! The shell's job in T-08.1 is to load a run and say what it is. Everything
//! about that is here, and none of it knows what a window or a GPU is: a run
//! is two byte buffers in, a labelled list of rows out. Where the bytes came
//! from — a directory, a pair of dropped files, an HTTP fetch — is the
//! caller's problem, which is what lets the same code serve the browser and
//! the native build (ADR-0006).

use std::fmt;

use termocline_format::{RunHeader, RunReadError, RunReader};

/// Seconds in a day, for reporting a run's model time in the unit its output
/// cadence is chosen in (`steady-trades.toml` writes a frame a day).
const SECONDS_PER_DAY: f64 = 86_400.0;

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
/// re-read — the bytes a user dropped or a fetch returned are the only copy —
/// and Epic 09 draws them. T-08.1 only counts them.
#[derive(Debug, Clone)]
pub struct LoadedRun {
    /// Where the run came from, for a reader checking they opened the one they
    /// meant to.
    source: String,
    /// Everything the frames do not say about themselves.
    header: RunHeader,
    /// The run's encoded frames, undecoded.
    frames: Vec<u8>,
}

impl LoadedRun {
    /// Decode `bytes`, which came from `source`.
    ///
    /// The header goes through [`RunReader`] rather than `serde_json` directly,
    /// so a run written by a format version this build does not read is
    /// refused here rather than mislabelled in the UI.
    ///
    /// # Errors
    /// The errors of [`RunReader::new`]: a header that is not valid JSON for a
    /// [`RunHeader`], or one from an unsupported [format version].
    ///
    /// [format version]: termocline_format::FORMAT_VERSION
    pub fn from_bytes(source: impl Into<String>, bytes: &RunBytes) -> Result<Self, RunReadError> {
        let reader = RunReader::new(bytes.header.as_slice(), bytes.frames.as_slice())?;
        let header = reader.header().clone();
        Ok(Self {
            source: source.into(),
            header,
            frames: bytes.frames.clone(),
        })
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

    /// The run's encoded frames, for the renderer of Epic 09.
    #[must_use]
    pub fn frame_bytes(&self) -> &[u8] {
        &self.frames
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
                format!("{} m", quantity(params.mean_depth_m)),
            ),
            MetadataRow::new(
                "Reduced gravity g'",
                format!("{} m s^-2", quantity(params.reduced_gravity_m_per_s2)),
            ),
            MetadataRow::new(
                "Coriolis gradient β",
                format!("{} m^-1 s^-1", quantity(params.beta_per_m_per_s)),
            ),
            MetadataRow::new(
                "Rayleigh damping r",
                format!("{} s^-1", quantity(params.rayleigh_damping_per_s)),
            ),
            MetadataRow::new(
                "Reference density ρ₀",
                format!("{} kg m^-3", quantity(params.reference_density_kg_per_m3)),
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

/// Below this magnitude a value is written in scientific notation.
///
/// β is 2.3e-11 m^-1 s^-1 and r is 1e-7 s^-1: written out in full they are a
/// row of zeros a reader has to count, and miscounting them is the kind of
/// error the panel exists to catch.
const SCIENTIFIC_BELOW: f64 = 1e-3;

/// A physical quantity, in the notation that keeps its magnitude readable.
fn quantity(value: f64) -> String {
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
