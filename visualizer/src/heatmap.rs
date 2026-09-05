//! The basin map: one frame's thermocline depth anomaly as a colour-mapped
//! image.
//!
//! This is the whole of T-08.2 that has a value in it, and none of it knows
//! what a GPU is — a [`Heatmap`] is a buffer of RGB triples that the shell
//! uploads as a texture. That is what lets the acceptance criterion, a visual
//! smoke test, be asserted in `tests/heatmap.rs` instead of looked at, and it
//! is what lets the same code draw the basin in a browser and natively
//! ([ADR-0006]).
//!
//! # The colour scale
//!
//! `h` is an anomaly, not a total depth (`CONTEXT.md`): it is signed, and the
//! reading that matters is which side of zero a place is on and by how much.
//! So the scale is **diverging and symmetric about zero** — zero is pinned to
//! a neutral colour whatever the frame contains, and the two halves reach
//! equally far — rather than a sequential ramp over the frame's min and max,
//! which would put the neutral colour wherever the data happened to start.
//!
//! The half-range is the run's, not the frame's, so the same colour means the
//! same anomaly in every frame of a run and a decaying anomaly is seen to
//! decay.
//!
//! The colours are ColorBrewer's 11-class `RdBu`, reversed. Three reasons for
//! that scheme in particular:
//!
//! - It is **published and colour-blind safe** (Brewer, Harrower &
//!   Pennsylvania State University, <https://colorbrewer2.org>); the red–blue
//!   axis stays separable under the common dichromacies, which a red–green
//!   diverging scheme does not.
//! - Its middle class is a near-white that reads as *no anomaly* against
//!   either a light or a dark panel, so the zero contour is visible rather
//!   than inferred.
//! - Warm-for-deep matches how the field is read physically: the deep western
//!   thermocline is the warm pool, the shallow eastern one is the cold tongue
//!   (`CONTEXT.md`, *Thermocline tilt*). A reader who knows the ocean and a
//!   reader who only knows the colour bar reach the same conclusion.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use termocline_format::{FormatError, Frame, GridSpec, Variable};

/// ColorBrewer's 11-class `RdBu`, reversed so the ramp runs from the most
/// negative anomaly to the most positive: blue for a shallow thermocline,
/// near-white for none, red for a deep one.
///
/// Transcribed from the published table at <https://colorbrewer2.org>. The
/// classes are anchors of a continuous ramp here rather than eleven discrete
/// bands: `h` is a continuous field, and banding it would draw contours the
/// model never produced.
const RD_BU_RAMP: [[u8; 3]; 11] = [
    [5, 48, 97],
    [33, 102, 172],
    [67, 147, 195],
    [146, 197, 222],
    [209, 229, 240],
    [247, 247, 247],
    [253, 219, 199],
    [244, 165, 130],
    [214, 96, 77],
    [178, 24, 43],
    [103, 0, 31],
];

/// What a point with no usable value is drawn in.
///
/// Off the ramp on purpose. A `NaN` in `h` means the integration diverged, and
/// the one reading that must not be available is "undisturbed ocean" — which
/// is exactly what the neutral middle of the ramp would say.
const NO_VALUE: [u8; 3] = [0, 0, 0];

/// Channels in one colour: red, green and blue.
const CHANNELS: usize = 3;

/// A diverging colour scale for a signed field, pinned neutral at zero.
///
/// The scale is stated in the field's own unit — metres of thermocline depth
/// anomaly — so the shell can label a colour bar with the numbers a reader
/// would quote.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DivergingScale {
    /// How far the scale reaches either side of zero, in metres. Never
    /// negative; zero for a field that is everywhere zero.
    half_range_m: f64,
}

impl DivergingScale {
    /// The scale that covers both this one and `other`.
    ///
    /// How a run-wide scale is built out of its frames: a frame at a time,
    /// while the frames are being walked for other reasons anyway.
    #[must_use]
    pub fn widened(self, other: Self) -> Self {
        Self {
            half_range_m: self.half_range_m.max(other.half_range_m),
        }
    }

    /// The scale that just covers `values_m`, symmetrically about zero.
    ///
    /// The half-range is the largest magnitude among the finite values, so no
    /// value in the field is clamped and zero stays neutral. Non-finite values
    /// are ignored rather than propagated: one `NaN` would otherwise leave the
    /// whole frame with no scale at all.
    #[must_use]
    pub fn symmetric_over(values_m: &[f64]) -> Self {
        let half_range_m = values_m
            .iter()
            .filter(|value| value.is_finite())
            .fold(0.0_f64, |largest, value| largest.max(value.abs()));
        Self { half_range_m }
    }

    /// How far the scale reaches either side of zero, in metres.
    #[must_use]
    pub const fn half_range_m(&self) -> f64 {
        self.half_range_m
    }

    /// The scale itself, as `samples` RGB triples running from its shallow
    /// end to its deep one: the colour bar that says what a map's colours
    /// mean.
    ///
    /// Here rather than in the shell so that every colour the visualizer draws
    /// comes from one place, and so the bar can be checked without a GPU. One
    /// sample is the neutral middle, which is all a bar that narrow could
    /// honestly show.
    #[must_use]
    pub fn bar_rgb(&self, samples: usize) -> Vec<u8> {
        let mut rgb = Vec::with_capacity(samples * CHANNELS);
        for sample in 0..samples {
            #[allow(clippy::cast_precision_loss)]
            let fraction = sample as f64 / (samples.saturating_sub(1).max(1)) as f64;
            rgb.extend_from_slice(&self.color(self.half_range_m * (2.0 * fraction - 1.0)));
        }
        rgb
    }

    /// The colour `value_m` is drawn in.
    ///
    /// Values past the ends of the scale take the end colours. That is
    /// clamping, and it is inherent to a colour bar rather than a silent
    /// substitution: the bar the shell draws states the range, so a saturated
    /// region reads as "at least this much".
    #[must_use]
    pub fn color(&self, value_m: f64) -> [u8; 3] {
        if !value_m.is_finite() {
            return NO_VALUE;
        }
        // A field that is everywhere zero has no range to normalize by, and
        // every point of it is exactly the neutral value anyway.
        if self.half_range_m == 0.0 {
            return ramp_at(0.5);
        }
        ramp_at((value_m / self.half_range_m + 1.0) / 2.0)
    }
}

/// The ramp colour at `position`, where 0 is its most negative end, 0.5 its
/// neutral middle and 1 its most positive end.
fn ramp_at(position: f64) -> [u8; 3] {
    #[allow(clippy::cast_precision_loss)]
    let last = (RD_BU_RAMP.len() - 1) as f64;
    let anchor = position.clamp(0.0, 1.0) * last;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let below = anchor.floor() as usize;
    let above = (below + 1).min(RD_BU_RAMP.len() - 1);
    #[allow(clippy::cast_precision_loss)]
    let fraction = anchor - below as f64;
    let (low, high) = (RD_BU_RAMP[below], RD_BU_RAMP[above]);
    [0, 1, 2].map(|channel| lerp_channel(low[channel], high[channel], fraction))
}

/// One channel of a colour interpolated `fraction` of the way from `low` to
/// `high`, quantized back to eight bits.
fn lerp_channel(low: u8, high: u8, fraction: f64) -> u8 {
    let low = f64::from(low);
    let value = fraction.mul_add(f64::from(high) - low, low);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let rounded = value.round().clamp(0.0, 255.0) as u8;
    rounded
}

/// One frame's thermocline depth anomaly, colour-mapped, ready to upload.
///
/// One pixel per cell of the basin, in reading order: the first row is the
/// northernmost and the first column the westernmost. The field itself is
/// row-major with `j` increasing northward, so the rows are reversed on the
/// way in — a map drawn straight from the buffer would be upside down, and an
/// upside-down basin looks perfectly plausible.
#[derive(Debug, Clone)]
pub struct Heatmap {
    /// Cells along x, west to east.
    width: usize,
    /// Cells along y, north to south as drawn.
    height: usize,
    /// `width * height` RGB triples, row-major from the northwest corner.
    rgb: Vec<u8>,
    /// The scale the colours came from, for the colour bar beside them.
    scale: DivergingScale,
}

impl Heatmap {
    /// The map of `frame`'s thermocline depth anomaly over `grid`, on `scale`.
    ///
    /// The scale is the caller's because it belongs to the run rather than to
    /// this frame: [`crate::LoadedRun::anomaly_scale`] covers every frame, so
    /// a colour means the same anomaly wherever in the run it is seen.
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] if `frame` does not fit `grid` — the frame
    /// of one run against the header of another.
    pub fn of_frame(
        grid: GridSpec,
        frame: &Frame,
        scale: DivergingScale,
    ) -> Result<Self, FormatError> {
        frame.validate(&grid)?;
        let (width, height) = grid
            .grid()
            .field_shape(Variable::ThermoclineDepthAnomaly.staggering());
        let h_m = frame.h();
        let mut rgb = Vec::with_capacity(width * height * CHANNELS);
        for row in 0..height {
            // Row 0 of the image is the northernmost, which is the last row of
            // the field.
            let j = height - 1 - row;
            for i in 0..width {
                rgb.extend_from_slice(&scale.color(h_m[j * width + i]));
            }
        }
        Ok(Self {
            width,
            height,
            rgb,
            scale,
        })
    }

    /// Cells along x, west to east.
    #[must_use]
    pub const fn width(&self) -> usize {
        self.width
    }

    /// Cells along y, north to south as drawn.
    #[must_use]
    pub const fn height(&self) -> usize {
        self.height
    }

    /// The image, as `width * height` RGB triples row-major from the northwest
    /// corner.
    #[must_use]
    pub fn rgb(&self) -> &[u8] {
        &self.rgb
    }

    /// The colour at column `x` and row `y`, or `None` outside the basin.
    #[must_use]
    pub fn pixel(&self, x: usize, y: usize) -> Option<[u8; 3]> {
        if x >= self.width || y >= self.height {
            return None;
        }
        let start = (y * self.width + x) * CHANNELS;
        Some([self.rgb[start], self.rgb[start + 1], self.rgb[start + 2]])
    }

    /// The scale these colours came from.
    #[must_use]
    pub const fn scale(&self) -> &DivergingScale {
        &self.scale
    }
}
