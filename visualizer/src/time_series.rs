//! The point time series: `h(t)` at one cell of the basin, across the whole
//! run.
//!
//! The map draws every place at one instant; this draws one place at every
//! instant. It is the "Niño-index" reading of a run — pick a spot, watch the
//! thermocline there rise and fall — and it is the view in which a wave is
//! seen to *arrive*: a westerly wind burst launches an eastward Kelvin wave
//! (`CONTEXT.md`), and a point near the eastern boundary deepens a basin
//! crossing later than the perturbation that caused it.
//!
//! Like [`crate::heatmap`], [`crate::wind`] and [`crate::cross_section`], none
//! of this knows what a GPU is: a [`PointSeries`] is a list of samples with a
//! time, an anomaly and a place on a unit rectangle, so the acceptance
//! criterion is asserted in `tests/time_series.rs` rather than looked at, and
//! the same code draws the chart in a browser and natively ([ADR-0006]).
//!
//! # What this view costs, and when
//!
//! Every other view of a run reads **one** frame: the map, the overlay and the
//! cross-section are all of the frame the scrubber names, and
//! [`crate::LoadedRun::frame`] reaches it without touching any other. A time
//! series is the transpose of that — one cell of *every* frame — so the
//! indexed lookup buys it nothing, and a series rebuilt per repaint would walk
//! all 731 frames of the scenario run sixty times a second.
//!
//! So the series is built **once per point the reader picks**, and nothing
//! else rebuilds it:
//!
//! - Scrubbing does not. The series is of the whole run; the chosen frame only
//!   moves the marker drawn on it, which is arithmetic on a sample already in
//!   hand.
//! - Playing does not, for the same reason — and playback is the path that
//!   would have hurt most, since it changes the frame every repaint.
//! - Toggling the chart does not: the series is held by the panel, not by the
//!   chart.
//! - Re-clicking the cell already selected does not, because it names the same
//!   cell ([`BasinPoint`] compares by index).
//!
//! What is left is one walk of the run per *click*, which is the work the
//! click is for. `crate::app`'s tests pin that by name.
//!
//! # The vertical axis is this series', not the run's
//!
//! The map, the colour bar and the cross-section all share the run-wide
//! [`DivergingScale`], so that one colour or one height means one anomaly
//! everywhere. This chart deliberately does not, and the label says so.
//!
//! The run-wide half-range is set by the largest anomaly *anywhere* in the
//! run, which for a wind-burst run is the piled-up western wall. Drawn on that
//! axis, an eastern point's arrival — the thing this view exists to show —
//! would be a wiggle a few pixels tall. A time series is read for its shape in
//! time, not for its height against another point's, so the axis is this
//! series' own range and the reader is told what it reaches.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use termocline_format::{FormatError, Frame, GridSpec, Variable};

use crate::chart::{axis_fraction, longitude_at};
use crate::{DivergingScale, LoadedRun};

/// One cell of the basin, as a reader picks it off the map.
///
/// It carries the shape of the field it was picked out of as well as the two
/// indices, so that a point picked off one run's map can be *refused* by
/// another rather than silently reading a different place
/// ([`PointSeries::at_point`]).
#[derive(Debug, Clone, Copy)]
pub struct BasinPoint {
    /// Column into the thermocline-depth field, counting east from the western
    /// wall.
    column: usize,
    /// Row into the thermocline-depth field, counting **north** from the
    /// southern wall — the field's own order, not the map's, which is drawn
    /// from the north down (`crate::heatmap`).
    row: usize,
    /// Longitude of this cell's centre, in degrees east of the prime meridian,
    /// folded into `[-180, 180)` as [`BasinExtent`] states its bounds.
    longitude_deg_east: f64,
    /// Latitude of this cell's centre, in degrees north of the equator.
    latitude_deg_north: f64,
    /// Columns and rows of the cell-centred field this cell is one of.
    field_shape: (usize, usize),
}

/// By index and field shape alone: the longitude and the latitude are a pure
/// function of those, so comparing them would only add two float comparisons
/// to the answer the indices already gave.
///
/// Comparing by index is what stops a drag across one cell of the map from
/// walking the run once per mouse move: two clicks a pixel apart inside one
/// cell name the same cell.
impl PartialEq for BasinPoint {
    fn eq(&self, other: &Self) -> bool {
        self.column == other.column
            && self.row == other.row
            && self.field_shape == other.field_shape
    }
}

impl Eq for BasinPoint {}

impl BasinPoint {
    /// The cell under a click `east` of the way across the basin map and
    /// `down` it, both as fractions of the map's own rectangle.
    ///
    /// `down` rather than north-up because that is the map's order: row 0 of
    /// the image is the northernmost (`crate::heatmap`), and this is the one
    /// place that flip is undone.
    ///
    /// `None` for a click outside the map, or for a fraction that is not a
    /// number: a click off the basin selects nothing rather than the nearest
    /// edge cell, because a reader who missed the map did not mean the coast.
    #[must_use]
    pub fn at_map_fraction(grid: GridSpec, east: f64, down: f64) -> Option<Self> {
        // The field's own shape, from the staggering `h` declares, rather than
        // the grid's cell counts: the two agree for a cell-centred field, and
        // taking them from one place is what keeps them agreeing.
        let (width, height) = grid
            .grid()
            .field_shape(Variable::ThermoclineDepthAnomaly.staggering());
        let column = cell_index(east, width)?;
        let row_from_north = cell_index(down, height)?;
        let row = height - 1 - row_from_north;
        let extent = grid.extent();
        #[allow(clippy::cast_precision_loss)]
        let x_fraction = (column as f64 + 0.5) / width as f64;
        #[allow(clippy::cast_precision_loss)]
        let y_fraction = (row as f64 + 0.5) / height as f64;
        Some(Self {
            column,
            row,
            longitude_deg_east: longitude_at(extent, x_fraction),
            latitude_deg_north: (extent.north_deg_north - extent.south_deg_north)
                .mul_add(y_fraction, extent.south_deg_north),
            field_shape: (width, height),
        })
    }

    /// Column into the field, counting east from the western wall.
    #[must_use]
    pub const fn column(&self) -> usize {
        self.column
    }

    /// Row into the field, counting north from the southern wall.
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }

    /// Longitude of this cell's centre, in degrees east of the prime meridian:
    /// negative west of it, as the basin's bounds are written.
    #[must_use]
    pub const fn longitude_deg_east(&self) -> f64 {
        self.longitude_deg_east
    }

    /// Latitude of this cell's centre, in degrees north of the equator.
    #[must_use]
    pub const fn latitude_deg_north(&self) -> f64 {
        self.latitude_deg_north
    }

    /// Columns and rows of the cell-centred field this cell is one of: the
    /// shape of the basin it was picked off, so a panel can place it back on
    /// the map it came from.
    #[must_use]
    pub const fn field_shape(&self) -> (usize, usize) {
        self.field_shape
    }

    /// Where this cell's value sits in its field.
    const fn offset(&self) -> usize {
        self.row * self.field_shape.0 + self.column
    }
}

/// One sample of a series: what the run held at the chosen cell in one frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeriesSample {
    /// The frame's model time, in seconds since the start of the run.
    t_s: f64,
    /// Thermocline depth anomaly `h` at the chosen cell, in metres. Positive is
    /// deeper than the mean depth `H` (`CONTEXT.md`).
    h_m: f64,
    /// Mixed-layer SST anomaly `T'` there, in kelvin, or `None` for a frame
    /// that carries none.
    ///
    /// `Option` and not zero: an uncoupled run has no SST to report, and zero
    /// kelvin of anomaly is a claim about the ocean that such a run never made
    /// (`termocline_format::Frame`).
    sst_anomaly_k: Option<f64>,
}

impl SeriesSample {
    /// The frame's model time, in seconds since the start of the run.
    #[must_use]
    pub const fn t_s(&self) -> f64 {
        self.t_s
    }

    /// Thermocline depth anomaly `h` here, in metres.
    #[must_use]
    pub const fn h_m(&self) -> f64 {
        self.h_m
    }

    /// Mixed-layer SST anomaly `T'` here, in kelvin, or `None` if this frame
    /// carries none.
    #[must_use]
    pub const fn sst_anomaly_k(&self) -> Option<f64> {
        self.sst_anomaly_k
    }
}

/// How far an SST axis reaches either side of zero, in kelvin.
///
/// The kelvin twin of [`DivergingScale`], which is stated in metres. `T'` is an
/// anomaly like `h` is (`CONTEXT.md`, *SST anomaly*), so the axis is symmetric
/// about zero for the same reason: the reading that matters is which side of
/// zero it is on and by how much.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SstScale {
    /// How far the axis reaches either side of zero, in kelvin. Never
    /// negative; zero for a series that is everywhere zero.
    half_range_k: f64,
}

impl SstScale {
    /// The scale that just covers `values_k`, symmetrically about zero.
    ///
    /// Non-finite values are ignored rather than propagated: one `NaN` would
    /// otherwise leave the whole series with no axis at all.
    #[must_use]
    fn symmetric_over(values_k: impl Iterator<Item = f64>) -> Self {
        let half_range_k = values_k
            .filter(|value| value.is_finite())
            .fold(0.0_f64, |largest, value| largest.max(value.abs()));
        Self { half_range_k }
    }

    /// How far the axis reaches either side of zero, in kelvin.
    #[must_use]
    pub const fn half_range_k(&self) -> f64 {
        self.half_range_k
    }
}

/// One cell's thermocline depth anomaly through the whole run.
#[derive(Debug, Clone)]
pub struct PointSeries {
    /// The cell this is a series of.
    point: BasinPoint,
    /// One sample per frame of the run, in the order the run wrote them.
    samples: Vec<SeriesSample>,
    /// The scale the `h` axis is drawn on: this series', not the run's (see
    /// the module docs).
    scale: DivergingScale,
    /// The scale the `T'` axis is drawn on, or `None` for a run that carries no
    /// SST anomaly at this cell.
    sst_scale: Option<SstScale>,
}

impl PointSeries {
    /// The series at `point` through every frame of `run`.
    ///
    /// This is the one place in the visualizer that walks a whole run, and it
    /// is why the caller must hold on to the result: see the module docs for
    /// what must not rebuild it.
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] if `point` is not a cell of `run`'s basin —
    /// a point picked off the map of one run and asked of another — or if a
    /// frame of `run` carries a field too short to hold it.
    pub fn at_point(run: &LoadedRun, point: BasinPoint) -> Result<Self, FormatError> {
        let grid = run.header().grid;
        let shape = grid
            .grid()
            .field_shape(Variable::ThermoclineDepthAnomaly.staggering());
        // The point's own basin against this run's, before a single frame is
        // read: two basins of different shapes have a cell at any given index
        // pair, and reading the wrong one would be a plausible-looking series
        // of somewhere else entirely.
        if point.field_shape() != shape {
            return Err(FormatError::FieldShape {
                variable: Variable::ThermoclineDepthAnomaly,
                expected: grid.field_len(Variable::ThermoclineDepthAnomaly),
                actual: point.field_shape().0 * point.field_shape().1,
            });
        }
        let offset = point.offset();
        let mut samples = Vec::with_capacity(usize::try_from(run.frame_count()).unwrap_or(0));
        for index in 0..run.frame_count() {
            let frame = run
                .frame(index)
                .expect("a run reports only the frames it holds");
            samples.push(sample_at(&frame, offset)?);
        }
        let scale = DivergingScale::symmetric_over(
            &samples.iter().map(SeriesSample::h_m).collect::<Vec<f64>>(),
        );
        // `None` unless *some* frame reported an SST anomaly here: a run that
        // never coupled SST gets no second axis at all, rather than an axis
        // drawn over values it did not produce.
        let sst_scale = samples
            .iter()
            .any(|sample| sample.sst_anomaly_k.is_some())
            .then(|| {
                SstScale::symmetric_over(samples.iter().filter_map(SeriesSample::sst_anomaly_k))
            });
        Ok(Self {
            point,
            samples,
            scale,
            sst_scale,
        })
    }

    /// The cell this is a series of.
    #[must_use]
    pub const fn point(&self) -> BasinPoint {
        self.point
    }

    /// The samples, in the order the run wrote its frames.
    #[must_use]
    pub fn samples(&self) -> &[SeriesSample] {
        &self.samples
    }

    /// The scale the `h` axis is drawn on: this series' own range, for the
    /// reason in the module docs.
    #[must_use]
    pub const fn scale(&self) -> DivergingScale {
        self.scale
    }

    /// The scale the `T'` axis is drawn on, or `None` for a run that carries no
    /// SST anomaly here.
    ///
    /// `None` is the honest answer for an uncoupled run: there is no SST to
    /// draw, and a flat line at zero would say the ocean was at its
    /// climatological temperature throughout, which the run never claimed.
    #[must_use]
    pub const fn sst_scale(&self) -> Option<SstScale> {
        self.sst_scale
    }

    /// Whether this run reports an SST anomaly at this cell at all.
    #[must_use]
    pub const fn carries_sst_anomaly(&self) -> bool {
        self.sst_scale.is_some()
    }

    /// Model time from the first sample to the last, in seconds. Zero for a
    /// series of one sample, which spans no time at all.
    #[must_use]
    pub fn span_s(&self) -> f64 {
        match (self.samples.first(), self.samples.last()) {
            (Some(first), Some(last)) => last.t_s - first.t_s,
            _ => 0.0,
        }
    }

    /// Where `sample`'s thermocline depth anomaly goes on a unit rectangle:
    /// `(east, down)` from its north-west corner, or `None` for an anomaly that
    /// is not a number.
    ///
    /// The same convention [`crate::CrossSection::plot_position`] uses — `y`
    /// down, because that is how a panel is laid out — so a
    /// deeper-than-average anomaly sits above the middle of the chart, where
    /// the map beside it draws it warm ([`crate::chart::axis_fraction`]).
    #[must_use]
    pub fn plot_position(&self, sample: &SeriesSample) -> Option<(f64, f64)> {
        Some((
            self.time_fraction(sample.t_s),
            axis_fraction(sample.h_m, self.scale.half_range_m())?,
        ))
    }

    /// Where `sample`'s SST anomaly goes on the same unit rectangle, or `None`
    /// if this frame carries none, this run carries none, or the value is not a
    /// number.
    #[must_use]
    pub fn sst_plot_position(&self, sample: &SeriesSample) -> Option<(f64, f64)> {
        let half_range_k = self.sst_scale?.half_range_k();
        Some((
            self.time_fraction(sample.t_s),
            axis_fraction(sample.sst_anomaly_k?, half_range_k)?,
        ))
    }

    /// How far across the chart the instant `t_s` sits, as a fraction of the
    /// run's span.
    ///
    /// Clamped to the chart: an instant outside the run — the frame time of
    /// another run — belongs at the end it ran past, not off the axis. A run
    /// that spans no time at all puts everything in the middle, which is the
    /// whole of the axis it has.
    #[must_use]
    pub fn time_fraction(&self, t_s: f64) -> f64 {
        let span_s = self.span_s();
        if span_s <= 0.0 || !t_s.is_finite() {
            return 0.5;
        }
        let first_s = self.samples.first().map_or(0.0, SeriesSample::t_s);
        ((t_s - first_s) / span_s).clamp(0.0, 1.0)
    }
}

/// The sample `frame` holds at `offset` into its cell-centred fields.
///
/// # Errors
/// [`FormatError::FieldShape`] if a field the frame *does* carry is too short
/// to hold `offset`. A carried-but-short `T'` is an error and not an absence:
/// reporting it as "this run has no SST" would be exactly the lie the
/// [`Option`] exists to avoid, told about a run that does have one.
fn sample_at(frame: &Frame, offset: usize) -> Result<SeriesSample, FormatError> {
    let too_short = |variable, field: &[f64]| FormatError::FieldShape {
        variable,
        expected: offset + 1,
        actual: field.len(),
    };
    let h = frame.h();
    let h_m = *h
        .get(offset)
        .ok_or_else(|| too_short(Variable::ThermoclineDepthAnomaly, h))?;
    // The SST anomaly is cell-centred like `h`
    // (`termocline_format::Variable::staggering`), so one offset serves both.
    // A frame that carries none stays `None` rather than becoming a zero.
    let sst_anomaly_k = match frame.sst_anomaly_k() {
        None => None,
        Some(field) => Some(
            *field
                .get(offset)
                .ok_or_else(|| too_short(Variable::SstAnomaly, field))?,
        ),
    };
    Ok(SeriesSample {
        t_s: frame.t_s(),
        h_m,
        sst_anomaly_k,
    })
}

/// The index of the cell `fraction` of the way along an axis of `cells` cells,
/// or `None` for a fraction outside the axis or one that is not a number.
fn cell_index(fraction: f64, cells: usize) -> Option<usize> {
    if !fraction.is_finite() || !(0.0..1.0).contains(&fraction) || cells == 0 {
        return None;
    }
    #[allow(clippy::cast_precision_loss)]
    let scaled = fraction * cells as f64;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    // `floor` then a bound: a fraction a hair under one can still round up to
    // `cells` in floating point, and that would index past the field.
    let index = (scaled.floor() as usize).min(cells - 1);
    Some(index)
}
