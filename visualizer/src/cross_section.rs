//! The equatorial cross-section: one frame's `h` along the equator, as a line
//! against longitude.
//!
//! The basin map draws the whole field and leaves the reader to judge a slope
//! from two colours; this draws the slope itself. It is the thermocline tilt
//! (`CONTEXT.md`) as a number one can read off an axis — deep in the west,
//! shallow in the east, and the single zero crossing between them — and it is
//! where a Kelvin pulse is seen to travel rather than inferred from a
//! reddening patch.
//!
//! Like [`crate::heatmap`] and [`crate::wind`], none of this knows what a GPU
//! is: a [`CrossSection`] is a list of points with a longitude, an anomaly and
//! a place on a unit rectangle, so the acceptance criterion — that the line of
//! a known equilibrium frame matches T-07.4's analytic tilt — is asserted in
//! `tests/cross_section.rs` rather than looked at, and the same code draws the
//! chart in a browser and natively ([ADR-0006]).
//!
//! # Why the vertical axis is the run's
//!
//! The scale is [`crate::DivergingScale`], the run's and not the frame's, for
//! the reason the heatmap's colours are: a chart rescaled frame by frame draws
//! every frame at full height, and a tilt that collapses over a run — which is
//! El Niño (`CONTEXT.md`, *ENSO*) — would look like a tilt that never moved.
//!
//! # Where "the equator" is
//!
//! `h` lives at cell centres (ADR-0003), and a basin laid out symmetrically
//! about the equator in an even number of rows has no cell centre on it: the
//! scenario basin's two nearest sit at 0.25°S and 0.25°N. Taking either alone
//! would put the section off the axis of the waveguide the model is about, so
//! the rows nearest the equator are averaged — one of them where a row sits on
//! the equator, two where it falls between them.
//! [`CrossSection::latitude_deg_north`] and [`CrossSection::rows_averaged`] say
//! which happened, so the chart can label itself with the latitude it actually
//! read.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use termocline_format::{FormatError, Frame, GridSpec, Variable};

use crate::DivergingScale;

/// A full turn of longitude, in degrees.
///
/// The basin crosses the antimeridian (`CONTEXT.md`, *Basin*), so a longitude
/// accumulated eastward from the western wall has to be folded back into the
/// degrees-east convention [`termocline_format::BasinExtent`] states its bounds
/// in. Named for the same reason `engine/src/basin.rs` names it: it is the
/// modulus a zonal span is measured in, not a magic number.
const FULL_TURN_DEG: f64 = 360.0;

/// Half a turn of longitude, in degrees: the fold point of the degrees-east
/// convention, which runs from 180°W to 180°E.
const HALF_TURN_DEG: f64 = FULL_TURN_DEG / 2.0;

/// How far two rows' distances from the equator may differ and still count as
/// equal, as a fraction of the row spacing.
///
/// A basin symmetric about the equator puts its two innermost rows the same
/// distance either side of it, and that distance is computed from bounds
/// written in decimal degrees — which are not exact in binary, so the two
/// distances can land a few ulp apart. This bound is many orders of magnitude
/// looser than that slack and many orders tighter than the half-row offset
/// that would mean the basin genuinely is not symmetric.
const SAME_DISTANCE_FRACTION: f64 = 1e-9;

/// One point of the section: the anomaly at one cell of the equatorial row,
/// and where on the chart it goes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrossSectionPoint {
    /// How far east of the basin's western wall this point sits, as a fraction
    /// of the basin's width. Strictly between 0 and 1: the point is a cell
    /// centre, not a wall.
    x_fraction: f64,
    /// The longitude of that cell centre, in degrees east of the prime
    /// meridian, folded into `[-180, 180)` as
    /// [`termocline_format::BasinExtent`] states its bounds.
    longitude_deg_east: f64,
    /// The thermocline depth anomaly here, in metres. Positive is deeper than
    /// the mean depth `H` (`CONTEXT.md`).
    h_m: f64,
}

impl CrossSectionPoint {
    /// How far east of the basin's western wall this point sits, as a fraction
    /// of the basin's width.
    ///
    /// The axis position rather than the longitude, because the longitude
    /// wraps through the antimeridian and an axis that wrapped with it would
    /// draw the basin folded over itself.
    #[must_use]
    pub const fn x_fraction(&self) -> f64 {
        self.x_fraction
    }

    /// The longitude of this cell centre, in degrees east of the prime
    /// meridian: negative west of it, as the basin's bounds are written.
    #[must_use]
    pub const fn longitude_deg_east(&self) -> f64 {
        self.longitude_deg_east
    }

    /// The thermocline depth anomaly here, in metres. Positive is a
    /// deeper-than-average thermocline (`CONTEXT.md`).
    #[must_use]
    pub const fn h_m(&self) -> f64 {
        self.h_m
    }
}

/// One frame's thermocline depth anomaly along the equator.
#[derive(Debug, Clone)]
pub struct CrossSection {
    /// One point per cell of the basin's zonal axis, west to east.
    points: Vec<CrossSectionPoint>,
    /// The latitude the section was actually read at, in degrees north.
    latitude_deg_north: f64,
    /// How many rows of the field were averaged to get it: one where a row
    /// sits on the equator, two where it falls between them.
    rows_averaged: usize,
    /// The scale the vertical axis is drawn on, shared with the basin map.
    scale: DivergingScale,
}

impl CrossSection {
    /// The equatorial section of `frame`'s thermocline depth anomaly over
    /// `grid`, on `scale`.
    ///
    /// The scale is the caller's because it belongs to the run rather than to
    /// this frame: [`crate::LoadedRun::anomaly_scale`] covers every frame, so
    /// the same height means the same anomaly in every frame of the run, and
    /// the chart and the map under it agree.
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
        let (width, _height) = grid
            .grid()
            .field_shape(Variable::ThermoclineDepthAnomaly.staggering());
        let extent = grid.extent();
        let rows = equatorial_rows(grid);
        let h_m = frame.h();
        #[allow(clippy::cast_precision_loss)]
        let rows_per_point = rows.len() as f64;
        let mut points = Vec::with_capacity(width);
        for i in 0..width {
            let sum_m: f64 = rows.iter().map(|&j| h_m[j * width + i]).sum();
            #[allow(clippy::cast_precision_loss)]
            let x_fraction = (i as f64 + 0.5) / width as f64;
            points.push(CrossSectionPoint {
                x_fraction,
                longitude_deg_east: longitude_at(
                    extent.west_deg_east,
                    extent.east_deg_east,
                    x_fraction,
                ),
                h_m: sum_m / rows_per_point,
            });
        }
        Ok(Self {
            latitude_deg_north: mean_latitude_deg_north(grid, &rows),
            rows_averaged: rows.len(),
            points,
            scale,
        })
    }

    /// The points of the line, west to east.
    #[must_use]
    pub fn points(&self) -> &[CrossSectionPoint] {
        &self.points
    }

    /// The latitude the section was read at, in degrees north. Zero for a
    /// basin laid out symmetrically about the equator, which every scenario's
    /// is (`CONTEXT.md`, *Basin*).
    #[must_use]
    pub const fn latitude_deg_north(&self) -> f64 {
        self.latitude_deg_north
    }

    /// How many rows of the field were averaged: one where a cell centre sits
    /// on the equator, two where the equator falls between two of them.
    #[must_use]
    pub const fn rows_averaged(&self) -> usize {
        self.rows_averaged
    }

    /// The scale the vertical axis is drawn on: the run's, so the chart and
    /// the colour bar beside it say the same thing.
    #[must_use]
    pub const fn scale(&self) -> DivergingScale {
        self.scale
    }

    /// Where `point` goes on a unit rectangle: `(east, down)` from its
    /// north-west corner, or `None` for an anomaly that is not a number.
    ///
    /// The same convention [`crate::WindArrow`] places arrows in — `y` down,
    /// because that is how a panel is laid out — so a deeper-than-average
    /// anomaly sits *above* the middle of the chart, where the map beside it
    /// draws it warm.
    ///
    /// An anomaly past the ends of the scale is clamped to them. That cannot
    /// happen for a scale built over the run being drawn; it can for a caller
    /// that mixed two runs, and a point on the frame is a better answer than
    /// one drawn off it. A non-finite anomaly gets no position at all: the
    /// integration diverged there, and a line drawn through the gap would
    /// claim a value the run never produced.
    #[must_use]
    pub fn plot_position(&self, point: &CrossSectionPoint) -> Option<(f64, f64)> {
        if !point.h_m.is_finite() {
            return None;
        }
        let half_range_m = self.scale.half_range_m();
        // A run that is everywhere zero has no range to normalize by, and
        // every point of it is exactly on the zero line anyway.
        let above_zero = if half_range_m == 0.0 {
            0.0
        } else {
            (point.h_m / half_range_m).clamp(-1.0, 1.0)
        };
        Some((point.x_fraction, 0.5 - above_zero / 2.0))
    }
}

/// The rows of a cell-centred field nearest the equator: one, or the two that
/// straddle it.
fn equatorial_rows(grid: GridSpec) -> Vec<usize> {
    let mut nearest = Vec::new();
    let mut smallest_deg = f64::INFINITY;
    let same_deg = row_spacing_deg(grid).abs() * SAME_DISTANCE_FRACTION;
    for j in 0..grid.ny() {
        let distance_deg = latitude_of_row_deg_north(grid, j).abs();
        if distance_deg < smallest_deg - same_deg {
            smallest_deg = distance_deg;
            nearest.clear();
        }
        if distance_deg <= smallest_deg + same_deg {
            smallest_deg = smallest_deg.min(distance_deg);
            nearest.push(j);
        }
    }
    nearest
}

/// The mean latitude of `rows`, in degrees north: the latitude the section was
/// read at.
fn mean_latitude_deg_north(grid: GridSpec, rows: &[usize]) -> f64 {
    let sum_deg: f64 = rows
        .iter()
        .map(|&j| latitude_of_row_deg_north(grid, j))
        .sum();
    #[allow(clippy::cast_precision_loss)]
    let count = rows.len() as f64;
    sum_deg / count
}

/// The latitude of the centres of row `j`, in degrees north.
fn latitude_of_row_deg_north(grid: GridSpec, j: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let offset = j as f64 + 0.5;
    offset.mul_add(row_spacing_deg(grid), grid.extent().south_deg_north)
}

/// The meridional size of one cell, in degrees.
fn row_spacing_deg(grid: GridSpec) -> f64 {
    let extent = grid.extent();
    #[allow(clippy::cast_precision_loss)]
    let rows = grid.ny() as f64;
    (extent.north_deg_north - extent.south_deg_north) / rows
}

/// The longitude `fraction` of the way east across a basin running from
/// `west_deg_east` to `east_deg_east`, in degrees east of the prime meridian.
///
/// The span is taken eastward around the globe — the basin crosses the
/// antimeridian, so the eastern bound is numerically the smaller of the two —
/// and the result is folded back into `[-180, 180)`.
fn longitude_at(west_deg_east: f64, east_deg_east: f64, fraction: f64) -> f64 {
    let span_deg = (east_deg_east - west_deg_east).rem_euclid(FULL_TURN_DEG);
    let absolute_deg = span_deg.mul_add(fraction, west_deg_east);
    (absolute_deg + HALF_TURN_DEG).rem_euclid(FULL_TURN_DEG) - HALF_TURN_DEG
}
