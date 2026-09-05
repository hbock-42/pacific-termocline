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

use termocline_format::{BasinExtent, FormatError, Frame, GridSpec, Variable};

use crate::chart::{axis_fraction, longitude_at};
use crate::DivergingScale;

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
        // Both extents come from the staggering `h` declares, rather than from
        // the cell counts: the two agree for a cell-centred field, and taking
        // them from one place is what keeps them agreeing.
        let (width, height) = grid
            .grid()
            .field_shape(Variable::ThermoclineDepthAnomaly.staggering());
        let axis = MeridionalAxis::of(height, grid.extent());
        let rows = axis.rows_nearest_the_equator();
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
                longitude_deg_east: longitude_at(grid.extent(), x_fraction),
                h_m: sum_m / rows_per_point,
            });
        }
        Ok(Self {
            latitude_deg_north: axis.mean_latitude_deg_north(&rows),
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
    /// An anomaly past the ends of the scale is clamped to them, and a
    /// non-finite one gets no position at all — see
    /// [`crate::chart::axis_fraction`], which the point time series places its
    /// samples with too.
    #[must_use]
    pub fn plot_position(&self, point: &CrossSectionPoint) -> Option<(f64, f64)> {
        Some((
            point.x_fraction,
            axis_fraction(point.h_m, self.scale.half_range_m())?,
        ))
    }
}

/// The basin's meridional axis, as a cell-centred field sees it: how many rows
/// it has, and where each one sits.
///
/// One type rather than four functions each re-deriving the row spacing from a
/// [`GridSpec`]. The rows are the field's own, not the grid's cell count, so
/// nothing here has to assume the two agree.
#[derive(Debug, Clone, Copy)]
struct MeridionalAxis {
    /// Rows of cell centres on the axis.
    rows: usize,
    /// The basin's southern boundary, in degrees north.
    south_deg_north: f64,
    /// The meridional size of one cell, in degrees.
    spacing_deg: f64,
}

impl MeridionalAxis {
    /// The axis of a field `rows` rows tall over `extent`.
    fn of(rows: usize, extent: BasinExtent) -> Self {
        #[allow(clippy::cast_precision_loss)]
        let spacing_deg = (extent.north_deg_north - extent.south_deg_north) / rows as f64;
        Self {
            rows,
            south_deg_north: extent.south_deg_north,
            spacing_deg,
        }
    }

    /// The latitude of the centres of row `j`, in degrees north.
    fn latitude_deg_north(self, j: usize) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let offset = j as f64 + 0.5;
        offset.mul_add(self.spacing_deg, self.south_deg_north)
    }

    /// The rows nearest the equator: one where a row sits on it, or the two
    /// that straddle it.
    fn rows_nearest_the_equator(self) -> Vec<usize> {
        let mut nearest = Vec::new();
        let mut smallest_deg = f64::INFINITY;
        let same_deg = self.spacing_deg.abs() * SAME_DISTANCE_FRACTION;
        for j in 0..self.rows {
            let distance_deg = self.latitude_deg_north(j).abs();
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

    /// The mean latitude of `rows`, in degrees north: the latitude a section
    /// averaged over them was read at.
    fn mean_latitude_deg_north(self, rows: &[usize]) -> f64 {
        let sum_deg: f64 = rows.iter().map(|&j| self.latitude_deg_north(j)).sum();
        #[allow(clippy::cast_precision_loss)]
        let count = rows.len() as f64;
        sum_deg / count
    }
}
