//! The arithmetic every chart of the basin shares.
//!
//! Two pieces, both of them things more than one view needs and neither of
//! them anything a view should re-derive: where a longitude falls across a
//! basin that crosses the antimeridian, and where a value falls on an axis
//! drawn symmetrically about zero.
//!
//! Like the views that use it, none of this knows what a GPU is.

use termocline_format::BasinExtent;

/// A full turn of longitude, in degrees.
///
/// The basin crosses the antimeridian (`CONTEXT.md`, *Basin*), so a longitude
/// accumulated eastward from the western wall has to be folded back into the
/// degrees-east convention [`BasinExtent`] states its bounds in. Named for the
/// same reason `engine/src/basin.rs` names it: it is the modulus a zonal span
/// is measured in, not a magic number.
const FULL_TURN_DEG: f64 = 360.0;

/// Half a turn of longitude, in degrees: the fold point of the degrees-east
/// convention, which runs from 180°W to 180°E.
const HALF_TURN_DEG: f64 = FULL_TURN_DEG / 2.0;

/// The longitude `fraction` of the way east across `extent`, in degrees east of
/// the prime meridian.
///
/// The span is taken eastward around the globe — the basin crosses the
/// antimeridian, so the eastern bound is numerically the smaller of the two —
/// and the result is folded back into `[-180, 180)`.
pub(crate) fn longitude_at(extent: BasinExtent, fraction: f64) -> f64 {
    let span_deg = (extent.east_deg_east - extent.west_deg_east).rem_euclid(FULL_TURN_DEG);
    let absolute_deg = span_deg.mul_add(fraction, extent.west_deg_east);
    (absolute_deg + HALF_TURN_DEG).rem_euclid(FULL_TURN_DEG) - HALF_TURN_DEG
}

/// How far *down* a unit rectangle `value` sits on an axis reaching
/// `half_range` either side of zero, or `None` for a value that is not a
/// number.
///
/// `value` and `half_range` are in the axis's own unit — metres of thermocline
/// depth anomaly for one line, kelvin of SST anomaly for another — and the
/// answer is dimensionless, which is exactly why one function serves both.
///
/// `y` is down, because that is how a panel is laid out, so a positive value
/// sits *above* the middle: a deeper-than-average thermocline is drawn high,
/// where the basin map beside it draws it warm.
///
/// A value past the ends of the axis is clamped to them. That cannot happen
/// for an axis built over the data being drawn; it can for a caller that mixed
/// two runs, and a point on the chart is a better answer than one drawn off
/// it. A non-finite value gets no position at all: the integration diverged
/// there, and a line drawn through the gap would claim a value the run never
/// produced.
pub(crate) fn axis_fraction(value: f64, half_range: f64) -> Option<f64> {
    if !value.is_finite() {
        return None;
    }
    // A field that is everywhere zero has no range to normalize by, and every
    // value of it is exactly on the zero line anyway.
    let above_zero = if half_range == 0.0 {
        0.0
    } else {
        (value / half_range).clamp(-1.0, 1.0)
    };
    Some(0.5 - above_zero / 2.0)
}
