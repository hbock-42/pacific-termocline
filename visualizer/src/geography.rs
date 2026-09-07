//! Where the basin is, in the words an atlas uses: longitudes as degrees east
//! and west, and the two coasts the ocean ends at.
//!
//! Every view in this crate already knows where a point of the basin *is* — a
//! fraction of the way east, and the longitude that fraction lands on
//! ([`crate::chart::longitude_at`]). What none of them had was an axis a
//! reader can read: a picture labelled "cell 160" says nothing to a visitor,
//! and one labelled "180°" says the middle of the Pacific.
//!
//! Everything here is derived from the [`BasinExtent`] the run's header
//! declares, and nothing is hardcoded to the scenario basin. The basin is a
//! scenario parameter (`CONTEXT.md`, *Basin*: roughly 120°E–80°W, "exact
//! truncation is a scenario parameter"), so an axis or a coast name pinned to
//! one particular truncation would be wrong for the next scenario that moves
//! it — and wrong silently, which is the failure mode a labelled picture
//! exists to prevent.
//!
//! Like the views it labels, none of this knows what a GPU is: the ticks and
//! the coasts are values, and `tests/teaching_labels.rs` asserts them
//! ([ADR-0006]).
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use termocline_format::BasinExtent;

/// A full turn of longitude, in degrees — the modulus a zonal span is measured
/// in, as [`crate::chart`] names it for the same reason.
const FULL_TURN_DEG: f64 = 360.0;

/// Half a turn of longitude: the fold point of the degrees-east convention,
/// which runs from 180°W to 180°E.
const HALF_TURN_DEG: f64 = FULL_TURN_DEG / 2.0;

/// The tick spacings an axis of longitude may use, in degrees, coarsest last.
///
/// The round numbers a printed atlas rules its meridians at. The first one
/// that fits the basin in [`MAX_TICKS`] is the one used, so a wide basin gets a
/// coarse axis and a narrow one a fine axis without either being told which.
const TICK_STEPS_DEG: [f64; 5] = [10.0, 20.0, 30.0, 45.0, 60.0];

/// The most ticks an axis of longitude is drawn with.
///
/// The basin is drawn far wider than it is tall and the labels sit under it in
/// a single row: nine of them across the scenario basin's 160° is a label
/// every 20°, which is about as close as two labels can sit before they touch
/// at the width a browser window gives the map.
const MAX_TICKS: usize = 9;

/// Which end of the basin a coast is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Wall {
    /// The western boundary: where the trade winds pile the warm water up.
    West,
    /// The eastern boundary: where the thermocline is shallowest.
    East,
}

/// A stretch of land at the equator, and the longitudes it spans.
///
/// A span rather than a point with a tolerance: a coast is a piece of the
/// world with a width, and asking whether a wall of the basin stands against
/// it is asking whether the wall's longitude is inside that width. There is
/// then no tolerance to justify — only geography.
struct Landfall {
    /// What a reader calls it.
    name: &'static str,
    /// The western end of the land, in degrees east of the prime meridian.
    west_deg_east: f64,
    /// Its eastern end, the same way.
    east_deg_east: f64,
}

/// The land the equatorial Pacific runs between, at the latitudes this model
/// is about.
///
/// Approximate equatorial longitudes of the two continental margins: the
/// maritime continent from Sumatra's west coast (~95°E) to the eastern tip of
/// New Guinea (~151°E), and South America from the coast of Ecuador (~81°W) to
/// the Brazilian coast (~35°W). They bound the basin `CONTEXT.md`, *Basin*
/// describes, and a wall inside one of them is standing on that shore however
/// the scenario truncates the domain.
const LANDFALLS: [Landfall; 2] = [
    Landfall {
        name: "Indonesia / New Guinea",
        west_deg_east: 95.0,
        east_deg_east: 151.0,
    },
    Landfall {
        name: "South America",
        west_deg_east: -81.0,
        east_deg_east: -35.0,
    },
];

/// One end of the basin: which wall it is, where it is, and what land stands
/// there.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coast {
    /// Which end of the basin this is.
    wall: Wall,
    /// The longitude of the wall, in degrees east of the prime meridian,
    /// folded into `[-180, 180)` as [`BasinExtent`] states its bounds.
    longitude_deg_east: f64,
    /// The land the wall stands against, or `None` where the scenario has
    /// truncated the basin in open water.
    land: Option<&'static str>,
}

impl Coast {
    /// Which end of the basin this coast is.
    #[must_use]
    pub const fn wall(&self) -> Wall {
        self.wall
    }

    /// The longitude of the wall, in degrees east of the prime meridian.
    #[must_use]
    pub const fn longitude_deg_east(&self) -> f64 {
        self.longitude_deg_east
    }

    /// The land this wall stands against, or `None` where the basin was
    /// truncated in open water and there is no coast to name.
    #[must_use]
    pub const fn land(&self) -> Option<&'static str> {
        self.land
    }

    /// The wall as a label: the coast's name where it has one, and where it
    /// is, always.
    ///
    /// The longitude is on it either way, so the name and the axis under the
    /// picture cannot come to say different things.
    #[must_use]
    pub fn label(&self) -> String {
        let place = longitude_text(self.longitude_deg_east);
        match self.land {
            Some(name) => format!("{name} ({place})"),
            None => {
                let side = match self.wall {
                    Wall::West => "western",
                    Wall::East => "eastern",
                };
                format!("the {side} edge of the basin ({place})")
            }
        }
    }
}

/// The two ends of `extent`, west first.
#[must_use]
pub fn coasts(extent: BasinExtent) -> [Coast; 2] {
    [
        coast(Wall::West, extent.west_deg_east),
        coast(Wall::East, extent.east_deg_east),
    ]
}

/// The coast at `longitude_deg_east`, on the given end of the basin.
fn coast(wall: Wall, longitude_deg_east: f64) -> Coast {
    let folded = fold(longitude_deg_east);
    Coast {
        wall,
        longitude_deg_east: folded,
        land: LANDFALLS
            .iter()
            .find(|land| (land.west_deg_east..=land.east_deg_east).contains(&folded))
            .map(|land| land.name),
    }
}

/// One labelled meridian of the longitude axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LongitudeTick {
    /// How far east of the basin's western wall this meridian sits, as a
    /// fraction of the basin's width.
    x_fraction: f64,
    /// Its longitude, in degrees east of the prime meridian, folded into
    /// `[-180, 180)`.
    longitude_deg_east: f64,
}

impl LongitudeTick {
    /// How far east of the western wall this meridian sits, as a fraction of
    /// the basin's width.
    ///
    /// The axis position rather than the longitude, for the reason
    /// [`crate::CrossSectionPoint::x_fraction`] gives: the longitude wraps
    /// through the antimeridian and an axis that wrapped with it would draw
    /// the basin folded over itself.
    #[must_use]
    pub const fn x_fraction(&self) -> f64 {
        self.x_fraction
    }

    /// Its longitude, in degrees east of the prime meridian.
    #[must_use]
    pub const fn longitude_deg_east(&self) -> f64 {
        self.longitude_deg_east
    }

    /// The meridian as an atlas writes it: `140°E`, `180°`, `120°W` — a
    /// magnitude and a hemisphere, never a minus sign.
    #[must_use]
    pub fn label(&self) -> String {
        longitude_text(self.longitude_deg_east)
    }
}

/// The meridians of `extent`, west to east.
///
/// Ruled at round longitudes rather than at the basin's own bounds: a reader
/// finds 180° on the axis and knows where the middle of the Pacific is, which
/// a tick at 120°E and another at 80°W does not tell them. Where the basin is
/// too narrow for any round meridian to fall inside it, the two walls are the
/// ticks — an axis with the ends labelled is still an axis.
#[must_use]
pub fn longitude_ticks(extent: BasinExtent) -> Vec<LongitudeTick> {
    let west_deg_east = extent.west_deg_east;
    let span_deg = (extent.east_deg_east - west_deg_east).rem_euclid(FULL_TURN_DEG);
    let step_deg = TICK_STEPS_DEG
        .iter()
        .copied()
        .find(|&step| tick_count(west_deg_east, span_deg, step) <= MAX_TICKS);
    let ticks: Vec<LongitudeTick> = step_deg
        .map(|step| {
            let first_deg = first_tick_deg(west_deg_east, step);
            (0..tick_count(west_deg_east, span_deg, step))
                .map(|index| {
                    #[allow(clippy::cast_precision_loss)]
                    let absolute_deg = (index as f64).mul_add(step, first_deg);
                    tick(west_deg_east, span_deg, absolute_deg)
                })
                .collect()
        })
        .unwrap_or_default();
    if ticks.len() >= 2 {
        return ticks;
    }
    [west_deg_east, west_deg_east + span_deg]
        .into_iter()
        .map(|absolute_deg| tick(west_deg_east, span_deg, absolute_deg))
        .collect()
}

/// The tick at `absolute_deg`, an unfolded longitude east of the western wall.
fn tick(west_deg_east: f64, span_deg: f64, absolute_deg: f64) -> LongitudeTick {
    LongitudeTick {
        x_fraction: if span_deg == 0.0 {
            0.0
        } else {
            (absolute_deg - west_deg_east) / span_deg
        },
        longitude_deg_east: fold(absolute_deg),
    }
}

/// How many meridians a spacing of `step_deg` puts across the basin.
fn tick_count(west_deg_east: f64, span_deg: f64, step_deg: f64) -> usize {
    let first = first_tick_deg(west_deg_east, step_deg);
    if first > west_deg_east + span_deg {
        return 0;
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let count = ((west_deg_east + span_deg - first) / step_deg).floor() as usize + 1;
    count
}

/// The first multiple of `step_deg` at or east of the basin's western wall.
fn first_tick_deg(west_deg_east: f64, step_deg: f64) -> f64 {
    (west_deg_east / step_deg).ceil() * step_deg
}

/// `deg_east` folded into the `[-180, 180)` the basin's bounds are written in.
fn fold(deg_east: f64) -> f64 {
    (deg_east + HALF_TURN_DEG).rem_euclid(FULL_TURN_DEG) - HALF_TURN_DEG
}

/// A longitude as an atlas writes it: a magnitude and a hemisphere.
///
/// The prime meridian and the antimeridian are in neither hemisphere and get
/// neither letter, which is how they are printed and also how a reader tells
/// the middle of the basin from a place 180° away from it.
pub(crate) fn longitude_text(deg_east: f64) -> String {
    let folded = fold(deg_east);
    let magnitude = degrees_text(folded.abs());
    if folded == 0.0 {
        return format!("{magnitude}°");
    }
    // A wall exactly on the antimeridian folds to −180, which is the same
    // meridian as +180 and is written without a hemisphere either way.
    if folded.abs() == HALF_TURN_DEG {
        return format!("{}°", degrees_text(HALF_TURN_DEG));
    }
    let hemisphere = if folded < 0.0 { 'W' } else { 'E' };
    format!("{magnitude}°{hemisphere}")
}

/// A number of degrees, to a tenth, with a trailing `.0` dropped: `120`,
/// `82.5`. A whole-degree meridian is written whole, as an atlas writes it.
fn degrees_text(magnitude_deg: f64) -> String {
    let text = format!("{magnitude_deg:.1}");
    match text.strip_suffix(".0") {
        Some(whole) => whole.to_owned(),
        None => text,
    }
}

/// Where a view of the equator was read, in words: "along the equator" for the
/// basin every scenario declares, and the latitude itself for one laid out
/// some other way.
///
/// Shared by the cross-section and the side view because they are the same
/// reading: the side view is the section's own data drawn as an ocean
/// (`crate::side_view`), so the two must not be able to name different places.
pub(crate) fn latitude_phrase(latitude_deg_north: f64) -> String {
    if latitude_deg_north == 0.0 {
        "along the equator".to_owned()
    } else {
        format!(
            "along {:.2}°{}",
            latitude_deg_north.abs(),
            if latitude_deg_north < 0.0 { 'S' } else { 'N' }
        )
    }
}
