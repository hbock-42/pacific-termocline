//! The wind-stress overlay: one frame's `τx, τy` as arrows over the basin map.
//!
//! The alizés are what drives the response the heatmap draws, and the two are
//! only readable together — a thermocline that has stopped tilting means one
//! thing under steady trades and another under a wind burst. So the overlay is
//! a *layer*, not a second picture: its arrows are laid out in the same cell
//! coordinates [`crate::Heatmap`] draws in, tail at a cell centre, and the
//! shell paints them over the map it already uploaded.
//!
//! Like the heatmap, none of this knows what a GPU is. A [`WindOverlay`] is a
//! list of [`WindArrow`]s with positions and lengths in cells, so what a reader
//! would check by eye — arrows pointing west along the equator — is checked by
//! index in `tests/wind_overlay.rs`, and the same code draws the layer in a
//! browser and natively ([ADR-0006]).
//!
//! # Why the arrows are geometry rather than pixels
//!
//! Nothing here touches the heatmap's buffer. That is what makes the second
//! acceptance criterion of T-09.1 — toggling the overlay does not affect the
//! map under it — structural rather than something to be careful about: the
//! map is a texture built from `h` alone, the arrows are shapes drawn on top of
//! it, and hiding a shape cannot reach a texture.
//!
//! # Where an arrow's stress comes from
//!
//! `τx` lives on the cells' east/west faces and `τy` on their north/south
//! faces (ADR-0003), because each forces the current that lives there. Neither
//! is defined at a cell centre, which is where an arrow has to be drawn if it
//! is to sit on the map's own grid. So each component is averaged from the
//! faces the cell spans, using the point counts [`Staggering::extra_points`]
//! states, rather than index arithmetic that assumes which axis is staggered.
//! The averaging reads the frame's buffers where they lie rather than
//! materializing a field of its own: a run is walked once per load and once per
//! drawn frame, and a browser tab is short of exactly that (ADR-0006).
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use termocline_format::{FormatError, Frame, GridSpec, Variable};
use termocline_grid::{Grid, Staggering};

/// How far the scale reaches when there is no wind at all.
const CALM_PA: f64 = 0.0;

/// Cells between neighbouring arrows, along both axes.
///
/// The stress is a field: one arrow per cell would be 32 000 of them over the
/// control basin, which is ink rather than information. Twelve of that basin's
/// half-degree cells is six degrees of arc.
pub const ARROW_SPACING_CELLS: usize = 12;

/// Length in cells of an arrow drawing the strongest stress in the run.
///
/// Shorter than [`ARROW_SPACING_CELLS`], so that even a basin of arrows at full
/// length does not run into itself.
pub const MAX_ARROW_LENGTH_CELLS: f64 = 9.0;

/// Below this length, in cells, an arrow is not drawn at all.
///
/// An arrow shorter than the cell it sits on has no readable direction: it is a
/// dot, and a scatter of dots over the map reads as data the run does not
/// carry. The trades fall off as a Gaussian, so without this the quiet
/// two-thirds of the basin would be stippled with marks standing for stresses
/// of order 10^-17 Pa.
pub const MIN_ARROW_LENGTH_CELLS: f64 = 1.0;

/// The stress magnitude an arrow of full length stands for.
///
/// The run's, not the frame's, for the reason [`crate::DivergingScale`] is the
/// run's: an overlay rescaled frame by frame would draw every frame at full
/// length, and a wind burst relaxing back to the trades would look like a wind
/// that never moved.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StressScale {
    /// The largest stress magnitude the scale covers, in pascals (N m⁻²).
    /// Never negative; zero for an ocean under no wind at all.
    max_magnitude_pa: f64,
}

impl StressScale {
    /// The scale of an ocean under no wind: nothing to draw, and no length to
    /// draw it at.
    ///
    /// The seed a run-wide scale is folded from, one frame at a time.
    #[must_use]
    pub const fn calm() -> Self {
        Self {
            max_magnitude_pa: CALM_PA,
        }
    }

    /// The scale that just covers the stress `frame` carries over `grid`.
    ///
    /// Non-finite stresses are ignored rather than propagated: one `NaN` would
    /// otherwise leave the whole run with no scale, and so no overlay at all.
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] if `frame` does not fit `grid` — the frame
    /// of one run against the header of another.
    pub fn covering(grid: GridSpec, frame: &Frame) -> Result<Self, FormatError> {
        let stress = CellCentreStress::of_frame(grid, frame)?;
        let max_magnitude_pa = stress
            .every_cell()
            .map(|stress| stress.magnitude_pa())
            .filter(|magnitude| magnitude.is_finite())
            .fold(CALM_PA, f64::max);
        Ok(Self { max_magnitude_pa })
    }

    /// The scale that covers both this one and `other`.
    #[must_use]
    pub fn widened(self, other: Self) -> Self {
        Self {
            max_magnitude_pa: self.max_magnitude_pa.max(other.max_magnitude_pa),
        }
    }

    /// The largest stress magnitude the scale covers, in pascals (N m⁻²), so
    /// the shell can label the overlay with a number a reader would quote.
    #[must_use]
    pub const fn max_magnitude_pa(&self) -> f64 {
        self.max_magnitude_pa
    }

    /// How long an arrow for `stress` is, as a fraction of the longest the
    /// overlay draws, or `None` if it is not worth drawing at all.
    ///
    /// Nothing is drawn for a calm cell, for one whose stress is non-finite —
    /// which has no direction to point in — or for one whose arrow would come
    /// out shorter than [`MIN_ARROW_LENGTH_CELLS`].
    ///
    /// A stress past the end of the scale is drawn at full length. That cannot
    /// happen for a scale built over the run being drawn; it can for a caller
    /// that mixed two runs, and a saturated arrow is a better answer than one
    /// drawn across half the basin.
    fn drawable_length_fraction(self, stress: Stress) -> Option<f64> {
        let magnitude_pa = stress.magnitude_pa();
        if !magnitude_pa.is_finite() || self.max_magnitude_pa == CALM_PA {
            return None;
        }
        let fraction = (magnitude_pa / self.max_magnitude_pa).min(1.0);
        (fraction * MAX_ARROW_LENGTH_CELLS >= MIN_ARROW_LENGTH_CELLS).then_some(fraction)
    }
}

/// One arrow of the overlay: the wind stress at one cell centre, and where to
/// draw it.
///
/// Positions are in cells of the basin map, measured from its northwest corner
/// with `y` increasing southward — the same coordinates [`crate::Heatmap`]'s
/// pixels are in, so the shell places an arrow by scaling the map's drawn
/// rectangle and nothing else.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindArrow {
    /// Cells east of the map's western edge, at the arrow's tail.
    x_cells: f64,
    /// Cells south of the map's northern edge, at the arrow's tail.
    y_cells: f64,
    /// The stress the arrow draws.
    stress: Stress,
    /// This arrow's length, as a fraction of the longest the overlay draws.
    /// Strictly positive: an arrow with nothing to say is not built.
    length_fraction: f64,
}

impl WindArrow {
    /// Where the arrow starts: the centre of the cell whose stress it draws,
    /// as `(east, south)` in cells from the map's northwest corner.
    #[must_use]
    pub const fn tail_cells(&self) -> (f64, f64) {
        (self.x_cells, self.y_cells)
    }

    /// Where the arrow ends, in the same cell coordinates as its tail.
    ///
    /// The tip is displaced in the direction the stress pushes the ocean: east
    /// for `τx > 0`, and *up the image* for `τy > 0`, because the field's `y`
    /// increases northward and an image's increases downward. How far is
    /// [`MAX_ARROW_LENGTH_CELLS`] scaled by this arrow's share of the run's
    /// stress scale.
    #[must_use]
    pub fn tip_cells(&self) -> (f64, f64) {
        let per_pascal = self.length_cells() / self.stress.magnitude_pa();
        (
            self.stress.tau_x_pa.mul_add(per_pascal, self.x_cells),
            (-self.stress.tau_y_pa).mul_add(per_pascal, self.y_cells),
        )
    }

    /// How long the arrow is drawn, in cells. Never below
    /// [`MIN_ARROW_LENGTH_CELLS`]: a shorter one is not built.
    #[must_use]
    pub fn length_cells(&self) -> f64 {
        self.length_fraction * MAX_ARROW_LENGTH_CELLS
    }

    /// Zonal stress here, in pascals. Easterly stress is negative
    /// (`CONTEXT.md`, *Wind stress*).
    #[must_use]
    pub const fn tau_x_pa(&self) -> f64 {
        self.stress.tau_x_pa
    }

    /// Meridional stress here, in pascals. Northward is positive.
    #[must_use]
    pub const fn tau_y_pa(&self) -> f64 {
        self.stress.tau_y_pa
    }

    /// The magnitude of the stress here, `√(τx² + τy²)`, in pascals.
    #[must_use]
    pub fn magnitude_pa(&self) -> f64 {
        self.stress.magnitude_pa()
    }

    /// This arrow's length as a fraction of the longest the overlay draws:
    /// the stress here against [`StressScale::max_magnitude_pa`].
    #[must_use]
    pub const fn length_fraction(&self) -> f64 {
        self.length_fraction
    }
}

/// One frame's wind stress, as arrows over the basin map.
#[derive(Debug, Clone)]
pub struct WindOverlay {
    /// The arrows, in reading order from the northwest corner.
    arrows: Vec<WindArrow>,
    /// The scale their lengths came from, for the legend beside them.
    scale: StressScale,
}

impl WindOverlay {
    /// The arrows of `frame`'s wind stress over `grid`, on `scale`.
    ///
    /// The stress is a field, so there is one value per cell and 32 000 of them
    /// for the control basin — far more arrows than a reader can separate, and
    /// enough ink to hide the map underneath. One cell in
    /// [`ARROW_SPACING_CELLS`] carries an arrow instead.
    ///
    /// Which cells is anchored on the *middle* of the basin rather than on its
    /// northwestern corner, so that a basin laid out symmetrically about the
    /// equator — which every scenario's is (`CONTEXT.md`, *Basin*) — gets a row
    /// of arrows along the equator, where the trades are strongest and where
    /// the response they drive is. Anchoring at a corner instead leaves the
    /// equator between two rows, and the strongest forcing in the run undrawn.
    ///
    /// A cell gets no arrow when its stress is calm, non-finite, or too weak to
    /// draw ([`MIN_ARROW_LENGTH_CELLS`]).
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] if `frame` does not fit `grid`.
    pub fn of_frame(
        grid: GridSpec,
        frame: &Frame,
        scale: StressScale,
    ) -> Result<Self, FormatError> {
        let stress = CellCentreStress::of_frame(grid, frame)?;
        let (width, height) = (stress.width(), stress.height());
        let mut arrows = Vec::new();
        // Row 0 of the map is the northernmost, which is the last row of the
        // field — the same flip the heatmap makes.
        for row in arrow_line(height) {
            let j = height - 1 - row;
            for i in arrow_line(width) {
                let stress = stress.at(i, j);
                let Some(length_fraction) = scale.drawable_length_fraction(stress) else {
                    continue;
                };
                #[allow(clippy::cast_precision_loss)]
                let (x_cells, y_cells) = (i as f64 + 0.5, row as f64 + 0.5);
                arrows.push(WindArrow {
                    x_cells,
                    y_cells,
                    stress,
                    length_fraction,
                });
            }
        }
        Ok(Self { arrows, scale })
    }

    /// The arrows to draw, in reading order from the map's northwest corner.
    #[must_use]
    pub fn arrows(&self) -> &[WindArrow] {
        &self.arrows
    }

    /// The scale these lengths came from.
    #[must_use]
    pub const fn scale(&self) -> StressScale {
        self.scale
    }
}

/// The wind stress at one point: the two components, together, because
/// nothing reads one without the other and the magnitude is of the pair.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Stress {
    /// Zonal stress, in pascals. Easterly — the alizés — is negative.
    tau_x_pa: f64,
    /// Meridional stress, in pascals. Northward is positive.
    tau_y_pa: f64,
}

impl Stress {
    /// `√(τx² + τy²)`, in pascals.
    fn magnitude_pa(self) -> f64 {
        self.tau_x_pa.hypot(self.tau_y_pa)
    }
}

/// A frame's wind stress, read off the C-grid faces it is stored on at the
/// cell centres an arrow can be drawn at.
///
/// It borrows the frame rather than building a field of its own: the whole of
/// a run is walked this way at load to find the stress scale, and one frame of
/// it again for every frame drawn.
struct CellCentreStress<'a> {
    /// The basin, in cells.
    cells: Grid,
    /// `τx`, on the cells' east/west faces.
    tau_x: FaceField<'a>,
    /// `τy`, on the cells' north/south faces.
    tau_y: FaceField<'a>,
}

impl<'a> CellCentreStress<'a> {
    /// A reader for `frame`'s stress over `grid`.
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] if `frame` does not fit `grid`.
    fn of_frame(grid: GridSpec, frame: &'a Frame) -> Result<Self, FormatError> {
        frame.validate(&grid)?;
        let cells = grid.grid();
        Ok(Self {
            cells,
            tau_x: FaceField::of(cells, frame, Variable::ZonalWindStress),
            tau_y: FaceField::of(cells, frame, Variable::MeridionalWindStress),
        })
    }

    /// Cells along x.
    const fn width(&self) -> usize {
        self.cells.nx()
    }

    /// Cells along y.
    const fn height(&self) -> usize {
        self.cells.ny()
    }

    /// The stress at the centre of cell `(i, j)`, with `j` counted northward
    /// from the southern edge of the basin.
    fn at(&self, i: usize, j: usize) -> Stress {
        Stress {
            tau_x_pa: self.tau_x.at_cell(i, j),
            tau_y_pa: self.tau_y.at_cell(i, j),
        }
    }

    /// The stress at the centre of every cell of the basin.
    fn every_cell(&self) -> impl Iterator<Item = Stress> + '_ {
        (0..self.height()).flat_map(move |j| (0..self.width()).map(move |i| self.at(i, j)))
    }
}

/// One component of the forcing, as it is stored: on the C-grid faces the
/// current it drives lives on, and read back at cell centres.
///
/// A field carries one extra line of points on the axis it is staggered on
/// ([`Staggering::extra_points`]), so a cell spans two of its points along that
/// axis and one along the other. Averaging over exactly the points the
/// staggering claims is why this is one type rather than one per component:
/// nothing here has to know which axis `τx` is staggered on.
struct FaceField<'a> {
    /// The values, row-major with `j` increasing northward.
    values: &'a [f64],
    /// Points along x, which is the row stride.
    points_x: usize,
    /// Extra points beyond the cell count, as `(along x, along y)`.
    extra: (usize, usize),
}

impl<'a> FaceField<'a> {
    /// `variable`'s field of `frame`, which must already have been validated
    /// against the grid `cells` describes.
    ///
    /// # Panics
    /// If `frame` does not carry `variable`. Only the two wind-stress
    /// components are read here and every frame carries both
    /// ([`Variable::LINEAR_CORE`]), so an absent one means this code asked for
    /// the wrong variable rather than that the run is short of a field.
    fn of(cells: Grid, frame: &'a Frame, variable: Variable) -> Self {
        let values = frame
            .field(variable)
            .expect("every frame carries the variables of the linear core");
        Self::at(cells, values, variable.staggering())
    }

    /// `values`, read as a field staggered at `staggering` over `cells`.
    fn at(cells: Grid, values: &'a [f64], staggering: Staggering) -> Self {
        let (points_x, _) = cells.field_shape(staggering);
        Self {
            values,
            points_x,
            extra: staggering.extra_points(),
        }
    }

    /// The value at the centre of cell `(i, j)`: the mean of the points the
    /// cell spans.
    fn at_cell(&self, i: usize, j: usize) -> f64 {
        let (extra_x, extra_y) = self.extra;
        let mut sum = 0.0;
        for offset_y in 0..=extra_y {
            for offset_x in 0..=extra_x {
                sum += self.values[(j + offset_y) * self.points_x + (i + offset_x)];
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let points_per_cell = ((extra_x + 1) * (extra_y + 1)) as f64;
        sum / points_per_cell
    }
}

/// The lines of cells that carry arrows along an axis `extent` cells long.
///
/// Anchored on the middle of the axis rather than on its start, so that the
/// pattern is symmetric about the basin's centre — which for every scenario's
/// basin puts a line on the equator.
fn arrow_line(extent: usize) -> impl Iterator<Item = usize> {
    (extent / 2 % ARROW_SPACING_CELLS..extent).step_by(ARROW_SPACING_CELLS)
}
