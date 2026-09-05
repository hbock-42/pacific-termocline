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
//! faces the cell spans, using the point counts
//! [`termocline_grid::Staggering::extra_points`]
//! states, rather than index arithmetic that assumes which axis is staggered.
//!
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

use termocline_format::{FormatError, Frame, GridSpec, Variable};

/// How far the scale reaches when there is no wind at all.
const CALM_PA: f64 = 0.0;

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
            .magnitudes_pa()
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

    /// How long an arrow for `magnitude_pa` is, as a fraction of the longest
    /// the overlay draws.
    ///
    /// A stress past the end of the scale is drawn at full length. That cannot
    /// happen for a scale built over the run being drawn; it can for a caller
    /// that mixed two runs, and a saturated arrow is a better answer than one
    /// drawn across half the basin.
    fn length_fraction(self, magnitude_pa: f64) -> f64 {
        if self.max_magnitude_pa == CALM_PA {
            return 0.0;
        }
        (magnitude_pa / self.max_magnitude_pa).min(1.0)
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
    /// Zonal stress here, in pascals. Easterly — the alizés — is negative.
    tau_x_pa: f64,
    /// Meridional stress here, in pascals. Northward is positive.
    tau_y_pa: f64,
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

    /// Where the arrow ends, when an arrow at the full extent of the scale is
    /// `max_length_cells` long.
    ///
    /// The tip is displaced in the direction the stress pushes the ocean: east
    /// for `τx > 0`, and *up the image* for `τy > 0`, because the field's `y`
    /// increases northward and an image's increases downward.
    #[must_use]
    pub fn tip_cells(&self, max_length_cells: f64) -> (f64, f64) {
        let magnitude_pa = self.magnitude_pa();
        let length_cells = self.length_fraction * max_length_cells;
        (
            self.tau_x_pa
                .mul_add(length_cells / magnitude_pa, self.x_cells),
            (-self.tau_y_pa).mul_add(length_cells / magnitude_pa, self.y_cells),
        )
    }

    /// Zonal stress here, in pascals. Easterly stress is negative
    /// (`CONTEXT.md`, *Wind stress*).
    #[must_use]
    pub const fn tau_x_pa(&self) -> f64 {
        self.tau_x_pa
    }

    /// Meridional stress here, in pascals. Northward is positive.
    #[must_use]
    pub const fn tau_y_pa(&self) -> f64 {
        self.tau_y_pa
    }

    /// The magnitude of the stress here, `√(τx² + τy²)`, in pascals.
    #[must_use]
    pub fn magnitude_pa(&self) -> f64 {
        self.tau_x_pa.hypot(self.tau_y_pa)
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
    /// The arrows of `frame`'s wind stress over `grid`, on `scale`, one every
    /// `spacing_cells` cells along each axis.
    ///
    /// The stress is a field, so there is one value per cell and 32 000 of them
    /// for the control basin — far more arrows than a reader can separate, and
    /// enough ink to hide the map underneath. `spacing_cells` is how many are
    /// drawn; the first is half a spacing in from the northwest corner, so the
    /// outermost arrows sit inside the map rather than on its rim.
    ///
    /// Cells whose stress is calm or non-finite get no arrow. A zero-length
    /// arrow is a dot, and a scatter of dots over the map reads as data; a
    /// non-finite stress has no direction to point in at all.
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] if `frame` does not fit `grid`.
    ///
    /// # Panics
    /// If `spacing_cells` is zero. The spacing is the shell's own layout
    /// constant rather than anything a run carries, so a zero means this code
    /// is wrong, not that the input is.
    pub fn of_frame(
        grid: GridSpec,
        frame: &Frame,
        scale: StressScale,
        spacing_cells: usize,
    ) -> Result<Self, FormatError> {
        assert!(
            spacing_cells > 0,
            "arrows cannot be spaced zero cells apart"
        );
        let stress = CellCentreStress::of_frame(grid, frame)?;
        let (width, height) = (stress.width, stress.height);
        let mut arrows = Vec::new();
        // Row 0 of the map is the northernmost, which is the last row of the
        // field — the same flip the heatmap makes.
        for row in (spacing_cells / 2..height).step_by(spacing_cells) {
            let j = height - 1 - row;
            for i in (spacing_cells / 2..width).step_by(spacing_cells) {
                let (tau_x_pa, tau_y_pa) = stress.at(i, j);
                let magnitude_pa = tau_x_pa.hypot(tau_y_pa);
                if !magnitude_pa.is_finite() || magnitude_pa == CALM_PA {
                    continue;
                }
                let length_fraction = scale.length_fraction(magnitude_pa);
                if length_fraction == 0.0 {
                    continue;
                }
                #[allow(clippy::cast_precision_loss)]
                let (x_cells, y_cells) = (i as f64 + 0.5, row as f64 + 0.5);
                arrows.push(WindArrow {
                    x_cells,
                    y_cells,
                    tau_x_pa,
                    tau_y_pa,
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

/// A frame's wind stress moved off the C-grid faces it is stored on and onto
/// the cell centres an arrow can be drawn at.
struct CellCentreStress {
    /// Cells along x.
    width: usize,
    /// Cells along y.
    height: usize,
    /// `τx` at each cell centre, row-major with `j` increasing northward.
    tau_x_pa: Vec<f64>,
    /// `τy` at each cell centre, in the same order.
    tau_y_pa: Vec<f64>,
}

impl CellCentreStress {
    /// `frame`'s stress at the centre of every cell of `grid`.
    ///
    /// # Errors
    /// [`FormatError::FieldShape`] if `frame` does not fit `grid`.
    fn of_frame(grid: GridSpec, frame: &Frame) -> Result<Self, FormatError> {
        frame.validate(&grid)?;
        let cells = grid.grid();
        Ok(Self {
            width: cells.nx(),
            height: cells.ny(),
            tau_x_pa: at_cell_centres(grid, frame, Variable::ZonalWindStress),
            tau_y_pa: at_cell_centres(grid, frame, Variable::MeridionalWindStress),
        })
    }

    /// The stress `(τx, τy)` at the centre of cell `(i, j)`, with `j` counted
    /// northward from the southern edge of the basin.
    fn at(&self, i: usize, j: usize) -> (f64, f64) {
        let offset = j * self.width + i;
        (self.tau_x_pa[offset], self.tau_y_pa[offset])
    }

    /// The magnitude of the stress at every cell centre, in pascals.
    fn magnitudes_pa(&self) -> impl Iterator<Item = f64> + '_ {
        self.tau_x_pa
            .iter()
            .zip(&self.tau_y_pa)
            .map(|(tau_x_pa, tau_y_pa)| tau_x_pa.hypot(*tau_y_pa))
    }
}

/// `variable`'s field, averaged from wherever it is staggered onto the cell
/// centres, row-major with `j` increasing northward.
///
/// A field carries one extra line of points on the axis it is staggered on
/// ([`termocline_grid::Staggering::extra_points`]), so a cell spans two of its
/// points along that axis and one along the other. Averaging over exactly the
/// points the staggering claims is why this is one function rather than one per
/// component: nothing here has to know which axis `τx` is staggered on.
///
/// `frame` must already have been validated against `grid`.
fn at_cell_centres(grid: GridSpec, frame: &Frame, variable: Variable) -> Vec<f64> {
    let cells = grid.grid();
    let (width, height) = (cells.nx(), cells.ny());
    let staggering = variable.staggering();
    let (points_x, _) = cells.field_shape(staggering);
    let (extra_x, extra_y) = staggering.extra_points();
    let values = frame.field(variable);
    #[allow(clippy::cast_precision_loss)]
    let points_per_cell = ((extra_x + 1) * (extra_y + 1)) as f64;
    let mut centres = Vec::with_capacity(width * height);
    for j in 0..height {
        for i in 0..width {
            let mut sum = 0.0;
            for offset_y in 0..=extra_y {
                for offset_x in 0..=extra_x {
                    sum += values[(j + offset_y) * points_x + (i + offset_x)];
                }
            }
            centres.push(sum / points_per_cell);
        }
    }
    centres
}
