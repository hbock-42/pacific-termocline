//! Shared 2D field and grid types: what a grid cell is, and where each
//! prognostic variable sits on it.
//!
//! This crate is deliberately physics-free — it holds data structures and
//! indexing math only, so the engine, the tests and (later) the visualizer
//! share one definition instead of three. The staggering it encodes is the
//! Arakawa C-grid chosen in [ADR-0003]: `h` at cell centers, `u` at cell
//! east/west faces, `v` at cell north/south faces. Solver code addresses those
//! positions through [`Staggering`] rather than through raw `+1`/`-1` index
//! arithmetic.
//!
//! There is no metric geometry here — no cell spacing in metres, no basin
//! origin. Where the basin sits and how wide a cell is are scenario
//! parameters (see `CONTEXT.md`, *Basin*), and they arrive with the physics in
//! Epic 01.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

pub mod sweep;

use std::fmt;

/// One of the grid's two axes, named so an error can say which one it means.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Axis {
    /// The zonal axis, indexed by `i`, increasing eastward.
    X,
    /// The meridional axis, indexed by `j`, increasing northward.
    Y,
}

impl fmt::Display for Axis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::X => "x",
            Self::Y => "y",
        })
    }
}

/// Where a variable sits relative to a grid cell on the Arakawa C-grid.
///
/// The three variants are the three positions the scheme uses; the named
/// constants [`H_STAGGERING`], [`U_STAGGERING`] and [`V_STAGGERING`] bind each
/// prognostic variable of the 1.5-layer model to one of them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Staggering {
    /// The center of the cell, where the thermocline depth anomaly `h` lives.
    CellCenter,
    /// The cell's western edge, where the zonal current anomaly `u` lives.
    /// Indexing runs `0..=nx`, so the eastern basin boundary is addressable.
    EastWestFace,
    /// The cell's southern edge, where the meridional current anomaly `v`
    /// lives. Indexing runs `0..=ny`, so the northern boundary is addressable.
    NorthSouthFace,
}

/// Staggering of the thermocline depth anomaly `h`, per ADR-0003.
pub const H_STAGGERING: Staggering = Staggering::CellCenter;
/// Staggering of the zonal current anomaly `u`, per ADR-0003.
pub const U_STAGGERING: Staggering = Staggering::EastWestFace;
/// Staggering of the meridional current anomaly `v`, per ADR-0003.
pub const V_STAGGERING: Staggering = Staggering::NorthSouthFace;

/// Offset of a cell center from the cell's southwest corner, in cell widths.
const HALF_CELL: f64 = 0.5;
/// Offset of a face value along the axis it is staggered on, in cell widths:
/// the face is the corner itself.
const ON_FACE: f64 = 0.0;

impl Staggering {
    /// Offset of this position from the cell's southwest corner, as
    /// `(x, y)` fractions of a cell width and a cell height.
    ///
    /// Dimensionless on purpose: converting to a length needs a cell spacing,
    /// which is a scenario parameter this crate does not carry.
    #[must_use]
    pub const fn offset_in_cells(self) -> (f64, f64) {
        match self {
            Self::CellCenter => (HALF_CELL, HALF_CELL),
            Self::EastWestFace => (ON_FACE, HALF_CELL),
            Self::NorthSouthFace => (HALF_CELL, ON_FACE),
        }
    }

    /// Extra points a field at this position carries beyond the cell count, as
    /// `(along x, along y)`. A face field needs one more line of points than
    /// there are cells on the axis it is staggered on, so that both basin
    /// boundaries are represented.
    #[must_use]
    pub const fn extra_points(self) -> (usize, usize) {
        match self {
            Self::CellCenter => (0, 0),
            Self::EastWestFace => (1, 0),
            Self::NorthSouthFace => (0, 1),
        }
    }
}

/// Why a grid or a field could not be built.
///
/// These describe invalid *input* — a scenario asking for a degenerate basin —
/// so they are returned rather than panicked, and each names the offending
/// value and the bound it violated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GridError {
    /// A cell count was zero: a basin needs at least one cell on each axis.
    EmptyExtent {
        /// The axis whose extent was zero.
        axis: Axis,
    },
    /// A backing buffer did not hold exactly `nx * ny` values.
    ShapeMismatch {
        /// Points expected, `nx * ny`.
        expected: usize,
        /// Values actually supplied.
        actual: usize,
    },
}

impl fmt::Display for GridError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyExtent { axis } => {
                write!(f, "n{axis} is 0; a grid needs at least 1 cell on each axis")
            }
            Self::ShapeMismatch { expected, actual } => {
                write!(f, "expected {expected} values for this shape, got {actual}")
            }
        }
    }
}

impl std::error::Error for GridError {}

/// A dense 2D array of values, one per point, stored row-major.
///
/// "Row-major" means `x` varies fastest: the flat offset of `(i, j)` is
/// `j * nx + i`. Nothing here knows what the values mean — a `Field2D` holds
/// `h` in metres just as happily as it holds a mask of booleans.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Field2D<T> {
    nx: usize,
    ny: usize,
    values: Vec<T>,
}

impl<T: Clone> Field2D<T> {
    /// A field of `nx` by `ny` points, every one set to `value`.
    ///
    /// # Errors
    /// [`GridError::EmptyExtent`] if either extent is zero.
    pub fn filled(nx: usize, ny: usize, value: T) -> Result<Self, GridError> {
        check_extents(nx, ny)?;
        Ok(Self {
            nx,
            ny,
            values: vec![value; nx * ny],
        })
    }
}

impl<T> Field2D<T> {
    /// A field of `nx` by `ny` points wrapping an existing row-major buffer.
    ///
    /// # Errors
    /// [`GridError::EmptyExtent`] if either extent is zero, or
    /// [`GridError::ShapeMismatch`] if `values` does not hold exactly
    /// `nx * ny` entries.
    pub fn from_vec(nx: usize, ny: usize, values: Vec<T>) -> Result<Self, GridError> {
        check_extents(nx, ny)?;
        let expected = nx * ny;
        if values.len() != expected {
            return Err(GridError::ShapeMismatch {
                expected,
                actual: values.len(),
            });
        }
        Ok(Self { nx, ny, values })
    }

    /// Number of points along x.
    #[must_use]
    pub const fn nx(&self) -> usize {
        self.nx
    }

    /// Number of points along y.
    #[must_use]
    pub const fn ny(&self) -> usize {
        self.ny
    }

    /// Total number of points, `nx * ny`.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.nx * self.ny
    }

    /// Whether the field holds no points — always `false`, since both extents
    /// are checked to be non-zero at construction.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Flat offset of the point `(i, j)`, or `None` if it lies outside the
    /// field.
    #[must_use]
    pub const fn flat_index(&self, i: usize, j: usize) -> Option<usize> {
        if i >= self.nx || j >= self.ny {
            return None;
        }
        Some(j * self.nx + i)
    }

    /// The `(i, j)` this flat offset addresses — the inverse of
    /// [`Field2D::flat_index`] — or `None` if the offset is past the end.
    #[must_use]
    pub const fn cell_of(&self, flat: usize) -> Option<(usize, usize)> {
        if flat >= self.len() {
            return None;
        }
        Some((flat % self.nx, flat / self.nx))
    }

    /// The value at `(i, j)`, or `None` if it lies outside the field.
    #[must_use]
    pub fn get(&self, i: usize, j: usize) -> Option<&T> {
        self.flat_index(i, j).map(|flat| &self.values[flat])
    }

    /// Mutable access to the value at `(i, j)`, or `None` if it lies outside
    /// the field.
    pub fn get_mut(&mut self, i: usize, j: usize) -> Option<&mut T> {
        self.flat_index(i, j).map(|flat| &mut self.values[flat])
    }

    /// The `j`th row of the field: its `nx` values, contiguous.
    ///
    /// The read side of [`sweep::write_rows`], which hands a kernel one output
    /// row at a time and leaves it to reach for the input rows it needs.
    ///
    /// # Panics
    /// If `j` is past the last row. A kernel indexes rows the field's own
    /// shape defines, so an out-of-range row means the calling code is wrong,
    /// which is what panics are for.
    #[must_use]
    pub fn row(&self, j: usize) -> &[T] {
        assert!(
            j < self.ny,
            "row {j} is outside a field of {} rows",
            self.ny
        );
        &self.values[j * self.nx..][..self.nx]
    }

    /// The backing buffer, row-major.
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        &self.values
    }

    /// The backing buffer, row-major and mutable. Time-stepping writes through
    /// this rather than reallocating per step.
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.values
    }
}

/// The shape of a rectangular basin in cells, and the shapes the staggered
/// fields over it must have.
///
/// The grid carries no physics — it does not know that `h` is a depth anomaly,
/// only that `h` lives at cell centers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Grid {
    nx: usize,
    ny: usize,
}

impl Grid {
    /// A basin of `nx` by `ny` cells.
    ///
    /// # Errors
    /// [`GridError::EmptyExtent`] if either cell count is zero.
    pub fn new(nx: usize, ny: usize) -> Result<Self, GridError> {
        check_extents(nx, ny)?;
        Ok(Self { nx, ny })
    }

    /// Number of cells along x.
    #[must_use]
    pub const fn nx(&self) -> usize {
        self.nx
    }

    /// Number of cells along y.
    #[must_use]
    pub const fn ny(&self) -> usize {
        self.ny
    }

    /// Shape a field at this staggering must have to cover the basin, as
    /// `(points along x, points along y)`.
    #[must_use]
    pub const fn field_shape(&self, staggering: Staggering) -> (usize, usize) {
        let (extra_x, extra_y) = staggering.extra_points();
        (self.nx + extra_x, self.ny + extra_y)
    }

    /// A field covering the basin at this staggering, every point set to
    /// `value`.
    #[must_use]
    pub fn allocate<T: Clone>(&self, staggering: Staggering, value: T) -> Field2D<T> {
        let (nx, ny) = self.field_shape(staggering);
        Field2D::filled(nx, ny, value)
            .expect("a grid's extents are non-zero, so its field shapes are too")
    }
}

fn check_extents(nx: usize, ny: usize) -> Result<(), GridError> {
    if nx == 0 {
        return Err(GridError::EmptyExtent { axis: Axis::X });
    }
    if ny == 0 {
        return Err(GridError::EmptyExtent { axis: Axis::Y });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejection_messages_name_the_offending_extent() {
        // Invalid input must be actionable, per CODING_STANDARDS.md: the
        // message states which axis and what bound it violated.
        let err = Grid::new(4, 0).expect_err("a basin needs at least one row");
        assert_eq!(err, GridError::EmptyExtent { axis: Axis::Y });
        let message = err.to_string();
        assert!(message.contains("ny is 0"), "{message}");

        let err = Field2D::from_vec(2, 2, vec![0.0_f64; 3]).expect_err("3 values do not fill 2x2");
        assert_eq!(
            err,
            GridError::ShapeMismatch {
                expected: 4,
                actual: 3
            }
        );
    }

    #[test]
    fn the_mutable_buffer_is_the_same_storage_as_the_indexed_view() {
        // Time stepping writes through the slice and reads back by cell; the
        // two views must not drift apart.
        let mut field = Field2D::filled(3, 2, 0.0_f64).expect("valid shape");
        field.as_mut_slice()[4] = 2.0;
        assert_eq!(field.get(1, 1), Some(&2.0));
    }
}
