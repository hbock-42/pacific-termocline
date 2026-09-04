//! Why a header or a frame could not be built.

use std::fmt;

use termocline_grid::GridError;

use crate::Variable;

/// Why a header or a frame could not be built from the values supplied.
///
/// Every variant describes invalid *input* — a scenario asking for a
/// degenerate basin, or field data that does not cover the grid it claims to
/// be on — so it is returned rather than panicked, and names both the
/// offending value and the bound it violated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FormatError {
    /// The grid the header describes is not a valid basin.
    Grid(GridError),
    /// A frame's field did not carry one value per point of its staggered
    /// position on the grid.
    FieldShape {
        /// The variable whose field was the wrong length.
        variable: Variable,
        /// Values the grid asks for at this variable's staggering.
        expected: usize,
        /// Values actually supplied.
        actual: usize,
    },
}

impl fmt::Display for FormatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Grid(err) => err.fmt(f),
            Self::FieldShape {
                variable,
                expected,
                actual,
            } => write!(
                f,
                "field {}: expected {expected} values for this grid, got {actual}",
                variable.symbol()
            ),
        }
    }
}

impl std::error::Error for FormatError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Grid(err) => Some(err),
            Self::FieldShape { .. } => None,
        }
    }
}

impl From<GridError> for FormatError {
    fn from(err: GridError) -> Self {
        Self::Grid(err)
    }
}
