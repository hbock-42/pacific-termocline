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
    /// A frame and its header disagreed about which variables the run wrote:
    /// the frame carried a variable the header does not declare, or lacked one
    /// it does.
    ///
    /// The header's variable list is what a reader indexes a run by, so a
    /// frame that does not match it would be read under the wrong labels — or,
    /// worse, would offer a field the run never announced.
    UndeclaredVariable {
        /// The variable the two disagree about.
        variable: Variable,
        /// Whether the header declares it.
        declared: bool,
    },
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
            Self::UndeclaredVariable {
                variable,
                declared: true,
            } => write!(
                f,
                "the header lists {} among the run's variables and the frame does not carry it",
                variable.symbol()
            ),
            Self::UndeclaredVariable {
                variable,
                declared: false,
            } => write!(
                f,
                "the frame carries {} and the header does not list it among the run's variables",
                variable.symbol()
            ),
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
            Self::UndeclaredVariable { .. } | Self::FieldShape { .. } => None,
        }
    }
}

impl From<GridError> for FormatError {
    fn from(err: GridError) -> Self {
        Self::Grid(err)
    }
}
