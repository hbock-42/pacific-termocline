//! The run header: everything a reader needs before it opens the frames.

use serde::{Deserialize, Serialize};
use termocline_grid::Grid;

use crate::{FormatError, Variable, VariableSpec, FORMAT_VERSION};

/// Where the basin sits on the globe, in degrees.
///
/// The equatorial Pacific crosses the antimeridian, so `east_deg_east` is
/// numerically smaller than `west_deg_east` for the default basin
/// (120°E to 80°W is `120.0` to `-80.0`); the pair is a west-to-east span, not
/// a min and a max.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BasinExtent {
    /// Western boundary, degrees east of the prime meridian.
    pub west_deg_east: f64,
    /// Eastern boundary, degrees east of the prime meridian.
    pub east_deg_east: f64,
    /// Southern boundary, degrees north of the equator.
    pub south_deg_north: f64,
    /// Northern boundary, degrees north of the equator.
    pub north_deg_north: f64,
}

impl BasinExtent {
    /// The basin spanning `west` to `east` and `south` to `north`, in degrees.
    #[must_use]
    pub const fn new(
        west_deg_east: f64,
        east_deg_east: f64,
        south_deg_north: f64,
        north_deg_north: f64,
    ) -> Self {
        Self {
            west_deg_east,
            east_deg_east,
            south_deg_north,
            north_deg_north,
        }
    }
}

/// The grid a run was computed on: how many cells, and what stretch of ocean
/// they cover.
///
/// A `GridSpec` always describes a basin with at least one cell on each axis,
/// however it was obtained: reading one back from a header runs the same check
/// [`GridSpec::new`] does, so a truncated or hand-edited file fails at the
/// `serde` call with a message rather than panicking later on.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "GridSpecFields")]
pub struct GridSpec {
    nx: usize,
    ny: usize,
    extent: BasinExtent,
}

/// The wire shape of a [`GridSpec`], deserialized before the cell counts have
/// been checked. Serialization goes straight from `GridSpec`, whose derive
/// produces these same three fields.
#[derive(Deserialize)]
struct GridSpecFields {
    nx: usize,
    ny: usize,
    extent: BasinExtent,
}

impl TryFrom<GridSpecFields> for GridSpec {
    type Error = FormatError;

    fn try_from(fields: GridSpecFields) -> Result<Self, Self::Error> {
        Self::new(fields.nx, fields.ny, fields.extent)
    }
}

impl GridSpec {
    /// A basin of `nx` by `ny` cells covering `extent`.
    ///
    /// # Errors
    /// [`FormatError::Grid`] if either cell count is zero.
    pub fn new(nx: usize, ny: usize, extent: BasinExtent) -> Result<Self, FormatError> {
        Grid::new(nx, ny)?;
        Ok(Self { nx, ny, extent })
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

    /// Where the basin sits on the globe.
    #[must_use]
    pub const fn extent(&self) -> BasinExtent {
        self.extent
    }

    /// The cell geometry this spec describes, for code that wants to index
    /// into a field rather than merely size one.
    #[must_use]
    pub fn grid(&self) -> Grid {
        Grid::new(self.nx, self.ny)
            .expect("a GridSpec is only constructible from non-zero cell counts")
    }

    /// How many values a frame's field for `variable` must carry, given where
    /// that variable sits on the C-grid.
    #[must_use]
    pub fn field_len(&self, variable: Variable) -> usize {
        let (nx, ny) = self.grid().field_shape(variable.staggering());
        nx * ny
    }
}

/// How often a run wrote a frame, and how many it wrote.
///
/// The interval is the output cadence, not the solver's timestep: long runs
/// write a decimated series (see T-05.2).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OutputTiming {
    /// Number of frames in the companion frame file.
    pub frame_count: u64,
    /// Model time between consecutive frames, in seconds.
    pub interval_s: f64,
}

/// The physical parameters a run was integrated with.
///
/// These are scenario inputs, not constants of the model: two runs of the same
/// build differ here and nowhere else in this struct.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PhysicalParams {
    /// Mean thermocline depth `H`, in metres — the resting thickness of the
    /// upper layer, which frame values of `h` are anomalies from.
    pub mean_depth_m: f64,
    /// Reduced gravity `g'`, in m s^-2.
    pub reduced_gravity_m_per_s2: f64,
    /// Beta-plane Coriolis gradient `β`, in m^-1 s^-1.
    pub beta_per_m_per_s: f64,
    /// Rayleigh damping coefficient `r`, in s^-1.
    pub rayleigh_damping_per_s: f64,
    /// Reference seawater density `ρ₀`, in kg m^-3.
    pub reference_density_kg_per_m3: f64,
}

/// The header of a run: written once, in JSON, alongside the binary frames.
///
/// It is deliberately self-describing — version, grid, parameters, variable
/// list with units, and output cadence — so a reader never guesses at the
/// shape or meaning of the frames beside it (ADR-0004).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RunHeader {
    /// Version of the on-disk format this run was written with.
    pub format_version: u32,
    /// The grid the run was computed on.
    pub grid: GridSpec,
    /// The physical parameters the run was integrated with.
    pub physical_params: PhysicalParams,
    /// Free text naming the scenario, for a human reading the header.
    pub scenario_description: String,
    /// The variables each frame carries, in frame order, with their units.
    pub variables: Vec<VariableSpec>,
    /// How many frames were written, and how far apart in model time.
    pub output: OutputTiming,
}

impl RunHeader {
    /// The header for a run on `grid` with `physical_params`, stamped with the
    /// current [`FORMAT_VERSION`] and the variable list of the linear core.
    ///
    /// A run that couples SST says so with [`RunHeader::with_sst_anomaly`];
    /// the list is what a reader indexes a run's frames by, so it names the
    /// variables that run actually wrote and not the ones the format knows
    /// how to write.
    #[must_use]
    pub fn new(
        grid: GridSpec,
        physical_params: PhysicalParams,
        scenario_description: impl Into<String>,
        output: OutputTiming,
    ) -> Self {
        Self {
            format_version: FORMAT_VERSION,
            grid,
            physical_params,
            scenario_description: scenario_description.into(),
            variables: Variable::LINEAR_CORE.map(VariableSpec::of).to_vec(),
            output,
        }
    }

    /// The same header, declaring that the run's frames also carry the
    /// mixed-layer SST anomaly `T'` of the Epic 12 coupling.
    ///
    /// Appended after the linear core rather than inserted among it, so that
    /// the first five entries of a coupled run's list are exactly an
    /// uncoupled run's — the extension is additive on the page as well as in
    /// the equations.
    #[must_use]
    pub fn with_sst_anomaly(mut self) -> Self {
        if !self.carries(Variable::SstAnomaly) {
            self.variables.push(VariableSpec::of(Variable::SstAnomaly));
        }
        self
    }

    /// Whether this run's frames carry `variable`.
    ///
    /// Asked of the header rather than of a frame, because it is a fact about
    /// the run: every frame of a run carries the same variables, and a reader
    /// needs to know which before it decodes one.
    #[must_use]
    pub fn carries(&self, variable: Variable) -> bool {
        self.variables.iter().any(|spec| spec.variable == variable)
    }
}
