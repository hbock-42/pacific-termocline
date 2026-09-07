//! What a scenario looks like on disk, and how it becomes engine types.
//!
//! `CONTEXT.md` defines a *scenario* as "a complete, runnable specification of
//! one simulation: grid, physical parameters, wind-forcing description, and
//! run length. The engine's unit of input." This module is that definition
//! written down twice: [`ScenarioConfig`] is the TOML shape, a plain `serde`
//! record with no invariants, and [`Scenario`] is the validated result — a
//! [`Basin`], a [`PhysicalParams`], an ordered list of wind forcings and an
//! [`OutputSchedule`], each already through the constructor that checks it.
//!
//! Keeping the two apart is what makes the error messages good. Every
//! constructor in the engine already refuses its own bad input by name
//! ([`PhysicalParamsError`], [`WindStressError`], [`CflError`], …), so
//! [`ScenarioConfig::build`] re-implements almost nothing: it wires the file's
//! numbers into those constructors and wraps whichever one objected. A
//! malformed file, an unknown forcing type or an unstable timestep is
//! therefore a [`ScenarioError`] naming the offending value — never a panic
//! (CODING_STANDARDS.md § *Correctness and failure*).
//!
//! # Everything is checked before anything runs
//!
//! [`ScenarioConfig::build`] is the pre-flight check of T-06.3, and its
//! contract is that an accepted scenario is a scenario the run will not stop
//! on: no bound is left for the solver, the writer or the allocator to
//! discover partway through a long run. Two of them are here for that reason
//! rather than because this is where they are computed — the rotation bound on
//! `dt_s`, which [`Solver::new`](crate::Solver::new) also enforces, and the
//! memory budget on the grid, which is the one bound this module owns
//! outright because no single constructor can see both the cell count and what
//! a run does with it.
//!
//! # Why this is not in `termocline-format`
//!
//! ADR-0004 makes `termocline-format` the one place a file format is defined,
//! but the format it is about is the engine's *output* — the header and frames
//! the visualizer reads, the contract between two processes that never share a
//! type. A scenario is *input*, read by the engine alone, and its whole job is
//! to name engine concepts: the [`WindStress`] implementations of `forcing`,
//! the [`PhysicalParams`] of `params`, the [`Basin`] of `basin`. Putting it in
//! `termocline-format` would either invert the dependency — the format crate
//! is explicitly free of simulation logic, and both `engine` and `visualizer`
//! depend on it — or duplicate every constructor there and re-validate in two
//! places. So the scenario format lives here, in the crate whose vocabulary it
//! speaks, and ADR-0004's rule keeps its scope: one place per format, and this
//! is the only place this one is defined.
//!
//! # The file
//!
//! ```toml
//! [basin]                      # every key optional; the whole section too
//! western_longitude_deg = 120.0   # 120°E
//! eastern_longitude_deg = -80.0   # 80°W, counted eastward across the dateline
//! southern_latitude_deg = -25.0
//! northern_latitude_deg = 25.0
//! resolution_deg = 0.5            # cell size, both axes
//!
//! [physics]
//! reduced_gravity_m_per_s2 = 0.06
//! mean_thermocline_depth_m = 150.0
//! rayleigh_damping_per_s = 1.0e-7
//! # beta_per_m_per_s            — optional, defaults to EQUATORIAL_BETA_PER_M_PER_S
//! # reference_density_kg_per_m3 — optional, defaults to SEAWATER_REFERENCE_DENSITY_KG_PER_M3
//!
//! [run]
//! dt_s = 3600.0
//! total_steps = 17520
//! output_every_n_steps = 24
//!
//! [[wind]]                     # zero or more, summed in the order written
//! type = "steady_trade_winds"
//! equatorial_zonal_stress_pa = -0.05
//! meridional_decay_scale_m = 361000.0
//! ```
//!
//! The sketch above is a tour, not the specification: every field, its unit,
//! its default and the bound behind it are in
//! `docs/scenario-config-reference.md`, which
//! `engine/tests/scenario_config_reference.rs` holds to this module field by
//! field. Three worked examples live in `engine/scenarios/`, one per scenario
//! of `docs/planning/01-scientific-model.md`.

use std::fmt;
#[cfg(feature = "fs")]
use std::fs;
#[cfg(feature = "fs")]
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use termocline_numerics::{check_timestep, CflError, WaveSpeed};

use crate::basin::{Basin, BasinBounds, BasinBoundsError};
use crate::coriolis::BetaPlane;
use crate::forcing::{
    CompositeWind, SeasonalTradeWinds, SteadyTradeWinds, TimeDependence, WindBurstAnomaly,
    WindStress, WindStressError,
};
use crate::params::{
    PhysicalParams, PhysicalParamsError, EQUATORIAL_BETA_PER_M_PER_S,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
};
use crate::run_writer::{OutputSchedule, OutputScheduleError};
use crate::solver::{check_rotation_timestep, RotationLimitError};
use crate::sst::{SstParams, SstParamsError, DEFAULT_SURFACE_DRAG_PER_S};
use crate::wind_response::{
    WindResponseError, WindResponseParams, DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M,
};

/// Why a scenario could not be read.
///
/// All but one variant describe invalid *input* — a file that is not a
/// scenario, or a scenario the engine cannot run — so they are returned rather
/// than panicked; [`ScenarioError::Unwritable`] is the exception, and reports
/// a scenario the engine holds that TOML cannot express. The variants that
/// wrap another crate's error do so precisely to keep that error's message,
/// which already names the offending value and the bound it violated; this
/// type only says which part of the file it came from.
#[derive(Debug)]
pub enum ScenarioError {
    /// The scenario file could not be read from disk.
    ///
    /// Only a build with the `fs` feature can reach a file to fail to read
    /// (ADR-0012); parsing a scenario from text cannot produce this.
    #[cfg(feature = "fs")]
    Unreadable {
        /// The path that was asked for.
        path: PathBuf,
        /// What the filesystem said.
        source: std::io::Error,
    },
    /// The file was not valid TOML, was missing a section, named a forcing
    /// that does not exist, or carried a key the format does not define.
    Malformed(toml::de::Error),
    /// A valid scenario could not be written back out as TOML.
    Unwritable(toml::ser::Error),
    /// `[basin]` described something that is not a basin: a boundary that is
    /// not a position on the planet, or a resolution that does not divide it.
    Basin(BasinBoundsError),
    /// `[basin]` described a grid too fine for the engine to hold a run on.
    ///
    /// A perfectly well-formed basin can still ask for more cells than the
    /// machine has memory for, and the engine says so before it starts
    /// allocating rather than letting the allocator end the run partway
    /// through it.
    BasinTooLarge {
        /// Cells east–west.
        nx: usize,
        /// Cells north–south.
        ny: usize,
        /// Bytes of solver state a run over this grid would hold resident.
        resident_bytes: u64,
    },
    /// `[physics]` asked for an unphysical ocean.
    PhysicalParams(PhysicalParamsError),
    /// A `[[wind]]` entry described a forcing that cannot exist.
    Wind(WindStressError),
    /// `[run]` asked for a timestep or an output cadence that is not a
    /// schedule.
    Schedule(OutputScheduleError),
    /// `[run]` asked for a timestep the grid cannot carry stably.
    Cfl(CflError),
    /// `[run]` asked for a timestep longer than the basin's rotation allows
    /// (ADR-0007).
    Rotation(RotationLimitError),
    /// `[sst]` asked for a mixed layer that cannot exist.
    Sst(SstParamsError),
    /// `[sst]` asked for an atmospheric wind response that cannot exist.
    WindResponse(WindResponseError),
}

impl fmt::Display for ScenarioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            #[cfg(feature = "fs")]
            Self::Unreadable { path, source } => {
                write!(
                    f,
                    "could not read the scenario {}: {source}",
                    path.display()
                )
            }
            Self::Malformed(source) => write!(f, "this is not a scenario: {source}"),
            Self::Unwritable(source) => write!(f, "could not write the scenario: {source}"),
            Self::Basin(source) => write!(f, "[basin]: {source}"),
            Self::BasinTooLarge {
                nx,
                ny,
                resident_bytes,
            } => write!(
                f,
                "[basin]: the basin is {nx} × {ny} cells, whose solver state alone would be \
                 {} — more than the {} this build will start a run with; coarsen \
                 resolution_deg, or bring the basin's boundaries closer together",
                gibibytes(*resident_bytes),
                gibibytes(MAX_RESIDENT_STATE_BYTES)
            ),
            Self::PhysicalParams(source) => write!(f, "[physics]: {source}"),
            Self::Wind(source) => write!(f, "[[wind]]: {source}"),
            Self::Schedule(source) => write!(f, "[run]: {source}"),
            Self::Cfl(source) => write!(f, "[run]: {source}"),
            Self::Rotation(source) => write!(f, "[run]: {source}"),
            Self::Sst(source) => write!(f, "[sst]: {source}"),
            Self::WindResponse(source) => write!(f, "[sst]: {source}"),
        }
    }
}

/// `bytes` as a number of gibibytes, to a tenth: the unit the memory budget is
/// stated in, and the one a reader compares against the machine in front of
/// them.
fn gibibytes(bytes: u64) -> String {
    format!("{:.1} GiB", bytes as f64 / BYTES_PER_GIB as f64)
}

/// Bytes in a gibibyte.
const BYTES_PER_GIB: u64 = 1 << 30;

/// Bytes of solver state a run holds resident per grid cell.
///
/// Counted off what a run holds alive from before its first step to after its
/// last (`run.rs` § *What is allocated, and when*): the state and RK4's five
/// stage buffers are six [`OceanState`](crate::OceanState)s of three fields
/// each, and the two [`WindStressField`](crate::WindStressField)s carry two
/// components apiece — 18 + 4 = 22 `f64` a cell, rounded up to 24 for the
/// extra row and column the staggered `u` and `v` fields carry. At 8 bytes a
/// `f64` that is 192 bytes a cell.
///
/// It is an estimate of the resident set, not a measurement of it: what the
/// budget below is protecting is the difference between a basin that fits in
/// memory and one that is three orders of magnitude past it, and no plausible
/// error in the buffer count moves that line.
const RESIDENT_BYTES_PER_CELL: u64 = 24 * 8;

/// Bytes of solver state a run with the Epic 12 SST coupling holds resident
/// per grid cell.
///
/// The linear core's [`RESIDENT_BYTES_PER_CELL`] plus what the coupling adds:
/// the SST anomaly itself in the state and in RK4's five stage buffers (6
/// fields), and the [`SstTerm`](crate::SstTerm)'s eight — the two mixed-layer
/// velocity components, the two interpolated stress components, the two
/// divergence halves, the upwelling, and the zonal current on cell centers.
/// That is 14 more `f64` a cell, plus the two components of the extra
/// [`WindStressField`](crate::WindStressField) a coupled forcing sums the
/// prescribed winds and the atmospheric response into (T-12.2) — 16, rounded
/// up to 20 for the extra rows and columns the five staggered ones carry, on
/// the same reasoning as the count above.
///
/// Counted separately rather than folded into one worst case so that a
/// scenario of the validated linear model is held to the budget it has always
/// been held to: turning the coupling on is what costs the memory, and it is
/// the only scenario that should pay for it.
const COUPLED_RESIDENT_BYTES_PER_CELL: u64 = RESIDENT_BYTES_PER_CELL + 20 * 8;

/// The largest resident solver state this build will start a run with.
///
/// A project policy rather than a measured constant, in the same sense as
/// `CFL_SAFETY_FACTOR`: 2 GiB is comfortably inside any machine this project
/// is developed or run on, and it admits 11.2 million cells — 350 times the
/// 320 × 100 of the default Pacific, or that basin at 0.03°, which is far
/// finer than the deformation radius the model resolves. A scenario past it is
/// not a scenario with an ambitious grid, it is a scenario with a mistyped
/// `resolution_deg`.
///
/// The limit is refused, never silently coarsened
/// (CODING_STANDARDS.md § *No silent clamping*).
const MAX_RESIDENT_STATE_BYTES: u64 = 2 * BYTES_PER_GIB;

impl std::error::Error for ScenarioError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            #[cfg(feature = "fs")]
            Self::Unreadable { source, .. } => Some(source),
            Self::Malformed(source) => Some(source),
            Self::Unwritable(source) => Some(source),
            Self::Basin(source) => Some(source),
            Self::PhysicalParams(source) => Some(source),
            Self::Wind(source) => Some(source),
            Self::Schedule(source) => Some(source),
            Self::Cfl(source) => Some(source),
            Self::Rotation(source) => Some(source),
            Self::Sst(source) => Some(source),
            Self::WindResponse(source) => Some(source),
            Self::BasinTooLarge { .. } => None,
        }
    }
}

impl From<toml::de::Error> for ScenarioError {
    fn from(source: toml::de::Error) -> Self {
        Self::Malformed(source)
    }
}

impl From<BasinBoundsError> for ScenarioError {
    fn from(source: BasinBoundsError) -> Self {
        Self::Basin(source)
    }
}

impl From<PhysicalParamsError> for ScenarioError {
    fn from(source: PhysicalParamsError) -> Self {
        Self::PhysicalParams(source)
    }
}

impl From<WindStressError> for ScenarioError {
    fn from(source: WindStressError) -> Self {
        Self::Wind(source)
    }
}

impl From<OutputScheduleError> for ScenarioError {
    fn from(source: OutputScheduleError) -> Self {
        Self::Schedule(source)
    }
}

impl From<CflError> for ScenarioError {
    fn from(source: CflError) -> Self {
        Self::Cfl(source)
    }
}

impl From<RotationLimitError> for ScenarioError {
    fn from(source: RotationLimitError) -> Self {
        Self::Rotation(source)
    }
}

impl From<SstParamsError> for ScenarioError {
    fn from(source: SstParamsError) -> Self {
        Self::Sst(source)
    }
}

impl From<WindResponseError> for ScenarioError {
    fn from(source: WindResponseError) -> Self {
        Self::WindResponse(source)
    }
}

/// The `[basin]` section: which part of the ocean the scenario runs on, in
/// degrees, and how finely it is cut into cells.
///
/// Every key is optional and defaults to the equatorial Pacific of
/// `CONTEXT.md` (*Basin*), so a scenario states a bound only when it means
/// something other than the basin this project is about.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct BasinSection {
    /// Western boundary, in degrees east. Omitted means
    /// [`PACIFIC_WESTERN_LONGITUDE_DEG`](crate::basin::PACIFIC_WESTERN_LONGITUDE_DEG).
    pub western_longitude_deg: f64,
    /// Eastern boundary, in degrees east — counted eastward from the western
    /// one, so `-80.0` and `280.0` are the same meridian and a basin may
    /// cross the dateline. Omitted means [`PACIFIC_EASTERN_LONGITUDE_DEG`](crate::basin::PACIFIC_EASTERN_LONGITUDE_DEG).
    pub eastern_longitude_deg: f64,
    /// Southern boundary, in degrees north. Omitted means
    /// [`PACIFIC_SOUTHERN_LATITUDE_DEG`](crate::basin::PACIFIC_SOUTHERN_LATITUDE_DEG).
    pub southern_latitude_deg: f64,
    /// Northern boundary, in degrees north. Omitted means
    /// [`PACIFIC_NORTHERN_LATITUDE_DEG`](crate::basin::PACIFIC_NORTHERN_LATITUDE_DEG).
    pub northern_latitude_deg: f64,
    /// Cell size, in degrees, on both axes. Omitted means
    /// [`PACIFIC_RESOLUTION_DEG`](crate::basin::PACIFIC_RESOLUTION_DEG).
    ///
    /// One number, not two: on the equatorial beta-plane a degree of longitude
    /// and a degree of latitude are the same degree of arc, so a basin stated
    /// in degrees has square cells. A scenario wanting a coarser zonal step
    /// than meridional one cannot say so here, which is deliberate — an
    /// anisotropic grid is a numerical decision, and it would arrive with the
    /// ADR that justifies it rather than as a second key.
    pub resolution_deg: f64,
}

impl Default for BasinSection {
    /// The Pacific, read off [`BasinBounds::pacific`] rather than restated
    /// here so the file format and the basin cannot drift apart about what the
    /// default basin is.
    fn default() -> Self {
        Self::of(BasinBounds::pacific())
    }
}

impl BasinSection {
    /// The section that states `bounds`.
    #[must_use]
    pub fn of(bounds: BasinBounds) -> Self {
        Self {
            western_longitude_deg: bounds.western_longitude_deg(),
            eastern_longitude_deg: bounds.eastern_longitude_deg(),
            southern_latitude_deg: bounds.southern_latitude_deg(),
            northern_latitude_deg: bounds.northern_latitude_deg(),
            resolution_deg: bounds.resolution_deg(),
        }
    }

    /// The [`BasinBounds`] this section describes.
    ///
    /// # Errors
    /// [`ScenarioError::Basin`], naming the boundary or the resolution that is
    /// not one.
    pub fn bounds(&self) -> Result<BasinBounds, ScenarioError> {
        Ok(BasinBounds::new(
            self.western_longitude_deg,
            self.eastern_longitude_deg,
            self.southern_latitude_deg,
            self.northern_latitude_deg,
            self.resolution_deg,
        )?)
    }
}

/// The `[physics]` section: the constants of the scenario's ocean.
///
/// `β` and `ρ₀` are properties of the planet rather than of the experiment, so
/// they default to the named constants of [`crate::params`] and a scenario
/// only states them when it deliberately varies them.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhysicsSection {
    /// Reduced gravity `g'`, in m/s².
    pub reduced_gravity_m_per_s2: f64,
    /// Mean thermocline depth `H`, in metres.
    pub mean_thermocline_depth_m: f64,
    /// Rayleigh damping coefficient `r`, in s⁻¹.
    pub rayleigh_damping_per_s: f64,
    /// Meridional gradient of the Coriolis parameter `β`, in m⁻¹s⁻¹. Omitted
    /// means [`EQUATORIAL_BETA_PER_M_PER_S`].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub beta_per_m_per_s: Option<f64>,
    /// Reference seawater density `ρ₀`, in kg/m³. Omitted means
    /// [`SEAWATER_REFERENCE_DENSITY_KG_PER_M3`].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_density_kg_per_m3: Option<f64>,
}

impl PhysicsSection {
    /// The [`PhysicalParams`] this section describes.
    ///
    /// # Errors
    /// [`ScenarioError::PhysicalParams`], naming the first parameter outside
    /// its bound.
    pub fn build(&self) -> Result<PhysicalParams, ScenarioError> {
        Ok(PhysicalParams::new(
            self.reduced_gravity_m_per_s2,
            self.mean_thermocline_depth_m,
            self.rayleigh_damping_per_s,
            self.beta_per_m_per_s.unwrap_or(EQUATORIAL_BETA_PER_M_PER_S),
            self.reference_density_kg_per_m3
                .unwrap_or(SEAWATER_REFERENCE_DENSITY_KG_PER_M3),
        )?)
    }
}

/// The `[run]` section: how long the run is and how often it is saved.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunSection {
    /// Length of one solver step, in seconds.
    pub dt_s: f64,
    /// Steps the run takes from its initial state — the run length, in steps
    /// rather than in seconds, because that is what makes the frame count
    /// exact.
    pub total_steps: u64,
    /// Steps between saved frames.
    pub output_every_n_steps: u64,
}

impl RunSection {
    /// The [`OutputSchedule`] this section describes.
    ///
    /// # Errors
    /// [`ScenarioError::Schedule`] for a non-positive timestep or a cadence of
    /// zero.
    pub fn build(&self) -> Result<OutputSchedule, ScenarioError> {
        Ok(OutputSchedule::new(
            self.dt_s,
            self.total_steps,
            self.output_every_n_steps,
        )?)
    }
}

/// One `[[wind]]` entry: which forcing, and with what parameters.
///
/// The `type` key names the [`WindStress`] implementation, so an unknown
/// forcing is a `serde` unknown-variant error that lists the ones that do
/// exist. [`CompositeWind`] is not among them: a composite is the *list* of
/// entries, not an entry, so nesting is not expressible and does not need to
/// be — the components sum, and a sum is flat.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum WindSection {
    /// [`SteadyTradeWinds`]: the control scenario.
    SteadyTradeWinds {
        /// Zonal stress `τ₀` on the equator, in Pa. Strictly negative.
        equatorial_zonal_stress_pa: f64,
        /// Meridional decay scale `Ly`, in metres. Omitted means a field with
        /// no meridional structure at all — the `Ly → ∞` limit, not a large
        /// `Ly`.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        meridional_decay_scale_m: Option<f64>,
    },
    /// [`SeasonalTradeWinds`]: the same field scaled by an annual harmonic.
    SeasonalTradeWinds {
        /// Zonal stress `τ₀` on the equator, in Pa, before modulation.
        equatorial_zonal_stress_pa: f64,
        /// Meridional decay scale `Ly`, in metres. Omitted means no
        /// meridional structure.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        meridional_decay_scale_m: Option<f64>,
        /// Relative amplitude `a` of the annual harmonic, dimensionless, in
        /// `[0, 1]`.
        relative_amplitude: f64,
        /// The instant the alizés are strongest, in seconds into the run.
        peak_time_s: f64,
    },
    /// [`WindBurstAnomaly`]: a westerly anomaly stacked on whatever else the
    /// list carries.
    WindBurstAnomaly {
        /// Peak zonal stress `τ_burst`, in Pa. Strictly positive.
        peak_zonal_stress_pa: f64,
        /// Zonal centre `x₀` of the burst, in metres.
        center_x_m: f64,
        /// Zonal `e`-folding scale `Lx`, in metres.
        zonal_scale_m: f64,
        /// Meridional `e`-folding scale `Ly`, in metres, about the equator.
        meridional_scale_m: f64,
        /// Instant `t₀` of the burst's peak, in seconds into the run.
        peak_time_s: f64,
        /// Temporal `e`-folding scale `Lt`, in seconds.
        duration_s: f64,
    },
}

impl WindSection {
    /// The forcing this entry describes.
    ///
    /// # Errors
    /// [`ScenarioError::Wind`], naming the parameter the forcing refused and
    /// the bound it violated.
    pub fn build(&self) -> Result<ScenarioWind, ScenarioError> {
        Ok(match *self {
            Self::SteadyTradeWinds {
                equatorial_zonal_stress_pa,
                meridional_decay_scale_m,
            } => ScenarioWind::Steady(steady(
                equatorial_zonal_stress_pa,
                meridional_decay_scale_m,
            )?),
            Self::SeasonalTradeWinds {
                equatorial_zonal_stress_pa,
                meridional_decay_scale_m,
                relative_amplitude,
                peak_time_s,
            } => ScenarioWind::Seasonal(SeasonalTradeWinds::new(
                steady(equatorial_zonal_stress_pa, meridional_decay_scale_m)?,
                relative_amplitude,
                peak_time_s,
            )?),
            Self::WindBurstAnomaly {
                peak_zonal_stress_pa,
                center_x_m,
                zonal_scale_m,
                meridional_scale_m,
                peak_time_s,
                duration_s,
            } => ScenarioWind::Burst(WindBurstAnomaly::new(
                peak_zonal_stress_pa,
                center_x_m,
                zonal_scale_m,
                meridional_scale_m,
                peak_time_s,
                duration_s,
            )?),
        })
    }
}

/// The steady field a `[[wind]]` entry describes, with or without a decay
/// scale. Shared by the two entries that carry one.
fn steady(
    equatorial_zonal_stress_pa: f64,
    meridional_decay_scale_m: Option<f64>,
) -> Result<SteadyTradeWinds, WindStressError> {
    match meridional_decay_scale_m {
        None => SteadyTradeWinds::uniform(equatorial_zonal_stress_pa),
        Some(scale_m) => {
            SteadyTradeWinds::with_meridional_decay(equatorial_zonal_stress_pa, scale_m)
        }
    }
}

/// One validated wind forcing of a scenario: whichever [`WindStress`] its
/// `[[wind]]` entry named.
///
/// A closed enum rather than a `Box<dyn WindStress>` because a scenario is
/// worth inspecting — the CLI reports it, the tests assert on it — and a trait
/// object would have thrown away exactly the parameters a reader wants to see.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ScenarioWind {
    /// The control scenario.
    Steady(SteadyTradeWinds),
    /// The seasonal cycle.
    Seasonal(SeasonalTradeWinds),
    /// A westerly wind burst.
    Burst(WindBurstAnomaly),
}

impl WindStress for ScenarioWind {
    fn stress(&self, x_m: f64, y_m: f64, t_s: f64) -> (f64, f64) {
        match self {
            Self::Steady(wind) => wind.stress(x_m, y_m, t_s),
            Self::Seasonal(wind) => wind.stress(x_m, y_m, t_s),
            Self::Burst(wind) => wind.stress(x_m, y_m, t_s),
        }
    }

    /// Whatever the wind the scenario named says. The enum is a dispatcher and
    /// has no opinion of its own — least of all here, where an opinion would
    /// be a promise about a forcing it did not write.
    fn time_dependence(&self) -> TimeDependence {
        match self {
            Self::Steady(wind) => wind.time_dependence(),
            Self::Seasonal(wind) => wind.time_dependence(),
            Self::Burst(wind) => wind.time_dependence(),
        }
    }
}

/// The `[sst]` section: the Epic 12 mixed-layer coupling, and the switch that
/// turns it on.
///
/// Omitting the section entirely is what keeps a scenario the validated linear
/// model of Epics 01-07 — `CONTEXT.md` puts the SST anomaly outside the ocean
/// core, and so does this format. A scenario that carries the section gets the
/// fourth prognostic variable `T'` and the equation of [`crate::sst`]; one
/// that does not gets exactly the three-variable run it always got.
///
/// The section is the config option the ticket's deliverable asks for, and it
/// is opt-in rather than a boolean flag beside a set of always-present
/// parameters: there is no defensible default mixed layer, and a `enabled =
/// false` sitting above five numbers nobody chose would be a scenario claiming
/// more than it means.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SstSection {
    /// Mixed-layer depth `H_m`, in metres.
    pub mixed_layer_depth_m: f64,
    /// Zonal gradient of the mean SST, `∂T̄/∂x`, in K/m.
    pub mean_zonal_sst_gradient_k_per_m: f64,
    /// Sensitivity `γ = ∂T_sub/∂h` of the entrained water to the thermocline
    /// depth anomaly, in K/m.
    pub subsurface_temperature_sensitivity_k_per_m: f64,
    /// Thermal damping `ε_T` of an SST anomaly, in s⁻¹.
    pub thermal_damping_per_s: f64,
    /// Rayleigh drag `r_s` of the wind-driven surface layer, in s⁻¹. Omitted
    /// means [`DEFAULT_SURFACE_DRAG_PER_S`].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub surface_drag_per_s: Option<f64>,
    /// Feedback strength `μ` of the T-12.2 atmospheric wind response, in Pa/K.
    /// Omitted means zero — the prescribed-wind model T-12.1 left, in which
    /// the coupling runs one way and the alizés are whatever the `[[wind]]`
    /// entries say.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wind_feedback_strength_pa_per_k: Option<f64>,
    /// Meridional scale `L_a` of that response, in metres. Omitted means
    /// [`DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M`].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wind_response_meridional_scale_m: Option<f64>,
}

impl SstSection {
    /// The validated mixed-layer parameters this section describes.
    ///
    /// # Errors
    /// [`ScenarioError::Sst`], wrapping the [`SstParamsError`] that names the
    /// offending value and the bound it violated.
    pub fn build(&self) -> Result<SstParams, ScenarioError> {
        Ok(SstParams::new(
            self.mixed_layer_depth_m,
            self.surface_drag_per_s
                .unwrap_or(DEFAULT_SURFACE_DRAG_PER_S),
            self.mean_zonal_sst_gradient_k_per_m,
            self.subsurface_temperature_sensitivity_k_per_m,
            self.thermal_damping_per_s,
        )?)
    }

    /// The validated atmospheric wind response this section describes.
    ///
    /// Every coupled scenario has one, because "no feedback" is a strength of
    /// zero rather than an absent object: the loop is open at `μ = 0` and the
    /// run is the one T-12.1 validated, bit for bit, which is exactly the
    /// claim T-12.2's acceptance criterion makes.
    ///
    /// # Errors
    /// [`ScenarioError::WindResponse`], wrapping the [`WindResponseError`]
    /// that names the offending value and the bound it violated.
    pub fn build_wind_response(&self) -> Result<WindResponseParams, ScenarioError> {
        Ok(WindResponseParams::new(
            self.wind_feedback_strength_pa_per_k.unwrap_or(0.0),
            self.wind_response_meridional_scale_m
                .unwrap_or(DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M),
        )?)
    }
}

/// A scenario as it is written in a file: five sections, no invariants.
///
/// This is the `serde` record and nothing more. It round-trips through TOML,
/// so a run can record the scenario that produced it, and it becomes a
/// runnable [`Scenario`] only through [`ScenarioConfig::build`], which is
/// where every bound is checked.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScenarioConfig {
    /// The `[basin]` section. Omitted entirely means the Pacific.
    #[serde(default)]
    pub basin: BasinSection,
    /// The `[physics]` section.
    pub physics: PhysicsSection,
    /// The `[run]` section.
    pub run: RunSection,
    /// The `[[wind]]` entries, in the order they are summed. An empty list is
    /// a calm ocean, which is the undriven limit of
    /// `docs/planning/01-scientific-model.md` rather than a mistake.
    #[serde(default)]
    pub wind: Vec<WindSection>,
    /// The `[sst]` section. Omitted is the validated linear model of
    /// Epics 01-07; present switches on the Epic 12 mixed-layer coupling.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sst: Option<SstSection>,
}

impl ScenarioConfig {
    /// Parse a scenario from the text of a TOML file.
    ///
    /// # Errors
    /// [`ScenarioError::Malformed`] for a file that is not valid TOML, is
    /// missing a section, names a forcing that does not exist, or carries a
    /// key this format does not define — the last because silently ignoring a
    /// misspelled parameter would run a scenario nobody asked for.
    pub fn from_toml(source: &str) -> Result<Self, ScenarioError> {
        Ok(toml::from_str(source)?)
    }

    /// This scenario as the text of a TOML file.
    ///
    /// # Errors
    /// [`ScenarioError::Unwritable`] if the config holds a value TOML cannot
    /// represent, such as a non-finite float.
    pub fn to_toml(&self) -> Result<String, ScenarioError> {
        toml::to_string(self).map_err(ScenarioError::Unwritable)
    }

    /// Validate every section and produce the runnable [`Scenario`].
    ///
    /// This is the whole pre-flight check of T-06.3: when it returns `Ok`,
    /// every bound the engine holds a scenario to has been cleared, so a file
    /// the loader accepts is a file the run will not stop on. That is why the
    /// two timestep bounds are here rather than left to the solver, and why
    /// the grid-size budget is checked before anything allocates: the point of
    /// validating up front is that nothing downstream gets to discover a new
    /// objection partway through a long run.
    ///
    /// The order is the one `docs/scenario-config-reference.md` § *Errors*
    /// documents, and the first failure is the one reported. The two bounds on
    /// `dt_s` come last because each needs more than `[run]`: the gravity-wave
    /// CFL bound needs the wave speed `[physics]` implies and the spacing
    /// `[basin]` implies, and the rotation bound of [ADR-0007] needs `β` and
    /// how far from the equator `[basin]` reaches. Neither ever shortens the
    /// timestep (CODING_STANDARDS.md § *No silent clamping*).
    ///
    /// # Errors
    /// A [`ScenarioError`] naming the section that was wrong and, through the
    /// error it wraps, the value and the bound.
    ///
    /// [ADR-0007]: ../../docs/planning/adr/0007-rotation-timestep-bound.md
    pub fn build(&self) -> Result<Scenario, ScenarioError> {
        let bounds = self.basin.bounds()?;
        let sst_params = self.sst.map(|sst| sst.build()).transpose()?;
        let wind_response_params = self.sst.map(|sst| sst.build_wind_response()).transpose()?;
        check_grid_fits_in_memory(bounds.nx(), bounds.ny(), sst_params.is_some())?;
        let basin = bounds.basin();
        let physical_params = self.physics.build()?;
        let schedule = self.run.build()?;
        let winds = self
            .wind
            .iter()
            .map(WindSection::build)
            .collect::<Result<Vec<_>, _>>()?;

        let wave_speed = WaveSpeed::new(physical_params.kelvin_wave_speed_m_per_s())?;
        check_timestep(schedule.dt_s(), basin.spacing(), wave_speed)?;
        let plane = BetaPlane::of_basin(physical_params, basin);
        check_rotation_timestep(schedule.dt_s(), basin.grid(), plane)?;

        Ok(Scenario {
            bounds,
            basin,
            physical_params,
            winds,
            output_schedule: schedule,
            sst_params,
            wind_response_params,
        })
    }
}

/// Refuse an `nx` × `ny` grid whose run would not fit in
/// [`MAX_RESIDENT_STATE_BYTES`].
///
/// The count is done in `u64` and saturates rather than wrapping, so a grid
/// far past the budget is refused for being past it rather than accepted for
/// having overflowed — on a 32-bit target `nx · ny · 192` does not fit in a
/// `usize`.
///
/// # Errors
/// [`ScenarioError::BasinTooLarge`], naming the grid and what it would cost.
fn check_grid_fits_in_memory(nx: usize, ny: usize, couples_sst: bool) -> Result<(), ScenarioError> {
    let cells = (nx as u64).saturating_mul(ny as u64);
    let resident_bytes = cells.saturating_mul(if couples_sst {
        COUPLED_RESIDENT_BYTES_PER_CELL
    } else {
        RESIDENT_BYTES_PER_CELL
    });
    if resident_bytes > MAX_RESIDENT_STATE_BYTES {
        return Err(ScenarioError::BasinTooLarge {
            nx,
            ny,
            resident_bytes,
        });
    }
    Ok(())
}

/// A scenario the engine can run: every value already through the constructor
/// that checks it.
///
/// The engine's unit of input (`CONTEXT.md`, *Scenario*). Holding the winds as
/// an ordered `Vec` rather than as an assembled [`CompositeWind`] keeps the
/// scenario inspectable and keeps the sum's order — and so the run's
/// floating-point result — fixed by the file (CODING_STANDARDS.md §
/// *Correctness and failure*).
#[derive(Debug, Clone, PartialEq)]
pub struct Scenario {
    bounds: BasinBounds,
    basin: Basin,
    physical_params: PhysicalParams,
    winds: Vec<ScenarioWind>,
    output_schedule: OutputSchedule,
    sst_params: Option<SstParams>,
    wind_response_params: Option<WindResponseParams>,
}

impl Scenario {
    /// Read and validate a scenario from a TOML file on disk.
    ///
    /// # Errors
    /// [`ScenarioError::Unreadable`] naming `path` if the file cannot be read,
    /// otherwise whatever [`ScenarioConfig::from_toml`] or
    /// [`ScenarioConfig::build`] objected to.
    #[cfg(feature = "fs")]
    pub fn load(path: &Path) -> Result<Self, ScenarioError> {
        let source = fs::read_to_string(path).map_err(|source| ScenarioError::Unreadable {
            path: path.to_path_buf(),
            source,
        })?;
        Self::from_toml(&source)
    }

    /// Parse and validate a scenario from the text of a TOML file.
    ///
    /// # Errors
    /// Whatever [`ScenarioConfig::from_toml`] or [`ScenarioConfig::build`]
    /// objected to.
    pub fn from_toml(source: &str) -> Result<Self, ScenarioError> {
        ScenarioConfig::from_toml(source)?.build()
    }

    /// Where this scenario's basin is on the globe, in degrees, and how
    /// finely it is cut into cells.
    ///
    /// Kept beside the [`Basin`] rather than derived from it because the two
    /// answer different questions and only one of them can answer this one: a
    /// `Basin` is the truncation in metres the equations are solved on, with
    /// `x` measured east from its own western wall, and it has forgotten which
    /// meridians that wall lies between. A run's header records the basin in
    /// degrees, so the degrees have to survive the build.
    #[must_use]
    pub const fn bounds(&self) -> BasinBounds {
        self.bounds
    }

    /// Where this scenario's basin is, and how big its cells are.
    #[must_use]
    pub const fn basin(&self) -> Basin {
        self.basin
    }

    /// The constants of this scenario's ocean.
    #[must_use]
    pub const fn physical_params(&self) -> PhysicalParams {
        self.physical_params
    }

    /// The wind forcings of this scenario, in the order they are summed.
    #[must_use]
    pub fn winds(&self) -> &[ScenarioWind] {
        &self.winds
    }

    /// The forcings as the one [`WindStress`] the solver reads: their sum, in
    /// file order.
    ///
    /// An empty `[[wind]]` list gives an empty composite, which is calm.
    ///
    /// This allocates one box per component, so it is a once-per-run call:
    /// assemble the composite before the first step and reuse it, never inside
    /// the time-stepping loop (CODING_STANDARDS.md § *Performance*).
    #[must_use]
    pub fn wind(&self) -> CompositeWind {
        self.winds
            .iter()
            .fold(CompositeWind::new(), |composite, wind| {
                composite.with(*wind)
            })
    }

    /// How long this run is and how often it is saved.
    #[must_use]
    pub const fn output_schedule(&self) -> OutputSchedule {
        self.output_schedule
    }

    /// The Epic 12 mixed-layer coupling this scenario asked for, or `None` for
    /// the validated linear model of Epics 01-07.
    ///
    /// The one place a run learns whether it integrates `T'`: a `None` here is
    /// what makes the run allocate three prognostic fields and step the
    /// three-variable right-hand side it always did.
    #[must_use]
    pub const fn sst_params(&self) -> Option<SstParams> {
        self.sst_params
    }

    /// The T-12.2 atmospheric wind response of this scenario, or `None` when
    /// there is no `[sst]` section to respond to.
    ///
    /// `Some` with a zero feedback strength is the open loop — the answer for
    /// every coupled scenario written before this ticket, and the one that
    /// runs exactly as it did.
    #[must_use]
    pub const fn wind_response_params(&self) -> Option<WindResponseParams> {
        self.wind_response_params
    }
}
