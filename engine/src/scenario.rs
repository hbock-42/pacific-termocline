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
//! [`ScenarioConfig::build`] does not re-implement a single bound: it wires
//! the file's numbers into those constructors and wraps whichever one
//! objected. A malformed file, an unknown forcing type or an unstable
//! timestep is therefore a [`ScenarioError`] naming the offending value —
//! never a panic (CODING_STANDARDS.md § *Correctness and failure*).
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
use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use termocline_numerics::{check_timestep, CflError, WaveSpeed};

use crate::basin::{Basin, BasinBounds, BasinBoundsError};
use crate::forcing::{
    CompositeWind, SeasonalTradeWinds, SteadyTradeWinds, WindBurstAnomaly, WindStress,
    WindStressError,
};
use crate::params::{
    PhysicalParams, PhysicalParamsError, EQUATORIAL_BETA_PER_M_PER_S,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
};
use crate::run_writer::{OutputSchedule, OutputScheduleError};

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
    /// `[physics]` asked for an unphysical ocean.
    PhysicalParams(PhysicalParamsError),
    /// A `[[wind]]` entry described a forcing that cannot exist.
    Wind(WindStressError),
    /// `[run]` asked for a timestep or an output cadence that is not a
    /// schedule.
    Schedule(OutputScheduleError),
    /// `[run]` asked for a timestep the grid cannot carry stably.
    Cfl(CflError),
}

impl fmt::Display for ScenarioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
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
            Self::PhysicalParams(source) => write!(f, "[physics]: {source}"),
            Self::Wind(source) => write!(f, "[[wind]]: {source}"),
            Self::Schedule(source) => write!(f, "[run]: {source}"),
            Self::Cfl(source) => write!(f, "[run]: {source}"),
        }
    }
}

impl std::error::Error for ScenarioError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unreadable { source, .. } => Some(source),
            Self::Malformed(source) => Some(source),
            Self::Unwritable(source) => Some(source),
            Self::Basin(source) => Some(source),
            Self::PhysicalParams(source) => Some(source),
            Self::Wind(source) => Some(source),
            Self::Schedule(source) => Some(source),
            Self::Cfl(source) => Some(source),
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

    /// The [`Basin`] this section describes, in metres.
    ///
    /// # Errors
    /// Whatever [`BasinSection::bounds`] objected to.
    pub fn build(&self) -> Result<Basin, ScenarioError> {
        Ok(self.bounds()?.basin())
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
}

/// A scenario as it is written in a file: four sections, no invariants.
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
    /// The CFL check comes last, because it needs the wave speed the physics
    /// section implies and the spacing the basin section implies; it refuses
    /// an unstable timestep rather than shortening it
    /// (CODING_STANDARDS.md § *No silent clamping*).
    ///
    /// # Errors
    /// A [`ScenarioError`] naming the section that was wrong and, through the
    /// error it wraps, the value and the bound.
    pub fn build(&self) -> Result<Scenario, ScenarioError> {
        let basin = self.basin.build()?;
        let physical_params = self.physics.build()?;
        let schedule = self.run.build()?;
        let winds = self
            .wind
            .iter()
            .map(WindSection::build)
            .collect::<Result<Vec<_>, _>>()?;

        let wave_speed = WaveSpeed::new(physical_params.kelvin_wave_speed_m_per_s())?;
        check_timestep(schedule.dt_s(), basin.spacing(), wave_speed)?;

        Ok(Scenario {
            basin,
            physical_params,
            winds,
            output_schedule: schedule,
        })
    }
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
    basin: Basin,
    physical_params: PhysicalParams,
    winds: Vec<ScenarioWind>,
    output_schedule: OutputSchedule,
}

impl Scenario {
    /// Read and validate a scenario from a TOML file on disk.
    ///
    /// # Errors
    /// [`ScenarioError::Unreadable`] naming `path` if the file cannot be read,
    /// otherwise whatever [`ScenarioConfig::from_toml`] or
    /// [`ScenarioConfig::build`] objected to.
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
}
