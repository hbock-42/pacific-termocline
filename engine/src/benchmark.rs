//! The workloads the `cargo bench` suite measures, and the counts it reports
//! throughput against.
//!
//! This module holds no timing code. Criterion owns the measurement — how many
//! samples, how long a warm-up, what statistics — and the benchmarks in
//! `benches/` own the calls being timed. What lives here is the *input*: which
//! grids the suite runs on, which scenario it runs, how long the short run is,
//! and the analytic state the right-hand side is evaluated over. It is a
//! module of the library rather than a `benches/common/` helper because a
//! benchmark input is exactly the kind of thing that quietly drifts —
//! `tests/benchmark_workloads.rs` asserts on everything below, and it can only
//! do that if the definitions are importable.
//!
//! # What the suite measures, and what it does not
//!
//! Two workloads, one per benchmark in the ticket:
//!
//! - **The right-hand side** — one [`ShallowWaterRhs::evaluate`] over a whole
//!   basin: the pressure gradient, the surface stress, the Rayleigh damping
//!   and the continuity divergence of Epic 02, and nothing else. This is the
//!   innermost thing a run does, called four times per step by RK4, and the
//!   figure it reports is **grid cells per second** — [`grid_cells`] over the
//!   measured duration.
//! - **A short scenario run** — [`run_scenario`] end to end at
//!   [`SHORT_RUN_STEPS`] steps: scenario build, solver construction, the time
//!   loop with its forcing re-sampled at every RK4 stage, and the run
//!   directory written. The figure it reports is **timesteps per second** —
//!   [`timesteps`] over the measured duration. Grid cells per second for a run
//!   is the product of the two: `steps/s × grid_cells`.
//!
//!   Its output cadence is set so the run writes exactly two frames, the
//!   initial state and the final one. The measurement therefore *includes*
//!   the filesystem, because a run does, but it includes a constant two frames
//!   of it rather than a number that grows with the run length: the variable
//!   part of the timing is the time loop. A benchmark whose variance came from
//!   the page cache would not be able to see the changes Epic 10's later
//!   tickets make.
//!
//! Neither workload is a *scaling* study of the timestep: `dt` is the control
//! scenario's one hour at every resolution, so that the two grids differ in
//! the grid alone. That is deliberate — the CFL bound would let the 1.0° grid
//! take a longer step, and a suite whose coarse case also stepped differently
//! could not attribute a change to either.
//!
//! # Why the control scenario, and why only two knobs
//!
//! Every workload is [`engine/scenarios/steady-trades.toml`] — the project's
//! control scenario, embedded at compile time — with two fields overridden:
//! the basin resolution and the run length. Its `g'`, `H`, `r`, `β`, `ρ₀`,
//! timestep and wind are the ones a real run uses, so a figure from this suite
//! is a figure about the simulation rather than about a benchmark-only
//! configuration. `tests/benchmark_workloads.rs` holds the two to that.
//!
//! # Reproducibility
//!
//! Nothing here is random, timed, or read from the environment. The initial
//! state is a closed-form analytic field, the scenario is compiled in, and the
//! run is the engine's own deterministic one (CODING_STANDARDS.md
//! § *Correctness and failure*), so two runs of a workload write byte-identical
//! output and two evaluations of the right-hand side produce bit-identical
//! tendencies. That is what makes a difference between two measurements
//! attributable to the code rather than to the input.
//!
//! [`engine/scenarios/steady-trades.toml`]: ../../scenarios/steady-trades.toml
//! [`grid_cells`]: BenchmarkWorkload::grid_cells
//! [`timesteps`]: BenchmarkWorkload::timesteps

use std::fs;
use std::path::{Path, PathBuf};
use std::process;

use termocline_grid::{Field2D, Grid, Staggering, H_STAGGERING, U_STAGGERING};

use crate::basin::Basin;
use crate::forcing::WindStressField;
use crate::run::{run_scenario, RunError, RunReport};
use crate::scenario::{Scenario, ScenarioConfig};
use crate::shallow_water::ShallowWaterRhs;
use crate::state::OceanState;

/// The control scenario, embedded so a benchmark reads no files and depends on
/// no working directory.
const CONTROL_SCENARIO_TOML: &str = include_str!("../scenarios/steady-trades.toml");

/// Steps the short-run benchmark takes: ten days of model time at the control
/// scenario's one-hour timestep.
///
/// Long enough that the time loop, not the construction in front of it,
/// dominates the measurement; short enough that criterion can sample it
/// repeatedly. Ten days is also long enough for a Kelvin wave to cross about a
/// seventh of the basin at `c = 3 m/s`, so the run is doing recognisable
/// physics rather than sitting in its first transient.
pub const SHORT_RUN_STEPS: u64 = 240;

/// Amplitude `A` of the benchmark state's thermocline depth anomaly, in
/// metres.
///
/// The same 10 m pulse the Kelvin wave validation of T-07.1 launches: a
/// realistic first-baroclinic anomaly against a mean thermocline depth of
/// 150 m, and small enough that the linear equations it is fed to are the
/// right ones (`docs/planning/01-scientific-model.md`).
pub const BENCHMARK_ANOMALY_AMPLITUDE_M: f64 = 10.0;

/// Zonal width `W` of the benchmark state's anomaly, in metres.
///
/// A pulse of 1 500 km, as in T-07.1: wide enough to be resolved at the
/// coarsest grid in the suite (1.0° ≈ 111 km cells, so ~14 cells across) and
/// narrow enough to sit well inside a basin 17 800 km across without touching
/// either meridional wall.
pub const BENCHMARK_ANOMALY_WIDTH_M: f64 = 1.5e6;

/// One benchmark input: the control scenario at a stated resolution, cut down
/// to a short run.
///
/// A `BenchmarkWorkload` is a description, not a fixture — it holds two
/// numbers, and every method below builds what it describes on demand. That is
/// what lets a benchmark construct fresh inputs outside criterion's timing
/// loop and a test construct the same ones without a benchmark.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BenchmarkWorkload {
    /// Cell size, in degrees, on both axes.
    resolution_deg: f64,
    /// Steps the short run takes.
    total_steps: u64,
}

/// The resolutions the suite measures at: the control scenario's own 0.5°
/// grid, and the 1.0° grid at a quarter of its cells.
///
/// Two, because the ticket asks for "a couple of representative grid
/// resolutions" and because two points are what turn a duration into a scaling
/// law — the question the profiling and parallelisation tickets behind this
/// one need answered is how cost grows with the basin, and one grid cannot say.
/// They are a factor of two apart in each direction, so a cost linear in the
/// cell count shows up as a factor of four.
///
/// Finer than 0.5° is deliberately absent: a 0.25° run is four times the work
/// again, which is a benchmark nobody runs on a laptop and therefore a
/// benchmark nobody runs.
pub const BENCHMARK_WORKLOADS: [BenchmarkWorkload; 2] = [
    BenchmarkWorkload {
        resolution_deg: 1.0,
        total_steps: SHORT_RUN_STEPS,
    },
    BenchmarkWorkload {
        resolution_deg: 0.5,
        total_steps: SHORT_RUN_STEPS,
    },
];

impl BenchmarkWorkload {
    /// Cell size of this workload's basin, in degrees.
    #[must_use]
    pub const fn resolution_deg(&self) -> f64 {
        self.resolution_deg
    }

    /// Steps this workload's short run takes — the count the timesteps-per-
    /// second figure is measured against.
    #[must_use]
    pub const fn timesteps(&self) -> u64 {
        self.total_steps
    }

    /// This workload's basin, in cells.
    #[must_use]
    pub fn grid(&self) -> Grid {
        self.scenario().basin().grid()
    }

    /// Cells in this workload's basin — the count the grid-cells-per-second
    /// figure is measured against.
    ///
    /// The cell *centres*, `nx · ny`, where `h` and the divergence live. The
    /// staggered velocity fields carry one extra column and one extra row on
    /// their own faces, which is the C-grid's bookkeeping rather than a second
    /// grid, and counting them would report a throughput no reader could
    /// compare against a resolution.
    #[must_use]
    pub fn grid_cells(&self) -> u64 {
        let grid = self.grid();
        (grid.nx() * grid.ny()) as u64
    }

    /// How this workload names itself in a criterion report: its grid shape,
    /// as `320x100`.
    ///
    /// Derived from the grid rather than stored beside it, so a benchmark id
    /// cannot come to say one thing while the workload runs another.
    #[must_use]
    pub fn label(&self) -> String {
        let grid = self.grid();
        format!("{}x{}", grid.nx(), grid.ny())
    }

    /// The scenario this workload runs: the control scenario at this
    /// resolution and length.
    ///
    /// The output cadence is the whole run, so the run writes exactly two
    /// frames — the initial state and the final one — and the filesystem is a
    /// constant rather than a term that grows with the measurement.
    ///
    /// # Panics
    /// If the embedded control scenario is not a scenario the engine will run
    /// at this resolution. Both are compiled in, so that is a statement about
    /// this file rather than about anything a user did — which is what panics
    /// are for (CODING_STANDARDS.md § *Correctness and failure*).
    #[must_use]
    pub fn scenario(&self) -> Scenario {
        let mut config = ScenarioConfig::from_toml(CONTROL_SCENARIO_TOML)
            .expect("the control scenario is embedded from this repository and parses");
        config.basin.resolution_deg = self.resolution_deg;
        config.run.total_steps = self.total_steps;
        config.run.output_every_n_steps = self.total_steps;
        config
            .build()
            .expect("the control scenario is valid at every resolution the suite measures")
    }

    /// Run this workload's scenario into `directory` — the call the short-run
    /// benchmark times.
    ///
    /// # Errors
    /// The errors of [`run_scenario`]: in this suite, only a run directory
    /// that could not be written, since the scenario itself is compiled in and
    /// validated by [`BenchmarkWorkload::scenario`].
    pub fn run_into(&self, directory: &Path) -> Result<RunReport, RunError> {
        run_scenario(&self.scenario(), "benchmark", directory)
    }

    /// A right-hand-side evaluator for this workload's basin — the object the
    /// right-hand-side benchmark calls [`ShallowWaterRhs::evaluate`] on.
    ///
    /// Built once, outside the timing loop, exactly as a run builds it once
    /// per run (CODING_STANDARDS.md § *Performance*): what is measured is an
    /// evaluation, not an allocation.
    #[must_use]
    pub fn rhs_evaluator(&self) -> ShallowWaterRhs {
        let scenario = self.scenario();
        let basin = scenario.basin();
        ShallowWaterRhs::new(basin.grid(), basin.spacing(), scenario.physical_params())
    }

    /// The state the right-hand side is evaluated over: an equatorially
    /// trapped Kelvin structure.
    ///
    /// ```text
    /// h(x, y) = A · exp(−((x − x_c)/W)²) · exp(−η²/2)      η = y/Le
    /// u(x, y) = (c/H) · h(x, y)                            v ≡ 0
    /// ```
    ///
    /// with `Le = √(c/β)` the equatorial deformation radius and `c = √(g'H)`
    /// the Kelvin wave speed (CONTEXT.md, *Kelvin wave*). `u = (c/H)·h` is the
    /// Kelvin balance, which is what makes this a state a run could actually
    /// be in rather than an arbitrary field.
    ///
    /// It matters that it is not an ocean at rest. The evaluator has no
    /// data-dependent branches, so the values do not change *which* arithmetic
    /// is done — but a field of exact zeros is not what a machine spends its
    /// time on in a real run, and a benchmark should not be the one place the
    /// hardware sees an easier problem than the simulation does. Each field is
    /// evaluated at its own C-grid position, so `h` at cell centres and `u` on
    /// east–west faces are the same analytic function sampled where each
    /// variable lives (CODING_STANDARDS.md § *Scope guards*).
    #[must_use]
    pub fn benchmark_state(&self) -> OceanState {
        let scenario = self.scenario();
        let basin = scenario.basin();
        let params = scenario.physical_params();
        let mut state = OceanState::at_rest(basin.grid());

        let wave_speed_m_per_s = params.kelvin_wave_speed_m_per_s();
        let deformation_radius_m = (wave_speed_m_per_s / params.beta_per_m_per_s()).sqrt();
        let centre_x_m = basin.western_edge_x_m() + basin.zonal_extent_m() / 2.0;
        let anomaly_m = |x_m: f64, y_m: f64| {
            let zonal = ((x_m - centre_x_m) / BENCHMARK_ANOMALY_WIDTH_M).powi(2);
            let meridional = (y_m / deformation_radius_m).powi(2) / 2.0;
            BENCHMARK_ANOMALY_AMPLITUDE_M * (-zonal - meridional).exp()
        };

        fill_field(state.h_mut(), basin, H_STAGGERING, anomaly_m);
        let current_scale_per_s = wave_speed_m_per_s / params.mean_thermocline_depth_m();
        fill_field(state.u_mut(), basin, U_STAGGERING, |x_m, y_m| {
            current_scale_per_s * anomaly_m(x_m, y_m)
        });
        state
    }

    /// The surface stress the right-hand side is evaluated under: this
    /// workload's own wind, sampled at `t = 0`.
    ///
    /// The scenario's forcing rather than a calm field, because the stress
    /// term is part of what the benchmark is meant to cover.
    #[must_use]
    pub fn wind_stress(&self) -> WindStressField {
        let scenario = self.scenario();
        WindStressField::sampled(scenario.basin(), &scenario.wind(), 0.0)
    }

    /// A tendency buffer for this workload's basin — the `&mut` argument of
    /// [`ShallowWaterRhs::evaluate`].
    ///
    /// Every point of it is written by an evaluation, so one buffer serves
    /// every iteration of the benchmark and no allocation lands inside the
    /// timing loop.
    #[must_use]
    pub fn tendency_buffer(&self) -> OceanState {
        OceanState::at_rest(self.grid())
    }
}

/// A directory for one workload's output, removed when it is dropped.
///
/// The short-run benchmark writes a real run, so it needs somewhere to write
/// it, and so does every test that runs a workload. It lives here rather than
/// beside `tests/common::ScratchDir` for the same reason the workloads do: a
/// `benches/` target cannot reach into a test target's modules, and one
/// definition both can import is better than two that drift.
///
/// Two workloads must not share one. A frame file left by a larger grid would
/// be truncated by the next run rather than replaced, and the run would be
/// measured writing over it.
#[derive(Debug)]
pub struct BenchmarkOutputDir {
    path: PathBuf,
}

impl BenchmarkOutputDir {
    /// A fresh empty directory, labelled `label` and carrying the process id
    /// so a leftover one says what wrote it and two processes cannot collide.
    ///
    /// # Panics
    /// If the system temp directory cannot be written, which is a broken
    /// environment rather than a failed measurement.
    #[must_use]
    pub fn new(label: &str) -> Self {
        let path = std::env::temp_dir().join(format!("termocline-bench-{label}-{}", process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).expect("the system temp directory is writable");
        Self { path }
    }

    /// Where the directory is.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for BenchmarkOutputDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

/// Write `value_at(x, y)` into every point of `field`, at the C-grid positions
/// `staggering` puts them.
fn fill_field(
    field: &mut Field2D<f64>,
    basin: Basin,
    staggering: Staggering,
    value_at: impl Fn(f64, f64) -> f64,
) {
    let (nx, ny) = (field.nx(), field.ny());
    for j in 0..ny {
        let y_m = basin.y_of_row_m(staggering, j);
        for i in 0..nx {
            let x_m = basin.x_of_column_m(staggering, i);
            *field
                .get_mut(i, j)
                .expect("i and j are inside the field's own shape") = value_at(x_m, y_m);
        }
    }
}
