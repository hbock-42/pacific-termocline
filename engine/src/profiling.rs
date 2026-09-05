//! T-10.2 — where a timestep's time goes, measured rather than guessed.
//!
//! [`docs/benchmarks.md`] says how fast the engine is; this module says *what
//! it is spending that on*. Epic 10's remaining tickets — a rayon-parallel
//! inner loop (T-10.3) and an `f32` field layout (T-10.4) — must cite a
//! profile before they touch anything (CODING_STANDARDS.md § *Performance*),
//! and this is the instrument that produces it. The findings it produced are
//! written up in [`docs/performance-notes.md`].
//!
//! Nothing here is part of a run. The solver, the right-hand side and the
//! integrator are untouched by this ticket: a profiler that changed the code
//! it measures would be measuring itself.
//!
//! # Two resolutions, and why both
//!
//! A profile of a step can be read at two depths, and neither is sufficient
//! alone.
//!
//! - [`StepProfiler`] splits a whole step into the five [`StepPhase`]s a step
//!   is made of, using the engine's *own* [`ShallowWaterRhs`],
//!   [`CoriolisTerm`], [`NoNormalFlow`] and [`Rk4`]. It is exact and it adds
//!   up: the four instrumented phases are timed directly and
//!   [`StepPhase::StageAlgebra`] is the residual, so the five shares sum to
//!   one by construction and no work can hide between them.
//!   `tests/step_profile.rs` pins the profiled step to
//!   [`Solver::step_forced_by`](crate::Solver::step_forced_by) bit for bit, so
//!   what is being timed is the step a run actually takes. It holds its
//!   forcing across steps, as [`run_scenario`](crate::run_scenario) does, so
//!   the wind phase is charged what a run pays rather than what a solver
//!   handed a fresh wind each call pays.
//!
//! - [`TermProfiler`] splits the two evaluator phases further, into the
//!   fourteen [`RhsTerm`] array kernels the right-hand side and the Coriolis
//!   term are built from. This is the level the optimisation tickets need —
//!   "the right-hand side is 40% of a step" does not tell you which loop to
//!   parallelise. It calls the same kernels the evaluators call, in the same
//!   order, over the same buffers, and `tests/step_profile.rs` asserts that
//!   its tendency is bit-identical to the evaluators' own.
//!
//! # What a timing decomposition cannot see
//!
//! Two caveats, stated here because a profile read without them is a profile
//! read wrong.
//!
//! - **The clock is not free.** A step reads [`Instant::now`]
//!   [`CLOCK_READS_PER_STEP`] times, so a decomposition charges the
//!   instrumented phases that many clock reads that a real step does not
//!   pay, and credits the residual with the same amount. On the machine the
//!   note was taken on that is well under a microsecond against a step of
//!   milliseconds — [`clock_read_cost`] measures it, and the example prints
//!   it beside the table so a reader can check the ratio for themselves rather
//!   than trusting this sentence.
//!
//! - **Splitting a computation can change it.** Timing a kernel on its own
//!   puts a barrier where the optimiser might otherwise have fused two loops
//!   or kept a value in a register across them, so [`TermProfiler`]'s terms
//!   may not sum to [`StepPhase::ShallowWaterTerms`] plus
//!   [`StepPhase::Coriolis`]. That gap is a real quantity — it is what
//!   inlining was worth — and it is why the note reads the term table as
//!   *proportions within the evaluators* and the phase table as absolute
//!   shares of a step, rather than pretending the two are one table.
//!
//! Both are why the note behind this module also carries a sampled profile
//! from an external sampler, which perturbs nothing: two instruments that
//! disagree are a finding, and two that agree are a conclusion.
//!
//! [`docs/benchmarks.md`]: ../../docs/benchmarks.md
//! [`docs/performance-notes.md`]: ../../docs/performance-notes.md

use std::time::{Duration, Instant};

use termocline_grid::{Field2D, Grid, H_STAGGERING, U_STAGGERING, V_STAGGERING};
use termocline_numerics::{check_timestep, CGridOperators, Spacing, WaveSpeed};

use crate::boundary::NoNormalFlow;
use crate::coriolis::{accumulate_rows, BetaPlane, CoriolisTerm};
use crate::forcing::{CompositeWind, WindForcing, WindStressField};
use crate::integrator::Rk4;
use crate::params::PhysicalParams;
use crate::scenario::Scenario;
use crate::shallow_water::{
    subtract_damping, turn_gradient_into_acceleration, write_continuity, ShallowWaterRhs,
};
use crate::solver::{check_rotation_timestep, SolverError};
use crate::state::OceanState;

/// One of the five pieces a timestep is made of.
///
/// The split is by *what owns the work*, not by what it costs: the wind is
/// sampled by the forcing, the shallow-water terms by [`ShallowWaterRhs`],
/// rotation by [`CoriolisTerm`], the walls by [`NoNormalFlow`], and everything
/// left over is [`Rk4`] moving states around. That is the granularity at which
/// an optimisation ticket picks its target.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StepPhase {
    /// Asking the forcing for the stress at the stage's instant, and sampling
    /// it onto the C-grid when the field in hand is not already that field.
    ///
    /// Once per RK4 stage until T-10.5; since then, once per instant the
    /// forcing has not already been asked about, which for the control
    /// scenario is once per run (`docs/performance-notes.md`).
    WindStressSampling,
    /// [`ShallowWaterRhs::evaluate`]: pressure gradient, surface stress,
    /// Rayleigh damping and the continuity divergence.
    ShallowWaterTerms,
    /// [`CoriolisTerm::add_to_tendency`]: the two four-point interpolations
    /// and the `±f·` accumulations on top of them.
    Coriolis,
    /// [`NoNormalFlow`] on the incoming state and on each stage's tendency.
    BoundaryCondition,
    /// What [`Rk4`] itself does: assigning the stage state and accumulating
    /// the four weighted stages into the result.
    ///
    /// Measured as the residual — the step's own duration less the four phases
    /// above — so that the five shares sum to one and unattributed work cannot
    /// vanish.
    StageAlgebra,
}

/// Every [`StepPhase`], in the order a step performs them.
pub const STEP_PHASES: [StepPhase; 5] = [
    StepPhase::WindStressSampling,
    StepPhase::ShallowWaterTerms,
    StepPhase::Coriolis,
    StepPhase::BoundaryCondition,
    StepPhase::StageAlgebra,
];

/// The four phases a [`StepProfiler`] times directly; the fifth is their
/// residual.
const TIMED_PHASES: usize = STEP_PHASES.len() - 1;

/// How many times [`StepProfiler::step`] reads the clock, per step.
///
/// The instrument's whole footprint, and the number a reader multiplies by
/// [`clock_read_cost`] to decide whether a phase table is worth reading. It is
/// stated here, beside [`StepProfiler::step`], and *derived* rather than
/// counted by eye at each place that quotes it — a phase added to the step
/// changes this number, and a hand-copied one would go quietly stale.
pub const CLOCK_READS_PER_STEP: u32 = STEP_CLOCK_READS + RK4_STAGES * STAGE_CLOCK_READS;

/// Stages of the RK4 tableau, each of which reads the clock
/// [`STAGE_CLOCK_READS`] times.
const RK4_STAGES: u32 = 4;

/// Clock reads inside one RK4 stage: one mark before the wind is sampled, and
/// one after each of the stage's four phases.
const STAGE_CLOCK_READS: u32 = 1 + TIMED_PHASES as u32;

/// Clock reads a step performs outside its stages: the pair around the step
/// itself, and the pair around the boundary condition applied to the incoming
/// state.
const STEP_CLOCK_READS: u32 = 4;

impl StepPhase {
    /// How this phase names itself in a report.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::WindStressSampling => "wind stress sampling",
            Self::ShallowWaterTerms => "shallow-water terms",
            Self::Coriolis => "coriolis",
            Self::BoundaryCondition => "boundary condition",
            Self::StageAlgebra => "rk4 stage algebra",
        }
    }

    /// Position of this phase in [`STEP_PHASES`], which is also its slot in a
    /// [`StepProfile`].
    const fn index(self) -> usize {
        match self {
            Self::WindStressSampling => 0,
            Self::ShallowWaterTerms => 1,
            Self::Coriolis => 2,
            Self::BoundaryCondition => 3,
            Self::StageAlgebra => 4,
        }
    }
}

/// What a run of [`StepProfiler`] measured: the wall time of some number of
/// whole steps, and how it divided among the [`STEP_PHASES`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepProfile {
    /// Steps the profile covers.
    steps: u64,
    /// Wall time of those steps, end to end.
    total: Duration,
    /// Time in each phase, indexed by [`StepPhase::index`].
    phases: [Duration; STEP_PHASES.len()],
    /// How many times each phase was charged, indexed the same way.
    charges: [u32; STEP_PHASES.len()],
}

impl StepProfile {
    /// Steps this profile covers.
    #[must_use]
    pub const fn steps(&self) -> u64 {
        self.steps
    }

    /// Wall time of those steps, end to end.
    #[must_use]
    pub const fn total(&self) -> Duration {
        self.total
    }

    /// Mean wall time of one step.
    #[must_use]
    pub fn per_step(&self) -> Duration {
        mean_over(self.total, self.steps)
    }

    /// Time spent in `phase` across the whole profile.
    #[must_use]
    pub const fn phase(&self, phase: StepPhase) -> Duration {
        self.phases[phase.index()]
    }

    /// Mean time one step spent in `phase`.
    #[must_use]
    pub fn per_step_in(&self, phase: StepPhase) -> Duration {
        mean_over(self.phase(phase), self.steps)
    }

    /// How many times `phase` was charged over the whole profile.
    ///
    /// The count rather than the duration, so that a test can check the
    /// decomposition is complete without asserting on a clock: a step charges
    /// each of the four evaluator phases once per RK4 stage, the boundary
    /// condition once more for the state it brings in, and
    /// [`StepPhase::StageAlgebra`] once, as the residual of the step. A phase
    /// the profiler forgot to charge shows up here as a zero on any machine,
    /// where a zero *duration* would only show up on a fast one.
    #[must_use]
    pub const fn charges(&self, phase: StepPhase) -> u32 {
        self.charges[phase.index()]
    }

    /// The fraction of a step spent in `phase`, between 0 and 1.
    ///
    /// The five shares sum to one:
    /// [`StepPhase::StageAlgebra`] is the residual of the other four, so a
    /// share is a division of the measured step rather than of a subtotal.
    #[must_use]
    pub fn share(&self, phase: StepPhase) -> f64 {
        let total = self.total.as_secs_f64();
        if total == 0.0 {
            return 0.0;
        }
        self.phase(phase).as_secs_f64() / total
    }
}

/// A step of one scenario, taken with every phase timed.
///
/// It is a [`Solver`](crate::Solver) with a stopwatch between the pieces: the
/// same [`ShallowWaterRhs`], the same [`CoriolisTerm`], the same
/// [`NoNormalFlow`] and the same [`Rk4`], composed in the same order, with the
/// same buffers allocated once at construction
/// (CODING_STANDARDS.md § *Performance*). It is written out here rather than
/// wrapped around `Solver` because the phase boundaries are *inside* the
/// closure `Solver` hands to the integrator, and a caller cannot reach in
/// there. `tests/step_profile.rs` is what holds the two to being the same
/// step.
#[derive(Debug)]
pub struct StepProfiler {
    /// The scenario's wind and the field it is sampled into, re-sampled at
    /// each stage exactly as a run does — which since T-10.5 means *when the
    /// stage asks about an instant the field is not already the field of*.
    forcing: WindForcing<CompositeWind>,
    /// Length of one step, in seconds.
    dt_s: f64,
    /// The pressure-gradient, continuity, surface-stress and damping terms.
    rhs: ShallowWaterRhs,
    /// The beta-plane rotation terms.
    coriolis: CoriolisTerm,
    /// The integrator and its stage buffers.
    integrator: Rk4<OceanState>,
    /// The state being advanced.
    state: OceanState,
    /// Model time of `state`, in seconds.
    t_s: f64,
    /// Time accumulated in each directly timed phase.
    timed: [Duration; TIMED_PHASES],
    /// How many times each directly timed phase has been charged.
    charged: [u32; TIMED_PHASES],
    /// Wall time of the whole steps taken so far.
    total: Duration,
    /// Steps taken so far.
    steps: u64,
}

impl StepProfiler {
    /// A profiler for `scenario`, starting from the ocean at rest at `t = 0`.
    ///
    /// # Errors
    /// The errors [`Solver::new`](crate::Solver::new) returns, and for the
    /// same reasons: a timestep past the gravity-wave CFL bound or past the
    /// basin's rotation bound is a scenario that cannot be run, so it is not
    /// one that can be profiled either. Both bounds are checked here rather
    /// than inherited, exactly as the scenario loader checks them
    /// (`solver.rs`, *Two bounds on the timestep, not one*).
    pub fn new(scenario: &Scenario) -> Result<Self, SolverError> {
        let basin = scenario.basin();
        let grid = basin.grid();
        let spacing = basin.spacing();
        let params = scenario.physical_params();
        let dt_s = scenario.output_schedule().dt_s();

        let wave_speed = WaveSpeed::new(params.kelvin_wave_speed_m_per_s())
            .expect("physical parameters are validated positive, so `√(g'·H)` is too");
        check_timestep(dt_s, spacing, wave_speed)?;
        let plane = BetaPlane::of_basin(params, basin);
        check_rotation_timestep(dt_s, grid, plane)?;

        Ok(Self {
            forcing: WindForcing::new(basin, scenario.wind()),
            dt_s,
            rhs: ShallowWaterRhs::new(grid, spacing, params),
            coriolis: CoriolisTerm::new(grid, spacing, plane),
            integrator: Rk4::new(&OceanState::at_rest(grid)),
            state: OceanState::at_rest(grid),
            t_s: 0.0,
            timed: [Duration::ZERO; TIMED_PHASES],
            charged: [0; TIMED_PHASES],
            total: Duration::ZERO,
            steps: 0,
        })
    }

    /// The state as far as the profiler has advanced it.
    ///
    /// The reason a profile is trustworthy: a test can compare this against
    /// the state a [`Solver`](crate::Solver) reaches from the same scenario
    /// and see that the instrumented step is the real one.
    #[must_use]
    pub const fn state(&self) -> &OceanState {
        &self.state
    }

    /// Take one step, charging each phase the time it took.
    pub fn step(&mut self) {
        let Self {
            forcing,
            dt_s,
            rhs,
            coriolis,
            integrator,
            state,
            t_s,
            timed,
            charged,
            total,
            steps,
        } = self;

        // Charged to the whole step, so the residual phase is the step less
        // everything named — including this clock read.
        let step_started = Instant::now();
        let mut this_step = [Duration::ZERO; TIMED_PHASES];
        let mut charges = [0_u32; TIMED_PHASES];

        let boundary_started = Instant::now();
        NoNormalFlow::apply_to_state(state);
        this_step[StepPhase::BoundaryCondition.index()] += boundary_started.elapsed();
        charges[StepPhase::BoundaryCondition.index()] += 1;

        integrator.step(
            state,
            *t_s,
            *dt_s,
            &mut |now: &OceanState, stage_t_s: f64, tendency: &mut OceanState| {
                // The order is the solver's own: the shallow-water evaluator
                // writes every point, the Coriolis term adds to it, and the
                // boundary condition has the last word.
                let sampled = Instant::now();
                let stage_stress = forcing.at(stage_t_s);
                let evaluated = Instant::now();
                rhs.evaluate(now, stage_stress, tendency);
                let rotated = Instant::now();
                coriolis.add_to_tendency(now, tendency);
                let bounded = Instant::now();
                NoNormalFlow::apply_to_tendency(tendency);
                let finished = Instant::now();

                this_step[StepPhase::WindStressSampling.index()] += evaluated - sampled;
                this_step[StepPhase::ShallowWaterTerms.index()] += rotated - evaluated;
                this_step[StepPhase::Coriolis.index()] += bounded - rotated;
                this_step[StepPhase::BoundaryCondition.index()] += finished - bounded;
                for charge in &mut charges {
                    *charge += 1;
                }
            },
        );

        *total += step_started.elapsed();
        for (accumulated, this) in timed.iter_mut().zip(this_step) {
            *accumulated += this;
        }
        for (accumulated, this) in charged.iter_mut().zip(charges) {
            *accumulated += this;
        }
        *t_s += *dt_s;
        *steps += 1;
    }

    /// Take `steps` steps and report how their time divided.
    ///
    /// The profile covers every step this profiler has ever taken, so a caller
    /// that wants warm-up excluded takes it on a profiler it then throws away.
    pub fn profile(&mut self, steps: u64) -> StepProfile {
        for _ in 0..steps {
            self.step();
        }
        self.profile_so_far()
    }

    /// What has been measured so far, without taking another step.
    #[must_use]
    pub fn profile_so_far(&self) -> StepProfile {
        let mut phases = [Duration::ZERO; STEP_PHASES.len()];
        let mut charges = [0_u32; STEP_PHASES.len()];
        let mut named = Duration::ZERO;
        for (phase, measured) in phases.iter_mut().zip(self.timed) {
            *phase = measured;
            named += measured;
        }
        charges[..TIMED_PHASES].copy_from_slice(&self.charged);
        // The residual is charged once per step, being the one thing a step
        // computes rather than measures.
        charges[StepPhase::StageAlgebra.index()] = u32::try_from(self.steps).unwrap_or(u32::MAX);
        // Saturating: the residual is a difference of two clocks, and a step
        // whose named phases measured longer than the step around them is a
        // measurement at the resolution of the clock rather than negative work.
        phases[StepPhase::StageAlgebra.index()] = self.total.saturating_sub(named);
        StepProfile {
            steps: self.steps,
            total: self.total,
            phases,
            charges,
        }
    }
}

/// One of the fourteen array kernels the right-hand side and the Coriolis term
/// are built from.
///
/// In the order they are performed. Each is a pass over one or two basin-sized
/// fields, which is why this is the level a parallelisation or a narrower
/// float type acts at.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RhsTerm {
    /// `∂h/∂x` from cell centres onto east/west faces.
    ZonalPressureGradient,
    /// `∂h/∂y` from cell centres onto north/south faces.
    MeridionalPressureGradient,
    /// `−g'·∂h/∂x + τx/(ρ₀·H)` in place on the zonal tendency.
    ZonalPressureAndStress,
    /// `−g'·∂h/∂y + τy/(ρ₀·H)` in place on the meridional tendency.
    MeridionalPressureAndStress,
    /// `−r·u` on the zonal tendency.
    ZonalDamping,
    /// `−r·v` on the meridional tendency.
    MeridionalDamping,
    /// `∂u/∂x` from east/west faces onto cell centres.
    ZonalDivergence,
    /// `∂v/∂y` from north/south faces onto cell centres.
    MeridionalDivergence,
    /// `−H·(∂u/∂x + ∂v/∂y)` onto the thermocline tendency.
    Continuity,
    /// `−r·h` on the thermocline tendency.
    ThermoclineDamping,
    /// `v` interpolated onto the east/west faces.
    MeridionalVelocityOntoZonalFaces,
    /// `u` interpolated onto the north/south faces.
    ZonalVelocityOntoMeridionalFaces,
    /// `+f·v` accumulated onto the zonal tendency.
    ZonalRotation,
    /// `−f·u` accumulated onto the meridional tendency.
    MeridionalRotation,
}

/// Every [`RhsTerm`], in the order an evaluation performs them.
pub const RHS_TERMS: [RhsTerm; 14] = [
    RhsTerm::ZonalPressureGradient,
    RhsTerm::MeridionalPressureGradient,
    RhsTerm::ZonalPressureAndStress,
    RhsTerm::MeridionalPressureAndStress,
    RhsTerm::ZonalDamping,
    RhsTerm::MeridionalDamping,
    RhsTerm::ZonalDivergence,
    RhsTerm::MeridionalDivergence,
    RhsTerm::Continuity,
    RhsTerm::ThermoclineDamping,
    RhsTerm::MeridionalVelocityOntoZonalFaces,
    RhsTerm::ZonalVelocityOntoMeridionalFaces,
    RhsTerm::ZonalRotation,
    RhsTerm::MeridionalRotation,
];

impl RhsTerm {
    /// How this term names itself in a report.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::ZonalPressureGradient => "d(h)/dx  centre -> u face",
            Self::MeridionalPressureGradient => "d(h)/dy  centre -> v face",
            Self::ZonalPressureAndStress => "-g'.dh/dx + taux/(rho.H)",
            Self::MeridionalPressureAndStress => "-g'.dh/dy + tauy/(rho.H)",
            Self::ZonalDamping => "-r.u",
            Self::MeridionalDamping => "-r.v",
            Self::ZonalDivergence => "d(u)/dx  u face -> centre",
            Self::MeridionalDivergence => "d(v)/dy  v face -> centre",
            Self::Continuity => "-H.(du/dx + dv/dy)",
            Self::ThermoclineDamping => "-r.h",
            Self::MeridionalVelocityOntoZonalFaces => "v -> u faces (4-point)",
            Self::ZonalVelocityOntoMeridionalFaces => "u -> v faces (4-point)",
            Self::ZonalRotation => "+f.v",
            Self::MeridionalRotation => "-f.u",
        }
    }

    /// Which [`StepPhase`] this term is part of, so a term table can be read
    /// against the phase table above it.
    #[must_use]
    pub const fn phase(self) -> StepPhase {
        match self {
            Self::MeridionalVelocityOntoZonalFaces
            | Self::ZonalVelocityOntoMeridionalFaces
            | Self::ZonalRotation
            | Self::MeridionalRotation => StepPhase::Coriolis,
            _ => StepPhase::ShallowWaterTerms,
        }
    }

    /// Position of this term in [`RHS_TERMS`], which is also its slot in a
    /// [`TermProfile`].
    const fn index(self) -> usize {
        match self {
            Self::ZonalPressureGradient => 0,
            Self::MeridionalPressureGradient => 1,
            Self::ZonalPressureAndStress => 2,
            Self::MeridionalPressureAndStress => 3,
            Self::ZonalDamping => 4,
            Self::MeridionalDamping => 5,
            Self::ZonalDivergence => 6,
            Self::MeridionalDivergence => 7,
            Self::Continuity => 8,
            Self::ThermoclineDamping => 9,
            Self::MeridionalVelocityOntoZonalFaces => 10,
            Self::ZonalVelocityOntoMeridionalFaces => 11,
            Self::ZonalRotation => 12,
            Self::MeridionalRotation => 13,
        }
    }
}

/// What a run of [`TermProfiler`] measured: how the time inside the two
/// evaluators divided among the [`RHS_TERMS`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TermProfile {
    /// Evaluations the profile covers. A step performs four.
    evaluations: u64,
    /// Time in each term, indexed by [`RhsTerm::index`].
    terms: [Duration; RHS_TERMS.len()],
    /// How many times each term was charged, indexed the same way.
    charges: [u32; RHS_TERMS.len()],
}

impl TermProfile {
    /// Evaluations this profile covers.
    #[must_use]
    pub const fn evaluations(&self) -> u64 {
        self.evaluations
    }

    /// Time spent in `term` across the whole profile.
    #[must_use]
    pub const fn term(&self, term: RhsTerm) -> Duration {
        self.terms[term.index()]
    }

    /// Mean time one evaluation spent in `term`.
    #[must_use]
    pub fn per_evaluation(&self, term: RhsTerm) -> Duration {
        mean_over(self.term(term), self.evaluations)
    }

    /// Mean time one evaluation spent in every term together.
    #[must_use]
    pub fn per_evaluation_total(&self) -> Duration {
        mean_over(self.total(), self.evaluations)
    }

    /// How many times `term` was charged over the whole profile.
    ///
    /// One per evaluation, for every term, which is what makes a missing row
    /// visible without asserting on a clock — see
    /// [`StepProfile::charges`].
    #[must_use]
    pub const fn charges(&self, term: RhsTerm) -> u32 {
        self.charges[term.index()]
    }

    /// Time spent in every term together.
    ///
    /// Not the same as the two evaluator phases of a [`StepProfile`]: timing
    /// the kernels one at a time gives the optimiser less to work with, so
    /// this total runs longer than the evaluators do when left alone. The
    /// difference is what fusing them across a call boundary was worth, and a
    /// report states it rather than hiding it.
    #[must_use]
    pub fn total(&self) -> Duration {
        self.terms.iter().sum()
    }

    /// The fraction of the measured kernel time spent in `term`, between 0 and
    /// 1.
    #[must_use]
    pub fn share(&self, term: RhsTerm) -> f64 {
        let total = self.total().as_secs_f64();
        if total == 0.0 {
            return 0.0;
        }
        self.term(term).as_secs_f64() / total
    }
}

/// One right-hand-side evaluation of a scenario, term by term.
///
/// It calls the same kernels [`ShallowWaterRhs::evaluate`] and
/// [`CoriolisTerm::add_to_tendency`] call — the `termocline-numerics`
/// operators and the crate-private array loops beside them — in the same order
/// over the same buffers, so its tendency is theirs.
/// `tests/step_profile.rs` asserts exactly that, bit for bit, which is what
/// makes a term table a statement about the engine rather than about this
/// module.
#[derive(Debug)]
pub struct TermProfiler {
    /// Physical parameters `(g', H, r, β, ρ₀)` of the ocean being evaluated.
    params: PhysicalParams,
    /// The beta-plane the rotation terms read `f = β·y` from.
    plane: BetaPlane,
    /// The C-grid derivative and interpolation operators at this spacing.
    operators: CGridOperators,
    /// `∂u/∂x` at cell centres, in s⁻¹. Scratch, rewritten every evaluation.
    zonal_divergence_per_s: Field2D<f64>,
    /// `∂v/∂y` at cell centres, in s⁻¹. Scratch, rewritten every evaluation.
    meridional_divergence_per_s: Field2D<f64>,
    /// `v` on the east/west faces. Scratch, rewritten every evaluation.
    v_on_u_faces: Field2D<f64>,
    /// `u` on the north/south faces. Scratch, rewritten every evaluation.
    u_on_v_faces: Field2D<f64>,
    /// Time accumulated in each term.
    terms: [Duration; RHS_TERMS.len()],
    /// How many times each term has been charged.
    charges: [u32; RHS_TERMS.len()],
    /// Evaluations performed so far.
    evaluations: u64,
}

impl TermProfiler {
    /// A term profiler for `grid` at `spacing`, for an ocean with `params` on
    /// `plane`.
    ///
    /// Private: every caller has a [`Scenario`] in hand, and
    /// [`TermProfiler::of_scenario`] is the constructor that cannot place a
    /// profiler on a different beta-plane from the run it is profiling.
    fn new(grid: Grid, spacing: Spacing, params: PhysicalParams, plane: BetaPlane) -> Self {
        Self {
            params,
            plane,
            operators: CGridOperators::new(grid, spacing),
            zonal_divergence_per_s: grid.allocate(H_STAGGERING, 0.0),
            meridional_divergence_per_s: grid.allocate(H_STAGGERING, 0.0),
            v_on_u_faces: grid.allocate(U_STAGGERING, 0.0),
            u_on_v_faces: grid.allocate(V_STAGGERING, 0.0),
            terms: [Duration::ZERO; RHS_TERMS.len()],
            charges: [0; RHS_TERMS.len()],
            evaluations: 0,
        }
    }

    /// A term profiler for `scenario`'s basin, parameters and beta-plane.
    #[must_use]
    pub fn of_scenario(scenario: &Scenario) -> Self {
        let basin = scenario.basin();
        let params = scenario.physical_params();
        Self::new(
            basin.grid(),
            basin.spacing(),
            params,
            BetaPlane::of_basin(params, basin),
        )
    }

    /// Write the full tendency of `state` under `wind_stress` into `tendency`,
    /// charging each kernel the time it took.
    ///
    /// The full tendency: the shallow-water terms *and* rotation, which is what
    /// an RK4 stage evaluates. The boundary condition is not applied — that is
    /// [`StepPhase::BoundaryCondition`], and it belongs to the step rather than
    /// to the evaluators.
    ///
    /// # Panics
    /// If `state`, `wind_stress` or `tendency` covers a different basin from
    /// the one this profiler was built for — the panic the operators raise,
    /// for the reason [`ShallowWaterRhs::evaluate`] gives.
    pub fn evaluate(
        &mut self,
        state: &OceanState,
        wind_stress: &WindStressField,
        tendency: &mut OceanState,
    ) {
        let minus_g_prime_m_per_s2 = -self.params.reduced_gravity_m_per_s2();
        let layer_mass_kg_per_m2 =
            self.params.reference_density_kg_per_m3() * self.params.mean_thermocline_depth_m();
        let damping_per_s = self.params.rayleigh_damping_per_s();
        let minus_mean_depth_m = -self.params.mean_thermocline_depth_m();

        let mut charge = |term: RhsTerm, started: Instant| {
            self.terms[term.index()] += started.elapsed();
            self.charges[term.index()] += 1;
            Instant::now()
        };

        let mut mark = Instant::now();
        self.operators
            .ddx_center_to_face(state.h(), tendency.u_mut());
        mark = charge(RhsTerm::ZonalPressureGradient, mark);

        self.operators
            .ddy_center_to_face(state.h(), tendency.v_mut());
        mark = charge(RhsTerm::MeridionalPressureGradient, mark);

        turn_gradient_into_acceleration(
            tendency.u_mut(),
            minus_g_prime_m_per_s2,
            wind_stress.tau_x_pa(),
            layer_mass_kg_per_m2,
        );
        mark = charge(RhsTerm::ZonalPressureAndStress, mark);

        turn_gradient_into_acceleration(
            tendency.v_mut(),
            minus_g_prime_m_per_s2,
            wind_stress.tau_y_pa(),
            layer_mass_kg_per_m2,
        );
        mark = charge(RhsTerm::MeridionalPressureAndStress, mark);

        subtract_damping(tendency.u_mut(), state.u(), damping_per_s);
        mark = charge(RhsTerm::ZonalDamping, mark);

        subtract_damping(tendency.v_mut(), state.v(), damping_per_s);
        mark = charge(RhsTerm::MeridionalDamping, mark);

        self.operators
            .ddx_face_to_center(state.u(), &mut self.zonal_divergence_per_s);
        mark = charge(RhsTerm::ZonalDivergence, mark);

        self.operators
            .ddy_face_to_center(state.v(), &mut self.meridional_divergence_per_s);
        mark = charge(RhsTerm::MeridionalDivergence, mark);

        write_continuity(
            tendency.h_mut(),
            &self.zonal_divergence_per_s,
            &self.meridional_divergence_per_s,
            minus_mean_depth_m,
        );
        mark = charge(RhsTerm::Continuity, mark);

        subtract_damping(tendency.h_mut(), state.h(), damping_per_s);
        mark = charge(RhsTerm::ThermoclineDamping, mark);

        self.operators
            .face_y_to_face_x(state.v(), &mut self.v_on_u_faces);
        mark = charge(RhsTerm::MeridionalVelocityOntoZonalFaces, mark);

        self.operators
            .face_x_to_face_y(state.u(), &mut self.u_on_v_faces);
        mark = charge(RhsTerm::ZonalVelocityOntoMeridionalFaces, mark);

        let plane = self.plane;
        accumulate_rows(tendency.u_mut(), &self.v_on_u_faces, |j| {
            plane.coriolis_at_row_per_s(U_STAGGERING, j)
        });
        mark = charge(RhsTerm::ZonalRotation, mark);

        accumulate_rows(tendency.v_mut(), &self.u_on_v_faces, |j| {
            -plane.coriolis_at_row_per_s(V_STAGGERING, j)
        });
        let _ = charge(RhsTerm::MeridionalRotation, mark);

        self.evaluations += 1;
    }

    /// What has been measured so far.
    #[must_use]
    pub const fn profile(&self) -> TermProfile {
        TermProfile {
            evaluations: self.evaluations,
            terms: self.terms,
            charges: self.charges,
        }
    }
}

/// `total` divided by `count`, or zero if there is nothing to divide.
///
/// One definition, because every mean in this module is a duration over a
/// count of steps, evaluations or charges, and `Duration` divides only by a
/// `u32`.
fn mean_over(total: Duration, count: u64) -> Duration {
    total
        .checked_div(u32::try_from(count).unwrap_or(u32::MAX))
        .unwrap_or_default()
}

/// The cost of one [`Instant::now`], measured by taking `reads` of them back
/// to back.
///
/// A decomposition charges the phases it names for the clock reads that
/// delimit them, so a reader needs this number to know whether the
/// decomposition is worth reading at all: [`CLOCK_READS_PER_STEP`] of these
/// against the duration of one step is the instrument's own footprint.
///
/// # Panics
/// If `reads` is zero, which is a caller asking for a mean of nothing.
#[must_use]
pub fn clock_read_cost(reads: u32) -> Duration {
    assert!(
        reads > 0,
        "the cost of a clock read needs at least one read"
    );
    // One discarded pass first. The wanted figure is what a clock read costs
    // during a profile — inside a loop the processor has already spun up for —
    // and the first reads of a cold process measure the ramp instead, which on
    // a laptop is a factor of two or three. Whatever ramp is left after the
    // discarded pass makes this an over-estimate of the instrument's footprint,
    // which is the safe direction for a number a reader uses to decide whether
    // to trust a table.
    for _ in 0..reads {
        std::hint::black_box(Instant::now());
    }
    let started = Instant::now();
    for _ in 0..reads {
        // `black_box` so the loop cannot be folded away: what is wanted is the
        // cost of asking the operating system, not of a constant.
        std::hint::black_box(Instant::now());
    }
    started.elapsed() / reads
}

#[cfg(test)]
mod tests {
    use super::{CLOCK_READS_PER_STEP, RHS_TERMS, STEP_PHASES};

    #[test]
    fn a_profiled_step_reads_the_clock_twenty_four_times() {
        // The constant is arithmetic over three others; this is the number it
        // comes to, written out once so that a change to any of them has to be
        // acknowledged here rather than silently changing what the example and
        // `docs/performance-notes.md` say the instrument costs.
        assert_eq!(CLOCK_READS_PER_STEP, 24);
    }

    #[test]
    fn every_phase_sits_at_its_own_index_in_the_phase_list() {
        // A phase's index is where a profile keeps its measurement, and
        // `STEP_PHASES` is the order a report walks. Two statements of one
        // ordering can disagree, and a disagreement would put a measurement in
        // another row's slot rather than fail.
        for (position, phase) in STEP_PHASES.into_iter().enumerate() {
            assert_eq!(phase.index(), position, "{phase:?} is listed out of order");
        }
    }

    #[test]
    fn every_term_sits_at_its_own_index_in_the_term_list() {
        for (position, term) in RHS_TERMS.into_iter().enumerate() {
            assert_eq!(term.index(), position, "{term:?} is listed out of order");
        }
    }
}
