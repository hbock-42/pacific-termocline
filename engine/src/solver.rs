//! The time loop's single step: the whole right-hand side, advanced by RK4.
//!
//! This is where Epic 02's terms meet Epic 01's integrator. [`Solver`] holds
//! the three pieces a step needs — the pressure-gradient, continuity, surface
//! stress and Rayleigh damping terms of [`ShallowWaterRhs`], the beta-plane
//! rotation of [`CoriolisTerm`], and an [`Rk4`] sized for the basin — and
//! composes them into the `state_{n+1} = rk4_step(state_n, dt, forcing_fn)`
//! that [ADR-0003] describes:
//!
//! ```text
//! ∂u/∂t = +f·v − g'·∂h/∂x + τx/(ρ₀·H) − r·u
//! ∂v/∂t = −f·u − g'·∂h/∂y + τy/(ρ₀·H) − r·v
//! ∂h/∂t =       −H·(∂u/∂x + ∂v/∂y)    − r·h
//! ```
//!
//! Two ways in exist, and they differ only in how the forcing arrives.
//! [`Solver::step`] takes a [`WindStressField`] already sampled onto the
//! C-grid, which is what the Epic 02 term tests want; [`Solver::step_forced_by`]
//! takes a [`WindStress`] and a [`Basin`] and re-samples it at each stage,
//! which is what a scenario wants.
//!
//! The two halves stay separate types rather than merging into one evaluator:
//! rotation is the only term that needs the basin's position on the
//! beta-plane, and keeping it out of [`ShallowWaterRhs`] keeps the terms that
//! do not care about latitude testable without one. Composing them is this
//! module's whole job, and the order is fixed by their contracts — the
//! shallow-water evaluator *writes* every point of the tendency, the Coriolis
//! term *adds* to it, and [`NoNormalFlow`] overrules both at the walls.
//!
//! # The closed basin
//!
//! [`NoNormalFlow`] is applied here rather than inside either evaluator,
//! because it is a statement about the *system being integrated* rather than
//! about any one term: the normal velocity at a coast is not a degree of
//! freedom, so no term may accelerate it. A step therefore puts the incoming
//! state on the condition once, and holds the walls at rest at each of RK4's
//! four stages; the `boundary` module explains why both, and why per stage
//! rather than per step.
//!
//! # What the solver owns, and why
//!
//! A `Solver` is built once per run and stepped many times. Everything a step
//! writes into — the four RK4 stage buffers, the stage state, the divergence
//! scratch, the two Coriolis interpolation buffers — is allocated here and
//! reused, so a whole simulation allocates them exactly once
//! (CODING_STANDARDS.md § Performance).
//!
//! The timestep is part of that construction rather than an argument to each
//! step, because it is what makes a solver valid: [`Solver::new`] refuses a
//! `dt` the scheme cannot take instead of shortening it
//! (CODING_STANDARDS.md § No silent clamping), and a solver that exists is one
//! whose steps are stable.
//!
//! # Two bounds on the timestep, not one
//!
//! The CFL bound of T-01.3 is derived for the gravity-wave terms alone: it
//! keeps `c·κ_max·dt` inside RK4's stability region. Rotation is a second
//! oscillation in the same system — the momentum pair `u̇ = +f·v`,
//! `v̇ = −f·u` has eigenvalues `±i·f` — and it carries its own limit,
//! `|f|·dt ≤ 2√2`, which the gravity-wave bound knows nothing about.
//!
//! The two are independent, and either can be the binding one. On a basin
//! resolved finely enough for its Kelvin waves the gravity-wave bound wins by
//! a wide margin, which is the case `termocline-numerics` assumes when it
//! calls rotation something the safety factor absorbs. On a coarse basin
//! reaching far from the equator it is the other way round: at 625 km cells
//! over ±2500 km the CFL bound admits a 32-hour step while the inertial
//! period at the northern wall is 30 hours, and a run at that step amplifies
//! by a factor of 70 per step at the wall. This is the first module that
//! sees both terms at once, so it is the first that can check both, and it
//! refuses a timestep either bound rejects. [ADR-0007] records that decision
//! and the alternatives weighed against it.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md
//! [ADR-0007]: ../../docs/planning/adr/0007-rotation-timestep-bound.md

use termocline_grid::Grid;

use crate::basin::Basin;

use crate::boundary::NoNormalFlow;
use crate::coriolis::{BetaPlane, CoriolisTerm};
use crate::forcing::{WindStress, WindStressField};
use crate::integrator::Rk4;
use crate::params::PhysicalParams;
use crate::shallow_water::ShallowWaterRhs;
use crate::state::OceanState;
use termocline_numerics::{
    check_timestep, CflError, Spacing, WaveSpeed, CFL_SAFETY_FACTOR, RK4_IMAGINARY_AXIS_LIMIT,
};

use std::fmt;

/// Why a solver could not be built.
///
/// Both variants describe invalid *scenario input* — a run asking for a
/// timestep the scheme cannot take — so they are returned rather than
/// panicked, and each names the value it rejected and the bound it violated.
#[derive(Debug, Clone, PartialEq)]
pub enum SolverError {
    /// The timestep was rejected by the gravity-wave CFL bound of T-01.3.
    Cfl(CflError),
    /// The timestep was longer than the rotation of the beta-plane allows.
    ///
    /// The run is refused rather than quietly shortened, per
    /// CODING_STANDARDS.md § *No silent clamping*.
    TimestepExceedsRotationLimit {
        /// The timestep asked for, in seconds.
        requested_s: f64,
        /// The largest timestep this basin's rotation allows, in seconds.
        max_stable_s: f64,
        /// The largest `|f| = |β·y|` in the basin, in s⁻¹ — the rotation rate
        /// the bound is derived from.
        largest_coriolis_per_s: f64,
    },
}

impl fmt::Display for SolverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cfl(error) => error.fmt(f),
            Self::TimestepExceedsRotationLimit {
                requested_s,
                max_stable_s,
                largest_coriolis_per_s,
            } => write!(
                f,
                "dt is {requested_s} s, but the Coriolis parameter reaches \
                 {largest_coriolis_per_s} s⁻¹ at this basin's meridional boundary, whose \
                 inertial oscillation RK4 can only follow up to {max_stable_s} s; the run would \
                 go unstable there. Set dt to at most {max_stable_s} s, or bring the basin's \
                 boundaries closer to the equator"
            ),
        }
    }
}

impl std::error::Error for SolverError {}

impl From<CflError> for SolverError {
    fn from(error: CflError) -> Self {
        Self::Cfl(error)
    }
}

/// The longest timestep, in seconds, at which RK4 stays stable on an
/// oscillation of angular rate `largest_coriolis_per_s`.
///
/// The rotation pair `u̇ = +f·v`, `v̇ = −f·u` has eigenvalues `±i·f`, so this
/// is `RK4_IMAGINARY_AXIS_LIMIT / |f|` with the same
/// [`CFL_SAFETY_FACTOR`] margin the gravity-wave bound holds back — the two
/// bounds are the same stability region read against two different
/// oscillations, and there is no reason to trust one of them closer to its
/// boundary than the other.
///
/// A basin whose boundaries lie exactly on the equator has `f ≡ 0` and no
/// rotation limit at all, which is what the infinite bound means.
#[must_use]
fn max_stable_dt_for_rotation(largest_coriolis_per_s: f64) -> f64 {
    CFL_SAFETY_FACTOR * RK4_IMAGINARY_AXIS_LIMIT / largest_coriolis_per_s
}

/// One basin's time loop: the full right-hand side, stepped by RK4 at a fixed,
/// CFL-checked timestep.
#[derive(Debug, Clone)]
pub struct Solver {
    /// Length of one step, in seconds. Checked against the CFL bound at
    /// construction and constant for the life of the solver.
    dt_s: f64,
    /// The pressure-gradient, continuity, surface-stress and damping terms.
    rhs: ShallowWaterRhs,
    /// The beta-plane rotation terms, added on top of them.
    coriolis: CoriolisTerm,
    /// The integrator and its stage buffers, sized for this basin.
    integrator: Rk4<OceanState>,
    /// The surface stress of the current RK4 stage, re-sampled in place by
    /// [`Solver::step_forced_by`]. Allocated here so that a run driven by a
    /// [`WindStress`] allocates its forcing once rather than once per stage
    /// (CODING_STANDARDS.md § Performance).
    stage_stress: WindStressField,
}

impl Solver {
    /// A solver for `grid` at `spacing`, for an ocean with `params` on `plane`,
    /// stepping `dt_s` seconds at a time.
    ///
    /// # Errors
    /// [`SolverError::Cfl`] if `dt_s` is not a finite, positive duration or is
    /// longer than the gravity-wave CFL maximum for this spacing and the
    /// Kelvin wave speed `√(g'·H)` implied by `params`, and
    /// [`SolverError::TimestepExceedsRotationLimit`] if it is longer than the
    /// basin's rotation allows. The timestep is never adjusted: an unstable
    /// one is the scenario's error to fix.
    pub fn new(
        grid: Grid,
        spacing: Spacing,
        params: PhysicalParams,
        plane: BetaPlane,
        dt_s: f64,
    ) -> Result<Self, SolverError> {
        let wave_speed = WaveSpeed::new(params.kelvin_wave_speed_m_per_s())
            .expect("physical parameters are validated positive, so `√(g'·H)` is too");
        check_timestep(dt_s, spacing, wave_speed)?;

        let largest_coriolis_per_s = plane.largest_coriolis_magnitude_per_s(grid);
        let max_stable_s = max_stable_dt_for_rotation(largest_coriolis_per_s);
        if dt_s > max_stable_s {
            return Err(SolverError::TimestepExceedsRotationLimit {
                requested_s: dt_s,
                max_stable_s,
                largest_coriolis_per_s,
            });
        }

        Ok(Self {
            dt_s,
            rhs: ShallowWaterRhs::new(grid, spacing, params),
            coriolis: CoriolisTerm::new(grid, spacing, plane),
            integrator: Rk4::new(&OceanState::at_rest(grid)),
            stage_stress: WindStressField::calm(grid),
        })
    }

    /// Length of one step, in seconds.
    #[must_use]
    pub const fn dt_s(&self) -> f64 {
        self.dt_s
    }

    /// Advance `state` from time `t_s` to `t_s + dt`, in place.
    ///
    /// `wind_stress_at(t)` supplies the surface stress at a given time in
    /// seconds; it is called once per RK4 stage, at the tableau's nodes `t`,
    /// `t + dt/2`, `t + dt/2` and `t + dt`. A steady scenario samples its
    /// [`WindStress`](crate::WindStress) once and hands back the same field at
    /// every stage — `|_t| &stress` — but the argument is a function of time
    /// because the seasonal and burst scenarios of Epic 03 are, and a step
    /// that sampled a varying stress once would integrate the wrong forcing.
    ///
    /// `state` is brought onto the closed basin's boundary condition on the
    /// way in: a normal velocity at a wall is not a degree of freedom, so a
    /// state handed in carrying one has it set to zero rather than integrated
    /// (see [`NoNormalFlow`]). That is stated here because it is a write the
    /// caller did not ask for — everything the step does *after* it preserves
    /// the condition exactly, so only a state that arrives off it is changed.
    ///
    /// # Panics
    /// If `state` or a returned wind stress covers a different basin from the
    /// one this solver was built for. A shape mismatch means the calling code
    /// is wrong, which is what panics are for (CODING_STANDARDS.md
    /// § Correctness and failure).
    pub fn step<'w, W>(&mut self, state: &mut OceanState, t_s: f64, wind_stress_at: W)
    where
        W: Fn(f64) -> &'w WindStressField,
    {
        // Destructured so the closure below can borrow the two evaluators while
        // the integrator is borrowed mutably: they are disjoint fields.
        let Self {
            dt_s,
            rhs,
            coriolis,
            integrator,
            stage_stress: _,
        } = self;
        NoNormalFlow::apply_to_state(state);
        integrator.step(
            state,
            t_s,
            *dt_s,
            &mut |now: &OceanState, stage_t_s: f64, tendency: &mut OceanState| {
                // Order matters: the shallow-water evaluator writes every point
                // of the tendency, the Coriolis term adds to what it wrote, and
                // the boundary condition has the last word over both.
                rhs.evaluate(now, wind_stress_at(stage_t_s), tendency);
                coriolis.add_to_tendency(now, tendency);
                NoNormalFlow::apply_to_tendency(tendency);
            },
        );
    }

    /// Advance `state` from time `t_s` to `t_s + dt` under `wind`, in place.
    ///
    /// The form a scenario uses: it takes the [`WindStress`] itself rather
    /// than a field already sampled from one, and re-samples it over `basin`
    /// at each of RK4's four stage times. That is what makes the forcing a
    /// *function of time* all the way through the integration — a steady
    /// scenario gets the same field four times, a seasonal or burst scenario
    /// gets the four the tableau's nodes actually ask for.
    ///
    /// Sampling writes into a buffer the solver owns, so a whole run allocates
    /// its forcing exactly once (CODING_STANDARDS.md § Performance).
    ///
    /// `state` is brought onto the boundary condition on the way in, exactly
    /// as in [`Solver::step`].
    ///
    /// # Panics
    /// If `state` or `basin` covers a different grid from the one this solver
    /// was built for. A shape mismatch means the calling code is wrong, which
    /// is what panics are for (CODING_STANDARDS.md § Correctness and failure).
    pub fn step_forced_by<W: WindStress + ?Sized>(
        &mut self,
        state: &mut OceanState,
        t_s: f64,
        basin: Basin,
        wind: &W,
    ) {
        let Self {
            dt_s,
            rhs,
            coriolis,
            integrator,
            stage_stress,
        } = self;
        NoNormalFlow::apply_to_state(state);
        integrator.step(
            state,
            t_s,
            *dt_s,
            &mut |now: &OceanState, stage_t_s: f64, tendency: &mut OceanState| {
                stage_stress.sample(basin, wind, stage_t_s);
                rhs.evaluate(now, stage_stress, tendency);
                coriolis.add_to_tendency(now, tendency);
                NoNormalFlow::apply_to_tendency(tendency);
            },
        );
    }
}

/// `state` advanced by one step of `dt_s` seconds from time zero, as a freshly
/// allocated state.
///
/// The named deliverable of T-02.5, and the convenient form for a test or a
/// one-off step. It allocates a solver and a state per call, so a time loop
/// builds one [`Solver`] and steps it instead (CODING_STANDARDS.md
/// § Performance).
///
/// `wind_stress_at` is [`Solver::step`]'s, and is sampled at the same four
/// stage times — measured from zero, since a one-off step has no run behind it
/// to place it in.
///
/// # Errors
/// The errors of [`Solver::new`]: a timestep that is not a finite, positive
/// duration, or one past either the gravity-wave or the rotation bound.
pub fn step<'w, W>(
    state: &OceanState,
    dt_s: f64,
    params: PhysicalParams,
    spacing: Spacing,
    plane: BetaPlane,
    wind_stress_at: W,
) -> Result<OceanState, SolverError>
where
    W: Fn(f64) -> &'w WindStressField,
{
    let mut solver = Solver::new(state.grid(), spacing, params, plane, dt_s)?;
    let mut advanced = state.clone();
    solver.step(&mut advanced, 0.0, wind_stress_at);
    Ok(advanced)
}
