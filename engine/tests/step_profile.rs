//! Acceptance tests for T-10.2 — the instrument behind
//! `docs/performance-notes.md`.
//!
//! The ticket's acceptance criterion is that the findings are backed by a
//! profile rather than by guesswork. A test cannot assert that a note is
//! honest, and it certainly cannot assert a duration — a wall-clock threshold
//! in a test suite is the flaky measurement `docs/benchmarks.md` exists to
//! avoid. What *is* testable is the only thing that makes a profile worth
//! believing: that the code being timed is the code the engine runs.
//!
//! So these tests pin the two profilers of `engine::profiling` to the engine
//! itself.
//!
//! - The instrumented step must reach the same state as
//!   [`Solver::step_forced_by`], bit for bit, over enough steps that a
//!   difference in any term would have shown. A phase table taken from a step
//!   that was not the solver's step would describe a program nobody runs.
//! - The term-by-term evaluation must produce the same tendency as
//!   [`ShallowWaterRhs::evaluate`] followed by
//!   [`CoriolisTerm::add_to_tendency`], bit for bit. Its fourteen kernels are
//!   the ones the evaluators call, and this is what says so.
//! - The two decompositions must be exhaustive: every phase and every term is
//!   charged the number of times a step charges it, and the phase shares sum
//!   to one, so no work can sit outside the table.
//!
//! Bit-for-bit rather than within a tolerance, deliberately. The profilers
//! perform the same operations in the same order over the same `f64` fields,
//! so the results are identical or the profiler is doing something else; a
//! tolerance here would be a place for a real divergence to hide.
//!
//! **Nothing here asserts a duration**, not even a positive one. `cargo test`
//! is the CI gate (AGENTS.md § *The CI gate*) and CI runs on a shared runner,
//! where a wall-clock threshold fails for reasons that have nothing to do with
//! the change — the reasoning `docs/benchmarks.md` § *What is not measured*
//! gives for keeping the benchmark job off the gate. A phase that the profiler
//! forgot to charge is caught by its charge *count*, which is four or five per
//! step on every machine, where a zero duration would only appear on a fast
//! one.

use engine::benchmark::BenchmarkWorkload;
use engine::benchmark::BENCHMARK_WORKLOADS;
use engine::profiling::{
    RhsTerm, StepPhase, StepProfiler, TermProfiler, CLOCK_READS_PER_STEP, RHS_TERMS, STEP_PHASES,
};
use engine::{BetaPlane, CoriolisTerm, OceanState, Scenario, ShallowWaterRhs, Solver};

/// The workload the profile in `docs/performance-notes.md` is taken on: the
/// control scenario at its own 0.5° resolution, the finer of the two the
/// benchmark suite measures.
fn control_workload() -> BenchmarkWorkload {
    *BENCHMARK_WORKLOADS
        .last()
        .expect("the benchmark suite defines at least one workload")
}

/// That workload's scenario.
fn control_workload_scenario() -> Scenario {
    control_workload().scenario()
}

/// Stages of the RK4 tableau, and so how many times a step evaluates the right
/// hand side and charges each evaluator phase.
const RK4_STAGES: u32 = 4;

/// Steps the comparison against the solver takes.
///
/// Enough that a difference in any one term would have propagated across the
/// basin and shown up: the shallow-water terms couple `h`, `u` and `v` at
/// every step, so an error in one of the fourteen kernels reaches all three
/// within a handful of steps and grows. Small enough that the test costs a
/// fraction of a second at 320 × 100.
const COMPARED_STEPS: u64 = 8;

#[test]
fn the_profiled_step_is_the_solvers_own_step() {
    // The whole basis of the phase table: if this differs, the profile
    // describes a program nobody runs.
    let scenario = control_workload_scenario();
    let basin = scenario.basin();
    let params = scenario.physical_params();
    let schedule = scenario.output_schedule();
    let wind = scenario.wind();

    let mut solver = Solver::new(
        basin.grid(),
        basin.spacing(),
        params,
        BetaPlane::of_basin(params, basin),
        schedule.dt_s(),
    )
    .expect("the control scenario's timestep clears both bounds");
    let mut solver_state = OceanState::at_rest(basin.grid());

    let mut profiler =
        StepProfiler::new(&scenario).expect("the profiler accepts what the solver accepts");

    for step in 0..COMPARED_STEPS {
        let t_s = schedule.model_time_at_step(step);
        solver.step_forced_by(&mut solver_state, t_s, basin, &wind);
        profiler.step();
        assert_eq!(
            profiler.state(),
            &solver_state,
            "the instrumented step diverged from the solver's step at step {}",
            step + 1
        );
    }
}

#[test]
fn every_phase_of_a_step_is_reported_and_the_shares_sum_to_one() {
    // A decomposition with a hole in it is worse than none: the missing work
    // would be silently credited to whichever phase a reader assumed.
    // `StageAlgebra` is the residual of the step, which is what makes this
    // exact rather than approximate.
    let scenario = control_workload_scenario();
    let mut profiler =
        StepProfiler::new(&scenario).expect("the profiler accepts what the solver accepts");
    let profile = profiler.profile(2);

    assert_eq!(profile.steps(), 2);

    let shares: f64 = STEP_PHASES.iter().map(|phase| profile.share(*phase)).sum();
    // Exact to a few ulps, not to a tolerance on the physics: the five
    // durations are constructed to sum to the total, so the only error here is
    // the rounding of five f64 divisions. 1e-12 is far above that and far
    // below any real gap.
    assert!(
        (shares - 1.0).abs() < 1e-12,
        "the five phase shares summed to {shares}, so a step has unattributed work"
    );

    let named: std::time::Duration = STEP_PHASES.iter().map(|phase| profile.phase(*phase)).sum();
    assert_eq!(
        named,
        profile.total(),
        "the phases must partition the step exactly"
    );
}

#[test]
fn each_phase_is_charged_once_for_every_time_a_step_performs_it() {
    // The count, not the duration: a phase the profiler forgot to charge is a
    // missing column of the table, and this is what catches it without
    // depending on how fast the machine under CI happens to be.
    //
    // The three evaluator phases are charged once per RK4 stage, because that
    // is how often a step evaluates the right-hand side. The boundary
    // condition is charged once more than that, for the incoming state a step
    // puts on the condition before integrating it (`Solver::step_forced_by`).
    // The stage algebra is the residual of the step, so it is charged once per
    // step.
    const STEPS: u64 = 3;
    let scenario = control_workload_scenario();
    let mut profiler =
        StepProfiler::new(&scenario).expect("the profiler accepts what the solver accepts");
    let profile = profiler.profile(STEPS);

    let steps = u32::try_from(STEPS).expect("three steps fit in a u32");
    for (phase, expected) in [
        (StepPhase::WindStressSampling, steps * RK4_STAGES),
        (StepPhase::ShallowWaterTerms, steps * RK4_STAGES),
        (StepPhase::Coriolis, steps * RK4_STAGES),
        (StepPhase::BoundaryCondition, steps * (RK4_STAGES + 1)),
        (StepPhase::StageAlgebra, steps),
    ] {
        assert_eq!(
            profile.charges(phase),
            expected,
            "phase {phase:?} was charged the wrong number of times, so the decomposition does \
             not describe what a step does"
        );
    }
}

#[test]
fn the_term_by_term_evaluation_is_the_evaluators_own() {
    // The basis of the term table, and the reason the optimisation tickets may
    // act on it: the fourteen kernels are the ones the two evaluators call, in
    // the order they call them.
    let workload = control_workload();
    let scenario = workload.scenario();
    let basin = scenario.basin();
    let params = scenario.physical_params();
    let plane = BetaPlane::of_basin(params, basin);

    // The benchmark suite's own state and stress: an equatorial Kelvin
    // structure under the control scenario's wind, so the comparison is made
    // on a state a run could actually be in rather than on an ocean at rest,
    // where several terms would be identically zero and could not disagree.
    let state = workload.benchmark_state();
    let wind_stress = workload.wind_stress();

    let mut rhs = ShallowWaterRhs::new(basin.grid(), basin.spacing(), params);
    let mut coriolis = CoriolisTerm::new(basin.grid(), basin.spacing(), plane);
    let mut expected = OceanState::at_rest(basin.grid());
    rhs.evaluate(&state, &wind_stress, &mut expected);
    coriolis.add_to_tendency(&state, &mut expected);

    let mut profiler = TermProfiler::of_scenario(&scenario);
    let mut measured = OceanState::at_rest(basin.grid());
    profiler.evaluate(&state, &wind_stress, &mut measured);

    assert_eq!(
        measured, expected,
        "the term-by-term evaluation is not the evaluators' own arithmetic"
    );
}

#[test]
fn a_term_profile_reuses_its_buffers_across_evaluations() {
    // The tendency buffer is overwritten in full by each evaluation, exactly
    // as an RK4 stage requires (`ShallowWaterRhs::evaluate`). A profiler that
    // let one evaluation leak into the next would report the second one's
    // kernels running over different data from the first's.
    let workload = control_workload();
    let scenario = workload.scenario();
    let state = workload.benchmark_state();
    let wind_stress = workload.wind_stress();

    let mut profiler = TermProfiler::of_scenario(&scenario);
    let mut first = OceanState::at_rest(scenario.basin().grid());
    profiler.evaluate(&state, &wind_stress, &mut first);
    let mut second = first.clone();
    profiler.evaluate(&state, &wind_stress, &mut second);

    assert_eq!(
        first, second,
        "two evaluations of one state produced different tendencies"
    );
    assert_eq!(profiler.profile().evaluations(), 2);
}

#[test]
fn every_term_is_reported_and_belongs_to_one_evaluator_phase() {
    // The term table has to be exhaustive in the same sense the phase table
    // is: fourteen distinct kernels, each charged once per evaluation, each
    // attributed to the evaluator it came from so the two tables can be read
    // against each other.
    const EVALUATIONS: u32 = 3;
    let workload = control_workload();
    let scenario = workload.scenario();
    let state = workload.benchmark_state();
    let wind_stress = workload.wind_stress();

    let mut profiler = TermProfiler::of_scenario(&scenario);
    let mut tendency = OceanState::at_rest(scenario.basin().grid());
    for _ in 0..EVALUATIONS {
        profiler.evaluate(&state, &wind_stress, &mut tendency);
    }
    let profile = profiler.profile();

    for term in RHS_TERMS {
        assert_eq!(
            profile.charges(term),
            EVALUATIONS,
            "term {term:?} was charged the wrong number of times, so the decomposition has a \
             missing or duplicated row"
        );
        assert!(
            matches!(
                term.phase(),
                StepPhase::ShallowWaterTerms | StepPhase::Coriolis
            ),
            "term {:?} claims to belong to a phase that evaluates no terms",
            term
        );
    }

    let shares: f64 = RHS_TERMS.iter().map(|term| profile.share(*term)).sum();
    assert!(
        (shares - 1.0).abs() < 1e-12,
        "the fourteen term shares summed to {shares}"
    );

    // Every label distinct: two rows of a table reading the same is a table a
    // reader cannot act on.
    let mut labels: Vec<&str> = RHS_TERMS.iter().map(|term| term.label()).collect();
    labels.sort_unstable();
    let count = labels.len();
    labels.dedup();
    assert_eq!(labels.len(), count, "two terms report the same label");
}

#[test]
fn the_four_coriolis_terms_are_the_ones_the_rotation_phase_owns() {
    // The split between the two evaluators is what lets the term table be read
    // against the phase table: ten shallow-water kernels, four rotation ones.
    let rotation = RHS_TERMS
        .iter()
        .filter(|term| term.phase() == StepPhase::Coriolis)
        .count();
    assert_eq!(rotation, 4);
    assert_eq!(RHS_TERMS.len() - rotation, 10);

    // Named individually, so that a term moved between evaluators has to be
    // moved here too rather than silently rebalancing the table.
    for term in [
        RhsTerm::MeridionalVelocityOntoZonalFaces,
        RhsTerm::ZonalVelocityOntoMeridionalFaces,
        RhsTerm::ZonalRotation,
        RhsTerm::MeridionalRotation,
    ] {
        assert_eq!(term.phase(), StepPhase::Coriolis);
    }
}

#[test]
fn the_instruments_footprint_is_the_clock_reads_a_step_actually_performs() {
    // `CLOCK_READS_PER_STEP` is what the example and
    // `docs/performance-notes.md` multiply by the cost of a clock read to
    // state the instrument's footprint. It is only honest if it is the number
    // of reads a step performs, so it is checked against the phase structure
    // it is derived from rather than against a copy of the number.
    //
    // Five marks per RK4 stage — one before the wind is sampled and one after
    // each of the four phases — plus the two pairs a step takes outside its
    // stages, around the step itself and around the incoming boundary
    // condition.
    let per_stage = 1 + u32::try_from(STEP_PHASES.len() - 1).expect("five phases fit in a u32");
    assert_eq!(CLOCK_READS_PER_STEP, 4 + RK4_STAGES * per_stage);
}
