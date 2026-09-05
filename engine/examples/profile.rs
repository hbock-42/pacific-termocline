//! T-10.2 — the driver behind `docs/performance-notes.md`.
//!
//! Two modes, because the note rests on two instruments:
//!
//! ```sh
//! cargo run --release --example profile                # the tables
//! cargo run --release --example profile -- spin 20     # 20 s of hot path
//! ```
//!
//! **The tables.** `phases` (the default) prints, for each grid in
//! `BENCHMARK_WORKLOADS`, how a timestep divides among the five
//! [`StepPhase`]s and how the two evaluator phases divide among the fourteen
//! [`RhsTerm`] kernels. Both come from `engine::profiling`, whose module
//! documentation says what they can and cannot see, and
//! `engine/tests/step_profile.rs` is what pins them to the engine's own step.
//!
//! **The spin.** `spin <seconds>` takes uninstrumented steps through the real
//! [`Solver`] until the time is up, and does nothing else. It exists so that
//! an external sampler — `/usr/bin/sample` on macOS, `perf record` on Linux —
//! can be pointed at a process whose stacks are the hot path and nothing else.
//! No clock is read between phases in that mode, so a sampled profile taken
//! against it is free of the instrument's own footprint, which is exactly the
//! independent check the tables need.
//!
//! It is an example rather than a subcommand of the `termocline` binary
//! because it is a tool for developing the engine, not for running a
//! simulation: nobody profiling a laptop wants it in `termocline --help`.

use std::process::ExitCode;
use std::time::{Duration, Instant};

use engine::benchmark::BENCHMARK_WORKLOADS;
use engine::profiling::{
    clock_read_cost, StepPhase, StepProfiler, TermProfiler, CLOCK_READS_PER_STEP, RHS_TERMS,
    STEP_PHASES,
};
use engine::{BetaPlane, OceanState, Solver};

/// Steps each workload is profiled over.
///
/// A tenth of the benchmark suite's short run, which at 0.5° is a few seconds
/// of profiling — long enough that the per-phase totals are many thousands of
/// clock ticks each, short enough to sit inside an edit-run loop.
const PROFILED_STEPS: u64 = 24;

/// Right-hand-side evaluations one step performs — RK4's four stages — so that
/// the whole-evaluator figure printed beside the term table is per evaluation
/// like the terms are, rather than per step.
const RK4_STAGES_PER_STEP: u32 = 4;

/// How long each workload is stepped before its clock starts.
///
/// A *duration* rather than a step count, because two different things need
/// settling and only one of them scales with the grid. The first step of a run
/// touches every buffer for the first time and pays the page faults for the
/// whole basin; that is over in a handful of steps. The processor's clock is
/// the other, and it takes about a tenth of a second of sustained work to come
/// up — long enough that a fixed step count warms the 0.5° grid and leaves the
/// 1.0° grid, four times cheaper per step, still on a cold core. Measured that
/// way the coarse grid appears to cost 2.2 times the fine grid's *time per
/// cell*, which is an artefact of the processor rather than a property of the
/// engine. A quarter of a second removes it: the two grids then agree on time
/// per cell to within 2%.
const WARMUP: Duration = Duration::from_millis(250);

/// Right-hand-side evaluations the term table is measured over.
///
/// Four per step, so this is the same amount of evaluator work as
/// [`PROFILED_STEPS`] steps.
const PROFILED_EVALUATIONS: u64 = 4 * PROFILED_STEPS;

fn main() -> ExitCode {
    let mut arguments = std::env::args().skip(1);
    match arguments.next().as_deref() {
        None | Some("phases") => {
            print_tables();
            ExitCode::SUCCESS
        }
        Some("spin") => match arguments.next().map(|seconds| seconds.parse::<f64>()) {
            Some(Ok(seconds)) if seconds.is_finite() && seconds > 0.0 => {
                spin(Duration::from_secs_f64(seconds));
                ExitCode::SUCCESS
            }
            _ => {
                eprintln!("usage: profile spin <seconds>, with a positive number of seconds");
                ExitCode::FAILURE
            }
        },
        Some(other) => {
            eprintln!(
                "unknown mode {other:?}; expected `phases` (the default) or `spin <seconds>`"
            );
            ExitCode::FAILURE
        }
    }
}

/// Print the phase and term tables for every benchmark workload.
fn print_tables() {
    // The instrument's own footprint, printed first so the tables below are
    // read against it rather than trusted.
    let clock = clock_read_cost(4096);
    println!("clock: one Instant::now costs {clock:?};");
    println!(
        "       a profiled step reads it {CLOCK_READS_PER_STEP} times, for {:?} of overhead.",
        clock * CLOCK_READS_PER_STEP
    );

    for workload in BENCHMARK_WORKLOADS {
        let scenario = workload.scenario();
        let mut profiler =
            StepProfiler::new(&scenario).expect("a benchmark workload is a runnable scenario");
        // The warm-up steps are taken and then discarded by rebuilding the
        // profiler: a profile covers every step its profiler has ever taken,
        // so the only way to exclude them is not to have them in it.
        let warming = Instant::now();
        while warming.elapsed() < WARMUP {
            profiler.step();
        }
        let mut profiler =
            StepProfiler::new(&scenario).expect("a benchmark workload is a runnable scenario");
        let profile = profiler.profile(PROFILED_STEPS);

        println!();
        println!(
            "== {} cells, {} steps, {:?} per step ==",
            workload.label(),
            profile.steps(),
            profile.per_step()
        );
        println!("   phase                       share    per step");
        for phase in STEP_PHASES {
            println!(
                "   {:<24} {:>7.1}%   {:>10?}",
                phase.label(),
                100.0 * profile.share(phase),
                profile.per_step_in(phase)
            );
        }

        let mut terms = TermProfiler::of_scenario(&scenario);
        let state = workload.benchmark_state();
        let wind_stress = workload.wind_stress();
        let mut tendency = workload.tendency_buffer();
        for _ in 0..PROFILED_EVALUATIONS {
            terms.evaluate(&state, &wind_stress, &mut tendency);
        }
        let term_profile = terms.profile();

        println!();
        println!(
            "   the two evaluator phases are {:.1}% of a step; term by term, over {} \
             evaluations:",
            100.0
                * (profile.share(StepPhase::ShallowWaterTerms)
                    + profile.share(StepPhase::Coriolis)),
            term_profile.evaluations()
        );
        println!("   kernel                       phase    share    per evaluation");
        for term in RHS_TERMS {
            println!(
                "   {:<26} {:>7}  {:>6.1}%   {:>10?}",
                term.label(),
                match term.phase() {
                    StepPhase::Coriolis => "rot",
                    _ => "sw",
                },
                100.0 * term_profile.share(term),
                term_profile.per_evaluation(term)
            );
        }
        println!(
            "   kernels timed one at a time total {:?} per evaluation, against {:?} for the \
             two evaluators left whole.",
            term_profile.per_evaluation_total(),
            (profile.per_step_in(StepPhase::ShallowWaterTerms)
                + profile.per_step_in(StepPhase::Coriolis))
                / RK4_STAGES_PER_STEP
        );
    }
}

/// Take uninstrumented steps of the finest benchmark workload for `duration`,
/// so that an external sampler has a hot path to sample.
///
/// The real [`Solver`], with no clock between its phases: what a sampler sees
/// here is the step a run takes.
fn spin(duration: Duration) {
    let workload = BENCHMARK_WORKLOADS
        .last()
        .expect("the benchmark suite defines at least one workload");
    let scenario = workload.scenario();
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
    .expect("a benchmark workload is a runnable scenario");
    let mut state = OceanState::at_rest(basin.grid());

    println!(
        "spinning {} for {duration:?}; pid {}",
        workload.label(),
        std::process::id()
    );
    let started = Instant::now();
    let mut steps: u64 = 0;
    while started.elapsed() < duration {
        // The clock is read once per batch rather than once per step, so the
        // loop around the step is negligible against the step even at the
        // coarser grid.
        for _ in 0..32 {
            solver.step_forced_by(&mut state, steps as f64 * schedule.dt_s(), basin, &wind);
            steps += 1;
        }
    }
    println!(
        "{steps} steps in {:?} — {:.0} steps/s",
        started.elapsed(),
        steps as f64 / started.elapsed().as_secs_f64()
    );
}
