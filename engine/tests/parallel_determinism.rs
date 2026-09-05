//! Acceptance tests for T-10.3 — that threading the array sweeps does not
//! change a single bit of a run.
//!
//! The ticket's acceptance criteria are two, and this file is the second:
//! *"existing validation tests (Epic 07) still pass bit-for-bit … parallelism
//! must not change results beyond floating-point summation-order noise."* The
//! first — a meaningful speed-up — is a measurement rather than an assertion,
//! and it is in `docs/performance-notes.md`, for the reason
//! `docs/benchmarks.md` gives for keeping wall-clock thresholds out of the
//! gate.
//!
//! **Bit for bit, not within summation-order noise.** The criterion allows the
//! looser bound, and this suite does not take it, because the right-hand side
//! has no summation order to disturb: it contains no reduction of any kind.
//! Every one of its fourteen kernels is a map — each output point is written
//! once, from a fixed handful of input points, in a fixed order — so splitting
//! the sweep across rows moves *which thread* performs an operation and
//! nothing else. A tolerance here would be a place for a real divergence to
//! hide, which is `tests/step_profile.rs`'s and
//! `tests/wind_stress_cache.rs`'s reasoning for comparing the same way.
//!
//! The thread counts are chosen to break a scheme that only worked by
//! accident. Ten is the machine's; two is the smallest real split; **three is
//! the important one**, because it divides neither the 50 rows of the coarse
//! benchmark basin nor the 100 of the fine one, so a sweep whose last chunk
//! were handled differently from the rest would show here and nowhere else.
//! One thread exercises the pool's own path against no pool at all.
//!
//! `CODING_STANDARDS.md` § *Correctness and failure* is what this is holding:
//! identical scenario in, byte-identical output — which has to survive a
//! machine deciding for itself how many workers to wake.
//!
//! What this file does **not** check is that the parallel kernels compute the
//! right thing at all: a kernel wrong in the same way on every thread agrees
//! with itself, and every reference stepper in this suite calls the same
//! evaluators. What catches that is Epic 07 — the wave, dispersion and
//! conservation validations, which compare against analytic results rather
//! than against another run — and they pass unchanged. The two are
//! complementary and neither is sufficient: Epic 07 would not notice a
//! thread-count dependence that stayed inside its tolerances, and this file
//! would not notice a kernel that is uniformly wrong.

use engine::benchmark::BENCHMARK_WORKLOADS;
use engine::{BetaPlane, CoriolisTerm, OceanState, Scenario, ShallowWaterRhs, Solver, WindForcing};

/// Thread counts every result below must be identical across.
///
/// Three is deliberately coprime to both benchmark basins' row counts; see
/// the module documentation.
const THREAD_COUNTS: [usize; 4] = [1, 2, 3, 10];

/// Steps a compared run takes.
///
/// Long enough that a divergence in any one kernel has been fed back through
/// the whole state several times over — a Kelvin wave crosses a few cells in
/// this many hours — and short enough to run four times in a gate.
const COMPARED_STEPS: u64 = 24;

/// Run `body` on a rayon pool of exactly `threads` workers.
///
/// `install` runs the closure *on* that pool, so every `write_rows` inside it
/// splits across those workers rather than across the global pool's.
fn with_threads<T: Send>(threads: usize, body: impl FnOnce() -> T + Send) -> T {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("a pool of a stated size can be built")
        .install(body)
}

/// The two benchmark workloads' scenarios: the control scenario at 1.0° and at
/// 0.5°, which are the basins the note's measurements are taken on.
fn compared_scenarios() -> Vec<Scenario> {
    BENCHMARK_WORKLOADS
        .iter()
        .map(engine::benchmark::BenchmarkWorkload::scenario)
        .collect()
}

/// One right-hand-side evaluation of `scenario`'s benchmark state, as the
/// three fields of the tendency it writes.
fn one_evaluation(scenario: &Scenario, workload_index: usize) -> OceanState {
    let workload = BENCHMARK_WORKLOADS[workload_index];
    let basin = scenario.basin();
    let mut rhs = ShallowWaterRhs::new(basin.grid(), basin.spacing(), scenario.physical_params());
    let mut coriolis = CoriolisTerm::new(
        basin.grid(),
        basin.spacing(),
        BetaPlane::of_basin(scenario.physical_params(), basin),
    );
    let state = workload.benchmark_state();
    let stress = workload.wind_stress();
    let mut tendency = OceanState::at_rest(basin.grid());
    rhs.evaluate(&state, &stress, &mut tendency);
    coriolis.add_to_tendency(&state, &mut tendency);
    tendency
}

/// `COMPARED_STEPS` steps of `scenario` from its benchmark state, as the state
/// they reach.
fn a_short_run(scenario: &Scenario, workload_index: usize) -> OceanState {
    let workload = BENCHMARK_WORKLOADS[workload_index];
    let basin = scenario.basin();
    let dt_s = scenario.output_schedule().dt_s();
    let mut solver = Solver::new(
        basin.grid(),
        basin.spacing(),
        scenario.physical_params(),
        BetaPlane::of_basin(scenario.physical_params(), basin),
        dt_s,
    )
    .expect("a benchmark scenario is a runnable one");
    let mut forcing = WindForcing::new(basin, scenario.wind());
    let mut state = workload.benchmark_state();
    for step in 0..COMPARED_STEPS {
        solver.step_with_forcing(&mut state, step as f64 * dt_s, &mut forcing);
    }
    state
}

#[test]
fn a_right_hand_side_evaluation_is_the_same_bits_on_any_number_of_threads() {
    for (index, scenario) in compared_scenarios().iter().enumerate() {
        let reference = with_threads(1, || one_evaluation(scenario, index));
        for threads in THREAD_COUNTS {
            let evaluated = with_threads(threads, || one_evaluation(scenario, index));
            assert_eq!(
                evaluated,
                reference,
                "the tendency over {} differs on {threads} threads; the right-hand side has no \
                 reduction, so a row split cannot change a value and this is a real divergence",
                scenario.basin().grid().nx(),
            );
        }
    }
}

#[test]
fn a_run_reaches_the_same_bits_on_any_number_of_threads() {
    for (index, scenario) in compared_scenarios().iter().enumerate() {
        let reference = with_threads(1, || a_short_run(scenario, index));
        for threads in THREAD_COUNTS {
            let reached = with_threads(threads, || a_short_run(scenario, index));
            assert_eq!(
                reached,
                reference,
                "{COMPARED_STEPS} steps over {} reach a different state on {threads} threads; \
                 CODING_STANDARDS.md requires identical scenario in, byte-identical output",
                scenario.basin().grid().nx(),
            );
        }
    }
}
