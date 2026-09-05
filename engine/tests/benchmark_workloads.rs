//! Acceptance tests for T-10.1 — the workloads the `cargo bench` suite
//! measures, and the figures it divides their timings by.
//!
//! The benchmarks themselves live in `benches/`, where `cargo test` cannot
//! assert on them: criterion reports a distribution, not a value, and a
//! wall-clock number is not something a test can pin down. What *is* testable
//! is everything the reported figure is built out of — that each workload is a
//! scenario the engine will actually run, that it is the control scenario's
//! physics at a stated resolution rather than a benchmark-only invention, that
//! two runs of one workload produce byte-identical output, and that the
//! element counts criterion divides by are the work the benchmark performed.
//! Those are the acceptance criteria of the ticket ("benchmarks run
//! reproducibly and report timestep-per-second / grid-cells-per-second
//! figures"); the timings are the report, and these are the guarantees the
//! report rests on.
//!
//! Nothing here asserts a duration. A benchmark harness whose *tests* were
//! timing-sensitive would be exactly the flaky measurement this ticket exists
//! to avoid.

use std::fs;

use engine::benchmark::{
    BenchmarkWorkload, BENCHMARK_ANOMALY_AMPLITUDE_M, BENCHMARK_ANOMALY_WIDTH_M,
    BENCHMARK_WORKLOADS, SHORT_RUN_STEPS,
};
use engine::{BetaPlane, Scenario, ShallowWaterRhs, Solver, FRAME_FILE_NAME, HEADER_FILE_NAME};

use engine::benchmark::BenchmarkOutputDir;

/// The control scenario the benchmark workloads are cut down from, as it sits
/// on disk. Read through the filesystem rather than through
/// [`engine::benchmark`], so that the assertion below compares two independent
/// readings of it rather than one reading with itself.
fn control_scenario() -> Scenario {
    Scenario::load(std::path::Path::new("scenarios/steady-trades.toml"))
        .expect("the control scenario ships with the engine and is valid")
}

/// The zonal span of the equatorial Pacific basin, in degrees: 120°E to 80°W
/// (CONTEXT.md, *Basin*).
const PACIFIC_ZONAL_SPAN_DEG: f64 = 160.0;

/// The meridional span of the equatorial Pacific basin, in degrees: 25°S to
/// 25°N (CONTEXT.md, *Basin*).
const PACIFIC_MERIDIONAL_SPAN_DEG: f64 = 50.0;

#[test]
fn every_workload_is_a_scenario_the_engine_will_run() {
    // A benchmark whose scenario the engine refuses does not measure anything:
    // both timestep bounds have to be clear at every resolution in the suite,
    // or the harness reports on a run that cannot happen.
    for workload in BENCHMARK_WORKLOADS {
        let scenario = workload.scenario();
        let basin = scenario.basin();
        let params = scenario.physical_params();
        let plane = BetaPlane::of_basin(params, basin);
        Solver::new(
            basin.grid(),
            basin.spacing(),
            params,
            plane,
            scenario.output_schedule().dt_s(),
        )
        .unwrap_or_else(|error| {
            panic!(
                "the {} workload asks for a timestep the scheme refuses: {error}",
                workload.label()
            )
        });
    }
}

#[test]
fn each_workload_covers_the_grid_its_resolution_implies() {
    // The expected cell counts come from the basin's span in degrees
    // (CONTEXT.md, *Basin*) divided by the workload's resolution, not from
    // running the code: a 1.0° cut of 160° × 50° is 160 × 50 cells, and a 0.5°
    // cut is 320 × 100.
    for workload in BENCHMARK_WORKLOADS {
        let expected_nx = (PACIFIC_ZONAL_SPAN_DEG / workload.resolution_deg()) as usize;
        let expected_ny = (PACIFIC_MERIDIONAL_SPAN_DEG / workload.resolution_deg()) as usize;
        let grid = workload.grid();

        assert_eq!(
            (grid.nx(), grid.ny()),
            (expected_nx, expected_ny),
            "the {} workload's grid is not the {}° cut of the Pacific basin",
            workload.label(),
            workload.resolution_deg()
        );
        assert_eq!(
            workload.grid_cells(),
            (expected_nx * expected_ny) as u64,
            "the cell count criterion divides by is not this workload's grid"
        );
    }
}

#[test]
fn the_suite_measures_more_than_one_resolution() {
    // "A couple of representative grid resolutions": a single point tells you
    // a duration, two tell you how the engine scales with the basin, which is
    // what the optimisation tickets behind this one need.
    assert!(
        BENCHMARK_WORKLOADS.len() >= 2,
        "a benchmark suite at one resolution cannot show how cost scales with the grid"
    );
    let mut resolutions: Vec<f64> = BENCHMARK_WORKLOADS
        .iter()
        .map(BenchmarkWorkload::resolution_deg)
        .collect();
    resolutions.dedup();
    assert_eq!(
        resolutions.len(),
        BENCHMARK_WORKLOADS.len(),
        "two workloads at the same resolution measure the same thing twice"
    );
}

#[test]
fn a_workload_is_the_control_scenario_at_another_resolution() {
    // The suite measures the physics the project actually runs. If a benchmark
    // could drift into its own g', H, r or wind, a later optimisation would be
    // measured against a scenario nobody simulates.
    let control = control_scenario();
    for workload in BENCHMARK_WORKLOADS {
        let scenario = workload.scenario();
        assert_eq!(
            scenario.physical_params(),
            control.physical_params(),
            "the {} workload's ocean is not the control scenario's",
            workload.label()
        );
        assert_eq!(
            scenario.winds(),
            control.winds(),
            "the {} workload's forcing is not the control scenario's",
            workload.label()
        );
        assert_eq!(
            scenario.output_schedule().dt_s(),
            control.output_schedule().dt_s(),
            "the {} workload does not step at the control scenario's timestep",
            workload.label()
        );
    }
}

#[test]
fn the_reference_workload_runs_the_control_scenarios_own_grid() {
    // One workload has to be the basin the project's own scenarios are written
    // for; otherwise every figure the suite reports is about a grid nobody
    // runs.
    let control_grid = control_scenario().basin().grid();
    assert!(
        BENCHMARK_WORKLOADS
            .iter()
            .any(|workload| workload.grid() == control_grid),
        "no workload runs the control scenario's {}x{} grid",
        control_grid.nx(),
        control_grid.ny()
    );
}

#[test]
fn the_short_run_is_short_but_not_a_single_step() {
    // The run benchmark trades length for turnaround, and both directions have
    // a floor: a run of a handful of steps measures the solver's construction
    // rather than its loop, and a run of the control scenario's full length
    // (17 520 steps) cannot be sampled repeatedly.
    let control_steps = control_scenario().output_schedule().total_steps();
    for workload in BENCHMARK_WORKLOADS {
        let steps = workload.timesteps();
        assert_eq!(
            steps,
            SHORT_RUN_STEPS,
            "the {} workload does not take the number of steps the suite states",
            workload.label()
        );
        assert!(
            steps >= 100,
            "a run of {steps} steps is dominated by everything that happens before the first one"
        );
        assert!(
            steps < control_steps,
            "the benchmark run is not shorter than the {control_steps}-step scenario it is cut \
             down from"
        );
        assert_eq!(
            workload.scenario().output_schedule().total_steps(),
            steps,
            "the step count criterion divides by is not the number of steps the run takes"
        );
    }
}

#[test]
fn the_run_benchmark_reports_the_steps_it_actually_took() {
    // The timestep-per-second figure is the workload's step count over the
    // measured duration, so the run has to take exactly that many steps —
    // criterion divides by a number this test is what ties to reality.
    for workload in BENCHMARK_WORKLOADS {
        let directory = BenchmarkOutputDir::new(&format!("steps-taken-{}", workload.label()));
        let report = workload
            .run_into(directory.path())
            .expect("a benchmark workload is a runnable scenario");

        assert_eq!(
            report.steps_taken(),
            workload.timesteps(),
            "the {} workload took a different number of steps from the one criterion divides by",
            workload.label()
        );
    }
}

#[test]
fn two_runs_of_one_workload_write_byte_identical_output() {
    // "Benchmarks run reproducibly": the same workload run twice is the same
    // arithmetic in the same order (CODING_STANDARDS.md § *Correctness and
    // failure*), so any difference a later measurement shows is the change
    // being measured rather than the workload wandering.
    let workload = BENCHMARK_WORKLOADS[0];
    let first = BenchmarkOutputDir::new("reproducible-a");
    let second = BenchmarkOutputDir::new("reproducible-b");

    workload.run_into(first.path()).expect("the run succeeds");
    workload.run_into(second.path()).expect("the run succeeds");

    for file in [HEADER_FILE_NAME, FRAME_FILE_NAME] {
        assert_eq!(
            fs::read(first.path().join(file)).expect("the run wrote its header and frames"),
            fs::read(second.path().join(file)).expect("the run wrote its header and frames"),
            "two runs of the {} workload disagree about {file}",
            workload.label()
        );
    }
}

#[test]
fn the_right_hand_side_benchmark_evaluates_the_same_tendency_every_time() {
    // The other half of reproducibility: the RHS benchmark's inputs are built
    // from an analytic initial state rather than from a random or a
    // time-dependent one, so every iteration of the benchmark evaluates the
    // identical arithmetic.
    for workload in BENCHMARK_WORKLOADS {
        let state = workload.benchmark_state();
        let stress = workload.wind_stress();

        let mut rhs = workload.rhs_evaluator();
        let mut first = workload.tendency_buffer();
        rhs.evaluate(&state, &stress, &mut first);

        let mut second_rhs = workload.rhs_evaluator();
        let mut second = workload.tendency_buffer();
        second_rhs.evaluate(
            &workload.benchmark_state(),
            &workload.wind_stress(),
            &mut second,
        );

        assert_eq!(
            first.h().as_slice(),
            second.h().as_slice(),
            "the {} workload's ∂h/∂t is not reproducible",
            workload.label()
        );
        assert_eq!(first.u().as_slice(), second.u().as_slice());
        assert_eq!(first.v().as_slice(), second.v().as_slice());
    }
}

#[test]
fn the_benchmark_state_is_the_analytic_kelvin_structure_not_an_ocean_at_rest() {
    // A state at rest has a tendency of zeros, which is not the arithmetic a
    // run does. The benchmark state is the equatorial Kelvin structure of
    // CONTEXT.md (*Kelvin wave*): `h` a Gaussian of amplitude `A` and zonal
    // width `W` times `e^{−η²/2}` with `η = y/Le`, and `u = (c/H)·h` — the
    // balance that makes it a state a run could be in rather than an arbitrary
    // field.
    //
    // Both peaks are asserted against `A` and `(c/H)·A`, which are analytic
    // constants of the structure rather than anything read out of a run. The
    // *lower* bound is the sampling loss: a grid samples the Gaussian at cell
    // positions, and the worst case puts the analytic peak half a cell away in
    // each direction, so the largest sampled value cannot fall below
    // `A · exp(−(dx/2W)² − (dy/2Le)²/2)`. The upper bound is `A` itself, which
    // no sample can exceed.
    for workload in BENCHMARK_WORKLOADS {
        let scenario = workload.scenario();
        let params = scenario.physical_params();
        let basin = scenario.basin();
        let state = workload.benchmark_state();

        let deformation_radius_m =
            (params.kelvin_wave_speed_m_per_s() / params.beta_per_m_per_s()).sqrt();
        let zonal_loss = (basin.spacing().dx_m() / (2.0 * BENCHMARK_ANOMALY_WIDTH_M)).powi(2);
        let meridional_loss = (basin.spacing().dy_m() / (2.0 * deformation_radius_m)).powi(2) / 2.0;
        let smallest_sampled_fraction = (-zonal_loss - meridional_loss).exp();

        let largest_h_m = largest_magnitude(state.h().as_slice());
        assert!(
            largest_h_m <= BENCHMARK_ANOMALY_AMPLITUDE_M
                && largest_h_m >= BENCHMARK_ANOMALY_AMPLITUDE_M * smallest_sampled_fraction,
            "the {} workload's thermocline anomaly peaks at {largest_h_m} m, not at the \
             {BENCHMARK_ANOMALY_AMPLITUDE_M} m Gaussian the benchmark state is defined as",
            workload.label()
        );

        // The Kelvin balance: the same structure scaled by `c/H`, so its peak
        // is the peak of `h` scaled by `c/H` under the identical sampling.
        let balanced_current_m_per_s = BENCHMARK_ANOMALY_AMPLITUDE_M
            * params.kelvin_wave_speed_m_per_s()
            / params.mean_thermocline_depth_m();
        let largest_u_m_per_s = largest_magnitude(state.u().as_slice());
        assert!(
            largest_u_m_per_s <= balanced_current_m_per_s
                && largest_u_m_per_s >= balanced_current_m_per_s * smallest_sampled_fraction,
            "the {} workload's zonal current peaks at {largest_u_m_per_s} m/s, which is not in \
             Kelvin balance with a {BENCHMARK_ANOMALY_AMPLITUDE_M} m anomaly \
             ({balanced_current_m_per_s} m/s)",
            workload.label()
        );
    }
}

/// The largest magnitude in `values` — the peak of a field, whichever sign it
/// is on.
fn largest_magnitude(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0_f64, |largest, value| largest.max(value.abs()))
}

#[test]
fn the_right_hand_side_benchmark_evaluates_over_the_grid_it_reports() {
    // The grid-cells-per-second figure is the workload's cell count over the
    // measured duration. The cell count is the *centres* — `h` and the
    // divergence live there — and the staggered velocity fields carry the one
    // extra column and row the C-grid puts on their own faces, which is not a
    // second grid.
    for workload in BENCHMARK_WORKLOADS {
        let state = workload.benchmark_state();
        let grid = workload.grid();
        assert_eq!(state.grid(), grid);
        assert_eq!(workload.grid_cells(), (grid.nx() * grid.ny()) as u64);
        assert_eq!(state.h().len(), workload.grid_cells() as usize);
    }
}

#[test]
fn a_workloads_evaluator_and_buffers_fit_its_own_state() {
    // `ShallowWaterRhs::evaluate` panics on a shape mismatch, so a benchmark
    // built out of mismatched pieces would fail as a panic mid-measurement
    // rather than as a test.
    for workload in BENCHMARK_WORKLOADS {
        let scenario = workload.scenario();
        let mut rhs = workload.rhs_evaluator();
        let state = workload.benchmark_state();
        let stress = workload.wind_stress();
        let mut tendency = workload.tendency_buffer();
        rhs.evaluate(&state, &stress, &mut tendency);

        // And the evaluator is the one the scenario's own solver would build,
        // not a differently parameterised one.
        let mut expected_rhs = ShallowWaterRhs::new(
            scenario.basin().grid(),
            scenario.basin().spacing(),
            scenario.physical_params(),
        );
        let mut expected_tendency = workload.tendency_buffer();
        expected_rhs.evaluate(&state, &stress, &mut expected_tendency);
        assert_eq!(tendency.h().as_slice(), expected_tendency.h().as_slice());
        assert_eq!(tendency.u().as_slice(), expected_tendency.u().as_slice());
        assert_eq!(tendency.v().as_slice(), expected_tendency.v().as_slice());
    }
}
