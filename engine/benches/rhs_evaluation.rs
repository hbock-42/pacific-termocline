//! T-10.1 — how fast the engine evaluates the shallow-water right-hand side.
//!
//! One [`ShallowWaterRhs::evaluate`] over a whole basin: the pressure
//! gradient, the surface stress, the Rayleigh damping and the continuity
//! divergence of Epic 02. This is the innermost thing a run does — RK4 calls
//! it four times per timestep — so it is where a change to the hot path shows
//! up first and least diluted.
//!
//! **The figure this reports is grid cells per second.** Criterion is given
//! the workload's cell count as the element count of one iteration, so its
//! `thrpt` line is cells per second directly, and the two resolutions are
//! comparable to each other despite differing in size by a factor of four.
//!
//! Everything the evaluation touches — the evaluator, the state, the stress
//! field and the tendency buffer — is built before the timing loop, which is
//! how a run holds them too (CODING_STANDARDS.md § *Performance*). What is
//! measured is one evaluation and nothing around it.
//!
//! The inputs come from `engine::benchmark`, which is where what they are and
//! why is written down, and which `tests/benchmark_workloads.rs` holds to it.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

use engine::benchmark::BENCHMARK_WORKLOADS;

fn rhs_evaluation(criterion: &mut Criterion) {
    // Criterion's throughput line says `elem/s` whatever the element is, and a
    // reader of a CI job summary has no other way to find out which.
    println!(
        "rhs_evaluation: one element is one grid cell, so `thrpt` reads as grid cells per second."
    );
    let mut group = criterion.benchmark_group("rhs_evaluation");
    for workload in BENCHMARK_WORKLOADS {
        let mut rhs = workload.rhs_evaluator();
        let state = workload.benchmark_state();
        let wind_stress = workload.wind_stress();
        let mut tendency = workload.tendency_buffer();

        // One "element" is one grid cell, so criterion's throughput line reads
        // as grid cells per second.
        group.throughput(Throughput::Elements(workload.grid_cells()));
        group.bench_function(BenchmarkId::from_parameter(workload.label()), |bencher| {
            bencher.iter(|| {
                rhs.evaluate(
                    black_box(&state),
                    black_box(&wind_stress),
                    black_box(&mut tendency),
                );
            });
        });
    }
    group.finish();
}

criterion_group!(benches, rhs_evaluation);
criterion_main!(benches);
