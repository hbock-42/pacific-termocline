//! T-10.1 — how fast the engine runs a scenario end to end.
//!
//! A whole short run through the same entry point the `run` command uses:
//! scenario build, solver construction, the time loop with its wind re-sampled
//! at each of RK4's four stages, and the run directory written. The
//! right-hand-side benchmark beside this one isolates the hot path; this one
//! is the number a user feels, and it is what says whether a speed-up in the
//! inner loop survives everything around it.
//!
//! **The figure this reports is timesteps per second.** Criterion is given the
//! run's step count as the element count of one iteration, so its `thrpt` line
//! is steps per second directly. Grid cells per second for a run is that
//! figure times the workload's cell count, which the right-hand-side benchmark
//! reports on its own.
//!
//! The run writes exactly two frames — the initial state and the final one —
//! so the filesystem is a constant term of the measurement rather than one
//! that grows with the run length. It is written into one directory per
//! workload, created before the timing loop and truncated by each run, so what
//! varies between iterations is the arithmetic rather than the number of files
//! on disk.
//!
//! The workloads — and the output directory they are written into — come from
//! `engine::benchmark`, which is where what they are and why is written down,
//! and which `tests/benchmark_workloads.rs` holds to it.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

use engine::benchmark::{BenchmarkOutputDir, BENCHMARK_WORKLOADS};

fn scenario_run(criterion: &mut Criterion) {
    // Criterion's throughput line says `elem/s` whatever the element is, and a
    // reader of a CI job summary has no other way to find out which. Printing
    // the legend beside the numbers is the cheapest place to say it.
    println!(
        "scenario_run: one element is one timestep, so `thrpt` reads as timesteps per second. \
         Grid cells per second for a run is that figure times the workload's cell count; the \
         rhs_evaluation benchmark reports cells per second directly."
    );
    let mut group = criterion.benchmark_group("scenario_run");
    // A run of a couple of hundred steps is far too long for criterion's
    // default hundred samples; ten is its minimum and enough for the
    // dispersion of a workload this deterministic.
    group.sample_size(10);
    for workload in BENCHMARK_WORKLOADS {
        let directory = BenchmarkOutputDir::new(&workload.label());

        // One "element" is one timestep, so criterion's throughput line reads
        // as timesteps per second.
        group.throughput(Throughput::Elements(workload.timesteps()));
        group.bench_function(BenchmarkId::from_parameter(workload.label()), |bencher| {
            bencher.iter(|| {
                workload
                    .run_into(directory.path())
                    .expect("a benchmark workload is a runnable scenario")
            });
        });
    }
    group.finish();
}

criterion_group!(benches, scenario_run);
criterion_main!(benches);
