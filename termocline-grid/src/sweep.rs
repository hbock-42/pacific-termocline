//! Row-blocked sweeps over a [`Field2D`], run across threads where the field
//! is big enough to pay for it.
//!
//! Every array kernel in the numerical core writes one output point from a
//! fixed handful of input points, and `docs/performance-notes.md` (T-10.2)
//! found fourteen of them, flat, with no hot kernel among them. There is
//! therefore nothing to parallelise *inside* a kernel; what T-10.3
//! parallelises is the sweep the kernels share, once, here.
//!
//! # Why a row is the unit
//!
//! A [`Field2D`] is row-major, so a row is a contiguous slice and rows are
//! disjoint. Splitting a sweep by rows therefore needs no synchronisation and,
//! more importantly, **changes no arithmetic**: each output point is written
//! by exactly one call, from the same inputs, in the same order, whatever
//! thread runs it. There is no reduction anywhere in the right-hand side — no
//! sum whose order a work-stealing scheduler could permute — which is what
//! lets the engine keep the bit-for-bit determinism CODING_STANDARDS.md
//! § *Correctness and failure* requires of it. A sweep of the same field
//! produces the same bytes on one thread and on ten, and
//! `engine/tests/parallel_determinism.rs` is what holds it to that.
//!
//! # Why there is a threshold
//!
//! Handing a sweep to a thread pool costs something, and a field small enough
//! makes that cost the whole measurement. [`SERIAL_POINT_LIMIT`] is where the
//! sweep stays on the calling thread instead.
//!
//! Neither it nor [`MIN_POINTS_PER_TASK`] is a measured optimum, and neither
//! claims to be: they are a guard chosen before the measurement, because the
//! measurement is what this ticket exists to take. What can be said about them
//! without one is that they cannot change a *result*, only its cost — both
//! paths perform identical arithmetic in identical order — so an arbitrary
//! value here is a performance question and never a correctness one. What the
//! measurement then said is in `docs/performance-notes.md` § *After T-10.3*,
//! and it is that no value of either rescues the idea.

use rayon::prelude::*;

use crate::Field2D;

/// Points at or below which a sweep runs on the calling thread rather than on
/// the thread pool.
///
/// A guard, not an optimum: it is set below the smallest basin the benchmark
/// suite measures (160 × 50, 8 000 cells) and above the few-hundred-point
/// basins the unit tests integrate, so that a suite of tiny grids does not
/// spend its time in a thread pool. Where the real crossover is, is a question
/// for the measurement rather than for this constant — and
/// `docs/performance-notes.md` § *After T-10.3* is where it was answered.
pub const SERIAL_POINT_LIMIT: usize = 2_048;

/// Fewest points a single worker is handed.
///
/// Rayon splits until every task is at least this big, so a 100-row field does
/// not become 100 tasks of one row each. Like [`SERIAL_POINT_LIMIT`] it is a
/// starting value rather than a measured one; the note records what happened
/// when it was varied, including at the coarsest split there is — one task per
/// worker per sweep — which was no better.
const MIN_POINTS_PER_TASK: usize = 1_024;

/// Call `write_row(j, row)` for every row `j` of `field`, in parallel across
/// rows once the field is larger than [`SERIAL_POINT_LIMIT`].
///
/// `row` is the `j`th row's `nx` values, contiguous and mutable; the closure
/// owns it exclusively, and reads whatever input rows its kernel needs through
/// the shared references it captures. Rows are disjoint, so nothing here
/// synchronises and nothing accumulates across rows.
pub fn write_rows<T, F>(field: &mut Field2D<T>, write_row: F)
where
    T: Send,
    F: Fn(usize, &mut [T]) + Send + Sync,
{
    let points_per_row = field.nx();
    if field.len() <= SERIAL_POINT_LIMIT {
        for (j, row) in field
            .as_mut_slice()
            .chunks_exact_mut(points_per_row)
            .enumerate()
        {
            write_row(j, row);
        }
        return;
    }
    // `div_ceil(..).max(1)` rather than a bare division: a field wider than
    // `MIN_POINTS_PER_TASK` still gives every worker at least one row.
    let min_rows_per_task = MIN_POINTS_PER_TASK.div_ceil(points_per_row).max(1);
    field
        .as_mut_slice()
        .par_chunks_exact_mut(points_per_row)
        .with_min_len(min_rows_per_task)
        .enumerate()
        .for_each(|(j, row)| write_row(j, row));
}
