//! What halving the width of a field buys a sweep over it, against the size of
//! the field (T-10.4).
//!
//! `docs/performance-notes.md` concludes that the engine's right-hand side is
//! **memory-bandwidth-bound**: fourteen flat array kernels, none hot, each
//! costing about what the traffic it streams costs rather than what its
//! handful of flops costs. T-10.3 measured the corollary — adding cores to a
//! saturated bus makes a step slower — and left the other corollary open: if
//! traffic is what is being paid for, storing the fields as `f32` halves it,
//! and should buy about `2x` on the phases that stream.
//!
//! *Should* is the word this example exists to replace. It holds the kernel
//! fixed and varies the width: it sweeps `c = a + b` over three arrays for one
//! flop — the shape of every kernel in the term table — at the engine's own
//! 320 x 100 and at sizes up to two hundred times larger, in `f32` and in
//! `f64`, and reports nanoseconds per point for each.
//!
//! It is a **ceiling**, not a speed-up. The narrowing it measures is the best
//! an `f32` field layout could do to a phase that does nothing but stream, on
//! a kernel with no boundary handling, no interpolation and no scratch. What
//! an engine built on it would actually get is at most this, and the accuracy
//! question — `engine/tests/f32_field_storage.rs` — is a separate one that
//! this example says nothing about.
//!
//! ```sh
//! cargo run --release -p termocline-grid --example width_scaling
//! ```

use std::time::{Duration, Instant};

use termocline_grid::Field2D;

/// Shapes swept, from the engine's 0.5 degree basin upward by factors of four.
///
/// The same ladder `docs/performance-notes.md` reports T-10.3's thread scaling
/// on, so the two tables are read against each other: the first row is the
/// basin the benchmark suite runs, and the rest say whether what happens there
/// is a property of the machine's caches or of the kernel.
const SHAPES: [(usize, usize); 5] = [
    (320, 100),
    (640, 200),
    (1_280, 400),
    (2_560, 800),
    (5_120, 1_600),
];

/// Seconds each measurement is repeated for, and the same again as warm-up.
const SECONDS_PER_MEASUREMENT: f64 = 0.5;

fn main() {
    println!("width_scaling: c = a + b over a whole field, f64 against f32");
    println!(
        "   {:>14}  {:>12}  {:>12}  {:>8}",
        "shape", "f64", "f32", "narrowing"
    );
    for (nx, ny) in SHAPES {
        let wide = per_point_f64(nx, ny);
        let narrow = per_point_f32(nx, ny);
        println!(
            "   {:>14}  {:>10.4} ns  {:>10.4} ns  {:>7.2}x",
            format!("{nx}x{ny}"),
            wide,
            narrow,
            wide / narrow
        );
    }
}

/// Nanoseconds per point of `c = a + b` over an `nx` by `ny` field of `f64`.
fn per_point_f64(nx: usize, ny: usize) -> f64 {
    let a = Field2D::filled(nx, ny, 1.0_f64).expect("a positive shape");
    let b = Field2D::filled(nx, ny, 2.0_f64).expect("a positive shape");
    let mut c = Field2D::filled(nx, ny, 0.0_f64).expect("a positive shape");
    per_point(c.len(), || {
        for (value, (left, right)) in c
            .as_mut_slice()
            .iter_mut()
            .zip(a.as_slice().iter().zip(b.as_slice()))
        {
            *value = left + right;
        }
    })
}

/// Nanoseconds per point of `c = a + b` over an `nx` by `ny` field of `f32` —
/// the same sweep over the same number of points, streaming half the bytes.
fn per_point_f32(nx: usize, ny: usize) -> f64 {
    let a = Field2D::filled(nx, ny, 1.0_f32).expect("a positive shape");
    let b = Field2D::filled(nx, ny, 2.0_f32).expect("a positive shape");
    let mut c = Field2D::filled(nx, ny, 0.0_f32).expect("a positive shape");
    per_point(c.len(), || {
        for (value, (left, right)) in c
            .as_mut_slice()
            .iter_mut()
            .zip(a.as_slice().iter().zip(b.as_slice()))
        {
            *value = left + right;
        }
    })
}

/// Nanoseconds per point of one `sweep` over `points`, repeated for
/// [`SECONDS_PER_MEASUREMENT`] after a warm-up of the same length.
///
/// A fixed *duration* rather than a fixed count, for the reason
/// `engine/examples/profile.rs` gives: a fixed count warms the large shapes
/// and leaves the small ones on a processor that has not ramped its clock.
fn per_point(points: usize, mut sweep: impl FnMut()) -> f64 {
    let budget = Duration::from_secs_f64(SECONDS_PER_MEASUREMENT);
    let warm_up = Instant::now();
    while warm_up.elapsed() < budget {
        sweep();
    }
    let mut sweeps = 0_u64;
    let started = Instant::now();
    while started.elapsed() < budget {
        sweep();
        sweeps += 1;
    }
    started.elapsed().as_secs_f64() * 1e9 / (sweeps as f64 * points as f64)
}
