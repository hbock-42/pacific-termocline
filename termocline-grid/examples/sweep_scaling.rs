//! How a row-split sweep scales with threads, against the size of the field
//! it sweeps (T-10.3).
//!
//! The engine's own measurement is in `docs/performance-notes.md`: threading
//! the fourteen array kernels of the right-hand side makes a step slower, and
//! worse with every core added. That measurement cannot say *why* on its own,
//! because two explanations fit it — the sweeps are too short to amortise a
//! fork and a join, or the computation is bandwidth-bound and extra cores are
//! extra claimants on one bus. The two call for different things next, so it
//! is worth separating them.
//!
//! This example separates them by holding the kernel fixed and growing the
//! field. It sweeps `c = a + b` over rows — three basin-sized arrays streamed
//! for one flop, which is the shape of every kernel the term table lists — at
//! the engine's own 320 x 100 and at sizes up to two hundred times larger, and
//! reports nanoseconds per point for a plain serial loop and for
//! [`sweep::write_rows`](termocline_grid::sweep::write_rows).
//!
//! If the parallel sweep beats the serial one once the field is big enough,
//! the cost at basin size is the fork and the join, and a design that
//! synchronised less often could still win. If it never beats it at any size,
//! the memory system is the limit and no arrangement of threads will help.
//!
//! ```sh
//! for n in 1 2 4 10; do RAYON_NUM_THREADS=$n cargo run --release -p termocline-grid --example sweep_scaling; done
//! ```

use std::time::{Duration, Instant};

use termocline_grid::sweep::write_rows;
use termocline_grid::Field2D;

/// Shapes swept, from the engine's 0.5 degree basin upward by factors of four.
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
    println!(
        "sweep_scaling: c = a + b over rows, {} rayon threads",
        rayon::current_num_threads()
    );
    println!(
        "   {:>14}  {:>11}  {:>11}  {:>8}",
        "shape", "serial", "rayon", "speed-up"
    );
    for (nx, ny) in SHAPES {
        let a = Field2D::filled(nx, ny, 1.0_f64).expect("a positive shape");
        let b = Field2D::filled(nx, ny, 2.0_f64).expect("a positive shape");
        let mut c = Field2D::filled(nx, ny, 0.0_f64).expect("a positive shape");

        let serial = per_point(&mut c, |field| {
            let points_per_row = field.nx();
            for (j, row) in field
                .as_mut_slice()
                .chunks_exact_mut(points_per_row)
                .enumerate()
            {
                let (left, right) = (a.row(j), b.row(j));
                for (value, (left, right)) in row.iter_mut().zip(left.iter().zip(right)) {
                    *value = left + right;
                }
            }
        });
        let parallel = per_point(&mut c, |field| {
            write_rows(field, |j, row| {
                let (left, right) = (a.row(j), b.row(j));
                for (value, (left, right)) in row.iter_mut().zip(left.iter().zip(right)) {
                    *value = left + right;
                }
            });
        });
        println!(
            "   {:>14}  {:>9.3} ns  {:>9.3} ns  {:>7.2}x",
            format!("{nx}x{ny}"),
            serial,
            parallel,
            serial / parallel
        );
    }
}

/// Nanoseconds per point of one `sweep` of `field`, repeated for
/// [`SECONDS_PER_MEASUREMENT`] after a warm-up of the same length.
///
/// A fixed *duration* rather than a fixed count, for the reason
/// `engine/examples/profile.rs` gives: a fixed count warms the large shapes
/// and leaves the small ones on a processor that has not ramped its clock.
fn per_point(field: &mut Field2D<f64>, mut sweep: impl FnMut(&mut Field2D<f64>)) -> f64 {
    let budget = Duration::from_secs_f64(SECONDS_PER_MEASUREMENT);
    let warm_up = Instant::now();
    while warm_up.elapsed() < budget {
        sweep(field);
    }
    let mut sweeps = 0_u64;
    let started = Instant::now();
    while started.elapsed() < budget {
        sweep(field);
        sweeps += 1;
    }
    started.elapsed().as_secs_f64() * 1e9 / (sweeps as f64 * field.len() as f64)
}
