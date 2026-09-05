//! The memory half of T-05.3's acceptance criteria: reading a run must not
//! cost memory proportional to the run's length.
//!
//! This is its own test binary because it installs a counting global
//! allocator, and a peak-allocation measurement is only meaningful if nothing
//! else is allocating alongside it. It holds one test for the same reason.
//!
//! # What is measured, and against what
//!
//! Two runs of the same grid are read to the end, one [`SHORT_FRAMES`] frames
//! long and one [`LONG_FRAMES`], and the peak live heap during each iteration
//! is recorded. Neither number is compared against a threshold picked by
//! running the code:
//!
//! - the long run's peak must stay below the bytes its frames occupy, which is
//!   what a reader that decoded the run into a `Vec<Frame>` would need. That
//!   bound is arithmetic on the fixture: frame count times encoded frame size.
//! - the two peaks must agree to within [`SLACK_FRAMES`] frames' worth of
//!   heap, which is the criterion itself — a thousandfold longer run must not
//!   cost more memory.
//!
//! Neither byte source is a buffer holding the whole run: [`RepeatingFrames`]
//! synthesizes the frame stream one block at a time, exactly as an HTTP body
//! arriving in a browser would (ADR-0006). A fixture that materialized 8192
//! frames up front would swamp the measurement with its own allocation.

use std::alloc::{GlobalAlloc, Layout, System};
use std::io::{Cursor, Read};
use std::sync::atomic::{AtomicUsize, Ordering};

use termocline_format::{
    frame_encoding, BasinExtent, Frame, GridSpec, OutputTiming, PhysicalParams, RunHeader,
    RunReader, Variable,
};

/// Live heap bytes, and the high-water mark since it was last reset.
static LIVE_BYTES: AtomicUsize = AtomicUsize::new(0);
static PEAK_BYTES: AtomicUsize = AtomicUsize::new(0);

/// The system allocator, keeping a running total of live bytes.
struct Counting;

// SAFETY-free by construction: every method forwards to `System` unchanged and
// only adds two atomic counters, so the allocator's contract is `System`'s.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            record_allocation(layout.size());
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        LIVE_BYTES.fetch_sub(layout.size(), Ordering::Relaxed);
        unsafe { System.dealloc(pointer, layout) };
    }

    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let new_pointer = unsafe { System.realloc(pointer, layout, new_size) };
        if !new_pointer.is_null() {
            LIVE_BYTES.fetch_sub(layout.size(), Ordering::Relaxed);
            record_allocation(new_size);
        }
        new_pointer
    }
}

fn record_allocation(size: usize) {
    let live = LIVE_BYTES.fetch_add(size, Ordering::Relaxed) + size;
    PEAK_BYTES.fetch_max(live, Ordering::Relaxed);
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

/// Cells of the test basin. Large enough that one frame is tens of kilobytes,
/// so the measurement is dominated by frame data rather than by bookkeeping.
const NX: usize = 40;
const NY: usize = 20;

/// The short run: long enough that the reader has reached its steady state.
const SHORT_FRAMES: u64 = 8;
/// The long run, a thousandfold longer. At this grid it is a few hundred
/// megabytes of frames — more than the peak this test allows by three orders
/// of magnitude, so a reader that buffered the run could not pass.
const LONG_FRAMES: u64 = 8_192;
/// Frames' worth of heap the two peaks may differ by. The reader holds one
/// decoded frame at a time and grows the next frame's buffers beside it, so a
/// couple of frames of headroom covers the transient without covering any
/// growth that scales with run length.
const SLACK_FRAMES: usize = 3;

/// Model time between frames, in seconds: one day of output cadence.
const INTERVAL_S: f64 = 86_400.0;

fn grid() -> GridSpec {
    // The equatorial Pacific basin of CONTEXT.md: 120°E-80°W, 25°S-25°N.
    GridSpec::new(NX, NY, BasinExtent::new(120.0, -80.0, -25.0, 25.0))
        .expect("a 40x20 basin is a valid grid")
}

fn header(frame_count: u64) -> RunHeader {
    // Scenario values, not physical constants: nothing here is integrated.
    let params = PhysicalParams {
        mean_depth_m: 150.0,
        reduced_gravity_m_per_s2: 0.05,
        beta_per_m_per_s: 2.28e-11,
        rayleigh_damping_per_s: 1.0e-6,
        reference_density_kg_per_m3: 1025.0,
    };
    RunHeader::new(
        grid(),
        params,
        "a long run, synthesized frame by frame",
        OutputTiming {
            frame_count,
            interval_s: INTERVAL_S,
        },
    )
}

/// One encoded frame, the block the run is built out of.
fn encoded_frame() -> Vec<u8> {
    let g = grid();
    let field = |variable: Variable| vec![1.0_f64; g.field_len(variable)];
    let frame = Frame::new(
        0.0,
        &g,
        field(Variable::ThermoclineDepthAnomaly),
        field(Variable::ZonalCurrentAnomaly),
        field(Variable::MeridionalCurrentAnomaly),
        field(Variable::ZonalWindStress),
        field(Variable::MeridionalWindStress),
    )
    .expect("every field is built at the length the grid asks for");
    bincode::serde::encode_to_vec(&frame, frame_encoding()).expect("the frame encodes")
}

/// A frame stream that never exists in full: one encoded frame, repeated
/// `remaining` times, handed out in whatever sized bites the reader asks for.
struct RepeatingFrames {
    block: Vec<u8>,
    position: usize,
    remaining: u64,
}

impl RepeatingFrames {
    fn new(block: Vec<u8>, frames: u64) -> Self {
        Self {
            block,
            position: 0,
            remaining: frames,
        }
    }
}

impl Read for RepeatingFrames {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        if self.remaining == 0 {
            return Ok(0);
        }
        let available = &self.block[self.position..];
        let taken = available.len().min(buffer.len());
        buffer[..taken].copy_from_slice(&available[..taken]);
        self.position += taken;
        if self.position == self.block.len() {
            self.position = 0;
            self.remaining -= 1;
        }
        Ok(taken)
    }
}

/// Peak heap bytes allocated while reading a `frame_count`-frame run to the
/// end, over and above what was already live when the reading started.
fn peak_bytes_reading(frame_count: u64, block: &[u8]) -> usize {
    let header_bytes = serde_json::to_vec(&header(frame_count)).expect("the header serializes");
    let frames = RepeatingFrames::new(block.to_vec(), frame_count);
    let reader = RunReader::new(Cursor::new(header_bytes), frames)
        .expect("the synthesized run has a readable header");

    let baseline = LIVE_BYTES.load(Ordering::Relaxed);
    PEAK_BYTES.store(baseline, Ordering::Relaxed);

    // The frames are consumed rather than collected: a `Vec<Frame>` here would
    // measure the test's appetite, not the reader's.
    let mut read = 0_u64;
    let mut checksum = 0.0_f64;
    for frame in reader {
        let frame = frame.expect("every synthesized frame decodes");
        checksum += frame.h()[0];
        read += 1;
    }
    assert_eq!(read, frame_count, "the whole run was read");
    assert_eq!(checksum, frame_count as f64, "every frame was looked at");

    PEAK_BYTES.load(Ordering::Relaxed) - baseline
}

#[test]
fn reading_a_run_costs_memory_that_does_not_grow_with_its_length() {
    let block = encoded_frame();
    let frame_bytes = block.len();

    let short_peak = peak_bytes_reading(SHORT_FRAMES, &block);
    let long_peak = peak_bytes_reading(LONG_FRAMES, &block);

    // A reader that decoded the run before returning it would need every
    // frame's bytes at once; the lazy one needs a bounded few.
    let whole_run_bytes = LONG_FRAMES as usize * frame_bytes;
    assert!(
        long_peak < whole_run_bytes,
        "reading {LONG_FRAMES} frames peaked at {long_peak} bytes, \
         which is not less than the {whole_run_bytes} bytes of the run itself"
    );

    // The criterion itself: a run 1024 times longer costs the same memory.
    let slack = SLACK_FRAMES * frame_bytes;
    assert!(
        long_peak <= short_peak + slack,
        "reading {LONG_FRAMES} frames peaked at {long_peak} bytes against \
         {short_peak} for {SHORT_FRAMES} frames, more than the {slack} bytes \
         of headroom a bounded reader is allowed"
    );
}
