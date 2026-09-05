# ADR-0008: A counting global allocator, in one test binary

## Status
Accepted. Records the workspace's only use of `unsafe`, which
[CODING_STANDARDS.md](../../../CODING_STANDARDS.md) § *Correctness and failure*
("No `unsafe` without an ADR") requires an ADR for.

## Context
T-05.3's `RunReader` iterates a run lazily so that reading it does not put the
whole run in memory — load-bearing under [ADR-0006](0006-web-visualizer.md),
where a run sits in a browser tab rather than on local disk. The ticket asks
for that property to be *verified*: "memory usage doesn't scale with total run
length (verified by a test with a deliberately long run and a memory-bound
assertion, or at minimum a design review note if a hard memory test proves
impractical)."

A memory-bound assertion needs a measurement, and safe Rust has no way to ask
how much heap a piece of code is holding. The only portable instrument is a
custom `GlobalAlloc`, whose `impl` is `unsafe` — not because it does anything
dangerous here, but because implementing that trait is an unsafe contract
whatever the body does.

## Decision

**`termocline-format/tests/reader_memory.rs` installs a counting global
allocator, and it is the only `unsafe` in the workspace.** The allocator
forwards every method to `System` unchanged and adds two atomic counters; the
test reads the peak live heap while iterating a run of 8 frames and again at
8192, and asserts the two agree.

It is confined to that one test binary, which is why the test is a file of its
own: a global allocator is process-wide, and a peak-allocation measurement
means nothing if other tests are allocating beside it.

## Considered options

- **Bound the bytes pulled from the byte source instead**, with a counting
  `Read`. Safe, and it does prove the reader decodes lazily — but it cannot
  see retention. A reader that decoded one frame at a time and pushed every
  one into an internal `Vec` would pass it while growing without bound, which
  is the exact failure the criterion is about.
- **Read the process's resident set size.** Needs `libc` (so, `unsafe` anyway)
  or `/proc`, which macOS does not have; and RSS is too noisy for a bound this
  tight.
- **Skip the measurement and write the design note the ticket allows.** The
  ticket's own fallback, and the honest choice if the test were impractical.
  It is not: the test is thirty lines of forwarding and runs in three seconds.
  A note asserts nothing, and this is the property a browser tab runs out of
  memory over.
- **A crate such as `stats_alloc`.** The same `unsafe` behind someone else's
  name, plus a dependency, and it would still need this ADR.

## Consequences
- The workspace has exactly one `unsafe` block region, in test code, reviewed
  here. Any further use — in a crate or in another test — needs its own ADR;
  this one is not a general licence.
- `cargo clippy -- -D warnings` and the standards review both keep working as
  the tripwire: a second `unsafe` appearing without an ADR is a review defect,
  not a precedent.
- The measurement covers the reader, not the byte source it is handed. A
  browser that holds an entire run in a `Vec<u8>` before reading it pays for
  that buffer regardless; bounding *that* is Epic 08's problem, and streaming
  the fetch is what the `Read`-only reader leaves open.
