# ADR-0011: A reader keeps the old frame layout rather than refusing or migrating old runs

## Status
Accepted

## Context
[ADR-0004](0004-data-interchange-format.md) versions the interchange format
from day one so that a change to it "doesn't silently break old runs or
require a 'does this file work with this build' guessing game". T-05.4 is the
first change to actually spend that version: the mixed-layer SST anomaly `T'`
of the Epic 12 coupling joins the frame, and `format_version` goes from 1 to
2.

That makes concrete a question the version field only deferred: what should a
build that reads version 2 do when handed a version 1 run? Version 1 runs are
not *wrong* about anything. A run of the validated linear core of Epics 01–07
is a complete, correct run; the only thing version 2 added to it is a field
that run never had, because the scenario never asked for the coupling. And
runs are archives — the validation results of Epic 07 and the profiles of
Epic 10 are directories on disk that nothing re-runs cheaply.

The frame layouts are genuinely different, not merely differently populated.
`bincode` writes an `Option` as a one-byte tag followed by the payload, so a
version 2 frame has one byte a version 1 frame does not — and decoding a
version 1 file with the version 2 layout would read the *next* frame's first
byte as that tag and desynchronize everything after it. There is no way to
read old bytes with new code by accident, and no way to read them by luck.

## Options considered
1. **Refuse anything but the current version.** What the reader did before
   this change, when there was only one version. Simple, and the failure is
   loud. But it turns every format change into a flag day for every archived
   run, and it refuses files that are not defective.
2. **Migrate old runs on open** — read version 1, write version 2 beside it,
   or rewrite in place. The reader then has one layout to think about. But
   `T'` does not exist in a version 1 run, so the migration would have to
   invent one; and a reader that rewrites the user's archive is doing
   something a reader should not do, particularly the browser one of
   [ADR-0006](0006-web-visualizer.md), which has no filesystem to write to.
3. **Read a range of versions, each with its own frame layout.** The reader
   keeps the version 1 layout as a decoder, converts what it decodes into the
   current `Frame`, and reports the absent `T'` as absent.

## Decision
Option 3. `termocline-format` publishes `OLDEST_READABLE_FORMAT_VERSION`
alongside `FORMAT_VERSION`; the reader maps a header's version to a frame
layout and refuses anything outside the range by name. Writers only ever write
`FORMAT_VERSION`, so there is exactly one layout being produced and a bounded
set being consumed.

A version 1 run's `T'` reads back as `None`, never as zeros. This is the part
that is a decision rather than a detail: a buffer of zeros would round-trip
perfectly and would state, in a unit the header tells a reader to believe, that
the entire basin sat at exactly its climatological temperature — a physical
claim about a run that never made one. `None` says the run has no `T'` to
report, which is true, and it forces every consumer to say what it does about
that instead of plotting a fabricated field.

## Consequences
- Every archived run stays readable. A version 1 run opens, inspects, plots
  and scrubs exactly as it did, minus a variable it never had.
- The reader carries one decoder per readable version. That is a real cost and
  it grows: each future version adds a layout that must be kept, and tested,
  for as long as its runs are worth reading. The cost is bounded by dropping
  a version from the range deliberately — which is then a decision with a
  date, not an accident of a struct edit.
- `Frame::field` and `Frame::sst_anomaly_k` return `Option`, so absence is in
  the type and a consumer cannot read a zero it has no way to tell from a
  measurement.
- The range is not open-ended in the other direction: a version *newer* than
  this build's is still refused, because its frames may have any layout at all.
