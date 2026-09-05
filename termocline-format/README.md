# termocline-format

The on-disk contract between the engine and the visualizer: a JSON header plus
a sequence of binary frames, defined once and shared by both sides.

Depends on neither simulation logic nor UI code — that independence is the
point (see [ADR-0004](../docs/planning/adr/0004-data-interchange-format.md)).

- `RunHeader` — written once per run, as JSON: format version, grid and basin
  extent, physical parameters, scenario description, the variable list with
  units, and the output cadence. Self-describing on purpose, so a reader never
  guesses at the shape or meaning of the frames beside it.
- `Frame` — one saved timestep: model time plus `h`, `u`, `v`, `τx` and `τy`,
  each a flat row-major buffer sized by where that variable sits on the
  C-grid (`GridSpec::field_len`), and — for a run whose scenario asked for the
  Epic 12 coupling — the mixed-layer SST anomaly `T'` in kelvin at cell
  centers.

`T'` is an `Option`, not a field that happens to be zero. A run of the linear
core never integrated one, and a buffer of zeros would claim the ocean sat at
exactly its climatological temperature; absence says the run has no `T'` to
report, which is the true statement about it.

The writer lands in T-05.2 and the reader in T-05.3; this crate holds the
types and nothing else.

## Versions

The header carries `format_version`, and this crate carries the two ends of the
range a build reads:

| Version | Frames |
| --- | --- |
| 1 | `t`, `h`, `u`, `v`, `tau_x`, `tau_y` |
| 2 (`FORMAT_VERSION`) | the above, then an optional `sst` |

Writers only ever write `FORMAT_VERSION`; readers accept anything from
`OLDEST_READABLE_FORMAT_VERSION` up to it, decoding each version with its own
frame layout, and refuse anything else by name. A version 1 run therefore still
opens, with its `T'` absent — see
[ADR-0011](../docs/planning/adr/0011-reading-runs-from-older-format-versions.md).
