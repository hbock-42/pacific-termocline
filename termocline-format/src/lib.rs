//! The on-disk contract between the engine and the visualizer.
//!
//! Per [ADR-0001] the two programs communicate only through files, and per
//! [ADR-0004] this crate is the single definition of that format: a JSON
//! [`RunHeader`], written once, plus a sequence of `bincode`-encoded
//! [`Frame`]s, one per saved timestep. It depends on neither simulation logic
//! nor UI code, so both sides can share it without sharing anything else.
//!
//! # Reading the header losslessly
//!
//! `serde_json`'s default number parser is fast rather than exact: it can land
//! one ULP away from the `f64` that was written, which would silently perturb
//! a run's physical parameters. Any crate that reads a [`RunHeader`] back must
//! enable `serde_json`'s `float_roundtrip` feature; the round-trip test here
//! fails without it, and cargo's feature unification then carries it to every
//! reader in the workspace.
//!
//! # One place the format is defined
//!
//! The types here are data and nothing else: they carry no solver and no UI.
//! What they *do* carry, alongside the types, is the rest of what a run is on
//! disk — the two file names of [`HEADER_FILE_NAME`] and [`FRAME_FILE_NAME`],
//! and the `bincode` configuration of [`frame_encoding`]. None of that is
//! recoverable from the bytes, so a writer (T-05.2) and a reader (T-05.3) that
//! each chose their own would disagree silently. ADR-0004 asks for exactly one
//! place the format is defined; this crate is it, and that has to include the
//! choices the `serde` types alone do not express.
//!
//! # Reading a run without a filesystem
//!
//! [`RunReader`] is defined over byte sources rather than paths, per
//! [ADR-0006]: the visualizer runs in a browser, where a run arrives by file
//! selection, drag-and-drop or HTTP fetch and there is nothing to open. The
//! path-taking `RunReader::open` is a native convenience behind this crate's
//! default `fs` feature; turn default features off for `wasm32`.
//!
//! [ADR-0001]: ../../docs/planning/adr/0001-engine-visualizer-split.md
//! [ADR-0004]: ../../docs/planning/adr/0004-data-interchange-format.md
//! [ADR-0006]: ../../docs/planning/adr/0006-web-visualizer.md

mod error;
mod frame;
mod header;
mod reader;
mod variable;

pub use error::FormatError;

/// Re-exported so a caller of [`encode_frame`] can name what it refused
/// without depending on `bincode` itself: the encoding is this crate's choice
/// (ADR-0004), and so is the error it fails with.
pub use bincode::error::EncodeError;
pub use frame::{encode_frame, Frame};
pub use header::{BasinExtent, GridSpec, OutputTiming, PhysicalParams, RunHeader};
pub use reader::{decode_frame, RunReadError, RunReader};
pub use variable::{Variable, VariableSpec};

/// Version of the on-disk format, written into every run's header so a reader
/// can tell whether it understands a file rather than guessing.
///
/// Version 2 added the optional mixed-layer SST anomaly `T'` to the frame
/// (T-05.4). Version 1 runs are still read — see
/// [`OLDEST_READABLE_FORMAT_VERSION`].
pub const FORMAT_VERSION: u32 = 2;

/// The oldest format version this build reads.
///
/// Runs are archives: a version 1 run is a complete, valid run of the linear
/// core, and the only thing version 2 added to it is a field it never had. So
/// the reader keeps version 1's frame layout and reports the `T'` of such a
/// run as absent, rather than refusing an archive that is not wrong about
/// anything (ADR-0011). Writers only ever write [`FORMAT_VERSION`].
pub const OLDEST_READABLE_FORMAT_VERSION: u32 = 1;

/// Name of the JSON [`RunHeader`] inside a run directory.
///
/// A run on a filesystem is a directory of two files; on the web it is two
/// byte sources (ADR-0006) and these names are what a fetch path or a pair of
/// dropped files is matched against.
pub const HEADER_FILE_NAME: &str = "header.json";

/// Name of the binary [`Frame`] sequence inside a run directory.
pub const FRAME_FILE_NAME: &str = "frames.bin";

/// The `bincode` configuration every [`Frame`] is encoded and decoded with.
///
/// Half of the frame encoding lives in the `serde` derive on [`Frame`] and the
/// other half lives here: nothing in the bytes records which configuration
/// wrote them, so a reader that picks a different one decodes garbage without
/// noticing. It is a function rather than a `const` because
/// `bincode::config::standard` is not one.
#[must_use]
pub fn frame_encoding() -> bincode::config::Configuration {
    bincode::config::standard()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_version_is_stable() {
        // A change here is a breaking change to every archived run; it should
        // be a deliberate edit accompanied by a migration note, not a drift.
        assert_eq!(FORMAT_VERSION, 2);
        // And every version from the oldest readable one up to it must be a
        // version this build actually has a frame layout for.
        assert_eq!(OLDEST_READABLE_FORMAT_VERSION, 1);
    }

    #[test]
    fn every_variable_is_described_exactly_once() {
        // The header's variable list is what a reader indexes frames by, so a
        // duplicated or missing entry would silently mislabel a field.
        let symbols: Vec<&str> = Variable::ALL.iter().map(|v| v.symbol()).collect();
        assert_eq!(symbols, ["h", "u", "v", "tau_x", "tau_y", "sst"]);
        // The linear core is the first five of them, in the same order: a
        // coupled run's variable list extends an uncoupled one rather than
        // rearranging it.
        let core: Vec<&str> = Variable::LINEAR_CORE.iter().map(|v| v.symbol()).collect();
        assert_eq!(core, symbols[..core.len()]);
    }
}
