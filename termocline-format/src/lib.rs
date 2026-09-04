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
//! The types here are data and nothing else: they carry no solver, and they do
//! not choose an encoding. The header is JSON and the frames are `bincode`
//! because the *writer* (T-05.2) and *reader* (T-05.3) say so; every type in
//! this crate is plain `serde`, and round-trips losslessly through either.
//!
//! [ADR-0001]: ../../docs/planning/adr/0001-engine-visualizer-split.md
//! [ADR-0004]: ../../docs/planning/adr/0004-data-interchange-format.md

mod error;
mod frame;
mod header;
mod variable;

pub use error::FormatError;
pub use frame::Frame;
pub use header::{BasinExtent, GridSpec, OutputTiming, PhysicalParams, RunHeader};
pub use variable::{Variable, VariableSpec};

/// Version of the on-disk format, written into every run's header so a reader
/// can tell whether it understands a file rather than guessing.
pub const FORMAT_VERSION: u32 = 1;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_version_is_stable() {
        // A change here is a breaking change to every archived run; it should
        // be a deliberate edit accompanied by a migration note, not a drift.
        assert_eq!(FORMAT_VERSION, 1);
    }

    #[test]
    fn every_variable_is_described_exactly_once() {
        // The header's variable list is what a reader indexes frames by, so a
        // duplicated or missing entry would silently mislabel a field.
        let symbols: Vec<&str> = Variable::ALL.iter().map(|v| v.symbol()).collect();
        assert_eq!(symbols, ["h", "u", "v", "tau_x", "tau_y"]);
    }
}
