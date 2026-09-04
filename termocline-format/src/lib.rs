//! The on-disk contract between the engine and the visualizer.
//!
//! Per [ADR-0001] the two programs communicate only through files, and per
//! [ADR-0004] this crate is the single definition of that format: a JSON
//! header plus a sequence of binary frames. It depends on neither simulation
//! logic nor UI code, so both sides can share it without sharing anything
//! else.
//!
//! The format itself lands in Epic 05; this crate is currently a placeholder.
//!
//! [ADR-0001]: ../../docs/planning/adr/0001-engine-visualizer-split.md
//! [ADR-0004]: ../../docs/planning/adr/0004-data-interchange-format.md

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
}

#[cfg(test)]
mod gate_verification {
    /// Deliberately failing: proves the ruleset blocks a red PR. Removed with
    /// this throwaway branch.
    #[test]
    fn this_must_fail() {
        assert_eq!(2 + 2, 5, "intentional failure for gate verification");
    }
}
