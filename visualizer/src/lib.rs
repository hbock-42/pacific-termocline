//! Rendering for the engine's output. Reads files through
//! [`termocline_format`] and never links against the simulation code, so
//! either side can be reimplemented without touching the other (ADR-0001).
//!
//! Rendering lands in Epics 08–09; this crate is currently a placeholder.

/// Re-exported so the visualizer and engine agree on one format version.
pub use termocline_format::FORMAT_VERSION;

#[cfg(test)]
mod tests {
    #[test]
    fn workspace_links_the_format_crate() {
        assert_eq!(crate::FORMAT_VERSION, termocline_format::FORMAT_VERSION);
    }
}
