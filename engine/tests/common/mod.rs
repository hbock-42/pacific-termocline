//! What more than one of the engine's integration tests needs.
//!
//! Kept here rather than copied into each test binary: `tests/*.rs` are
//! separate crates, so a helper two of them need is otherwise two helpers that
//! drift apart.

use std::fs;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU32, Ordering};

/// A directory of one's own under the system temp directory, removed when the
/// test that made it ends.
///
/// Tests that write real runs must not share one: `cargo test` runs them in
/// parallel threads of one process, and two runs into one directory would each
/// see the other's files. The name carries the caller's label, the process id
/// and a per-process counter, which is unique across both.
pub struct ScratchDir {
    path: PathBuf,
}

impl ScratchDir {
    /// A fresh empty directory, labelled `ticket` and `name` so a leftover one
    /// says which test left it.
    ///
    /// # Panics
    /// If the system temp directory cannot be written, which is a broken
    /// environment rather than a failed test.
    pub fn new(ticket: &str, name: &str) -> Self {
        static COUNTER: AtomicU32 = AtomicU32::new(0);
        let unique = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "termocline-{ticket}-{name}-{}-{unique}",
            process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).expect("the system temp directory is writable");
        Self { path }
    }

    /// Where the directory is.
    pub fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for ScratchDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}
