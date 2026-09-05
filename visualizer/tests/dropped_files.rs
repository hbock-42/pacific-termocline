//! A run reaching the shell by drag-and-drop.
//!
//! ADR-0006: a browser has no filesystem, so on the web a run arrives as its
//! two files dropped onto the window. They arrive one drop at a time and in
//! whatever order the user drags them, so the shell holds what it has been
//! given until the pair is complete.

use termocline_format::{FRAME_FILE_NAME, HEADER_FILE_NAME};
use visualizer::PendingRun;

#[test]
fn a_run_is_complete_only_once_both_of_its_files_have_arrived() {
    let mut pending = PendingRun::default();
    assert!(pending.offer(HEADER_FILE_NAME, b"{}".to_vec()));
    assert_eq!(pending.take_run(), None, "a header alone is not a run");
    assert_eq!(pending.still_needed(), vec![FRAME_FILE_NAME]);

    assert!(pending.offer(FRAME_FILE_NAME, b"frames".to_vec()));
    let run = pending.take_run().expect("both files have arrived");
    assert_eq!(run.header, b"{}");
    assert_eq!(run.frames, b"frames");
}

#[test]
fn the_two_files_may_be_dropped_in_either_order() {
    let mut pending = PendingRun::default();
    pending.offer(FRAME_FILE_NAME, b"frames".to_vec());
    assert_eq!(pending.still_needed(), vec![HEADER_FILE_NAME]);
    pending.offer(HEADER_FILE_NAME, b"{}".to_vec());
    assert!(pending.take_run().is_some());
}

#[test]
fn taking_a_run_clears_what_was_dropped() {
    // Otherwise the frames of one run would pair with the header of the next.
    let mut pending = PendingRun::default();
    pending.offer(HEADER_FILE_NAME, b"{}".to_vec());
    pending.offer(FRAME_FILE_NAME, b"frames".to_vec());
    assert!(pending.take_run().is_some());
    assert_eq!(pending.take_run(), None);
    assert_eq!(pending.still_needed().len(), 2);
}

#[test]
fn a_file_is_recognised_by_its_name_not_its_path() {
    // Native drops carry a full path; web drops carry a bare file name.
    let mut pending = PendingRun::default();
    assert!(pending.offer("/tmp/run-demo/header.json", b"{}".to_vec()));
    assert!(pending.offer("C:\\runs\\demo\\frames.bin", b"frames".to_vec()));
    assert!(pending.take_run().is_some());
}

#[test]
fn a_file_that_is_not_part_of_a_run_is_refused() {
    let mut pending = PendingRun::default();
    assert!(!pending.offer("notes.txt", b"hello".to_vec()));
    assert!(
        !pending.offer("Header.json", b"{}".to_vec()),
        "names are exact"
    );
    assert_eq!(pending.still_needed().len(), 2);
}

#[test]
fn a_second_drop_of_the_same_file_replaces_the_first() {
    // Dropping the right header after the wrong one should fix the run rather
    // than leave the first one in place.
    let mut pending = PendingRun::default();
    pending.offer(HEADER_FILE_NAME, b"first".to_vec());
    pending.offer(HEADER_FILE_NAME, b"second".to_vec());
    pending.offer(FRAME_FILE_NAME, Vec::new());
    assert_eq!(pending.take_run().expect("complete").header, b"second");
}
