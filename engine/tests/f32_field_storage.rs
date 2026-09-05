//! The engine ships at `f64`, and T-10.4's probe is off (T-10.4).
//!
//! `engine/src/precision.rs` carries an instrument that stores every
//! prognostic state at `f32`, so that Epic 07's validation suite can be re-run
//! at the width T-10.4 proposes without a single one of its tolerances being
//! touched. `docs/performance-notes.md` reports what that measured and what
//! was decided from it.
//!
//! An instrument that could be left switched on by accident would be worse
//! than no instrument: every budget in Epic 07 is derived against `f64`'s unit
//! round-off, and a suite silently running at `2²⁹` times that would be
//! asserting something nobody wrote down. So the probe is a `--cfg` flag rather
//! than a cargo feature — CI runs `cargo test --workspace --all-features`,
//! which would switch a feature on and cannot reach a `--cfg` — and this file
//! is the guard that the default is the default. `cargo test` proves the
//! shipped engine is the `f64` one;
//! `RUSTFLAGS="--cfg f32_storage_probe" cargo test -p engine --no-fail-fast`
//! is the measurement.
//!
//! CODING_STANDARDS.md § *Physical quantities* is the rule this holds:
//! **`f64` throughout the solver.**

use engine::StorageWidth;
#[cfg(not(f32_storage_probe))]
use engine::FIELD_STORAGE;

/// The width the engine stores its fields at, asserted rather than assumed.
///
/// Compiled only into an unprobed build, which is the only build whose passing
/// means anything: a probe run is the measurement rather than the engine, so
/// this assertion would be a claim about the wrong thing there. What says a
/// probe run is a probe run is that the rest of the suite goes red — the
/// finding `docs/performance-notes.md` § *After T-10.4* records.
#[test]
#[cfg(not(f32_storage_probe))]
fn the_shipped_engine_stores_its_fields_at_f64() {
    assert_eq!(
        FIELD_STORAGE,
        StorageWidth::F64,
        "the engine is storing its fields at {FIELD_STORAGE:?}, but every error budget in Epic \
         07 is derived against f64's unit round-off. See docs/performance-notes.md, T-10.4"
    );
}

/// The two widths differ by the factor the whole ticket turns on: `f32`'s unit
/// round-off is `2²⁹` times `f64`'s.
///
/// Written out because it is the number every argument in
/// `docs/performance-notes.md` § *After T-10.4* is scaled by, and because it is the
/// one thing about the two widths that no measurement can move.
#[test]
fn the_narrow_width_rounds_off_two_to_the_twenty_ninth_more_coarsely() {
    let ratio = StorageWidth::F32.unit_round_off() / StorageWidth::F64.unit_round_off();
    assert_eq!(ratio, 2.0_f64.powi(29));
}

/// Rounding to `f64` is the identity, on values that are not `f32`s.
///
/// This is what makes an unprobed build the engine rather than a perturbation
/// of it, so it is asserted on a value chosen to have bits an `f32` could not
/// hold: `√2` is irrational, so its `f64` and `f32` roundings differ.
#[test]
fn rounding_to_the_wide_width_changes_nothing() {
    let irrational = std::f64::consts::SQRT_2;
    assert_eq!(StorageWidth::F64.round(irrational), irrational);
    assert_ne!(StorageWidth::F32.round(irrational), irrational);
}
