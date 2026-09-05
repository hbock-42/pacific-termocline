//! The width the engine stores its fields at, and the probe that measures what
//! narrowing it would cost the answer (T-10.4).
//!
//! `docs/performance-notes.md` makes the case *for* narrowing. After T-10.5
//! the shallow-water right-hand side, the beta-plane rotation and the RK4
//! stage algebra are 97% of a timestep, and all three are
//! memory-bandwidth-bound: fourteen flat kernels paying for the traffic they
//! stream rather than for the arithmetic they do. T-10.3 measured the
//! corollary that adding cores to a saturated bus makes a step *slower*, and
//! left this one standing — halving the width of every field halves the
//! traffic, and `termocline-grid/examples/width_scaling.rs` measures what that
//! is worth on the engine's own basin.
//!
//! What it cannot say is what it costs. An `f32` carries about 7.2 decimal
//! digits against an `f64`'s 15.9, and Epic 07 does not assert that the engine
//! is roughly right: it asserts derived error budgets, the tightest of them
//! (T-07.4) to four significant figures with a margin of 1.0006x. So the
//! question is not an opinion about precision, it is a measurement — and this
//! module is the instrument that takes it.
//!
//! # The engine ships at `f64`
//!
//! [`FIELD_STORAGE`] is [`StorageWidth::F64`] in every build the project
//! produces, which is CODING_STANDARDS.md § *Physical quantities* — `f64`
//! throughout the solver — and `engine/tests/f32_field_storage.rs` guards it.
//! The `f32-storage-probe` feature is off by default, is never enabled by any
//! target in this repository, and exists so that the *validation suite itself*
//! can be re-run at the narrower width.
//!
//! # What the probe emulates, and why it is exact
//!
//! The ticket's proposal is `f32` **field storage** with `f64`
//! **accumulation**: the bulk grid data is narrow, the arithmetic over it is
//! not. That combination has an exact emulation in `f64`:
//!
//! - a value read from a narrow field widens to `f64` losslessly, because
//!   every `f32` is an `f64`;
//! - the arithmetic between reads is `f64` in both worlds, so it is the same
//!   arithmetic;
//! - a value written back to a narrow field rounds to the nearest `f32`.
//!
//! An `f64` field whose every value is exactly representable as an `f32`
//! therefore holds precisely the bits an `f32` field would, provided every
//! *store* rounds. [`narrow_stored_state`] is that store, and with the feature
//! off it is the identity — not "cheap", but literally nothing, since the body
//! is compiled out. The engine under test is the engine that ships.
//!
//! # What the probe does not narrow, and which way that points
//!
//! It rounds where a *state* is stored: the prognostic state itself, RK4's
//! stage state, and the four stage tendencies. A real `f32` field layout would
//! also narrow storage this probe leaves wide — the divergence scratch of
//! [`ShallowWaterRhs`](crate::ShallowWaterRhs), the interpolation buffers of
//! [`CoriolisTerm`](crate::CoriolisTerm), the sampled wind stress of
//! [`WindStressField`](crate::WindStressField), and the gradient written into
//! the tendency before being turned into an acceleration in place.
//!
//! Every one of those is a rounding the probe does not do, so what it measures
//! is a **lower bound** on the error a narrowed engine would make. That is the
//! useful direction, and it is the direction the conclusion in
//! `docs/performance-notes.md` is written in: a lower bound that already
//! crowds a validated budget settles the question, and no amount of care in
//! the implementation recovers the margin.
//!
//! # Running it
//!
//! ```sh
//! cargo test -p engine --features f32-storage-probe
//! ```
//!
//! That is Epic 07's suite — every tolerance exactly as it was derived, not
//! one of them touched — over an engine whose fields are stored at `f32`.

use crate::state::OceanState;
use termocline_grid::Field2D;

/// The width a field is stored at between one kernel and the next.
///
/// Two variants because T-10.4 is a comparison: [`StorageWidth::F64`] is the
/// engine as it is, and [`StorageWidth::F32`] is the engine the ticket
/// proposes. Running the same validation suite at both is the measurement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StorageWidth {
    /// The width the solver stores at: a double, ~15.9 decimal digits.
    F64,
    /// The width T-10.4 proposes for the bulk grid data: a single, ~7.2
    /// decimal digits.
    F32,
}

impl StorageWidth {
    /// The unit round-off of this width — half the gap between 1.0 and its
    /// successor, and the largest relative error a single store can introduce.
    ///
    /// Rust's `EPSILON` is the *gap*, so the round-off is half of it. The
    /// ratio between the two widths is `2²⁹`, and that factor is the whole of
    /// what this ticket is asking Epic 07 to absorb.
    #[must_use]
    pub const fn unit_round_off(self) -> f64 {
        match self {
            Self::F64 => 0.5 * f64::EPSILON,
            Self::F32 => 0.5 * f32::EPSILON as f64,
        }
    }

    /// `value` as it would come back out of a field stored at this width.
    #[must_use]
    pub fn round(self, value: f64) -> f64 {
        match self {
            Self::F64 => value,
            Self::F32 => f64::from(value as f32),
        }
    }

    /// Round every point of `field` to this width, in place.
    ///
    /// Only [`narrow_stored_state`] calls this, and only the probe build
    /// compiles that call, so without the feature there is nothing left to
    /// reach it — which is the point rather than an oversight.
    #[cfg_attr(not(feature = "f32-storage-probe"), allow(dead_code))]
    fn round_field(self, field: &mut Field2D<f64>) {
        if self == Self::F64 {
            return;
        }
        for value in field.as_mut_slice() {
            *value = self.round(*value);
        }
    }
}

/// The width this build stores its fields at.
///
/// [`StorageWidth::F64`] in every build the project ships; the
/// `f32-storage-probe` feature is what a measurement turns on, and nothing
/// else does.
pub const FIELD_STORAGE: StorageWidth = if cfg!(feature = "f32-storage-probe") {
    StorageWidth::F32
} else {
    StorageWidth::F64
};

/// Round a state that has just been written to [`FIELD_STORAGE`].
///
/// Called wherever the time loop *stores* a state — after the accumulation in
/// [`StateVector::add_scaled`](crate::StateVector::add_scaled), and after a
/// tendency is finished — so that with the probe on, a state held between two
/// kernels holds exactly the bits an `f32` field would.
///
/// With the probe off this is not a branch that is cheap to take: the body is
/// `cfg`-ed away entirely, so the shipped engine has no rounding in it at all.
#[inline]
pub fn narrow_stored_state(state: &mut OceanState) {
    #[cfg(feature = "f32-storage-probe")]
    {
        FIELD_STORAGE.round_field(state.h_mut());
        FIELD_STORAGE.round_field(state.u_mut());
        FIELD_STORAGE.round_field(state.v_mut());
        if let Some(sst_anomaly_k) = state.sst_anomaly_k_mut() {
            FIELD_STORAGE.round_field(sst_anomaly_k);
        }
    }
    #[cfg(not(feature = "f32-storage-probe"))]
    {
        let _ = state;
    }
}
