//! The CFL bound on the timestep, and the check a run passes before it starts.
//!
//! An explicit scheme carries a stability limit: a timestep past it does not
//! merely lose accuracy, it amplifies the grid-scale mode every step until the
//! run is numerically meaningless. This module turns the cell spacing and the
//! fastest signal speed into the largest timestep the scheme fixed in
//! [ADR-0003] can take, so the engine can refuse an unsafe one up front rather
//! than write a file full of garbage.
//!
//! # Where the bound comes from
//!
//! A plane wave `exp(i(kx + ly))` sampled on the C-grid turns the centred
//! difference `(f[i+1] − f[i]) / dx` into a multiplication by
//! `2i·sin(k·dx/2)/dx`. The linear gravity-wave operator therefore has purely
//! imaginary eigenvalues `±i·c·κ`, and the grid-scale mode (`k·dx = l·dy = π`)
//! maximises `κ` at
//!
//! ```text
//! κ_max = 2·√(1/dx² + 1/dy²)
//! ```
//!
//! Classic RK4's stability region meets the imaginary axis at `|λ·dt| ≤ 2√2`
//! ([`RK4_IMAGINARY_AXIS_LIMIT`]), so the scheme is stable while
//!
//! ```text
//! dt ≤ 2√2 / (c · κ_max)
//! ```
//!
//! and [`max_stable_dt`] returns [`CFL_SAFETY_FACTOR`] times that.
//!
//! Both axes enter through `κ_max`: on an anisotropic grid the bound is
//! stricter than the smaller spacing alone would suggest, because the fastest
//! mode is the diagonal one.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use crate::Spacing;
use std::fmt;

/// Where RK4's stability region crosses the imaginary axis, as `|λ·dt|`.
///
/// The classic four-stage method's amplification factor is the degree-4
/// truncation of the exponential, `R(z) = 1 + z + z²/2 + z³/6 + z⁴/24`, and
/// `|R(iθ)| ≤ 1` exactly on `θ ∈ [0, 2√2]` — the standard result for RK4
/// (Hairer & Wanner, *Solving Ordinary Differential Equations I*, § II.2), and
/// the reason ADR-0003's choice of RK4 buys a longer timestep than a
/// forward-Euler scheme, whose stability region touches the imaginary axis
/// only at the origin.
pub const RK4_IMAGINARY_AXIS_LIMIT: f64 = 2.0 * std::f64::consts::SQRT_2;

/// Dimensionless margin held back from the raw stability bound.
///
/// The bound above is derived for the gravity-wave terms alone. A run also
/// carries Rayleigh damping and the wind forcing, which move the eigenvalues
/// off the pure imaginary axis and can nudge a marginally-stable step over the
/// boundary. `0.8` keeps the timestep a fifth clear of it — the customary
/// margin in ocean models, and small enough that a run is bounded by physics
/// rather than by the margin.
///
/// The margin does **not** absorb the Coriolis term. Rotation is a second
/// oscillation with a frequency of its own, `|f| = β·|y|`, and its stability
/// limit involves neither the wave speed nor the cell spacing, so no fixed
/// factor on this bound can cover it. It is a separate bound, enforced by the
/// engine's solver where both terms are visible at once; see
/// [ADR-0007](../../docs/planning/adr/0007-rotation-timestep-bound.md).
///
/// It is a project policy number chosen in T-01.3, not a measured physical
/// constant and not a value taken from the literature — there is nothing to
/// cite but this comment. A scenario that wants more margin asks for a smaller
/// `dt`, which is always allowed.
pub const CFL_SAFETY_FACTOR: f64 = 0.8;

/// Why a timestep or a wave speed could not be accepted.
///
/// These all describe invalid *scenario input* rather than a broken invariant,
/// so they are returned rather than panicked, and each names the offending
/// value and the bound it violated.
#[derive(Debug, Clone, PartialEq)]
pub enum CflError {
    /// The fastest wave speed was zero, negative, or not a finite number.
    WaveSpeedNotPositive {
        /// The value supplied, in metres per second.
        value_m_per_s: f64,
    },
    /// The timestep was zero, negative, or not a finite number.
    TimestepNotPositive {
        /// The value supplied, in seconds.
        value_s: f64,
    },
    /// The timestep was longer than the CFL-stable maximum for this grid and
    /// wave speed. The run is refused rather than quietly shortened, per
    /// CODING_STANDARDS.md § *No silent clamping*.
    TimestepExceedsCfl {
        /// The timestep asked for, in seconds.
        requested_s: f64,
        /// The largest timestep this grid and wave speed allow, in seconds.
        max_stable_s: f64,
    },
}

impl fmt::Display for CflError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WaveSpeedNotPositive { value_m_per_s } => write!(
                f,
                "wave speed is {value_m_per_s} m/s; the fastest wave speed must be finite and \
                 greater than 0"
            ),
            Self::TimestepNotPositive { value_s } => write!(
                f,
                "dt is {value_s} s; the timestep must be finite and greater than 0"
            ),
            Self::TimestepExceedsCfl {
                requested_s,
                max_stable_s,
            } => write!(
                f,
                "dt is {requested_s} s, past the CFL-stable maximum of {} s for this grid \
                 spacing and wave speed; the run would go unstable. Set dt to at most {} s, or \
                 coarsen the grid",
                suggestable_bound_s(*max_stable_s),
                suggestable_bound_s(*max_stable_s)
            ),
        }
    }
}

impl std::error::Error for CflError {}

/// Steps per second the suggested bound is rounded to — a millisecond.
///
/// The exact bound is a quotient of square roots, so it prints as something
/// like `40000.00000000001 s`. Telling a user to "set dt to at most
/// 40000.00000000001 s" is not actionable, so the message rounds *down* to
/// the millisecond: the suggestion stays inside the stability region, which
/// rounding to nearest would not guarantee.
const SUGGESTION_STEPS_PER_S: f64 = 1_000.0;

/// The CFL bound rounded down to the millisecond, for a message a user can act
/// on. The un-rounded bound stays in [`CflError::TimestepExceedsCfl`] for
/// callers that want the exact number.
fn suggestable_bound_s(max_stable_s: f64) -> f64 {
    (max_stable_s * SUGGESTION_STEPS_PER_S).floor() / SUGGESTION_STEPS_PER_S
}

/// The fastest signal speed the grid has to carry, in metres per second.
///
/// For the 1.5-layer model that is the Kelvin wave speed `c = √(g'·H)` (see
/// `CONTEXT.md`), but nothing here depends on where the number came from —
/// this crate stays physics-free. Validating the speed on the way in is what
/// lets [`max_stable_dt`] be total: there is no such thing as a wave speed for
/// which the bound is undefined.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WaveSpeed {
    m_per_s: f64,
}

impl WaveSpeed {
    /// A wave speed of `m_per_s` metres per second.
    ///
    /// # Errors
    /// [`CflError::WaveSpeedNotPositive`] if the speed is zero, negative or
    /// not finite. Zero is rejected rather than treated as "no limit": it
    /// would make the CFL bound infinite and silently disable the check.
    pub fn new(m_per_s: f64) -> Result<Self, CflError> {
        if !m_per_s.is_finite() || m_per_s <= 0.0 {
            return Err(CflError::WaveSpeedNotPositive {
                value_m_per_s: m_per_s,
            });
        }
        Ok(Self { m_per_s })
    }

    /// The speed, in metres per second.
    #[must_use]
    pub const fn m_per_s(self) -> f64 {
        self.m_per_s
    }
}

/// The longest timestep, in seconds, that keeps the scheme stable on this grid
/// at this wave speed — safety factor included.
///
/// Both arguments are validated on construction, so the bound is always a
/// finite, positive number of seconds and no `Result` is needed. See the
/// module comment for the derivation, and [`CFL_SAFETY_FACTOR`] for the
/// margin held back from the raw bound.
#[must_use]
pub fn max_stable_dt(grid_spacing: Spacing, wave_speed: WaveSpeed) -> f64 {
    let dx_m = grid_spacing.dx_m();
    let dy_m = grid_spacing.dy_m();
    // κ_max = 2·√(1/dx² + 1/dy²), the grid-scale mode of the centred-difference
    // gravity-wave operator.
    let fastest_wavenumber_per_m = 2.0 * (1.0 / (dx_m * dx_m) + 1.0 / (dy_m * dy_m)).sqrt();
    CFL_SAFETY_FACTOR * RK4_IMAGINARY_AXIS_LIMIT / (wave_speed.m_per_s() * fastest_wavenumber_per_m)
}

/// Accept `dt_s` for a run on this grid at this wave speed, or say why not.
///
/// The engine calls this before the first step of a run. It never adjusts the
/// timestep: an unsafe `dt` is the scenario's error to fix, and substituting a
/// safe value silently would hand back a run nobody asked for.
///
/// # Errors
/// [`CflError::TimestepNotPositive`] if `dt_s` is not a finite, positive
/// duration, or [`CflError::TimestepExceedsCfl`] if it is longer than
/// [`max_stable_dt`].
pub fn check_timestep(
    dt_s: f64,
    grid_spacing: Spacing,
    wave_speed: WaveSpeed,
) -> Result<(), CflError> {
    if !dt_s.is_finite() || dt_s <= 0.0 {
        return Err(CflError::TimestepNotPositive { value_s: dt_s });
    }
    let max_stable_s = max_stable_dt(grid_spacing, wave_speed);
    if dt_s > max_stable_s {
        return Err(CflError::TimestepExceedsCfl {
            requested_s: dt_s,
            max_stable_s,
        });
    }
    Ok(())
}
