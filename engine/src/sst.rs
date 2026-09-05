//! The mixed-layer SST anomaly `T'`, its equation, and the upwelling that
//! couples it to the thermocline.
//!
//! `CONTEXT.md` introduces the *SST anomaly* as belonging to the Epic 12
//! coupling extension and not to the linear ocean core, and this module is
//! that boundary made structural: nothing here is reachable from a run that
//! did not ask for it, and no term of the shallow-water equations reads `T'`.
//! The core stays exactly the model Epic 07 validated.
//!
//! # The equation
//!
//! The linearized mixed-layer temperature equation of the intermediate coupled
//! models (Zebiak & Cane, *Mon. Wea. Rev.* 115, 1987, § 2b; Battisti,
//! *J. Atmos. Sci.* 45, 1988), in the anomaly variables this project already
//! names:
//!
//! ```text
//! ∂T'/∂t = −u'·∂T̄/∂x + (w⁺/H_m)·(γ·h − T') − ε_T·T'
//! ```
//!
//! - `−u'·∂T̄/∂x` — the anomalous zonal current advecting the *mean* SST
//!   gradient. This is linear: `∂T̄/∂x` is a prescribed constant of the
//!   scenario, not a field the run evolves, so the term carries no product of
//!   two anomalies and the scope guard of CODING_STANDARDS.md § *Scope guards*
//!   is not touched. The `u'·∂T'/∂x` that would break it is deliberately
//!   absent.
//! - `(w⁺/H_m)·(γ·h − T')` — entrainment. Water is drawn into the mixed layer
//!   from just beneath it at the upwelling rate `w`, and the water down there
//!   is warmer when the thermocline is deeper, by `γ = ∂T_sub/∂h`. **This is
//!   the coupling to `h`** the ticket is about, and it is the ocean half of
//!   the Bjerknes feedback (`CONTEXT.md`): a flatter thermocline in the east
//!   means a shallower `h`, colder entrained water, and a cooler cold tongue.
//!   Only *upwelling* entrains, so the rate is `w⁺ = max(w, 0)`: when the
//!   surface layer converges, mixed-layer water leaves downward and brings
//!   nothing back. The clamp is on `w`, which is a function of the prescribed
//!   wind alone, so the term stays linear in the prognostic variables.
//! - `−ε_T·T'` — the surface heat flux relaxing an anomaly back to the
//!   climatology, the SST equation's counterpart of the Rayleigh damping the
//!   core already carries.
//!
//! # Why the advection is zonal only
//!
//! There is no `−v'·∂T̄/∂y` term, and that is a decision rather than an
//! omission. The mean SST of the equatorial Pacific is a *maximum* near the
//! equator, so `∂T̄/∂y` changes sign across it and vanishes on it — the one
//! row where this model does most of its work. A single prescribed constant,
//! which is what the zonal gradient legitimately is (the warm pool falls away
//! to the cold tongue almost uniformly along the equator), would therefore be
//! the wrong shape for the meridional one: it would advect heat across the
//! equator in a direction the real ocean does not. Representing it honestly
//! needs a `T̄(y)` profile rather than a number, which is a bigger change than
//! this ticket's equation and is not what makes the Bjerknes loop close.
//!
//! # The upwelling is implied by the wind, not prescribed
//!
//! `w` is not a parameter. It is diagnosed at each evaluation from the wind
//! stress the run is being forced by, through the steady Rayleigh-drag balance
//! of a wind-driven surface layer ([`SurfaceLayer`]):
//!
//! ```text
//! r_s·u_ml − f·v_ml = τx/(ρ₀·H_m)
//! r_s·v_ml + f·u_ml = τy/(ρ₀·H_m)
//! ```
//!
//! whose solution is [`mixed_layer_velocity_m_per_s`], and then
//! `w = H_m·(∂u_ml/∂x + ∂v_ml/∂y)` — what the surface layer's divergence must
//! be fed from below. The drag `r_s` is what makes this finite at the equator,
//! where the classical Ekman transport `−τx/(ρ₀·f)` is singular and where all
//! of the interesting upwelling is; far from the equator, `|f| ≫ r_s`, the
//! solution returns to Ekman's to within `(r_s/f)²`.
//!
//! Under the alizés (`τx < 0`) this gives poleward flow on both flanks of the
//! equator, a surface divergence there, and upward `w` — the equatorial
//! upwelling the cold tongue is made of.
//!
//! # Where it sits, and what it does not do
//!
//! `T'` is a cell-centered field, beside `h`, because the entrainment term
//! multiplies the two together and a term should not have to interpolate its
//! own state. The zonal current arrives from the east/west faces through the
//! C-grid average of `termocline-numerics`, as everywhere else in the engine.
//!
//! [`SstTerm`] *adds* to a tendency, in the same contract
//! [`CoriolisTerm`](crate::CoriolisTerm) keeps: the shallow-water evaluator
//! writes every point of the tendency first — including zeroing the SST rate,
//! which none of its terms contributes to — and this adds on top. The
//! no-normal-flow condition does not apply: `T'` is a scalar at a cell center,
//! not a velocity through a wall.
//!
//! Nothing here narrows the timestep. The rates this equation carries are
//! `w/H_m + ε_T`, of order `10⁻⁶ s⁻¹` at the parameters above — a decay
//! timescale of weeks against a gravity-wave CFL bound of hours, so the
//! stability of a step is still decided entirely by the two bounds
//! [`Solver::new`](crate::Solver::new) already checks.

use std::fmt;

use termocline_grid::{Field2D, Grid, Staggering, H_STAGGERING, U_STAGGERING, V_STAGGERING};
use termocline_numerics::{CGridOperators, Spacing};

use crate::boundary::NoNormalFlow;
use crate::coriolis::{row_of, BetaPlane};
use crate::forcing::WindStressField;
use crate::params::PhysicalParams;
use crate::state::OceanState;

/// Rayleigh drag `r_s` of the wind-driven surface layer, in s⁻¹, when a
/// scenario does not state one.
///
/// The inverse of two days, the value Zebiak & Cane (*Mon. Wea. Rev.* 115,
/// 1987, § 2b) use for the frictional surface layer of the equatorial Pacific.
/// It sets the width of the equatorial band over which the Ekman singularity
/// is resolved — `|βy| = r_s` at about 250 km — and, with it, the strength of
/// the equatorial upwelling.
pub const DEFAULT_SURFACE_DRAG_PER_S: f64 = 1.0 / (2.0 * 86_400.0);

/// Why a set of mixed-layer parameters was rejected.
///
/// These describe invalid *scenario input* — a `[sst]` section asking for a
/// mixed layer that cannot exist — so they are returned rather than panicked,
/// and each names the offending parameter and the value it carried
/// (CODING_STANDARDS.md § *Correctness and failure*).
#[derive(Debug, Clone, PartialEq)]
pub enum SstParamsError {
    /// A parameter that must be strictly positive and finite was not.
    NotPositive {
        /// Name of the parameter, matching its accessor.
        parameter: &'static str,
        /// The value supplied, in the unit the parameter's name states.
        value: f64,
    },
    /// A parameter that must be non-negative and finite was negative (or not a
    /// number).
    Negative {
        /// Name of the parameter, matching its accessor.
        parameter: &'static str,
        /// The value supplied, in the unit the parameter's name states.
        value: f64,
    },
    /// A parameter that may take either sign was not a finite number.
    NotFinite {
        /// Name of the parameter, matching its accessor.
        parameter: &'static str,
        /// The value supplied, in the unit the parameter's name states.
        value: f64,
    },
}

impl fmt::Display for SstParamsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotPositive { parameter, value } => write!(
                f,
                "{parameter} is {value}; it must be finite and greater than 0"
            ),
            Self::Negative { parameter, value } => write!(
                f,
                "{parameter} is {value}; it must be finite and at least 0"
            ),
            Self::NotFinite { parameter, value } => {
                write!(f, "{parameter} is {value}; it must be a finite number")
            }
        }
    }
}

impl std::error::Error for SstParamsError {}

/// The constants of one scenario's mixed layer: everything the SST anomaly
/// equation needs beyond the ocean parameters the core already carries.
///
/// Constructed once per run and read from every right-hand-side evaluation, so
/// it is `Copy` and validated at the boundary rather than at each use — the
/// same shape as [`PhysicalParams`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SstParams {
    /// Mixed-layer depth `H_m`, in metres. A *total* thickness, like `H` and
    /// unlike the anomalies the model solves for.
    mixed_layer_depth_m: f64,
    /// Rayleigh drag `r_s` of the surface layer, in s⁻¹.
    surface_drag_per_s: f64,
    /// Zonal gradient of the mean SST, `∂T̄/∂x`, in K/m. Negative in the
    /// equatorial Pacific, where the ocean cools from the warm pool eastward
    /// to the cold tongue.
    mean_zonal_sst_gradient_k_per_m: f64,
    /// Sensitivity `γ = ∂T_sub/∂h` of the entrained water's temperature to the
    /// thermocline depth anomaly, in K/m.
    subsurface_temperature_sensitivity_k_per_m: f64,
    /// Thermal damping `ε_T` of an SST anomaly, in s⁻¹.
    thermal_damping_per_s: f64,
}

impl SstParams {
    /// The mixed-layer parameter set, in SI units.
    ///
    /// `∂T̄/∂x` may take either sign — it is a gradient, and which way the
    /// ocean warms is the scenario's to say — so it is only checked finite.
    /// `γ` and `ε_T` may be zero, which switches off entrainment's sensitivity
    /// to the thermocline and the surface heat flux respectively, but never
    /// negative: a negative `γ` would make a deeper thermocline colder, and a
    /// negative `ε_T` would amplify an anomaly instead of relaxing it. `H_m`
    /// and `r_s` must be strictly positive, since both are divided by.
    ///
    /// # Errors
    /// An [`SstParamsError`] naming the first parameter that failed and the
    /// value it carried.
    pub fn new(
        mixed_layer_depth_m: f64,
        surface_drag_per_s: f64,
        mean_zonal_sst_gradient_k_per_m: f64,
        subsurface_temperature_sensitivity_k_per_m: f64,
        thermal_damping_per_s: f64,
    ) -> Result<Self, SstParamsError> {
        check_positive("mixed_layer_depth_m", mixed_layer_depth_m)?;
        check_positive("surface_drag_per_s", surface_drag_per_s)?;
        check_finite(
            "mean_zonal_sst_gradient_k_per_m",
            mean_zonal_sst_gradient_k_per_m,
        )?;
        check_non_negative(
            "subsurface_temperature_sensitivity_k_per_m",
            subsurface_temperature_sensitivity_k_per_m,
        )?;
        check_non_negative("thermal_damping_per_s", thermal_damping_per_s)?;
        Ok(Self {
            mixed_layer_depth_m,
            surface_drag_per_s,
            mean_zonal_sst_gradient_k_per_m,
            subsurface_temperature_sensitivity_k_per_m,
            thermal_damping_per_s,
        })
    }

    /// Mixed-layer depth `H_m`, in metres.
    #[must_use]
    pub const fn mixed_layer_depth_m(self) -> f64 {
        self.mixed_layer_depth_m
    }

    /// Rayleigh drag `r_s` of the surface layer, in s⁻¹.
    #[must_use]
    pub const fn surface_drag_per_s(self) -> f64 {
        self.surface_drag_per_s
    }

    /// Zonal gradient of the mean SST, `∂T̄/∂x`, in K/m.
    #[must_use]
    pub const fn mean_zonal_sst_gradient_k_per_m(self) -> f64 {
        self.mean_zonal_sst_gradient_k_per_m
    }

    /// Sensitivity `γ = ∂T_sub/∂h` of the entrained water to the thermocline
    /// depth anomaly, in K/m.
    #[must_use]
    pub const fn subsurface_temperature_sensitivity_k_per_m(self) -> f64 {
        self.subsurface_temperature_sensitivity_k_per_m
    }

    /// Thermal damping `ε_T` of an SST anomaly, in s⁻¹.
    #[must_use]
    pub const fn thermal_damping_per_s(self) -> f64 {
        self.thermal_damping_per_s
    }
}

fn check_positive(parameter: &'static str, value: f64) -> Result<(), SstParamsError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(SstParamsError::NotPositive { parameter, value });
    }
    Ok(())
}

fn check_non_negative(parameter: &'static str, value: f64) -> Result<(), SstParamsError> {
    if !value.is_finite() || value < 0.0 {
        return Err(SstParamsError::Negative { parameter, value });
    }
    Ok(())
}

fn check_finite(parameter: &'static str, value: f64) -> Result<(), SstParamsError> {
    if !value.is_finite() {
        return Err(SstParamsError::NotFinite { parameter, value });
    }
    Ok(())
}

/// The steady wind-driven velocity `(u_ml, v_ml)` of the mixed layer, in m/s,
/// under a stress of `(τx, τy)` pascals where the Coriolis parameter is
/// `coriolis_per_s`.
///
/// The solution of the surface layer's steady momentum balance,
///
/// ```text
/// r_s·u_ml − f·v_ml = τx/(ρ₀·H_m)
/// r_s·v_ml + f·u_ml = τy/(ρ₀·H_m)
/// ```
///
/// namely
///
/// ```text
/// u_ml = (r_s·τx + f·τy) / (ρ₀·H_m·(r_s² + f²))
/// v_ml = (r_s·τy − f·τx) / (ρ₀·H_m·(r_s² + f²))
/// ```
///
/// `layer_mass_kg_per_m2` is `ρ₀·H_m`, the mass of the mixed layer per unit
/// area — the same grouping the momentum equations divide a stress by, spelled
/// out once here rather than reassembled at each point.
///
/// Where `|f| ≫ r_s` the transport `H_m·v_ml` returns to the classical Ekman
/// (1905) transport `−τx/(ρ₀·f)`, to a relative accuracy of `(r_s/f)²`. Where
/// `f → 0` it does not blow up: the drag takes over, and the flow is downwind
/// at `τ/(ρ₀·H_m·r_s)`. That is the whole reason the drag is here — the
/// equator is where this model's upwelling is.
///
/// A free function rather than a method because it is pure: it is the closed
/// form a test compares the field against, and it carries no grid.
#[must_use]
pub fn mixed_layer_velocity_m_per_s(
    tau_x_pa: f64,
    tau_y_pa: f64,
    coriolis_per_s: f64,
    layer_mass_kg_per_m2: f64,
    surface_drag_per_s: f64,
) -> (f64, f64) {
    let denominator = layer_mass_kg_per_m2
        * (surface_drag_per_s * surface_drag_per_s + coriolis_per_s * coriolis_per_s);
    (
        (surface_drag_per_s * tau_x_pa + coriolis_per_s * tau_y_pa) / denominator,
        (surface_drag_per_s * tau_y_pa - coriolis_per_s * tau_x_pa) / denominator,
    )
}

/// The wind-driven surface layer of one basin, and the upwelling its
/// divergence implies.
///
/// Built once per run and re-diagnosed at every right-hand-side evaluation, so
/// it allocates nothing per call: every buffer below is allocated here and
/// reused across steps (CODING_STANDARDS.md § *Performance*).
///
/// The two stress components live on different faces of the C-grid — `τx` with
/// `u`, `τy` with `v` — but each mixed-layer velocity component needs both, so
/// each is interpolated onto the other's faces first. A C-grid interpolation is
/// undefined on a wall face, which has cells on one side only, so the
/// mixed-layer flow would be reading an unwritten number exactly at the coast.
/// It does not: the flow is brought onto the closed basin's boundary condition
/// before its divergence is taken
/// ([`NoNormalFlow::apply_to_surface_flow`]). A coast is a coast to the
/// surface layer as much as to the currents the solver integrates, so this is
/// the physics rather than a patch over the interpolation — and it is what
/// keeps a fictitious upwelling from appearing along the perimeter under any
/// wind with a meridional component.
#[derive(Debug, Clone)]
pub struct SurfaceLayer {
    /// Where the basin sits on the beta-plane, which is what `f` at each row
    /// comes from.
    plane: BetaPlane,
    /// The C-grid derivative and interpolation operators at this spacing.
    operators: CGridOperators,
    /// The mixed-layer constants this layer reads `H_m` and `r_s` from. Held
    /// whole rather than copied field by field: they arrive together, they are
    /// validated together, and a second copy of two of them is a place for the
    /// two to disagree.
    params: SstParams,
    /// Mass of the mixed layer per unit area, `ρ₀·H_m`, in kg/m². Derived from
    /// `params` and the ocean's `ρ₀` once, because it is what a stress is
    /// divided by at every point of every evaluation.
    layer_mass_kg_per_m2: f64,
    /// `τy` interpolated onto the east/west faces, where `u_ml` needs it.
    tau_y_on_u_faces_pa: Field2D<f64>,
    /// `τx` interpolated onto the north/south faces, where `v_ml` needs it.
    tau_x_on_v_faces_pa: Field2D<f64>,
    /// Zonal mixed-layer velocity `u_ml`, in m/s, at east/west faces.
    u_mixed_layer_m_per_s: Field2D<f64>,
    /// Meridional mixed-layer velocity `v_ml`, in m/s, at north/south faces.
    v_mixed_layer_m_per_s: Field2D<f64>,
    /// `∂u_ml/∂x` at cell centers, in s⁻¹. Scratch.
    zonal_divergence_per_s: Field2D<f64>,
    /// `∂v_ml/∂y` at cell centers, in s⁻¹. Scratch.
    meridional_divergence_per_s: Field2D<f64>,
    /// Upwelling `w` out of the mixed layer's base, in m/s, at cell centers.
    /// Positive is upward.
    upwelling_m_per_s: Field2D<f64>,
}

impl SurfaceLayer {
    /// The surface layer of `grid` at `spacing`, on `plane`, for an ocean with
    /// `params` and a mixed layer with `sst`.
    #[must_use]
    pub fn new(
        grid: Grid,
        spacing: Spacing,
        plane: BetaPlane,
        params: PhysicalParams,
        sst: SstParams,
    ) -> Self {
        Self {
            plane,
            operators: CGridOperators::new(grid, spacing),
            params: sst,
            layer_mass_kg_per_m2: params.reference_density_kg_per_m3() * sst.mixed_layer_depth_m(),
            tau_y_on_u_faces_pa: grid.allocate(U_STAGGERING, 0.0),
            tau_x_on_v_faces_pa: grid.allocate(V_STAGGERING, 0.0),
            u_mixed_layer_m_per_s: grid.allocate(U_STAGGERING, 0.0),
            v_mixed_layer_m_per_s: grid.allocate(V_STAGGERING, 0.0),
            zonal_divergence_per_s: grid.allocate(H_STAGGERING, 0.0),
            meridional_divergence_per_s: grid.allocate(H_STAGGERING, 0.0),
            upwelling_m_per_s: grid.allocate(H_STAGGERING, 0.0),
        }
    }

    /// Recompute the mixed-layer flow and the upwelling it implies under
    /// `wind_stress`.
    ///
    /// Every point of [`SurfaceLayer::upwelling_m_per_s`] is written, so the
    /// same buffer serves every stage of every step without carrying one
    /// stage's values into the next.
    ///
    /// # Panics
    /// If `wind_stress` covers a different basin from the one this layer was
    /// built for. A shape mismatch means the calling code is wrong, which is
    /// what panics are for (CODING_STANDARDS.md § *Correctness and failure*).
    pub fn diagnose(&mut self, wind_stress: &WindStressField) {
        assert!(
            wind_stress.grid() == self.operators.grid(),
            "wind stress covers {:?}, but this surface layer was built for {:?}",
            wind_stress.grid(),
            self.operators.grid()
        );

        self.operators
            .face_y_to_face_x(wind_stress.tau_y_pa(), &mut self.tau_y_on_u_faces_pa);
        self.operators
            .face_x_to_face_y(wind_stress.tau_x_pa(), &mut self.tau_x_on_v_faces_pa);

        let (mass, drag) = (self.layer_mass_kg_per_m2, self.params.surface_drag_per_s());
        let plane = self.plane;
        write_rows(
            &mut self.u_mixed_layer_m_per_s,
            wind_stress.tau_x_pa(),
            &self.tau_y_on_u_faces_pa,
            |j| plane.coriolis_at_row_per_s(U_STAGGERING, j),
            |tau_x, tau_y, coriolis| {
                mixed_layer_velocity_m_per_s(tau_x, tau_y, coriolis, mass, drag).0
            },
        );
        write_rows(
            &mut self.v_mixed_layer_m_per_s,
            &self.tau_x_on_v_faces_pa,
            wind_stress.tau_y_pa(),
            |j| plane.coriolis_at_row_per_s(V_STAGGERING, j),
            |tau_x, tau_y, coriolis| {
                mixed_layer_velocity_m_per_s(tau_x, tau_y, coriolis, mass, drag).1
            },
        );

        // The coast, before the divergence: no wind-driven flow through a wall,
        // and so no upwelling read off an interpolation that was never defined
        // there.
        NoNormalFlow::apply_to_surface_flow(
            &mut self.u_mixed_layer_m_per_s,
            &mut self.v_mixed_layer_m_per_s,
        );

        self.operators.ddx_face_to_center(
            &self.u_mixed_layer_m_per_s,
            &mut self.zonal_divergence_per_s,
        );
        self.operators.ddy_face_to_center(
            &self.v_mixed_layer_m_per_s,
            &mut self.meridional_divergence_per_s,
        );

        // What diverges out of the surface layer horizontally has to arrive
        // through its base, so the upward velocity there is the depth-
        // integrated divergence.
        let depth_m = self.params.mixed_layer_depth_m();
        let divergence = self
            .zonal_divergence_per_s
            .as_slice()
            .iter()
            .zip(self.meridional_divergence_per_s.as_slice());
        for (upwelling, (zonal, meridional)) in self
            .upwelling_m_per_s
            .as_mut_slice()
            .iter_mut()
            .zip(divergence)
        {
            *upwelling = depth_m * (zonal + meridional);
        }
    }

    /// Zonal wind-driven flow `u_ml` of the mixed layer, in m/s, at east/west
    /// faces, as of the last [`SurfaceLayer::diagnose`]. Zero on the western
    /// and eastern walls, which is the closed basin's boundary condition and
    /// not a gap.
    #[must_use]
    pub const fn zonal_flow_m_per_s(&self) -> &Field2D<f64> {
        &self.u_mixed_layer_m_per_s
    }

    /// Meridional wind-driven flow `v_ml` of the mixed layer, in m/s, at
    /// north/south faces. The twin of [`SurfaceLayer::zonal_flow_m_per_s`],
    /// held at zero on the southern and northern walls.
    #[must_use]
    pub const fn meridional_flow_m_per_s(&self) -> &Field2D<f64> {
        &self.v_mixed_layer_m_per_s
    }

    /// Upwelling `w` out of the mixed layer's base, in m/s, at cell centers,
    /// as of the last [`SurfaceLayer::diagnose`]. Positive is upward.
    #[must_use]
    pub const fn upwelling_m_per_s(&self) -> &Field2D<f64> {
        &self.upwelling_m_per_s
    }
}

/// Write `component(along, across, f(j))` into every point of `out`, with the
/// Coriolis parameter evaluated once per row because `f` depends on `y` alone.
///
/// All three fields are at the same staggering — the two stress components
/// have already been brought onto `out`'s faces — so this is a pointwise loop,
/// as in [`crate::coriolis::accumulate_rows`], and it reads its companion rows
/// through the same [`row_of`] that module owns.
fn write_rows(
    out: &mut Field2D<f64>,
    zonal_stress_pa: &Field2D<f64>,
    meridional_stress_pa: &Field2D<f64>,
    coriolis_per_s: impl Fn(usize) -> f64,
    component: impl Fn(f64, f64, f64) -> f64,
) {
    let points_per_row = out.nx();
    for (j, row) in out
        .as_mut_slice()
        .chunks_exact_mut(points_per_row)
        .enumerate()
    {
        let coriolis = coriolis_per_s(j);
        let zonal = row_of(zonal_stress_pa, j);
        let meridional = row_of(meridional_stress_pa, j);
        for ((value, tau_x), tau_y) in row.iter_mut().zip(zonal).zip(meridional) {
            *value = component(*tau_x, *tau_y, coriolis);
        }
    }
}

/// The SST anomaly equation's contribution to a tendency, over one basin.
///
/// Built once per run and applied at every right-hand-side evaluation, so it
/// allocates nothing per call: it owns a [`SurfaceLayer`] and the one extra
/// interpolation buffer the advection term needs
/// (CODING_STANDARDS.md § *Performance*).
#[derive(Debug, Clone)]
pub struct SstTerm {
    grid: Grid,
    params: SstParams,
    operators: CGridOperators,
    surface_layer: SurfaceLayer,
    /// The zonal current anomaly `u` interpolated onto the cell centers, in
    /// m/s, where `T'` lives. Reused every call; never read before written.
    zonal_current_on_centers_m_per_s: Field2D<f64>,
}

impl SstTerm {
    /// The SST term for `grid` at `spacing`, on `plane`, for an ocean with
    /// `params` and a mixed layer with `sst`.
    #[must_use]
    pub fn new(
        grid: Grid,
        spacing: Spacing,
        plane: BetaPlane,
        params: PhysicalParams,
        sst: SstParams,
    ) -> Self {
        Self {
            grid,
            params: sst,
            operators: CGridOperators::new(grid, spacing),
            surface_layer: SurfaceLayer::new(grid, spacing, plane, params, sst),
            zonal_current_on_centers_m_per_s: grid.allocate(H_STAGGERING, 0.0),
        }
    }

    /// The mixed-layer parameters this term was built with.
    #[must_use]
    pub const fn params(&self) -> SstParams {
        self.params
    }

    /// The surface layer this term diagnoses its upwelling from, as of the
    /// last [`SstTerm::add_to_tendency`].
    #[must_use]
    pub const fn surface_layer(&self) -> &SurfaceLayer {
        &self.surface_layer
    }

    /// Add `∂T'/∂t` to `tendency`, for `state` under `wind_stress`.
    ///
    /// Added rather than assigned, in the same contract
    /// [`CoriolisTerm`](crate::CoriolisTerm) keeps: the shallow-water
    /// evaluator has already written every point of the tendency, the SST rate
    /// included — it zeroes that one, since no term of the shallow-water
    /// equations contributes to it — and this adds the whole SST equation on
    /// top. Nothing here writes the `h`, `u` or `v` rates, which is what makes
    /// the extension additive.
    ///
    /// # Panics
    /// If `state`, `wind_stress` or `tendency` covers a different basin from
    /// the one this term was built for, or if either state is missing the SST
    /// anomaly a coupled run is integrating. Both mean the calling code is
    /// wrong, which is what panics are for (CODING_STANDARDS.md
    /// § *Correctness and failure*).
    pub fn add_to_tendency(
        &mut self,
        state: &OceanState,
        wind_stress: &WindStressField,
        tendency: &mut OceanState,
    ) {
        self.check_grid("state", state.grid());
        self.check_grid("wind stress", wind_stress.grid());
        self.check_grid("tendency", tendency.grid());

        self.surface_layer.diagnose(wind_stress);
        self.operators
            .face_to_center_x(state.u(), &mut self.zonal_current_on_centers_m_per_s);

        let SstParams {
            mixed_layer_depth_m,
            mean_zonal_sst_gradient_k_per_m,
            subsurface_temperature_sensitivity_k_per_m,
            thermal_damping_per_s,
            ..
        } = self.params;
        let anomaly_k = state
            .sst_anomaly_k()
            .expect("a coupled state carries an SST anomaly");
        let rate_k_per_s = tendency
            .sst_anomaly_k_mut()
            .expect("a coupled tendency carries an SST anomaly rate");

        let inputs = anomaly_k
            .as_slice()
            .iter()
            .zip(state.h().as_slice())
            .zip(self.zonal_current_on_centers_m_per_s.as_slice())
            .zip(self.surface_layer.upwelling_m_per_s().as_slice());
        for (rate, (((anomaly, thermocline_m), current_m_per_s), upwelling_m_per_s)) in
            rate_k_per_s.as_mut_slice().iter_mut().zip(inputs)
        {
            // Only upwelling entrains: when the surface layer converges,
            // mixed-layer water leaves downward and brings nothing back up.
            let entrainment_per_s = upwelling_m_per_s.max(0.0) / mixed_layer_depth_m;
            let subsurface_anomaly_k = subsurface_temperature_sensitivity_k_per_m * thermocline_m;
            *rate += -current_m_per_s * mean_zonal_sst_gradient_k_per_m
                + entrainment_per_s * (subsurface_anomaly_k - anomaly)
                - thermal_damping_per_s * anomaly;
        }
    }

    /// Panic unless `grid` is the basin this term was built for.
    fn check_grid(&self, role: &str, grid: Grid) {
        assert!(
            grid == self.grid,
            "{role} covers {grid:?}, but this SST term was built for {:?}",
            self.grid
        );
    }
}

/// Staggering of the SST anomaly `T'`: cell centers, beside `h`.
///
/// Named here rather than in `termocline-grid` because that crate carries the
/// staggering of the *linear core*, which `CONTEXT.md` is explicit `T'` is not
/// part of. The choice is the entrainment term's: it multiplies `T'` by `h`
/// point for point, and a variable should not have to interpolate the field it
/// is coupled to.
pub const SST_STAGGERING: Staggering = H_STAGGERING;
