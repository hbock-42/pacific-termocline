//! The right-hand side of the linear 1.5-layer shallow-water equations.
//!
//! The equations are the linear ones of `docs/planning/01-scientific-model.md`.
//! This module owns three of their terms — the pressure gradient `−g'·∇h` in
//! the momentum equations, the divergence `−H·(∂u/∂x + ∂v/∂y)` in the `h`
//! equation, and the Rayleigh damping `−r` on all three variables — plus the
//! surface stress that makes the wind-stress argument of
//! [`shallow_water_rhs`] mean something:
//!
//! ```text
//! ∂u/∂t = −g'·∂h/∂x + τx/(ρ₀·H) − r·u
//! ∂v/∂t = −g'·∂h/∂y + τy/(ρ₀·H) − r·v
//! ∂h/∂t = −H·(∂u/∂x + ∂v/∂y)    − r·h
//! ```
//!
//! The Coriolis terms (T-02.2) fold into the same evaluation later, and the
//! Epic 12 SST equation is added on top of the tendency this writes
//! ([`crate::sst`]); nothing else belongs here. In particular the v1 core is linear
//! (CODING_STANDARDS.md § Scope guards): the divergence carries the *mean*
//! depth `H`, never the total `H + h`, and there is no advection.
//!
//! Damping every prognostic variable at the same rate `r` is what makes the
//! decay a property of the whole basin rather than of one variable: it turns
//! the system into `ẋ = (L − r)·x`, where `L` is the conservative
//! pressure-gradient and continuity pair. Wherever `L` is skew in the discrete
//! energy `E = (g'/2)·Σh² + (H/2)·Σ(u² + v²)` it contributes nothing to it, so
//! `Ė = −2·r·E` and an unforced perturbation relaxes monotonically to the rest
//! state however the waves redistribute it — the acceptance criterion of
//! T-02.4.
//!
//! That skewness is a summation by parts whose boundary term is `h·u` at the
//! walls, so the identity is exact only while the wall faces carry no
//! velocity. This module cannot promise that on its own: it zeroes the
//! *pressure gradient* on a wall, which leaves a wall velocity to decay at
//! `−r` rather than forcing it to zero, and nothing here stops a stress
//! applied at the wall from starting one. What does is the boundary condition
//! of T-04.2: [`NoNormalFlow`](crate::NoNormalFlow) holds the wall velocities
//! at exactly zero at every RK4 stage, so the boundary term vanishes for all
//! time and the budget closes however the basin is forced.
//!
//! The spatial derivatives come from `termocline-numerics`, which is where the
//! C-grid neighbour arithmetic lives. Two consequences are worth stating,
//! because they are contracts rather than accidents:
//!
//! - A center→face difference is undefined on the basin's four boundary faces,
//!   which have a cell on one side only; the operators write zero there. So on
//!   a wall the acceleration this module produces is the wind stress alone,
//!   and the solver's boundary condition is what discards it.
//! - A face→center difference is defined at every cell, so `∂h/∂t` has no such
//!   gap.

use termocline_grid::{Field2D, Grid, H_STAGGERING};
use termocline_numerics::{CGridOperators, Spacing};

use crate::forcing::WindStressField;
use crate::params::PhysicalParams;
use crate::state::OceanState;

/// The right-hand side of the shallow-water equations for one basin, one cell
/// spacing and one parameter set, holding the scratch space it needs.
///
/// One evaluator is built per run and reused for every stage of every step, so
/// a whole simulation allocates its intermediate fields exactly once
/// (CODING_STANDARDS.md § Performance). [`shallow_water_rhs`] is the
/// allocating convenience wrapper for tests and one-off evaluation.
#[derive(Debug, Clone)]
pub struct ShallowWaterRhs {
    /// Physical parameters `(g', H, r, β, ρ₀)` of the ocean being integrated.
    params: PhysicalParams,
    /// The C-grid derivative operators at this basin's cell spacing.
    operators: CGridOperators,
    /// `∂u/∂x` at cell centers, in s⁻¹. Scratch, rewritten every evaluation.
    zonal_divergence_per_s: Field2D<f64>,
    /// `∂v/∂y` at cell centers, in s⁻¹. Scratch, rewritten every evaluation.
    meridional_divergence_per_s: Field2D<f64>,
}

impl ShallowWaterRhs {
    /// An evaluator for `grid` at `spacing`, for an ocean with `params`.
    #[must_use]
    pub fn new(grid: Grid, spacing: Spacing, params: PhysicalParams) -> Self {
        Self {
            params,
            operators: CGridOperators::new(grid, spacing),
            zonal_divergence_per_s: grid.allocate(H_STAGGERING, 0.0),
            meridional_divergence_per_s: grid.allocate(H_STAGGERING, 0.0),
        }
    }

    /// Write the time derivative of `state` under `wind_stress` into
    /// `tendency`.
    ///
    /// `tendency` is an [`OceanState`] read as rates: `∂h/∂t` in m/s, `∂u/∂t`
    /// and `∂v/∂t` in m/s². Every point of it is written, so the same buffer
    /// can be reused across RK4 stages without carrying a stage's values into
    /// the next.
    ///
    /// # Panics
    /// If `state`, `wind_stress` or `tendency` covers a different basin from
    /// the one this evaluator was built for. A shape mismatch means the
    /// calling code is wrong, which is what panics are for
    /// (CODING_STANDARDS.md § Correctness and failure).
    pub fn evaluate(
        &mut self,
        state: &OceanState,
        wind_stress: &WindStressField,
        tendency: &mut OceanState,
    ) {
        self.check_grid("state", state.grid());
        self.check_grid("wind stress", wind_stress.grid());
        self.check_grid("tendency", tendency.grid());
        assert!(
            state.couples_sst() == tendency.couples_sst(),
            "state and tendency disagree about the Epic 12 SST anomaly; a coupled state is \
             integrated with a coupled tendency or neither is"
        );

        // Momentum: the pressure gradient lands straight on the velocity
        // faces, so it is written into the tendency and then turned into an
        // acceleration in place.
        self.operators
            .ddx_center_to_face(state.h(), tendency.u_mut());
        self.operators
            .ddy_center_to_face(state.h(), tendency.v_mut());
        let minus_g_prime_m_per_s2 = -self.params.reduced_gravity_m_per_s2();
        // Mass of the upper layer per unit area, in kg/m²: the `ρ₀·H` a
        // surface stress in Pa is divided by to become an acceleration.
        let layer_mass_kg_per_m2 =
            self.params.reference_density_kg_per_m3() * self.params.mean_thermocline_depth_m();
        turn_gradient_into_acceleration(
            tendency.u_mut(),
            minus_g_prime_m_per_s2,
            wind_stress.tau_x_pa(),
            layer_mass_kg_per_m2,
        );
        turn_gradient_into_acceleration(
            tendency.v_mut(),
            minus_g_prime_m_per_s2,
            wind_stress.tau_y_pa(),
            layer_mass_kg_per_m2,
        );
        let damping_per_s = self.params.rayleigh_damping_per_s();
        subtract_damping(tendency.u_mut(), state.u(), damping_per_s);
        subtract_damping(tendency.v_mut(), state.v(), damping_per_s);

        // Continuity: both halves of the divergence land on the cell centers,
        // where `h` lives.
        self.operators
            .ddx_face_to_center(state.u(), &mut self.zonal_divergence_per_s);
        self.operators
            .ddy_face_to_center(state.v(), &mut self.meridional_divergence_per_s);
        let minus_mean_depth_m = -self.params.mean_thermocline_depth_m();
        let thickness_rate = tendency.h_mut().as_mut_slice().iter_mut();
        let divergence = self
            .zonal_divergence_per_s
            .as_slice()
            .iter()
            .zip(self.meridional_divergence_per_s.as_slice());
        for (rate, (zonal, meridional)) in thickness_rate.zip(divergence) {
            *rate = minus_mean_depth_m * (zonal + meridional);
        }
        subtract_damping(tendency.h_mut(), state.h(), damping_per_s);

        // No term of the shallow-water equations touches the Epic 12 SST
        // anomaly, and this evaluator promises to write *every* point of the
        // tendency — so on a coupled run it writes the one rate it does not
        // contribute to as the zero it contributes, leaving
        // [`SstTerm`](crate::SstTerm) to add the SST equation on top exactly
        // as the Coriolis term adds rotation.
        if let Some(rate_k_per_s) = tendency.sst_anomaly_k_mut() {
            rate_k_per_s.as_mut_slice().fill(0.0);
        }
    }

    /// Panic unless `grid` is the basin this evaluator was built for.
    fn check_grid(&self, role: &str, grid: Grid) {
        assert!(
            grid == self.operators.grid(),
            "{role} covers {grid:?}, but this right-hand side was built for {:?}",
            self.operators.grid()
        );
    }
}

/// Overwrite a gradient of `h`, in m/m, with the acceleration it produces:
/// `−g'·∂h/∂· + τ/(ρ₀·H)`, in m/s².
///
/// In-place because the gradient operators write straight into the tendency
/// buffer, which is where the acceleration has to end up; `field` therefore
/// arrives holding a gradient and leaves holding an acceleration.
///
/// `pub(crate)` so that [`crate::profiling`] can time this kernel on its own
/// rather than re-implementing it: a per-term profile of a re-implementation
/// would be a profile of the profiler.
pub(crate) fn turn_gradient_into_acceleration(
    field: &mut Field2D<f64>,
    minus_g_prime_m_per_s2: f64,
    stress_pa: &Field2D<f64>,
    layer_mass_kg_per_m2: f64,
) {
    for (value, stress) in field.as_mut_slice().iter_mut().zip(stress_pa.as_slice()) {
        *value = minus_g_prime_m_per_s2 * *value + stress / layer_mass_kg_per_m2;
    }
}

/// Subtract the Rayleigh damping `r·anomaly` from a tendency that already
/// holds the rest of its equation's terms.
///
/// Pointwise, and applied at every point of the field including the basin's
/// walls: damping is a property of the water at a point rather than of a
/// difference stencil, so unlike the pressure gradient it has no boundary gap.
/// `anomaly` and `tendency` are at the same staggering — this is `−r·u` on the
/// east/west faces, `−r·v` on the north/south faces, `−r·h` at cell centers —
/// so there is nothing to interpolate.
///
/// `pub(crate)` for the same reason as
/// [`turn_gradient_into_acceleration`]: [`crate::profiling`] times the kernel
/// itself, not a copy of it.
pub(crate) fn subtract_damping(
    tendency: &mut Field2D<f64>,
    anomaly: &Field2D<f64>,
    rayleigh_damping_per_s: f64,
) {
    for (rate, value) in tendency
        .as_mut_slice()
        .iter_mut()
        .zip(anomaly.as_slice().iter())
    {
        *rate -= rayleigh_damping_per_s * value;
    }
}

/// The time derivative of `state` under `wind_stress`, as a freshly allocated
/// tendency.
///
/// The named deliverable of T-02.3 and T-02.4, and the convenient form for a
/// test or a one-off evaluation. It allocates a state and an evaluator per call, so the
/// time loop uses [`ShallowWaterRhs::evaluate`] instead
/// (CODING_STANDARDS.md § Performance).
///
/// # Panics
/// If `wind_stress` covers a different basin from `state`.
#[must_use]
pub fn shallow_water_rhs(
    state: &OceanState,
    params: PhysicalParams,
    spacing: Spacing,
    wind_stress: &WindStressField,
) -> OceanState {
    let mut evaluator = ShallowWaterRhs::new(state.grid(), spacing, params);
    let mut tendency = OceanState::at_rest(state.grid());
    evaluator.evaluate(state, wind_stress, &mut tendency);
    tendency
}
