//! The prognostic state of the 1.5-layer reduced-gravity ocean.
//!
//! [`OceanState`] holds the three prognostic variables of the shallow-water
//! core — the thermocline depth anomaly `h` and the current anomalies `u` and
//! `v` — each allocated at its own position on the Arakawa C-grid fixed in
//! [ADR-0003]. The staggering comes from `termocline-grid`'s named constants
//! rather than from index arithmetic here, so this module never has to know
//! that a face field carries one extra line of points.
//!
//! A fourth, *optional* field rides along: the mixed-layer SST anomaly `T'`
//! of the Epic 12 coupling extension. It is an `Option` rather than a field
//! that happens to be zero because `CONTEXT.md` is explicit that `T'` is "not
//! part of the linear ocean core": a run of the validated Epics 01-07 model
//! holds three fields, allocates three fields, and integrates three fields,
//! exactly as it did before the extension existed. Only a scenario that asks
//! for the coupling gets the fourth, through
//! [`OceanState::at_rest_with_sst_anomaly`]. See [`crate::sst`] for the
//! equation it obeys.
//!
//! The same type serves as a *tendency*: a right-hand-side evaluation (T-02.3)
//! produces an `OceanState` whose components are `∂h/∂t`, `∂u/∂t` and `∂v/∂t`
//! in SI-per-second, and [`StateVector`] is what lets [`crate::Rk4`] combine
//! those with a state.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use crate::integrator::StateVector;
use crate::params::PhysicalParams;
use crate::sst::SST_STAGGERING;
use termocline_grid::{Field2D, Grid, H_STAGGERING, U_STAGGERING, V_STAGGERING};

/// Value every prognostic variable takes in the undisturbed ocean.
///
/// All three are anomalies (`CONTEXT.md`), so the resting ocean is exactly
/// zero rather than some reference profile: at rest the thermocline sits at
/// its mean depth `H` and nothing is moving.
const AT_REST: f64 = 0.0;

/// The prognostic state of the upper layer over one basin.
///
/// The three fields are anomalies in SI units: `h` in metres, `u` and `v` in
/// m/s. `h` is a departure from the mean thermocline depth `H`, **not** a
/// total depth — the total is `H + h`, which
/// [`OceanState::total_thermocline_depth_m`] is there to spell out.
#[derive(Debug, Clone, PartialEq)]
pub struct OceanState {
    /// Shape of the basin the three fields cover.
    grid: Grid,
    /// Thermocline depth anomaly `h`, in metres, at cell centers.
    h: Field2D<f64>,
    /// Zonal current anomaly `u`, in m/s, at east/west faces.
    u: Field2D<f64>,
    /// Meridional current anomaly `v`, in m/s, at north/south faces.
    v: Field2D<f64>,
    /// Mixed-layer SST anomaly `T'`, in kelvin, at cell centers — `None` for
    /// the uncoupled linear core, which is what every run before Epic 12 is.
    sst_anomaly_k: Option<Field2D<f64>>,
}

impl OceanState {
    /// The ocean at rest over `grid`: every anomaly exactly zero.
    ///
    /// This is the initial condition of every unforced run and the state a
    /// damped, unforced basin relaxes back to.
    #[must_use]
    pub fn at_rest(grid: Grid) -> Self {
        Self {
            grid,
            h: grid.allocate(H_STAGGERING, AT_REST),
            u: grid.allocate(U_STAGGERING, AT_REST),
            v: grid.allocate(V_STAGGERING, AT_REST),
            sst_anomaly_k: None,
        }
    }

    /// The ocean at rest over `grid`, carrying the mixed-layer SST anomaly
    /// `T'` of the Epic 12 coupling as a fourth prognostic variable.
    ///
    /// The initial condition of a coupled run, and the only way a state comes
    /// to hold `T'`. `T'` is an anomaly like the other three, so at rest it is
    /// exactly zero: the mixed layer sits at its climatological temperature.
    #[must_use]
    pub fn at_rest_with_sst_anomaly(grid: Grid) -> Self {
        Self {
            sst_anomaly_k: Some(grid.allocate(SST_STAGGERING, AT_REST)),
            ..Self::at_rest(grid)
        }
    }

    /// Shape of the basin this state covers.
    #[must_use]
    pub const fn grid(&self) -> Grid {
        self.grid
    }

    /// Thermocline depth anomaly `h`, in metres, at cell centers.
    #[must_use]
    pub const fn h(&self) -> &Field2D<f64> {
        &self.h
    }

    /// Zonal current anomaly `u`, in m/s, at east/west faces.
    #[must_use]
    pub const fn u(&self) -> &Field2D<f64> {
        &self.u
    }

    /// Meridional current anomaly `v`, in m/s, at north/south faces.
    #[must_use]
    pub const fn v(&self) -> &Field2D<f64> {
        &self.v
    }

    /// Mixed-layer SST anomaly `T'`, in kelvin, at cell centers, or `None` if
    /// this state is the uncoupled linear core's.
    #[must_use]
    pub fn sst_anomaly_k(&self) -> Option<&Field2D<f64>> {
        self.sst_anomaly_k.as_ref()
    }

    /// Whether this state carries the Epic 12 SST anomaly.
    #[must_use]
    pub fn couples_sst(&self) -> bool {
        self.sst_anomaly_k.is_some()
    }

    /// Mutable access to the mixed-layer SST anomaly `T'`, in kelvin, or
    /// `None` if this state is the uncoupled linear core's.
    pub fn sst_anomaly_k_mut(&mut self) -> Option<&mut Field2D<f64>> {
        self.sst_anomaly_k.as_mut()
    }

    /// Mutable access to the thermocline depth anomaly `h`, in metres.
    pub fn h_mut(&mut self) -> &mut Field2D<f64> {
        &mut self.h
    }

    /// Mutable access to the zonal current anomaly `u`, in m/s.
    pub fn u_mut(&mut self) -> &mut Field2D<f64> {
        &mut self.u
    }

    /// Mutable access to the meridional current anomaly `v`, in m/s.
    pub fn v_mut(&mut self) -> &mut Field2D<f64> {
        &mut self.v
    }

    /// Total thermocline depth `H + h` at the cell center `(i, j)`, in metres,
    /// or `None` if that cell lies outside the basin.
    ///
    /// The one place the anomaly is turned into an absolute depth. Nothing in
    /// the solver needs this — the equations are written in `h` — but output
    /// and validation do, and doing the addition here keeps the distinction
    /// `CONTEXT.md` draws from being re-derived at each call site.
    #[must_use]
    pub fn total_thermocline_depth_m(
        &self,
        params: PhysicalParams,
        i: usize,
        j: usize,
    ) -> Option<f64> {
        self.h
            .get(i, j)
            .map(|anomaly_m| params.mean_thermocline_depth_m() + anomaly_m)
    }

    /// The prognostic variables in a fixed order, for the operations that
    /// treat the state as one vector rather than as separate fields.
    ///
    /// The three of the linear core first and the optional SST anomaly last,
    /// so that a coupled state's `h`, `u` and `v` are combined in exactly the
    /// order and exactly the associativity an uncoupled one's are — which is
    /// what makes the extension additive down to the last bit.
    fn components(&self) -> impl Iterator<Item = &Field2D<f64>> {
        [&self.h, &self.u, &self.v]
            .into_iter()
            .chain(self.sst_anomaly_k.as_ref())
    }

    /// The prognostic variables in the same order, mutably.
    fn components_mut(&mut self) -> impl Iterator<Item = &mut Field2D<f64>> {
        [&mut self.h, &mut self.u, &mut self.v]
            .into_iter()
            .chain(self.sst_anomaly_k.as_mut())
    }

    /// Round every prognostic field to `width`, in place — the *store* that
    /// T-10.4's probe narrows (`crate::precision`).
    ///
    /// Compiled only into a probe build, so a shipped engine does not merely
    /// skip the rounding, it does not carry the code that would do it.
    ///
    /// It lives here rather than in `precision` so that there is one list of
    /// what a state's components are: [`OceanState::components_mut`], whose
    /// order the extension's bit-exactness already depends on. A probe holding
    /// its own copy of that list would silently stop narrowing a fifth field
    /// the day one is added, and would say nothing about it.
    #[cfg(f32_storage_probe)]
    pub(crate) fn round_components_to(&mut self, width: crate::precision::StorageWidth) {
        for component in self.components_mut() {
            width.round_field(component);
        }
    }

    /// Panic unless `other` covers the same basin as `self` and agrees with it
    /// about whether the SST anomaly is being integrated.
    ///
    /// An `OceanState`'s shape is not in its type, so [`StateVector`] requires
    /// a panic rather than a truncation here: two states over different basins
    /// mean the calling code is wrong, and silently combining the overlap
    /// would be exactly the silent clamping CODING_STANDARDS.md forbids.
    fn check_same_grid(&self, other: &Self) {
        assert!(
            self.grid == other.grid,
            "ocean states cover different basins: this grid is {:?}, the other is {:?}",
            self.grid,
            other.grid
        );
        assert!(
            self.couples_sst() == other.couples_sst(),
            "ocean states disagree about the SST anomaly: this one {} it, the other {} it",
            carries(self.couples_sst()),
            carries(other.couples_sst())
        );
    }
}

impl StateVector for OceanState {
    fn assign(&mut self, source: &Self) {
        self.check_same_grid(source);
        for (target, values) in self.components_mut().zip(source.components()) {
            target.as_mut_slice().copy_from_slice(values.as_slice());
        }
    }

    fn add_scaled(&mut self, factor: f64, other: &Self) {
        self.check_same_grid(other);
        for (target, values) in self.components_mut().zip(other.components()) {
            for (value, term) in target.as_mut_slice().iter_mut().zip(values.as_slice()) {
                *value += factor * term;
            }
        }
        // The accumulation above is the `f64` arithmetic of the scheme; this
        // is the *store* that follows it, and it is where T-10.4's probe
        // narrows a field. Nothing happens here in a shipped build — the body
        // is compiled out — so the width question is asked without the engine
        // carrying an answer to it (`crate::precision`).
        crate::precision::narrow_stored_state(self);
    }
}

/// "carries" or "does not carry", so that the mismatch panic above reads as a
/// sentence rather than as two booleans.
fn carries(couples_sst: bool) -> &'static str {
    if couples_sst {
        "carries"
    } else {
        "does not carry"
    }
}
