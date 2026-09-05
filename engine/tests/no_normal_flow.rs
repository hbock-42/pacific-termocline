//! Acceptance tests for T-04.2 — no-normal-flow at the closed basin's walls.
//!
//! The basin of `CONTEXT.md` is "closed on all four boundaries", and this is
//! what closed means on the Arakawa C-grid of [ADR-0003]: the zonal current
//! anomaly `u` is exactly zero on the western and eastern boundary faces (the
//! `u` columns `i = 0` and `i = nx`), and the meridional anomaly `v` is exactly
//! zero on the southern and northern ones (the `v` rows `j = 0` and `j = ny`).
//! Those four lines of points are not degrees of freedom of the closed basin;
//! they are where the coast is.
//!
//! # Why the tests below force the walls on purpose
//!
//! Every forced case here applies a surface stress *at* the wall faces, which
//! is the one thing that can start a wall flow. Nothing else can: the C-grid
//! pressure-gradient operators of T-01.1 leave the boundary faces at zero, the
//! Coriolis term interpolates a velocity that is itself zero there, and
//! Rayleigh damping only decays what is already moving. A test forced by a
//! stress that vanishes at the walls would therefore pass on an engine with no
//! boundary condition at all — as T-03.1 found, and worked around with a
//! sampling rule that this ticket replaces.
//!
//! So the fields below are built with
//! [`WindStressField::uniform`](engine::WindStressField::uniform), which
//! carries its stress at *every* face, walls included, and the scenario-driven
//! case uses a test wind that is deliberately non-zero at the coast. What the
//! boundary condition has to do is hold the walls at rest anyway.
//!
//! # Where the expected values come from
//!
//! Three of the four checks are exact statements about the discretisation
//! rather than measurements of it, so they are asserted at exactly zero or to
//! round-off:
//!
//! - **The walls.** `u = 0` and `v = 0` there are the boundary condition
//!   itself. RK4 forms every stage as `state + a·k`, so if both the state and
//!   every stage tendency carry an exact `0.0` on a wall face, so does the
//!   result — `0.0 + a·dt·0.0` is exactly `0.0` in IEEE-754. Nothing here is
//!   allowed a tolerance.
//! - **A stress at a wall does nothing.** With the wall velocities held, the
//!   `τ/(ρ₀·H)` a wall face receives is discarded at every stage, so a run
//!   forced at the walls must produce *bit-identical* interior fields to the
//!   same run whose wall stress was zeroed by hand.
//! - **Volume.** Summing the discrete continuity equation over the basin
//!   telescopes: `Σ_cells ∂h/∂t = −(H/Δx)·Σ_j (u[nx,j] − u[0,j])
//!   − (H/Δy)·Σ_i (v[i,ny] − v[i,0])`, the divergence theorem written on the
//!   C-grid. With no normal flow and no damping every term on the right is
//!   zero, so an undamped closed basin conserves `Σh` exactly, whatever the
//!   wind does to it. This is the check that distinguishes a condition applied
//!   at each RK4 *stage* from one applied only after the completed step: the
//!   latter leaves each intermediate stage with a wall velocity that leaks
//!   volume through the coast.
//!
//! [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md

use std::f64::consts::TAU;

use engine::{
    max_stable_dt, Basin, BetaPlane, Grid, OceanState, PhysicalParams, Solver, Spacing, WaveSpeed,
    WindStress, WindStressField, U_STAGGERING, V_STAGGERING,
};

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s². Standard value for the 1.5-layer model (Gill, *Atmosphere–Ocean
/// Dynamics*, ch. 11; Cane & Sarachik 1981).
const PACIFIC_REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H` of the equatorial Pacific, in metres — the
/// canonical 150 m upper layer of the same 1.5-layer configuration.
const PACIFIC_MEAN_DEPTH_M: f64 = 150.0;
/// The equatorial beta-plane gradient, in m⁻¹s⁻¹ — `CONTEXT.md`, *Beta-plane*.
const BETA_PER_M_PER_S: f64 = engine::EQUATORIAL_BETA_PER_M_PER_S;
/// Reference seawater density `ρ₀`, in kg/m³ — `CONTEXT.md` and Gill, appendix 3.
const REFERENCE_DENSITY_KG_PER_M3: f64 = engine::SEAWATER_REFERENCE_DENSITY_KG_PER_M3;

/// Rayleigh damping `r` of the forced runs, in s⁻¹: an `e`-folding time of
/// about 11.6 days. Far stronger than the equatorial Pacific's own damping,
/// for the reason `rayleigh_damping.rs` spells out — a decay has to be visible
/// inside a run of CFL-admissible steps.
const STRONG_DAMPING_PER_S: f64 = 1.0e-6;
/// Rayleigh damping of the volume-budget run, in s⁻¹. The telescoping identity
/// this file's header derives is a statement about the continuity equation
/// alone; `−r·h` is a sink that would hide a leak behind a decay, so the
/// budget is checked on the undamped basin where `Σh` is exactly conserved.
const UNDAMPED_PER_S: f64 = 0.0;

/// Zonal extent of the test basin, in metres — the order of the equatorial
/// Pacific's width (`CONTEXT.md`, *Basin*). Deliberately different from
/// [`BASIN_LY_M`] so an x/y swap cannot pass.
const BASIN_LX_M: f64 = 1.0e7;
/// Meridional extent of the test basin, in metres: an equatorial channel
/// reaching ±500 km, about 1.4 equatorial deformation radii either side of the
/// equator. Narrow for the reason `time_stepping.rs` gives at length — it
/// leaves the gravity-wave CFL bound the binding one rather than the rotation
/// bound of ADR-0007.
const BASIN_LY_M: f64 = 1.0e6;
/// Cells along x of the test basin.
const BASIN_NX: usize = 16;
/// Cells along y. Different from [`BASIN_NX`] so an x/y swap cannot pass.
const BASIN_NY: usize = 8;

/// Zonal wind stress of the forced cases, in Pa. Easterly trade-wind stress is
/// `τx < 0` (`CONTEXT.md`, *Wind stress*).
const TRADE_WIND_STRESS_X_PA: f64 = -0.05;
/// Meridional wind stress of the forced cases, in Pa. Non-zero so the
/// southern and northern walls are forced too, and different in magnitude from
/// [`TRADE_WIND_STRESS_X_PA`] so an x/y swap cannot pass.
const TRADE_WIND_STRESS_Y_PA: f64 = 0.02;

/// Steps every forced run takes. At the CFL-safe timestep of this basin
/// (≈3×10⁴ s) that is of the order of a year of simulated time — "many steps
/// under active wind forcing", not the initial condition.
const FORCED_RUN_STEPS: usize = 500;

/// A velocity a wall face must never carry, in m/s. The boundary condition is
/// exact, so this is `0.0` and not a tolerance.
const AT_REST_M_PER_S: f64 = 0.0;

/// The equatorial-Pacific parameter set at a given Rayleigh damping.
fn pacific_params(rayleigh_damping_per_s: f64) -> PhysicalParams {
    PhysicalParams::new(
        PACIFIC_REDUCED_GRAVITY_M_PER_S2,
        PACIFIC_MEAN_DEPTH_M,
        rayleigh_damping_per_s,
        BETA_PER_M_PER_S,
        REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the published equatorial-Pacific parameters are physical")
}

/// The test basin: [`BASIN_NX`] by [`BASIN_NY`] cells spanning [`BASIN_LX_M`]
/// by [`BASIN_LY_M`], with the equator through its middle.
fn test_basin() -> Basin {
    let grid = Grid::new(BASIN_NX, BASIN_NY).expect("extents are non-zero");
    let spacing = Spacing::new(BASIN_LX_M / BASIN_NX as f64, BASIN_LY_M / BASIN_NY as f64)
        .expect("a basin spanned by whole cells has positive spacing");
    Basin::new(grid, spacing, 0.0, -0.5 * BASIN_LY_M).expect("both edges are finite positions")
}

/// A solver for `basin` at the CFL-safe maximum timestep, and that timestep.
fn solver_for(basin: Basin, params: PhysicalParams) -> (Solver, f64) {
    let wave_speed =
        WaveSpeed::new(params.kelvin_wave_speed_m_per_s()).expect("a positive wave speed");
    let dt_s = max_stable_dt(basin.spacing(), wave_speed);
    let plane = BetaPlane::centered_on_equator(params, basin.spacing(), basin.grid());
    let solver = Solver::new(basin.grid(), basin.spacing(), params, plane, dt_s)
        .unwrap_or_else(|error| panic!("the test's own timestep must be admissible: {error}"));
    (solver, dt_s)
}

/// Panic unless every one of the basin's four wall lines is exactly at rest.
///
/// `when` names the moment, so a failure says which step opened the coast.
fn assert_walls_at_rest(state: &OceanState, when: &str) {
    let (nx, ny) = (state.grid().nx(), state.grid().ny());
    for j in 0..state.u().ny() {
        for (wall, i) in [("western", 0), ("eastern", nx)] {
            let u_m_per_s = *state.u().get(i, j).expect("a wall face of the u field");
            assert_eq!(
                u_m_per_s, AT_REST_M_PER_S,
                "{when}: the {wall} wall passes water at row {j}: u = {u_m_per_s} m/s"
            );
        }
    }
    for i in 0..state.v().nx() {
        for (wall, j) in [("southern", 0), ("northern", ny)] {
            let v_m_per_s = *state.v().get(i, j).expect("a wall face of the v field");
            assert_eq!(
                v_m_per_s, AT_REST_M_PER_S,
                "{when}: the {wall} wall passes water at column {i}: v = {v_m_per_s} m/s"
            );
        }
    }
}

/// The basin's total volume anomaly, in m³: `Σ h·Δx·Δy` over the cell centers.
///
/// The quantity the discrete continuity equation conserves in a closed,
/// undamped basin (this file's header).
fn volume_anomaly_m3(state: &OceanState, spacing: Spacing) -> f64 {
    let cell_area_m2 = spacing.dx_m() * spacing.dy_m();
    state.h().as_slice().iter().sum::<f64>() * cell_area_m2
}

/// The sum of `|h|·Δx·Δy` over the basin, in m³ — the scale the volume
/// anomaly's round-off is measured against, since cancellation is what the
/// budget is testing.
fn volume_scale_m3(state: &OceanState, spacing: Spacing) -> f64 {
    let cell_area_m2 = spacing.dx_m() * spacing.dy_m();
    state
        .h()
        .as_slice()
        .iter()
        .map(|h_m| h_m.abs())
        .sum::<f64>()
        * cell_area_m2
}

/// The largest `|h|` anywhere in the basin, in metres.
fn peak_anomaly_m(state: &OceanState) -> f64 {
    state
        .h()
        .as_slice()
        .iter()
        .fold(0.0_f64, |peak, h_m| peak.max(h_m.abs()))
}

/// A test double: a spatially uniform stress in *both* components that
/// strengthens and weakens with time.
///
/// Not a physical scenario — the alizés are easterly and have no meridional
/// component — but a `WindStress` whose value at the coast is deliberately
/// non-zero in both components, which is what makes the scenario-driven path
/// ([`Solver::step_forced_by`], which re-samples the trait at every RK4 stage)
/// test all four walls rather than only the zonal pair.
struct GustsAtEveryFace {
    /// Period over which the stress swings, in seconds.
    period_s: f64,
}

impl WindStress for GustsAtEveryFace {
    fn stress(&self, _x_m: f64, _y_m: f64, t_s: f64) -> (f64, f64) {
        let phase = TAU * t_s / self.period_s;
        (
            TRADE_WIND_STRESS_X_PA * phase.cos(),
            TRADE_WIND_STRESS_Y_PA * phase.sin(),
        )
    }
}

/// A test double: the same uniform stress in both components, but exactly zero
/// on the basin's coastline.
///
/// The control for [`a_stress_applied_at_a_wall_changes_nothing_in_the_basin`].
/// A wall face sits *on* an edge of the basin — the `τx` faces at
/// `x = west` and `x = west + nx·Δx`, the `τy` faces at `y = south` and
/// `y = south + ny·Δy` — and those positions are compared exactly because
/// [`Basin::x_of_column_m`] computes both sides from the same expression.
struct CalmAtTheCoast {
    /// The basin whose coastline is being avoided.
    basin: Basin,
}

impl WindStress for CalmAtTheCoast {
    fn stress(&self, x_m: f64, y_m: f64, _t_s: f64) -> (f64, f64) {
        let west_m = self.basin.x_of_column_m(U_STAGGERING, 0);
        let east_m = self
            .basin
            .x_of_column_m(U_STAGGERING, self.basin.grid().nx());
        let south_m = self.basin.y_of_row_m(V_STAGGERING, 0);
        let north_m = self.basin.y_of_row_m(V_STAGGERING, self.basin.grid().ny());
        let tau_x_pa = if x_m == west_m || x_m == east_m {
            0.0
        } else {
            TRADE_WIND_STRESS_X_PA
        };
        let tau_y_pa = if y_m == south_m || y_m == north_m {
            0.0
        } else {
            TRADE_WIND_STRESS_Y_PA
        };
        (tau_x_pa, tau_y_pa)
    }
}

/// A test double: the same uniform stress everywhere, coastline included.
struct GaleAtEveryFace;

impl WindStress for GaleAtEveryFace {
    fn stress(&self, _x_m: f64, _y_m: f64, _t_s: f64) -> (f64, f64) {
        (TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA)
    }
}

/// A test double: a stress that doubles from the basin's western edge to its
/// eastern one, and from its southern edge to its northern one.
///
/// The asymmetry is the point, and it is what
/// [`an_undamped_closed_basin_conserves_its_volume_anomaly`] needs. Under a
/// *uniform* stress the western and eastern walls accelerate identically, so a
/// basin whose coasts pass water lets exactly as much in at one wall as it
/// lets out at the other and conserves its volume anomaly by accident; the
/// budget would then pass with no boundary condition at all. Under a tilted
/// stress the two fluxes differ, and the budget can fail.
struct GaleTiltedAcrossTheBasin {
    /// The basin the tilt is measured across.
    basin: Basin,
}

impl WindStress for GaleTiltedAcrossTheBasin {
    fn stress(&self, x_m: f64, y_m: f64, _t_s: f64) -> (f64, f64) {
        let eastward = (x_m - self.basin.western_edge_x_m()) / self.basin.zonal_extent_m();
        let northward = (y_m - self.basin.southern_edge_y_m()) / self.basin.meridional_extent_m();
        (
            TRADE_WIND_STRESS_X_PA * (1.0 + eastward),
            TRADE_WIND_STRESS_Y_PA * (1.0 + northward),
        )
    }
}

/// Run `basin` from rest for [`FORCED_RUN_STEPS`] steps under a fixed stress
/// field, and return the final state.
fn run_under(basin: Basin, params: PhysicalParams, wind: &WindStressField) -> OceanState {
    let (mut solver, dt_s) = solver_for(basin, params);
    let mut state = OceanState::at_rest(basin.grid());
    for step in 0..FORCED_RUN_STEPS {
        solver.step(&mut state, step as f64 * dt_s, |_t| wind);
    }
    state
}

#[test]
fn wall_faces_carry_no_flow_after_many_steps_of_wind_forced_at_the_coast() {
    // The ticket's acceptance criterion: no flow ever appears at a boundary
    // face regardless of forcing, checked after every one of many steps rather
    // than at t = 0. The stress is applied at every face, walls included, so
    // only a boundary condition can keep this true.
    let basin = test_basin();
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let (mut solver, dt_s) = solver_for(basin, params);
    let wind =
        WindStressField::uniform(basin.grid(), TRADE_WIND_STRESS_X_PA, TRADE_WIND_STRESS_Y_PA);

    let mut state = OceanState::at_rest(basin.grid());
    assert_walls_at_rest(&state, "at rest");
    for step in 0..FORCED_RUN_STEPS {
        solver.step(&mut state, step as f64 * dt_s, |_t| &wind);
        assert_walls_at_rest(&state, &format!("after step {}", step + 1));
    }

    // The run has to have done something, or the criterion above is vacuous.
    // The steady damped balance of `wind_forcing.rs` puts this basin's
    // thermocline anomaly at about 17 m at its eastern and western ends under
    // this stress; 1 m is that floored to the order below, so the check is
    // about the run being alive rather than about the tilt's magnitude.
    let peak_m = peak_anomaly_m(&state);
    assert!(
        peak_m > 1.0,
        "the forced run barely moved the thermocline: peak |h| = {peak_m} m"
    );
}

#[test]
fn wall_faces_carry_no_flow_under_a_scenario_wind_resampled_at_every_stage() {
    // The other way into the solver: a `WindStress` re-sampled at each RK4
    // stage, which is the path a scenario takes. The stress is non-zero at the
    // coast in both components and changes sign within the run, so a wall face
    // is forced in both directions.
    let basin = test_basin();
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let (mut solver, dt_s) = solver_for(basin, params);
    let wind = GustsAtEveryFace {
        period_s: FORCED_RUN_STEPS as f64 * dt_s / 4.0,
    };

    let mut state = OceanState::at_rest(basin.grid());
    for step in 0..FORCED_RUN_STEPS {
        solver.step_forced_by(&mut state, step as f64 * dt_s, basin, &wind);
        assert_walls_at_rest(&state, &format!("after step {}", step + 1));
    }
}

#[test]
fn a_state_that_arrives_with_flow_at_the_walls_is_brought_to_rest() {
    // "No flow ever appears at a boundary face" has to hold for the state the
    // solver is handed, not only for the one it builds from rest: a hand-built
    // initial condition, or one restored from a run, must be brought onto the
    // boundary condition rather than integrated off it.
    let basin = test_basin();
    let params = pacific_params(STRONG_DAMPING_PER_S);
    let (mut solver, _) = solver_for(basin, params);
    let calm = WindStressField::calm(basin.grid());

    let mut state = OceanState::at_rest(basin.grid());
    for value in state.u_mut().as_mut_slice() {
        *value = 1.0;
    }
    for value in state.v_mut().as_mut_slice() {
        *value = -1.0;
    }

    solver.step(&mut state, 0.0, |_t| &calm);

    assert_walls_at_rest(&state, "after one step from a state with flow at the walls");
}

#[test]
fn a_stress_applied_at_a_wall_changes_nothing_in_the_basin() {
    // The boundary owns the wall velocities outright, so the `τ/(ρ₀·H)` a wall
    // face receives is discarded at every stage. Two runs that differ only in
    // their coastal stress must therefore agree bit for bit — which is what
    // makes the sampling rule T-03.1 put in `WindStressField::sample`
    // redundant rather than load-bearing, and this ticket free to drop it.
    let basin = test_basin();
    let params = pacific_params(STRONG_DAMPING_PER_S);

    let forced_at_the_coast = WindStressField::sampled(basin, &GaleAtEveryFace, 0.0);
    let calm_at_the_coast = WindStressField::sampled(basin, &CalmAtTheCoast { basin }, 0.0);
    assert_ne!(
        forced_at_the_coast, calm_at_the_coast,
        "the two fields must actually differ at the coast, or the comparison is vacuous"
    );

    let forced = run_under(basin, params, &forced_at_the_coast);
    let sheltered = run_under(basin, params, &calm_at_the_coast);

    assert_eq!(
        forced, sheltered,
        "a stress at the wall must not reach the interior"
    );
}

#[test]
fn an_undamped_closed_basin_conserves_its_volume_anomaly() {
    // The divergence theorem on the C-grid (this file's header): with no
    // normal flow at any wall and no damping, `Σ h·Δx·Δy` is conserved
    // exactly, however violently the wind stirs the interior. It starts at
    // zero, so it must stay at zero.
    //
    // This is the check that says *where* the condition has to be applied. RK4
    // evaluates the continuity equation at four stage states per step; a
    // condition applied only to the completed step would leave each of those
    // stages with a wall velocity of order `dt·τ/(ρ₀·H)`, and the volume those
    // stages leak through the coast does not cancel.
    let basin = test_basin();
    let params = pacific_params(UNDAMPED_PER_S);
    let wind = WindStressField::sampled(basin, &GaleTiltedAcrossTheBasin { basin }, 0.0);
    let (nx, ny) = (basin.grid().nx(), basin.grid().ny());
    assert_ne!(
        wind.tau_x_pa().get(0, 0),
        wind.tau_x_pa().get(nx, 0),
        "the two zonal walls must be forced differently, or the budget cannot fail"
    );
    assert_ne!(
        wind.tau_y_pa().get(0, 0),
        wind.tau_y_pa().get(0, ny),
        "the two meridional walls must be forced differently, or the budget cannot fail"
    );

    let state = run_under(basin, params, &wind);

    let drift_m3 = volume_anomaly_m3(&state, basin.spacing()).abs();
    let scale_m3 = volume_scale_m3(&state, basin.spacing());
    // The identity is exact in exact arithmetic, so the only admissible
    // residual is round-off: 500 steps of four stages each, every stage
    // accumulating a few ulps of `f64` (ε ≈ 2.2×10⁻¹⁶) into a sum over 128
    // cells, is of order 10⁻¹² of the volume being cancelled. That is the
    // bound, taken relative to the scale of what cancels.
    let tolerance_m3 = 1.0e-12 * scale_m3;
    assert!(
        drift_m3 <= tolerance_m3,
        "a closed, undamped basin leaked {drift_m3} m³ of volume through its coasts, \
         against a bound of {tolerance_m3} m³ on a total |h| volume of {scale_m3} m³"
    );
}
