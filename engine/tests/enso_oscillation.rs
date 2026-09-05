//! T-12.3 — does the closed Bjerknes loop oscillate, and does it oscillate at
//! the period the delayed-oscillator theory predicts?
//!
//! T-12.1 gave the mixed layer an equation, T-12.2 gave the atmosphere an
//! answer to it, and this file asks the question the two were built for: run
//! the fully coupled basin for decades and watch an equatorial-Pacific SST
//! index. `docs/enso-oscillation-report.md` is the written record of what these
//! runs measured; this file is the suite the record is made of.
//!
//! # What is predicted, before anything is run
//!
//! The delayed oscillator (Suarez & Schopf, *J. Atmos. Sci.* 45, 1988;
//! Battisti & Hirst, *J. Atmos. Sci.* 46, 1989) reduces the coupled basin to
//!
//! ```text
//! dT/dt = a·T − b·T(t − δ)
//! ```
//!
//! — an instantaneous positive Bjerknes feedback `a`, and a negative one
//! delayed by `δ`, the time an off-equatorial Rossby wave takes to reach the
//! western wall and return as a reflected Kelvin wave. Two consequences of
//! that equation are what this file tests, and neither of them is a number
//! read out of the engine:
//!
//! - **A threshold.** Substituting `T = e^{(σ+iω)t}` gives `σ + iω = a −
//!   b·e^{−(σ+iω)δ}`, so at the margin `cos(ωδ) = a/b` and `ω = b·sin(ωδ)`.
//!   Both `a` and `b` are proportional to the feedback strength `μ` of
//!   [`WindResponseParams`], while the damping that opposes them is not, so
//!   there is a critical `μ` below which a perturbation decays and above which
//!   it grows.
//! - **A period set by wave transit, not by heat.** `δ` is `(3·L_R + L_K)/c`
//!   for a Rossby leg `L_R` and a Kelvin leg `L_K`, because the gravest Rossby
//!   mode travels at `c/3` and the Kelvin wave at `c` — the two speeds
//!   `docs/validation-report.md` already validated. Every path through the
//!   basin is some fraction of `L/c`, so `δ ∝ L/c` and, since the phase
//!   condition above fixes `ωδ`, the **period is proportional to `L/c`** with
//!   a coefficient of order one. This is the falsifiable claim: change the
//!   basin's width or its wave speed, and the period must follow. A period
//!   set by the mixed layer's thermal relaxation (125 days), by the ocean's
//!   Rayleigh damping (2.5 years), or by a numerical artefact of the grid
//!   would not move at all.
//!
//! # What is measured
//!
//! Each experiment holds the basin under steady alizés until the coupled state
//! settles ([`SPIN_UP_YEARS`]), adds a uniform 1 K warm anomaly to the mixed
//! layer — an idealised El Niño — and records the eastern equatorial SST index
//! for [`WINDOW_YEARS`] more, as its departure from the value the spin-up left
//! it at. [`IndexRecord`] then reports three things about that series: its
//! amplitude, its exponential growth rate, and the period of its dominant
//! spectral peak.
//!
//! # The headline result, stated here because it is a negative one
//!
//! The model oscillates, the oscillation has a threshold in `μ`, and its
//! period tracks `L/c` across a factor of two in each — all three predictions
//! hold. The period itself is about `5·L/c`, which for this basin is close to
//! one year and therefore **short of the observed 2–7 year ENSO band**, the
//! range the ticket's acceptance criterion names.
//! [`the_period_falls_short_of_the_observed_enso_band`] records that
//! discrepancy as an assertion rather than as prose, so that a change to the
//! physics which fixes it fails this file loudly.
//! `docs/enso-oscillation-report.md` § *Why the period is short* is the
//! diagnosis.

use std::f64::consts::TAU;
use std::sync::OnceLock;

use engine::sst::{SstParams, DEFAULT_SURFACE_DRAG_PER_S};
use engine::wind_response::{CoupledWind, SstWindResponse, WindResponseParams};
use engine::{
    run_scenario, Basin, BasinBounds, BetaPlane, CompositeWind, OceanState, PhysicalParams,
    RunReader, ScenarioConfig, Solver, SteadyTradeWinds, EQUATORIAL_BETA_PER_M_PER_S, H_STAGGERING,
    SEAWATER_REFERENCE_DENSITY_KG_PER_M3, TROPICAL_YEAR_S,
};

mod common;

use common::ScratchDir;

/// This file's ticket, which labels the directories it leaves in the system
/// temp directory.
const TICKET: &str = "t123";

// ---------------------------------------------------------------------------
// The ocean, the mixed layer and the atmosphere these experiments run in
// ---------------------------------------------------------------------------

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s² (Gill, *Atmosphere–Ocean Dynamics*, ch. 11; Cane & Sarachik 1981) —
/// the same value `docs/validation-report.md` states for every Epic 07 suite.
const REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H`, in metres — the canonical 150 m upper layer of
/// the same 1.5-layer configuration.
const MEAN_THERMOCLINE_DEPTH_M: f64 = 150.0;
/// Seconds in a day.
const SECONDS_PER_DAY: f64 = 86_400.0;

/// Rayleigh damping timescale `1/r` of the ocean dynamics, in tropical years.
///
/// Zebiak & Cane (*Mon. Wea. Rev.* 115, 1987, § 2a) damp their equatorial
/// ocean model on 2.5 years, and Battisti (*J. Atmos. Sci.* 45, 1988) — the
/// study that first found the delayed oscillation in that model — inherits it.
/// It is far weaker than the 100-day damping the steady-state validations of
/// Epic 07 use, and it has to be: a Rossby wave needs `3L/c ≈ 0.6 years` to
/// cross this basin, so a damping that removed it in 100 days would remove the
/// delayed branch of the feedback along with it, and the loop this file is
/// about would not close.
const OCEAN_DAMPING_TIMESCALE_YEARS: f64 = 2.5;
/// Rayleigh damping `r`, in s⁻¹.
const RAYLEIGH_DAMPING_PER_S: f64 = 1.0 / (OCEAN_DAMPING_TIMESCALE_YEARS * TROPICAL_YEAR_S);

/// Mixed-layer depth `H_m`, in metres (Zebiak & Cane 1987, § 2b) — the value
/// the T-12.1 and T-12.2 suites already run at.
const MIXED_LAYER_DEPTH_M: f64 = 50.0;
/// Zonal gradient of the mean SST, `∂T̄/∂x`, in K/m: about 7 K of warm-pool to
/// cold-tongue contrast across this basin's 17 800 km, negative because the
/// equatorial Pacific cools eastward.
const MEAN_ZONAL_SST_GRADIENT_K_PER_M: f64 = -4.0e-7;
/// Sensitivity `γ = ∂T_sub/∂h` of the entrained water to the thermocline depth
/// anomaly, in K/m (Zebiak & Cane 1987, § 2c).
const SUBSURFACE_SENSITIVITY_K_PER_M: f64 = 0.1;
/// Thermal damping `ε_T` of an SST anomaly, in s⁻¹: a 125-day relaxation
/// (Zebiak & Cane 1987, § 2b).
const THERMAL_DAMPING_PER_S: f64 = 1.0 / (125.0 * SECONDS_PER_DAY);

/// Equatorial trade-wind stress, in Pa. Easterly, so negative (`CONTEXT.md`,
/// *Wind stress*), and the same 0.05 Pa the T-12.2 suite drives the basin
/// with.
const TRADE_WIND_STRESS_PA: f64 = -0.05;

// ---------------------------------------------------------------------------
// The basin, and how long each experiment runs
// ---------------------------------------------------------------------------

/// Western boundary of every basin here, in degrees east — the Maritime
/// Continent edge of the Pacific (`CONTEXT.md`, *Basin*).
const WESTERN_LONGITUDE_DEG: f64 = 120.0;
/// Eastern boundary of the reference basin, in degrees east: 80°W, the South
/// American coast.
const PACIFIC_EASTERN_LONGITUDE_DEG: f64 = 280.0;
/// Southern and northern boundaries, in degrees north.
const SOUTHERN_LATITUDE_DEG: f64 = -25.0;
const NORTHERN_LATITUDE_DEG: f64 = 25.0;

/// Cell size, in degrees.
///
/// Two degrees is 222.6 km, which puts 1.55 cells in the 345 km deformation
/// radius of `docs/validation-report.md`. That is coarse beside the Epic 07
/// wave suites, which resolve `Le` seven to fourteen times over — and
/// deliberately so: what is measured here is a basin-crossing *time*, an
/// integral over the whole waveguide rather than the shape of one wave's
/// meridional profile, and every experiment is compared against another on the
/// same grid. The tolerance each comparison carries says what that costs.
const RESOLUTION_DEG: f64 = 2.0;

/// Timestep, in seconds.
///
/// Bounded by rotation rather than by the gravity waves: `|f| = β·y` reaches
/// 6.4 × 10⁻⁵ s⁻¹ at 25°, whose inertial oscillation RK4 follows up to 35 400 s
/// (ADR-0007), where the CFL bound of a 222.6 km cell at `c = 2.74 m/s` is
/// 65 000 s. Halving it to 12 000 s moves the measured period by less than one
/// part in 10⁴, so nothing below is timestep-limited.
const DT_S: f64 = 30_000.0;

/// Years the basin is held under steady alizés before it is perturbed.
///
/// Eight `e`-foldings of [`OCEAN_DAMPING_TIMESCALE_YEARS`], which leaves
/// 3 × 10⁻⁴ of the switch-on transient behind.
const SPIN_UP_YEARS: f64 = 20.0;
/// Years the index is recorded for after the perturbation.
///
/// The whole record is what [`IndexRecord::growth_rate_per_year`] reads,
/// because the perturbation is largest and cleanest right after it is applied.
const RECORD_YEARS: f64 = SETTLING_YEARS + WINDOW_YEARS;
/// Years of the record that belong to the settling of the step.
///
/// The 1 K anomaly is a step, and a step is every mode of the coupled basin at
/// once; the amplitude and the period are properties of the one that outlives
/// the others. Twelve years is five `e`-foldings of
/// [`OCEAN_DAMPING_TIMESCALE_YEARS`] and thirty-five of the mixed layer's
/// 125-day relaxation, so what the settled window opens on is that mode and
/// not the step.
const SETTLING_YEARS: f64 = 12.0;
/// Years of the record the amplitude and the period are read from.
///
/// Thirty basin-crossing times of the reference ocean, which is what makes a
/// spectral peak sharp enough to locate, and four times the slow edge of the
/// observed ENSO band, so that a period anywhere in that band would be seen.
const WINDOW_YEARS: f64 = 32.0;
/// Samples of the index taken across [`WINDOW_YEARS`].
///
/// Two thousand, so the shortest period this file resolves is 0.03 years and
/// the ones it measures are sampled fifty times per cycle.
const SAMPLES: usize = 2_000;

/// Amplitude of the warm anomaly the mixed layer is perturbed with, in kelvin.
///
/// One kelvin, uniform: the scale of an observed El Niño's SST anomaly, and
/// large enough that a decay of ten orders of magnitude is still above the
/// rounding of the sums that form the index.
const PERTURBATION_K: f64 = 1.0;

/// Half-width of the SST index's meridional band, in degrees of latitude.
///
/// Five degrees either side of the equator, the band of the observed Niño
/// indices, and about 1.6 deformation radii — wide enough to hold the
/// equatorial waveguide the coupling lives on.
const INDEX_HALF_WIDTH_DEG: f64 = 5.0;
/// Fraction of the basin, measured from the eastern wall, the index averages
/// over.
///
/// The eastern third. For the reference basin that is 227°E–280°E, which is
/// most of the observed Niño-3 region (150°W–90°W); stated as a fraction
/// rather than as two longitudes so that it follows the basin when
/// [`the_period_scales_with_the_basin_crossing_time`] makes the basin
/// narrower.
const INDEX_EASTERN_FRACTION: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// One experiment
// ---------------------------------------------------------------------------

/// One coupled run: an ocean, a basin, and how hard the atmosphere answers.
///
/// The three fields are exactly what the tests below vary. Everything else —
/// the mixed layer, the trades, the grid, the timestep, the run length — is
/// the same in every experiment, so a difference between two of them is a
/// difference in one of these.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Experiment {
    /// Reduced gravity `g'`, in m/s². Sets `c = √(g'·H)`, and nothing else the
    /// coupling reads.
    reduced_gravity_m_per_s2: f64,
    /// Eastern boundary, in degrees east. Sets the basin width `L`.
    eastern_longitude_deg: f64,
    /// Feedback strength `μ`, in Pa/K.
    feedback_strength_pa_per_k: f64,
}

impl Experiment {
    /// The reference experiment: the equatorial Pacific of `CONTEXT.md` at
    /// `feedback_strength_pa_per_k`.
    const fn pacific(feedback_strength_pa_per_k: f64) -> Self {
        Self {
            reduced_gravity_m_per_s2: REDUCED_GRAVITY_M_PER_S2,
            eastern_longitude_deg: PACIFIC_EASTERN_LONGITUDE_DEG,
            feedback_strength_pa_per_k,
        }
    }

    /// The same experiment in a basin reaching only to `eastern_longitude_deg`.
    const fn narrowed_to(self, eastern_longitude_deg: f64) -> Self {
        Self {
            eastern_longitude_deg,
            ..self
        }
    }

    /// The same experiment in an ocean of reduced gravity
    /// `reduced_gravity_m_per_s2`, and so of wave speed `√(g'·H)`.
    const fn with_reduced_gravity(self, reduced_gravity_m_per_s2: f64) -> Self {
        Self {
            reduced_gravity_m_per_s2,
            ..self
        }
    }

    /// The basin this experiment runs in.
    fn basin(self) -> Basin {
        BasinBounds::new(
            WESTERN_LONGITUDE_DEG,
            self.eastern_longitude_deg,
            SOUTHERN_LATITUDE_DEG,
            NORTHERN_LATITUDE_DEG,
            RESOLUTION_DEG,
        )
        .expect("the experiment's boundaries are a whole number of cells apart")
        .basin()
    }

    /// The ocean parameters.
    fn params(self) -> PhysicalParams {
        PhysicalParams::new(
            self.reduced_gravity_m_per_s2,
            MEAN_THERMOCLINE_DEPTH_M,
            RAYLEIGH_DAMPING_PER_S,
            EQUATORIAL_BETA_PER_M_PER_S,
            SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
        )
        .expect("the experiment's ocean is a physical one")
    }

    /// The time a Kelvin wave takes to cross this basin, `L/c`, in tropical
    /// years.
    ///
    /// Written from the definitions of `CONTEXT.md` — `c = √(g'·H)` and the
    /// basin's own width — and never asked of the engine, because it is the
    /// quantity the measured period is predicted from.
    fn basin_crossing_years(self) -> f64 {
        let basin = self.basin();
        let wave_speed_m_per_s = (self.reduced_gravity_m_per_s2 * MEAN_THERMOCLINE_DEPTH_M).sqrt();
        basin.zonal_extent_m() / wave_speed_m_per_s / TROPICAL_YEAR_S
    }

    /// Run it, and report the eastern SST index it produced.
    fn record(self) -> IndexRecord {
        let basin = self.basin();
        let grid = basin.grid();
        let params = self.params();
        let plane = BetaPlane::of_basin(params, basin);
        let mut solver =
            Solver::coupled_to_sst(grid, basin.spacing(), params, plane, DT_S, sst_params())
                .expect("the timestep is inside both stability bounds");
        let index = EasternSstIndex::over(basin);

        let mut state = OceanState::at_rest_with_sst_anomaly(grid);
        let mut forcing = self.forcing();
        let spin_up_steps = (SPIN_UP_YEARS * TROPICAL_YEAR_S / DT_S) as u64;
        for step in 0..spin_up_steps {
            solver.step_with_forcing(&mut state, step as f64 * DT_S, &mut forcing);
        }
        let settled_index_k = index.of(&state);

        for value_k in state
            .sst_anomaly_k_mut()
            .expect("a coupled run's state carries `T'`")
            .as_mut_slice()
        {
            *value_k += PERTURBATION_K;
        }

        let record_steps = (RECORD_YEARS * TROPICAL_YEAR_S / DT_S) as u64;
        let window_steps = (WINDOW_YEARS * TROPICAL_YEAR_S / DT_S) as u64;
        let sample_every = (window_steps / SAMPLES as u64).max(1);
        let mut departures_k = Vec::with_capacity((record_steps / sample_every) as usize + 1);
        for step in 0..=record_steps {
            if step % sample_every == 0 {
                departures_k.push(index.of(&state) - settled_index_k);
            }
            if step < record_steps {
                solver.step_with_forcing(
                    &mut state,
                    (spin_up_steps + step) as f64 * DT_S,
                    &mut forcing,
                );
            }
        }
        let sample_interval_years = sample_every as f64 * DT_S / TROPICAL_YEAR_S;
        IndexRecord {
            settled_from: (SETTLING_YEARS / sample_interval_years) as usize,
            sample_interval_years,
            departures_k,
        }
    }

    /// The forcing of this experiment: steady alizés, plus the atmosphere's
    /// answer to the SST anomaly at this feedback strength.
    fn forcing(self) -> CoupledWind<CompositeWind> {
        let basin = self.basin();
        let response = WindResponseParams::new(
            self.feedback_strength_pa_per_k,
            engine::DEFAULT_WIND_RESPONSE_MERIDIONAL_SCALE_M,
        )
        .expect("a non-negative strength and the default scale are a valid response");
        CoupledWind::new(
            basin,
            CompositeWind::new()
                .with(SteadyTradeWinds::uniform(TRADE_WIND_STRESS_PA).expect("an easterly stress")),
            SstWindResponse::new(basin, response),
        )
    }
}

/// The mixed layer every experiment carries.
fn sst_params() -> SstParams {
    SstParams::new(
        MIXED_LAYER_DEPTH_M,
        DEFAULT_SURFACE_DRAG_PER_S,
        MEAN_ZONAL_SST_GRADIENT_K_PER_M,
        SUBSURFACE_SENSITIVITY_K_PER_M,
        THERMAL_DAMPING_PER_S,
    )
    .expect("these are the standard Zebiak-Cane mixed-layer parameters")
}

// ---------------------------------------------------------------------------
// The index
// ---------------------------------------------------------------------------

/// The eastern equatorial SST index of one basin: `T'` averaged over the
/// eastern [`INDEX_EASTERN_FRACTION`] of it, within
/// [`INDEX_HALF_WIDTH_DEG`] of the equator.
///
/// The cells are found once and reused, so recording an index costs one pass
/// over them rather than a search per sample.
struct EasternSstIndex {
    /// Cell-center columns inside the box.
    columns: Vec<usize>,
    /// Cell-center rows inside the box.
    rows: Vec<usize>,
}

impl EasternSstIndex {
    /// The index box of `basin`.
    ///
    /// # Panics
    /// If the box is empty, which would mean a basin too small or too coarse
    /// to hold one — a broken experiment rather than a failed measurement.
    fn over(basin: Basin) -> Self {
        let grid = basin.grid();
        let western_edge_of_box_m = basin.zonal_extent_m() * (1.0 - INDEX_EASTERN_FRACTION);
        let half_width_m = INDEX_HALF_WIDTH_DEG * engine::basin::METRES_PER_DEGREE_OF_ARC;
        let columns: Vec<usize> = (0..grid.nx())
            .filter(|&i| basin.x_of_column_m(H_STAGGERING, i) >= western_edge_of_box_m)
            .collect();
        let rows: Vec<usize> = (0..grid.ny())
            .filter(|&j| basin.y_of_row_m(H_STAGGERING, j).abs() <= half_width_m)
            .collect();
        assert!(
            !columns.is_empty() && !rows.is_empty(),
            "the index box of this basin holds no cells: {} columns by {} rows",
            columns.len(),
            rows.len()
        );
        Self { columns, rows }
    }

    /// The index of `state`, in kelvin.
    ///
    /// # Panics
    /// If `state` is not a coupled state; an index of a run with no `T'` is
    /// the caller asking the wrong question.
    fn of(&self, state: &OceanState) -> f64 {
        self.of_field(
            state
                .sst_anomaly_k()
                .expect("a coupled run's state carries `T'`")
                .as_slice(),
            state.grid().nx(),
        )
    }

    /// The index of an SST anomaly laid out row-major in `nx` columns —
    /// the shape both an [`OceanState`]'s field and a written frame's slice
    /// have.
    fn of_field(&self, sst_anomaly_k: &[f64], nx: usize) -> f64 {
        // Rows summed in row order and then added, so the result is fixed by
        // the box rather than by an iteration order
        // (CODING_STANDARDS.md § *Correctness and failure*).
        let total_k: f64 = self
            .rows
            .iter()
            .map(|&j| {
                self.columns
                    .iter()
                    .map(|&i| sst_anomaly_k[j * nx + i])
                    .sum::<f64>()
            })
            .sum();
        total_k / (self.rows.len() * self.columns.len()) as f64
    }
}

// ---------------------------------------------------------------------------
// What a recorded index says
// ---------------------------------------------------------------------------

/// The SST index of one experiment over [`RECORD_YEARS`] after the
/// perturbation, as its departure from the value the spin-up left it at.
#[derive(Debug, Clone)]
struct IndexRecord {
    /// Model time between two samples, in tropical years.
    sample_interval_years: f64,
    /// Index of the first sample past [`SETTLING_YEARS`] — where the settled
    /// window [`IndexRecord::settled`] starts.
    settled_from: usize,
    /// The departures themselves, in kelvin, evenly spaced in model time.
    departures_k: Vec<f64>,
}

impl IndexRecord {
    /// The part of the record the step has settled out of: the last
    /// [`WINDOW_YEARS`] of it.
    ///
    /// The amplitude and the period are read from here, because both are
    /// properties of the surviving mode rather than of the perturbation that
    /// excited it. The growth rate is not: it is read from the whole record,
    /// where the mode is largest and furthest above the residue of the
    /// spin-up.
    fn settled(&self) -> &[f64] {
        &self.departures_k[self.settled_from..]
    }

    /// Half the peak-to-trough range of the settled window, in kelvin.
    fn amplitude_k(&self) -> f64 {
        let (min, max) = self
            .settled()
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &value| {
                (lo.min(value), hi.max(value))
            });
        (max - min) / 2.0
    }

    /// The exponential growth rate of the record, in reciprocal years.
    ///
    /// The logarithm of the ratio of the root-mean-square departure over the
    /// second half of the window to the first, divided by the time between
    /// their midpoints. For a mode `e^{σt}·cos(ωt + φ)` sampled over many
    /// cycles this is `σ` — the oscillation contributes the same mean square
    /// to both halves — and for a mode that has saturated it is zero.
    ///
    /// The whole record, not [`IndexRecord::settled`]: a strongly damped mode
    /// has already fallen to the residue of the spin-up by the time the
    /// settled window opens, and a rate measured there would be that residue's
    /// and not the mode's.
    fn growth_rate_per_year(&self) -> f64 {
        let half = self.departures_k.len() / 2;
        let rms = |samples: &[f64]| {
            (samples.iter().map(|value| value * value).sum::<f64>() / samples.len() as f64).sqrt()
        };
        let early = rms(&self.departures_k[..half]);
        let late = rms(&self.departures_k[half..]);
        let separation_years = half as f64 * self.sample_interval_years;
        (late / early).ln() / separation_years
    }

    /// The period of the record's strongest spectral line, in tropical years.
    ///
    /// The mean is removed and the discrete Fourier power `|Σ y·e^{−2πift}|²`
    /// evaluated on a grid of frequencies four times finer than the record's
    /// own resolution `1/T`; the largest is refined by fitting a parabola to it
    /// and its two neighbours, which is the standard sub-bin peak estimate for
    /// a spectral line. Only periods between [`Self::MIN_PERIOD_YEARS`] and
    /// [`Self::MAX_PERIOD_YEARS`] are considered.
    ///
    /// # Panics
    /// If the record is too short to hold a whole period of the longest
    /// candidate, which is a mis-configured experiment rather than a failed
    /// measurement.
    fn dominant_period_years(&self) -> f64 {
        let samples = self.settled();
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let span_years = (samples.len() - 1) as f64 * self.sample_interval_years;
        assert!(
            span_years >= 2.0 * Self::MAX_PERIOD_YEARS,
            "a {span_years} year record cannot resolve a {} year period",
            Self::MAX_PERIOD_YEARS
        );
        // Four times the record's own frequency resolution `1/T`: enough to
        // put a peak within a bin of its true position, which the parabola
        // then refines.
        let frequency_step_per_year = 1.0 / (4.0 * span_years);
        let first = (1.0 / (Self::MAX_PERIOD_YEARS * frequency_step_per_year)).ceil() as usize;
        let last = (1.0 / (Self::MIN_PERIOD_YEARS * frequency_step_per_year)).floor() as usize;
        let power_at = |bin: usize| {
            let frequency_per_year = bin as f64 * frequency_step_per_year;
            let (mut real, mut imaginary) = (0.0, 0.0);
            for (index, value) in samples.iter().enumerate() {
                let phase = TAU * frequency_per_year * index as f64 * self.sample_interval_years;
                real += (value - mean) * phase.cos();
                imaginary -= (value - mean) * phase.sin();
            }
            real * real + imaginary * imaginary
        };
        let mut peak_bin = first;
        let mut peak_power = f64::NEG_INFINITY;
        for bin in first..=last {
            let power = power_at(bin);
            if power > peak_power {
                peak_power = power;
                peak_bin = bin;
            }
        }
        let (below, above) = (power_at(peak_bin - 1), power_at(peak_bin + 1));
        let curvature = below - 2.0 * peak_power + above;
        let offset_bins = if curvature == 0.0 {
            0.0
        } else {
            0.5 * (below - above) / curvature
        };
        1.0 / ((peak_bin as f64 + offset_bins) * frequency_step_per_year)
    }

    /// Shortest period the search considers, in years — six times the
    /// timestep-resolved gravity-wave ringing is not what is being looked for,
    /// and 0.2 years is below every basin crossing time in this file.
    const MIN_PERIOD_YEARS: f64 = 0.2;
    /// Longest period the search considers, in years. A quarter of
    /// [`WINDOW_YEARS`], so that the slowest line it can report still has four
    /// whole cycles in the record — and, at eight years, it is past the slow
    /// edge of the observed ENSO band, so a model that did oscillate there
    /// would be seen to.
    const MAX_PERIOD_YEARS: f64 = WINDOW_YEARS / 4.0;
}

// ---------------------------------------------------------------------------
// The feedback strengths the tests are run at
// ---------------------------------------------------------------------------

/// Feedback strengths `μ`, in Pa/K, the sensitivity study sweeps.
///
/// Nothing about the physics picks these; they bracket the threshold the
/// theory says exists, which is what
/// [`the_growth_rate_rises_with_the_feedback_strength`] and
/// [`the_perturbation_decays_below_the_threshold_and_grows_above_it`] are
/// about. `docs/enso-oscillation-report.md` tabulates what each one measured.
const SUBCRITICAL_STRENGTHS_PA_PER_K: [f64; 3] = [0.02, 0.04, 0.06];
/// The feedback strength the oscillation is characterised at, in Pa/K — above
/// the threshold the sweep brackets, and the strength every period below is
/// measured at.
const SUPERCRITICAL_STRENGTH_PA_PER_K: f64 = 0.08;
/// The open loop: `μ = 0` is the prescribed-wind model of T-12.1.
const OPEN_LOOP_STRENGTH_PA_PER_K: f64 = 0.0;

/// The reference oscillating run, shared by every test that reads its period.
fn reference_record() -> &'static IndexRecord {
    static RECORD: OnceLock<IndexRecord> = OnceLock::new();
    RECORD.get_or_init(|| Experiment::pacific(SUPERCRITICAL_STRENGTH_PA_PER_K).record())
}

/// The subcritical sweep, shared by the two tests that read it.
fn subcritical_records() -> &'static [IndexRecord; 3] {
    static RECORDS: OnceLock<[IndexRecord; 3]> = OnceLock::new();
    RECORDS.get_or_init(|| {
        SUBCRITICAL_STRENGTHS_PA_PER_K.map(|strength| Experiment::pacific(strength).record())
    })
}

// ---------------------------------------------------------------------------
// The threshold — the sensitivity half of the acceptance criteria
// ---------------------------------------------------------------------------

#[test]
fn an_open_loop_run_returns_to_the_state_the_trades_hold_it_at() {
    // `μ = 0` is the model of T-12.1: `T'` reads the ocean and nothing reads
    // `T'`. With no feedback the only terms acting on the anomaly are the
    // entrainment and the thermal damping, both of which are relaxations, so
    // the 1 K perturbation must decay and the basin must return to the state
    // the trades hold it at. The delayed oscillator has nothing to oscillate
    // with here: `a = b = 0`.
    let record = Experiment::pacific(OPEN_LOOP_STRENGTH_PA_PER_K).record();

    // Nothing in the open loop relaxes more slowly than the ocean's Rayleigh
    // damping at 2.5 years, so by the time the window opens —
    // [`SETTLING_YEARS`] after the perturbation, five `e`-foldings of it — at
    // most `e^{−4.8} = 8×10⁻³` of the 1 K anomaly is left, and it goes on
    // decaying through the window. The bound is that factor, loosened to
    // 10⁻² so that it states "the perturbation is gone" and not a decay rate
    // the equation does not promise.
    let decayed_bound_k = PERTURBATION_K * 1.0e-2;
    assert!(
        record.amplitude_k() < decayed_bound_k,
        "the open loop still carries {} K of the {PERTURBATION_K} K perturbation \
         {SETTLING_YEARS} years on, more than the {decayed_bound_k} K a 2.5-year relaxation \
         leaves",
        record.amplitude_k()
    );
}

#[test]
fn the_growth_rate_rises_with_the_feedback_strength() {
    // The delayed oscillator's `a` and `b` are both proportional to `μ` and
    // the damping opposing them is not, so `σ` is an increasing function of
    // `μ`. This is the statement that makes "there is a threshold" more than
    // an accident of the two strengths
    // [`the_perturbation_decays_below_the_threshold_and_grows_above_it`]
    // happens to pick: the rate has to climb monotonically towards it.
    //
    // No tolerance: the assertion is on the order of three numbers.
    let records = subcritical_records();
    for (window, strengths) in records
        .windows(2)
        .zip(SUBCRITICAL_STRENGTHS_PA_PER_K.windows(2))
    {
        let (weaker, stronger) = (
            window[0].growth_rate_per_year(),
            window[1].growth_rate_per_year(),
        );
        assert!(
            stronger > weaker,
            "raising μ from {} to {} Pa/K moved the growth rate from {weaker} to {stronger} \
             per year, the wrong way",
            strengths[0],
            strengths[1]
        );
    }
    // And every one of them is a decay: the sweep is the subcritical side of
    // the threshold, which is what makes the comparison below a bracket.
    for (record, strength) in records.iter().zip(SUBCRITICAL_STRENGTHS_PA_PER_K) {
        assert!(
            record.growth_rate_per_year() < 0.0,
            "μ = {strength} Pa/K was expected to be subcritical, but its perturbation grew at \
             {} per year",
            record.growth_rate_per_year()
        );
    }
}

#[test]
fn the_perturbation_decays_below_the_threshold_and_grows_above_it() {
    // The acceptance criterion's sensitivity clause, and the delayed
    // oscillator's instability threshold: below a critical feedback strength a
    // perturbation dies, above it the coupled mode grows. The comparison is
    // made against the perturbation the run was given rather than against any
    // measured amplitude — a run that ends up smaller than what it started
    // with decayed, and one that ends up larger grew, and neither statement
    // needs a tolerance.
    // The middle of the sweep rather than its top: the strength closest to
    // the threshold decays so slowly that the window still holds part of the
    // step it was perturbed with, and "smaller than the perturbation" would
    // then be a statement about the settling time rather than about the sign
    // of the growth rate — which is what
    // [`the_growth_rate_rises_with_the_feedback_strength`] asserts on
    // directly, and where the tight bracket on the threshold lives.
    const DECAYING_STRENGTH_INDEX: usize = 1;
    let decaying = &subcritical_records()[DECAYING_STRENGTH_INDEX];
    assert!(
        decaying.amplitude_k() < PERTURBATION_K,
        "μ = {} Pa/K left {} K of a {PERTURBATION_K} K perturbation, so it is not below the \
         threshold",
        SUBCRITICAL_STRENGTHS_PA_PER_K[DECAYING_STRENGTH_INDEX],
        decaying.amplitude_k()
    );
    assert!(
        reference_record().amplitude_k() > PERTURBATION_K,
        "μ = {SUPERCRITICAL_STRENGTH_PA_PER_K} Pa/K left only {} K of a {PERTURBATION_K} K \
         perturbation, so it is not above the threshold",
        reference_record().amplitude_k()
    );
}

// ---------------------------------------------------------------------------
// The oscillation itself
// ---------------------------------------------------------------------------

#[test]
fn the_supercritical_run_settles_into_a_self_sustained_oscillation() {
    // Above the threshold the linearised mode grows exponentially; what stops
    // it is the one nonlinearity the coupled model has, the `w⁺ = max(w, 0)`
    // clamp of `crate::sst` acting on a wind that now depends on the state.
    // The result is a limit cycle rather than a runaway, and the test of that
    // is that the amplitude stops changing: an unsaturated mode at this
    // strength grows by orders of magnitude over the window, so a second half
    // within a factor of two of the first is a cycle and not growth.
    let settled = reference_record().settled();
    let half = settled.len() / 2;
    let peak = |samples: &[f64]| samples.iter().fold(0.0_f64, |a, &b| a.max(b.abs()));
    let early = peak(&settled[..half]);
    let late = peak(&settled[half..]);
    assert!(
        early > 0.0 && (0.5..=2.0).contains(&(late / early)),
        "the oscillation's peak went from {early} K to {late} K across the window, which is \
         not a settled cycle"
    );
}

#[test]
fn the_period_scales_with_the_basin_crossing_time() {
    // The delayed oscillator's central claim, and the only one about the
    // period that can be made before the model is run: the delay is
    // `(3·L_R + L_K)/c` — Rossby legs at `c/3`, Kelvin legs at `c`, the two
    // speeds `docs/validation-report.md` validated — so it is proportional to
    // the basin crossing time `L/c`, and the phase condition `cos(ωδ) = a/b`
    // makes the period proportional to it too.
    //
    // Three oceans are compared against the reference: one basin three
    // quarters as wide, and two with the wave speed multiplied by `√2` and by
    // 2. Every other parameter — the mixed layer's 125-day relaxation, the
    // ocean's 2.5-year damping, the atmosphere's width, the grid — is held
    // fixed, so the alternative hypotheses this discriminates against (that
    // the period is a thermal relaxation time, a damping time, or an artefact
    // fixed by the grid) all predict a ratio of exactly one.
    //
    // The tolerance is 25% on the ratio. It is not a precision claim: the
    // period is `4δ` only in the purely-delayed limit, and `a/b` is not held
    // fixed when the basin or the wave speed changes, so theory gives a
    // scaling with a coefficient of order one rather than an equality. What
    // the band has to do is separate the prediction from its alternatives, and
    // the smallest change it is asked to see is a factor of 0.71 — nearly
    // three times the band away from the ratio of one the alternatives
    // predict. It is also the order of the largest numerical term in the
    // comparison: the meridional truncation `(Δy/Le)²` of the Epic 07 budget
    // is 0.42 on this grid at `c = 2.74 m/s` and 0.21 at `c = 5.48 m/s`, and
    // the 0.21 difference between them does not cancel between two oceans.
    const RATIO_TOLERANCE: f64 = 0.25;

    let reference = Experiment::pacific(SUPERCRITICAL_STRENGTH_PA_PER_K);
    let reference_period_years = reference_record().dominant_period_years();
    let comparisons = [
        (
            "a basin narrowed to 120° of longitude",
            reference.narrowed_to(240.0),
        ),
        (
            "an ocean with `c` multiplied by √2",
            reference.with_reduced_gravity(2.0 * REDUCED_GRAVITY_M_PER_S2),
        ),
        (
            "an ocean with `c` doubled",
            reference.with_reduced_gravity(4.0 * REDUCED_GRAVITY_M_PER_S2),
        ),
    ];

    for (description, experiment) in comparisons {
        let predicted_ratio = experiment.basin_crossing_years() / reference.basin_crossing_years();
        let measured_ratio = experiment.record().dominant_period_years() / reference_period_years;
        let relative_error = (measured_ratio - predicted_ratio).abs() / predicted_ratio;
        assert!(
            relative_error < RATIO_TOLERANCE,
            "{description} crosses in {predicted_ratio} of the reference basin's `L/c`, but \
             oscillates at {measured_ratio} of its period — a relative departure of \
             {relative_error}"
        );
    }
}

#[test]
fn the_period_falls_short_of_the_observed_enso_band() {
    // **This test records a negative result.** The ticket's acceptance
    // criterion asks for a period in the observed ENSO range of roughly two to
    // seven years (`CONTEXT.md`, *ENSO*; the range the issue names). The model
    // oscillates, the oscillation has the threshold theory predicts, and its
    // period scales with `L/c` exactly as the delayed oscillator says it must
    // — but the coefficient puts it at about `5·L/c`, which for this basin is
    // near one year, below the band.
    //
    // The assertion is written the way the measurement came out rather than
    // the way the criterion asked, so that the discrepancy is in the suite and
    // not only in prose: physics that lengthens the period into the band will
    // fail this test and send whoever made the change to
    // `docs/enso-oscillation-report.md` § *Why the period is short*, which
    // gives the diagnosis and what would be needed to close the gap.
    //
    // Both bounds come from observation, not from this model: two years is the
    // fast edge of the observed ENSO band, and the basin crossing time `L/c`
    // is the shortest timescale the delayed oscillator can build a period out
    // of at all.
    const ENSO_BAND_FASTEST_YEARS: f64 = 2.0;

    let period_years = reference_record().dominant_period_years();
    let crossing_years =
        Experiment::pacific(SUPERCRITICAL_STRENGTH_PA_PER_K).basin_crossing_years();
    assert!(
        period_years > crossing_years,
        "the oscillation's period is {period_years} years, shorter than the {crossing_years} \
         years a Kelvin wave needs to cross the basin — no delayed oscillator can be that fast"
    );
    assert!(
        period_years < ENSO_BAND_FASTEST_YEARS,
        "the oscillation's period is {period_years} years, which is inside the observed ENSO \
         band this model was not able to reach; if that is right, \
         docs/enso-oscillation-report.md needs rewriting rather than this assertion relaxing"
    );
}

// ---------------------------------------------------------------------------
// The same claim, from what a run actually wrote to disk
// ---------------------------------------------------------------------------

/// The coupled scenario the written-run test drives, as the text of a scenario
/// file.
///
/// The same ocean, mixed layer, atmosphere and basin as [`Experiment`], stated
/// in the config format so that what is exercised is the whole path a user
/// takes: `[sst]` section, coupled solver, frame format, reader.
fn coupled_scenario_toml(total_steps: u64, output_every_n_steps: u64) -> String {
    format!(
        "[basin]\n\
         western_longitude_deg = {WESTERN_LONGITUDE_DEG}\n\
         eastern_longitude_deg = {PACIFIC_EASTERN_LONGITUDE_DEG}\n\
         southern_latitude_deg = {SOUTHERN_LATITUDE_DEG}\n\
         northern_latitude_deg = {NORTHERN_LATITUDE_DEG}\n\
         resolution_deg = {RESOLUTION_DEG}\n\
         \n\
         [physics]\n\
         reduced_gravity_m_per_s2 = {REDUCED_GRAVITY_M_PER_S2}\n\
         mean_thermocline_depth_m = {MEAN_THERMOCLINE_DEPTH_M}\n\
         rayleigh_damping_per_s = {RAYLEIGH_DAMPING_PER_S}\n\
         \n\
         [run]\n\
         dt_s = {DT_S}\n\
         total_steps = {total_steps}\n\
         output_every_n_steps = {output_every_n_steps}\n\
         \n\
         [sst]\n\
         mixed_layer_depth_m = {MIXED_LAYER_DEPTH_M}\n\
         mean_zonal_sst_gradient_k_per_m = {MEAN_ZONAL_SST_GRADIENT_K_PER_M}\n\
         subsurface_temperature_sensitivity_k_per_m = {SUBSURFACE_SENSITIVITY_K_PER_M}\n\
         thermal_damping_per_s = {THERMAL_DAMPING_PER_S}\n\
         wind_feedback_strength_pa_per_k = {SUPERCRITICAL_STRENGTH_PA_PER_K}\n\
         \n\
         [[wind]]\n\
         type = \"steady_trade_winds\"\n\
         equatorial_zonal_stress_pa = {TRADE_WIND_STRESS_PA}\n"
    )
}

#[test]
fn the_written_run_carries_the_oscillating_sst_index() {
    // The claim of this ticket has to be visible in what a run *writes*, not
    // only in a state a test holds in memory: T-05.4 put `T'` in the frame
    // format precisely so that an SST index could be read back off disk. So
    // this drives the whole path — a `[sst]` scenario, `run_scenario`, the
    // frame file, `RunReader` — and looks for the oscillation in the frames.
    //
    // Unlike the experiments above, nothing perturbs this run: it starts at
    // rest, and the mode grows out of the switch-on of the alizés. That is the
    // scenario a user would write.
    let directory = ScratchDir::new(TICKET, "written-oscillation");
    let total_steps = (RUN_YEARS * TROPICAL_YEAR_S / DT_S) as u64;
    let source = coupled_scenario_toml(total_steps, FRAME_EVERY_N_STEPS);
    let config = ScenarioConfig::from_toml(&source).expect("the coupled scenario parses");
    let scenario = config.build().expect("the coupled scenario is runnable");
    let run_directory = directory.path().join("run");
    run_scenario(&scenario, "enso-oscillation", &run_directory).expect("the coupled scenario runs");

    let reader = RunReader::open(&run_directory).expect("the run directory was written");
    let index = EasternSstIndex::over(scenario.basin());
    let nx = scenario.basin().grid().nx();
    let series_k: Vec<f64> = reader
        .map(|frame| {
            let frame = frame.expect("every frame of a run just written decodes");
            index.of_field(
                frame
                    .sst_anomaly_k()
                    .expect("a coupled run's frames carry `T'`"),
                nx,
            )
        })
        .collect();

    // The last third of the run, by which the mode has emerged from the
    // switch-on transient — the same window the in-memory experiments read,
    // for the same reason.
    let tail = &series_k[series_k.len() * 2 / 3..];
    let mean_k = tail.iter().sum::<f64>() / tail.len() as f64;
    let crossings = tail
        .windows(2)
        .filter(|pair| (pair[0] - mean_k) * (pair[1] - mean_k) < 0.0)
        .count();
    // Four crossings of the mean is two whole cycles: enough that what is seen
    // is an oscillation rather than a single overshoot on the way to a steady
    // state. The tail is 15 years long and the period is near one, so a
    // settled cycle crosses about thirty times.
    assert!(
        crossings >= 4,
        "the written run's SST index crossed its mean {crossings} times over the last third \
         of the run, which is not an oscillation"
    );
}

/// Years the written run integrates for.
///
/// Longer than the in-memory experiments' spin-up plus window, because this
/// one is not perturbed: the mode has to grow out of the switch-on transient
/// on its own.
const RUN_YEARS: f64 = 45.0;
/// Steps between the frames the written run saves.
///
/// About 0.19 years, which samples the near-annual oscillation five times a
/// cycle — enough to count its crossings — while keeping the run directory
/// under a hundred megabytes.
const FRAME_EVERY_N_STEPS: u64 = 200;
