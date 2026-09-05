//! The acceptance criteria of T-05.4, at the level of a whole run: a coupled
//! scenario's mixed-layer SST anomaly `T'` reaches disk, comes back the way it
//! went in, and is announced by the header and by `inspect`; an uncoupled
//! scenario's frames say it has none.
//!
//! The format's own half of the ticket — the version 1 archive, the byte cost
//! of an absent field — is in `termocline-format/tests/sst_anomaly.rs`. What
//! is here is the join: the writer, the solver's state, the reader and the
//! `inspect` command, on runs the engine actually produced.
//!
//! There is no tolerance in this file. Writing a field is a transcription, not
//! an approximation, so the round trip is compared as exact IEEE-754 bit
//! patterns; a value that came back merely close would mean the run on disk is
//! not the run that was integrated.

use engine::{
    inspect_run, run_scenario, BetaPlane, Grid, OceanState, OutputSchedule, PhysicalParams,
    RunReader, RunWriteError, RunWriter, Scenario, Solver, Spacing, WindStressField,
    EQUATORIAL_BETA_PER_M_PER_S, H_STAGGERING, SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
};
use termocline_format::{BasinExtent, GridSpec, RunHeader, Variable};

mod common;

use common::ScratchDir;

/// This file's ticket, which labels the directories it leaves in the system
/// temp directory.
const TICKET: &str = "t054";

// --- A run written straight from a state -------------------------------

/// Reduced gravity `g'` of the equatorial Pacific's first baroclinic mode, in
/// m/s² (Gill, *Atmosphere–Ocean Dynamics*, ch. 11).
const REDUCED_GRAVITY_M_PER_S2: f64 = 0.05;
/// Mean thermocline depth `H`, in metres — the canonical 150 m upper layer of
/// the same 1.5-layer configuration.
const MEAN_DEPTH_M: f64 = 150.0;
/// Rayleigh damping `r`, in s⁻¹: an `e`-folding time of about 11.6 days.
const DAMPING_PER_S: f64 = 1.0e-6;

/// Cells along x and y of the hand-written basin. Different from one another
/// so a transposed field cannot pass.
const NX: usize = 6;
const NY: usize = 4;
/// Zonal and meridional extent of that basin, in metres — the order of the
/// equatorial Pacific's width and an equatorial channel reaching ±500 km.
const BASIN_LX_M: f64 = 1.0e7;
const BASIN_LY_M: f64 = 1.0e6;

/// Timestep of the hand-written run, in seconds: 15 minutes, well inside the
/// CFL bound of this basin.
const DT_S: f64 = 900.0;
const TOTAL_STEPS: u64 = 4;
const EVERY_N_STEPS: u64 = 2;

fn params() -> PhysicalParams {
    PhysicalParams::new(
        REDUCED_GRAVITY_M_PER_S2,
        MEAN_DEPTH_M,
        DAMPING_PER_S,
        EQUATORIAL_BETA_PER_M_PER_S,
        SEAWATER_REFERENCE_DENSITY_KG_PER_M3,
    )
    .expect("the Pacific parameter set is physical")
}

fn basin() -> (Grid, Spacing) {
    let grid = Grid::new(NX, NY).expect("the test basin has cells on both axes");
    let spacing = Spacing::new(BASIN_LX_M / NX as f64, BASIN_LY_M / NY as f64)
        .expect("a basin spanned by whole cells has positive spacing");
    (grid, spacing)
}

fn grid_spec() -> GridSpec {
    GridSpec::new(NX, NY, BasinExtent::new(120.0, -80.0, -5.0, 5.0))
        .expect("the test basin has cells on both axes")
}

fn header() -> RunHeader {
    let schedule = OutputSchedule::new(DT_S, TOTAL_STEPS, EVERY_N_STEPS)
        .expect("a positive timestep and a non-zero cadence are a valid schedule");
    RunHeader::new(
        grid_spec(),
        params().into(),
        "a hand-written coupled run",
        schedule.timing(),
    )
}

/// Values chosen to stress the encoder rather than to be physical: a negative
/// zero, a subnormal and the extremes of f64 all have bit patterns a lossy
/// encoder would not return unchanged.
fn awkward_values(len: usize) -> Vec<f64> {
    let stressors = [
        -0.0_f64,
        f64::MIN_POSITIVE / 2.0,
        f64::MAX,
        f64::MIN,
        f64::EPSILON,
        -1.234_567_890_123_456_7e-9,
    ];
    (0..len).map(|i| stressors[i % stressors.len()]).collect()
}

fn assert_bit_identical(context: &str, written: &[f64], read: &[f64]) {
    assert_eq!(written.len(), read.len(), "{context}: length changed");
    for (i, (w, r)) in written.iter().zip(read).enumerate() {
        assert_eq!(
            w.to_bits(),
            r.to_bits(),
            "{context}[{i}]: {w} came back as {r}"
        );
    }
}

#[test]
fn a_coupled_run_round_trips_the_anomaly_the_state_held_bit_for_bit() {
    // The ticket's first criterion, through the engine's own writer and
    // reader: the `T'` a coupled state holds is the `T'` a reader gets back,
    // to the bit.
    let (grid, _) = basin();
    let scratch = ScratchDir::new(TICKET, "round-trip");
    let header = header().with_sst_anomaly();

    let mut state = OceanState::at_rest_with_sst_anomaly(grid);
    let written =
        awkward_values(grid.field_shape(H_STAGGERING).0 * grid.field_shape(H_STAGGERING).1);
    state
        .sst_anomaly_k_mut()
        .expect("a coupled state carries the anomaly")
        .as_mut_slice()
        .copy_from_slice(&written);
    let wind = WindStressField::uniform(grid, -0.05, 0.02);

    let mut writer =
        RunWriter::create(scratch.path(), &header).expect("the scratch directory is writable");
    for frame in 0..header.output.frame_count {
        writer
            .append(frame as f64 * DT_S * EVERY_N_STEPS as f64, &state, &wind)
            .expect("the state covers the basin the header describes");
    }
    writer.finish().expect("the schedule's frames were written");

    let mut reader = RunReader::open(scratch.path()).expect("the run was written whole");
    assert!(reader.header().carries(Variable::SstAnomaly));
    let frames: Vec<_> = reader
        .by_ref()
        .collect::<Result<Vec<_>, _>>()
        .expect("every frame decodes");
    assert_eq!(frames.len() as u64, header.output.frame_count);
    for (index, frame) in frames.iter().enumerate() {
        let read = frame
            .sst_anomaly_k()
            .expect("a coupled run's frame carries its anomaly");
        assert_bit_identical(&format!("frame {index} sst"), &written, read);
        assert_eq!(
            frame.field(Variable::SstAnomaly),
            Some(read),
            "frame {index}: the variable and the named accessor are one field"
        );
    }
}

#[test]
fn a_state_that_disagrees_with_its_header_about_the_coupling_is_refused_by_name() {
    // The header's variable list is what a reader believes about every frame
    // beside it, so a run cannot promise a `T'` it does not integrate or write
    // one it never promised. Invalid input, so a Result rather than a panic
    // (CODING_STANDARDS.md).
    let (grid, _) = basin();
    let wind = WindStressField::uniform(grid, -0.05, 0.02);

    let mut sink = Vec::new();
    let mut writer = RunWriter::new(Vec::new(), &mut sink, &header())
        .expect("the header is writable into a buffer");
    let error = writer
        .append(0.0, &OceanState::at_rest_with_sst_anomaly(grid), &wind)
        .expect_err("an uncoupled header cannot carry a coupled state's anomaly");
    assert!(
        matches!(
            error,
            RunWriteError::SstAnomalyMismatch {
                header_carries: false
            }
        ),
        "{error:?}"
    );

    let mut sink = Vec::new();
    let mut writer = RunWriter::new(Vec::new(), &mut sink, &header().with_sst_anomaly())
        .expect("the header is writable into a buffer");
    let error = writer
        .append(0.0, &OceanState::at_rest(grid), &wind)
        .expect_err("a coupled header has no anomaly to write from an uncoupled state");
    assert!(
        matches!(
            error,
            RunWriteError::SstAnomalyMismatch {
                header_carries: true
            }
        ),
        "{error:?}"
    );
}

// --- Runs produced by a scenario ---------------------------------------

/// A scenario file over a coarse Pacific, with an optional `[sst]` section
/// appended. Two years of daily frames would be a slow test, so the run is
/// short; what is under test is what reaches disk, not what the ocean does.
fn scenario_toml(sst_section: &str) -> String {
    format!(
        "[basin]\n\
         resolution_deg = 2.0\n\
         \n\
         [physics]\n\
         reduced_gravity_m_per_s2 = 0.06\n\
         mean_thermocline_depth_m = 150.0\n\
         rayleigh_damping_per_s = 1e-7\n\
         \n\
         [run]\n\
         dt_s = 3600.0\n\
         total_steps = 48\n\
         output_every_n_steps = 12\n\
         \n\
         [[wind]]\n\
         type = \"steady_trade_winds\"\n\
         equatorial_zonal_stress_pa = -0.05\n\
         meridional_decay_scale_m = 361000.0\n\
         {sst_section}"
    )
}

/// The `[sst]` section of Zebiak & Cane (*Mon. Wea. Rev.* 115, 1987, § 2b):
/// a 50 m mixed layer, the equatorial Pacific's eastward cooling, γ = 0.1, and
/// a 125-day thermal relaxation.
const SST_SECTION: &str = "\n[sst]\n\
     mixed_layer_depth_m = 50.0\n\
     mean_zonal_sst_gradient_k_per_m = -4e-7\n\
     subsurface_temperature_sensitivity_k_per_m = 0.1\n\
     thermal_damping_per_s = 9.26e-8\n";

fn run_into(scratch: &ScratchDir, sst_section: &str) {
    let scenario = Scenario::from_toml(&scenario_toml(sst_section)).expect("a valid scenario");
    run_scenario(&scenario, "coarse trades", scratch.path()).expect("the run completes");
}

#[test]
fn a_coupled_scenario_writes_its_anomaly_into_every_frame() {
    // The deliverable: `T'` written to frames when the coupling is enabled,
    // and readable through `RunReader`. The field is one value per cell, and
    // it is not everywhere zero — the alizés upwell, so a run that reported
    // zeros would be reporting a field it did not integrate.
    let scratch = ScratchDir::new(TICKET, "coupled-scenario");
    run_into(&scratch, SST_SECTION);

    let mut reader = RunReader::open(scratch.path()).expect("the run was written whole");
    let header = reader.header().clone();
    assert!(header.carries(Variable::SstAnomaly));
    let cells = header.grid.field_len(Variable::SstAnomaly);
    assert_eq!(cells, header.grid.nx() * header.grid.ny());

    let frames: Vec<_> = reader
        .by_ref()
        .collect::<Result<Vec<_>, _>>()
        .expect("every frame decodes");
    assert_eq!(frames.len() as u64, header.output.frame_count);
    for (index, frame) in frames.iter().enumerate() {
        let sst = frame
            .sst_anomaly_k()
            .unwrap_or_else(|| panic!("frame {index} of a coupled run carries its anomaly"));
        assert_eq!(sst.len(), cells, "frame {index}");
        assert!(
            sst.iter().all(|value| value.is_finite()),
            "frame {index}: the anomaly must be a temperature, not a NaN"
        );
    }
    // Frame 0 is the initial condition, so its anomaly is identically zero;
    // the last frame has been forced for two days and cannot be.
    let first = frames.first().expect("the run wrote frames");
    assert!(first
        .sst_anomaly_k()
        .expect("frame 0 carries its anomaly")
        .iter()
        .all(|value| *value == 0.0));
    let last = frames.last().expect("the run wrote frames");
    assert!(
        last.sst_anomaly_k()
            .expect("the last frame carries its anomaly")
            .iter()
            .any(|value| *value != 0.0),
        "a forced coupled run must write an anomaly it actually integrated"
    );
}

#[test]
fn an_uncoupled_scenario_writes_no_anomaly_at_all() {
    // "An uncoupled run does not pay for a field it does not have", and says
    // so rather than writing a basin of zeros that a reader could not tell
    // from a measurement.
    let scratch = ScratchDir::new(TICKET, "uncoupled-scenario");
    run_into(&scratch, "");

    let mut reader = RunReader::open(scratch.path()).expect("the run was written whole");
    let header = reader.header().clone();
    assert!(!header.carries(Variable::SstAnomaly));
    let symbols: Vec<&str> = header.variables.iter().map(|v| v.symbol.as_str()).collect();
    assert_eq!(symbols, ["h", "u", "v", "tau_x", "tau_y"]);

    for (index, frame) in reader.by_ref().enumerate() {
        let frame = frame.expect("every frame decodes");
        assert_eq!(frame.sst_anomaly_k(), None, "frame {index}");
    }
}

#[test]
fn inspect_lists_the_anomaly_and_its_unit() {
    // The header is self-describing (ADR-0004) and `inspect` is what a human
    // reads it with, so the new variable appears there with its unit — kelvin,
    // because `T'` is a temperature difference — rather than only in the JSON.
    let coupled = ScratchDir::new(TICKET, "inspect-coupled");
    run_into(&coupled, SST_SECTION);
    let summary = inspect_run(coupled.path()).expect("the run was written whole");
    assert!(
        summary.contains(
            "variables: h [m], u [m s^-1], v [m s^-1], tau_x [N m^-2], tau_y [N m^-2], sst [K]"
        ),
        "{summary}"
    );

    let uncoupled = ScratchDir::new(TICKET, "inspect-uncoupled");
    run_into(&uncoupled, "");
    let summary = inspect_run(uncoupled.path()).expect("the run was written whole");
    assert!(
        summary
            .contains("variables: h [m], u [m s^-1], v [m s^-1], tau_x [N m^-2], tau_y [N m^-2]\n"),
        "{summary}"
    );
    assert!(!summary.contains("sst"), "{summary}");
}

#[test]
fn the_solver_and_the_header_agree_about_the_coupling_by_construction() {
    // The switch is read in one place: a scenario with an `[sst]` section
    // allocates the fourth field *and* declares it, and one without does
    // neither. A run cannot end up half-coupled.
    for (section, expected) in [("", false), (SST_SECTION, true)] {
        let scenario = Scenario::from_toml(&scenario_toml(section)).expect("a valid scenario");
        assert_eq!(scenario.sst_params().is_some(), expected);

        let scratch = ScratchDir::new(TICKET, "switch");
        run_scenario(&scenario, "switch", scratch.path()).expect("the run completes");
        let reader = RunReader::open(scratch.path()).expect("the run was written whole");
        assert_eq!(reader.header().carries(Variable::SstAnomaly), expected);
    }
}

/// A one-step coupled solver, kept so this file's hand-written run and the
/// scenario runs are checked against the same allocation rule: a coupled
/// solver gets a coupled state, and the writer refuses any other pairing.
#[test]
fn a_coupled_solver_asks_for_the_state_the_writer_will_accept() {
    let (grid, spacing) = basin();
    let plane = BetaPlane::centered_on_equator(params(), spacing, grid);
    let solver = Solver::new(grid, spacing, params(), plane, DT_S)
        .expect("a 15-minute step is inside this basin's CFL bound");
    assert!(!solver.couples_sst());
    assert!(!OceanState::at_rest(grid).couples_sst());
    assert!(OceanState::at_rest_with_sst_anomaly(grid).couples_sst());
}
