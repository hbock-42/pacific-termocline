//! Reporting a run's header to a terminal.
//!
//! `inspect` is the lightweight companion to `run`: it answers "what is in
//! this directory?" without a visualizer and without decoding a single frame.
//! Everything it prints comes from the run's JSON header, which [ADR-0004]
//! makes self-describing precisely so a reader never has to guess — grid,
//! physical parameters, scenario, variable list with units, and output
//! cadence.
//!
//! # Why only the header
//!
//! [`RunReader`](termocline_format::RunReader) decodes the header eagerly and
//! the frames lazily, so opening a run and never asking for a frame reads only
//! the header file. That is what keeps the command cheap on a run of any
//! length, and it is also what makes it useful on a run a crash cut short: the
//! header is on disk from the run's first moment, so a run whose frames are
//! missing or truncated still describes its scenario.
//!
//! # Why the numbers are printed the way they are
//!
//! Every quantity is printed in the SI unit the header records it in — no
//! conversion to days or kilometres, because a summary that silently rescaled
//! a parameter would misreport the run it is being used to check. Each value
//! is rendered in its shortest form that round-trips back to the same `f64`
//! (`{:?}` formatting), so what reaches the terminal names the number the run
//! was integrated with exactly, rather than a rounded neighbour of it.
//!
//! [ADR-0004]: ../../docs/planning/adr/0004-data-interchange-format.md

use std::fmt::Write as _;
use std::path::Path;

use termocline_format::{RunHeader, RunReadError, RunReader};

/// The header of the run in `directory` as a human-readable summary, the first
/// line naming the directory it was read from.
///
/// Only [`termocline_format::HEADER_FILE_NAME`] is read: the frames are never
/// decoded, so a run cut short still reports its scenario.
///
/// # Errors
/// The errors of [`RunReader::open`]: [`RunReadError::Open`] if the run's
/// header could not be opened, [`RunReadError::Header`] if it is not the JSON
/// a header is, and [`RunReadError::UnsupportedVersion`] if the run was
/// written by a format version this build does not read.
pub fn inspect_run(directory: &Path) -> Result<String, RunReadError> {
    let reader = RunReader::open(directory)?;
    Ok(format!(
        "run: {}\n{}",
        directory.display(),
        render_header(reader.header())
    ))
}

/// `header` as the lines [`inspect_run`] prints below the run's directory.
///
/// The rendering is a plain, line-per-field transcription rather than a table:
/// it is meant to be read in a terminal and grepped, and every field of
/// [`RunHeader`] appears exactly once.
#[must_use]
pub fn render_header(header: &RunHeader) -> String {
    let grid = header.grid;
    let extent = grid.extent();
    let params = header.physical_params;

    let mut out = String::new();
    // Writing into a `String` cannot fail, so the results are discarded rather
    // than propagated: there is no error here to report.
    let _ = writeln!(out, "format version: {}", header.format_version);
    let _ = writeln!(out, "scenario: {}", header.scenario_description);
    let _ = writeln!(out, "grid: {} x {} cells", grid.nx(), grid.ny());
    let _ = writeln!(
        out,
        "basin extent: {} to {} degrees east, {} to {} degrees north",
        number(extent.west_deg_east),
        number(extent.east_deg_east),
        number(extent.south_deg_north),
        number(extent.north_deg_north),
    );
    let _ = writeln!(
        out,
        "mean thermocline depth H = {} m",
        number(params.mean_depth_m)
    );
    let _ = writeln!(
        out,
        "reduced gravity g' = {} m s^-2",
        number(params.reduced_gravity_m_per_s2)
    );
    let _ = writeln!(out, "beta = {} m^-1 s^-1", number(params.beta_per_m_per_s));
    let _ = writeln!(
        out,
        "Rayleigh damping r = {} s^-1",
        number(params.rayleigh_damping_per_s)
    );
    let _ = writeln!(
        out,
        "reference density rho_0 = {} kg m^-3",
        number(params.reference_density_kg_per_m3)
    );
    let _ = writeln!(
        out,
        "frames: {}, one every {} s",
        header.output.frame_count,
        number(header.output.interval_s)
    );
    let variables = header
        .variables
        .iter()
        .map(|spec| format!("{} [{}]", spec.symbol, spec.unit))
        .collect::<Vec<_>>()
        .join(", ");
    let _ = writeln!(out, "variables: {variables}");
    out
}

/// `value` in its shortest text form that reads back as the same `f64`.
///
/// A run's parameters are the identity of the run, so the summary states them
/// exactly: `2.3e-11` rather than a fixed number of decimal places that would
/// print `β` as zero.
fn number(value: f64) -> String {
    format!("{value:?}")
}
