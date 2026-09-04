//! The variables a run writes, and what a reader needs to know about each.

use serde::{Deserialize, Serialize};
use termocline_grid::{Staggering, H_STAGGERING, U_STAGGERING, V_STAGGERING};

/// A field a run writes once per frame.
///
/// The set is closed: it is exactly the state and forcing of the 1.5-layer
/// model (`docs/planning/01-scientific-model.md`), and adding a member is a
/// format change, hence a [`crate::FORMAT_VERSION`] bump.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Variable {
    /// Thermocline depth anomaly `h` — a departure from the mean depth `H`,
    /// never a total depth.
    ThermoclineDepthAnomaly,
    /// Zonal (eastward) current anomaly `u`.
    ZonalCurrentAnomaly,
    /// Meridional (northward) current anomaly `v`.
    MeridionalCurrentAnomaly,
    /// Zonal component of the wind stress forcing, `τx`. Easterly trade-wind
    /// stress is negative.
    ZonalWindStress,
    /// Meridional component of the wind stress forcing, `τy`.
    MeridionalWindStress,
}

impl Variable {
    /// Every variable a frame carries, in the order the header lists them.
    pub const ALL: [Self; 5] = [
        Self::ThermoclineDepthAnomaly,
        Self::ZonalCurrentAnomaly,
        Self::MeridionalCurrentAnomaly,
        Self::ZonalWindStress,
        Self::MeridionalWindStress,
    ];

    /// The symbol this variable is written with in `CONTEXT.md`, unaccented so
    /// it doubles as the field name a non-Rust reader looks for.
    #[must_use]
    pub const fn symbol(self) -> &'static str {
        match self {
            Self::ThermoclineDepthAnomaly => "h",
            Self::ZonalCurrentAnomaly => "u",
            Self::MeridionalCurrentAnomaly => "v",
            Self::ZonalWindStress => "tau_x",
            Self::MeridionalWindStress => "tau_y",
        }
    }

    /// The SI unit of the values, so a reader never has to assume one.
    #[must_use]
    pub const fn unit(self) -> &'static str {
        match self {
            Self::ThermoclineDepthAnomaly => "m",
            Self::ZonalCurrentAnomaly | Self::MeridionalCurrentAnomaly => "m s^-1",
            Self::ZonalWindStress | Self::MeridionalWindStress => "N m^-2",
        }
    }

    /// Where this variable sits on the Arakawa C-grid of [ADR-0003].
    ///
    /// Each wind-stress component sits where the current it forces sits: `τx`
    /// enters the `u` equation, `τy` the `v` equation.
    ///
    /// [ADR-0003]: ../../docs/planning/adr/0003-numerical-scheme.md
    #[must_use]
    pub const fn staggering(self) -> Staggering {
        match self {
            Self::ThermoclineDepthAnomaly => H_STAGGERING,
            Self::ZonalCurrentAnomaly | Self::ZonalWindStress => U_STAGGERING,
            Self::MeridionalCurrentAnomaly | Self::MeridionalWindStress => V_STAGGERING,
        }
    }
}

/// One entry of the header's variable list: what a frame's field means, and in
/// what unit, spelled out for a reader that does not share these Rust types.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct VariableSpec {
    /// The variable this entry describes.
    pub variable: Variable,
    /// Its symbol, as written in `CONTEXT.md`.
    pub symbol: String,
    /// Its SI unit.
    pub unit: String,
}

impl VariableSpec {
    /// The header entry describing `variable`.
    #[must_use]
    pub fn of(variable: Variable) -> Self {
        Self {
            variable,
            symbol: variable.symbol().to_owned(),
            unit: variable.unit().to_owned(),
        }
    }
}
