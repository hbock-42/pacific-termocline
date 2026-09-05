//! T-11.2 — the scenario config reference documents the format the engine
//! actually accepts, and keeps doing so.
//!
//! `docs/scenario-config-reference.md` is prose about a `serde` record, and
//! prose does not fail to compile when a field is added, renamed or made
//! optional. These tests are the CI check the ticket asks for: they read
//! `engine/src/scenario.rs` as text, extract the fields the format really
//! defines, and compare that set against the fields the reference really
//! documents. A field on one side and not the other fails the build.
//!
//! Neither side is generated from the other. The source is the truth about
//! what the engine accepts; the document is a hand-written claim about it; the
//! test is the equality of the two, which is the only thing that can go stale
//! silently.

use std::collections::BTreeSet;

use engine::ScenarioConfig;

/// The format's definition, read as text rather than through `serde`, because
/// a `Deserialize` impl cannot be asked what field names it accepts.
const SOURCE: &str = include_str!("../src/scenario.rs");

/// The document under test.
const REFERENCE: &str = include_str!("../../docs/scenario-config-reference.md");

/// The four top-level items of the format. The `[[wind]]` entries are not
/// listed here: an entry's fields depend on its `type`, so the guard finds the
/// variants in the enum itself rather than trusting a list that a new forcing
/// would not update.
const SECTION_ITEMS: &[(&str, &str)] = &[
    ("pub struct ScenarioConfig {", "ScenarioConfig"),
    ("pub struct BasinSection {", "BasinSection"),
    ("pub struct PhysicsSection {", "PhysicsSection"),
    ("pub struct RunSection {", "RunSection"),
];

#[test]
fn every_field_of_the_format_is_documented() {
    for (marker, body, _) in documented_items() {
        let in_source = fields_of(&body);
        let in_reference: BTreeSet<String> = documented_rows(&marker)
            .into_iter()
            .map(|(name, _)| name)
            .collect();

        assert_eq!(
            in_source, in_reference,
            "docs/scenario-config-reference.md and engine/src/scenario.rs disagree about \
             `{marker}`: the source defines {in_source:?}, the reference documents \
             {in_reference:?}. Update the table under `<!-- fields: {marker} -->`."
        );
    }
}

#[test]
fn every_field_is_documented_as_required_or_optional_to_match_the_source() {
    for (marker, body, item_default) in documented_items() {
        for (name, row) in documented_rows(&marker) {
            // What makes a field omittable is `#[serde(default)]`, not its
            // type: an `Option<f64>` without one is still a key the file has
            // to carry. The attribute counts whether it sits on the field or
            // on the item — `#[serde(default)]` on a struct makes every one of
            // its fields omittable at once, which is how `[basin]` states a
            // default basin.
            let optional_in_source = item_default || is_optional(&body, &name);
            let documented_optional = row.to_lowercase().contains("optional");
            let documented_required = row.to_lowercase().contains("required");

            assert!(
                documented_optional ^ documented_required,
                "`{marker}.{name}` must be documented as exactly one of required or optional; \
                 its row reads: {row}"
            );
            assert_eq!(
                optional_in_source,
                documented_optional,
                "`{marker}.{name}` is documented as {}, but engine/src/scenario.rs declares it \
                 {} `#[serde(default)]`, which is what decides whether the file may omit it",
                if documented_optional {
                    "optional"
                } else {
                    "required"
                },
                if optional_in_source {
                    "with"
                } else {
                    "without"
                }
            );
        }
    }
}

#[test]
fn every_wind_type_tag_is_documented() {
    // The tag is what a scenario author types; `serde(rename_all =
    // "snake_case")` derives it from the variant name, so the reference has to
    // spell out the derived form rather than the Rust one.
    for tag in wind_variants().iter().map(|variant| snake_case(variant)) {
        assert!(
            REFERENCE.contains(&format!("type = \"{tag}\"")),
            "docs/scenario-config-reference.md never shows `type = \"{tag}\"`, so a reader \
             cannot discover the forcing"
        );
    }
}

#[test]
fn every_worked_example_in_the_reference_is_a_scenario_the_engine_accepts() {
    let examples = marked_toml_blocks();
    assert!(
        !examples.is_empty(),
        "docs/scenario-config-reference.md carries no `<!-- scenario -->` example, so nothing \
         proves its syntax is the syntax the engine parses"
    );
    for (index, example) in examples.iter().enumerate() {
        let config = ScenarioConfig::from_toml(example).unwrap_or_else(|error| {
            panic!("`<!-- scenario -->` example {index} does not parse: {error}\n{example}")
        });
        config.build().unwrap_or_else(|error| {
            panic!("`<!-- scenario -->` example {index} does not validate: {error}\n{example}")
        });
    }
}

/// Every item the reference has to carry a table for, as its marker and its
/// body in the source: the four sections, plus one per `WindSection` variant.
///
/// The variants are discovered rather than listed, so that adding a forcing to
/// the enum and forgetting to document it fails here — which is exactly the
/// staleness this guard exists to catch.
fn documented_items() -> Vec<(String, String, bool)> {
    let mut items: Vec<(String, String, bool)> = SECTION_ITEMS
        .iter()
        .map(|(header, marker)| {
            (
                (*marker).to_string(),
                item_body(SOURCE, header),
                item_carries_default(SOURCE, header),
            )
        })
        .collect();
    let wind = item_body(SOURCE, "pub enum WindSection {");
    let wind_default = item_carries_default(SOURCE, "pub enum WindSection {");
    for variant in wind_variants() {
        let body = item_body(&wind, &format!("{variant} {{"));
        items.push((format!("WindSection::{variant}"), body, wind_default));
    }
    items
}

/// Whether the item `header` opens carries `#[serde(default)]` on itself,
/// which makes every one of its fields omittable without any of them saying
/// so.
///
/// The attributes are the `#[…]` lines immediately above the declaration, so
/// this reads back from the header to the first line that is not one.
fn item_carries_default(source: &str, header: &str) -> bool {
    let start = source
        .find(header)
        .unwrap_or_else(|| panic!("engine/src/scenario.rs no longer declares `{header}`"));
    source[..start]
        .lines()
        .rev()
        .take_while(|line| line.trim_start().starts_with("#["))
        .any(|line| line.contains("serde(") && line.contains("default"))
}

/// The variants of `WindSection`, in declaration order: the `[[wind]]` types
/// the format defines.
fn wind_variants() -> Vec<String> {
    item_body(SOURCE, "pub enum WindSection {")
        .lines()
        .filter_map(|line| {
            let name = line.trim().strip_suffix(" {")?;
            let is_variant = name
                .chars()
                .next()
                .is_some_and(|first| first.is_ascii_uppercase())
                && name.chars().all(char::is_alphanumeric);
            is_variant.then(|| name.to_string())
        })
        .collect()
}

/// A variant name as `serde(rename_all = "snake_case")` writes it.
fn snake_case(variant: &str) -> String {
    let mut tag = String::new();
    for (index, character) in variant.char_indices() {
        if character.is_ascii_uppercase() && index > 0 {
            tag.push('_');
        }
        tag.push(character.to_ascii_lowercase());
    }
    tag
}

/// The field names one struct or enum variant of the format defines.
fn fields_of(body: &str) -> BTreeSet<String> {
    body.lines().filter_map(field_name).collect()
}

/// Whether the file may leave `name` out.
///
/// That is decided by `#[serde(default)]` on the declaration and by nothing
/// else — an `Option<f64>` without one is still a key the file has to carry,
/// and a plain `f64` with one is not.
fn is_optional(body: &str, name: &str) -> bool {
    let mut carries_default = false;
    for line in body.lines() {
        let line = line.trim();
        if let Some(field) = field_name(line) {
            if field == name {
                return carries_default;
            }
            carries_default = false;
        } else if line.starts_with("#[serde(") {
            carries_default = line.contains("default");
        } else if !line.starts_with("///") && !line.is_empty() {
            carries_default = false;
        }
    }
    panic!("no field `{name}` in this item")
}

/// The field a line declares, if it declares one.
///
/// A field line is `name: Type,` with an optional `pub` and optional
/// indentation; doc comments, attributes and variant headers are none of
/// those.
fn field_name(line: &str) -> Option<String> {
    let line = line.trim().strip_prefix("pub ").unwrap_or(line.trim());
    let (name, _) = line.split_once(':')?;
    let is_field = !name.is_empty()
        && name
            .chars()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_');
    is_field.then(|| name.to_string())
}

/// The text between the braces of the item `header` opens.
///
/// Brace-matched rather than indentation-matched, so a nested variant body is
/// found correctly. None of these items contains a brace inside a string
/// literal, which is the one case this would get wrong.
fn item_body(source: &str, header: &str) -> String {
    let start = source
        .find(header)
        .unwrap_or_else(|| panic!("engine/src/scenario.rs no longer declares `{header}`"));
    let open = start
        + source[start..]
            .find('{')
            .expect("an item header ends with an opening brace");
    let mut depth = 0usize;
    for (offset, character) in source[open..].char_indices() {
        match character {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    return source[open + 1..open + offset].to_string();
                }
            }
            _ => {}
        }
    }
    panic!("`{header}` is never closed");
}

/// The `| \`field\` | … |` rows of the table the reference marks with `marker`,
/// each as its field name and the whole row.
fn documented_rows(marker: &str) -> Vec<(String, String)> {
    let opening = format!("<!-- fields: {marker} -->");
    let start = REFERENCE.find(&opening).unwrap_or_else(|| {
        panic!("docs/scenario-config-reference.md carries no `{opening}` marker")
    });
    let region = &REFERENCE[start + opening.len()..];
    let end = region.find("<!-- end fields -->").unwrap_or_else(|| {
        panic!("the `{opening}` table is never closed by `<!-- end fields -->`")
    });

    region[..end]
        .lines()
        .filter_map(|line| {
            let cell = line.trim().strip_prefix("| `")?;
            let (name, _) = cell.split_once('`')?;
            Some((name.to_string(), line.trim().to_string()))
        })
        .collect()
}

/// Every fenced TOML block the reference marks as a complete scenario.
fn marked_toml_blocks() -> Vec<String> {
    REFERENCE
        .split("<!-- scenario -->")
        .skip(1)
        .filter_map(|rest| {
            let fence = rest.find("```toml")? + "```toml".len();
            let end = rest[fence..].find("```")? + fence;
            Some(rest[fence..end].to_string())
        })
        .collect()
}
