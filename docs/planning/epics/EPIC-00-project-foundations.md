# Epic 00 — Project Foundations

## Goal
Stand up the repository structure, tooling, and CI so that every later epic
lands in a working, testable, lint-clean workspace from its first MR.

## Scope
Cargo workspace layout, shared crates, CI, licensing/contribution basics.
No physics, no rendering.

## Out of scope
Anything involving the equations, the grid, or rendering — those start in
Epic 01 / Epic 08.

---

### T-00.1: Cargo workspace skeleton
- **Description:** Create the root `Cargo.toml` workspace with member
  crates: `engine/` (binary + lib), `termocline-format/` (lib, per
  ADR-0004), `visualizer/` (binary + lib), and a placeholder `README.md` in
  each crate stating its purpose.
- **Deliverable:** `cargo build` and `cargo test` succeed on an empty
  workspace.
- **Acceptance criteria:**
  - Workspace builds with stable Rust, no crate has physics/UI code yet.
  - `.gitignore` covers `target/`, IDE files, OS cruft.
- **Depends on:** none.

### T-00.2: CI pipeline
- **Description:** GitHub Actions workflow running `cargo fmt --check`,
  `cargo clippy -- -D warnings`, and `cargo test --workspace` on every push
  and PR.
- **Deliverable:** `.github/workflows/ci.yml`.
- **Acceptance criteria:**
  - CI is green on the empty-but-building workspace from T-00.1.
  - Clippy warnings fail the build (not just format).
- **Depends on:** T-00.1.

### T-00.3: License and contribution basics
- **Description:** Add a license file (user to confirm which — default
  MIT/Apache-2.0 dual, standard for Rust projects, if no preference given)
  and a minimal `CONTRIBUTING.md` describing the epic/MR workflow this repo
  follows (link back to `docs/planning/`).
- **Deliverable:** `LICENSE-MIT`, `LICENSE-APACHE`, `CONTRIBUTING.md`.
- **Acceptance criteria:** files present, referenced from root `README.md`.
- **Depends on:** T-00.1.

### T-00.4: Shared numeric/grid types crate stub
- **Description:** Create `termocline-grid`, a small crate for the 2D
  field/grid types (`Field2D<T>`, grid indexing, C-grid staggering
  conventions per ADR-0003) with no physics — just data structures and
  indexing math, so both the engine and (later) tests can depend on a single
  definition of "what a grid cell is."
- **Deliverable:** `termocline-grid` crate with `Field2D` and grid-geometry
  types, unit-tested for indexing correctness only (no physics yet).
- **Acceptance criteria:**
  - Indexing round-trips correctly (`(i, j) -> flat -> (i, j)`).
  - C-grid staggering offsets for `h`, `u`, `v` are explicit named
    constants/types, not magic numbers, so Epic 01/02 code reads clearly.
- **Depends on:** T-00.1.
