# Domain Docs

How the engineering skills should consume this repo's domain documentation when exploring the codebase.

## Before exploring, read these

- **`CONTEXT.md`** at the repo root: the domain glossary. Physics terms carry
  their symbols; `h` is a depth *anomaly*, not a total depth, and mistaking the
  two is the error this glossary exists to prevent.
- **`docs/planning/adr/`**: read ADRs that touch the area you're about to work in. This repo keeps its ADRs under `docs/planning/` alongside the vision/scope and scientific-model docs, not at the more common `docs/adr/` — don't create a second ADR directory at `docs/adr/`, and don't propose moving these.
- **`docs/planning/00-vision-and-scope.md`** and **`docs/planning/01-scientific-model.md`**: background on what's being built and the physics being modeled. Useful context even though they aren't `CONTEXT.md`/ADR-shaped.

`CONTEXT.md` is a glossary and nothing else — no implementation detail, no spec. Numerical-method vocabulary (C-grid, RK4, CFL) belongs to ADR-0003, not here. Extend it through `/domain-modeling` (reached via `/grill-with-docs` and `/improve-codebase-architecture`) as terms get resolved, pointing at `docs/planning/adr/` for ADRs, not the default `docs/adr/`.

## File structure

Single-context repo (this repo, and most repos):

```
/
├── CONTEXT.md                          ← domain glossary
├── CODING_STANDARDS.md                 ← read by /code-review's Standards axis
├── docs/planning/adr/                  ← existing ADRs live here
│   ├── 0001-engine-visualizer-split.md
│   ├── 0002-visualizer-language-choice.md
│   ├── 0003-numerical-scheme.md
│   └── 0004-data-interchange-format.md
└── engine/, visualizer/                ← not created yet; repo is in planning phase
```

This repo has no monorepo signals (no workspace config, single `engine/` + planned `visualizer/`), so multi-context layout (`CONTEXT-MAP.md` + per-context `CONTEXT.md`) doesn't apply. Revisit if that changes.

## Use the glossary's vocabulary

When your output names a domain concept (in an issue title, a refactor proposal, a hypothesis, a test name), use the term as defined in `CONTEXT.md`. Don't drift to synonyms the glossary explicitly avoids.

If the concept you need isn't in the glossary yet, that's a signal: either you're inventing language the project doesn't use (reconsider) or there's a real gap (note it for `/domain-modeling`).

## Flag ADR conflicts

If your output contradicts an existing ADR in `docs/planning/adr/`, surface it explicitly rather than silently overriding:

> _Contradicts ADR-0002 (visualizer language choice), but worth reopening because…_
