# Backlog — Epics

Each epic file lists its MRs (merge-request-sized tickets) with a
description, deliverable, acceptance criteria, and dependencies. Rough
delivery order (later epics may reorder once the engine is minimally
running, but dependencies below are hard):

| Epic | Title | Depends on |
|---|---|---|
| [00](EPIC-00-project-foundations.md) | Project Foundations | — |
| [01](EPIC-01-numerical-core.md) | Numerical Core | 00 |
| [02](EPIC-02-ocean-dynamics.md) | Ocean Dynamics (Shallow-Water Core) | 01 |
| [03](EPIC-03-wind-forcing-alizes.md) | Wind Forcing (the Alizés) | 02 |
| [04](EPIC-04-boundaries-and-basin.md) | Boundaries and the Real Basin | 03 |
| [05](EPIC-05-data-io-serialization.md) | Data I/O & Serialization | 00, 02 |
| [06](EPIC-06-engine-cli-and-scenarios.md) | Engine CLI & Scenario Runner | 01, 03, 04, 05 |
| [07](EPIC-07-scientific-validation.md) | Scientific Validation | 04, 06 |
| [08](EPIC-08-visualizer-foundation.md) | Visualizer Foundation | 05, 06 |
| [09](EPIC-09-visualization-features.md) | Visualization Features | 08 |
| [10](EPIC-10-performance.md) | Performance | 06 |
| [11](EPIC-11-documentation.md) | Documentation | 06, 08 |
| [12](EPIC-12-sst-enso-coupling.md) | SST & ENSO Coupling (stretch) | 07 |

Epics 00–07 form the **engine v1** (a scientifically validated,
wind-forced, linear thermocline model). Epics 08–09 form the **visualizer
v1**. Epic 10 (performance) and Epic 11 (documentation) run alongside once
there's something working to profile/document. Epic 12 is an explicit
stretch goal, not part of the initial delivery target.
