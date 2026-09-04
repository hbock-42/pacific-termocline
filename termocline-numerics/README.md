# termocline-numerics

The finite-difference operators the solver runs on: `∂/∂x` and `∂/∂y` mapping
Arakawa C-grid cell centers to faces and back, plus the matching two-point
interpolations.

Where [`termocline-grid`](../termocline-grid) says *where* each variable sits,
this crate says how to differentiate and interpolate between those positions.
It is the metric half of the numerical core — the first place a cell width in
metres appears — and it is still physics-free: nothing here knows that `h` is a
thermocline depth anomaly, only that `h` lives at cell centers.

The staggering comes from
[ADR-0003](../docs/planning/adr/0003-numerical-scheme.md). Because a C-grid
face sits exactly midway between two cell centers and a cell center exactly
midway between two faces, every operator here is a centred stencil over one
cell width, and therefore second-order accurate.

Operators write through a caller-owned output field so that time stepping
allocates nothing per step.
